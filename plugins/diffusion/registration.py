##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

from __future__ import division

import numpy
from scipy import ndimage

import medipy.base
from medipy.diffusion.utils import (spectral_decomposition, compose_spectral, 
    log_transformation, exp_transformation)
import medipy.medimax.recalage
import medipy.io

##########################
# Transformations fields #
##########################

def load_trf(filename, image, athr=500) :
    """ Load a .trf transformation field. Field values that are higher than
        ``athr`` are set to 0.
    """
    ux = medipy.base.Image((1), numpy.float32)
    uy = medipy.base.Image((1), numpy.float32)
    uz = medipy.base.Image((1), numpy.float32)
    imSource = medipy.base.Image((1), numpy.float32,
        image.spacing, image.origin, image.direction)
    medipy.medimax.recalage.LoadTrfFile(str(filename),ux,uy,uz,imSource)

    field = medipy.base.Image(
        ux.shape+(3,), numpy.single, data_type="vector")
    field[...,0] = ux.data
    field[...,1] = uy.data
    field[...,2] = uz.data

    field.data[numpy.where(field.data>athr)] = 0.0

    return field

def invert_trf(fname_aff,fname_nl,fname_inv_aff,fname_inv_nl, fname_inv_affnl, model_wrap) :

    dpthref,hghtref,wdthref = model_wrap.shape
    dzref,dyref,dxref = model_wrap.spacing
    dzref = float(dzref)
    dyref = float(dyref)
    dxref = float(dxref)
    medipy.medimax.recalage.InvertTransfo3d(str(fname_aff), str(fname_inv_aff), 
        wdthref, hghtref, dpthref, dxref, dyref, dzref, 0.01)
    medipy.medimax.recalage.InvertTransfo3d(str(fname_nl), str(fname_inv_nl), 
        wdthref, hghtref, dpthref, dxref, dyref, dzref, 0.01)
    medipy.medimax.recalage.CombineTransfo3d(
        str(fname_inv_nl), str(fname_inv_aff), str(fname_inv_affnl), 5)

def invert_trf_rigid(fname_rigid,fname_inv_rigid, model_wrap) :

    dpthref,hghtref,wdthref = model_wrap.shape
    dzref,dyref,dxref = model_wrap.spacing
    dzref = float(dzref)
    dyref = float(dyref)
    dxref = float(dxref)
    medipy.medimax.recalage.InvertTransfo3d(
        str(fname_rigid), str(fname_inv_rigid), 
        wdthref, hghtref, dpthref, dxref, dyref, dzref, 0.01)

############
# Fiber
############

def interpolate_field(points,field):
    """ Interpolate (linear) Image instance at arbitrary points in world space.
        The resampling is done with scipy.ndimage.
    """
    
    points_voxel = numpy.transpose(points)

    comps = []
    for i in range(3):
        im_comp = field[...,i]
        comp = ndimage.interpolation.map_coordinates(
            im_comp, points_voxel, order=0, prefilter=False)
        comps.insert(0,comp)

    return numpy.asarray(comps).T

def register_multi(points,field, moving, fixed):
    """
    Register a track.
    """
    
    indices_moving = [moving.physical_to_index(p) for p in points]
    transfo = interpolate_field(indices_moving,field)
    indices_fixed = indices_moving+transfo
    return [fixed.index_to_physical(i) for i in indices_fixed]

def register_fibers(fibers,ftrf,model1,model2):
    field_backward_inv = load_trf(ftrf,model1,100.0)
    fout = [] 
    for f in fibers :
        f_registered = register_multi(f,field_backward_inv, model2, model1)
        fout.append(f_registered)
    return fout

#############
# Voxel
#############

def apply_tensor_trf(model_ref,model_wrap,ftrf) :
    """ Apply the deformation field stored in file ``ftrf`` to ``model_wrap``,
        using the grid of ``model_ref``.
    
        <gui>
            <item name="model_ref" type="Image" label="Reference tensor image" />
            <item name="model_wrap" type="Image" label="Tensor image to wrap" />
            <item name="ftrf" type="File" label="Deformation field file" />
            <item name="registered" type="Image" initializer="output=True" 
                  role="return" label="Registered tensor image"/>
        </gui>
    """
    # apply transfo on tensor image
    field_backward = load_trf(ftrf,model_wrap)
    log_model_wrap = log_transformation(model_wrap)
    log_model_out = interpolation_tensor_trf(log_model_wrap,model_ref,ftrf)
    model_out = exp_transformation(log_model_out)
    model_out = ppd_tensor_trf(model_out,field_backward)
    
    # creation masks
    
    # MASK :take into account voxels out of the brain so sometimes initialized at 0
    mask_t = medipy.base.Image(shape=model_wrap.shape, dti="tensor_2", value=0, dtype=model_wrap.dtype)
    mask_t.copy_information(model_wrap)
    mask_t[numpy.nonzero(model_wrap)] = 1
    
    mask_out = medipy.base.Image(shape=model_ref.shape, data_type="scalar", image_type="normal", value=0, dtype=model_ref.dtype)
    mask_out.copy_information(model_ref)
    
    mask = medipy.base.Image(shape=model_wrap.shape, data_type="scalar", image_type="normal", value=0, dtype=model_wrap.dtype)
    mask.copy_information(model_wrap)
    mask[:,:,:] = mask_t[:,:,:,0]
    
    medipy.medimax.recalage.recalage_interface.ApplyTransfo3d_GUI(mask, ftrf, mask_out, 'Nearest')
    indices = numpy.nonzero(mask_out)
    
    # apply mask on registered tensor image
    model_registered = medipy.base.Image(shape=model_ref.shape, dti="tensor_2", value=0, dtype=model_ref.dtype)
    model_registered.copy_information(model_ref)
    
    model_registered[indices] = model_out[indices]
    
    return model_registered

def ppd_tensor_trf(model_wrap,F) :
    """ Return the transform of ``model_wrap`` by the field ``F`` using the 
        Preservation of the Principal Direction (PPD) method. 
    """
    F[numpy.isnan(F)] = 0.0

    # approximate vector field F by local affine transformations A
    delV = jacobian_from_trf_field(F)
    identity = numpy.zeros(delV.shape,dtype=numpy.single)
    identity[...,0,0] = 1
    identity[...,1,1] = 1
    identity[...,2,2] = 1
    A = identity + delV

    # inverse backward field
    # F: ref->wrap
    # delV: ref->wrap
    # I1(p1)-A->I2(p2)
    # T1(p1)-A-1->I1(p1)
    A = _invert_itk(numpy.ascontiguousarray(A))

    # apply PPD
    return ppd(model_wrap,A)

def ppd(dt6,F):
    """ Implementation of Alexander's preservation of principle direction 
        algorithm. 
    """

    imgEigVal,imgEigVec = spectral_decomposition(dt6)
    imgEigVec.data = imgEigVec.data.reshape(imgEigVec.shape+(3,3))

    # use F to find our n1, n2 
    dimDelV = F.shape

    # n1 = F*v1 
    n1 = numpy.zeros((dimDelV[0],dimDelV[1],dimDelV[2],3),dtype=numpy.single)
    n1[...,0] = numpy.sum(F[...,0,i]*imgEigVec[...,2,i] for i in range(3))
    n1[...,1] = numpy.sum(F[...,1,i]*imgEigVec[...,2,i] for i in range(3))
    n1[...,2] = numpy.sum(F[...,2,i]*imgEigVec[...,2,i] for i in range(3))
    # norm(n1)
    normN = numpy.sqrt(n1[...,0]**2 + n1[...,1]**2 + n1[...,2]**2)
    normN = normN + (normN == 0)
    n1[...,0] = n1[...,0]/normN
    n1[...,1] = n1[...,1]/normN
    n1[...,2] = n1[...,2]/normN

    # n2 = F*v2
    n2 = numpy.zeros((dimDelV[0],dimDelV[1],dimDelV[2],3),dtype=numpy.single)
    n2[...,0] = numpy.sum(F[...,0,i]*imgEigVec[...,1,i] for i in range(3))
    n2[...,1] = numpy.sum(F[...,1,i]*imgEigVec[...,1,i] for i in range(3))
    n2[...,2] = numpy.sum(F[...,2,i]*imgEigVec[...,1,i] for i in range(3))

    # Projecting n2 onto n1-n2 plane: P(n2) = Pn2 = n2 - (n2*n1')*n1
    Pn2 = numpy.zeros((dimDelV[0],dimDelV[1],dimDelV[2],3),dtype=numpy.single)
    # temp = (n2*n1)*n1, so Pn2 = n2 - temp
    tempTemp = n1[...,0]*n2[...,0] + n1[...,1]*n2[...,1] + n1[...,2]*n2[...,2]
    Pn2 [...,0] = tempTemp * n1[...,0]
    Pn2 [...,1] = tempTemp * n1[...,1]
    Pn2 [...,2] = tempTemp * n1[...,2]
    Pn2 = n2 - Pn2
    normN = numpy.sqrt(Pn2[...,0]**2 + Pn2[...,1]**2 + Pn2[...,2]**2)
    normN = normN + (normN==0)
    Pn2[...,0] = Pn2[...,0]/normN
    Pn2[...,1] = Pn2[...,1]/normN
    Pn2[...,2] = Pn2[...,2]/normN

    # Computing n3 in order to have an orthogonal basis
    n3 = numpy.cross(n1,Pn2)
    normN = numpy.sqrt(n3[...,0]**2 + n3[...,1]**2 + n3[...,2]**2)
    normN = normN + (normN==0)
    n3[...,0] = n3[...,0]/normN
    n3[...,1] = n3[...,1]/normN
    n3[...,2] = n3[...,2]/normN

    imgEigVec[...,2,:] = n1
    imgEigVec[...,1,:] = Pn2
    imgEigVec[...,0,:] = n3

    return medipy.base.Image(data=compose_spectral(imgEigVec,imgEigVal), 
        spacing=dt6.spacing, origin=dt6.origin, direction=dt6.direction, 
        data_type="vector", image_type="tensor_2")

def interpolation_tensor_trf(model,model_ref,ftrf):
    """ Linear interpolation of tensor model
    """
    nb_of_components = model._get_number_of_components()
    output = medipy.base.Image(model_ref.shape+(nb_of_components,), numpy.single,
        model_ref.spacing, model_ref.origin, model_ref.direction, 
        data_type="vector", image_type="tensor_2")

    for i in range(nb_of_components):
        array_in = medipy.base.Image(data=numpy.ascontiguousarray(model[...,i]),
            spacing=model.spacing, origin=model.origin, direction=model.direction,
            data_type="scalar")
        array_out = medipy.base.Image((1), numpy.single, data_type="scalar")	
        medipy.medimax.recalage.recalage.ApplyTransfo3d(array_in,str(ftrf),array_out,1)
        output[...,i] = array_out.data

    return output

def jacobian_from_trf_field(trf_field):
    """ Approximates Jacobian of a vector field - each voxel has an associated 3x3 jacobian
    INPUT:
	    trf_field	[Z,Y,X,3] 	    Deformation field (dx,dy,dz) of voxel displacements 
    OUTPUT:
	    J		    [Z,Y,X,3,3]	    Jacobian array 
    """

    J = numpy.zeros(trf_field.shape+(3,3),dtype=numpy.single)

    # compute Jx <- dx
    g = _gradient_itk(trf_field[...,0])
    J[...,0,0] = g[...,0] #df1/dx
    J[...,0,1] = g[...,1] #df1/dy
    J[...,0,2] = g[...,2] #df1/dz
    # compute Jy <- dy
    g = _gradient_itk(trf_field[...,1])
    J[...,1,0] = g[...,0]
    J[...,1,1] = g[...,1]
    J[...,1,2] = g[...,2]
    # compute Jz <- dz
    g = _gradient_itk(trf_field[...,2])
    J[...,2,0] = g[...,0]
    J[...,2,1] = g[...,1]
    J[...,2,2] = g[...,2]

    return J

def _gradient_itk(im):
    """ Compute the gradient of an image
    """
    im_i = medipy.base.Image(data=numpy.ascontiguousarray(im), data_type="scalar")
    grad = medipy.base.Image(im.shape+(3,), numpy.single, data_type="vector")
    medipy.diffusion.itkutils.gradient(im_i,grad)

    return grad

def _invert_itk(M):
    """ Invert 3 by 3 Matrix at each voxel of M 
    """
    M = M.reshape(M.shape[:3]+(9,))
    M_i = medipy.base.Image(data=M, data_type="vector")
    medipy.diffusion.itkutils.dtiInvMatrix(M_i)
    M = M_i.data.reshape(M.shape[:3]+(3,3))

    return M
