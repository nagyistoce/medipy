##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

from __future__ import division
import multiprocessing
import itertools
import numpy
from scipy import ndimage
import medipy.medimax.recalage
from medipy.diffusion.utils import gradient_itk, invert_itk, spectral_decomposition, compose_spectral, log_transformation, exp_transformation
import medipy.base
import medipy.io
from medipy.diffusion.tensor import ls_SecondOrderSymmetricTensorEstimation_




############
# Fiber
############

def interpolate_field(points,field):
    """
    Interpolate (linear) Image instance at arbitrary points in world space   
    The resampling is done with scipy.ndimage.
    """
    
    #points_voxel = (points-origin)/spacing
    points_voxel = points.T

    comps = []
    for i in range(3):
        im_comp = field[:,:,:,i]
        comp = ndimage.interpolation.map_coordinates(im_comp, points_voxel, order=0, prefilter=False)
        comps.insert(0,comp)

    return numpy.asarray(comps).T

def apply_transfo(points,transfo):
    """
    Apply the transformation to the fiber points.
    """
    #transfo = transfo*spacing
    #points x y z in mm
    return points + transfo


def register_multi(points,field,spacing_wrap,origin_wrap,spacing,origin):
    """
    Register a track.
    """
    points_voxel = (points-origin_wrap)/spacing_wrap # in index
    transfo = interpolate_field(points_voxel,field)
    points_wrap_voxel = apply_transfo(points_voxel,transfo) # in index
    return points_wrap_voxel*spacing + origin # in physical

def register_multi_(points,field,model_wrap,model_fix):
    """
    Register a track.
    """
    points_voxel = numpy.asarray( [model_wrap.physical_to_index(point[::-1])[::-1] for point in points] )  # in index
    transfo = interpolate_field(points_voxel,field)
    points_wrap_voxel = apply_transfo(points_voxel,transfo) # in index
    return numpy.asarray( [model_fix.index_to_physical(point[::-1])[::-1] for point in points_wrap_voxel] ) # in physical

def register_star(params):
    return register_multi(*params)

def register_fibers_wrong(T,ftrf,model_wrap,model,inv=False) :

    if inv :
        field_backward = load_trf(ftrf,model,athr=100.0)
    else :
        field_backward = load_trf(ftrf,model_wrap,athr=100.0)

    pool = multiprocessing.Pool() # Create a group of CPUs to run on
    asyncResult = pool.map_async(register_star, [a for a in itertools.izip( T,itertools.repeat(field_backward),itertools.repeat(model_wrap.spacing),itertools.repeat(model_wrap.origin),itertools.repeat(model.spacing),itertools.repeat(model.origin) ) ])
    resultList = asyncResult.get()

    return resultList


def register_fibers(fibers,ftrf,model1,model2):
    field_backward_inv = load_trf(ftrf,model1,athr=100.0)
    fout = [] 
    for f in fibers :
        fout.append(register_multi(f,field_backward_inv,model2.spacing,model2.origin,model1.spacing,model1.origin))
    return fout

def register_fibers_(fibers,field_backward_inv,model1,model2):
    fout = [] 
    for f in fibers :
        fout.append(register_multi(f,field_backward_inv,model2.spacing,model2.origin,model1.spacing,model1.origin))
    return fout


def invert_trf(fname_aff,fname_nl,fname_inv_aff,fname_inv_nl, fname_inv_affnl, model_wrap) :

    dpthref,hghtref,wdthref = model_wrap.shape
    dzref,dyref,dxref = model_wrap.spacing
    dzref = float(dzref)
    dyref = float(dyref)
    dxref = float(dxref)
    medipy.medimax.recalage.InvertTransfo3d(str(fname_aff), str(fname_inv_aff), wdthref, hghtref, dpthref, dxref, dyref, dzref, 0.01)
    medipy.medimax.recalage.InvertTransfo3d(str(fname_nl), str(fname_inv_nl), wdthref, hghtref, dpthref, dxref, dyref, dzref, 0.01)
    medipy.medimax.recalage.CombineTransfo3d(str(fname_inv_nl), str(fname_inv_aff), str(fname_inv_affnl), 5)

def invert_trf_rigid(fname_rigid,fname_inv_rigid, model_wrap) :

    dpthref,hghtref,wdthref = model_wrap.shape
    dzref,dyref,dxref = model_wrap.spacing
    dzref = float(dzref)
    dyref = float(dyref)
    dxref = float(dxref)
    medipy.medimax.recalage.InvertTransfo3d(str(fname_rigid), str(fname_inv_rigid), wdthref, hghtref, dpthref, dxref, dyref, dzref, 0.01)


#############
# Voxel
#############

def apply_tensor_trf_gui(*args,**kwargs) :
    """ Interpolation + PPD reorientation from a .trf deformation field

    <gui>
        <item name="fmodel_ref" type="File" label="Reference diffusion image"/>
        <item name="fmodel_wrap" type="File" label="Diffusion Image to wrap"/>
        <item name="ftrf" type="File" label="Deformation field"/>
        <item name="tensor_ref" type="Image" initializer="output=True" role="return" label="Reference tensor"/>
        <item name="tensor_registered" type="Image" initializer="output=True" role="return" label="Registered tensor"/>
    </gui>
    """
    fmodel_ref = kwargs['fmodel_ref']
    fmodel_wrap = kwargs['fmodel_wrap']
    ftrf = kwargs['ftrf']

    images_wrap = medipy.io.io.load_serie(fmodel_wrap)
    model_wrap = ls_SecondOrderSymmetricTensorEstimation_(images_wrap)

    images_ref = medipy.io.io.load_serie(fmodel_ref)
    model_ref = ls_SecondOrderSymmetricTensorEstimation_(images_ref)
    tensor_ref = model_ref

    tensor_registered = apply_tensor_trf(model_ref,model_wrap,ftrf)
 
    return tensor_ref,tensor_registered


def apply_tensor_trf(model_ref,model_wrap,ftrf) :
    """ Interpolation + PPD reorientation from a .trf deformation field
    """
    field_backward = load_trf(ftrf,model_wrap)
    log_model_wrap = log_transformation(model_wrap)
    log_model_out = interpolation_tensor_trf(log_model_wrap,model_ref,ftrf)
    model_out = exp_transformation(log_model_out)
    model_out = ppd_tensor_trf(model_out,field_backward)
    return model_out

def load_trf(ftrf,model_wrap,athr=500) :
    """ Load a .trf transformation field
    """
    ux = medipy.base.Image(shape=(1), dtype=numpy.float32)
    uy = medipy.base.Image(shape=(1), dtype=numpy.float32)
    uz = medipy.base.Image(shape=(1), dtype=numpy.float32)

    imSource = medipy.base.Image(shape=(1), spacing=model_wrap.spacing, origin=model_wrap.origin, direction=model_wrap.direction, dtype=numpy.float32)
    medipy.medimax.recalage.LoadTrfFile(str(ftrf),ux,uy,uz,imSource)

    field = medipy.base.Image(data=numpy.zeros(ux.shape+(3,),dtype=numpy.single),data_type="vector")
    field[...,0] = ux.data
    field[...,1] = uy.data
    field[...,2] = uz.data

    field.data[numpy.where(field.data>athr)] = 0.0

    return field


def ppd_tensor_trf(model_wrap,F) :
    """ Reorients a DT image using the Preservation of the Principal Direction (PPD) method.
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
    A = invert_itk(numpy.ascontiguousarray(A))

    # apply PPD
    return ppd(model_wrap,A)

def ppd(dt6,F):
    """ Implementation of Alexander's preservation of principle direction algorithm. 
    """

    imgEigVal,imgEigVec = spectral_decomposition(dt6)
    imgEigVec.data = imgEigVec.data.reshape(imgEigVec.shape+(3,3))

    # use F to find our n1, n2 
    dimDelV = F.shape

    # n1 = F*v1 
    n1 = numpy.zeros((dimDelV[0],dimDelV[1],dimDelV[2],3),dtype=numpy.single)
    n1[...,0] = F[...,0,0]*imgEigVec[...,2,0] + F[...,0,1]*imgEigVec[...,2,1] + F[...,0,2]*imgEigVec[...,2,2]
    n1[...,1] = F[...,1,0]*imgEigVec[...,2,0] + F[...,1,1]*imgEigVec[...,2,1] + F[...,1,2]*imgEigVec[...,2,2] 
    n1[...,2] = F[...,2,0]*imgEigVec[...,2,0] + F[...,2,1]*imgEigVec[...,2,1] + F[...,2,2]*imgEigVec[...,2,2] 
    # norm(n1)
    normN = numpy.sqrt(n1[...,0]**2 + n1[...,1]**2 + n1[...,2]**2)
    normN = normN + (normN == 0)
    n1[...,0] = n1[...,0]/normN
    n1[...,1] = n1[...,1]/normN
    n1[...,2] = n1[...,2]/normN

    # n2 = F*v2
    n2 = numpy.zeros((dimDelV[0],dimDelV[1],dimDelV[2],3),dtype=numpy.single)
    n2[...,0] = F[...,0,0]*imgEigVec[...,1,0] + F[...,0,1]*imgEigVec[...,1,1] + F[...,0,2]*imgEigVec[...,1,2]
    n2[...,1] = F[...,1,0]*imgEigVec[...,1,0] + F[...,1,1]*imgEigVec[...,1,1] + F[...,1,2]*imgEigVec[...,1,2] 
    n2[...,2] = F[...,2,0]*imgEigVec[...,1,0] + F[...,2,1]*imgEigVec[...,1,1] + F[...,2,2]*imgEigVec[...,1,2] 

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

    return medipy.base.Image(data=compose_spectral(imgEigVec,imgEigVal), spacing=dt6.spacing, origin=dt6.origin, direction=dt6.direction, data_type="vector", image_type="tensor_2")

     

def interpolation_tensor_trf(model,model_ref,ftrf):
    """ Linear interpolation of tensor model
    """
    nb_of_components = model._get_number_of_components()
    output = medipy.base.Image(data=numpy.zeros(model_ref.shape+(nb_of_components,),dtype=numpy.single), spacing=model_ref.spacing, origin=model_ref.origin, direction=model_ref.direction, data_type="vector", image_type="tensor_2")

    for i in range(nb_of_components):
        array_in = medipy.base.Image(data=numpy.ascontiguousarray(model[...,i]),spacing=model.spacing,origin=model.origin,direction=model.direction,data_type="scalar")
        array_out = medipy.base.Image(shape=(1),dtype=numpy.single,data_type="scalar")	
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
    g = gradient_itk(trf_field[...,0])
    J[...,0,0] = g[...,0] #df1/dx
    J[...,0,1] = g[...,1] #df1/dy
    J[...,0,2] = g[...,2] #df1/dz
    # compute Jy <- dy
    g = gradient_itk(trf_field[...,1])
    J[...,1,0] = g[...,0]
    J[...,1,1] = g[...,1]
    J[...,1,2] = g[...,2]
    # compute Jz <- dz
    g = gradient_itk(trf_field[...,2])
    J[...,2,0] = g[...,0]
    J[...,2,1] = g[...,1]
    J[...,2,2] = g[...,2]

    return J


if __name__ == '__main__':

    fdata1 = "/home/grigis/data/SOM/01/data_fsl_eddy.nii.gz"
    fdata2 = "/home/grigis/data/SOM/02/data_fsl_eddy.nii.gz"
    ftrf = "/home/grigis/data/SOM/affnl.trf"

    images2 = medipy.io.io.load_serie(fdata2)
    model2 = ls_SecondOrderSymmetricTensorEstimation_(images2)

    imgEigVal,imgEigVec = spectral_decomposition(model2)
    imgEigVec.data = imgEigVec.data.reshape(imgEigVec.shape+(3,3))
    t2 = compose_spectral(imgEigVec,imgEigVal)

    images1 = medipy.io.io.load_serie(fdata1)
    model1 = ls_SecondOrderSymmetricTensorEstimation_(images1)

    model2_registered = apply_tensor_trf(model1,model2,ftrf)

