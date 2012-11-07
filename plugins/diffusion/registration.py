##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import numpy
import medipy.medimax.recalage
from utils import gradient_itk, invert_itk, spectral_decomposition, compose_spectral
import medipy.base

def apply_tensor_trf(model_ref,model_wrap,ftrf) :

    field_backward = load_trf(ftrf,model_wrap)
    model_out = interpolation_tensor_trf(model_wrap,model_ref,ftrf)
    model_out = ppd_tensor_trf(model_out,field_backward)
    return model_out

def load_trf(ftrf,model_wrap) :
    """ Load a .trf transformation field
    """
    ux = medipy.base.Image(shape=(1), dtype=numpy.float32)
    uy = medipy.base.Image(shape=(1), dtype=numpy.float32)
    uz = medipy.base.Image(shape=(1), dtype=numpy.float32)

    imSource = medipy.base.Image(shape=(1), spacing=model_wrap.spacing, origin=model_wrap.origin, direction=model_wrap.direction, dtype=numpy.float32)
    medipy.medimax.recalage.LoadTrfFile(ftrf,ux,uy,uz,imSource)

    field = numpy.zeros(ux.shape+(3,),dtype=numpy.single)
    field[:,:,:,0] = ux.data
    field[:,:,:,1] = uy.data
    field[:,:,:,2] = uz.data

    return field


def ppd_tensor_trf(model_wrap,F) :
    """ Reorients a DT image using the Preservation of the Principal Direction (PPD) method.
    """
    F[numpy.isnan(F)] = 0.0

    # approximate vector field F by local affine transformations A
    delV = jacobian_from_trf_field(F)
    identity = numpy.zeros(delV.shape,dtype=numpy.single)
    identity[:,:,:,0,0] = 1
    identity[:,:,:,1,1] = 1
    identity[:,:,:,2,2] = 1
    A = identity + delV

    # inverse backward field
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
    n1[:,:,:,0] = F[:,:,:,0,0]*imgEigVec[:,:,:,2,0] + F[:,:,:,0,1]*imgEigVec[:,:,:,2,1] + F[:,:,:,0,2]*imgEigVec[:,:,:,2,2]
    n1[:,:,:,1] = F[:,:,:,1,0]*imgEigVec[:,:,:,2,0] + F[:,:,:,1,1]*imgEigVec[:,:,:,2,1] + F[:,:,:,1,2]*imgEigVec[:,:,:,2,2] 
    n1[:,:,:,2] = F[:,:,:,2,0]*imgEigVec[:,:,:,2,0] + F[:,:,:,2,1]*imgEigVec[:,:,:,2,1] + F[:,:,:,2,2]*imgEigVec[:,:,:,2,2] 
    # norm(n1)
    normN = numpy.sqrt(n1[:,:,:,0]**2 + n1[:,:,:,1]**2 + n1[:,:,:,2]**2)
    normN = normN + (normN == 0)
    n1[:,:,:,0] = n1[:,:,:,0]/normN
    n1[:,:,:,1] = n1[:,:,:,1]/normN
    n1[:,:,:,2] = n1[:,:,:,2]/normN

    # n2 = F*v2
    n2 = numpy.zeros((dimDelV[0],dimDelV[1],dimDelV[2],3),dtype=numpy.single)
    n2[:,:,:,0] = F[:,:,:,0,0]*imgEigVec[:,:,:,1,0] + F[:,:,:,0,1]*imgEigVec[:,:,:,1,1] + F[:,:,:,0,2]*imgEigVec[:,:,:,1,2]
    n2[:,:,:,1] = F[:,:,:,1,0]*imgEigVec[:,:,:,1,0] + F[:,:,:,1,1]*imgEigVec[:,:,:,1,1] + F[:,:,:,1,2]*imgEigVec[:,:,:,1,2] 
    n2[:,:,:,2] = F[:,:,:,2,0]*imgEigVec[:,:,:,1,0] + F[:,:,:,2,1]*imgEigVec[:,:,:,1,1] + F[:,:,:,2,2]*imgEigVec[:,:,:,1,2] 

    # Projecting n2 onto n1-n2 plane: P(n2) = Pn2 = n2 - (n2*n1')*n1
    Pn2 = numpy.zeros((dimDelV[0],dimDelV[1],dimDelV[2],3),dtype=numpy.single)
    # temp = (n2*n1)*n1, so Pn2 = n2 - temp
    tempTemp = n1[:,:,:,0]*n2[:,:,:,0] + n1[:,:,:,1]*n2[:,:,:,1] + n1[:,:,:,2]*n2[:,:,:,2]
    Pn2 [:,:,:,0] = tempTemp * n1[:,:,:,0]
    Pn2 [:,:,:,1] = tempTemp * n1[:,:,:,1]
    Pn2 [:,:,:,2] = tempTemp * n1[:,:,:,2]
    Pn2 = n2 - Pn2
    normN = numpy.sqrt(Pn2[:,:,:,0]**2 + Pn2[:,:,:,1]**2 + Pn2[:,:,:,2]**2)
    normN = normN + (normN==0)
    Pn2[:,:,:,0] = Pn2[:,:,:,0]/normN
    Pn2[:,:,:,1] = Pn2[:,:,:,1]/normN
    Pn2[:,:,:,2] = Pn2[:,:,:,2]/normN

    # Computing n3 in order to have an orthogonal basis
    n3 = numpy.cross(n1,Pn2)
    normN = numpy.sqrt(n3[:,:,:,0]**2 + n3[:,:,:,1]**2 + n3[:,:,:,2]**2)
    normN = normN + (normN==0)
    n3[:,:,:,0] = n3[:,:,:,0]/normN
    n3[:,:,:,1] = n3[:,:,:,1]/normN
    n3[:,:,:,2] = n3[:,:,:,2]/normN

    imgEigVec[:,:,:,2,:] = n1
    imgEigVec[:,:,:,1,:] = Pn2
    imgEigVec[:,:,:,0,:] = n3

    return compose_spectral(imgEigVec,imgEigVal)

def interpolation_tensor_trf(model,model_ref,ftrf):
    """ Linear interpolation of tensor model
    """
    nb_of_components = model._get_number_of_components()
    output = medipy.base.Image(data=numpy.zeros(model_ref.shape+(nb_of_components,),dtype=numpy.single), spacing=model_ref.spacing, origin=model_ref.origin, direction=model_ref.direction, data_type="vector")

    for i in range(nb_of_components):
        array_in = medipy.base.Image(data=numpy.ascontiguousarray(model[:,:,:,i]),spacing=model.spacing,origin=model.origin,direction=model.direction,data_type="scalar")
        array_out = medipy.base.Image(shape=(1),dtype=numpy.single,data_type="scalar")	
        medipy.medimax.recalage.recalage.ApplyTransfo3d(array_in,str(ftrf),array_out,1)
        output[:,:,:,i] = array_out.data

    return output

def jacobian_from_trf_field(trf_field):
    """ Approximates Jacobian of a vector field - each voxel has an associated 3x3 jacobian
    INPUT:
	    trf_field	[Z,Y,X,3] 	    Deformation field (dx,dy,dz) of voxel displacements 
    OUTPUT:
	    J		    [Z,Y,X,3,3]	    Jacobian array 
    """

    dim = trf_field.shape
    J = numpy.zeros((dim[0],dim[1],dim[2],3,3),dtype=numpy.single)

    # compute Jx <- dx
    g = gradient_itk(trf_field[:,:,:,0])
    J[:,:,:,0,0] = g[:,:,:,0] #df/dx
    J[:,:,:,0,1] = g[:,:,:,1] #df/dy
    J[:,:,:,0,2] = g[:,:,:,2] #df/dz
    # compute Jy <- dy
    g = gradient_itk(trf_field[:,:,:,1])
    J[:,:,:,1,0] = g[:,:,:,0]
    J[:,:,:,1,1] = g[:,:,:,1]
    J[:,:,:,1,2] = g[:,:,:,2]
    # compute Jz <- dz
    g = gradient_itk(trf_field[:,:,:,2])
    J[:,:,:,2,0] = g[:,:,:,0]
    J[:,:,:,2,1] = g[:,:,:,1]
    J[:,:,:,2,2] = g[:,:,:,2]

    return J

if __name__ == '__main__':

    import medipy.io
    from medipy.diffusion.tensor import ls_SecondOrderSymmetricTensorEstimation_

    fdata1 = "/home/grigis/data/SOM/01/data_fsl_eddy.nii.gz"
    fdata2 = "/home/grigis/data/SOM/02/data_fsl_eddy.nii.gz"
    ftrf = "/home/grigis/data/SOM/affnl.trf"

    images2 = medipy.io.io.load_serie(fdata2)
    model2 = ls_SecondOrderSymmetricTensorEstimation_(images2)

    images1 = medipy.io.io.load_serie(fdata1)
    model1 = ls_SecondOrderSymmetricTensorEstimation_(images1)

    model2_registered = apply_tensor_trf(model1,model2,ftrf)

