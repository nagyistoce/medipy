##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################


""" Basic tools for DTI
"""

import numpy as np
import medipy.base

from spectral_analysis import spectral_analysis

def gradient_itk(im):
    """ Compute the gradient of an image
    """
    im_i = medipy.base.Image(data=np.ascontiguousarray(im), data_type="scalar")
    grad = medipy.base.Image(data=np.zeros(im.shape+(3,),dtype=np.single), data_type="vector")
    medipy.diffusion.gradient(im_i,grad)

    return grad

def invert_itk(M):
    """ Invert 3 by 3 Matrix at each voxel of M 
    """
    M = M.reshape(M.shape[:3]+(9,))
    M_i = medipy.base.Image(data=M, data_type="vector")
    medipy.diffusion.dtiInvMatrix(M_i)
    M = M_i.data.reshape(M.shape[:3]+(3,3))

    return M

def generate_image_sampling(image,step=(1,1,1)) :
    """ Generate seeds to init tractographu
    step is expressed in mm
    """

    spacing = image.spacing
    shape = image.shape*spacing
    Z,Y,X = np.mgrid[1:shape[0]-1:step[0], 1:shape[1]-1:step[1], 1:shape[2]-1:step[2]]
    X = X.flatten()
    Y = Y.flatten()
    Z = Z.flatten()
    seeds = []
    for i,j,k in zip(X,Y,Z) :
        seeds.append((i,j,k))
    return seeds

def length(xyz, constant_step=None):
    """ Euclidean length of track line in mm 
    """
    if constant_step==None :
        if xyz.shape[0] < 2 :
            return 0
        else :
            dists = np.sqrt((np.diff(xyz, axis=0)**2).sum(axis=1))
            return np.sum(dists)
    else :
        return (xyz.shape[0]-1)*constant_step



def voxel_parameters(tensor,w_size_plane=3,w_size_depth=3,mask=None):
    shape = tensor.shape
    mean = medipy.base.Image(data=np.zeros(shape[:3]+(6,),dtype=np.single),data_type="vector")
    var = medipy.base.Image(data=np.zeros(shape[:3]+(1,),dtype=np.single),data_type="vector")
    if mask==None :
        mask = medipy.base.Image(data=np.zeros((1,1,1),dtype=np.single))
        medipy.diffusion.dtiParamItk(tensor,mean,var,mask,w_size_plane,w_size_depth,False)
    else :
        medipy.diffusion.dtiParamItk(tensor,mean,var,mask,w_size_plane,w_size_depth,True)
    var.data = np.sqrt( var.data/6.0 ) # compute standard deviation

    return mean,var,w_size_plane*w_size_plane*w_size_depth


def spectral_decomposition(slice_tensor):
    shape = slice_tensor.shape

    eigVal = medipy.base.Image(data=np.zeros(shape+(3,),dtype=np.single),data_type="vector")
    eigVal_itk = medipy.itk.medipy_image_to_itk_image(eigVal, False)
    
    eigVec = medipy.base.Image(data=np.zeros(shape+(9,),dtype=np.single),data_type="vector")
    eigVec_itk = medipy.itk.medipy_image_to_itk_image(eigVec, False)
    
    itk_tensors = medipy.itk.medipy_image_to_itk_image(slice_tensor, False)
    spectral_analysis(itk_tensors,eigVal_itk,eigVec_itk)

    return eigVal,eigVec



def decompose_tensor(tensor):
    #outputs multiplicity as well so need to unique
    eigenvals, eigenvecs = np.linalg.eig(tensor)

    #need to sort the eigenvalues and associated eigenvectors
    order = eigenvals.argsort()[::-1]
    eigenvecs = eigenvecs[:, order]
    eigenvals = eigenvals[order]

    #Forcing negative eigenvalues to 0
    eigenvals = np.maximum(eigenvals, 0)

    return eigenvals, eigenvecs

def dti6to33(dt6):
    """ Full second order symmetric tensor from the six independent components.
        dt6 must be a numpy array, NOT a medipy.base.Image
    """
    
    dt33 = np.zeros(dt6.shape[:-1]+(3,3),dtype=dt6.dtype)
    dt33[...,0,0] = dt6[...,0]
    dt33[...,0,1] = dt6[...,1]
    dt33[...,0,2] = dt6[...,2]
    dt33[...,1,0] = dt6[...,1]
    dt33[...,1,1] = dt6[...,3]
    dt33[...,1,2] = dt6[...,4]
    dt33[...,2,0] = dt6[...,2]
    dt33[...,2,1] = dt6[...,4]
    dt33[...,2,2] = dt6[...,5]

    return dt33

def dti33to6(dt33):
    """ Six independent components from the full second order symmetric tensor.
        dt33 must be a numpy array, NOT a medipy.base.Image
    """
    
    dt6 = np.zeros(dt33.shape[:-2]+(6,),dtype=dt33.dtype)
    dt6[...,0] = dt33[...,0,0]
    dt6[...,1] = dt33[...,0,1]
    dt6[...,2] = dt33[...,0,2]
    dt6[...,3] = dt33[...,1,1]
    dt6[...,4] = dt33[...,1,2]
    dt6[...,5] = dt33[...,2,2]

    return dt6

def rotation33todt6(dt6,R):
    """ Apply a rotation R to a 3D numpy array (NOT a medipy.base.Image) 
        representing a second-order tensor by its six independant components.
    """

    dt6_r = np.zeros(dt6.shape,dtype=np.single)

    dt6_r[:,:,:,0] = R[0,0]*(R[0,0]*dt6[:,:,:,0] + R[0,1]*dt6[:,:,:,1] + R[0,2]*dt6[:,:,:,2])\
                   + R[0,1]*(R[0,0]*dt6[:,:,:,1] + R[0,1]*dt6[:,:,:,3] + R[0,2]*dt6[:,:,:,4])\
                   + R[0,2]*(R[0,0]*dt6[:,:,:,2] + R[0,1]*dt6[:,:,:,4] + R[0,2]*dt6[:,:,:,5])

    dt6_r[:,:,:,1] = R[1,0]*(R[0,0]*dt6[:,:,:,0] + R[0,1]*dt6[:,:,:,1] + R[0,2]*dt6[:,:,:,2])\
                   + R[1,1]*(R[0,0]*dt6[:,:,:,1] + R[0,1]*dt6[:,:,:,3] + R[0,2]*dt6[:,:,:,4])\
                   + R[1,2]*(R[0,0]*dt6[:,:,:,2] + R[0,1]*dt6[:,:,:,4] + R[0,2]*dt6[:,:,:,5])
     	
    dt6_r[:,:,:,2] = R[2,0]*(R[0,0]*dt6[:,:,:,0] + R[0,1]*dt6[:,:,:,1] + R[0,2]*dt6[:,:,:,2])\
                   + R[2,1]*(R[0,0]*dt6[:,:,:,1] + R[0,1]*dt6[:,:,:,3] + R[0,2]*dt6[:,:,:,4])\
                   + R[2,2]*(R[0,0]*dt6[:,:,:,2] + R[0,1]*dt6[:,:,:,4] + R[0,2]*dt6[:,:,:,5])

    dt6_r[:,:,:,3] = R[1,0]*(R[1,0]*dt6[:,:,:,0] + R[1,1]*dt6[:,:,:,1] + R[1,2]*dt6[:,:,:,2])\
                   + R[1,1]*(R[1,0]*dt6[:,:,:,1] + R[1,1]*dt6[:,:,:,3] + R[1,2]*dt6[:,:,:,4])\
                   + R[1,2]*(R[1,0]*dt6[:,:,:,2] + R[1,1]*dt6[:,:,:,4] + R[1,2]*dt6[:,:,:,5])

    dt6_r[:,:,:,4] = R[2,0]*(R[1,0]*dt6[:,:,:,0] + R[1,1]*dt6[:,:,:,1] + R[1,2]*dt6[:,:,:,2])\
                   + R[2,1]*(R[1,0]*dt6[:,:,:,1] + R[1,1]*dt6[:,:,:,3] + R[1,2]*dt6[:,:,:,4])\
                   + R[2,2]*(R[1,0]*dt6[:,:,:,2] + R[1,1]*dt6[:,:,:,4] + R[1,2]*dt6[:,:,:,5])

    dt6_r[:,:,:,5] = R[2,0]*(R[2,0]*dt6[:,:,:,0] + R[2,1]*dt6[:,:,:,1] + R[2,2]*dt6[:,:,:,2])\
                   + R[2,1]*(R[2,0]*dt6[:,:,:,1] + R[2,1]*dt6[:,:,:,3] + R[2,2]*dt6[:,:,:,4])\
                   + R[2,2]*(R[2,0]*dt6[:,:,:,2] + R[2,1]*dt6[:,:,:,4] + R[2,2]*dt6[:,:,:,5])

    return dt6_r


def compose_spectral(eigVec, eigVal):
    """ Recovers a DTI image in dt6 format [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz]
    """

    dz,dy,dx = eigVal.shape[:3]
    tensor = np.zeros((dz,dy,dx,6),dtype=np.single)

    tensor[:,:,:,0] = eigVec[:,:,:,2,0]*eigVal[:,:,:,2]*eigVec[:,:,:,2,0]\
                	+ eigVec[:,:,:,1,0]*eigVal[:,:,:,1]*eigVec[:,:,:,1,0]\
                	+ eigVec[:,:,:,0,0]*eigVal[:,:,:,0]*eigVec[:,:,:,0,0]
    tensor[:,:,:,3] = eigVec[:,:,:,2,1]*eigVal[:,:,:,2]*eigVec[:,:,:,2,1]\
                	+ eigVec[:,:,:,1,1]*eigVal[:,:,:,1]*eigVec[:,:,:,1,1]\
                	+ eigVec[:,:,:,0,1]*eigVal[:,:,:,0]*eigVec[:,:,:,0,1]
    tensor[:,:,:,5] = eigVec[:,:,:,2,2]*eigVal[:,:,:,2]*eigVec[:,:,:,2,2]\
                	+ eigVec[:,:,:,1,2]*eigVal[:,:,:,1]*eigVec[:,:,:,1,2]\
                	+ eigVec[:,:,:,0,2]*eigVal[:,:,:,0]*eigVec[:,:,:,0,2]
    tensor[:,:,:,1] = eigVec[:,:,:,2,0]*eigVal[:,:,:,2]*eigVec[:,:,:,2,1]\
                	+ eigVec[:,:,:,1,0]*eigVal[:,:,:,1]*eigVec[:,:,:,1,1]\
                	+ eigVec[:,:,:,0,0]*eigVal[:,:,:,0]*eigVec[:,:,:,0,1]
    tensor[:,:,:,2] = eigVec[:,:,:,2,0]*eigVal[:,:,:,2]*eigVec[:,:,:,2,2]\
                	+ eigVec[:,:,:,1,0]*eigVal[:,:,:,1]*eigVec[:,:,:,1,2]\
                	+ eigVec[:,:,:,0,0]*eigVal[:,:,:,0]*eigVec[:,:,:,2,2]
    tensor[:,:,:,4] = eigVec[:,:,:,2,1]*eigVal[:,:,:,2]*eigVec[:,:,:,2,2]\
               	    + eigVec[:,:,:,1,1]*eigVal[:,:,:,1]*eigVec[:,:,:,1,2]\
                	+ eigVec[:,:,:,0,1]*eigVal[:,:,:,0]*eigVec[:,:,:,0,2]	

    return tensor
