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
import os
import tempfile
import shutil

import medipy.io
import medipy.base
import medipy.segmentation
import medipy.fsl
import medipy.logic

import itkutils
import spectral_analysis

def get_diffusion_information(image) :
    """ Return the diffusion information from the ``image`` metadata. The 
        diffusion information is a dictionary and may contain the following keys :
        
        * ``"diffusion_bvalue"``
        * ``"diffusion_gradient_orientation"``
    """
    
    result = {}
    if "mr_diffusion_sequence" in image.metadata :
        mr_diffusion = image.metadata["mr_diffusion_sequence"][0]
        if "diffusion_bvalue" in mr_diffusion :
            result["diffusion_bvalue"] = mr_diffusion.diffusion_bvalue.value
        if "diffusion_gradient_direction_sequence" in mr_diffusion :
            diffusion_gradient_direction = mr_diffusion.\
                diffusion_gradient_direction_sequence.value[0]
            if "diffusion_gradient_orientation" in diffusion_gradient_direction :
                result["diffusion_gradient_orientation"] = \
                    diffusion_gradient_direction.diffusion_gradient_orientation.value
    
    return result

def baseline(images, b_threshold=0) :
    """ Return the baseline image

    <gui>
        <item name="images" type="ImageSerie" label="Input"/>
        <item name="b_threshold" type="Int" initializer="0" label="Value B for Baseline"/>
        <item name="output" type="Image" initializer="output=True" 
              role="return" label="Output"/>
    </gui>
        
    """
    b0 = None
    for image in images :
        result = get_diffusion_information(image)
        b_value = result["diffusion_bvalue"]
        if b_value <= b_threshold:
            b0 = image
            break
    
    if b0 == None:
        raise medipy.base.Exception("No reference image in diffusion sequence")
    
    return b0

def get_mask(images, b_threshold=0) :
    """ Return a mask computed from baseline image
    """
    b0 = baseline(images, b_threshold)
    
    skull_mask = medipy.segmentation.skull(b0)
    skull = medipy.logic.apply_mask(b0, skull_mask, 0, 0)
    
    directory = tempfile.mkdtemp()
    
    medipy.io.save(skull, os.path.join(directory, "skull.nii.gz"))
    bet = medipy.fsl.BET(
        os.path.join(directory, "skull.nii.gz"), os.path.join(directory, "bet.nii.gz"),
        intensity_threshold=0.2, create_brain_mask=True)
    bet()
    mask = medipy.io.load(bet.brain_mask)
    
    shutil.rmtree(directory)
    
    return mask

def log_transformation(dt6,epsi=1e-5) :
    """ Compute Log tensors
    """
    imgEigVal,imgEigVec = spectral_decomposition(dt6)
    imgEigVec.data = imgEigVec.data.reshape(imgEigVec.shape+(3,3))
    imgEigVal.data[np.where(imgEigVal.data<=epsi)] = epsi
    imgEigVal.data = np.log(imgEigVal.data)  
    return medipy.base.Image(data=compose_spectral(imgEigVec,imgEigVal),origin=dt6.origin,spacing=dt6.spacing,direction=dt6.direction,data_type="vector",image_type="tensor_2")

def exp_transformation(dt6) :
    """ Compute Exp tensor
    """
    imgEigVal,imgEigVec = spectral_decomposition(dt6)
    imgEigVec.data = imgEigVec.data.reshape(imgEigVec.shape+(3,3))
    imgEigVal.data = np.exp(imgEigVal.data)
    return medipy.base.Image(data=compose_spectral(imgEigVec,imgEigVal),origin=dt6.origin,spacing=dt6.spacing,direction=dt6.direction,data_type="vector",image_type="tensor_2")

def spectral_decomposition(image):
    """ Spectral decomposition of a DTI image. Return the eigenvalues and 
        eigenvectors of the tensor field.
    """
    
    shape = image.shape

    imgEigVal = medipy.base.Image(data=np.zeros(shape+(3,),dtype=image.dtype),
                               data_type="vector")
    imgEigVal_itk = medipy.itk.medipy_image_to_itk_image(imgEigVal, False)
    
    imgEigVec = medipy.base.Image(data=np.zeros(shape+(9,),dtype=image.dtype),
                               data_type="vector")
    imgEigVec_itk = medipy.itk.medipy_image_to_itk_image(imgEigVec, False)
    
    itk_tensors = medipy.itk.medipy_image_to_itk_image(image, False)
    spectral_analysis.spectral_analysis(itk_tensors,imgEigVal_itk,imgEigVec_itk)

    return imgEigVal,imgEigVec

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

def eigVec9to33(eigVec):
    """ Full eigen vectors from the 9 independent components.
        eigVec must be a numpy array, NOT a medipy.base.Image
        eigVec is composed of row vectors
    """
    
    eigVec33 = np.zeros(eigVec.shape[:-1]+(3,3),dtype=eigVec.dtype)
    eigVec33[...,0,0] = eigVec[...,0]   #V1x
    eigVec33[...,0,1] = eigVec[...,1]   #V1y
    eigVec33[...,0,2] = eigVec[...,2]   #V1z
    eigVec33[...,1,0] = eigVec[...,3]   #V2x
    eigVec33[...,1,1] = eigVec[...,4]   #V2y
    eigVec33[...,1,2] = eigVec[...,5]   #V2z
    eigVec33[...,2,0] = eigVec[...,6]   #V3x
    eigVec33[...,2,1] = eigVec[...,7]   #V3y
    eigVec33[...,2,2] = eigVec[...,8]   #V3z
    
    return eigVec33

def eigVec33to9(eigVec33):
    """ Nine independent components from the full eigen vectors.
        eigVec33 must be a numpy array, NOT a medipy.base.Image
        eigVec is composed of row vectors
    """
    
    eigVec = np.zeros(eigVec33.shape[:-2]+(9,),dtype=eigVec33.dtype)
    eigVec[...,0] = eigVec33[...,0,0]
    eigVec[...,1] = eigVec33[...,0,1]
    eigVec[...,2] = eigVec33[...,0,2]
    eigVec[...,3] = eigVec33[...,1,0]
    eigVec[...,4] = eigVec33[...,1,1]
    eigVec[...,5] = eigVec33[...,1,2]
    eigVec[...,6] = eigVec33[...,2,0]
    eigVec[...,7] = eigVec33[...,2,1]
    eigVec[...,8] = eigVec33[...,2,2]

    return eigVec

def rotation33todt6(dt6,R):
    """ Apply a rotation R to a 3D numpy array (NOT a medipy.base.Image) 
        representing a second-order tensor by its six independant components.
    """

    dt6_r = np.zeros(dt6.shape,dtype=np.single)

    dt6_r[...,0] = R[0,0]*(R[0,0]*dt6[...,0] + R[0,1]*dt6[...,1] + R[0,2]*dt6[...,2])\
                 + R[0,1]*(R[0,0]*dt6[...,1] + R[0,1]*dt6[...,3] + R[0,2]*dt6[...,4])\
                 + R[0,2]*(R[0,0]*dt6[...,2] + R[0,1]*dt6[...,4] + R[0,2]*dt6[...,5])

    dt6_r[...,1] = R[1,0]*(R[0,0]*dt6[...,0] + R[0,1]*dt6[...,1] + R[0,2]*dt6[...,2])\
                 + R[1,1]*(R[0,0]*dt6[...,1] + R[0,1]*dt6[...,3] + R[0,2]*dt6[...,4])\
                 + R[1,2]*(R[0,0]*dt6[...,2] + R[0,1]*dt6[...,4] + R[0,2]*dt6[...,5])
     	
    dt6_r[...,2] = R[2,0]*(R[0,0]*dt6[...,0] + R[0,1]*dt6[...,1] + R[0,2]*dt6[...,2])\
                 + R[2,1]*(R[0,0]*dt6[...,1] + R[0,1]*dt6[...,3] + R[0,2]*dt6[...,4])\
                 + R[2,2]*(R[0,0]*dt6[...,2] + R[0,1]*dt6[...,4] + R[0,2]*dt6[...,5])

    dt6_r[...,3] = R[1,0]*(R[1,0]*dt6[...,0] + R[1,1]*dt6[...,1] + R[1,2]*dt6[...,2])\
                 + R[1,1]*(R[1,0]*dt6[...,1] + R[1,1]*dt6[...,3] + R[1,2]*dt6[...,4])\
                 + R[1,2]*(R[1,0]*dt6[...,2] + R[1,1]*dt6[...,4] + R[1,2]*dt6[...,5])

    dt6_r[...,4] = R[2,0]*(R[1,0]*dt6[...,0] + R[1,1]*dt6[...,1] + R[1,2]*dt6[...,2])\
                 + R[2,1]*(R[1,0]*dt6[...,1] + R[1,1]*dt6[...,3] + R[1,2]*dt6[...,4])\
                 + R[2,2]*(R[1,0]*dt6[...,2] + R[1,1]*dt6[...,4] + R[1,2]*dt6[...,5])

    dt6_r[...,5] = R[2,0]*(R[2,0]*dt6[...,0] + R[2,1]*dt6[...,1] + R[2,2]*dt6[...,2])\
                 + R[2,1]*(R[2,0]*dt6[...,1] + R[2,1]*dt6[...,3] + R[2,2]*dt6[...,4])\
                 + R[2,2]*(R[2,0]*dt6[...,2] + R[2,1]*dt6[...,4] + R[2,2]*dt6[...,5])

    return dt6_r


def compose_spectral(imgEigVec, imgEigVal):
    """ Recovers the six independent components dt6 format [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz]
        (tensor is a numpy array, NOT a medipy.base.Image)
        imgEigVec and imgEigVal are medipy.base.Image
    """

    tensor = np.zeros(imgEigVal.shape+(6,),dtype=np.single)

    tensor[...,0] = imgEigVec[...,2,0]*imgEigVal[...,2]*imgEigVec[...,2,0]\
                  + imgEigVec[...,1,0]*imgEigVal[...,1]*imgEigVec[...,1,0]\
                  + imgEigVec[...,0,0]*imgEigVal[...,0]*imgEigVec[...,0,0]
    tensor[...,3] = imgEigVec[...,2,1]*imgEigVal[...,2]*imgEigVec[...,2,1]\
                  + imgEigVec[...,1,1]*imgEigVal[...,1]*imgEigVec[...,1,1]\
                  + imgEigVec[...,0,1]*imgEigVal[...,0]*imgEigVec[...,0,1]
    tensor[...,5] = imgEigVec[...,2,2]*imgEigVal[...,2]*imgEigVec[...,2,2]\
                  + imgEigVec[...,1,2]*imgEigVal[...,1]*imgEigVec[...,1,2]\
                  + imgEigVec[...,0,2]*imgEigVal[...,0]*imgEigVec[...,0,2]
    tensor[...,1] = imgEigVec[...,2,0]*imgEigVal[...,2]*imgEigVec[...,2,1]\
                  + imgEigVec[...,1,0]*imgEigVal[...,1]*imgEigVec[...,1,1]\
                  + imgEigVec[...,0,0]*imgEigVal[...,0]*imgEigVec[...,0,1]
    tensor[...,2] = imgEigVec[...,2,0]*imgEigVal[...,2]*imgEigVec[...,2,2]\
                  + imgEigVec[...,1,0]*imgEigVal[...,1]*imgEigVec[...,1,2]\
                  + imgEigVec[...,0,0]*imgEigVal[...,0]*imgEigVec[...,0,2]
    tensor[...,4] = imgEigVec[...,2,1]*imgEigVal[...,2]*imgEigVec[...,2,2]\
               	  + imgEigVec[...,1,1]*imgEigVal[...,1]*imgEigVec[...,1,2]\
                  + imgEigVec[...,0,1]*imgEigVal[...,0]*imgEigVec[...,0,2]	

    return tensor
