##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import itk
import medipy.itk
import numpy as np
from medipy.base import Image

def ls_SecondOrderSymmetricTensorEstimation(images):
    """ Least Square Second Order Symmetric Tensor Estimation.
    A diffusion serie is composed of a float reference image (first element in the list) and a set of float diffusion weighted images (on shell, i.e. one bval).

    <gui>
        <item name="images" type="Image" label="Input"/>
        <item name="output" type="Image" initializer="output=True" role="return" label="Output"/>
    </gui>
    """
    
    # We're in the same package as itkSecondOrderSymmetricTensorReconstructionFilter, so it has already been included in itk by __init__
    
    estimation_filter = itk.SecondOrderSymmetricTensorReconstructionFilter[itk.Image[itk.F,3], itk.VectorImage[itk.F,3]].New()
    estimation_filter.SetBVal(images[1].metadata["mr_diffusion_sequence"][0].diffusion_bvalue)
    for cnt,image in enumerate(images) :
        image.data = np.cast[np.single](image.data)
        itk_image = medipy.itk.medipy_image_to_itk_image(image, False)
        grad = image.metadata["mr_diffusion_sequence"][0].diffusion_gradient_direction_sequence[0].diffusion_gradient_orientation
        itk_grad = itk.Point[itk.F,3]()
        itk_grad[0] = float(grad[0])
        itk_grad[1] = float(grad[1])
        itk_grad[2] = float(grad[2])
        estimation_filter.SetInput(cnt,itk_image)
        estimation_filter.SetGradientDirection(cnt,itk_grad)
    itk_output = estimation_filter()[0]
    tensors = medipy.itk.itk_image_to_medipy_image(itk_output,None,True)
    tensors.image_type = "tensor_2"

    return tensors

