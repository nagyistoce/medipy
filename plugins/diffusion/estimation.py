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

def least_squares(images):
    """ Least Square Second Order Symmetric Tensor Estimation.
        A diffusion serie is composed of a float reference image (first element 
        in the list) and a set of float diffusion weighted images (on shell, 
        i.e. one bval).
        
        All images must have the same shape and the same dtype, and must 
        contain diffusion metadata.
    """
    
    # We're in the same package as itkSecondOrderSymmetricTensorReconstructionFilter, 
    # so it has already been included in itk by __init__
    
    PixelType = medipy.itk.dtype_to_itk[images[0].dtype.type]
    Dimension = images[0].ndim
    InputImage = itk.Image[PixelType, Dimension]
    OutputImage = itk.VectorImage[PixelType, Dimension]
    EstimationFilter = itk.SecondOrderSymmetricTensorReconstructionFilter[
        InputImage, OutputImage]
    
    estimation_filter = EstimationFilter.New()
    estimation_filter.SetBVal(images[1].metadata["mr_diffusion_sequence"][0].diffusion_bvalue)
    for cnt,image in enumerate(images) :
        itk_image = medipy.itk.medipy_image_to_itk_image(image, False)
        estimation_filter.SetInput(cnt,itk_image)
        
        gradient = image.metadata["mr_diffusion_sequence"][0].\
            diffusion_gradient_direction_sequence[0].diffusion_gradient_orientation
        estimation_filter.SetGradientDirection(cnt, [float(x) for x in gradient])
    
    itk_output = estimation_filter()[0]
    tensors = medipy.itk.itk_image_to_medipy_image(itk_output,None,True)
    tensors.image_type = "tensor_2"

    return tensors

