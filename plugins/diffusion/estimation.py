##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import itk
import medipy.itk
import medipy.base

import numpy as np

def least_squares(limages, mask=None, return_baseline=False):
    """ Least Square Second Order Symmetric Tensor Estimation.
        <gui>
            <item name="limages" type="ImageSerie" label="Input"/>
            <item name="mask" type="Image" initializer="may_be_empty=True, may_be_empty_checked=True" 
                  label="Mask"/>
            <item name="output" type="Image" initializer="output=True" 
                  role="return" label="Output"/>
        </gui>
    """
    
    images = limages
    
    Dimension = images[0].ndim
    PixelType = medipy.itk.dtype_to_itk[images[0].dtype.type]
    
    ScalarImage = itk.Image[itk.F, Dimension]
    InputImage = itk.Image[PixelType, Dimension]
    TensorsImage = itk.VectorImage[PixelType, Dimension]
    BaselineImage = ScalarImage
    
    mask_itk = None
    MaskImage = ScalarImage
    if mask :
        mask_itk = medipy.itk.medipy_image_to_itk_image(mask, False)
        MaskImage = mask_itk.__class__
    
    EstimationFilter = itk.LeastSquaresImageFilter[
        InputImage, TensorsImage, BaselineImage, MaskImage]
    
    estimation_filter = EstimationFilter.New()
    
    if mask :
        estimation_filter.SetMaskImage(mask_itk)
        
    for cnt,image in enumerate(images) :
        itk_image = medipy.itk.medipy_image_to_itk_image(image, False)
        estimation_filter.SetInput(cnt,itk_image)
        
        gradient = image.metadata["mr_diffusion_sequence"][0].\
            diffusion_gradient_direction_sequence.value[0].diffusion_gradient_orientation.value
        estimation_filter.SetGradientDirection(cnt, [float(x) for x in gradient])
        
        b_value = image.metadata["mr_diffusion_sequence"][0].diffusion_bvalue.value
        
        estimation_filter.SetBvalueAndGradientDirection(cnt, float(b_value), [float(x) for x in gradient])
        
    estimation_filter()
    
    itk_output_tensors = estimation_filter.GetTensorsImage()
    tensors = medipy.itk.itk_image_to_medipy_image(itk_output_tensors,None,True)
    tensors.image_type = "tensor_2"

    itk_output_baseline = estimation_filter.GetBaselineImage()
    baseline = medipy.itk.itk_image_to_medipy_image(itk_output_baseline,None,True)
    baseline.image_type = "scalar"
    
    if return_baseline:
        return tensors, baseline
    else:
        return tensors

def weighted_least_squares(limages, mask=None, nb_iter=5, return_baseline=False):
    """ Weighted Least Square Second Order Symmetric Tensor Estimation.
        See publication : Hongtu ZHU "Statistical Analysis of Diffusion Tensors in Diffusion-Weighted Magnetic Resonance Imaging Data",
        Journal of the American Statistical Association, 102:480, 1085-1102, 2007.
        A diffusion serie is composed of a float reference image (first element 
        in the list) and a set of float diffusion weighted images (on shell, 
        i.e. one bval).
        
        All images must have the same shape and the same dtype, and must 
        contain diffusion metadata.
        
        <gui>
            <item name="limages" type="ImageSerie" label="Input"/>
            <item name="mask" type="Image" initializer="may_be_empty=True, may_be_empty_checked=True" 
                  label="Mask"/>
            <item name="nb_iter" type="Int" initializer="5" label="Iteration Count"
                tooltip="Number of iteration of the WLS estimation"/>
            <item name="output" type="Image" initializer="output=True" 
                  role="return" label="Output"/>
        </gui>
    """
    
    images = limages

    Dimension = images[0].ndim
    PixelType = medipy.itk.dtype_to_itk[images[0].dtype.type]
    
    ScalarImage = itk.Image[itk.F, Dimension]
    InputImage = itk.Image[PixelType, Dimension]
    TensorsImage = itk.VectorImage[PixelType, Dimension]
    BaselineImage = ScalarImage
    
    mask_itk = None
    MaskImage = ScalarImage
    if mask :
        mask_itk = medipy.itk.medipy_image_to_itk_image(mask, False)
        MaskImage = mask_itk.__class__
    
    EstimationFilter = itk.WeightedLeastSquaresImageFilter[
        InputImage, TensorsImage, BaselineImage, MaskImage]
    
    estimation_filter = EstimationFilter.New()
    estimation_filter.SetIterationCount(nb_iter)
    
    if mask :
        estimation_filter.SetMaskImage(mask_itk)
    
    for cnt,image in enumerate(images) :
        itk_image = medipy.itk.medipy_image_to_itk_image(image, False)
        estimation_filter.SetInput(cnt,itk_image)
        
        gradient = image.metadata["mr_diffusion_sequence"][0].\
            diffusion_gradient_direction_sequence.value[0].diffusion_gradient_orientation.value
        estimation_filter.SetGradientDirection(cnt, [float(x) for x in gradient])

        b_value = image.metadata["mr_diffusion_sequence"][0].diffusion_bvalue.value
        
        estimation_filter.SetBvalueAndGradientDirection(cnt, float(b_value), [float(x) for x in gradient])

    estimation_filter()
    
    itk_output_tensors = estimation_filter.GetTensorsImage()
    tensors = medipy.itk.itk_image_to_medipy_image(itk_output_tensors,None,True)
    tensors.image_type = "tensor_2"

    itk_output_baseline = estimation_filter.GetBaselineImage()
    baseline = medipy.itk.itk_image_to_medipy_image(itk_output_baseline,None,True)
    baseline.image_type = "scalar"
    
    if return_baseline:
        return tensors, baseline
    else:
        return tensors
