##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import itk

import medipy.itk

def registration(fixed, moving, metric, transform, optimizer, interpolator):
    """ Compute the transformation between the fixed and the moving image.
    """
    
    fixed_itk = medipy.itk.medipy_image_to_itk_image(fixed, False)
    moving_itk = medipy.itk.medipy_image_to_itk_image(moving, False)
    
    registration = itk.ImageRegistrationMethod[fixed_itk, moving_itk].New(
        Metric=metric, Transform=transform, 
        InitialTransformParameters=transform.GetParameters(),
        Optimizer=optimizer, Interpolator=interpolator, 
        FixedImage=fixed_itk, MovingImage=moving_itk, 
        FixedImageRegion=fixed_itk.GetBufferedRegion())
    registration()
    
    final_transform = transform.__class__.New(
        Parameters=registration.GetLastTransformParameters(), 
        FixedParameters=transform.GetFixedParameters())
    
    return final_transform

def apply_transform(fixed, moving, transform, interpolator) :
    
    fixed_itk = medipy.itk.medipy_image_to_itk_image(fixed, False)
    moving_itk = medipy.itk.medipy_image_to_itk_image(moving, False)
    
    resample = itk.ResampleImageFilter[moving_itk, fixed_itk].New(
        Transform=transform, Input=moving_itk, 
        Size=fixed_itk.GetLargestPossibleRegion().GetSize(), 
        OutputOrigin=fixed_itk.GetOrigin(), OutputSpacing=fixed_itk.GetSpacing(),
        OutputDirection=fixed_itk.GetDirection(), DefaultPixelValue=0, 
        Interpolator=interpolator)
    resample()
    
    itk_output = resample[0]
    output = medipy.itk.itk_image_to_medipy_image(itk_output, None, True)
    return output