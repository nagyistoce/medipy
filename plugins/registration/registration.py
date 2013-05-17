##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import itk

import medipy.itk
from helpers import *

def registration_gui(fixed, moving, registration_type, optimization_type, optimization_strategy, nb_levels):
    """ Affine or Rigid registration
        
        <gui>
            <item name="fixed" type="Image" label="Fixed Image"/>
            <item name="moving" type="Image" label="Moving Image"/>
            <item name="registration_type" type="Enum" initializer="('Rigid', 'Affine')" label="Registration Type" />
            <item name="optimization_type" type="Enum" initializer="('Simplex', 'Gradient Descent')" label="Optimization Type" />
            <item name="optimization_strategy" type="Enum" initializer="('Normal', 'Multi Resolution')" label="Optimization Strategy" />
            <item name="nb_levels" type="Int" initializer="3" label="Number of Levels in the Pyramid" tooltip="Main parameter (1 &lt;= nb_levels &lt;= 5)"/>
            <item name="output" type="Image" initializer="output=True" role="return" label="Output"/>
        </gui>
    """

    if registration_type=="Rigid" :
        metric = mattes_mutual_information(fixed, moving,bins_count=64,pixels_fraction=0.5)
        transform = euler_3d_transform(True,fixed)
        scale = scales (transform)
        if optimization_type=="Gradient Descent" :
            optimizer = rsgd(transform, maximum_step_length=0.1, minimum_step_length=0.001, max_iterations=300, set_relaxation_factor=0.8, scales=scale)
        elif optimization_type=="Simplex" :
            optimizer = amoeba(transform,parameters_tolerance=0.1,function_tolerance=0.0001,scales=scale,initial_simplex_size=1e-4) 
        else :
            raise medipy.base.Exception("Wrong Optimization Type!")        
        interpolator = linear_interpolator(moving)
        if optimization_strategy=="Normal" :
            final_transform = registration(fixed, moving, metric, transform, optimizer, interpolator)
        elif optimization_strategy=="Multi Resolution" :
            final_transform = multi_resolution_registration(fixed, moving, metric, transform, optimizer, interpolator, levels=nb_levels)
        else :
            raise medipy.base.Exception("Wrong Optimization Strategy!")
        output = apply_transform(fixed, moving, final_transform, interpolator)

    elif registration_type=="Affine" :
        metric = mattes_mutual_information(fixed, moving)
        transform = affine_transform(True,fixed) 
        scale = scales (transform)
        if optimization_type=="Gradient Descent" :
            optimizer = rsgd(transform, maximum_step_length=0.2, minimum_step_length=0.001, max_iterations=300, set_relaxation_factor=0.8, scales=scale) 
        elif optimization_type=="Simplex" :
            optimizer = amoeba(transform,parameters_tolerance=0.1,function_tolerance=0.0001,scales=scale,initial_simplex_size=1e-4) 
        else :
            raise medipy.base.Exception("Wrong Optimization Type!")  
        interpolator = linear_interpolator(moving)
        if optimization_strategy=="Normal" :
            final_transform = registration(fixed, moving, metric, transform, optimizer, interpolator)
        elif optimization_strategy=="Multi Resolution" :
            final_transform = multi_resolution_registration(fixed, moving, metric, transform, optimizer, interpolator, levels=nb_levels)
        else :
            raise medipy.base.Exception("Wrong Optimization Strategy!")
        output = apply_transform(fixed, moving, final_transform, interpolator)

    else :
        raise medipy.base.Exception("Wrong Registration Type!")

    return output

def histogram_matching(fixed, moving, levels, match_points) :
    """ Histogram Matching
        
        <gui>
            <item name="fixed" type="Image" label="Fixed Image"/>
            <item name="moving" type="Image" label="Moving Image"/>
            <item name="levels" type="Int" initializer="1024" label="Number of Histogram Levels"/>
            <item name="match_points" type="Int" initializer="7" label="Number of Match Points"/>
            <item name="output" type="Image" initializer="output=True" role="return" label="Output"/>
        </gui>
    """

    fixed_itk = medipy.itk.medipy_image_to_itk_image(fixed, False)
    moving_itk = medipy.itk.medipy_image_to_itk_image(moving, False)

    matcher = itk.HistogramMatchingImageFilter[moving_itk, fixed_itk].New(
        Input=moving_itk, ReferenceImage=fixed_itk)
    matcher.SetNumberOfHistogramLevels(levels)
    matcher.SetNumberOfMatchPoints(match_points)
    matcher.ThresholdAtMeanIntensityOn()
    matcher()
    itk_output = matcher.GetOutput()
    output = medipy.itk.itk_image_to_medipy_image(itk_output, None, True)
    return output
  

def registration(fixed, moving, metric, transform, optimizer, interpolator):
    """ Compute the transformation between the fixed and the moving image.
    """
    
    fixed_itk = medipy.itk.medipy_image_to_itk_image(fixed, False)
    moving_itk = medipy.itk.medipy_image_to_itk_image(moving, False)
    
    registration = itk.ImageRegistrationMethod[fixed_itk, moving_itk].New(
        Metric=metric, Transform=transform, 
        InitialTransformParameters=list(transform.GetParameters()),
        Optimizer=optimizer, Interpolator=interpolator, 
        FixedImage=fixed_itk, MovingImage=moving_itk, 
        FixedImageRegion=fixed_itk.GetBufferedRegion())
    registration()
    
    final_transform = transform.__class__.New(
        Parameters=registration.GetLastTransformParameters(), 
        FixedParameters=transform.GetFixedParameters())
    
    return final_transform

def multi_resolution_registration(fixed, moving, metric, transform, optimizer, interpolator, levels):
    """ Compute the transformation in pyramids between the fixed and the moving image.
    """

    fixed_itk = medipy.itk.medipy_image_to_itk_image(fixed, False)
    moving_itk = medipy.itk.medipy_image_to_itk_image(moving, False)
    
    fixed_pyramid = itk.MultiResolutionPyramidImageFilter[fixed_itk, fixed_itk].New()
    moving_pyramid = itk.MultiResolutionPyramidImageFilter[moving_itk, moving_itk].New()

    registration = itk.MultiResolutionImageRegistrationMethod[fixed_itk, moving_itk].New(
        Metric=metric, Transform=transform, 
        InitialTransformParameters=list(transform.GetParameters()),
        Optimizer=optimizer, Interpolator=interpolator, 
        FixedImage=fixed_itk, MovingImage=moving_itk, 
        FixedImageRegion=fixed_itk.GetBufferedRegion(),
        FixedImagePyramid=fixed_pyramid, MovingImagePyramid=moving_pyramid,
        NumberOfLevels=levels)
    registration()
    
    print registration.GetLastTransformParameters()
    print transform.GetFixedParameters()
    final_transform = transform.__class__.New(
        Parameters=registration.GetLastTransformParameters(), 
        FixedParameters=transform.GetFixedParameters())
    
    return final_transform

def level_set_registration(fixed, moving, interpolator) :

    fixed_itk = medipy.itk.medipy_image_to_itk_image(fixed, False)
    moving_itk = medipy.itk.medipy_image_to_itk_image(moving, False)
    transformation_field_type = itk.Image[itk.Vector[itk.F,3],3]

    registration = itk.LevelSetMotionRegistrationFilter[fixed_itk, moving_itk, transformation_field_type].New()
    registration.SetFixedImage(fixed_itk)
    registration.SetMovingImage(moving_itk)
    registration.SetNumberOfIterations(300)
    registration.SetGradientSmoothingStandardDeviations(4)
    registration.Update()

    warper = itk.WarpImageFilter[moving_itk, moving_itk, transformation_field_type].New(
        Input=moving_itk,
        Interpolator=interpolator, OutputSpacing=fixed_itk.GetSpacing(),
        OutputOrigin=fixed_itk.GetOrigin(),DeformationField=registration.GetOutput())
    warper()

    itk_output = warper[0]
    output = medipy.itk.itk_image_to_medipy_image(itk_output, None, True)
    return output

    
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
