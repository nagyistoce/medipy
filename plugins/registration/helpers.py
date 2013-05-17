##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import itk
import numpy

import medipy.base
import medipy.itk

def mutual_information_histogram(fixed, moving):
    if isinstance(fixed, medipy.base.Image) :
        FixedImageType = medipy.itk.itk_image_type(fixed)
    else :
        FixedImageType = fixed
    
    if isinstance(moving, medipy.base.Image) :
        MovingImageType = medipy.itk.itk_image_type(moving)
    else :
        MovingImageType = moving
        
    # Not wrapped !
    metric = itk.MutualInformationHistogramImageToImageMetric[FixedImageType, MovingImageType].New()

def mattes_mutual_information(fixed, moving, bins_count=24, pixels_fraction=None, **kwargs) :
    """ Return an object of type itk.MattesMutualInformation
        
        Arguments :
          * fixed : instance of itk.Image or medipy.base.Image to determine
            the first template parameter
          * moving : instance of itk.Image or medipy.base.Image to determine
            the second template parameter
          * bins_count : number of bins in the histogram. Defaults to 50 (cf.
            ITK Software guide, section 8.5.2
          * pixels_fraction : fraction of the total pixel count to use as number
            of spatial samples. Defaults to 0.1 (cf. ITK Software guide, section
            8.5.2)
          * all_pixels : if True, all pixels are used as spatial samples
          * spatial_samples_count : absolute number of pixels to be used as
            spatial samples
        
        Only one of pixels_fraction, all_pixels and spatial_samples_count may
        be specified.
    """
    
    pixels_count = 0
    if isinstance(fixed, medipy.base.Image) :
        pixels_count = numpy.prod(fixed.shape)
        FixedImageType = medipy.itk.itk_image_type(fixed)
    else :
        pixels_count = numpy.prod(fixed.GetLargestPossibleRegion().GetSize())
        FixedImageType = fixed
    
    if isinstance(moving, medipy.base.Image) :
        pixels_count = max(pixels_count, numpy.prod(moving.shape))
        MovingImageType = medipy.itk.itk_image_type(moving)
    else :
        pixels_count = max(pixels_count, numpy.prod(moving.GetLargestPossibleRegion().GetSize()))
        MovingImageType = moving

    arguments = {}
    
    # Medimax : recalage/imx_mtch.c:erreur_IM : histograms have 250 bins. We use
    # ITK suggested value
    arguments["NumberOfHistogramBins"] = int(bins_count)
    arguments["UseAllPixels"] = True

    if pixels_fraction is not None :
        arguments["NumberOfSpatialSamples"] = int(pixels_fraction*pixels_count)

    #if pixels_fraction is None :
    #    if kwargs == {} :
    #        arguments["UseAllPixels"] = False
    #        arguments["NumberOfSpatialSamples"] = int(0.1*pixels_count)
    #    elif kwargs.get("all_pixels", False) :
    #        arguments["UseAllPixels"] = True
   #     elif "spatial_samples_count" in kwargs :
    #        arguments["UseAllPixels"] = False
    #        arguments["NumberOfSpatialSamples"] = kwargs["spatial_samples_count"]
    #    else :
    #        raise Exception("Incorrect arguments")
    #else :
    #    arguments["UseAllPixels"] = False
    #    arguments["NumberOfSpatialSamples"] = int(pixels_fraction*pixels_count)
    
    metric = itk.MattesMutualInformationImageToImageMetric[
        FixedImageType, MovingImageType].New(**arguments)

    return metric

def euler_3d_transform(initialize=False, fixed=None, moments=True) :
    """ Return an object of type itk.Euler3DTransform. If initialize is True,
        then the fixed and moving must be medipy.base.Image. The value of moments 
        indicates wether moments-based or geometry-based initialization is
        performed.
    """
    
    transform = itk.Euler3DTransform[itk.D].New()
    
    if initialize :
        fixed_itk = medipy.itk.medipy_image_to_itk_image(fixed, False)
        #moving_itk = medipy.itk.medipy_image_to_itk_image(moving, False)

        fixedIndex = fixed_itk.GetLargestPossibleRegion().GetIndex()
        fixedSize = fixed_itk.GetLargestPossibleRegion().GetSize()

        centerIndex = (int(fixedIndex[0] + fixedSize[0] / 2.0), \
                       int(fixedIndex[1] + fixedSize[1] / 2.0), \
                       int(fixedIndex[2] + fixedSize[2] / 2.0))

        rotationCenter = fixed_itk.TransformIndexToPhysicalPoint(centerIndex)

        transform.SetIdentity()
        transform.SetCenter(rotationCenter)
        
        #initial_transform = itk.VersorRigid3DTransform[itk.D].New()
        #initializer = itk.CenteredTransformInitializer[
        #    initial_transform, fixed_itk, moving_itk].New(
        #    Transform=initial_transform, FixedImage=fixed_itk, MovingImage=moving_itk) 
        #if moments :
        #    initializer.MomentsOn()
        #else :
        #    initializer.GeometryOn()
        #initializer.InitializeTransform()
        #transform.SetCenter(initial_transform.GetCenter())
        #transform.SetOffset(initial_transform.GetOffset())
    
    return transform

def affine_transform(initialize=False, fixed=None, moments=True) :
    """ Return an object of type itk.AffineTransform. If initialize is True,
        then the fixed and moving must be medipy.base.Image. The value of moments 
        indicates wether moments-based or geometry-based initialization is
        performed.
    """
    
    transform = itk.AffineTransform[itk.D, fixed.ndim].New()
    
    if initialize :
        fixed_itk = medipy.itk.medipy_image_to_itk_image(fixed, False)
        #moving_itk = medipy.itk.medipy_image_to_itk_image(moving, False)

        fixedIndex = fixed_itk.GetLargestPossibleRegion().GetIndex()
        fixedSize = fixed_itk.GetLargestPossibleRegion().GetSize()

        centerIndex = (int(fixedIndex[0] + fixedSize[0] / 2.0), \
                       int(fixedIndex[1] + fixedSize[1] / 2.0), \
                       int(fixedIndex[2] + fixedSize[2] / 2.0))

        rotationCenter = fixed_itk.TransformIndexToPhysicalPoint(centerIndex)

        transform.SetIdentity()
        transform.SetCenter(rotationCenter)

        #initial_transform = itk.VersorRigid3DTransform[itk.D].New()
        #initializer = itk.CenteredTransformInitializer[
        #    initial_transform, fixed_itk, moving_itk].New(
        #    Transform=initial_transform, FixedImage=fixed_itk, MovingImage=moving_itk) 
        #if moments :
        #    initializer.MomentsOn()
        #else :
        #    initializer.GeometryOn()
        #initializer.InitializeTransform()
        #transform.SetCenter(initial_transform.GetCenter())
        #transform.SetOffset(initial_transform.GetOffset())
    
    return transform

def rsgd(transform, maximum_step_length=0.2, minimum_step_length=0.001, max_iterations=300, set_relaxation_factor=0.8, scales=None) :
    """ Retrun an object of type itk.RegularStepGradientDescentOptimizer
    """
    optimizer = itk.RegularStepGradientDescentOptimizer.New()
    optimizer.MinimizeOn()
    optimizer.SetMaximumStepLength(maximum_step_length)
    optimizer.SetMinimumStepLength(minimum_step_length)
    optimizer.SetNumberOfIterations( max_iterations )
    optimizer.SetRelaxationFactor(set_relaxation_factor)

    if scales is not None :
        optimizer.SetScales(scales)

    return optimizer

def amoeba(transform, parameters_tolerance=0.1, function_tolerance=0.0001, max_iterations=300, scales=None, initial_simplex_size=None):
    """ Return an object of type itk.AmoebaOptimizer.
    
        If initial_simplex_size is not None, then transform must be specified
        to perform simplex initialization.
    """

    #
    # Iteration Observer
    #
    #def iterationUpdate():
    #    print optimizer.GetInitialSimplexDelta()
    #    print transform.GetParameters()
   
    optimizer = itk.AmoebaOptimizer.New()
    optimizer.MinimizeOn()
    # Medimax <-> Numerical Recipes in C
    # recalage/mtch_3d.c:get_facteur_precision
    # NORMAL : 1
    # PRECIS : 0.1
    # TRES_PRECIS : 0.01
    # PRECISION_MAX : 0.0001
    optimizer.SetParametersConvergenceTolerance(parameters_tolerance) # 1/10th pixel
    optimizer.SetFunctionConvergenceTolerance(function_tolerance)  # 0.001 bits
    optimizer.SetMaximumNumberOfIterations(max_iterations)
       
    if initial_simplex_size is not None :
        optimizer.AutomaticInitialSimplexOff()
        delta = transform.GetNumberOfParameters()*(initial_simplex_size,) # the initial size of the simplex (initial_simplex_size units in each of the parameters)
        print delta
        optimizer.SetInitialSimplexDelta(delta)
    else :
        optimizer.AutomaticInitialSimplexOn()

    if scales is not None :
        optimizer.SetScales(scales)

    #iterationCommand = itk.PyCommand.New()
    #iterationCommand.SetCommandCallable( iterationUpdate )
    #optimizer.AddObserver( itk.IterationEvent(), iterationCommand.GetPointer() )
    
    return optimizer

def scales(transform) :

    optimizerScales = [0,]*transform.GetNumberOfParameters()

    if transform.GetNameOfClass()=="Euler3DTransform" :
        translationScale = 1.0 / 100.0
        optimizerScales[0] = 1.0
        optimizerScales[1] = 1.0
        optimizerScales[2] = 1.0
        optimizerScales[3] = 1.0
        optimizerScales[4] = translationScale
        optimizerScales[5] = translationScale

    elif transform.GetNameOfClass()=="AffineTransform" :
        translationScale = 1.0 / 1000.0
        optimizerScales[0] = 1.0
        optimizerScales[1] = 1.0
        optimizerScales[2] = 1.0
        optimizerScales[3] = 1.0
        optimizerScales[4] = 1.0
        optimizerScales[5] = 1.0
        optimizerScales[6] = 1.0
        optimizerScales[7] = 1.0
        optimizerScales[8] = 1.0
        optimizerScales[9] = translationScale
        optimizerScales[10] = translationScale
        optimizerScales[11] = translationScale

    else :
        raise medipy.base.Exception("Uknown Transformation Type!")

    return optimizerScales

def linear_interpolator(moving):
    """ Return an object of type itk.LinearInterpolateImageFunction
        
        Arguments :
          * moving : instance of itk.Image or medipy.base.Image to determine
            the first template parameter
    """
    
    if isinstance(moving, medipy.base.Image) :
        MovingImageType = medipy.itk.itk_image_type(moving)
    else :
        MovingImageType = moving
    
    return itk.LinearInterpolateImageFunction[MovingImageType, itk.D].New()
