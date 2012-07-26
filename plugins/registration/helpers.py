import itk
import numpy

import medipy.base
import medipy.itk

def mattes_mutual_information(fixed, moving, bins_count=50, pixels_fraction=None,
                              **kwargs) :
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
    
    if pixels_fraction is None :
        if kwargs == {} :
            arguments["UseAllPixels"] = False
            arguments["NumberOfSpatialSamples"] = int(0.1*pixels_count)
        elif kwargs.get("all_pixels", False) :
            arguments["UseAllPixels"] = True
        elif "spatial_samples_count" in kwargs :
            arguments["UseAllPixels"] = False
            arguments["NumberOfSpatialSamples"] = kwargs["spatial_samples_count"]
        else :
            raise Exception("Incorrect arguments")
    else :
        arguments["UseAllPixels"] = False
        arguments["NumberOfSpatialSamples"] = int(pixels_fraction*pixels_count)
    
    metric = itk.MattesMutualInformationImageToImageMetric[
        FixedImageType, MovingImageType].New(**arguments)

    return metric

def euler_3d_transform(initialize=False, fixed=None, moving=None, moments=True) :
    """ Return an object of type itk.Euler3DTransform. If initialize is True,
        then the fixed and moving must be medipy.base.Image. The value of moments 
        indicates wether moments-based or geometry-based initialization is
        performed.
    """
    
    transform = itk.Euler3DTransform[itk.D].New()
    
    if initialize :
        fixed_itk = medipy.itk.medipy_image_to_itk_image(fixed, False)
        moving_itk = medipy.itk.medipy_image_to_itk_image(moving, False)
        
        initial_transform = itk.VersorRigid3DTransform[itk.D].New()
        initializer = itk.CenteredTransformInitializer[
            initial_transform, fixed_itk, moving_itk].New(
            Transform=initial_transform, FixedImage=fixed_itk, MovingImage=moving_itk) 
        if moments :
            initializer.MomentsOn()
        else :
            initializer.GeometryOn()
        initializer.InitializeTransform()
        transform.SetCenter(initial_transform.GetCenter())
    
    return transform

def amoeba(parameters_tolerance=0.25, function_tolerance=0.001, max_iterations=200,
           scales=None, initial_simplex_size=None, transform=None):
    """ Return an object of type itk.AmoebaOptimizer.
    
        If initial_simplex_size is not None, then transform must be specified
        to perform simplex initialization.
    """
    
    optimizer = itk.AmoebaOptimizer.New()
    
    optimizer.SetParametersConvergenceTolerance(parameters_tolerance)
    # Medimax <-> Numerical Recipes in C
    # recalage/mtch_3d.c:get_facteur_precision
    # NORMAL : 1
    # PRECIS : 0.1
    # TRES_PRECIS : 0.01
    # PRECISION_MAX : 0.0001
    optimizer.SetFunctionConvergenceTolerance(function_tolerance)
    
    optimizer.SetMaximumNumberOfIterations(max_iterations)
    
    if scales is not None :
        optimizer.SetScales(scales)
    
    if initial_simplex_size is not None :
        optimizer.AutomaticInitialSimplexOff()
        delta = transform.GetNumberOfParameters()*(initial_simplex_size,)
        optimizer.SetInitialSimplexDelta(delta)
    else :
        optimizer.AutomaticInitialSimplexOn()
    
    return optimizer

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