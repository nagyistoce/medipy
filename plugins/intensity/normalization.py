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
import medipy.base.exception
import medipy.itk
import medipy.logic

def mean_stdev_normalization(reference, image, mask_ref=None, mask_image=None):
    """ Return a normalized version of image, so that the mean and standard 
        deviation match those of reference.
        
        <gui>
            <item name="reference" type="Image" label="Reference"/>
            <item name="image" type="Image" label="Image"/>
            <item name="mask_ref" type="Image" initializer="may_be_empty=True" 
                  label="Mask of reference image" />
            <item name="mask_image" type="Image" initializer="may_be_empty=True" 
                  label="Mask of the image to normalize" />       
            <item name="output" type="Image" initializer="output=True"
                  role="return" label="Output"/>
        </gui>
    """
    
    if mask_ref :
        meanref=reference[numpy.nonzero(mask_ref)].mean()
        stdref=reference[numpy.nonzero(mask_ref)].std()
    else :
        meanref=reference.data.mean()
        stdref=reference.data.std()
        
    if mask_image :
        meanimage=image[numpy.nonzero(mask_image)].mean()
        stdimage=image[numpy.nonzero(mask_image)].std()
    else :
        meanimage=image.data.mean()
        stdimage=image.data.std()
        
    alpha = stdref/stdimage
    beta = meanref-meanimage*alpha
    
    data = alpha*image.data+beta
    output = medipy.base.Image(data=data)
    output.copy_information(image)
    
    return output

def one_parameter_linear_regression_normalization(src,ref):
    """ Return a normalized version of image, so that the mean and standard 
        deviation match those of reference.
        
        <gui>
            <item name="src" type="Image" label="Image to normalize"/>
            <item name="ref" type="Image" label="Reference"/>
            <item name="output" type="Image" initializer="output=True"
                role="return" label="Output"/>
        </gui>
    """
    
    if src.shape != ref.shape :
        raise medipy.base.exception.Exception("Images must have the same size") 
         
    alpha=numpy.sum(numpy.multiply(src,ref),dtype=float)/numpy.sum(numpy.multiply(src,src),dtype=float)
    data=alpha*src.data
    output = medipy.base.Image(data=data)
    output.copy_information(src)
    
    return output

def joint_histogram(
        fixed, moving, mask=None, mask_value=1, 
        bins_count_fixed=200, bins_count_moving=200, method="Nearest Neighbor"):
    
    """ Intensity normalization based on joint histogram.
    
        <gui>
            <item name="fixed" type="Image" label="Fixed"/>
            <item name="moving" type="Image" label="Moving"/>
            <item name="mask" type="Image" initializer="may_be_empty=True" 
                  label="Mask" />
            <item name="mask_value" type="Float" initializer="1" 
                  label="Mask value" />
            <item name="bins_count_fixed" type="Int" initializer="200"
                  label="Bins count (fixed)" />
            <item name="bins_count_moving" type="Int" initializer="200"
                  label="Bins count (moving)" />
            <item name="method" type="Enum" 
                  initializer='["Nearest Neighbor", "Linear"]'
                  label="Histogram interpolation" />
            <item name="output" type="Image" initializer="output=True"
                  role="return" label="Output"/>
        </gui>
    """
    
    fixed_itk = medipy.itk.medipy_image_to_itk_image(fixed, False)
    moving_itk = medipy.itk.medipy_image_to_itk_image(moving, False)    
    
    mask_itk = None
    if mask:
        mask_itk = medipy.itk.medipy_image_to_itk_image(mask, False)
    
    ##################################
    # 1. Compute the joint histogram #
    ##################################
    histogram_calculator = itk.JointHistogramCalculator[
        moving_itk, mask_itk or moving_itk].New(
        Image1=moving_itk, Image2=fixed_itk,
        BinsCount1=bins_count_moving, BinsCount2=bins_count_fixed)
    if mask:
        histogram_calculator.SetMask(mask_itk)
    
    # FIXME: should be in ctor
    if method == "Nearest Neighbor":
        histogram_calculator.SetMethodToNearestNeighbor() 
    elif method == "Linear":
        histogram_calculator.SetMethodToLinearInterpolation()
    
    histogram_calculator.Compute()
    histogram = histogram_calculator.GetHistogram()

    ####################################
    # 2. Compute the transfer function #
    ####################################
    transfer_function_calculator = itk.JointHistogramTransferFunctionCalculator[
        histogram].New(
        Histogram=histogram, ResampleFactor=10
    )
    transfer_function_calculator.Compute()
    transfer_function = transfer_function_calculator.GetTransferFunction()
    
    ########################################
    # 3. Adjust the moving image intensity #
    ########################################
    transform_intensity = itk.TransformIntensityImageFilter[
        moving_itk, moving_itk].New(
        Input=moving_itk, TransferFunction=transfer_function)
    output_itk = transform_intensity()[0]
    output = medipy.itk.itk_image_to_medipy_image(output_itk, None, True)
    
    return output
