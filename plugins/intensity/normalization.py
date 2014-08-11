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
        image1, image2, mask=None, mask_value=1, 
        bins_count1=100, bins_count2=100, method=0):
    
    """ Intensity normalization based on joint histogram.
    """
    
    image1_itk = medipy.itk.medipy_image_to_itk_image(image1, False)
    image2_itk = medipy.itk.medipy_image_to_itk_image(image2, False)    
    
    mask_itk = None
    if mask:
        mask_itk = medipy.itk.medipy_image_to_itk_image(mask, False)

    if mask:
        FilterType = itk.JointHistogramNormalizationFilter[image1_itk, mask_itk, image2_itk]
    else:
        FilterType = itk.JointHistogramNormalizationFilter[image1_itk, image1_itk, image2_itk]
    
    filter_ = FilterType.New(image1_itk, image2_itk,
        MaskValue=mask_value, BinsCount1=bins_count1, BinsCount2=bins_count2)
    if mask:
        filter_.SetMask(mask_itk)
    
    # FIXME: should be in ctor
    if method == 0:
        filter_.SetMethodToNearestNeighbor() 
    elif method == 1:
        filter_.SetMethodToLinearInterpolation()
    else:
        raise medipy.base.Exception("Invalid estimation method: {0}".format(method))
    
    output_itk = filter_()[0]
    output = medipy.itk.itk_image_to_medipy_image(output_itk, None, True)
    
    return output
