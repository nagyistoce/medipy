##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import itertools

import itk
import numpy

import medipy.itk

def segmentation(images, atlas, mask=None, flair_image=-1, iterations=5, 
                 display_outliers=True, outliers_criterion=0, threshold=0):
    
    ImageType = itk.Image[medipy.itk.dtype_to_itk[images[0].dtype.type], images[0].ndim]
    
    # Create the filter and run it.
    
    MaskType = ImageType
    if mask:
        MaskType = itk.Image[medipy.itk.dtype_to_itk[mask.dtype.type], mask.ndim]
    
    hmc_segmentation_filter = itk.HMCSegmentationImageFilter[
        ImageType, MaskType, ImageType].New(
        FlairImage=flair_image, Iterations=iterations, Modalities=len(images),
        DisplayOutliers=display_outliers, OutliersCriterion=outliers_criterion, 
        Threshold=threshold)
        
    inputs = []
    for (index, image) in enumerate(itertools.chain(images, atlas)):
        itk_input = medipy.itk.medipy_image_to_itk_image(image, False)
        inputs.append(itk_input)
        hmc_segmentation_filter.SetInput(index, itk_input)
        
    mask_itk = None
    if mask:
        mask_itk = medipy.itk.medipy_image_to_itk_image(mask, False)
        hmc_segmentation_filter.SetMaskImage(mask_itk)
        
    hmc_segmentation_filter()
    
    segmentation_itk = hmc_segmentation_filter.GetSegmentationImage()
    segmentation = medipy.itk.itk_image_to_medipy_image(segmentation_itk, None, True)
    
    outliers_itk = hmc_segmentation_filter.GetOutliersImage()
    outliers = medipy.itk.itk_image_to_medipy_image(outliers_itk, None, True)
    
    return (segmentation, outliers)
