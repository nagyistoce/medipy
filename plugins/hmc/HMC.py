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

import medipy.intensity
import medipy.logic
import medipy.segmentation

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

def post_process_outliers(segmentation, outliers, white_matter_class, min_size=6):
    """ Post-process the outliers by removing
        
        <gui>
            <item name="segmentation" type="Image" label="Segmentation" />
            <item name="outliers" type="Image" label="Outliers" />
            <item name="white_matter_class" type="Int" label="White matter class" />
            <item name="min_size" type="Float" label="Minimum size" initializer="6"/>
            <item name="result" type="Image" role="return"
                  initializer="output=True" label="Result"/>
        </gui>
    """
    
    # Only keep those in the white matter
    white_matter = medipy.intensity.binary_threshold(segmentation, 
        white_matter_class, white_matter_class, 1, 0)
    outliers = medipy.logic.apply_mask(outliers, white_matter, 0, 0)
    
    # Only keep those with more than 6 voxels
    connected_components = medipy.base.Image((1,1,1))
    medipy.segmentation.label_connected_components(outliers, connected_components)
    outliers = medipy.segmentation.order_connected_components(
        connected_components, min_size)
    
    return outliers
