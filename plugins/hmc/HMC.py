##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import numpy

import itk

import medipy.itk

def segmentation(images, atlas, mask=None, flair_image=-1, iterations=5, 
                 display_outliers=True, outliers_criterion=0, threshold=0):
    
    ImageType = itk.Image[medipy.itk.dtype_to_itk[images[0].dtype.type], images[0].ndim]
    
    # Pad the inputs (images, atlas and mask) so that they have a cubic shape,
    # with a size equal to a power of 2.
    
    dimension = numpy.max([image.shape for image in images+atlas])
    if mask:
        dimension = max(dimension, numpy.max(mask.shape))
    padded_size = 2**int(numpy.ceil(numpy.log2(dimension)))
    
    padded_images = [_pad_image(image, padded_size) for image in images]
    padded_atlas = [_pad_image(image, padded_size) for image in atlas]
    padded_mask = None
    if mask:
        padded_mask = _pad_image(mask, padded_size)
    
    # Create the filter and run it.
    
    MaskType = ImageType
    if mask:
        MaskType = itk.Image[medipy.itk.dtype_to_itk[
            padded_mask.dtype.type], padded_mask.ndim]
    
    hmc_segmentation_filter = itk.HMCSegmentationImageFilter[
        ImageType, MaskType, ImageType].New(
        FlairImage=flair_image, Iterations=iterations, Modalities=len(images),
        DisplayOutliers=display_outliers, OutliersCriterion=outliers_criterion, 
        Threshold=threshold)
    
    inputs = []
    for (index, image) in enumerate(padded_images):
        inputs.append(medipy.itk.medipy_image_to_itk_image(image, False))
        hmc_segmentation_filter.SetInput(index, inputs[-1])
    for (index, image) in enumerate(padded_atlas):
        inputs.append(medipy.itk.medipy_image_to_itk_image(image, False))
        hmc_segmentation_filter.SetInput(len(images)+index, inputs[-1])

    padded_mask_itk = None
    if mask:
        print "setting mask"
        padded_mask_itk = medipy.itk.medipy_image_to_itk_image(
            padded_mask, False)
        hmc_segmentation_filter.SetMaskImage(padded_mask_itk)
    
    hmc_segmentation_filter()
    
    segmentation_itk = hmc_segmentation_filter.GetSegmentationImage()
    segmentation = medipy.itk.itk_image_to_medipy_image(segmentation_itk, None, True)
    
    outliers_itk = hmc_segmentation_filter.GetOutliersImage()
    outliers = medipy.itk.itk_image_to_medipy_image(outliers_itk, None, True)
    
    # Un-pad the output images
    for image in segmentation, outliers:
        image.data = image.data[[slice(0, x) for x in images[0].shape]]
    
    return (segmentation, outliers)

def _pad_image(image, padded_size):
    """ Pad the given image with zeros so that the resulting shape is 
        imge.ndim*(padded_size,)
    """
    
    padded_image = medipy.base.Image(
        shape=image.ndim*(padded_size,), dtype=image.dtype, value=0)
    
    padded_image[[slice(0, x) for x in image.shape]]=image.data
    
    return padded_image
