##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import logging
import math
import sys

import itk
import numpy

from medipy.base import Image
import medipy.base
import medipy.itk

from medipy.io.dicom.dictionary import data_dictionary
from medipy.io.dicom.private_dictionaries import private_dictionaries
from medipy.io.dicom import Tag
import medipy.io.dicom.csa2
import medipy.io.dicom.normalize
import medipy.io.dicom.sort
import medipy.io.dicom.split
from vr import *

default_skipped_tags = set([
    Tag(0x0002,0x0003), # Media Storage SOP Instance UID
    Tag(0x0008,0x0018), # SOP Instance UID
    Tag(0x0020,0x0013), # Instance Number
    Tag(0x0020,0x0032), # Image Position (Patient)
    Tag(0x0020,0x1041), # Slice Location
    Tag(0x0020,0x9057), # In-Stack Position Number
    Tag(0x0028,0x0010), # Rows
    Tag(0x0028,0x0002), # Samples per Pixel
    Tag(0x0028,0x0011), # Columns
    Tag(0x0028,0x0106), # Smallest Image Pixel Value
    Tag(0x0028,0x0107), # Largest Image Pixel Value
    Tag(0x0028,0x1050), # Window Center
    Tag(0x0028,0x1051), # Window Center
])

def to_axis_aligned_space(image, destination):
    """ Transform the image to an axis-aligned frame of reference.
    """
    
    if destination is None:
        destination=medipy.base.coordinate_system.RAS
    
    source_aa = medipy.base.coordinate_system.best_fitting_axes_aligned_matrix(image.direction)
    if (source_aa==destination).all() :
        # No need to reorient
        return
    
    # Adjust the destination to the number of dimensions of image, add ones on
    # the diagonal if necessary
    destination_padded = medipy.base.array.reshape(destination,
        (image.ndim, image.ndim), "constant", False, value=0)
    for rank in range(image.ndim):
        if numpy.less_equal(max(medipy.base.coordinate_system.RAS.shape), rank).all():
            destination_padded[image.ndim-rank-1, image.ndim-rank-1] = 1
    
    # Convert from numpy to ITK coordinates (x,y,z) and structure
    destination_itk = numpy.fliplr(numpy.flipud(destination_padded))
    MatrixType = itk.Matrix[
        medipy.itk.dtype_to_itk[destination_itk.dtype.type], 
        destination_itk.shape[-1], destination_itk.shape[-2]]
    MatrixBridge = itk.MatrixBridge[MatrixType]
    destination_itk = MatrixBridge.GetMatrixFromArray(destination_itk)
    
    source = image.direction
    
    source_aa_itk = numpy.fliplr(numpy.flipud(source_aa))
    source_aa_itk = MatrixBridge.GetMatrixFromArray(source_aa_itk)
    
    itk_image = medipy.itk.medipy_image_to_itk_image(image, False)
    itk_image.SetDirection(source_aa_itk)
    
    orienter = itk.OrientImageFilter[itk_image, itk_image].New()
    orienter.UseImageDirectionOn()
    orienter.SetInput(itk_image)
    orienter.SetDesiredCoordinateDirection(destination_itk)
    orienter.Update()
    
    medipy.itk.itk_image_to_medipy_image(orienter.GetOutput(), image, True)
    
    # Restore real orientation : the ideal original direction D_o has been 
    # transformed to D_t : D_t = M.D_o, hence M = D_t.D_o^-1
    # Since directions are orthogonal matrices, we have :
    transformation = numpy.dot(image.direction, source_aa.T)
    image.direction = numpy.dot(transformation, source)

def to_axis_aligned_ras_space(image):
    """ Transform the image to the closest axis-aligned approximation of RAS
        (i.e. Nifti) space
    """
    
    return to_axis_aligned_space(image, medipy.base.coordinate_system.RAS)


def data(datasets):
    """ Build 3-D or 4-D data from datasets.
        
        All Data Sets must belong to the same stack and be geometrically sorted.
    """
    
    sample_dataset = datasets[0]
    
    if "MOSAIC" in sample_dataset.get("image_type", CS([])).value :
        image_csa = medipy.io.dicom.csa2.parse_csa(sample_dataset[0x0029,0x1010])
        number_of_tiles = image_csa["NumberOfImagesInMosaic"][0]
        mosaic_size = int(math.ceil(math.sqrt(number_of_tiles)))
        shape = (len(datasets), number_of_tiles, 
                 sample_dataset.rows.value/mosaic_size, sample_dataset.columns.value/mosaic_size)
    else :
        shape = (len(datasets),) + sample_dataset.pixel_array.shape[-2:]
    
    dtype = sample_dataset.pixel_array.dtype
    for dataset in datasets :
        if "rescale_slope" in dataset or "rescale_intercept" in dataset :
            dtype = numpy.float32
            break
    array = numpy.ndarray(shape, dtype=dtype)
    
    for index, dataset in enumerate(datasets) :
        # Get data
        if "MOSAIC" in dataset.get("image_type", CS([])).value :
            image_csa = medipy.io.dicom.csa2.parse_csa(dataset[0x0029,0x1010])
            number_of_tiles = image_csa["NumberOfImagesInMosaic"][0]
            mosaic_size = int(math.ceil(math.sqrt(number_of_tiles)))
            tile_size = (dataset.rows.value/mosaic_size, dataset.columns.value/mosaic_size)
            
            itk_pixel_array = medipy.itk.array_to_itk_image(dataset.pixel_array, False)
            
            TileImageType = itk.Image[medipy.itk.dtype_to_itk[dataset.pixel_array.dtype.type], 2]
            VolumeImageType = itk.Image[medipy.itk.dtype_to_itk[dataset.pixel_array.dtype.type], 3]
            assemble_tiles = itk.AssembleTilesImageFilter[TileImageType, VolumeImageType].New(
                Input=itk_pixel_array, 
                TileSize=tile_size, NumberOfTiles=number_of_tiles)
            itk_volume = assemble_tiles()[0]
            pixel_data = medipy.itk.itk_image_to_array(itk_volume, True)
        else : 
            pixel_data = dataset.pixel_array
        
        # Apply LUT transform
        if "rescale_slope" in dataset or "rescale_intercept" in dataset :
            rescale_slope = dataset.get("rescale_slope", DS(1.0)).value
            rescale_intercept = dataset.get("rescale_intercept", DS(0.0)).value
            pixel_data = pixel_data.astype(float)*rescale_slope+rescale_intercept
        
        # Insert into reconstructed data
        pixel_data = pixel_data.reshape(shape[1:])
        array[index] = pixel_data
    
    return array

def metadata(datasets, skipped_tags="default"):
    
    if skipped_tags == "default" :
        skipped_tags = default_skipped_tags
    
    special_processing = set([
        Tag(0x0020,0x0032), # Image Position (Patient) : keep only one
        Tag(0x0020,0x0037), # Image Orientation Patient : keep only one
        Tag(0x0028,0x0030), # Pixel Spacing : keep only one
        Tag(0x7fe0,0x0010), # Pixel Data : skip
    ])
    
    result = {}
    for dataset in datasets :
        for key, value in dataset.items() :
            value = value.value
            
            if key in skipped_tags or key in special_processing :
                continue
            
            if key.private :
                if key.element != 0x0010 :
                    private_creator = dataset.get(
                        (key.group, 0x0010), medipy.io.dicom.LO(None))
                    key = (private_creator, key)
                else :
                    # Private Creator : useless in here, skip it
                    continue
            
            if value not in result.setdefault(key, []) :
                result[key].append(value)
    
    # If all values are the same, replace list by value
    for key, value in result.items() :
        if len(result[key]) == 1 :
            result[key] = value[0]
    
    # Replace the tags in the keys by their names
    named_result = {}
    for key, value in result.items() :
        if isinstance(key, Tag) :
            if key in data_dictionary :
                name = data_dictionary[key][4]
            else :
                logging.warning("Unknown public tag : {0}".format(str(key)))
                name = str(key)
        else :
            private_creator, tag = key
            if private_creator.value in private_dictionaries :
                tag = "{0:04x}xx{1:02x}".format(tag.group, tag.element%0x100)
                name = private_dictionaries[private_creator.value].get(
                    tag, ("", "", "", "", str(tag)))[4]
            else :
                name = str(tag)
            # Make sure name is not empty
            name = name or str(tag)
        
        named_result[name] = value
    result = named_result
    
    # Origin
    origin = datasets[0].get("image_position_patient", FD((0,0,0))).value
    result["origin"] = tuple(reversed(origin))

    # Spacing
    spacing = datasets[0].get("pixel_spacing", FD((1.,1.))).value
        
    if len(datasets) >= 2 :
        slice_spacing = numpy.linalg.norm(numpy.subtract(
            datasets[0].get("image_position_patient", FD((0,0,0))).value, 
            datasets[1].get("image_position_patient", FD((0,0,0))).value))
    else :
        # Sane default
        slice_spacing = 1.0
    
    if slice_spacing == 0 :
        slice_spacing = 1.
    result["spacing"] = (slice_spacing,)+tuple(reversed(spacing))
    
    # Orientation
    orientation = datasets[0].get("image_orientation_patient", FD((1., 0., 0., 
                                                                0., 1., 0))).value
    # Use column vectors, cf. PS 3.3, C.7.6.2.1.1
    v1 = orientation[:3]
    v2 = orientation[3:]
    normal = numpy.cross(v1, v2)
    result["direction"] = numpy.transpose(numpy.asarray(
        (tuple(reversed(normal)), tuple(reversed(v2)), tuple(reversed(v1)))))
    
    return result 

def image(datasets, skipped_tags="default", do_sort=True, sort_function=None, 
          align_to_ras=False):
    """ Create an Image from the datasets. All datasets must belong to the same
        stack.
    """

    datasets = medipy.io.dicom.split.images(datasets)
    datasets = medipy.io.dicom.normalize.normalize(datasets)
    stacks = medipy.io.dicom.split.stacks(datasets)
    
    if len(stacks)>1 :
        raise medipy.base.Exception("All datasets must belong to the same stack")

    if do_sort :
        medipy.io.dicom.sort.sort(stacks[0], sort_function)
    
    dictionary = metadata(stacks[0], skipped_tags)
    origin, spacing, direction = dictionary["origin"], dictionary["spacing"], dictionary["direction"]
    del dictionary["origin"]
    del dictionary["spacing"]
    del dictionary["direction"]
    
    _image = Image(data=data(stacks[0]), 
                  origin=origin, spacing=spacing, direction=direction,
                  metadata=dictionary)
    
    if align_to_ras :
        to_axis_aligned_ras_space(_image)
    
    return _image
