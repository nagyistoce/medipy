##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import numpy

from medipy.denoising import median_filter
from medipy.intensity import otsu, binary_threshold
from medipy.morphology import fill_2d_cavities, grayscale_open
from medipy.segmentation.label import label_connected_components, largest_connected_component

def skull(input):
    """ Segment the skull.
        
        <gui>
            <item name="input" type="Image" label="Input"/>
            <item name="output" type="Image" initializer="output=True" 
                role="return" label="Output"/>
        </gui>
    """
    
    # Three-classes Otsu thresholding
    output = otsu(input, 3)
    # Remove small isolated components (usually background noise)
    output = grayscale_open(output, "ball", 1)
    # Binarize image : everything that is non-0 is set to 1
    output = binary_threshold(output, 0, 0, 0, 1)
    # Label image, using a cross-shaped neighborhood
    label_connected_components(output, output)
    # Keep the label with the largest volume (in voxels)
    largest_connected_component(output, output)
    # Fill cavities
    output = fill_2d_cavities(output, 0)
    output = fill_2d_cavities(output, 1)
    output = fill_2d_cavities(output, 2)
    # Apply a 3x3x3 median filter
    median_filter(output, 1, output)
    # Fill cavities
    output = fill_2d_cavities(output, 0)
    output = fill_2d_cavities(output, 1)
    output = fill_2d_cavities(output, 2)
    # Binarize again
    output = binary_threshold(output, 0, 0, 0, 1)
    # Label
    label_connected_components(output, output)
    # Keep the label with the largest volume (in voxels)
    largest_connected_component(output, output)
    
    ###################################
    # Keep the top 18 cm of the image #
    ###################################
    
    # Inferior-superior axis in voxel space and corresponding coordinate index
    is_axis = input.direction[0]
    is_coordinate_index = numpy.argmax(numpy.abs(is_axis))
    
    nonzero=numpy.nonzero(output)
    size = 180/input.spacing[is_coordinate_index] # 18 cm in voxels on the IS axis
    
    if is_axis[is_coordinate_index]>0:
        end = numpy.max(nonzero[is_coordinate_index])
        begin = max(0, int(numpy.round(end-size)))
    else:
        begin = numpy.min(nonzero[is_coordinate_index])
        end = min(output.shape[is_coordinate_index], int(numpy.round(begin+size)))
    
    s = [slice(0, x) for x in output.shape]
    
    # Zero the "lower" part of the image
    s[is_coordinate_index] = slice(0,begin)
    output[s] = 0
    
    # Zero the "upper" part of the image
    s[is_coordinate_index] = slice(end, output.shape[is_coordinate_index])
    output[s] = 0
    
    output.copy_information(input)
    return output
