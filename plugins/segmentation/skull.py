##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import numpy

from medipy.denoising import median_filter
from medipy.intensity import otsu, binary_threshold
from medipy.morphology import fill_2d_cavities
from medipy.segmentation.label import label_connected_components, largest_connected_component

def skull(input):
    """ Segment the skull, using an algorithm from Medimax
        
        <gui>
            <item name="input" type="Image" label="Input"/>
            <item name="output" type="Image" initializer="output=True" 
                role="return" label="Output"/>
        </gui>
    """
    
    # Three-classes Otsu thresholding
    output = otsu(input, 3)
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
    # Keep the top 18 cm of the image
    nonzero=numpy.nonzero(output)
    z_end = nonzero[-3][-1]
    z_begin = max(0, z_end-180/input.spacing[-3])
    
    if z_begin != 0 :
        # Zero everything below
        output[...,:z_begin,:,:]=0 
    
    output.copy_information(input)
    return output
