##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import itk
import numpy

import medipy.itk
import medipy.logic

def create_mask(input, background):
    """ Create a binary image by setting every non-background voxel to 1,
        and every background voxel to 0.
        
        <gui>
            <item name="input" type="Image" label="Input"/>
            <item name="background" type="Float" initializer="0" label="Background"/>
            <item name="output" type="Image" initializer="output = True"
                role="return" label="Output"/>
        </gui>
    """
    
    output = medipy.base.Image(input.shape, dtype=input.dtype, value=1)
    output.copy_information(input)
    numpy.not_equal(input, background, output.data)
    
    return output

def apply_mask(input, mask, background, outside):
    """ Apply the mask to the input image : the value of the output image is :
          * input(p) if input(p) != background
          * outside otherwise
    
        <gui>
            <item name="input" type="Image" label="Input"/>
            <item name="mask" type="Image" label="Mask"/>
            <item name="background" type="Float" initializer="0" label="Background"/>
            <item name="outside" type="Float" initializer="0" label="Outside"/>
            <item name="output" type="Image" initializer="output = True"
                role="return" label="Output"/>
        </gui>
    """
    
    itk_input = medipy.itk.medipy_image_to_itk_image(input, False)
    itk_mask = medipy.itk.medipy_image_to_itk_image(mask, False)
    mask_filter = itk.MaskWithValueImageFilter[itk_input, itk_mask, itk_input].New(
        itk_input, Mask=itk_mask, BackgroundValue=background, OutsideValue=outside)
    itk_output = mask_filter()[0]
    return medipy.itk.itk_image_to_medipy_image(itk_output, None, True)
