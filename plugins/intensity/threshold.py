##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import itk
import medipy.itk

def binary_threshold(input, min, max, inside, outside) :
    """ Return a binarized version of the input image.
    
        <gui> 
            <item name="input" type="Image" label="Image"/>
            <item name="min" type="Float" label="Minimum value"
                tooltip="Minimum value in threshold (included)"/>
            <item name="max" type="Float" label="Maximum value"
                tooltip="Maximum value in threshold (included)"/>
            <item name="inside" type="Float" initializer="1" label="Inside value"
                tooltip="Value taken by pixels inside the image"/>
            <item name="outside" type="Float" initializer="0" label="Outside value"
                tooltip="Value taken by pixels outside the image"/>
            <item name="output" type="Image" initializer="output=True"
                role="return" label="Output image"/>
        </gui>
    """
    itk_input = medipy.itk.medipy_image_to_itk_image(input, False)
    
    threshold_filter = itk.BinaryThresholdImageFilter[itk_input, itk_input].New(
        LowerThreshold=min, UpperThreshold=max, 
        InsideValue=inside, OutsideValue=outside)

    threshold_filter(itk_input)
    itk_output = threshold_filter[0]
    
    return medipy.itk.itk_image_to_medipy_image(itk_output, None, True)

def threshold(input, min, max, outside) :
    """ Return a thresholded version of the input image.
    
        <gui> 
            <item name="input" type="Image" label="Image"/>
            <item name="min" type="Float" label="Minimum value"
                tooltip="Minimum value in threshold (included)"/>
            <item name="max" type="Float" label="Maximum value"
                tooltip="Maximum value in threshold (included)"/>
            <item name="outside" type="Float" initializer="0" label="Outside value"
                tooltip="Value taken by pixels outside the image"/>
            <item name="output" type="Image" initializer="output=True"
                role="return" label="Output image"/>
        </gui>
    """
    itk_input = medipy.itk.medipy_image_to_itk_image(input, False)
    
    threshold_filter = itk.ThresholdImageFilter[itk_input].New(
        Lower=min, Upper=max, OutsideValue=outside)

    threshold_filter(itk_input)
    itk_output = threshold_filter[0]
    
    return medipy.itk.itk_image_to_medipy_image(itk_output, None, True)
