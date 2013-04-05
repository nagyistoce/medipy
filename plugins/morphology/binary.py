##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import itk
import medipy.itk

from structuring_element import name_to_structuring_element

def erode(input, erode_value, *args, **kwargs):
    """ Binary erosion of an image using a name of a structuring element and a
        radius, or a structuring element.
        
        <gui>
            <item name="input" type="Image" label="Input"/>
            <item name="erode_value" type="Float" initializer="1" label="Erode value"/>
            <item name="shape" type="Enum" initializer="('ball', 'box','cross')"
                label="Shape"/>
            <item name="radius" type="Int" initializer="1" label="Radius"/>
            <item name="output" type="Image" initializer="output=True" role="return"
                label="Output"/>
        </gui>
    """
    
    if len(args) == 1 :
        return erode_se(input, erode_value, *args)
    elif len(args) == 2 :
        return erode_shape_and_radius(input, erode_value, *args)
    elif len(args) == 0 :
        if "shape" in kwargs :
            return erode_shape_and_radius(input, erode_value, **kwargs)
        elif "structuring_element" in kwargs :
            return erode_se(input, erode_value, **kwargs)
        else :
            raise Exception("Incorrect parameters")
    else :
        raise Exception("Incorrect parameters")

def erode_shape_and_radius(input, erode_value, shape, radius):
    """ Binary erosion of an image using a name of a structuring element and a
        radius
    """
    
    structuring_element = name_to_structuring_element(shape, input.ndim, radius)
    return erode_se(input, erode_value, structuring_element)

def erode_se(input, erode_value, structuring_element):
    """ Binary erosion of an image using a structuring element.
    """
    
    itk_input = medipy.itk.medipy_image_to_itk_image(input, False)
    
    filter = itk.BinaryErodeImageFilter[itk_input, itk_input, structuring_element].New(
        Input=itk_input, ErodeValue=erode_value, BackgroundValue=0,
        Kernel=structuring_element)
    
    itk_output = filter()[0]
    output = medipy.itk.itk_image_to_medipy_image(itk_output, None, True)
    return output

def dilate(input, dilate_value, *args, **kwargs):
    """ Binary dilation of an image using a name of a structuring element and a
        radius, or a structuring element.
        
        <gui>
            <item name="input" type="Image" label="Input"/>
            <item name="dilate_value" type="Float" initializer="1" label="Dilate value"/>
            <item name="shape" type="Enum" initializer="('ball', 'box','cross')"
                label="Shape"/>
            <item name="radius" type="Int" initializer="1" label="Radius"/>
            <item name="output" type="Image" initializer="output=True" role="return"
                label="Output"/>
        </gui>
    """
    
    if len(args) == 1 :
        return dilate_se(input, dilate_value, *args)
    elif len(args) == 2 :
        return dilate_shape_and_radius(input, dilate_value, *args)
    elif len(args) == 0 :
        if "shape" in kwargs :
            return dilate_shape_and_radius(input, dilate_value, **kwargs)
        elif "structuring_element" in kwargs :
            return dilate_se(input, dilate_value, **kwargs)
        else :
            raise Exception("Incorrect parameters")
    else :
        raise Exception("Incorrect parameters")

def dilate_shape_and_radius(input, dilate_value, shape, radius):
    """ Binary dilation of an image using a name of a structuring element and a
        radius
    """
    structuring_element = name_to_structuring_element(shape, input.ndim, radius)
    return dilate_se(input, dilate_value, structuring_element)

def dilate_se(input, dilate_value, structuring_element):
    """ Binary dilation of an image using a structuring element.
    """
    
    itk_input = medipy.itk.medipy_image_to_itk_image(input, False)
    
    filter = itk.BinaryDilateImageFilter[itk_input, itk_input, structuring_element].New(
        Input=itk_input, DilateValue=dilate_value, BackgroundValue=0, 
        Kernel=structuring_element)
    
    itk_output = filter()[0]
    output = medipy.itk.itk_image_to_medipy_image(itk_output, None, True)
    return output