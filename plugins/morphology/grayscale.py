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

def erode(input, *args, **kwargs):
    """ Gray-scale erosion of an image using a name of a structuring element and a
        radius, or a structuring element.
        
        <gui>
            <item name="input" type="Image" label="Input"/>
            <item name="shape" type="Enum" initializer="('ball', 'box','cross')"
                label="Shape"/>
            <item name="radius" type="Int" initializer="1" label="Radius"/>
            <item name="output" type="Image" initializer="output=True" role="return"
                label="Output"/>
        </gui>
    """
    
    if len(args) == 1 :
        return erode_se(input, *args)
    elif len(args) == 2 :
        return erode_shape_and_radius(input, **args)
    elif len(args) == 0 :
        if "shape" in kwargs :
            return erode_shape_and_radius(input, **kwargs)
        elif "structuring_element" in kwargs :
            return erode_se(input, **kwargs)
        else :
            raise Exception("Incorrect parameters")
    else :
        raise Exception("Incorrect parameters")

def erode_shape_and_radius(input, shape, radius):
    """ Gray-scale erosion of an image using a name of a structuring element and a
        radius
    """
    
    structuring_element = name_to_structuring_element(shape, input.ndim, radius)
    return erode_se(input, structuring_element)

def erode_se(input, structuring_element):
    """ Gray-scale erosion of an image using a structuring element.
    """
    
    itk_input = medipy.itk.medipy_image_to_itk_image(input, False)
    
    filter = itk.GrayscaleErodeImageFilter[itk_input, itk_input, structuring_element].New(
        Input=itk_input, Kernel=structuring_element)
    
    itk_output = filter()[0]
    output = medipy.itk.itk_image_to_medipy_image(itk_output, None, True)
    return output

def dilate(input, *args, **kwargs):
    """ Gray-scale dilation of an image using a name of a structuring element and a
        radius, or a structuring element.
        
        <gui>
            <item name="input" type="Image" label="Input"/>
            <item name="shape" type="Enum" initializer="('ball', 'box','cross')"
                label="Shape"/>
            <item name="radius" type="Int" initializer="1" label="Radius"/>
            <item name="output" type="Image" initializer="output=True" role="return"
                label="Output"/>
        </gui>
    """
    
    if len(args) == 1 :
        return dilate_se(input, *args)
    elif len(args) == 2 :
        return dilate_shape_and_radius(input, **args)
    elif len(args) == 0 :
        if "shape" in kwargs :
            return dilate_shape_and_radius(input, **kwargs)
        elif "structuring_element" in kwargs :
            return dilate_se(input, **kwargs)
        else :
            raise Exception("Incorrect parameters")
    else :
        raise Exception("Incorrect parameters")

def dilate_shape_and_radius(input, shape, radius):
    """ Gray-scale dilation of an image using a name of a structuring element and a
        radius
    """
    structuring_element = name_to_structuring_element(shape, input.ndim, radius)
    return dilate_se(input, structuring_element)

def dilate_se(input, structuring_element):
    """ Gray-scale dilation of an image using a structuring element.
    """
    
    itk_input = medipy.itk.medipy_image_to_itk_image(input, False)
    
    filter = itk.GrayscaleDilateImageFilter[itk_input, itk_input, structuring_element].New(
        Input=itk_input, Kernel=structuring_element)
    
    itk_output = filter()[0]
    output = medipy.itk.itk_image_to_medipy_image(itk_output, None, True)
    return output

def open(input, *args, **kwargs):
    """ Gray-scale opening of an image using a name of a structuring element and a
        radius, or a structuring element.
        
        <gui>
            <item name="input" type="Image" label="Input"/>
            <item name="shape" type="Enum" initializer="('ball', 'box','cross')"
                label="Shape"/>
            <item name="radius" type="Int" initializer="1" label="Radius"/>
            <item name="output" type="Image" initializer="output=True" role="return"
                label="Output"/>
        </gui>
    """
    
    if len(args) == 1 :
        return open_se(input, *args)
    elif len(args) == 2 :
        return open_shape_and_radius(input, **args)
    elif len(args) == 0 :
        if "shape" in kwargs :
            return open_shape_and_radius(input, **kwargs)
        elif "structuring_element" in kwargs :
            return open_se(input, **kwargs)
        else :
            raise Exception("Incorrect parameters")
    else :
        raise Exception("Incorrect parameters")

def open_shape_and_radius(input, shape, radius):
    """ Gray-scale opening of an image using a name of a structuring element and a
        radius
    """
    structuring_element = name_to_structuring_element(shape, input.ndim, radius)
    return open_se(input, structuring_element)

def open_se(input, structuring_element):
    """ Gray-scale opening of an image using a structuring element.
    """
    
    itk_input = medipy.itk.medipy_image_to_itk_image(input, False)
    
    erode_filter = itk.GrayscaleErodeImageFilter[itk_input, itk_input, structuring_element].New(
        Input=itk_input, Kernel=structuring_element)
    dilate_filter = itk.GrayscaleDilateImageFilter[itk_input, itk_input, structuring_element].New(
        Input=erode_filter[0], Kernel=structuring_element)
    
    itk_output = dilate_filter()[0]
    output = medipy.itk.itk_image_to_medipy_image(itk_output, None, True)
    return output

def close(input, *args, **kwargs):
    """ Gray-scale closing of an image using a name of a structuring element and a
        radius, or a structuring element.
        
        <gui>
            <item name="input" type="Image" label="Input"/>
            <item name="shape" type="Enum" initializer="('ball', 'box','cross')"
                label="Shape"/>
            <item name="radius" type="Int" initializer="1" label="Radius"/>
            <item name="output" type="Image" initializer="output=True" role="return"
                label="Output"/>
        </gui>
    """
    
    if len(args) == 1 :
        return close_se(input, *args)
    elif len(args) == 2 :
        return close_shape_and_radius(input, **args)
    elif len(args) == 0 :
        if "shape" in kwargs :
            return close_shape_and_radius(input, **kwargs)
        elif "structuring_element" in kwargs :
            return close_se(input, **kwargs)
        else :
            raise Exception("Incorrect parameters")
    else :
        raise Exception("Incorrect parameters")

def close_shape_and_radius(input, shape, radius):
    """ Gray-scale closing of an image using a name of a structuring element and a
        radius
    """
    structuring_element = name_to_structuring_element(shape, input.ndim, radius)
    return close_se(input, structuring_element)

def close_se(input, structuring_element):
    """ Gray-scale closing of an image using a structuring element.
    """
    
    itk_input = medipy.itk.medipy_image_to_itk_image(input, False)
    
    dilate_filter = itk.GrayscaleDilateImageFilter[itk_input, itk_input, structuring_element].New(
        Input=itk_input, Kernel=structuring_element)
    erode_filter = itk.GrayscaleErodeImageFilter[itk_input, itk_input, structuring_element].New(
        Input=dilate_filter[0], Kernel=structuring_element)
    
    itk_output = erode_filter()[0]
    output = medipy.itk.itk_image_to_medipy_image(itk_output, None, True)
    return output