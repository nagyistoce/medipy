##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import itk
import medipy.itk

def addition(image, constant):
    """ Add constant to each input image voxel and return the result.
    
        <gui>
            <item name="image" type="Image" label="Image"/>
            <item name="constant" type="Float" label="Constant"/>
            <item name="result" type="Image" role="return"
                  initializer="output=True" label="Result"/>
        </gui>
    """
    
    return constant_operation(image, constant, itk.AddImageFilter)

def subtraction(image, constant):
    """ Subtract constant from each input image voxel and return the result.
    
        <gui>
            <item name="image" type="Image" label="Image"/>
            <item name="constant" type="Float" label="Constant"/>
            <item name="result" type="Image" role="return"
                  initializer="output=True" label="Result"/>
        </gui>
    """
    
    return constant_operation(image, constant, itk.SubtractImageFilter)

def multiplication(image, constant):
    """ Multiply by constant each input image voxel and return the result.
    
        <gui>
            <item name="image" type="Image" label="Image"/>
            <item name="constant" type="Float" label="Constant"/>
            <item name="result" type="Image" role="return"
                  initializer="output=True" label="Result"/>
        </gui>
    """
    
    return constant_operation(image, constant, itk.MultiplyImageFilter)

def division(image, constant):
    """ Divide by constant each input image voxel and return the result.
    
        <gui>
            <item name="image" type="Image" label="Image"/>
            <item name="constant" type="Float" label="Constant"/>
            <item name="result" type="Image" role="return"
                  initializer="output=True" label="Result"/>
        </gui>
    """
    
    return constant_operation(image, constant, itk.DivideImageFilter)

def constant_operation(image, constant, filter_class) :
    """ Perform a pixelwise operation using an ITK filter on the images,
        return the result
    """ 
    
    itk_image = medipy.itk.medipy_image_to_itk_image(image, False)
    filter_ = filter_class[itk_image, itk_image, itk_image].New(
        Input = itk_image, Constant = constant)
    itk_output = filter_()[0]
    
    output = medipy.itk.itk_image_to_medipy_image(itk_output, None, True)
    
    return output
