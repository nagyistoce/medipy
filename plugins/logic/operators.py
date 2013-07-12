##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import itk
import numpy
import medipy.itk

def pixelwise_and(image1, image2):
    """ Return the pixelwise value of the AND operator
    
        <gui>
            <item name="image1" type="Image" label="Image 1"/>
            <item name="image2" type="Image" label="Image 2"/>
            <item name="result" type="Image" role="return"
                  initializer="output=True" label="Result"/>
        </gui>
    """
    
    return pixelwise_operation(image1, image2, itk.AndImageFilter)

def pixelwise_or(image1, image2):
    """ Return the pixelwise value of the OR operator
    
        <gui>
            <item name="image1" type="Image" label="Image 1"/>
            <item name="image2" type="Image" label="Image 2"/>
            <item name="result" type="Image" role="return"
                  initializer="output=True" label="Result"/>
        </gui>
    """
    
    return pixelwise_operation(image1, image2, itk.OrImageFilter)

def pixelwise_xor(image1, image2):
    """ Return the pixelwise value of the XOR operator
    
        <gui>
            <item name="image1" type="Image" label="Image 1"/>
            <item name="image2" type="Image" label="Image 2"/>
            <item name="result" type="Image" role="return"
                  initializer="output=True" label="Result"/>
        </gui>
    """
    
    return pixelwise_operation(image1, image2, itk.XorImageFilter)

def pixelwise_not(image):
    """ Return the pixel-wise value of the NOT operator.
    
        <gui>
            <item name="image" type="Image" label="Image"/>
            <item name="result" type="Image" role="return"
                  initializer="output=True" label="Result"/>
        </gui>
    """
    
    itk_image = medipy.itk.medipy_image_to_itk_image(image, False)
    filter = itk.NotImageFilter[itk_image, itk_image].New(
        Input = itk_image)
    itk_output = filter()[0]
    
    output = medipy.itk.itk_image_to_medipy_image(itk_output, None, True)
    
    return output

def pixelwise_operation(image1, image2, filter_class):
    """ Perform a pixelwise operation using an ITK filter on the images,
        return the result
    """ 
    
    if image1.dtype == numpy.uint8 :
        itk_image1 = medipy.itk.medipy_image_to_itk_image(image1, False)
    else :
        image1_uint8 = image1.astype(numpy.uint8)
        itk_image1 = medipy.itk.medipy_image_to_itk_image(image1_uint8, True)
    
    if image2.dtype == numpy.uint8 :
        itk_image2= medipy.itk.medipy_image_to_itk_image(image2, False)
    else :
        image2_uint8 = image2.astype(numpy.uint8)
        itk_image2 = medipy.itk.medipy_image_to_itk_image(image2_uint8, True)
    
    filter = filter_class[itk_image1, itk_image2, itk_image1].New(
        Input1 = itk_image1, Input2 = itk_image2)
    itk_output = filter()[0]
    
    output = medipy.itk.itk_image_to_medipy_image(itk_output, None, True)
    if image1.dtype != numpy.uint8 :
        output.data = output.data.astype(image1.dtype)
    
    return output
