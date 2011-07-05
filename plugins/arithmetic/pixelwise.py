import itk
import medipy.itk

def addition(image1, image2):
    """ Add (pixel-wise) the two input images and return the result.
    
        <gui>
            <item name="image1" type="Image" label="Image 1"/>
            <item name="image2" type="Image" label="Image 2"/>
            <item name="result" type="Image" role="return"
                  initializer="output=True" label="Result"/>
        </gui>
    """
    
    return pixelwise_operation(image1, image2, itk.AddImageFilter)

def subtraction(image1, image2):
    """ Subtract (pixel-wise) the two input images and return the result.
    
        <gui>
            <item name="image1" type="Image" label="Image 1"/>
            <item name="image2" type="Image" label="Image 2"/>
            <item name="result" type="Image" role="return"
                  initializer="output=True" label="Result"/>
        </gui>
    """
    
    return pixelwise_operation(image1, image2, itk.SubtractImageFilter)

def multiplication(image1, image2):
    """ Multiply (pixel-wise) the two input images and return the result.
    
        <gui>
            <item name="image1" type="Image" label="Image 1"/>
            <item name="image2" type="Image" label="Image 2"/>
            <item name="result" type="Image" role="return"
                  initializer="output=True" label="Result"/>
        </gui>
    """
    
    return pixelwise_operation(image1, image2, itk.MultiplyImageFilter)

def division(image1, image2):
    """ Divide (pixel-wise) the two input images and return the result.
    
        <gui>
            <item name="image1" type="Image" label="Image 1"/>
            <item name="image2" type="Image" label="Image 2"/>
            <item name="result" type="Image" role="return"
                  initializer="output=True" label="Result"/>
        </gui>
    """
    
    return pixelwise_operation(image1, image2, itk.DivideImageFilter)

def absolute_value(image):
    """ Return the pixel-wise absolute value of image.
    
        <gui>
            <item name="image" type="Image" label="Image"/>
            <item name="result" type="Image" role="return"
                  initializer="output=True" label="Result"/>
        </gui>
    """
    
    itk_image = medipy.itk.medipy_image_to_itk_image(image, False)
    filter = itk.AbsImageFilter[itk_image, itk_image].New(
        Input = itk_image)
    itk_output = filter()[0]
    
    output = medipy.itk.itk_image_to_medipy_image(itk_output, None, True)
    
    return output

def operation(image1, image2, filter_class):
    """ Perform a pixelwise operation using an ITK filter on the images,
        return the result
    """ 
    
    itk_image_1 = medipy.itk.medipy_image_to_itk_image(image1, False)
    itk_image_2 = medipy.itk.medipy_image_to_itk_image(image2, False)
    filter = filter_class[itk_image_1, itk_image_2, itk_image_1].New(
        Input1 = itk_image_1, Input2 = itk_image_2)
    itk_output = filter()[0]
    
    output = medipy.itk.itk_image_to_medipy_image(itk_output, None, True)
    
    return output

#import itk
#import numpy
#
#import medipy.itk
#


#
#def absolute_value(input, output):
#    """ Compute the absolute value element_wise
#        
#        :gui:
#            input : Image
#                Input
#            output : Image : output=True
#                Output
#    """
#    if input.shape != output.shape :
#       output.data = numpy.ndarray(shape=input.shape, dtype=input.dtype)
#    output[:] = numpy.abs(input)
#    output.copy_information(input)
#    

