import itk
import medipy.itk

def crop(input, index, shape):
    """ Return a sub-image from the input image, starting at index and with 
        given shape.
        
        <gui>
            <item name="input" type="Image" label="Input"/>
            <item name="index" type="Coordinates" label="Start"/>
            <item name="shape" type="Array" initializer="type=int" label="Shape"/>
            <item name="output" type="Image" initializer="output=True" role="return"
                label="Output"/>
        </gui>
    """
    
    itk_input = medipy.itk.medipy_image_to_itk_image(input, False)
    
    itk_index = [x for x in reversed(index)]
    itk_shape = [x for x in reversed(shape)]
    
    # Use RegionOfInterestImageFilter since we wish to modify the output's origin
    region = itk_input.GetRequestedRegion().__class__(itk_index, itk_shape)
    filter = itk.RegionOfInterestImageFilter[itk_input, itk_input].New(
        Input=itk_input, RegionOfInterest=region)
    filter()
    itk_output = filter[0]
    
    output = medipy.itk.itk_image_to_medipy_image(itk_output, None, True)
    
    return output
