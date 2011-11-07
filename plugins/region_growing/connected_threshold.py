import itk
import medipy.itk
import region_growing

def connected_threshold(input, lower, upper, value, seeds):
    """ Grow a region from a seed (or set of seeds)
    """
    
    itk_input = medipy.itk.medipy_image_to_itk_image(input, False)
    
    filter = itk.ConnectedThresholdImageFilter[itk_input, itk_input].New(
        Lower=lower, Upper=upper, ReplaceValue=value,
        Input=itk_input)
    for seed in seeds :
        filter.AddSeed(list(reversed(seed)))
    filter()
    
    itk_output = filter[0]
    output = medipy.itk.itk_image_to_medipy_image(itk_output, None, False)
    
    return output

def connected_threshold_with_radius(input, lower, upper, radius, value, seeds):
    """ Grow a region from a seed (or set of seeds)
    """
    
    itk_input = medipy.itk.medipy_image_to_itk_image(input, False)
    
    filter = itk.ConnectedThresholdWithRadiusImageFilter[itk_input, itk_input].New(
        Lower=lower, Upper=upper, 
        Radius=radius, ReplaceValue=value, Input=itk_input)
    for seed in seeds :
        l = [int(x) for x in reversed(seed)]
        filter.AddSeed(l)
    filter()
    
    itk_output = filter[0]
    return medipy.itk.itk_image_to_medipy_image(itk_output, None, True)
    