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

def label_connected_components(input, output):
    """ Label connected components in an image.
        
        <gui>
            <item name="input" type="Image" label="Input"/>
            <item name="output" type="Image" initializer="output = True" 
                role="output" label="Output"/>
        </gui>
    
    """
    if input.dtype == numpy.uint16 :
        itk_input = medipy.itk.medipy_image_to_itk_image(input,False)
    else :
        input_uint16 = input.astype(numpy.uint16)
        itk_input = medipy.itk.medipy_image_to_itk_image(input_uint16, True)
    
    connected_component_filter = itk.ConnectedComponentImageFilter[itk_input, itk_input].New()
    connected_component_filter(itk_input)
    
    itk_output = connected_component_filter[0]
    
    if output.dtype == numpy.uint16 :
        medipy.itk.itk_image_to_medipy_image(itk_output, output, True)
    else :
        output.data = medipy.itk.itk_image_to_array(itk_output, True).astype(output.dtype)
        output.spacing = [x for x in reversed(itk_output.GetSpacing())]

def largest_connected_component(input, output):
    """ Get the largest connected component from a labelled image.
        
        <gui>
            <item name="input" type="Image" label="Input"/>
            <item name="output" type="Image" initializer="output = True" 
                role="output" label="Output"/>
        </gui>
    """
    
    if input.dtype == numpy.uint16 :
        itk_input = medipy.itk.medipy_image_to_itk_image(input,False)
    else :
        input_uint16 = input.astype(numpy.uint16)
        itk_input = medipy.itk.medipy_image_to_itk_image(input_uint16, True)
    
    relabel_component_filter = itk.RelabelComponentImageFilter[itk_input, itk_input].New()
    relabel_component_filter(itk_input)
    
    volumes = [relabel_component_filter.GetSizeOfObjectInPixels(1+i)
               for i in range(relabel_component_filter.GetNumberOfObjects())]
    
    largest_label = 1+volumes.index(max(volumes))
    
    threshold_filter = itk.BinaryThresholdImageFilter[itk_input, itk_input].New(
        LowerThreshold=largest_label, UpperThreshold=largest_label, 
        InsideValue=largest_label, OutsideValue=0)
    threshold_filter(relabel_component_filter.GetOutput())
    
    itk_output = threshold_filter[0]
    
    if output.dtype == numpy.uint16 :
        medipy.itk.itk_image_to_medipy_image(itk_output, output, True)
    else :
        output.data = medipy.itk.itk_image_to_array(itk_output, True).astype(output.dtype)
        output.spacing = [x for x in reversed(itk_output.GetSpacing())]