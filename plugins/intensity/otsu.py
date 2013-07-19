##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import itk
import medipy.itk 

def otsu_single_class(input, inside, outside):
    """ Single-class Otsu threshold
        
        <gui>
            <item name="input" type="Image" label="Input"/>
            <item name="inside" type="Float" initializer="1" label="Inside"
                tooltip="Value of the voxels in the foreground"/>
            <item name="outside" type="Float" initializer="0" label="Outside"
                tooltip="Value of the voxels in the background"/>
            <item name="output" type="Image" initializer="output=True"
                role="return" label="Output"/>
        </gui>
    """
    
    itk_input = medipy.itk.medipy_image_to_itk_image(input, False)
    # ITK considers that "inside" is between -inf and threshold, so we switch
    # inside and outside
    otsu_threshold = itk.OtsuThresholdImageFilter[itk_input, itk_input].New(
        InsideValue=outside, OutsideValue=inside, NumberOfHistogramBins=256)
    otsu_threshold(itk_input)
    
    itk_output = otsu_threshold[0]
    
    return medipy.itk.itk_image_to_medipy_image(itk_output, None, True)

def otsu_multiple_classes(input, classes):
    """ Multiple-classes Otsu thresholding
        
        <gui>
            <item name="input" type="Image" label="Input"/>
            <item name="classes" type="Int" initializer="2" label="Number of classes"/>
            <item name="output" type="Image" initializer="output=True"
                role="return" label="Output"/>
        </gui>
    """
    
    itk_input = medipy.itk.medipy_image_to_itk_image(input, False)
    otsu_threshold = itk.OtsuMultipleThresholdsImageFilter[itk_input, itk_input].New(
        NumberOfThresholds=classes, NumberOfHistogramBins=256)
    otsu_threshold(itk_input)
    
    itk_output = otsu_threshold[0]
    return medipy.itk.itk_image_to_medipy_image(itk_output, None, True)

def otsu(input, classes):
    """ Otsu thresholding
    
        <gui>
            <item name="input" type="Image" label="Input"/>
            <item name="classes" type="Int" initializer="1" label="Number of classes"/>
            <item name="output" type="Image" initializer="output=True"
                role="return" label="Output"/>
        </gui>
    """ 
    
    if classes == 1 :
        return otsu_single_class(input, 1, 0)
    elif classes > 1 :
        return otsu_multiple_classes(input, classes)
