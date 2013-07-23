##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import numpy
import scipy.ndimage

def median_filter(input, size, output):
    """ Applies a median filter of given size on input image.
    
        <gui>
            <item name="input" type="Image" label="Input image"/>
            <item name="size" type="Int" initializer="1" label="Size"/>
            <item name="output" type="Image" initializer="output=True"
                role="output" label="Output image"/>
        </gui>
    """ 
    
    if input.shape != output.shape :
       output.data = numpy.ndarray(shape=input.shape, dtype=input.dtype)
       output.copy_information(input)		
    
    output.data = scipy.ndimage.median_filter(input, size)
    
def gaussian_filter(input, sigma, output):
    """ Applies a median filter of given standard deviation on input image.
    
        <gui>
            <item name="input" type="Image" label="Input image"/>
            <item name="sigma" type="Float" initializer="1" label="St. dev."/>
            <item name="output" type="Image" initializer="output=True"
                role="output" label="Output image"/>
        </gui>
    """ 
    
    if input.shape != output.shape :
       output.data = numpy.ndarray(shape=input.shape, dtype=input.dtype)
       output.copy_information(input)
    
    output.data = scipy.ndimage.gaussian_filter(input, sigma)
