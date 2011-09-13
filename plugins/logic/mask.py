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

def create_mask(input, background):
    """ Create a binary image by setting every non-background voxel to 1,
        and every background voxel to 0.
        
        <gui>
            <item name="input" type="Image" label="Input"/>
            <item name="background" type="Float" initializer="0" label="Background"/>
            <item name="output" type="Image" initializer="output = True"
                role="return" label="Output"/>
        </gui>
    """
    
    output = medipy.base.Image(input.shape, dtype=input.dtype, value=1)
    output.copy_information(input)
    numpy.not_equal(input, background, output.data)
    
    return output
