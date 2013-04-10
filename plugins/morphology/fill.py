##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import numpy
import scipy.ndimage
import medipy.base

def fill_2d_cavities(input, axis):
    """ Fill 2D cavities on each slice of given image axis, using 6-connectivity.
        
        <gui>
            <item name="input" type="Image" label="Input"/>
            <item name="axis" type="Int" initializer="range = (0, len(${input}.shape)-1)"
                label="Axis"/>
            <item name="output" type="Image" initializer="output = True"
                role="return" label="Output"/>
        </gui>
    """
    
    shape = [3,3]
    shape.insert(axis, 1)
    structure = numpy.ones(shape, dtype=int)
    
    corners = [ [0, 0], [0, 2], [2, 0], [2, 2] ]
    
    for corner in corners :
        corner.insert(axis, 0)
        structure[tuple(corner)] = 0
    
    data = scipy.ndimage.binary_fill_holes(input.data, structure).astype(input.dtype)
    output = medipy.base.Image(data=data)
    output.copy_information(input)
    return output

def fill_3d_cavities(input):
    """ Fill 3D cavities using 6-connectivity.
        
        <gui>
            <item name="input" type="Image" label="Input"/>
            <item name="output" type="Image" initializer="output = True"
                role="return" label="Output"/>
        </gui>
    """
    
    shape = [3,3,3]
    structure = numpy.ones(shape, dtype=int)
    
    corners = [ (0, 0, 0), (0, 0, 2), (0, 2, 0), (0, 2, 2),
                (2, 0, 0), (2, 0, 2), (2, 2, 0), (2, 2, 2) ]
    
    for corner in corners :
        structure[corner] = 0
    
    data = scipy.ndimage.binary_fill_holes(input.data, structure).astype(input.dtype)
    output = medipy.base.Image(data=data)
    output.copy_information(input)
    return output