##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import numpy

import medipy.base.array
from tools import KeyboardTool

class MoveCursor(KeyboardTool):
    """ Move the cursor on the slice using the arrows keys and Page Up/Page Down.
    """
    
    def __init__(self):
        super(MoveCursor, self).__init__()
    
    def press(self, rwi, slice):
        
        image = slice.layers[0].image
        
        world_to_slice = medipy.base.array.reshape(slice.world_to_slice, 
            (image.ndim, image.ndim), "constant", False, value=0)
        # Add ones on the diagonal when necessary
        for rank in range(image.ndim) :
            if numpy.less_equal(slice.world_to_slice.shape, rank).all() : 
                world_to_slice[image.ndim-rank-1, image.ndim-rank-1] = 1.
        
        if slice.display_coordinates in ["physical", "nearest_axis_aligned"] :
            matrix = numpy.dot(world_to_slice, image.direction)
        else :
            matrix = world_to_slice
        
        if rwi.GetKeySym() == "Left" :
            motion = -matrix[-1]
        elif rwi.GetKeySym() == "Right" :
            motion = matrix[-1]
        elif rwi.GetKeySym() == "Up" :
            motion = matrix[-2]
        elif rwi.GetKeySym() == "Down" :
            motion = -matrix[-2]
        elif matrix.shape[0] >= 3 and rwi.GetKeySym() in ["Prior", "PageUp"] :
            motion = matrix[-3]
        elif matrix.shape[0] >= 3 and rwi.GetKeySym() in ["Next", "PageDown"] :
            motion = -matrix[-3]
        else :
            motion = [0]*image.ndim
    
        image_position = numpy.add(slice.cursor_index_position, motion)
        
        if image.is_inside(image_position) :
            slice.cursor_index_position = image_position
        
        rwi.Render()

class Zoom(KeyboardTool):
    """ Zoom in (or out, depending on factor)
    """
    
    def __init__(self, factor):
        
        super(Zoom, self).__init__()
        self.factor = factor
    
    def press(self, rwi, slice):
        
        slice.zoom *= self.factor
        rwi.Render()

class ToggleInterpolation(KeyboardTool):
    """ Toggle the interpolation of the displayed image.
    """
    
    def __init__(self):
        super(ToggleInterpolation, self).__init__()
    
    def press(self, rwi, slice):
        slice.interpolation = not slice.interpolation
        rwi.Render()

class ToggleScalarBarVisibility(KeyboardTool):
    """ Toggle the display of the scalar bar.
    """
    
    def __init__(self):
        super(ToggleScalarBarVisibility, self).__init__()
    
    def press(self, rwi, slice):
        slice.scalar_bar_visibility = not slice.scalar_bar_visibility
        rwi.Render()

class ToggleOrientationVisibility(KeyboardTool):
    """ Toggle the display of the orientations.
    """
    
    def __init__(self):
        super(ToggleOrientationVisibility, self).__init__()
    
    def press(self, rwi, slice):
        slice.orientation_visibility = not slice.orientation_visibility
        rwi.Render()

class ToggleCornerAnnotationsVisibility(KeyboardTool):
    """ Toggle the display of the corner annotations.
    """
    
    def __init__(self):
        super(ToggleCornerAnnotationsVisibility, self).__init__()
    
    def press(self, rwi, slice):
        slice.corner_annotations_visibility = not slice.corner_annotations_visibility
        rwi.Render()