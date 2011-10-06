##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import numpy
from tools import KeyboardTool

class MoveCursor(KeyboardTool):
    """ Move the cursor on the slice using the arrows keys and Page Up/Page Down.
    """
    
    motion = {
        "Left" : numpy.asarray((0,0,-1)),
        "Right" : numpy.asarray((0,0,1)),
        "Up" : numpy.asarray((0,1,0)),
        "Down" : numpy.asarray((0,-1,0)),
        "Prior" : numpy.asarray((1,0,0)),
        "PageUp" : numpy.asarray((1,0,0)),
        "Next" : numpy.asarray((-1,0,0)),
        "PageDown" : numpy.asarray((-1,0,0)),
    }
    
    def __init__(self):
        super(MoveCursor, self).__init__()
    
    def press(self, rwi, slice):
        
        image = slice.layers[0].image
        ndim = image.ndim
        
        motion = self.motion[rwi.GetKeySym()][-ndim:]
    
        slice_position = numpy.dot(slice.world_to_slice, 
                                   slice.cursor_index_position)
        slice_position = numpy.add(slice_position, motion)
        
        image_position = numpy.dot(slice.slice_to_world, slice_position)
        
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