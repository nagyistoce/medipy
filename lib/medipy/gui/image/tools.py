##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

""" This module contain tools related to slice interactions.
    
"""

import numpy
from vtk import vtkCoordinate

class Tool(object) :
    """ Base class for mouse and keyboard tools.
    """
    
    def __init__(self) :
        self._coordinate = vtkCoordinate()
        self._coordinate.SetCoordinateSystemToDisplay()
    
    def _display_to_slice(self, p, slice) :
        """ Return the slice position (numpy order) in the renderer space related
            to the display position given as (x,y).
        """
        
        self._coordinate.SetValue(*p)
            
        world = self._coordinate.GetComputedWorldValue(slice.renderer)
        return list(reversed(world))
    
    def _display_to_world(self, p, slice) :
        """ Return the world position (numpy order) in the renderer space related
            to the display position given as (x,y).
        """
        
        slice_position = self._display_to_slice(p, slice)
        
        # _display_to_slice always returns the same z coordinate, adjust it.
        if slice.display_coordinates == "physical":
            position = slice.cursor_physical_position
        else :
            position = slice.cursor_index_position
        if len(position) == 2 :
            position = numpy.hstack([0, position])
        slice_position[0] = numpy.dot(slice.world_to_slice, position)[0]
        
        # Convert to renderer's world coordinates
        world_position = numpy.dot(slice.slice_to_world, slice_position)
        
        return world_position
    
    def _display_to_image_physical(self, p, slice) :
        """ Return the physical position (numpy order) in the image space
            related to the display position given as (x,y).
        """
        
        world_position = self._display_to_world(p, slice)
        
        if slice.display_coordinates == "physical" :
            # renderer's world coordinates match image physical coordinates
            return world_position
        else :
            # renderer's world coordinates match image index coordinates
            image = slice.layers[0].image
            return image.origin+world_position*image.spacing
    
    def _display_to_image_index(self, p, slice) :
        """ Return the index position (numpy order) in the image space
            related to the display position given as (x,y).
        """
        
        world_position = self._display_to_world(p, slice)
        
        if slice.display_coordinates == "physical" :
            # renderer's world coordinates match image physical coordinates
            if slice.layers :
                origin = slice.layers[0].image.origin
                if len(origin) == 2 :
                    origin = numpy.hstack([0, origin])
                spacing = slice.layers[0].image.spacing
                if len(spacing) == 2 :
                    spacing = numpy.hstack([1, spacing])
                return (world_position-origin)/spacing
            else :
                return world_position
            
        else :
            # renderer's world coordinates match image index coordinates
            return world_position

class MouseTool(Tool):
    """ Base class for all mouse tools.
    """
    
    def __init__(self) :
        super(MouseTool, self).__init__()
    
    def select(self) :
        """ Action to perform when the user selects this tool """
        pass
    
    def deselect(self) :
        """ Action to perform when the user deselects this tool """ 
        pass
    
    def start_interaction(self, rwi, slice) :
        """ Action to perform when the user clicks on the mouse button """
        pass
    
    def dispatch_interaction(self, rwi, slice) :
        """ Action to perform when the user moves the mouse with the button 
            pressed 
        """
        pass
    
    def stop_interaction(self, rwi, slice) :
        """ Action to perform when the user releases the mouse button """
        pass

class KeyboardTool(Tool):
    def __init__(self) :
        super(KeyboardTool, self).__init__()
    
    def select(self) :
        """ Action to perform when the user selects this tool """
        pass
    
    def deselect(self) :
        """ Action to perform when the user deselects this tool """ 
        pass
    
    def press(self, rwi, slice):
        """ Action to perform when the user presses a key. """
        
        pass
    
    def release(self, rwi, slice):
        """ Action to perform when the user releases a key. """
        
        pass
