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
    
    def _display_to_image_physical(self, p, slice) :
        """ Return the physical position (numpy order) in the image space
            related to the display position given as (x,y).
        """
        
        self._coordinate.SetValue(*p)
        world_position = self._coordinate.GetComputedWorldValue(slice.renderer)
        return slice.layers[0].world_to_physical(world_position)
    
    def _display_to_image_index(self, p, slice) :
        """ Return the index position (numpy order, floating point coordinates) 
            in the image space related to the display position given as (x,y).
        """
        
        self._coordinate.SetValue(*p)
        world_position = self._coordinate.GetComputedWorldValue(slice.renderer)
        return slice.layers[0].world_to_index(world_position)
    
    def _display_to_image_int_index(self, p, slice) :
        """ Return the integer index position (numpy order) in the image space
            related to the display position given as (x,y).
        """
        
        index_position = self._display_to_image_index(p, slice)
        # Coordinates are rounded since the pixels are "considered to be the
        # rectangular region surrounding the pixel center holding the data 
        # value" (ITK Software Guide, 4.1.4, p. 40)
        index_position = numpy.round(index_position).astype(int)
        
        return index_position

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
