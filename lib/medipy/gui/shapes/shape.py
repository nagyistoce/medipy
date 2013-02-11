##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

from vtk import vtkActor

import medipy.base

class Shape(object):
    """ Base class for a simple geometric shape.
    
        This class (and its derived classes) provide a simple interface to
        render shapes of given color, position, filled status, ...
        
        Position is in mm, expressed as z,y,x
    """
    
    def __init__(self, position=None, color=None, size=None, filled=None):
        """ Create a shape.  
        
            The  following keywords can be passed to the constructor :
              * position
              * color
              * size
              * filled
        """
        
        self._actor = vtkActor()
        self._position = None
        self._color = None
        self._size = None
        self._filled = None
        
        if position is not None : 
            self._set_position(position)
        if color is not None : 
            self._set_color(color)
        if size is not None : 
            self._set_size(size)
        if filled is not None : 
            self._set_filled(filled)
    
    ##############
    # Properties #
    ##############
    
    def _get_actor(self):
        return self._actor
    
    def _get_position(self):
        """ World position of the shape, in numpy order.
        """
        
        return self._position
    
    def _set_position(self, position):
        self._position = position
        # numpy->vtk
        self._actor.SetPosition(*reversed(position))
    
    def _get_color(self):
        """ RGB color of the shape.
        """
        
        return self._color
    
    def _set_color(self, color):
        self._color = color
        self._actor.GetProperty().SetColor(color)
    
    def _get_size(self):
        """ World size of the shape.
        """
        
        return self._size
    
    def _set_size(self, size):
        self._size = size
        
        if size <= 0 :
            self._actor.VisibilityOff()
        else :
            self._actor.VisibilityOn()
    
    def _get_filled(self):
        """ Filled status of the shape.
        """
        
        return self._filled
    
    def _set_filled(self, filled):
        self._filled = filled
        
        if filled : 
            self._actor.GetProperty().SetRepresentationToSurface()
        else :
            self._actor.GetProperty().SetRepresentationToWireframe()
    
    # Define all "set" properties as lazy, so that derived classes can override
    # them
    actor = medipy.base.LateBindingProperty(_get_actor)
    position = medipy.base.LateBindingProperty(_get_position, _set_position)
    color = medipy.base.LateBindingProperty(_get_color, _set_color)
    size = medipy.base.LateBindingProperty(_get_size, _set_size)
    filled = medipy.base.LateBindingProperty(_get_filled, _set_filled)
    
