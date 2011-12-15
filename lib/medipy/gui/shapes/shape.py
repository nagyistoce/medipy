##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import logging
from vtk import vtkActor

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
    
    def _set_position(self, position):
        """ Set the position of the shape by moving the actor to specified  
            position
        """
        self._position = position
        # numpy->vtk
        self._actor.SetPosition(*reversed(position))
    
    def _set_color(self, color):
        """ Set the color of the shape by setting the color of the actor. """
        self._color = color
        self._actor.GetProperty().SetColor(color)
    
    def _set_size(self, size):
        self._size = size
        
        if size <= 0 :
            self._actor.VisibilityOff()
        else :
            self._actor.VisibilityOn()
    
    def _set_filled(self, filled):
        """ Set the filled status of the shape by rendering either as a surface
            or as wireframe.
        """
        self._filled = filled
        
        if filled : 
            self._actor.GetProperty().SetRepresentationToSurface()
        else :
            self._actor.GetProperty().SetRepresentationToWireframe()
    
    # Define all "set" properties as lazy, so that derived classes can override
    # them
    actor = property(lambda self : self._actor)
    position = property(lambda self : self._position, lambda o, x : o._set_position(x))
    color = property(lambda self : self._color, lambda o, x : o._set_color(x))
    size = property(lambda self : self._size, lambda o, x : o._set_size(x))
    filled = property(lambda self : self._filled, lambda o, x : o._set_filled(x))
    
