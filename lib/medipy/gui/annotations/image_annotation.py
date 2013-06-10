##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import numpy
from vtk import vtkTextActor

import medipy.base
import medipy.gui.shapes

class ImageAnnotation(object) :
    """ 2D representation of a base.ImageAnnotation
    """

    def __init__(self, annotation, layer) :
        
        ############################
        # Property-related members #
        ############################
        
        self._annotation = None
        self._layer = None
        self._slice_position_world = None
        self._text_actor = vtkTextActor()
        
        ###################
        # Private members #
        ###################
        self._shape = None
        self._renderer = None
        
        ##################
        # Initialization #
        ##################
        
        text_property = self._text_actor.GetTextProperty()
        text_property.ShadowOff()
        text_property.SetFontFamilyToArial()
        text_property.BoldOff()
        text_property.ItalicOff()
        text_property.SetFontSize(12)
        
        self._set_annotation(annotation)
        self._set_layer(layer)
    
    def _get_annotation(self) :
        """ Instance of base.Image
        """
        
        return self._annotation
    
    def _set_annotation(self, annotation) :
        self._annotation = annotation
        
        self._shape = self._build_shape(
            annotation.shape, annotation.position, annotation.color,
            annotation.size, annotation.filled)
        
        self._update_shape()
        self._update_label()
        
        for event in ["position", "color", "size", "filled", "label"] :
            self._annotation.add_observer(event, self._on_annotation_modified)
        self._annotation.add_observer("shape", self._on_annotation_shape_modified)
    
    def _get_layer(self) :
        
        return self._layer
    
    def _set_layer(self, layer) :
        
        if self._layer :
            self._layer.remove_observer("any", self._on_layer_modified)
        
        self._layer = layer
        self._layer.add_observer("any", self._on_layer_modified)
        
        self._update_shape()
        self._update_label()
    
    def _get_slice_position_world(self) :
        """ Position of the slice, in world VTK coordinates, VTK order.
        """
        
        return self._slice_position_world
    
    def _set_slice_position_world(self, position) :
        self._slice_position_world = position
        
        self._update_shape()
        self._update_label()
    
    def _get_shape_actor(self) :
        """ Actor containing the shape.
        """
        
        return self._shape.actor
    
    def _get_text_actor(self) :
        """ Text actor containing the label of the annotation.
        """
        
        return self._text_actor
    
    def _get_renderer(self):
        return self._renderer
    
    def _set_renderer(self, renderer):
        self._renderer = renderer
    
    annotation = property(_get_annotation, _set_annotation)
    layer = property(_get_layer, _set_layer)
    slice_position_world = property(_get_slice_position_world, _set_slice_position_world)
    shape_actor = property(_get_shape_actor)
    text_actor = property(_get_text_actor)
    renderer = property(_get_renderer, _set_renderer)
    
    #####################
    # Private interface #
    #####################
    
    def _build_shape(self, shape, *args, **kwargs) :
        """ Factory for the shape class
        """
        
        shapes = {
            medipy.base.ImageAnnotation.Shape.sphere : medipy.gui.shapes.Disk,
            medipy.base.ImageAnnotation.Shape.cube : medipy.gui.shapes.Square,
            medipy.base.ImageAnnotation.Shape.cross : medipy.gui.shapes.Cross,
        }
        return shapes[shape](*args, **kwargs)
    
    def _update_shape(self) :
        """ Update the shape to reflect the current state (annotation, 
            world-to-slice, display coordinates, origin, and spacing).
        """
    
        if None in [self._annotation, self._layer] :
            return
    
        altitude = self._shape.actor.GetPosition()[2]
        position = self._layer.physical_to_world(self._annotation.position)[::-1]
        
        shape_position = numpy.copy(position)
        shape_position[0] = altitude
        
        self._shape.position = shape_position
        self._shape.color = self.annotation.color
        
        # Apparent size of the shape is how high the annotation is above the 
        # slice plane.
        # TODO : should depend on the annotation shape
        if None not in [self._annotation.size, self.slice_position_world] :
            # position is in numpy order, slice_position_world in VTK order 
            distance = abs(position[0]-self.slice_position_world[-1])
            self._shape.size = max(0, self._annotation.size-distance)
        
        self._shape.filled = self.annotation.filled
    
    def _update_label(self) :
        """ Update the label to reflect the current state of self._shape.
        """
        
        # TODO : display_coordinates
        
        if None in [self._annotation, self._shape] :
            return
        
        if self._shape.size == 0 :
            self._text_actor.VisibilityOff()
        else :
            self._text_actor.VisibilityOn()
            self._text_actor.GetPositionCoordinate().SetCoordinateSystemToWorld()
            position = numpy.add(self._shape.position, self._shape.size)
            self._text_actor.SetPosition(position[-1], position[-2])
        
        self._text_actor.GetTextProperty().SetColor(self._annotation.color)
        if self._annotation.label :
            self._text_actor.SetVisibility(self.shape_actor.GetVisibility())
            self._text_actor.SetInput(self._annotation.label)
        else :
            self._text_actor.VisibilityOff()
    
    def _on_layer_modified(self, event) :
        """ Event handler called when the layer has changed.
        """
        
        self._update_shape()
        self._update_label()
    
    def _on_annotation_modified(self, event) :
        """ Event handler called when annotation has changed.
        """
        
        self._update_shape()
        self._update_label()

    def _on_annotation_shape_modified(self, event):
        """ Event handler called when the annotation shape is modified
        """
        
        altitude = self.shape_actor.GetPosition()[2]
        
        if self._renderer is not None :
            self._renderer.RemoveActor(self.shape_actor)
        
        self._shape = self._build_shape(
            self._annotation.shape, self._annotation.position, self._annotation.color,
            self._annotation.size, self._annotation.filled)
        
        actor_position = self.shape_actor.GetPosition()
        self.shape_actor.SetPosition(actor_position[0], actor_position[2], altitude)
        
        if self._renderer is not None :
            self._renderer.AddActor(self.shape_actor)
        
        self._update_shape()
        self._update_label()
