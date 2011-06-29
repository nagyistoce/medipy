##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import numpy
from vtk import vtkTextActor

import medipy.gui.shapes

class ImageAnnotation(object) :
    """ 2D representation of a base.ImageAnnotation
    """

    # Mapping of shape ID to shape name
    _annotation_id_to_shape = dict( 
        [ (getattr(medipy.base.ImageAnnotation.Shape, name), name)
            for name in dir(medipy.base.ImageAnnotation.Shape) if name[:2] != "__"
        ]
    )
    
    def __init__(self, annotation, world_to_slice, display_coordinates, origin, spacing) :
        
        ############################
        # Property-related members #
        ############################
        
        self._annotation = None
        self._world_to_slice = None
        self._display_coordinates = None
        self._origin = None
        self._spacing = None
        self._slice_position = None
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
        self._set_world_to_slice(world_to_slice)
        self._set_display_coordinates(display_coordinates)
        self._set_origin(origin)
        self._set_spacing(spacing)
    
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
    
    def _get_world_to_slice(self) :
        """ Transformation matrix from world space to slice space.
        """
        
        return self._world_to_slice
    
    def _set_world_to_slice(self, world_to_slice) :
        
        self._world_to_slice = world_to_slice
        
        self._update_shape()
        self._update_label()
    
    def _get_display_coordinates(self) :
        """ Display annotation using physical or index coordinates.
        """
        
        return self._display_coordinates
    
    def _set_display_coordinates(self, display_coordinates) :
        
        self._display_coordinates = display_coordinates
        self._update_shape()
        self._update_label()
    
    def _get_origin(self) :
        """ Origin of the image, used when displaying in index coordinates.
        """
        
        return self._origin
    
    def _set_origin(self, origin) :
        
        self._origin = origin
        self._update_shape()
        self._update_label()
    
    def _get_spacing(self) :
        """ Spacing of the image, used when displaying in index coordinates.
        """
        
        return self._spacing
    
    def _set_spacing(self, spacing) :
        
        self._spacing = spacing
        self._update_shape()
        self._update_label()
    
    def _get_slice_position(self) :
        """ Position of the slice place.
        """
        
        return self._slice_position
    
    def _set_slice_position(self, slice_position) :
        self._slice_position = slice_position
        
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
    world_to_slice = property(_get_world_to_slice, _set_world_to_slice)
    display_coordinates = property(_get_display_coordinates, _set_display_coordinates)
    origin = property(_get_origin, _set_origin)
    spacing = property(_get_spacing, _set_spacing)
    slice_position = property(_get_slice_position, _set_slice_position)
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
    
        # TODO : display_coordinates
    
        if None in [self._annotation, self._world_to_slice, self._slice_position] :
            return
    
        altitude = self._shape.actor.GetPosition()[2]
        if self._display_coordinates == "physical" :
            position = numpy.dot(self._world_to_slice, self._annotation.position)
        else :
            index_position = numpy.subtract(self._annotation.position, self._origin)
            index_position /= self._spacing
            position = numpy.dot(self._world_to_slice, index_position)
        position[0] = altitude
        
        self._shape.position = position
        self._shape.color = self.annotation.color
        
        # Apparent size of the shape is how high the annotation is above the 
        # slice plane.
        # TODO : should depend on the annotation shape
        p1 = numpy.dot(self._world_to_slice, self._annotation.position)
        p2 = numpy.dot(self._world_to_slice, self._slice_position)
        distance = abs(p1[0]-p2[0])
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
