##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import sys

import numpy
import scipy.spatial
from vtk import (vtkActor, vtkCornerAnnotation, vtkLineSource, 
    vtkPolyDataMapper, vtkRenderer, vtkScalarBarActor)

import medipy.base
import medipy.base.array
from medipy.base import ObservableList, PropertySynchronized
from medipy.gui import colormaps
from medipy.gui.colormap import Colormap
from medipy.gui.annotations import ImageAnnotation as GUIImageAnnotation
from contour_layer import ContourLayer
from image_layer import ImageLayer
from medipy.vtk import vtkOrientationAnnotation

import mouse_tools
import keyboard_tools 

class Slice(PropertySynchronized) :
    """ Synchronized representation of several layers, with keyboard and mouse
        interaction, and with optional annotations.
        
        Default keyboard interactions are :
          * left, right, up, down : move cursor one voxel
          * page up, page down : move slice position one voxel
          * "+" : zoom in
          * "-" : zoom out
          * "i" : toggle interpolation
          * "s" : toggle scalar bar visibility
          * "c" : toggle corner annotations visibility
          * "o" : toggle orientation annotations visibility
        
        Default mouse interactions are :
          * Left : cursor position
          * Middle : center under mouse when clicked, pan when dragged
          * Right : window and level
          * Mouse : zoom
    """
    
    # Altitudes of the different actors, layers will be place between 0 and 1
    _actors_altitudes = {
        "outline" : -0.02,
        "back_plane" : -0.1,
        "annotations" : 1.01,
        "cursor" : 2.1,
        "camera" : 10,
    }
    
    # Index of location for vtkOrientationAnnotation
    _orientation_annotation_index = {
        "up" : 1, "down" : 3, "left" : 2, "right" : 0
    }
    
    # Index of location for vtkCornerAnnotation
    _corner_annotation_index = {
        "up_left" : 2, "up_right" : 3, "down_left" : 0, "down_right" : 1
    }
    
    def __init__(self, world_to_slice, layers=None, annotations=None,
                 interpolation=False, display_coordinates="physical", 
                 scalar_bar_visibility = False, orientation_visibility=True,
                 corner_annotations_visibility=False) :
        
        layers = layers or []
        annotations = annotations or ObservableList()
        
        ############################
        # Property-related members #
        ############################
        self._interpolation = None
        self._display_coordinates = None
        self._scalar_bar_visibility = None
        self._orientation_visibility = None
        self._corner_annotations_visibility = None
        
        self._world_to_slice = None
        
        self._slice_to_world = None
        
        self._layers = []
        
        self._annotations = None
        self._gui_annotations = {}
        
        self._image_physical_position = None
        self._image_index_position = None
        
        self._cursor_physical_position = None
        self._cursor_index_position = None
        
        self._zoom = None
        
        self._mouse_tools = {}
        self._keyboard_tools = {}
        
        self._renderer = vtkRenderer()
        
        ###################
        # Private members #
        ###################
        
        # World-to-slice matrix, with rows and columns added or removed so that
        # it is 3x3.
        self._3d_world_to_slice = None
        self._3d_slice_to_world = None
        
        # Slice extent is the physical extent of all layers, 
        # given as (x_min, x_max, y_min, y_max)
        self._slice_extent = (-100, 100, -100, 100)
        
        # VTK objects
        self._scalar_bar_actor = vtkScalarBarActor()
        self._corner_annotation = vtkCornerAnnotation()
        self._orientation_annotation = vtkOrientationAnnotation()
        self._horizontal_line_source = vtkLineSource()
        self._vertical_line_source = vtkLineSource()
        self._horizontal_line_actor = vtkActor()
        self._vertical_line_actor = vtkActor()
        
        # Tools and interactions
        self._observer_tags = []
        self._active_source = None
        
        ##################
        # Initialization #
        ##################
        
        super(Slice, self).__init__([
            "world_to_slice", "interpolation", "display_coordinates", 
            "scalar_bar_visibility", "orientation_visibility", 
            "corner_annotations_visibility", "zoom", 
        ])
        self.add_allowed_event("cursor_position")
        self.add_allowed_event("image_position")
        self.add_allowed_event("center")
        
        # Configure camera
        camera = self._renderer.GetActiveCamera()
        camera.ParallelProjectionOn()
        camera.SetPosition(0, 0, self._actors_altitudes["camera"])
        camera.SetFocalPoint(0, 0, 0)
        
        # Create cursor objects
        for direction in ["horizontal", "vertical"] :
            line_source = getattr(self, "_%s_line_source"%direction)
            actor = getattr(self, "_%s_line_actor"%direction)
            
            mapper = vtkPolyDataMapper()
            mapper.SetInputConnection(line_source.GetOutputPort())
            actor.SetMapper(mapper)
            actor.PickableOff()
            actor.GetProperty().SetColor(1,0,0)
            self._renderer.AddActor(actor)
        
        # Create scalar bar (from vtkInria3D)
        self._scalar_bar_actor.GetLabelTextProperty().SetColor(1.0,1.0,1.0)
        self._scalar_bar_actor.GetTitleTextProperty().SetColor(1.0,1.0,1.0)
        self._scalar_bar_actor.GetLabelTextProperty().BoldOff()
        self._scalar_bar_actor.GetLabelTextProperty().ShadowOff()
        self._scalar_bar_actor.GetLabelTextProperty().ItalicOff()
        self._scalar_bar_actor.SetNumberOfLabels(3)
        self._scalar_bar_actor.GetLabelTextProperty().SetFontSize(8)
        self._scalar_bar_actor.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
        self._scalar_bar_actor.SetWidth(0.1)
        self._scalar_bar_actor.SetHeight(0.5)
        self._scalar_bar_actor.SetPosition(0.9,0.3)
        self._scalar_bar_actor.PickableOff()
        self._renderer.AddActor(self._scalar_bar_actor)
        
        # Setup text-annotation actors
        self._corner_annotation.SetNonlinearFontScaleFactor(0.3)
        self._renderer.AddActor(self._corner_annotation)
        self._orientation_annotation.SetNonlinearFontScaleFactor(0.25)
        self._renderer.AddActor(self._orientation_annotation)
        
        self._set_interpolation(interpolation)
        self._set_display_coordinates(display_coordinates)
        
        self._set_scalar_bar_visibility(scalar_bar_visibility)
        self._set_orientation_visibility(orientation_visibility)
        self._set_corner_annotations_visibility(corner_annotations_visibility)
        
        self._set_world_to_slice(world_to_slice)
        
        for layer in layers :
            self.append_layer(**layer)
        
        if annotations is not None :
            self._set_annotations(annotations)
        
        # Position slice at middle of layer 0
        self.reset_view()
    
        # Configure default tools
        self.set_mouse_button_tool("Left", mouse_tools.Select())
        self.set_mouse_button_tool("Middle", mouse_tools.Pan())
        self.set_mouse_button_tool("Right", mouse_tools.WindowLevel())
        self.set_wheel_tool("Forward", mouse_tools.Zoom(1.1))
        self.set_wheel_tool("Backward", mouse_tools.Zoom(1./1.1))
        self.set_keyboard_tool("Left", keyboard_tools.MoveCursor())
        self.set_keyboard_tool("Right", keyboard_tools.MoveCursor())
        self.set_keyboard_tool("Up", keyboard_tools.MoveCursor())
        self.set_keyboard_tool("Down", keyboard_tools.MoveCursor())
        self.set_keyboard_tool("Prior", keyboard_tools.MoveCursor())
        self.set_keyboard_tool("Next", keyboard_tools.MoveCursor())
        self.set_keyboard_tool("PageUp", keyboard_tools.MoveCursor())
        self.set_keyboard_tool("PageDown", keyboard_tools.MoveCursor())
        self.set_keyboard_tool("+", keyboard_tools.Zoom(1.1))
        self.set_keyboard_tool("-", keyboard_tools.Zoom(1./1.1))
        self.set_keyboard_tool("i", keyboard_tools.ToggleInterpolation())
        self.set_keyboard_tool("s", keyboard_tools.ToggleScalarBarVisibility())
        self.set_keyboard_tool("c", keyboard_tools.ToggleCornerAnnotationsVisibility())
        self.set_keyboard_tool("o", keyboard_tools.ToggleOrientationVisibility())
    
    def close(self):
        """ Remove all actors from renderer, prepare for destruction.
        """
    
        self.renderer.RemoveActor(self._horizontal_line_actor)
        self.renderer.RemoveActor(self._vertical_line_actor)
        self.renderer.RemoveActor(self._scalar_bar_actor)
        self.renderer.RemoveActor(self._orientation_annotation)
        self.renderer.RemoveActor(self._corner_annotation)
        
        for layer in self._layers :
            self.renderer.RemoveActor(layer.actor)
        
        for gui_annotation in self._gui_annotations.values() :
            self.renderer.RemoveActor(gui_annotation.shape_actor)
            self.renderer.RemoveActor(gui_annotation.text_actor)
    
    def append_layer(self, *args, **kwargs) :
        """ Append a new layer.
        """
        
        self.insert_layer(len(self._layers), *args, **kwargs)
    
    def insert_layer(self, index, image, colormap=None, opacity=1.0) :
        """ Insert a new layer at specified position. The colormap defaults to
            a gray colormap. If the colormap's display range is None, it 
            defaults to the image range.
        """
        
        if colormap is None :
            colormap = Colormap(colormaps["gray"], None, False, False, False)
        
        if colormap.display_range is None :
            colormap.display_range = (image.data.min(), image.data.max())
        
        # Find out which layer class we will use
        classes = {
            "spectroscopy" : ContourLayer
        }
        LayerClass = classes.get(image.image_type, ImageLayer)
        
        # Create the Layer and insert it in the list
        layer = LayerClass(
            self._world_to_slice, image, self._display_coordinates, colormap, opacity)
        self._layers.insert(index, layer)
        
        # Update the physical extent
        self._compute_extent()
        
        # The scalar bar will always reflect layer 0
        if index == 0 :
            self._scalar_bar_actor.SetLookupTable(layer.colormap.vtk_colormap)
        
        # Adjust layer w.r.t. the current state.
        self._update_layers_positions()
        if self._cursor_physical_position is not None : 
            layer.physical_position = self._cursor_physical_position
        if isinstance(layer, ImageLayer) :
            layer.actor.SetInterpolate(self._interpolation)
        
        # And finally add it to the renderer
        self._renderer.AddActor(layer.actor)
    
    def delete_layer(self, index) :
        """ Remove a layer from the list.
        """
        
        # Remove the actor, delete the list item, and update the other layers.
        self._renderer.RemoveActor(self._layers[index].actor)
        del self._layers[index]
        self._update_layers_positions()
    
    def get_layer_visibility(self, index):
        return self._layers[index].actor.GetVisibility()
    
    def set_layer_visibility(self, index, visibility):
        return self._layers[index].actor.SetVisibility(visibility)
    
    def center_on_physical_position(self, position):
        self._set_cursor_physical_position(position)
        self._set_image_physical_position(position)
        self.notify_observers("center")
    
    def center_on_index_position(self, position):
        self._set_cursor_index_position(position)
        self._set_image_index_position(position)
        self.notify_observers("center")
    
    def reset_view(self):
        if self._layers :
            image = self._layers[0].image
            center = numpy.divide(image.shape, 2.).round()
            self._set_cursor_index_position(center)
            self._set_image_index_position(center)
        else :
            self._set_cursor_physical_position((0,0,0))
            self._set_image_physical_position((0,0,0))
        
        self._set_zoom(1.0)
    
    def setup_orientation_annotation(self) :
        """ Update orientation annotation to reflect the image-to-slice
            projection
        """
        
        # Anatomical directions in LPS convention, numpy order
        directions_anatomical = {
            "L" : (0,0,+1),
            "R" : (0,0,-1),
            "P" : (0,+1,0),
            "A" : (0,-1,0),
            "I" : (-1,0,0),
            "S" : (+1,0,0),
        }
        
        # Index directions, numpy order
        directions_index = {
            "+x" : (0,0,+1),
            "-x" : (0,0,-1),
            "+y" : (0,+1,0),
            "-y" : (0,-1,0),
            "+z" : (-1,0,0),
            "-z" : (+1,0,0),
        }
        
        directions = (directions_anatomical 
                      if self.display_coordinates in ["physical", "nearest_axis_aligned"]
                      else directions_index)
        
        # Window locations
        locations = {
            "up" : (1,0),
            "down" : (-1,0),
            "left" : (0,-1),
            "right" : (0,1)
        }
        
        for location, p in locations.items() :
            matrix = self._3d_world_to_slice
            direction = numpy.dot(self._3d_slice_to_world, numpy.hstack((0, p)))
            
            # Find closest in-slice direction based on dot product
            closest = None
            max_distance = -1
            for name, d in directions.items() :
                distance = numpy.dot(d, direction)
                if distance > max_distance :
                    max_distance = distance
                    closest = name
            
            # Set text
            index = self._orientation_annotation_index[location]
            self._orientation_annotation.SetText(index, closest)
    
    def get_label(self, where) :
        if where in self._orientation_annotation_index.keys() :
            index = self._orientation_annotation_index[where]
            return self._orientation_annotation.GetText(index)
        else :
            index = self._corner_annotation_index[where]
            return self._corner_annotation.GetText(index)
    
    def set_label(self, where, label) :
        if where in self._orientation_annotation_index.keys() :
            index = self._orientation_annotation_index[where]
            self._orientation_annotation.SetText(index, label)
        else :
            index = self._corner_annotation_index[where]
            return self._corner_annotation.SetText(index, label)
    
    def setup_rwi(self, rwi) :
        rwi.SetInteractorStyle(None)
        self._observer_tags.append(rwi.AddObserver("LeftButtonPressEvent", self._start_interaction))
        self._observer_tags.append(rwi.AddObserver("MiddleButtonPressEvent", self._start_interaction))
        self._observer_tags.append(rwi.AddObserver("RightButtonPressEvent", self._start_interaction))
        
        self._observer_tags.append(rwi.AddObserver("LeftButtonReleaseEvent", self._stop_interaction))
        self._observer_tags.append(rwi.AddObserver("MiddleButtonReleaseEvent", self._stop_interaction))
        self._observer_tags.append(rwi.AddObserver("RightButtonReleaseEvent", self._stop_interaction))
        
        self._observer_tags.append(rwi.AddObserver("MouseWheelForwardEvent", self._dispatch_interaction))
        self._observer_tags.append(rwi.AddObserver("MouseWheelBackwardEvent", self._dispatch_interaction))
        
        self._observer_tags.append(rwi.AddObserver("MouseMoveEvent", self._dispatch_interaction))
        
        self._observer_tags.append(rwi.AddObserver("KeyPressEvent", self._key_press))
        self._observer_tags.append(rwi.AddObserver("KeyReleaseEvent", self._key_release))
        
        self._observer_tags.append(rwi.AddObserver("ConfigureEvent", self._window_resize))
        
        self._set_parallel_scale()
    
    def unset_rwi(self, rwi):
        for tag in self._observer_tags :
            rwi.RemoveObserver(tag)
    
    def get_mouse_button_tool(self, button):
        event_name = "%sButton"%button
        return self._mouse_tools.get(event_name, None)
    
    def set_mouse_button_tool(self, button, tool) :
        """ Set a tool associated with given button (Left, Middle, or Right),
            with an optional modifier (Shift or Control). Example : Right,
            ShiftLeft. Set tool to None to have no tool connected to the button.
        """
        
        event_name = "%sButton"%button
        
        if event_name in self._mouse_tools : 
            self._mouse_tools[event_name].deselect()
        
        if tool :
            tool.select()
            self._mouse_tools[event_name] = tool
        elif event_name in self._mouse_tools :
            del self._mouse_tools[event_name]
    
    def get_wheel_tool(self, direction) :
        """ Return the tool associated with a mouse wheel direction (Forward or 
            Backward)
        """
        
        event_name = "MouseWheel%s"%direction
        return self._mouse_tools.get(event_name, None)
    
    def set_wheel_tool(self, direction, tool) :
        """ Set a tool associated with a mouse wheel direction (Forward or 
            Backward)
        """
        
        event_name = "MouseWheel%s"%direction
        
        if event_name in self._mouse_tools : 
            self._mouse_tools[event_name].deselect()
        
        if tool :
            tool.select()
            self._mouse_tools[event_name] = tool
        elif event_name in self._mouse_tools :
            del self._mouse_tools[event_name]
    
    def get_keyboard_tool(self, key):
        return self._keyboard_tools[key]
    
    def set_keyboard_tool(self, key, tool):
        if key in self._keyboard_tools : 
            self._keyboard_tools[key].deselect()
        if tool :
            tool.select()
            self._keyboard_tools[key] = tool
        elif key in self._keyboard_tools :
            del self._keyboard_tools[key]
        
    ##############
    # Properties #
    ##############
    
    def _get_annotations(self):
        return self._annotations
    
    def _set_annotations(self, annotations):
    
        if self._annotations is not None :
            self._annotations.remove_observer("any", self._on_annotations_changed)
            for slice_annotation in self._gui_annotations.values() :
                self._renderer.RemoveActor(slice_annotation.shape_actor)
                self._renderer.RemoveActor(slice_annotation.text_actor)
        
        self._gui_annotations = {}
                
        self._annotations = annotations
        self._annotations.add_observer("any", self._on_annotations_changed)
        
        for annotation in annotations :
            gui_annotation = GUIImageAnnotation(annotation, self._layers[0])
            gui_annotation.slice_position = self._cursor_physical_position
            gui_annotation.renderer = self._renderer
            
            actor_position = gui_annotation.shape_actor.GetPosition()
            gui_annotation.shape_actor.SetPosition(
                actor_position[0], actor_position[1], 
                self._actors_altitudes["annotations"]
            )
            self._gui_annotations[annotation] = gui_annotation
            
            self.renderer.AddActor(gui_annotation.shape_actor)
            self.renderer.AddActor(gui_annotation.text_actor)
    
    def _get_interpolation(self) :
        """ Interpolate the displayed data or not.
        """
        
        return self._interpolation
    
    def _set_interpolation(self, interpolation) :
        self._interpolation = interpolation
        for layer in self._layers :
            if isinstance(layer, ImageLayer) :
                layer.actor.SetInterpolate(interpolation)
        
        self.notify_observers("interpolation")
    
    def _get_display_coordinates(self) :
        """ Display image using physical or index coordinates.
        """
        
        return self._display_coordinates
    
    def _set_display_coordinates(self, display_coordinates) :
        if display_coordinates not in ["physical", "nearest_axis_aligned", "index"] :
            raise medipy.base.Exception("Unknown display coordinates : %s"%(display_coordinates,))
        
        self._display_coordinates = display_coordinates
        
        for layer in self._layers :
            layer.display_coordinates = display_coordinates
        self._compute_extent()
        for annotation in self._gui_annotations.values() :
            annotation.display_coordinates = display_coordinates
        
        # Keep the same pixel under the cursor and centered in the view
        self._locked = True
        if self._cursor_index_position is not None :
            self._set_cursor_index_position(self._get_cursor_index_position())
        if self._image_index_position is not None :
            self._set_image_index_position(self._get_image_index_position())
        self._locked = False
    
    def _get_scalar_bar_visibility(self) :
        """ Visibility of the scalar bar.
        """
        
        return self._scalar_bar_visibility
    
    def _set_scalar_bar_visibility(self, scalar_bar_visibility) :
        self._scalar_bar_visibility = scalar_bar_visibility
        self._scalar_bar_actor.SetVisibility(scalar_bar_visibility)
        self.notify_observers("scalar_bar_visibility")
    
    def _get_orientation_visibility(self) :
        """ Visibility of the anatomical orientation informations.
        """
        
        return self._orientation_visibility
    
    def _set_orientation_visibility(self, orientation_visibility) :
        self._orientation_visibility = orientation_visibility
        self._orientation_annotation.SetVisibility(orientation_visibility)
        self.notify_observers("orientation_visibility")
    
    def _get_corner_annotations_visibility(self) :
        """ Visibility of the corner annotations.
        """
        
        return self._corner_annotations_visibility
    
    def _set_corner_annotations_visibility(self, corner_annotations_visibility) :
        self._corner_annotations_visibility = corner_annotations_visibility
        self._corner_annotation.SetVisibility(corner_annotations_visibility)
        self.notify_observers("corner_annotations_visibility")
    
    def _get_world_to_slice(self) :
        """ Transformation matrix for this slice.
        """
        
        return self._world_to_slice
    
    def _set_world_to_slice(self, world_to_slice) :
        self._world_to_slice = world_to_slice
        self._slice_to_world = numpy.linalg.inv(world_to_slice)
        
        # Get a 3D matrix
        self._3d_world_to_slice = medipy.base.array.reshape(world_to_slice,
            (3,3), "constant", False, value=0)
        # Add ones on the diagonal when necessary
        for rank in range(3) :
            if numpy.less_equal(max(world_to_slice.shape), rank).all() : 
                self._3d_world_to_slice[3-rank-1, 3-rank-1] = 1.
        self._3d_slice_to_world = numpy.linalg.inv(self._3d_world_to_slice)
        
        for layer in self._layers :
            layer.world_to_slice = world_to_slice

        self.setup_orientation_annotation()        
        self._compute_extent()
        
        # Keep the same pixel under the cursor and centered in the view
        self._locked = True
        if self._cursor_index_position is not None :
            self._set_cursor_index_position(self._get_cursor_index_position())
        if self._image_index_position is not None :
            self._set_image_index_position(self._get_image_index_position())
        self._locked = False
        
        self.notify_observers("world_to_slice")
    
    def _get_slice_to_world(self) :
        """ Inverse of the transformation matrix for this slice.
        """
        
        return self._slice_to_world
    
    def _get_layers(self) :
        """ List the slice Layers. For item insertion or deletion, the
            append_layer, insert_layer, and delete_layer functions must be used.
        """
        
        return self._layers
    
    def _get_image_physical_position(self) :
        """ Physical position of the voxel to be centered in the view. 
        """
        
        return self._image_physical_position
    
    def _set_image_physical_position(self, position) :
        
        self._image_physical_position = numpy.asarray(position)
        
        if self._layers :
            image = self._layers[0].image
            
            self._image_index_position = image.physical_to_index(position)
            
            position = medipy.base.array.reshape(position, (image.ndim,), 
                "constant", False, value=0)
            world_vtk = self._layers[0].physical_to_world(position)
            
        else :
            self._image_index_position = numpy.asarray(position)
            
            position = medipy.base.array.reshape(numpy.asarray(position), 
                (3,), "constant", False, value=0)
            world_vtk = numpy.dot(self._3d_world_to_slice, position)
        
        self._renderer.GetActiveCamera().SetPosition(
            world_vtk[0], world_vtk[1], self._actors_altitudes["camera"])
        self._renderer.GetActiveCamera().SetFocalPoint(
            world_vtk[0], world_vtk[1], 0)
        
        self.notify_observers("image_position")
    
    def _get_image_index_position(self) :
        """ Index position of the voxel to be centered in the view. At least one
            layer must be present for this to be defined.
        """
        
        return self._image_index_position
    
    def _set_image_index_position(self, position) :
        
        image = self._layers[0].image
        
        # Normalize the dimension of position w.r.t. layers[0]'s image
        position = medipy.base.array.reshape(
            position, (image.ndim,), "constant", False, value=0)
        
        physical_position = image.index_to_physical(position)
        self._set_image_physical_position(physical_position)
    
    def _get_cursor_physical_position(self) :
        """ Physical position of the cursor, rounded to the center of the
            nearest pixel of layer 0. If no layer is present, the position is 
            not rounded. 
        """
        return self._cursor_physical_position
    
    def _set_cursor_physical_position(self, position) :
        
        self._cursor_physical_position = numpy.asarray(position)
        
        if self._layers :
            image = self._layers[0].image
            
            self._cursor_index_position = image.physical_to_index(position)
            
            position = medipy.base.array.reshape(position, (image.ndim,), 
                "constant", False, value=0)
            world_vtk = self._layers[0].physical_to_world(position)
            
        else :
            self._cursor_index_position = numpy.asarray(position)
            
            position = medipy.base.array.reshape(numpy.asarray(position), 
                (3,), "constant", False, value=0)
            world_vtk = numpy.dot(self._3d_world_to_slice, position)
        
        extent = self._slice_extent
        
        z = self._actors_altitudes["cursor"]
        self._horizontal_line_source.SetPoint1(extent[0], world_vtk[1], z)
        self._horizontal_line_source.SetPoint2(extent[1], world_vtk[1], z)
        self._vertical_line_source.SetPoint1(world_vtk[0], extent[2], z)
        self._vertical_line_source.SetPoint2(world_vtk[0], extent[3], z)
        
        # Update layers
        for layer in self._layers :
            layer.physical_position = self._cursor_physical_position
        
        # Update annotations
        for gui_annotation in self._gui_annotations.values() :
            gui_annotation.slice_position = self._cursor_physical_position
        
        self.notify_observers("cursor_position")
        
    def _get_cursor_index_position(self) :
        """ Index position of the cursor. At least one
            layer must be present for this to be defined.
        """
        return self._cursor_index_position
    
    def _set_cursor_index_position(self, position) :
        if self._layers :
            image = self._layers[0].image
            index_position = medipy.base.array.reshape(numpy.asarray(position),
                (image.ndim,), "constant", False, value=0)
            physical_position = image.index_to_physical(index_position)
        else :
            index_position = position
        
        self._set_cursor_physical_position(physical_position)
    
    def _get_zoom(self) :
        """ Relative zoom value. A zoom of 1 fits the data in the viewport.
        """
        
        # TODO : make it absolute zoom value : a zoom of 1 displays one data
        # pixel in one viewport pixel.
        
        return self._zoom
    
    def _set_zoom(self, zoom) :
        self._zoom = zoom
        self._set_parallel_scale()
        self.notify_observers("zoom")
    
    def _get_renderer(self) :
        """ VTK renderer associated with the slice.
        """
        
        return self._renderer
    
    annotations = property(_get_annotations, _set_annotations)
    interpolation = property(_get_interpolation, _set_interpolation)
    display_coordinates = property(_get_display_coordinates, 
                                   _set_display_coordinates)
    scalar_bar_visibility = property(_get_scalar_bar_visibility, 
                                     _set_scalar_bar_visibility)
    orientation_visibility = property(_get_orientation_visibility,
                                      _set_orientation_visibility)
    corner_annotations_visibility = property(_get_corner_annotations_visibility,
                                             _set_corner_annotations_visibility)
    world_to_slice = property(_get_world_to_slice, _set_world_to_slice)
    slice_to_world = property(_get_slice_to_world)
    layers = property(_get_layers)
    image_physical_position = property(_get_image_physical_position, 
                                       _set_image_physical_position)
    image_index_position = property(_get_image_index_position, 
                                    _set_image_index_position)
    cursor_physical_position = property(_get_cursor_physical_position,
                                        _set_cursor_physical_position)
    cursor_index_position = property(_get_cursor_index_position,
                                     _set_cursor_index_position)
    zoom = property(_get_zoom, _set_zoom)
    renderer = property(_get_renderer)

    def _update_layers_positions(self) :
        for i, layer in enumerate(self._layers) :
            altitude = float(i)/float(len(self._layers))
            
            position = list(layer.actor.GetPosition())
            position[2] = altitude
            layer.actor.SetPosition(position)
    
    def _set_parallel_scale(self):
        parallel_scale = self._renderer.GetSize()[1]/2./self._zoom
        self._renderer.GetActiveCamera().SetParallelScale(parallel_scale)
    
    def _compute_extent(self) :
        if len(self._layers) == 0 :
            self._slice_extent = (-100, 100, -100, 100)
            return
        
        slice_extent = list(self._layers[0].actor.GetBounds())
        for layer in self.layers[1:] :
            bounds = layer.actor.GetBounds()
            for i, value in enumerate(bounds) :
                if i%2 == 0 :
                    slice_extent[i] = min(slice_extent[i], value)
                else :
                    slice_extent[i] = max(slice_extent[i], value)
        
        self._slice_extent = slice_extent
    
    def _start_interaction(self, rwi, event) :
        
        if not self._renderer.IsInViewport(*rwi.GetEventPosition()) :
            return
        
        if event.endswith("PressEvent") :
            source = event[:-len("PressEvent")]

        modifier = ""
        if rwi.GetControlKey() :
            modifier += "Control"
        if rwi.GetShiftKey() :
            modifier += "Shift"
        
        if modifier+source in self._mouse_tools :
            if not self._active_source :
                self._active_source = (source, modifier)
            source_name = self._active_source[1]+self._active_source[0]
            self._mouse_tools[source_name].start_interaction(rwi, self)
    
    def _stop_interaction(self, rwi, event) :
        if not self._renderer.IsInViewport(*rwi.GetEventPosition()) :
            return
        
        if event.endswith("ReleaseEvent") :
            source = event[:-len("ReleaseEvent")]
        
        if self._active_source is not None and source == self._active_source[0] :
            source_name = self._active_source[1]+self._active_source[0]
            self._mouse_tools[source_name].stop_interaction(rwi, self)
            self._active_source = None
    
    def _dispatch_interaction(self, rwi, event) :
        if not self._renderer.IsInViewport(*rwi.GetEventPosition()) :
            return
        
        if self._active_source is not None :
            source_name = self._active_source[1]+self._active_source[0]
            self._mouse_tools[source_name].dispatch_interaction(rwi, self)
        elif event.startswith("MouseWheel") :
            source = event[:-len("Event")]
            if source in self._mouse_tools :
                self._mouse_tools[source].dispatch_interaction(rwi, self)
    
    def _get_key_name(self, rwi):
        if rwi.GetKeyCode() == "\x00" : 
            key = rwi.GetKeySym()
        else :
            key = rwi.GetKeyCode()	
        if rwi.GetControlKey() :
            key = "Ctrl"+key
        if getattr(rwi, "GetAltKey") and rwi.GetAltKey() :
            key = "Alt"+key
        
        return key
    
    def _key_press(self, rwi, dummy) :
        if not self._renderer.IsInViewport(*rwi.GetEventPosition()) :
            return
        
        key = self._get_key_name(rwi)
        
        if key in self._keyboard_tools :
            self._keyboard_tools[key].press(rwi, self)
            
    def _key_release(self, rwi, dummy) :
        if not self._renderer.IsInViewport(*rwi.GetEventPosition()) :
            return
        
        key = self._get_key_name(rwi)
        
        if key in self._keyboard_tools :
            self._keyboard_tools[key].release(rwi, self)
    
    def _window_resize(self, rwi, dummy):
        self._set_parallel_scale()
        
    def _on_annotations_changed(self, event):
        
        # Remove deleted annotations
        for annotation in list(self._gui_annotations.keys()) :
            if annotation not in self._annotations :
                gui_annotation = self._gui_annotations[annotation]
                self._renderer.RemoveActor(gui_annotation.shape_actor)
                self._renderer.RemoveActor(gui_annotation.text_actor)
                del self._gui_annotations[annotation]
            
        # Add new annotations
        for annotation in self._annotations :
            if annotation not in self._gui_annotations :
                
                gui_annotation = GUIImageAnnotation(annotation, 
                    self._3d_world_to_slice, self._display_coordinates, 
                    self._layers[0].image.origin, self._layers[0].image.spacing)
                gui_annotation.renderer = self._renderer
                
                gui_annotation.slice_position = self._cursor_physical_position
                
                actor_position = gui_annotation.shape_actor.GetPosition()
                gui_annotation.shape_actor.SetPosition(
                    actor_position[0], actor_position[1], 
                    self._actors_altitudes["annotations"]
                )
                self._gui_annotations[annotation] = gui_annotation
                
                self.renderer.AddActor(gui_annotation.shape_actor)
                self.renderer.AddActor(gui_annotation.text_actor)
