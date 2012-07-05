##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import logging
import numpy
from vtk import vtkCornerAnnotation, vtkRenderer
import wx

import medipy.base
from medipy.base import ObservableList, PropertySynchronized
from medipy.gui import colormaps
from medipy.gui.colormap import Colormap
from medipy.gui.wxVTKRenderWindowInteractor import wxVTKRenderWindowInteractor

from slice import Slice

def get_informations(image):
    """ Return the informations to be displayed on image.
    """
    
    if image.number_of_layers == 0 :
        return {}
    
    informations = {}
    
    # Cursor position, voxels and mm
    if image.cursor_index_position is not None :
        index_position = tuple(reversed(image.cursor_index_position))
        physical_position = tuple(reversed(image.cursor_physical_position))
        
        if index_position is not None :
            informations["up_left"] = ("Position (voxels) : "+
                                       (" x ".join(len(index_position)*["%i"]))%index_position)
        if physical_position is not None :
            informations["up_left"] += "\n"
            informations["up_left"] += ("Position (mm) : "+
                                        (" x ".join(len(physical_position)*["%i"]))%physical_position)
        
        # Voxel value
        informations["up_left"] += "\n"
        if image.get_layer_image(0).is_inside(image.cursor_index_position) : 
            value = image.get_layer_image(0)[tuple(image.cursor_index_position)]
        else :
            value = ""
        informations["up_left"] += str(value)
        
    # Image size and spacing
    shape = tuple(reversed(image.get_layer_image(0).shape))
    spacing = tuple(reversed(image.get_layer_image(0).spacing))
    informations["down_left"] = ("Image size : "+
                                 (" x ".join(len(shape)*["%i"]))%shape)
    informations["down_left"] += "\n"
    informations["down_left"] += ("Voxel size (mm): "+
                                  (" x ".join(len(spacing)*["%.2f"]))%spacing)
    
    # Zoom
    if image.zoom is not None :
        informations["down_left"] += "\n"
        informations["down_left"] += "Zoom : %.2f %%"%(100.*image.zoom)
    
    return informations

class Image(wx.Panel, PropertySynchronized):
    
    _viewport = {
        "axial" : (0.0, 0.0, 0.5, 0.5),
        "coronal" : (0.0, 0.5, 0.5, 1.0),
        "sagittal" : (0.5, 0.5, 1.0, 1.0),
        "informations" : (0.5, 0.0, 1.0, 0.5)
    }
    
    # Index of location for vtkCornerAnnotation
    _corner_annotation_index = {
        "up_left" : 2, "up_right" : 3, "down_left" : 0, "down_right" : 1
    }
    
    def __init__(self, parent, slice_mode="multiplanar", layers=None,
                 annotations=None, interpolation=False,
                 display_coordinates="physical", scalar_bar_visibility = False,
                 orientation_visibility=True, corner_annotations_visibility=False,
                 convention="radiological", *args, **kwargs):
        
        if annotations is None :
            annotations = ObservableList()
        
        ##############
        # Properties #
        ##############
        
        self._slice_mode = None
        
        self._interpolation = None
        self._display_coordinates = None
        self._scalar_bar_visibility = None
        self._orientation_visibility = None
        self._convention = None
        
        self._annotations = None
        
        self._cursor_physical_position = None
        self._cursor_index_position = None
        
        self._center_physical_position = None
        self._center_index_position = None
        
        self._zoom = None
        self._mouse_tools = {}
        self._keyboard_tools = {}
        
        ###################
        # Private members #
        ###################
        self._rwi = None
        self._layers = ObservableList()
        self._slices = []
        self._slices_names = []
        self._informations_renderer = vtkRenderer()
        self._informations_renderer.SetViewport(*self._viewport["informations"])
        
        self._informations_corner_annotation = vtkCornerAnnotation()
        self._informations_corner_annotation.SetNonlinearFontScaleFactor(0.3)
        self._informations_renderer.AddActor(self._informations_corner_annotation)
        
        ##################
        # Initialization #
        ##################
        
        wx.Panel.__init__(self, parent, *args, **kwargs)
        self._rwi = wxVTKRenderWindowInteractor(self, wx.ID_ANY)
        self._rwi.Enable(1)
        sizer = wx.BoxSizer()
        sizer.Add(self._rwi, 1, wx.EXPAND|wx.ALL, 3)
        self.SetSizer(sizer)
        self._rwi.Bind(wx.EVT_LEFT_DOWN, self._on_button_down)
        self._rwi.Bind(wx.EVT_MIDDLE_DOWN, self._on_button_down)
        self._rwi.Bind(wx.EVT_RIGHT_DOWN, self._on_button_down)
        self._rwi.Bind(wx.EVT_MOUSEWHEEL, self._on_button_down)
        self.Bind(wx.EVT_CLOSE, self._on_close)
        
        PropertySynchronized.__init__(self, ["interpolation", 
            "display_coordinates", "scalar_bar_visibility", 
            "orientation_visibility", "zoom"])
        self.add_allowed_event("cursor_position")
        self.add_allowed_event("image_position")
        self.add_allowed_event("center")
        for event in ["data", "display_range", "cut_low", "cut_high", "zero_transparency"] :
            self.add_allowed_event("colormap_{0}".format(event))
#        self.add_allowed_event("layer_modified")
#        self.add_allowed_event("activated")

        self._set_interpolation(interpolation)
        self._set_display_coordinates(display_coordinates)
        self._set_scalar_bar_visibility(scalar_bar_visibility)
        self._set_orientation_visibility(orientation_visibility)
        self._set_convention(convention)
        for layer in layers or [] :
            self.append_layer(**layer)
        self._set_slice_mode(slice_mode)
        
        self._set_annotations(annotations)
        
        # Position image at middle of layer 0
        self.reset_view()
    
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
        
        self._layers.insert(index, {"image" : image, "colormap" : colormap, 
                                    "opacity" : opacity})
        
        for slice in self._slices :
            slice.insert_layer(index, image, colormap, opacity)
        
        if len(self._slices)>1 :
            for event in ["data", "display_range", "cut_low", "cut_high", "zero_transparency"] :
                self._slices[0].layers[index].colormap.add_observer(event, self._on_colormap_event) 
            for slice_index, slice in enumerate(self._slices) :
                layer = slice.layers[index]
                next_slice_layer = self._slices[(slice_index+1)%len(self._slices)].layers[index]
                layer.colormap.append_child(next_slice_layer.colormap)
    
    def delete_layer(self, index) :
        """ Remove a layer from the list.
        """
        
        for event in ["data", "display_range", "cut_low", "cut_high", "zero_transparency"] :
            self._slices[0].layers[index].colormap.remove_observer(event, self._on_colormap_event)
        
        del self._layers[index]
        for slice in self._slices :
            slice.delete_layer(index)
    
    def get_layer_visibility(self, index):
        return self._slices[0].get_layer_visibility(index)
    
    def set_layer_visibility(self, index, visibility):
        for slice in self._slices :
            slice.set_layer_visibility(index, visibility)
    
    def get_layer_image(self, index):
        return self._layers[index]["image"]
    
    def set_layer_image(self, index, image):
        self._layers[index]["image"] = image
        for slice in self._slices :
            slice.layers[index].image = image
    
    def get_layer_colormap(self, index):
        return self._layers[index]["colormap"]
    
    def set_layer_colormap(self, index, colormap):
        self._layers[index]["colormap"] = colormap
        for slice in self._slices :
            slice.layers[index].colormap = colormap
    
    def get_layer_opacity(self, index):
        return self._layers[index]["opacity"]
    
    def set_layer_opacity(self, index, opacity):
        self._layers[index]["opacity"] = opacity
        for slice in self._slices :
            slice.layers[index].opacity = opacity
    
    def render(self):
        self._rwi.Render()
    
    def reset_view(self):
        for slice in self._slices :
            slice.reset_view()
    
        self.notify_observers("center")
    
    def get_mouse_button_tool(self, button) :
        """ Return a triplet (ToolClass, args, kwargs)
        """
        
        return self._mouse_tools[button]
        
    def set_mouse_button_tool(self, button, tool_class, *args, **kwargs) :
        """ Set a tool associated with given button (Left, Middle, or Right),
            with an optional modifier (Shift or Control). Example : Right,
            ShiftLeft. Set tool to None to have no tool connected to the button.
        """
        
        self._mouse_tools[button] = (tool_class, args, kwargs)
        
        for slice in self._slices :
            if tool_class :
                tool = tool_class(*args, **kwargs)
            else :
                tool = None
            slice.set_mouse_button_tool(button, tool)
    
    def get_wheel_tool(self, direction) :
        """ Return the tool associated with a mouse wheel direction (Forward or 
            Backward)
        """
        
        event_name = "MouseWheel%s"%direction
        return self._mouse_tools[event_name]
    
    def set_wheel_tool(self, direction, tool_class, *args, **kwargs) :
        """ Set a tool associated with a mouse wheel direction (Forward or 
            Backward)
        """
        
        self._mouse_tools[direction] = (tool_class, args, kwargs)
        
        for slice in self._slices :
            if tool_class :
                tool = tool_class(*args, **kwargs)
            else :
                tool = None
            slice.set_wheel_tool(direction, tool)
    
    def get_keyboard_tool(self, key):
        return self._keyboard_tools[key]
    
    def set_keyboard_tool(self, key, tool_class, *args, **kwargs):
        self._keyboard_tools[key] = (tool_class, args, kwargs)
        
        for slice in self._slices :
            if tool_class :
                tool = tool_class(*args, **kwargs)
            else :
                tool = None
            slice.set_keyboard_tool(key, tool)
    
    def set_next_window_info(self, info):
        """ Remap the VTK (and underlying) windows. This has to be called when
            reparenting.
        """
       
        self._rwi.GetRenderWindow().SetNextWindowInfo(info)
        self._rwi.GetRenderWindow().WindowRemap()
    
    ##############
    # Properties #
    ##############
    
    def _get_layers(self):
        return self._layers
    
    def _get_slice_mode(self):
        return self._slice_mode
    
    def _set_slice_mode(self, slice_mode):
        if slice_mode not in ["axial", "coronal", "sagittal", "multiplanar"] :
            raise medipy.base.Exception("Unknown slice mode : %s"%(slice_mode,))
        
        old_slice_mode = self._slice_mode
        self._slice_mode = slice_mode
        
        for slice in self._slices :
            self._rwi.GetRenderWindow().RemoveRenderer(slice.renderer)
            slice.unset_rwi(self._rwi)
        
        if old_slice_mode == "multiplanar" and slice_mode != "multiplanar" :
            self._rwi.GetRenderWindow().RemoveRenderer(self._informations_renderer) 
        
        self._slices = []
        self._slices_names = []
        
        names = (["axial", "coronal", "sagittal"] if slice_mode == "multiplanar"
                 else [slice_mode])
        
        # Build and configure the Slice objects
        for name in names :
            world_to_slice = medipy.base.coordinate_system.slices[self._convention][name]
            slice = Slice(world_to_slice, self._layers, 
                self._annotations, self._interpolation, 
                self._display_coordinates,self._scalar_bar_visibility, 
                self._orientation_visibility, slice_mode != "multiplanar"
            )
        
            if slice_mode == "multiplanar" :
                slice.renderer.SetViewport(*self._viewport[name])
            else :
                slice.renderer.SetViewport(0., 0., 1., 1.)
            self._rwi.GetRenderWindow().AddRenderer(slice.renderer)
            slice.setup_rwi(self._rwi)
            self._slices.append(slice)
            self._slices_names.append(name)
        
        for index, slice in enumerate(self._slices) :
            for event in slice.allowed_events :
                if event in ["any", "cursor_position", "image_position", "center",
                             "corner_annotations_visibility", "world_to_slice"] :
                    continue
                slice.add_observer(event, self._on_slice_event)
            
            slice.add_observer("cursor_position", self._on_cursor_position)
            slice.add_observer("center", self._on_center)
            
            # Synchronize the Slices' colormaps
            # Do not use PropertySynchronized API for better event handling
            # and thus better performances
            if len(self._slices)>1 :
                next_slice = self._slices[(index+1)%len(self._slices)]
                
                for layer_index, layer in enumerate(slice.layers) :
                    for event in ["data", "display_range", "cut_low", "cut_high", "zero_transparency"] :
                        layer.colormap.add_observer(event, self._on_colormap_event)
                    next_slice_layer = next_slice.layers[layer_index]
                    layer.colormap.append_child(next_slice_layer.colormap)
        
        # Add a 4th renderer with image info in multiplanar mode, otherwise
        # display the corner annotations
        if slice_mode == "multiplanar" :
            self._rwi.GetRenderWindow().AddRenderer(self._informations_renderer)
        else :
            self._slices[0].corner_annotations_visibility = True
        
        for button, (class_, args, kwargs) in self._mouse_tools.items() :
            self.set_mouse_button_tool(button, class_, *args, **kwargs)
        for key, (class_, args, kwargs) in self._keyboard_tools.items() :
            self.set_keyboard_tool(key, class_, *args, **kwargs)
        
        self._update_informations()
        
    def _get_number_of_layers(self):
        return len(self._layers)      
    
    def _get_interpolation(self) :
        return self._interpolation
    
    def _set_interpolation(self, interpolation):
        self._set_slice_property("interpolation", interpolation)
    
    def _get_display_coordinates(self) :
        return self._display_coordinates
    
    def _set_display_coordinates(self, display_coordinates):
        if display_coordinates not in ["physical", "nearest_axis_aligned", "index"] :
            raise medipy.base.Exception("Unknown display coordinates : %s"%(display_coordinates,))
        
        self._set_slice_property("display_coordinates", display_coordinates)
    
    def _get_scalar_bar_visibility(self) :
        return self._scalar_bar_visibility
    
    def _set_scalar_bar_visibility(self, scalar_bar_visibility):
        self._set_slice_property("scalar_bar_visibility", scalar_bar_visibility)
    
    def _get_orientation_visibility(self) :
        return self._orientation_visibility
    
    def _set_orientation_visibility(self, orientation_visibility):
        self._set_slice_property("orientation_visibility", orientation_visibility)
    
    def _get_convention(self):
        """ Image viewing convention, can be either :
              * radiological (left-is-right)
              * neurological (left-is-left)
        """
        
        return self._convention
    
    def _set_convention(self, convention):
        if convention not in ["radiological", "neurological"] :
            raise medipy.base.Exception("Unknown viewing convention : {0}".format(
                convention))
        self._convention = convention
        
        for index, slice in enumerate(self._slices) :
            name = self._slices_names[index]
            world_to_slice = medipy.base.coordinate_system.slices[self._convention][name]
            slice.world_to_slice = world_to_slice
    
    def _get_annotations(self):
        return self._annotations
    
    def _set_annotations(self, annotations):
        
        if self._annotations is not None :
            self._annotations.remove_observer("any", self._on_annotations_changed)
            
        self._annotations = annotations
        self._annotations.add_observer("any", self._on_annotations_changed)
        
        for slice in self._slices :
            slice.annotations = annotations
    
    def _get_cursor_physical_position(self) :
        return self._cursor_physical_position
    
    def _set_cursor_physical_position(self, cursor_physical_position):
        self._cursor_physical_position = self.cursor_physical_position
        if self._layers :
            image = self._layers[0]["image"]
            self._cursor_index_position = (cursor_physical_position-image.origin)/image.spacing
        else :
            self._cursor_index_position = cursor_physical_position
        
        if self._slices :
            # All slices are synchronized
            self._slices[0].cursor_physical_position = cursor_physical_position
        
        self._rwi.Render()
        
        self.notify_observers("cursor_position")
        
    def _get_cursor_index_position(self) :
        return self._cursor_index_position
    
    def _set_cursor_index_position(self, cursor_index_position):
        if self._layers :
            image = self._layers[0]["image"]
            cursor_physical_position = image.origin+cursor_index_position*image.spacing
        else :
            cursor_physical_position = cursor_index_position
        
        self._set_cursor_physical_position(cursor_physical_position)
    
    def _get_center_physical_position(self):
        if self._center_physical_position is not None :
            return self._center_physical_position
        elif self._slices :
            self._center_physical_position = self._slices[0].image_physical_position
            return self._center_physical_position
        else :
            return None
    
    def _set_center_physical_position(self, position):
        self._center_physical_position = position
        
        for slice in self._slices :
            slice.remove_observer("center", self._on_center)
        for slice in self._slices :
            slice.center_on_physical_position(position)
        for slice in self._slices :
            slice.add_observer("center", self._on_center)
        
        self._cursor_physical_position = position
        if self._layers :
            image = self._layers[0]["image"]
            self._cursor_index_position = (position-image.origin)/image.spacing
        else :
            self._cursor_index_position = position
        self._center_index_position = self._cursor_index_position
            
        self.notify_observers("center")
    
    def _get_center_index_position(self):
        if self._center_index_position is not None :
            return self._center_index_position
        elif self._slices :
            return self._slices[0].image_index_position
        else :
            return None
    
    def _set_center_index_position(self, position):
        self._center_index_position = position
        
        for slice in self._slices :
            slice.remove_observer("center", self._on_center)
        for slice in self._slices :
            slice.center_on_index_position(position)
        for slice in self._slices :
            slice.add_observer("center", self._on_center)
        
        self._cursor_index_position = position
        if self._layers :
            image = self._layers[0]["image"]
            self._cursor_physical_position = image.origin+position*image.spacing
        else :
            self._cursor_physical_position = position
        self._center_physical_position = self._cursor_physical_position
    
    def _get_zoom(self) :
        return self._zoom
    
    def _set_zoom(self, zoom):
        self._set_slice_property("zoom", zoom)
    
    layers = property(_get_layers)
    slice_mode = property(_get_slice_mode, _set_slice_mode)
    number_of_layers = property(_get_number_of_layers) 
    interpolation = property(_get_interpolation, _set_interpolation)
    display_coordinates = property(_get_display_coordinates, _set_display_coordinates)
    scalar_bar_visibility = property(_get_scalar_bar_visibility, _set_scalar_bar_visibility)
    orientation_visibility = property(_get_orientation_visibility, _set_orientation_visibility)
    convention = property(_get_convention, _set_convention)
    annotations = property(_get_annotations, _set_annotations)
    cursor_physical_position = property(_get_cursor_physical_position, _set_cursor_physical_position)
    cursor_index_position = property(_get_cursor_index_position, _set_cursor_index_position)
    center_physical_position = property(_get_center_physical_position, _set_center_physical_position)
    center_index_position = property(_get_center_index_position, _set_center_index_position)
    zoom = property(_get_zoom, _set_zoom)
    
    #####################
    # Private interface #
    #####################
    
    def _set_slice_property(self, name, value):
        """ Set a property across all slices, and notify observers.
        """
        
        setattr(self, "_{0}".format(name), value)
        
        for slice in self._slices :
            slice.remove_observer(name, self._on_slice_event)
        for slice in self._slices :
            setattr(slice, name, value)
        for slice in self._slices :
            slice.add_observer(name, self._on_slice_event)
            
        self._rwi.Render()
        
        self.notify_observers(name)
    
    def _on_slice_event(self, event) :
        
        # Don't use _set_xxx to avoid spurious events
        value = getattr(event.object, event.event)
        setattr(self, "_%s"%event.event, value)
        
        for slice in self._slices :
            slice.remove_observer(event.event, self._on_slice_event)
        
        for slice in self._slices :
            if slice != event.object :
                setattr(slice, event.event, value)
        self.notify_observers(event.event)
        
        for slice in self._slices :
            slice.add_observer(event.event, self._on_slice_event)
        
        if event.event == "zoom" :
            self._update_informations()
    
    def _on_colormap_event(self, event):
        if self._slices :
            try :
                layers_colormaps = [x["colormap"] for x in self._layers]
                layer_index = layers_colormaps.index(event.object)
            except ValueError :
                logging.warning("No such colormap in layers {0}".format(self)) 
            else :
                self.notify_observers("colormap_{0}".format(event.event), 
                                      layer_index=layer_index)
    
    def _on_cursor_position(self, event):
        self._cursor_index_position = event.object.cursor_index_position
        self._cursor_physical_position = event.object.cursor_physical_position
        
        for slice in self._slices :
            slice.remove_observer(event.event, self._on_cursor_position)
        
        for slice in self._slices :
            if slice != event.object :
                slice.cursor_physical_position = self._cursor_physical_position
        self.notify_observers(event.event)
        
        for slice in self._slices :
            slice.add_observer(event.event, self._on_cursor_position)
        
        self._update_informations()
    
    def _on_annotations_changed(self, event):
        self.render()
    
    def _on_center(self, event):
        
        for slice in self._slices :
            slice.remove_observer(event.event, self._on_center)
        
        for slice in self._slices :
            if slice != event.object :
                slice.center_on_physical_position(self._cursor_physical_position)
        self.notify_observers(event.event)
        
        for slice in self._slices :
            slice.add_observer(event.event, self._on_center)
        
        self._update_informations()
    
    def _on_button_down(self,event):
        """ Propagate mouse clicks to parent.
        """
        
        event.Skip()
        new_event = event.Clone()
        new_event.SetEventObject(self)
        self.AddPendingEvent(new_event)
    
    def _on_close(self, event):
        for slice in self._slices :
            slice.close()
            slice.unset_rwi(self._rwi)
        self._rwi.Disable()
        self._rwi.Close()
        self.Destroy()
    
    def _update_informations(self) :
        informations = get_informations(self)

        if self._slice_mode == "multiplanar" :
            for where, label in informations.items() :
                index = self._corner_annotation_index[where]
                self._informations_corner_annotation.SetText(index, label)
        else :
            for where, label in informations.items() :
                self._slices[0].set_label(where, label)
