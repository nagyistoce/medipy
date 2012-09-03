##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import os

from vtk import (vtkDataReader, vtkPolyDataReader, vtkPolyDataWriter,
    vtkRectilinearGridReader, vtkStructuredGridReader, 
    vtkUnstructuredGridReader, vtkVRMLImporter,)
import wx

import medipy
import medipy.base
import medipy.gui.image
import medipy.gui.dicom
import medipy.gui.dicom.reconstruction

from medipy.gui.annotations.annotations_dialog import AnnotationsDialog
from medipy.gui.base import Application
import medipy.gui.control
from medipy.gui.image.cine_dialog import CineDialog
from medipy.gui.function_gui_builder import FunctionGUIBuilder
from medipy.gui.image.layers_dialog import LayersDialog
from medipy.gui.image.tools import MouseTool

from medipy.base import ObservableList, Object3D

from main_frame import MainFrame
import menu_builder

class MediPyApp(Application) :
    
    def __init__(self, *args, **kwargs):
        
        # Public interface
        self.images = ObservableList()
        self.gui_images = ObservableList()
        self.viewer_3ds = ObservableList()
        
        # Properties
        self._display_coordinates = None
        self._display_convention = None
        self._synchronize_images = {}
        
        # Private interface
        self._frame = None
        self._active_image_index = None
        self._full_screen = False
        self._cine_dialog = None
        self._layers_dialog = None
        self._annotations_dialog = None
        
        super(MediPyApp, self).__init__(*args, **kwargs)
    
    ####################
    # Public interface #
    ####################
        
    def insert_image(self, index, image, layers=None):  
        """ Image and layers must be three dimensional
        """
        
        layers = layers or []
        
        self._frame.Freeze()
        
        self.images.insert(index, image)
        
        if len(image.shape) == 2 or len(image.shape) == 3 and image.shape[0] == 1 :
            mode = "axial"
        else : 
            mode = "multiplanar"
        
        gui_image = medipy.gui.image.Image(
            self._frame.images_panel, mode, 
            [{"image" : image}] + [{"image" : l} for l in layers],
            image.annotations, display_coordinates=self._display_coordinates,
            convention=self._display_convention
        )
        gui_image.reset_view()
        gui_image.Bind(wx.EVT_LEFT_DOWN, self.OnImageClicked)
        gui_image.Bind(wx.EVT_MIDDLE_DOWN, self.OnImageClicked)
        gui_image.Bind(wx.EVT_RIGHT_DOWN, self.OnImageClicked)
        gui_image.Bind(wx.EVT_MOUSEWHEEL, self.OnImageClicked)
        
        if self.gui_images :
            reference_image = self.gui_images[0]
            if self._synchronize_images["cursor_position"] :
                gui_image.cursor_physical_position = reference_image.cursor_physical_position
            if self._synchronize_images["center"] :
                gui_image.center_physical_position = reference_image.center_physical_position
            if self._synchronize_images["zoom"] :
                gui_image.zoom = reference_image.zoom
            if self._synchronize_images["display_range"] :
                display_range = reference_image.get_layer_colormap(0).display_range
                gui_image.get_layer_colormap(0).display_range = display_range
        
        gui_image.add_observer("cursor_position", self._on_cursor_position)
        gui_image.add_observer("center", self._on_center)
        gui_image.add_observer("colormap_display_range", self._on_colormap_display_range)
        
#        gui_image.set_mouse_button_tool("Right", ContextMenuImageTool, self._frame, "Right")
        
        self.gui_images.insert(index, gui_image)
        
        self.active_image = gui_image
        
        self._frame.insert_image(index, gui_image)
        
        # Set which attributes are synchronized, set them to another image if
        # they are to be synchronized
        gui_image.synchronize_on_event["zoom"] = self._synchronize_images["zoom"]
            
        # Set children for synchronization :
        for image_index, gui_image in enumerate(self.gui_images) :
            next_image = self.gui_images[(image_index+1)%len(self.gui_images)]
            gui_image.delete_all_children()
            gui_image.append_child(next_image)
        
        self._frame.Thaw()
    
    def append_image(self, image, layers = None):
        return self.insert_image(len(self.gui_images), image, layers)
    
    def insert_viewer_3d(self, index, viewer_3d):
        self.viewer_3ds.insert(index, viewer_3d)
        viewer_3d.add_observer("close", self._on_viewer_3d_close)
    
    def append_viewer_3d(self, viewer_3d):
        return self.insert_viewer_3d(len(self.viewer_3ds), viewer_3d)
    
    def insert_layer(self, image_index, layer_index, layer_data):
        self.gui_images[image_index].insert_layer(self, layer_index, 
                                                  **layer_data)
    
    def append_layer(self, image_index, layer_data):
        self.insert_layer(image_index, 
                          self.gui_images[image_index].number_of_layers, 
                          layer_data)
    
    def set_cine_dialog_visibility(self, visibility) :
        if visibility :
            if self._active_image_index is not None :
                image = self.gui_images[self._active_image_index]
                if self._cine_dialog is None :
                    self._cine_dialog = CineDialog(self._frame, wx.ID_ANY)
                self._cine_dialog.image = image
                self._cine_dialog.Show()
            elif self._cine_dialog is not None :
                self._cine_dialog.Hide()
        else :
            self._cine_dialog.Hide()
    
    def toggle_synchronize_images(self, attribute) :
        self.set_synchronize_images(attribute,
                                    not self.get_synchronize_images(attribute))
    
    def toggle_full_screen(self):
        self._full_screen = not self._full_screen
        
        if self._full_screen : 
            self._frame.full_screen(self.active_image)
        else :
            self._frame.full_screen(None)      
    
    def close_image(self, gui_image):
        
        dialogs = ["cine", "layers", "annotations"]
        for name in dialogs :
            dialog = getattr(self, "_%s_dialog"%name)
            if dialog is not None :
                dialog.Destroy()
                setattr(self, "_%s_dialog"%name, None)
        
        gui_image_index = self.gui_images.index(gui_image)
        image_index = self.images.index(gui_image.get_layer_image(0))
        
        self._frame.delete_image(gui_image)
        gui_image.Close()
        del self.gui_images[gui_image_index]
        
        if image_index is not None :
            del self.images[image_index]
        
        if len(self.images) > 0 and self._active_image_index >= len(self.images) : 
            self._set_active_image(self.gui_images[-1])
        elif len(self.images) > 0 :
            self._set_active_image(self.gui_images[self._active_image_index])
        
        # Set children for synchronization :
        for image_index, gui_image in enumerate(self.gui_images) :
            next_image = self.gui_images[(image_index+1)%len(self.gui_images)]
            gui_image.delete_all_children()
            gui_image.append_child(next_image)
    
    def close_all_images(self) :
        while len(self.gui_images)>0 : 
            self.close_image(self.gui_images[0])
    
    def set_layers_dialog_visibility(self, visibility):
        if visibility :
            if self._active_image_index is not None :
                image = self.gui_images[self._active_image_index]
                if self._layers_dialog is None :
                    self._layers_dialog = LayersDialog(self._frame, wx.ID_ANY)
                self._layers_dialog.image = image
                self._layers_dialog.Show()
            elif self._layers_dialog is not None :
                self._layers_dialog.Hide()
        else :
            self._layers_dialog.Hide()
    
    def set_annotations_dialog_visibility(self, visibility):
        if visibility :
            if self._active_image_index is not None :
                image = self.gui_images[self._active_image_index]
                if self._annotations_dialog is None :
                    self._annotations_dialog = AnnotationsDialog(
                        self._frame, title="Annotations", 
                        style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
                self._annotations_dialog.image = image
                self._annotations_dialog.Show()
            elif self._annotations_dialog is not None :
                self._annotations_dialog.Hide()
        else :
            self._annotations_dialog.Hide()
    
    def reset_view(self):
        self.active_image.reset_view()
        self.active_image.render()
    
    def current_image_screenshot(self):
        # Set image to inactive, avoid green border
        self.active_image.active = False
        self.active_image.Refresh()
        self.active_image.Update()
        
        image = self._image_from_window(self.active_image)
        
        # Reset the active status
        self.active_image.active = True
        self.active_image.Refresh()
        self.active_image.Update()
        
        return image
    
    def all_images_screenshot(self):
        return self._image_from_window(self._frame.images_panel)
    
    def whole_window_screenshot(self):
        return self._image_from_window(self._frame)
    
    def execute_script(self, filename):
        execfile(filename, globals(), locals())
    
    def quit(self):
        self.close_all_images()
        for viewer in self.viewer_3ds :
            viewer.Close()
    
    def load_object_3d(self, path, viewer_3d):
        generic_reader = vtkDataReader()
        generic_reader.SetFileName(path)
        
        if generic_reader.OpenVTKFile() and generic_reader.ReadHeader() :
            if generic_reader.IsFileStructuredPoints() :
                raise Exception("Cannot read VTK structured points")
            elif generic_reader.IsFilePolyData() :
                reader = vtkPolyDataReader()
                reader.SetFileName(path)
                object_3d = Object3D(reader.GetOutput(), path)
            elif generic_reader.IsFileStructuredGrid() :
                reader = vtkStructuredGridReader()
                reader.SetFileName(path)
                object_3d = Object3D(reader.GetOutput(), path)
            elif generic_reader.IsFileUnstructuredGrid() :
                reader = vtkUnstructuredGridReader()
                reader.SetFileName(path)
                object_3d = Object3D(reader.GetOutput(), path)
            elif generic_reader.IsFileRectilinearGrid() :
                reader = vtkRectilinearGridReader()
                reader.SetFileName(path)
                object_3d = Object3D(reader.GetOutput(), path)
            else : 
                raise Exception("Cannot read VTK file containing type %i"%generic_reader.GetFileType())
            viewer_3d.objects_3d.append(object_3d)
        else :
            importer = vtkVRMLImporter()
            importer.SetFileName(path)
            importer.Update()
            actors = importer.GetRenderer().GetActors()
            number_of_actors = actors.GetNumberOfItems()
            actors.InitTraversal()
            for i in range(number_of_actors) :
                actor = actors.GetNextItem()
                object_3d = Object3D(path + ", %i"%i)
                object_3d.actor.SetProperty(actor.GetProperty())
                object_3d.actor.SetMapper(actor.GetMapper())
                viewer_3d.objects_3d.append(object_3d)
#        else :
#            raise Exception("Cannot load file %s of unknown type"%path)
    
    def save_object_3d(self, object_3d, path) :
        writer = vtkPolyDataWriter()
        writer.SetFileName(path)
        writer.SetFileTypeToBinary()
        writer.SetHeader(object_3d.name)
        writer.SetInput(object_3d.dataset)
        writer.Update()
    
    def get_synchronize_images(self, attribute):
        return self._synchronize_images[attribute]
    
    def set_synchronize_images(self, attribute, value) :
        self._synchronize_images[attribute] = value
        wx.ConfigBase_Get().WriteBool("SynchronizeImages_%s"%attribute, value)
        wx.ConfigBase_Get().Flush()
        
        if attribute == "zoom" :
            for image in self.gui_images :
                if hasattr(image, attribute) :
                    image.synchronize_on_event[attribute] = value
    
    ##############
    # PROPERTIES #
    ##############

    def _get_display_coordinates(self) :
        return self._display_coordinates
    
    def _set_display_coordinates(self, value) :
        self._display_coordinates = value
        for image in self.gui_images :
            image.display_coordinates = value
            image.render()
        wx.ConfigBase_Get().Write("DisplayCoordinates", self._display_coordinates)
        wx.ConfigBase_Get().Flush()
    
    def _get_display_convention(self) :
        return self._display_convention
    
    def _set_display_convention(self, value) :
        self._display_convention = value
        for image in self.gui_images :
            image.convention = value
            image.render()
        wx.ConfigBase_Get().Write("DisplayConvention", self._display_convention)
        wx.ConfigBase_Get().Flush()

    def _set_slices(self, slices):
        if slices != self.active_image.slice_mode :
            self.active_image.slice_mode = slices
            self.active_image.render()
        self._frame.slices = slices
    
    def _get_active_image(self):
        """ Active GUI image
        """
        
        if self._active_image_index != -1 :
            return self.gui_images[self._active_image_index]
        else :
            return None
    
    def _set_active_image(self, gui_image):
        try :
            self._active_image_index = self.gui_images.index(gui_image)
        except ValueError :
            # Image is not in list, do nothing
            pass

        # Display active image path in title
        url = self.images[self._active_image_index].metadata.get("loader", {}).get("url", "")
        self._frame.SetTitle("MediPy ({0})".format(url))

        for other_image in self.gui_images :
            if other_image == gui_image :
                other_image.SetBackgroundColour(wx.GREEN)
            else :
                other_image.SetBackgroundColour(wx.BLACK)
        self._frame.slices = gui_image.slice_mode
        if self._cine_dialog :
            self._cine_dialog.image = gui_image
        if self._layers_dialog :
            self._layers_dialog.image = gui_image
        if self._annotations_dialog is not None :
            self._annotations_dialog.image = gui_image
            self._annotations_dialog.annotations = gui_image.annotations
        
        if (isinstance(self._frame.current_ui, FunctionGUIBuilder) and 
            hasattr(self._frame.current_ui, "controls")) :
            for control in self._frame.current_ui.controls.values() :
                if isinstance(control, medipy.gui.control.Coordinates) :
                    control.image = gui_image
        elif hasattr(self._frame.current_ui, "image") :
            self._frame.current_ui.image = gui_image
    
    display_coordinates = property(_get_display_coordinates, 
                                   _set_display_coordinates)
    display_convention = property(_get_display_convention,
                                  _set_display_convention)
    slices = property(lambda self:self.active_image.slices, _set_slices)
    active_image = property(_get_active_image, _set_active_image)
    
    ##################
    # Event handlers #
    ##################
    def OnInit(self) :
        self.SetAppName("MediPy")
        
        for attribute in ["cursor_position", "center", "zoom", "display_range"] :
            config_entry = "SynchronizeImages_%s"%attribute
            value = wx.ConfigBase_Get().ReadBool(config_entry) or False
            self.set_synchronize_images(attribute, value)
        
        self.display_coordinates = (
            wx.ConfigBase_Get().Read("DisplayCoordinates") or "physical")
        self.display_convention = (
            wx.ConfigBase_Get().Read("DisplayConvention") or "radiological")
        
        if self.options.menu_file is not None :
            menu = menu_builder.from_file.build_menu(self.options.menu_file)
        else :
            menu = []
            for directory in medipy.Importer().plugins_path :
                menu.extend(menu_builder.from_api.build_menu(directory))
        
        self._frame = MainFrame(None, wx.ID_ANY, menu, 
                                title="MediPy", size=(1000,800))
        self._frame.application = self
        self._frame.Show()
        self.SetTopWindow(self._frame)
        
        return True
    
    def OnImageClicked(self, event):
        self.active_image = event.GetEventObject()
    
    def _on_cursor_position(self, event) :
        
        if not self._synchronize_images.get("cursor_position", False) :
            return
        
        for gui_image in self.gui_images :
            if gui_image != event.object :
                gui_image.remove_observer("cursor_position", self._on_cursor_position)
                gui_image.cursor_physical_position = event.object.cursor_physical_position
                gui_image.add_observer("cursor_position", self._on_cursor_position)
        
    def _on_center(self, event) :
        
        if not self._synchronize_images.get("center", False) :
            return
        
        for gui_image in self.gui_images :
            if gui_image != event.object :
                gui_image.remove_observer("center", self._on_center)
                gui_image.center_physical_position = event.object.center_physical_position
                gui_image.render()
                gui_image.add_observer("center", self._on_center)
    
    def _on_colormap_display_range(self, event):
        
        if not self._synchronize_images.get("display_range", False) :
            return
        
        for gui_image in self.gui_images :
            if gui_image != event.object :
                gui_image.remove_observer("colormap_display_range", self._on_colormap_display_range)
                source_colormap = event.object.get_layer_colormap(event.layer_index)
                destination_colormap = gui_image.get_layer_colormap(event.layer_index)
                destination_colormap.display_range = source_colormap.display_range
                gui_image.render()
                gui_image.add_observer("colormap_display_range", self._on_colormap_display_range)
    
    def _on_viewer_3d_close(self, event):
        viewer = event.object
        index = self.viewer_3ds.index(viewer)
        del self.viewer_3ds[index]
    
    #####################
    # Private interface #
    #####################
    
    def _image_from_window(self, window):
        # Force a render of the images to remove menu on screenshot
        for image in self.gui_images :
            image.render()
        
        # Take the screenshot by rendering to memory
        self._frame.Refresh()
        self._frame.Update()
        screen = wx.ClientDC(window)
        size = window.GetSize()
        image = wx.EmptyBitmap(size[0], size[1], -1)
        memory = wx.MemoryDC()
        memory.SelectObject(image)
        memory.Blit(0, 0, size.GetWidth(), size.GetHeight(), screen, 0, 0)
        return image
