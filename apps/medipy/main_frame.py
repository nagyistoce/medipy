##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import inspect
import os
import sys

import numpy
import wx
import wx.xrc

import medipy.base
import medipy.gui.annotations
import medipy.gui.base
import medipy.gui.function_gui_builder
import medipy.gui.image
import medipy.gui.io
import medipy.gui.image.mouse_tools
import medipy.gui.image.tools
from medipy.gui.viewer_3d_frame import Viewer3DFrame

import medipy.diffusion

import menu_builder
from medipy.base.observable import Observable

class SelectMainTool(medipy.gui.image.tools.KeyboardTool) :
    def __init__(self, frame):
        super(SelectMainTool, self).__init__()
        self.frame = frame
    
    def press(self, rwi, slice):
        key = rwi.GetKeyCode() or rwi.GetKeySym()
        tools = {
            "s" : medipy.gui.image.mouse_tools.Select,
            "p" : medipy.gui.image.mouse_tools.Pan,
            "c" : medipy.gui.image.mouse_tools.WindowLevel,
        }
        if key in tools :
            self.frame.set_image_tool(tools[key])

class MainFrame(medipy.gui.base.Frame):
    
    class UI(medipy.gui.base.UI):
        def __init__(self):
            self.menu_treectrl = None
            self.function_ui_panel = None
            self.function_ui_sizer = None
            
            self.layers_tab = None
            self.layers = None
            
            self.annotations_tab = None
            self.annotations = None
            
            self.images_panel = None
            self.image_grid = None
            self.images_sizer = None
            
            self.controls = [
                "menu_treectrl", "function_ui_panel",
                "layers_tab", 
                "annotations_tab",
                "images_panel"]
                
            medipy.gui.base.UI.__init__(self)
        
        def from_window(self, window, names):
            medipy.gui.base.UI.from_window(self, window, names)
            
            self.layers = medipy.gui.image.LayersPanel(self.layers_tab)
            layers_sizer = wx.BoxSizer()
            layers_sizer.Add(self.layers, 1, wx.EXPAND)
            self.layers_tab.SetSizer(layers_sizer)
            
            self.annotations = medipy.gui.annotations.AnnotationsPanel(self.annotations_tab)
            annotations_sizer = wx.BoxSizer()
            annotations_sizer.Add(self.annotations, 1, wx.EXPAND)
            self.annotations_tab.SetSizer(annotations_sizer)
            
            self.image_grid = medipy.gui.image.ImageGrid(self.images_panel)
            
            # Setup sizers
            self.function_ui_sizer = wx.BoxSizer(wx.VERTICAL)
            self.function_ui_panel.SetSizer(self.function_ui_sizer)
            
            sizer = wx.BoxSizer()
            sizer.Add(self.image_grid, 1, wx.EXPAND)
            self.images_panel.SetSizer(sizer)
    
    def __init__(self, left_menu, parent=None, *args, **kwargs):
        self.ui = MainFrame.UI()

        xrc_file = medipy.base.find_resource(
            os.path.join("resources", "gui", "medipy_frame.xrc"))
        medipy.gui.base.Frame.__init__(self, xrc_file, "medipy_frame", 
            [], self.ui, self.ui.controls, parent, *args, **kwargs)
        
        self._current_ui = None
        self._full_screen_frame = None
        self._full_screen_image_index = None
        self._full_screen_image = None
        self._menus_active_when_image_loaded = [
            "save_image_as", "save_image_serie_as", "close_image",
            "current_image_screenshot", "all_images_screenshot",
            "whole_window_screenshot"
        ]
        
        self.images = medipy.base.ObservableList()
        self._preferences = medipy.gui.base.Preferences(
            wx.GetApp().GetAppName(), wx.GetApp().GetVendorName())
        
        self._display_coordinates = None
        self._display_convention = None
        self._slices = None
        self._crosshair = None
        
        self._synchronization = {
            "cursor_position" : False,
            "center" : False,
            "zoom" : False,
            "display_range" : False
        }
        
        self._image_tool = None
        
        ##################
        # Initialize GUI #
        ##################
        
        self.SetTitle(wx.GetApp().GetAppName())
        
        # Bind events
        self.ui.menu_treectrl.Bind(wx.EVT_TREE_SEL_CHANGED, 
                                   self.OnMenuTreeCtrlSelChanged)
        self.ui.image_grid.add_observer("active", self._on_active_image)
        
        # Fill function menu
        menu_builder.fill_treectrl(left_menu, self.ui.menu_treectrl)
        
        # Disable item that should not be active when no image is loaded
        for item in self._menus_active_when_image_loaded :
            menu_item = self.GetMenuBar().FindItemById(wx.xrc.XRCID(item))
            if menu_item is not None :
                menu_item.Enable(False)
        
        # Update UI from preferences
        self.display_coordinates = self._preferences.get("Display/coordinates", "index")
        self.display_convention = self._preferences.get("Display/convention", "radiological")
        self.tensor2_display_mode = self._preferences.get("Display/tensor2", "principal_direction_voxel")
        self.slices = self._preferences.get("Display/slices", "axial")
        self.crosshair = self._preferences.get("Display/crosshair", "full")
        for attribute in self._synchronization :
            value = self._preferences.get(
                "Display/Synchronization/{0}".format(attribute), False)
            self.set_synchronization(attribute, value)
        tool = self._preferences.get("Display/Tools/left", "medipy.gui.image.mouse_tools.Select")
        if "." in tool :
            module, class_ = tool.rsplit(".", 1)
            module = sys.modules[module]
            tool = getattr(module, class_)
        else :
            tool = locals()[tool]
        self.set_image_tool(tool)


    ####################
    # Public interface #
    ####################

    def append_image(self, layers) :
        nb_elements = len(self.images)
        self.insert_image(nb_elements, layers)
    
    def insert_image(self, index, layers):
        if len(layers) == 1 :
            self.images.insert(index, layers[0]["image"])
        else :
            self.images.insert(index, [layer["image"] for layer in layers])
        
        self.ui.image_grid.insert(index, layers)
        
        image = self.ui.image_grid[index]
        image.set_mouse_button_tool("Left", 
            self._image_tool[0], *self._image_tool[1], **self._image_tool[2])
        image.set_mouse_button_tool("Middle", None)
        image.set_mouse_button_tool("Right", None)
        image.set_keyboard_tool("s", SelectMainTool, self)
        image.set_keyboard_tool("p", SelectMainTool, self)
        image.set_keyboard_tool("c", SelectMainTool, self)
        
        for index in range(len(image.layers)) :
            if image.get_layer_class(index) == medipy.diffusion.gui.Tensor2Layer :
                image.set_layer_property(index, "display_mode", self.tensor2_display_mode)
        
        # Update menu items
        for item in self._menus_active_when_image_loaded :
            menu_item = self.GetMenuBar().FindItemById(wx.xrc.XRCID(item))
            if menu_item is not None :
                menu_item.Enable(True)
    
    def delete_image(self, index):
        for name in ["layers", "annotations"] :
            panel = getattr(self.ui, name)
            panel.image = None
        self.ui.image_grid.delete(index)
        del self.images[index]
        
        # Update menu items
        for item in self._menus_active_when_image_loaded :
            menu_item = self.GetMenuBar().FindItemById(wx.xrc.XRCID(item))
            if menu_item is not None :
                menu_item.Enable(len(self.images)!=0)
    
#    def full_screen(self, gui_image):
#        if gui_image is not None :
#            # Create a top-level frame to display the image
#            self._full_screen_frame = wx.Frame(None)
#            
#            functions = {
#                (0, wx.WXK_ESCAPE) : lambda event : wx.GetApp().toggle_full_screen(),
#                (wx.ACCEL_ALT, ord("a")) : self.OnViewAxial,
#                (wx.ACCEL_ALT, ord("c")) : self.OnViewCoronal,
#                (wx.ACCEL_ALT, ord("s")) : self.OnViewSagittal,
#                (wx.ACCEL_ALT, ord("m")) : self.OnViewMultiplanar,
#                # TODO : not working ?
#                #(wx.ACCEL_CTRL, ord("0")) : self.OnResetView
#            }
#            
#            entries = []
#            for parameters, action in functions.items() :
#                command_id = wx.NewId()
#                self._full_screen_frame.Bind(wx.EVT_MENU, action, id=command_id)
#                entry = list(parameters)
#                entry.append(command_id)
#                entries.append(tuple(entry))
#            
#            accelerator_table = wx.AcceleratorTable(entries)
#            self._full_screen_frame.SetAcceleratorTable(accelerator_table)
#
#            children = [c.GetWindow() for c in self.ui.images_sizer.GetChildren()]
#            self._full_screen_image_index = children.index(gui_image)
#            self._full_screen_image = gui_image
#            
#            # Remove gui_image from sizer
#            self.ui.images_sizer.Detach(gui_image)
#            
#            # Reparent it to the full screen frame
#            self._reparent_gui_image(gui_image, self._full_screen_frame)
#            
#            # Show the frame full screen
#            self._full_screen_frame.ShowFullScreen(True)
#        elif self._full_screen_frame is not None :
#            self._full_screen_frame.Hide()
#            
#            self._reparent_gui_image(self._full_screen_image, self.ui.images_panel)
#            
#            self.ui.images_sizer.Insert(self._full_screen_image_index, self._full_screen_image, 1, wx.EXPAND)
#            
#            # This crude hack forces the Sizer to correctly redraw its children
#            size = self.GetSize()
#            size[0] -= 1
#            self.SetSize(size)
#            
#            self._full_screen_image_index = None
#            self._full_screen_image = None
#            
#            self._full_screen_frame.Close()
#            self._full_screen_frame = None
#            
#        else : 
#            # remove full screen, but we are not in full screen mode : do nothing
#            pass
    
    def get_synchronization(self, attribute):
        """ Return the synchronization state of an image attribute (may be
            cursor_position, center, zoom or display range).
        """
        
        return self._synchronization[attribute]
    
    def set_synchronization(self, attribute, value):
        """ Set the synchronization state of an image attribute (may be
            cursor_position, center, zoom or display range).
        """
        
        if attribute not in self._synchronization :
            raise medipy.base.Exception("Invalid synchronization attribute: {0}".format(attribute))
        self._synchronization[attribute] = value
        
        self.ui.image_grid.set_synchronization(attribute, value)
        
        # Update GUI
        "Display/Synchronization/{0}".format(attribute)
        
        # Update GUI
        item_id = wx.xrc.XRCID("sync_" + attribute)
        item = self.GetMenuBar().FindItemById(item_id)
        item.Check(value)
        # Update preferences
        self._preferences.set("Display/Synchronization/{0}".format(attribute), value)
    
    def reset_view(self):
        """ Reset the view of the active image. If synchronization is enabled, 
            other images will also reset.
        """
        
        for image in self.ui.image_grid :
            image.reset_view()
            image.render()
    
    def current_image_screenshot(self):
        image = self.ui.image_grid[self.ui.image_grid.active]
        return self._image_from_window(image)
    
    def all_images_screenshot(self):
        return self._image_from_window(self.ui.images_grid)
    
    def whole_window_screenshot(self):
        return self._image_from_window(self)
    
    def get_image_tool(self) :
        return self._image_tool
    
    def set_image_tool(self, tool_class, *args, **kwargs):
        self._image_tool = (tool_class, args, kwargs)
        for image in self.ui.image_grid :
            image.set_mouse_button_tool("Left", tool_class, *args, **kwargs)
        
        class_to_label = {
            medipy.gui.image.mouse_tools.Select : "Select",
            medipy.gui.image.mouse_tools.Pan : "Pan",
            medipy.gui.image.mouse_tools.WindowLevel : "Contrast"
        }
        
        for tool_id in self.tool_ids :
            label = self.GetToolBar().FindById(tool_id).GetLabel()
            if label == class_to_label[tool_class] :
                self.GetToolBar().ToggleTool(tool_id, True)
        
        entry = "{0}.{1}".format(tool_class.__module__, tool_class.__name__)
        self._preferences.set("Display/Tools/left", entry)
        
    ##############
    # Properties #
    ##############
    
    def _get_display_coordinates(self):
        return self._display_coordinates
    
    def _set_display_coordinates(self, value):
        self._display_coordinates = value
        self.ui.image_grid.display_coordinates = value
        
        # Update GUI
        item_id = wx.xrc.XRCID("display_coordinates_" + value)
        item = self.GetMenuBar().FindItemById(item_id)
        item.Check()
        # Update preferences
        self._preferences.set("Display/coordinates", value)
    
    def _get_display_convention(self) :
        return self._display_convention
    
    def _set_display_convention(self, value):
        self._display_convention = value
        self.ui.image_grid.convention = value
        
        # Update GUI
        item_id = wx.xrc.XRCID("display_convention_" + value)
        item = self.GetMenuBar().FindItemById(item_id)
        item.Check()
        # Update preferences
        self._preferences.set("Display/convention", value)
    
    def _get_tensor2_display_mode(self):
        return self._tensor2_display_mode
    
    def _set_tensor2_display_mode(self, value):
        if value not in ["principal_direction_voxel", 
                         "principal_direction_line", "ellipsoid"] :
            raise medipy.base.Exception(
                "Unknown display mode: {0}".format(repr(value)))
        self._tensor2_display_mode = value
        
        for image in self.ui.image_grid :
            for index in range(len(image.layers)) :
                render = False
                if image.get_layer_class(index) == medipy.diffusion.gui.Tensor2Layer :
                    image.set_layer_property(index, "display_mode", self.tensor2_display_mode)
                    render = True
                if render :
                    image.render()
        
        # Update GUI
        item_id = wx.xrc.XRCID("tensor2_display_" + value)
        item = self.GetMenuBar().FindItemById(item_id)
        item.Check()
        
        # Update preferences
        self._preferences.set("Display/tensor2", value)
    
    def _get_slices(self) :
        return self._slices
    
    def _set_slices(self, value):
        self._slices = value
        self.ui.image_grid.slice_mode = value
        
        # Update the tensor display mode
        for image in self.ui.image_grid :
            for index in range(len(image.layers)) :
                render = False
                if image.get_layer_class(index) == medipy.diffusion.gui.Tensor2Layer :
                    image.set_layer_property(index, "display_mode", self.tensor2_display_mode)
                    render = True
                if render :
                    image.render()
        
        # Update GUI
        item_id = wx.xrc.XRCID("view_" + value)
        item = self.GetMenuBar().FindItemById(item_id)
        item.Check()
        # Update preferences
        self._preferences.set("Display/slices", value)
    
    def _get_crosshair(self) :
        return self._crosshair

    def _set_crosshair(self, value) :
        if value not in ["full", "none"] :
            raise medipy.base.Exception("Unknown crosshair mode: {0!r}".format(value))
        
        self._crosshair = value
        self.ui.image_grid.crosshair = value
        
        # Update GUI
        item_id = wx.xrc.XRCID("crosshair_" + value)
        item = self.GetMenuBar().FindItemById(item_id)
        item.Check()
        # Update preferences
        self._preferences.set("Display/crosshair", value)
    
    def _get_current_ui(self):
        return self._current_ui
    
    display_coordinates = property(_get_display_coordinates, 
                                   _set_display_coordinates)
    display_convention = property(_get_display_convention,
                                  _set_display_convention)
    tensor2_display_mode = property(_get_tensor2_display_mode, _set_tensor2_display_mode)
    slices = property(_get_slices, _set_slices)
    crosshair = property(_get_crosshair, _set_crosshair)
    current_ui = property(_get_current_ui)
    
    ##################
    # Event handlers #
    ##################    
    
    def OnOpenImage(self, dummy):
        
        images = medipy.gui.io.load(self, multiple=True)
        
        for image in images :
            # TODO : use force_3d as parameter ?
            if image.ndim == 2 :
                image.data = image.data.reshape(1, *image.shape)
                image.origin = numpy.insert(image.origin, 0, 0)
                image.direction = numpy.insert(
                    numpy.insert(image.direction,0, [0,0], 0), 
                    0, [1,0,0],1)
                image.spacing = numpy.insert(image.spacing, 0, 1)
            self.append_image([{"image":image}])

    def OnOpenImageSerie(self, dummy):
        
        images = medipy.gui.io.load(self, multiple=True, load_all_images=True)
        images = [{"image":image} for image in images]

        if isinstance(images, list) and len(images)>0 :
            self.append_image(images)
            image = self.ui.image_grid[-1]
            for index in range(1, len(image.layers)) :
                image.set_layer_visibility(index, False)
                
    def OnOpenRaw(self, dummy):
        raw_dialog = medipy.gui.image.ImportRawDialog(self)
        if raw_dialog.ShowModal() == wx.ID_OK :
            self.append_image([{"image" : raw_dialog.image}])
        
    def OnOpenSpectro(self, dummy):
        from medipy.gui.image.spectro_dialog import SpectroDialog
        spectro_dialog = SpectroDialog(None, size = (500,220))
        spectro_dialog.ShowModal()
    
    def OnFromDirectory(self, dummy):
        images = medipy.gui.io.import_dicom_directory(self)
        if images :
            self.append_image([{"image":image} for image in images])

    def OnFromDirectorySerie(self, dummy):
        images = medipy.gui.io.import_serie_from_dicom_directory(self)
        images = [{"image":image} for image in images]

        if isinstance(images, list) and len(images)>0 :
            self.append_image(images)
            image = self.ui.image_grid[-1]
            for index in range(1, len(image.layers)) :
                image.set_layer_visibility(index, False)
    
    def OnSaveImageAs(self, dummy):
        image = self.images[self.ui.image_grid.active]
        medipy.gui.io.save(image)

    def OnSaveImageSerieAs(self, dummy):
        images = self.images[self.ui.image_grid.active]
        medipy.gui.io.save(images)
    
    def OnCloseImage(self, dummy):
        self.delete_image(self.ui.image_grid.active)
    
    def OnCloseAllImages(self, dummy):
        while self.images :
            self.delete_image(0)
    
    def OnQuit(self, dummy):
        self.Close()
    
#    def OnCineImage(self, dummy) :
#        wx.GetApp().set_cine_dialog_visibility(True)
    
    def OnSyncCursorPosition(self, dummy) :
        self.set_synchronization(
            "cursor_position",not self.get_synchronization("cursor_position"))
    
    def OnSyncCenter(self, dummy):
        self.set_synchronization("center",not self.get_synchronization("center"))
    
    def OnSyncZoom(self, dummy) :
        self.set_synchronization("zoom",not self.get_synchronization("zoom"))
    
    def OnSyncDisplayRange(self, dummy):
        self.set_synchronization(
            "display_range",not self.get_synchronization("display_range"))

    def OnTensor2DisplayPrincipalDirectionVoxel(self, dummy) :
        self.tensor2_display_mode = "principal_direction_voxel"
    
    def OnTensor2DisplayPrincipalDirectionLine(self, dummy) :
        self.tensor2_display_mode = "principal_direction_line"
    
    def OnTensor2DisplayEllipsoid(self, dummy) :
        self.tensor2_display_mode = "ellipsoid"
    
    def OnDisplayCoordinatesPhysical(self, dummy) :
        self.display_coordinates = "physical"
    
    def OnDisplayCoordinatesNearestAxisAligned(self, dummy) :
        self.display_coordinates = "nearest_axis_aligned"
    
    def OnDisplayCoordinatesIndex(self, dummy) :
        self.display_coordinates = "index"
    
    def OnDisplayConventionRadiological(self, dummy) :
        self.display_convention = "radiological"
    
    def OnDisplayConventionNeurological(self, dummy) :
        self.display_convention = "neurological"
    
    def OnSelect(self, dummy):
        self.set_image_tool(medipy.gui.image.mouse_tools.Select)
    
    def OnPan(self, dummy):
        self.set_image_tool(medipy.gui.image.mouse_tools.Pan)
    
    def OnContrast(self, dummy):
        self.set_image_tool(medipy.gui.image.mouse_tools.WindowLevel)
    
    def OnMenuTreeCtrlSelChanged(self, evt):
        # Destroy existing UI
        if self._current_ui is not None :
            if isinstance(self.current_ui, medipy.gui.function_gui_builder.FunctionGUIBuilder) :
                self.current_ui.remove_observer("new_image", self.on_function_gui_builder_new_image)
                self.current_ui.remove_observer("replace_image", self.on_function_gui_builder_new_image)
            
            if hasattr(self._current_ui, "Close") :
                self._current_ui.Close() 
            self.ui.function_ui_sizer.Clear(True)
            del self._current_ui
        
        # Create UI for current item data
        item_data = self.ui.menu_treectrl.GetPyData(evt.GetItem())
        if item_data is not None : 
            if inspect.isfunction(item_data) :
                function = item_data
                function_gui=medipy.gui.function_gui_builder.FunctionGUIBuilder(
                    self.ui.function_ui_panel, function, self.images, 
                    wx.GetApp().viewer_3ds, "function_parameters")
                if self.images :
                    for control in function_gui.controls.values() :
                        if isinstance(control, medipy.gui.control.Coordinates) :
                            control.image = self.ui.image_grid[self.ui.image_grid.active]
                
                function_gui.add_observer("new_image", self.on_function_gui_builder_new_image)
                function_gui.add_observer("replace_image", self.on_function_gui_builder_replace_image)
                
                self._current_ui = function_gui
                self.ui.function_ui_sizer.Add(function_gui.panel, 1, wx.EXPAND)
                self.ui.function_ui_sizer.Layout()
            else :
                ui = item_data(self.ui.function_ui_panel)
                if hasattr(ui, "image") and self.images :
                    ui.image = self.ui.image_grid[self.ui.image_grid.active]
                self._current_ui = ui
                self.ui.function_ui_sizer.Add(ui, 1, wx.EXPAND)
                self.ui.function_ui_sizer.Layout()
        else : 
            self._current_ui = None
    
    def OnCrosshairFull(self, dummy):
        self.crosshair = "full"
    
#    def OnCrosshairPartial(self, dummy):
#        pass
    
    def OnCrosshairNone(self, dummy):
        self.crosshair = "none"
    
    def OnViewAxial(self, dummy):
        self.slices = "axial"
    
    def OnViewCoronal(self, dummy):
        self.slices = "coronal"
    
    def OnViewSagittal(self, dummy):
        self.slices = "sagittal"
        
    def OnViewMultiplanar(self, dummy):
        self.slices = "multiplanar"
    
    def OnResetView(self, dummy):
        self.reset_view()
    
#    def OnFullScreen(self, dummy):
#        wx.GetApp().toggle_full_screen()
    
    def OnCurrentImageScreenshot(self, dummy):
        image = self.current_image_screenshot()
        self._save_image(image)
    
    def OnAllImagesScreenshot(self, dummy):
        image = self.all_images_screenshot()
        self._save_image(image)
    
    def OnWholeWindowScreenshot(self, dummy):
        image = self.whole_window_screenshot()
        self._save_image(image)
    
#    def OnExecuteScript(self, dummy):
#        path = wx.GetApp().preferences.get("IO/execute_script_path", "")
#        
#        file_dialog = wx.FileDialog(self, style=wx.FD_OPEN, defaultDir=path, wildcard="*.py")
#        if file_dialog.ShowModal() == wx.ID_OK :
#            wx.GetApp().execute_script(file_dialog.GetPath())
#            wx.GetApp().preferences.set("IO/execute_script_path", 
#                                        file_dialog.GetDirectory())
    
    def OnNewViewer3d(self, dummy):
        viewer_3d = Viewer3DFrame(None, medipy.base.ObservableList())
        viewer_3d.Show()
        wx.GetApp().append_viewer_3d(viewer_3d)
    
    def on_function_gui_builder_new_image(self, event):
        self.append_image([{"image" : event.image}])
    
    def on_function_gui_builder_replace_image(self, event):
        index = self.images.index(event.old)
        self.delete_image(index)
        self.insert_image(index, [{"image":event.new}])
    
    def _on_active_image(self, event):
        # Display active image path in title
        active_image = self.images[self.ui.image_grid.active]
        if isinstance(active_image, list) :
            active_image = active_image[0]
        url = active_image.metadata.get("loader", {}).get("url", "")
        self.SetTitle("{0} ({1})".format(wx.GetApp().GetAppName(), url))
        
        image = self.ui.image_grid[self.ui.image_grid.active]
        
        # Update Layers, Annotations
        for name in ["layers", "annotations"] :
            panel = getattr(self.ui, name)
            if panel.image != image :
                panel.image = image
        
        # Update Coordinate controls
        if (isinstance(self.current_ui, medipy.gui.function_gui_builder.FunctionGUIBuilder) and 
            hasattr(self.current_ui, "controls")) :
            for control in self.current_ui.controls.values() :
                if isinstance(control, medipy.gui.control.Coordinates) :
                    control.image = image
        elif hasattr(self.current_ui, "image") :
            self.current_ui.image = image
    
    #####################
    # Private interface #
    #####################
    
    def _reparent_gui_image(self, gui_image, parent) :
        orig_render = gui_image.render
        gui_image.render = lambda : None
        if sys.platform == "linux2" :
            gui_image.set_next_window_info(str(parent.GetHandle()))
        gui_image.Reparent(parent)
        wx.SafeYield()
        # Restore the original render.
        gui_image.render = orig_render
        gui_image.render()
    
    def _image_from_window(self, window):
        # Take the screenshot by rendering to memory
        self.Refresh()
        self.Update()
        screen = wx.ClientDC(window)
        size = window.GetSize()
        image = wx.EmptyBitmap(size[0], size[1], -1)
        memory = wx.MemoryDC()
        memory.SelectObject(image)
        memory.Blit(0, 0, size.GetWidth(), size.GetHeight(), screen, 0, 0)
        return image
    
    def _save_image(self, image):
        path = self._preferences.get("IO/save_image_path", "")
        dialog = wx.FileDialog(self, "Save screenshot", 
            defaultDir = path,
            style = wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT,
            wildcard = "PNG format|*.png|JPG format|*.jpg")
        if dialog.ShowModal() == wx.ID_OK :
            filter = dialog.GetWildcard().split("|")[2*dialog.GetFilterIndex()+1]
            extension = "." + filter.split(".")[1]
            if extension == ".png" : 
                format = wx.BITMAP_TYPE_PNG
            elif extension == ".jpg" : 
                format = wx.BITMAP_TYPE_JPEG
            else : 
                raise Exception("Unknown format : %s"%extension)
            
            if not dialog.GetPath().endswith(extension) :
                filename = dialog.GetPath() + extension
            else : 
                filename = dialog.GetPath()
            
            image.SaveFile(filename, format)
            
            self._preferences.set("IO/save_image_path", 
                                        dialog.GetDirectory())
