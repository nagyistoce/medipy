##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import inspect
import logging
import math
import os.path
import sys
import xml.dom.minidom

import numpy
import wx
import wx.xrc

from medipy.base import find_resource, ObservableList
import medipy.io 
import medipy.gui.io
from medipy.gui import xrc_wrapper
from medipy.gui.viewer_3d_frame import Viewer3DFrame
from medipy.gui.utilities import underscore_to_camel_case, remove_access_letter_from_menu

from medipy.gui.function_gui_builder import FunctionGUIBuilder
import menu_builder

class MainFrame(xrc_wrapper.Frame):
    def __init__(self, parent, id_, left_menu, *args, **kwargs):
        
        self._current_ui = None
        self._full_screen_frame = None
        self._full_screen_image_index = None
        self._full_screen_image = None
        self._menus_active_when_image_loaded = [
            "save_image_as", "close_image", 
            "view_axial", "view_coronal", "view_sagittal", "view_multiplanar",
        ]
        
        ##################
        # Initialize GUI #
        ##################
        
        filename = find_resource("resources/gui/medipy_frame.xrc")
        
        xml_document = xml.dom.minidom.parse(filename)
        
        resource = wx.xrc.XmlResource(filename)
        
        frame = resource.LoadFrame(parent, "medipy_frame")
        xrc_wrapper.Frame.__init__(self, frame, id_, *args, **kwargs)
        
        # Get controls
        controls = ["menu_treectrl", "function_ui_panel", "images_panel"]
        for control in controls :
            setattr(self, "_"+control, wx.xrc.XRCCTRL(self, control))
        
        # Setup sizers
        self._function_ui_sizer = wx.BoxSizer(wx.VERTICAL)
        self._function_ui_panel.SetSizer(self._function_ui_sizer)
        self._images_sizer = wx.GridSizer(cols=0, rows=0, vgap=0, hgap=0)
        self._images_panel.SetSizer(self._images_sizer)
        
        # Bind widgets events
        self._menu_treectrl.Bind(wx.EVT_TREE_SEL_CHANGED, 
                                 self.OnMenuTreeCtrlSelChanged)
        
        # Bind menu items
        for menu_item in self._find_menu_items(xml_document) :
            item_id = wx.xrc.XRCID(menu_item)
            
            function_name = menu_item[:-len("_menu_item")]
            function_name = "On" + underscore_to_camel_case(function_name)
            
            # Bind to handler function
            if hasattr(self, function_name) :
                self.Bind(wx.EVT_MENU, getattr(self, function_name), id=item_id)
            else : 
                menu_item = self.GetMenuBar().FindItemById(item_id)
                menu = menu_item.GetMenu()
                menu.Delete(item_id)
                logging.warning("No function named %s", function_name)
                if menu.GetMenuItemCount() == 0 :
                    index = self.GetMenuBar().FindMenu(menu.GetTitle())
                    if index != wx.NOT_FOUND : 
                        self.GetMenuBar().EnableTop(index, False)
        
        # Fill function menu
        menu_builder.fill_treectrl(left_menu, self._menu_treectrl)
        
        # Disable item that should not be active when no image is loaded
        for item in self._menus_active_when_image_loaded :
            menu_item = self.GetMenuBar().FindItemById(wx.xrc.XRCID(item + "_menu_item"))
            if menu_item is not None :
                menu_item.Enable(False)
        
        # Update UI from preferences
        for attribute in ["cursor_position", "center", "zoom", "display_range"] :
            name = "sync_%s_menu_item"%attribute
            item = self.GetMenuBar().FindItemById(wx.xrc.XRCID(name))
            value = wx.GetApp().get_synchronize_images(attribute)
            item.Check(value)

    ##############
    # Properties #
    ##############
    
    images_panel = property(lambda x:x._images_panel)

    ####################
    # Public interface #
    ####################

    def append_image(self, gui_image) :
        nb_elements = len(self._images_sizer.GetChildren())
        self.insert_image(nb_elements, gui_image)
    
    def insert_image(self, index, gui_image):
        self._images_sizer.Insert(index, gui_image, 1, wx.EXPAND)
        self._adjust_sizer()
        
        mode = gui_image.slice_mode
        item = self.GetMenuBar().FindItemById(wx.xrc.XRCID("view_" + mode +"_menu_item"))
        if item is not None :
            item.Check(True)
        
        if len(wx.GetApp().images) > 0 :
            for view_mode in self._menus_active_when_image_loaded :
                item = self.GetMenuBar().FindItemById(wx.xrc.XRCID(view_mode+"_menu_item"))
                if item is not None :
                    item.Enable(True)
    
    def delete_image(self, image):
        self._images_sizer.Detach(image)
        self._adjust_sizer()
        
        # This is called /before/ the image is deleted
        if len(wx.GetApp().images) == 1 :
            for item in self._menus_active_when_image_loaded :
                menu_item = self.GetMenuBar().FindItemById(wx.xrc.XRCID(item + "_menu_item"))
                if menu_item is not None :
                    menu_item.Enable(False)
    
    def full_screen(self, gui_image):
        
        if gui_image is not None :
            # Create a top-level frame to display the image
            self._full_screen_frame = wx.Frame(None)
            
            command_id = wx.NewId()
            function = lambda event : wx.GetApp().toggle_full_screen()
            self._full_screen_frame.Bind(wx.EVT_MENU, function, id=command_id)

            accelerator_table = wx.AcceleratorTable([
                (0, wx.WXK_ESCAPE, command_id)
            ])
            self._full_screen_frame.SetAcceleratorTable(accelerator_table)

            children = [c.GetWindow() for c in self._images_sizer.GetChildren()]
            self._full_screen_image_index = children.index(gui_image)
            self._full_screen_image = gui_image
            
            # Remove gui_image from sizer
            self._images_sizer.Detach(gui_image)
            
            # Reparent it to the full screen frame
            self._reparent_gui_image(gui_image, self._full_screen_frame)
            
            # Show the frame full screen
            self._full_screen_frame.ShowFullScreen(True)
        elif self._full_screen_frame is not None :
            self._full_screen_frame.Hide()
            
            self._reparent_gui_image(self._full_screen_image, self.images_panel)
            
            self._images_sizer.Insert(self._full_screen_image_index, self._full_screen_image, 1, wx.EXPAND)
            
            # This crude hack forces the Sizer to correctly redraw its children
            size = self.GetSize()
            size[0] -= 1
            self.SetSize(size)
            
            self._full_screen_image_index = None
            self._full_screen_image = None
            
            self._full_screen_frame.Close()
            self._full_screen_frame = None
            
        else : 
            # remove full screen, but we are not in full screen mode : do nothing
            pass
    
    ##############
    # Properties #
    ##############
    
    def _get_current_ui(self):
        return self._current_ui
    
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
            wx.GetApp().append_image(image)
                
    def OnOpenRaw(self, dummy):
        from medipy.gui.image.raw_dialog import RawDialog
        raw_dialog = RawDialog(None, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        raw_dialog.Show()
        
    def OnOpenSpectro(self, dummy):
        from medipy.gui.image.spectro_dialog import SpectroDialog
        spectro_dialog = SpectroDialog(None, size = (500,220))
        spectro_dialog.ShowModal()
    
    def OnFromDirectory(self, dummy):
        
        images = medipy.gui.io.import_dicom_directory(self)
        if images :
            wx.GetApp().append_image(images[0])
    
    def OnSaveImageAs(self, dummy):
        medipy.gui.io.save(wx.GetApp().active_image.get_layer_image(0))
    
    def OnCloseAllImages(self, dummy):
        wx.GetApp().close_all_images()
    
    def OnCineImage(self, dummy) :
        wx.GetApp().set_cine_dialog_visibility(True)
    
    def OnSyncCursorPosition(self, dummy) :
        wx.GetApp().toggle_synchronize_images("cursor_position")
    
    def OnSyncCenter(self, dummy):
        wx.GetApp().toggle_synchronize_images("center")
    
    def OnSyncZoom(self, dummy) :
        wx.GetApp().toggle_synchronize_images("zoom")
    
    def OnSyncDisplayRange(self, dummy):
        wx.GetApp().toggle_synchronize_images("display_range")
    
    def OnMenuTreeCtrlSelChanged(self, evt):
       
        # Destroy existing UI
        if self._current_ui is not None :
            if hasattr(self._current_ui, "Close") :
                self._current_ui.Close() 
            del self._current_ui
            self._function_ui_sizer.Clear(True)
        
        # Create UI for current item data
        item_data = self._menu_treectrl.GetPyData(evt.GetItem())
        if item_data is not None : 
            if inspect.isfunction(item_data) :
                function = item_data
                function_gui=FunctionGUIBuilder(self._function_ui_panel,function,
                                                wx.GetApp().images, wx.GetApp().viewer_3ds,
                                                "function_parameters")
                if wx.GetApp().images :
                    for control in function_gui.controls.values() :
                        if isinstance(control, medipy.gui.control.Coordinates) :
                            control.image = wx.GetApp().active_image
                self._current_ui = function_gui
                self._function_ui_sizer.Add(function_gui.panel, 1, wx.EXPAND)
                self._function_ui_sizer.Layout()
            else :
                ui = item_data(self._function_ui_panel)
                if hasattr(ui, "image") and wx.GetApp().images :
                    ui.image = wx.GetApp().active_image
                self._current_ui = ui
                self._function_ui_sizer.Add(ui, 1, wx.EXPAND)
                self._function_ui_sizer.Layout()
        else : 
            self._current_ui = None
    
    def OnLoadLayer(self, dummy):
        
        images = medipy.gui.io.load(self)
        if images :
            gui_image = wx.GetApp().active_image 
            gui_image.insert_layer(gui_image.number_of_layers, images[0], 
                                   zero_transparency=True)
    
    def OnCloseImage(self, dummy):
        wx.GetApp().close_image(wx.GetApp().active_image)
    
    def OnQuit(self, dummy):
        wx.GetApp().quit()
        self.Close()
    
    def OnLayers(self, dummy):
        wx.GetApp().set_layers_dialog_visibility(True)
    
    def OnAnnotations(self, dummy):
        wx.GetApp().set_annotations_dialog_visibility(True)
    
    def OnViewAxial(self, dummy):
        wx.GetApp().slices = "axial"
    
    def OnViewCoronal(self, dummy):
        wx.GetApp().slices = "coronal"
    
    def OnViewSagittal(self, dummy):
        wx.GetApp().slices = "sagittal"
        
    def OnViewMultiplanar(self, dummy):
        wx.GetApp().slices = "multiplanar"
    
    def OnResetView(self, dummy):
        wx.GetApp().reset_view()
    
    def OnFullScreen(self, dummy):
        wx.GetApp().toggle_full_screen()
    
    def OnCurrentImageScreenshot(self, dummy):
        image = wx.GetApp().current_image_screenshot()
        self._save_image(image)
    
    def OnAllImagesScreenshot(self, dummy):
        image = wx.GetApp().all_images_screenshot()
        self._save_image(image)
    
    def OnWholeWindowScreenshot(self, dummy):
        image = wx.GetApp().whole_window_screenshot()
        self._save_image(image)
    
    def OnExecuteScript(self, dummy):
        path = wx.ConfigBase_Get().Read("LoadImagePath")
        
        file_dialog = wx.FileDialog(self, style=wx.FD_OPEN, defaultDir=path, wildcard="*.py")
        if file_dialog.ShowModal() == wx.ID_OK :
            wx.GetApp().execute_script(file_dialog.GetPath())
        
            wx.ConfigBase_Get().Write("LoadImagePath", file_dialog.GetDirectory())
            wx.ConfigBase_Get().Flush()
    
    def OnNewViewer3d(self, dummy):
        viewer_3d = Viewer3DFrame(None, ObservableList())
        viewer_3d.Show()
        wx.GetApp().append_viewer_3d(viewer_3d)
    
    #####################
    # Private interface #
    #####################
    
    def _adjust_sizer(self):
        nb_objects = len(self._images_sizer.GetChildren())
        rows = max(int(math.ceil(math.sqrt(nb_objects))),1)
        self._images_sizer.SetRows(rows)
        self._images_sizer.SetCols(rows)
        self._images_sizer.Layout()
    
    def _set_slices(self, mode):
        item_id = wx.xrc.XRCID("view_" + mode + "_menu_item")
        item = self.GetMenuBar().FindItemById(item_id)
        item.Check()
    
    slices = property(fset=_set_slices)
    
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
    
    def _find_menu_items(self, xml_document):
        # Find all "object" elements having a "class" property of wxMenuItem
        menu_items = []
        for element in xml_document.getElementsByTagName("object") :
            if element.hasAttribute("class") and \
                    element.getAttribute("class") == "wxMenuItem" and \
                    element.hasAttribute("name") :
                menu_items.append(element.getAttribute("name"))
        
        return menu_items
    
    def _save_image(self, image):
        dialog = wx.FileDialog(self, "Save screenshot", 
            defaultDir = wx.ConfigBase_Get().Read("LoadImagePath"),
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
