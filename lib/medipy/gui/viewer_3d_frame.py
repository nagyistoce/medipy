##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import logging
import os
import re
import sys
import traceback
import xml.dom.minidom

import numpy

from vtk import vtkPlanes, vtkPoints, vtkPolyData, vtkMatrix4x4, vtkTransform

import wx
import wx.xrc

from medipy.base.coordinate_system import slices

import medipy.base
from medipy.base import Object3D, Observable, ObservableList
import medipy.io
import medipy.gui.xrc_wrapper
from medipy.gui import Colormap, get_colormap_from_name, Cine3dDialog
from medipy.gui.control import Image
from medipy.base import find_resource

from medipy.gui.image import ImageLayer
from medipy.gui.utilities import underscore_to_camel_case, remove_access_letter_from_menu

from medipy.visualization_3d.texture import texture_from_depth

import viewer_3d_tools

class Viewer3DFrame(medipy.gui.xrc_wrapper.Frame, Observable):  
    
    def __init__(self, parent, objects_3d, *args, **kwargs):
        
        Observable.__init__(self, ["close"])
        
        self._objects_3d = objects_3d
        self._cine_3d_dialog = None

        self._display_coordinates = None
        if "display_coordinates" in kwargs.keys() and "slicing_mode" in  kwargs.keys():
            self._set_display_coordinates(kwargs["display_coordinates"])
            if kwargs["slicing_mode"] in ["radiological","index"] :
                self.slicing_matrices = slices[kwargs["slicing_mode"]]
            else :
                 raise medipy.base.Exception("Unknown slicing mode : %s"%(kwargs["slicing_mode"],)) 
        else :
            self.slicing_matrices = slices["radiological"]
            self._set_display_coordinates("physical")
        
#        self._interaction_mode = "Normal"
#        self._motion_factor = 0.5
#        self._mouse_wheel_motion_factor = 10
#        # used by the LOD actors
#        self._desired_update_rate = 1
#        self._still_update_rate = 0.0001
#        
#        self._drawings = []
        
        filename = find_resource("resources/gui/viewer_3d_frame.xrc")
        resource = open(filename).read()
        pattern = r"<bitmap>(.*)</bitmap>"
        replacement = r"<bitmap>{0}/\1</bitmap>".format(os.path.dirname(filename))
        resource = re.sub(pattern, replacement, resource)
        xml_resource = xml.dom.minidom.parseString(resource)
        
        resource = wx.xrc.EmptyXmlResource()
        resource.InsertHandler(medipy.gui.xrc_wrapper.Viewer3DXMLHandler())
        resource.InsertHandler(medipy.gui.xrc_wrapper.ScrolledPanelXMLHandler())
        resource.InsertHandler(medipy.gui.xrc_wrapper.ControlImageXMLHandler())
        resource.LoadFromString(str(xml_resource.toxml()))
        
        frame = resource.LoadFrame(parent, "main_frame")
        medipy.gui.xrc_wrapper.Frame.__init__(self, frame, *args, **kwargs)
                
        # window controls
        for control in ["objects_check_list_box", "object_editor_panel", 
                        "inner_display_scrolled_panel",
                        "inner_clipping_planes_scrolled_panel",
                        "viewer_3d",
                        "arrays_choice" 
                       ] :
            setattr(self, "_"+control, wx.xrc.XRCCTRL(self, control))
        
        
        self._outer_display_panel = resource.LoadPanel(self._inner_display_scrolled_panel, "outer_display_panel")
        # display_panel controls
        for control in ["color_type_choice", "color_label", "color_button", "opacity_label",
                        "opacity_slider", "texture_image_label", "texture_image_chooser", 
                        "depth_label", "depth_slider",
                        "representation_combobox", "array_name_combobox", "array_name_label", "points_size_slider", "shading_combobox"
                        ] :
            setattr(self, "_"+control, wx.xrc.XRCCTRL(self._outer_display_panel, control))
        
        
        clipping_planes_panel = resource.LoadPanel(self._inner_clipping_planes_scrolled_panel, "outer_clipping_planes_panel")
        #clipping_planes_panel controls
        for control in ["inside_out_checkbox", "synchronize_checkbox", "axial_slider", "axial_choice", 
                        "coronal_slider", "coronal_choice", "sagittal_slider", "sagittal_choice", 
                        "clipping_image_chooser"
                        ] :
            setattr(self, "_"+control, wx.xrc.XRCCTRL(clipping_planes_panel, control))

        
        sizer = wx.BoxSizer()
        self._inner_display_scrolled_panel.SetSizer(sizer)
        sizer.Add(self._outer_display_panel, 1, wx.EXPAND)
        self._inner_display_scrolled_panel.Layout()
        self._inner_display_scrolled_panel.SetAutoLayout(1)
        self._inner_display_scrolled_panel.SetupScrolling()
        
        sizer = wx.BoxSizer()
        self._inner_clipping_planes_scrolled_panel.SetSizer(sizer)
        sizer.Add(clipping_planes_panel, 1, wx.EXPAND)
        self._inner_clipping_planes_scrolled_panel.Layout()
        self._inner_clipping_planes_scrolled_panel.SetAutoLayout(1)
        self._inner_clipping_planes_scrolled_panel.SetupScrolling()
        
        
        self._viewer_3d.objects_3d = objects_3d
        self._objects_3d.add_observer("any", self._on_objects_3d_modified)
        self._viewer_3d.add_observer("active_object", self._on_active_object)
        self._clipping_image_chooser.may_be_empty = True
        self._clipping_image_chooser.add_observer("any", self._on_clipping_image_changed)
        self._texture_image_chooser.may_be_empty = True
        self._texture_image_chooser.add_observer("any", self._on_texture_image_changed)
        
        wx.GetApp().images.add_observer("any", self._on_medipy_images_list_changed)


        menu_items = []
        tools = []
        for object in xml_resource.getElementsByTagName("object") :
            if object.hasAttribute("class") and \
                    object.getAttribute("class") == "wxMenuItem" and \
                    object.hasAttribute("name") :
                menu_items.append(object.getAttribute("name"))
            if object.hasAttribute("class") and \
                    object.getAttribute("class") == "tool" and \
                    object.hasAttribute("name") :
                tools.append(object.getAttribute("name"))
        
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        self._objects_check_list_box.Bind(wx.EVT_CHAR, self.OnObjectChar)
        self._objects_check_list_box.Bind(wx.EVT_LISTBOX_DCLICK, self.OnObjectDoubleClick)
        self._objects_check_list_box.Bind(wx.EVT_LISTBOX, self.OnObjectSelected)
        self._objects_check_list_box.Bind(wx.EVT_CHECKLISTBOX, self.OnObjectChecked)
        self._color_type_choice.Bind(wx.EVT_CHOICE, self.OnColorTypeChanged)
        self._color_button.Bind(wx.EVT_BUTTON, self.OnColorChanged)
        self._opacity_slider.Bind(wx.EVT_SLIDER, self.OnOpacityChanged)
        self._representation_combobox.Bind(wx.EVT_COMBOBOX, self.OnRepresentationChanged)
        self._array_name_combobox.Bind(wx.EVT_COMBOBOX, self.OnArrayNameChanged)
        self._points_size_slider.Bind(wx.EVT_SLIDER, self.OnPointsSizeChanged)
        self._shading_combobox.Bind(wx.EVT_COMBOBOX, self.OnShadingChanged)
        self._depth_slider.Bind(wx.EVT_SLIDER, self.OnDepthChanged)
        self._inside_out_checkbox.Bind(wx.EVT_CHECKBOX, self.OnInsideOutChanged)
        self._synchronize_checkbox.Bind(wx.EVT_CHECKBOX, self.OnSynchronizeClippingPlanes)
        
        for axis in ("axial", "coronal", "sagittal") :
            for widget in ("slider", "choice") :
                strg = axis + "_" + widget
                window_widget = getattr(self, "_" + strg)
                function_name = "On" + underscore_to_camel_case(strg) + "Changed"
                evt = "EVT_" + widget.upper()
                evt = getattr(wx, evt)
                window_widget.Bind(evt, getattr(self, function_name))
                
        for menu_item in menu_items :
            id = wx.xrc.XRCID(menu_item)
            item = self.GetMenuBar().FindItemById(id)
            label = item.GetLabel()
            
            item.SetText(label)
            
            function_name = menu_item[:-len("_menu_item")]
            function_name = "On" + underscore_to_camel_case(function_name)
            
            # Bind to handler function
            if hasattr(self, function_name) :
                self.Bind(wx.EVT_MENU, getattr(self, function_name), id=id)
            else : 
                menu_item = self.GetMenuBar().FindItemById(id)
                menu = menu_item.GetMenu()
                menu.Delete(id)
                logging.warning("No function called %s"%function_name)
                
            # In order to disable them, when it is necessary :
            if menu_item == "animation_menu_item" :
                self._animation_menu_item = self.GetMenuBar().FindItemById(id)
            if menu_item == "save_animation_menu_item" :
                self._save_animation_menu_item = self.GetMenuBar().FindItemById(id)
            # To operate on distance_menu_item ... (verify check)
            if menu_item == "distance_menu_item" :
                self._distance_menu_item = self.GetMenuBar().FindItemById(id)
                
        for tool in tools : 
            id = wx.xrc.XRCID(tool)
            function_name = "On" + underscore_to_camel_case(tool)
            self.Bind(wx.EVT_TOOL, getattr(self, function_name), id=id)
            
        self._enable()
        self._axial_slider.Enable(False)
        self._coronal_slider.Enable(False)
        self._sagittal_slider.Enable(False)
        self._inside_out_checkbox.Enable(False)
        if self._representation_combobox.GetValue() != "Points" :
            self._points_size_slider.Enable(False)
        else :
            self._points_size_slider.Enable(True)
        
        self._texture_image_label.Hide()
        self._texture_image_chooser.Hide()
        self._depth_label.Hide()
        self._depth_slider.Hide()
        
        self._clipping_image_chooser.choices = wx.GetApp().images
        self._texture_image_chooser.choices = wx.GetApp().images
        
        # self._memento save the state and properties of all of the 3d objects 
        # in order to reflect it on the editor panel and in the check_list_box
        self._memento = {}
    
        # initialization of the viewer_3d_frame tools
        self._viewer_3d.tools["LeftButtonShift"] = viewer_3d_tools.Rotate(self._viewer_3d)

    def _set_display_coordinates(self, display_coordinates) :
        if display_coordinates not in ["physical", "nearest_axis_aligned", "index"] :
            raise medipy.base.Exception("Unknown display coordinates : %s"%(display_coordinates,))
        
        self._display_coordinates = display_coordinates      
    
    
    def view_all(self):
        self._viewer_3d.view_all()
    
    def update_object_editor(self):
        index = self._objects_check_list_box.GetSelection()
        if index == wx.NOT_FOUND :
            return
        
        object = self._objects_3d[index]
        
#        self._arrays_choice.Clear()
#        self._arrays_choice.Append("Uniform material")
#        self._arrays_choice.Append("--- Point arrays ---")
#        for array in object.point_arrays() :
#            self._arrays_choice.Append(array)
#        self._arrays_choice.Append("--- Cell arrays ---")
#        for array in object.cell_arrays() :
#            self._arrays_choice.Append(array)
#        self._arrays_choice.Append("--- Field arrays ---")
#        for array in object.field_arrays() :
#            self._arrays_choice.Append(array)
         
        color = [int(255.*c) for c in object.color]
        self._color_button.SetBackgroundColour(color)
        
        opacity = int(object.opacity*100)
        self._opacity_slider.SetValue(opacity)
        
        representation = object.representation
        self._representation_combobox.SetValue(representation)
        
        points_size = object.points_size
        self._points_size_slider.SetValue(points_size)
        
        shading = object.shading
        self._shading_combobox.SetValue(shading)
        
        state = self._memento[object]
        coloration_type = state["coloration"]
        self._color_type_choice.Select(coloration_type)
        if coloration_type == 0 : #"Color"
            self._color_label.Show()
            self._color_button.Show()
            self._opacity_label.Show()
            self._opacity_slider.Show()
            self._texture_image_label.Hide()
            self._texture_image_chooser.Hide()
            self._depth_label.Hide()
            self._depth_slider.Hide()
            self._array_name_label.Hide()
            self._array_name_combobox.Hide()
        elif coloration_type == 1 : #"Pseudo-texture"
            self._color_label.Hide()
            self._color_button.Hide()
            self._opacity_label.Hide()
            self._opacity_slider.Hide()
            self._texture_image_label.Show()
            self._texture_image_chooser.Show()
            self._depth_label.Show()
            self._depth_slider.Show()
            self._array_name_label.Hide()
            self._array_name_combobox.Hide()
        elif coloration_type == 2 : #"Point array"
            self._color_label.Hide()
            self._color_button.Hide()
            self._opacity_label.Hide()
            self._opacity_slider.Hide()
            self._texture_image_label.Hide()
            self._texture_image_chooser.Hide()
            self._depth_label.Hide()
            self._depth_slider.Hide()
            self._array_name_label.Show()
            self._array_name_combobox.Show()
        
        
        self._texture_image_chooser.remove_observer("any", self._on_texture_image_changed)
        if state["texture_image"] :
            self._texture_image_chooser.may_be_empty_checked = False
        else :
            self._texture_image_chooser.may_be_empty_checked = True
        self._texture_image_chooser.value = state["texture_image"]
        self._texture_image_chooser.add_observer("any", self._on_texture_image_changed)
        
        
        self._depth_slider.SetValue(state["depth"])
        
        self._update_clipping_panel(object)
        
        
    def _update_clipping_panel(self, object):
        state = self._memento[object]
        
        self._inside_out_checkbox.SetValue(state["check_box"])
           
        self._axial_slider.SetMin(state["axial_slider_min"])
        self._axial_slider.SetMax(state["axial_slider_max"])
        self._axial_slider.SetValue(state["axial_slider_value"])
        self._axial_choice.Select(state["axial_choice"])
            
        self._coronal_slider.SetMin(state["coronal_slider_min"])
        self._coronal_slider.SetMax(state["coronal_slider_max"])
        self._coronal_slider.SetValue(state["coronal_slider_value"])
        self._coronal_choice.Select(state["coronal_choice"])
            
        self._sagittal_slider.SetMin(state["sagittal_slider_min"])
        self._sagittal_slider.SetMax(state["sagittal_slider_max"])
        self._sagittal_slider.SetValue(state["sagittal_slider_value"])
        self._sagittal_choice.Select(state["sagittal_choice"])
        
        self._manage_clipping_sliders(object)
        
        self._manage_images(state, object)
         
    
    def select_object(self, n):
        """ Select the nth object
        """
        self._objects_check_list_box.SetSelection(n)
        self.update_object_editor()
    
    def delete_object(self, n):
        """ Delete the nth object
        """
        if n in range(len(self._objects_3d)) :
            del self._objects_3d[n]
            self._update_objects_check_list_box()
            if len(self._objects_3d) > 0 :
                index = min(n, len(self._objects_3d)-1)
                self._objects_check_list_box.SetSelection(index)
            self.update_object_editor()
        self._viewer_3d.render_window_interactor.Render()
            
    def show_cine_3d_dialog(self):
        # Create only one dialog, raise it if it is hidden
        if self._cine_3d_dialog is not None : 
            self._cine_3d_dialog.Raise()
        else :
            self._cine_3d_dialog = Cine3dDialog(self)
            self._cine_3d_dialog.Show()
            self._cine_3d_dialog.Bind(wx.EVT_CLOSE, self.OnCloseCine3dDialog)
    
        
    ##############
    # Properties #
    ##############
    
    def _get_memento(self):
        return dict(self._memento)
    
    objects_3d = property(lambda x:x._objects_3d)
    memento = property(_get_memento)
    
    
    ##################
    # Event handlers #
    ##################
    
    def OnAnimation(self, event):
        self.show_cine_3d_dialog()
    
    def OnArrayChanged(self, event):
        selection = self._arrays_choice.GetStringSelection()
        if selection in [ "--- %s arrays ---"%(x) for x in ["Point", "Cell", "Field"] ] :
            return
        
        #point_start = self._arrays_choice.FindString("--- Point arrays ---")
        cell_start = self._arrays_choice.FindString("--- Cell arrays ---")
        field_start = self._arrays_choice.FindString("--- Field arrays ---")
        
        object = self._objects_3d[self._objects_check_list_box.GetSelection()]
        if selection == "Uniform material" :
            object.color_by_material()
        else :
            if self._arrays_choice.GetSelection() > field_start : 
                object.color_by_field_scalars(selection)
                array = object.get_field_array(selection)
            elif self._arrays_choice.GetSelection() > cell_start :
                object.color_by_cell_scalars(selection)
                array = object.get_cell_array(selection)
            else : 
                object.color_by_point_scalars(selection)
                array = object.get_point_array(selection)
            object.set_scalar_range(array.GetRange())
        self._viewer_3d.render_window_interactor.Render()
    
    def OnClose(self, event):
        if self._cine_3d_dialog is not None : 
            self._cine_3d_dialog.OnClose(event)
            #self._cine_3d_dialog.Destroy()
        
        # Delete all objects to avoid crash on close
        while self._objects_3d :
            self.delete_object(0)
        
        wx.CallAfter(lambda :self.notify_observers("close"))
        self.Destroy()
    
    def OnColorTypeChanged(self, event):
        index = self._objects_check_list_box.GetSelection()
        _object = self._objects_3d[index]
        state = self._memento[_object]
        image = state["texture_image"]
        
        choice = self._color_type_choice.GetStringSelection()
        if choice == "Color" :
            self._color_label.Show()
            self._color_button.Show()
            self._opacity_label.Show()
            self._opacity_slider.Show()
            self._texture_image_label.Hide()
            self._texture_image_chooser.Hide()
            self._depth_label.Hide()
            self._depth_slider.Hide()
            self._array_name_label.Hide()
            self._array_name_combobox.Hide()
            
            _object.color_by_material()
        elif choice == "Pseudo-texture" :
            self._color_label.Hide()
            self._color_button.Hide()
            self._opacity_label.Hide()
            self._opacity_slider.Hide()
            self._texture_image_label.Show()
            self._texture_image_chooser.Show()
            self._depth_label.Show()
            self._depth_slider.Show()
            self._array_name_label.Hide()
            self._array_name_combobox.Hide()
            
            if image :
                depth = self._depth_slider.GetValue()/100.
                self._manage_texture_from_depth(_object, image, depth)
                self._depth_slider.Enable()
            else :
                self._depth_slider.Disable()

        elif choice == "Point array" :
            self._color_label.Hide()
            self._color_button.Hide()
            self._opacity_label.Hide()
            self._opacity_slider.Hide()
            self._texture_image_label.Hide()
            self._texture_image_chooser.Hide()
            self._depth_label.Hide()
            self._depth_slider.Hide()
            self._array_name_label.Show()
            self._array_name_combobox.Show()

            self._array_name_combobox.Clear()
            array_names = _object.point_arrays()
            for name in array_names :
                self._array_name_combobox.Append(str(name)) 
                self._array_name_combobox.SetValue(str(name))
                _object.color_by_point_scalars(str(name))         
     
        self._save_object_state(_object)
        
        #self._outer_display_scrolled_panel.Layout()
        
        self._viewer_3d.render_window_interactor.Render()
    
    
    def OnColorChanged(self, event):
        self._update_color()
    
    def OnObjectSelected(self, event):
        self.update_object_editor()
    
    def OnObjectChecked(self, event):
        # The EVT_CHECKLISTBOX is fired before the selection is updated
           
        index = event.GetSelection()
        _object = self._objects_3d[index]
        state = self._memento[_object]
        
        _object.visibility = self._objects_check_list_box.IsChecked(index) 
        state["static_visibility"] = _object.visibility
            
        # hide/unhide the image_layers linked with the checked object
        for axis in ("axial", "coronal", "sagittal") :
            choice = getattr(self, "_"+axis+"_choice")
            choice = choice.GetStringSelection()
            
            if choice != "Inactive" :
                state[axis + "_image_layer"].actor.SetVisibility(_object.visibility)

        self._viewer_3d.render_window_interactor.Render()
        
    def OnObjectDoubleClick(self, event):
        index = self._objects_check_list_box.GetSelection()
        was_checked = self._objects_check_list_box.IsChecked(index)
        text = self._objects_check_list_box.GetString(index)
        new_name = wx.GetTextFromUser('New name for object %s'%text, 'Rename object', text)
        if new_name != '':dataset = vtkPolyData()
        points = vtkPoints()
        dataset.SetPoints(points)
        object = Object3D(dataset, "New drawing")
        object.show_vertices()
        object.actor.GetProperty().SetPointSize(5)
        self._objects_3d.append(object)
        self._objects_check_list_box.SetString(index, new_name)
        self._objects_3d[index].name = new_name
        # Setting the string seems to un-check the checkbox
        self._objects_check_list_box.Check(index, was_checked)
    
    def OnOpacityChanged(self, event):
        index = self._objects_check_list_box.GetSelection()
        object = self._objects_3d[index]
        
        opacity = self._opacity_slider.GetValue()/100.
        object.opacity = opacity
        self._viewer_3d.render_window_interactor.Render()
        
    
    def OnOpen(self, event):
        path = wx.ConfigBase_Get().Read("LoadObject3DPath")
        
        dialog = wx.FileDialog(self, 
                               style = wx.FD_DEFAULT_STYLE | wx.FD_MULTIPLE,
                               defaultDir=path)
        if dialog.ShowModal() == wx.ID_OK :
            wx.ConfigBase_Get().Write("LoadObject3DPath", dialog.GetDirectory())
            wx.ConfigBase_Get().Flush()
            
            for path in dialog.GetPaths() :
                try :
                    wx.GetApp().load_object_3d(path, self)
                    self.select_object(len(self._objects_3d)-1)
                except Exception, e :
                    exc_info = sys.exc_info()
                    logging.error("".join(traceback.format_exception(*exc_info)))
                    wx.MessageBox(str(e), "Could not load %s"%path)
    
    def OnSave(self, event):
        path = wx.ConfigBase_Get().Read("LoadObject3DPath")
        
        dialog = wx.FileDialog(self, 
                               style = wx.FD_DEFAULT_STYLE | wx.FD_SAVE,
                               defaultDir=path)
        if dialog.ShowModal() == wx.ID_OK :
            wx.ConfigBase_Get().Write("LoadObject3DPath", dialog.GetDirectory())
            wx.ConfigBase_Get().Flush()
            try : 
                index = self._objects_check_list_box.GetSelection()
                object = self._objects_3d[index]
                wx.GetApp().save_object_3d(object, dialog.GetPath())
            except Exception, e :
                    wx.MessageBox(str(e), "Could not save %s"%path)
    
    def OnViewAll(self, event):
        self._viewer_3d.view_all()
    
    
    def OnRotateTool(self, event):
        self._setup_tool(viewer_3d_tools.Rotate(self._viewer_3d))
    
    
    def OnPickTool(self, event):
        self._setup_tool(viewer_3d_tools.Pick(self._viewer_3d))
        
    
    def OnDistanceTool(self, event):
        self._setup_tool(viewer_3d_tools.Distance(self._viewer_3d))
    
    
    def OnObjectChar(self, event):
        if event.GetKeyCode() == wx.WXK_DELETE :
            index = self._objects_check_list_box.GetSelection()
            self.delete_object(index)
    
    
    def OnDelete(self, event):
        index = self._objects_check_list_box.GetSelection()
        self.delete_object(index)
    
    
    def OnNewDrawing(self, event):
        dataset = vtkPolyData()
        points = vtkPoints()
        dataset.SetPoints(points)
        object = Object3D(dataset, "New drawing")
        object.show_vertices()
        object.actor.GetProperty().SetPointSize(5)
        self._objects_3d.append(object)
        self._drawings.append(object)
            
        
    def OnRepresentationChanged(self, event):
        index = self._objects_check_list_box.GetSelection()
        repr = self._representation_combobox.GetValue()
        self.objects_3d[index].representation = repr
        if repr != "Surface" :
            self.objects_3d[index].actor.GetProperty().SetAmbient(1.)
        else :
            self.objects_3d[index].actor.GetProperty().SetAmbient(0.)
            
        if repr != "Points" :
            self._points_size_slider.Enable(False)
        else :
            self._points_size_slider.Enable(True)
        self._viewer_3d.render_window_interactor.Render()

    def OnArrayNameChanged(self, event):
        index = self._objects_check_list_box.GetSelection()
        _object = self._objects_3d[index]
        index = self._objects_check_list_box.GetSelection()
        repr = self._array_name_combobox.GetValue()
        _object.color_by_point_scalars(repr)
 
        self._save_object_state(_object)       
        self._viewer_3d.render_window_interactor.Render()
          
    def OnPointsSizeChanged(self, event):
        size = self._points_size_slider.GetValue()
        index = self._objects_check_list_box.GetSelection()
        self.objects_3d[index].points_size = size
        self._viewer_3d.render_window_interactor.Render()
        

    def OnShadingChanged(self, event):
        index = self._objects_check_list_box.GetSelection()
        shad = self._shading_combobox.GetValue()
        self.objects_3d[index].shading = shad
        self._viewer_3d.render_window_interactor.Render()
        
        
    def OnDepthChanged(self, event):    
        index = self._objects_check_list_box.GetSelection()
        _object = self._objects_3d[index]
        state = self._memento[_object]
        image = state["texture_image"]
        depth = self._depth_slider.GetValue()/100.
        self._manage_texture_from_depth(_object, image, depth)
        self._save_object_state(_object)
        
        self._viewer_3d.render_window_interactor.Render()
        
        
    def OnCloseCine3dDialog(self, event):
        self._cine_3d_dialog.OnClose(event)
        self._cine_3d_dialog.Destroy()
        self._cine_3d_dialog = None
        event.Skip()
        
    def OnSaveAnimation(self, event):
        print "Save animation : " + str(event)
        
        
    def OnInsideOutChanged(self, event):
        index = self._objects_check_list_box.GetSelection()
        object = self._objects_3d[index]
        object.inside_out = not self._inside_out_checkbox.GetValue()
        self._save_object_state(object)
        self._viewer_3d.render_window_interactor.Render()
        
    def OnAxialChoiceChanged(self, event):
        self._synchronized_cut("choice", "axial")
        self._viewer_3d.render_window_interactor.Render()
    
    def OnAxialSliderChanged(self, event):
        self._synchronized_cut("slider", "axial")
        self._viewer_3d.render_window_interactor.Render()
        
    def OnCoronalChoiceChanged(self, event):
        self._synchronized_cut("choice", "coronal")
        self._viewer_3d.render_window_interactor.Render()
        
    def OnCoronalSliderChanged(self, event):
        self._synchronized_cut("slider", "coronal")
        self._viewer_3d.render_window_interactor.Render()
        
    def OnSagittalChoiceChanged(self, event):
        self._synchronized_cut("choice", "sagittal")
        self._viewer_3d.render_window_interactor.Render()
        
    def OnSagittalSliderChanged(self, event):
        self._synchronized_cut("slider", "sagittal")
        self._viewer_3d.render_window_interactor.Render()
        
    def OnSynchronizeClippingPlanes(self, event):
        self._viewer_3d.render_window_interactor.Render()
        
        
        
    def _on_clipping_image_changed(self, event):
        index = self._objects_check_list_box.GetSelection()
        _object = self._objects_3d[index]
        state = self._memento[_object]
        
        _object.image = self._clipping_image_chooser.value
        if _object.image :
            (min, max) = _object.image.data.min(), _object.image.data.max()
            for axis in ("axial", "coronal", "sagittal") :
                self._display_image_layers(axis, _object)
                state[axis + "_image_layer"].image = self._clipping_image_chooser.value
                state[axis + "_image_layer"].display_range = (min, max)
                self._move_image_layers(axis, _object)
        else :
            for axis in ("axial", "coronal", "sagittal") :
                if state.has_key(axis+"_image_layer") :
                    state[axis + "_image_layer"].actor.SetVisibility(False)
        
        self._viewer_3d.render_window_interactor.Render()
        
        
    def _on_texture_image_changed(self, event):
        index = self._objects_check_list_box.GetSelection()
        _object = self._objects_3d[index]
        state = self._memento[_object]
        
        image = self._texture_image_chooser.value
        state["texture_image"] = image 
        
        if image :
            depth = self._depth_slider.GetValue()/100.
            self._manage_texture_from_depth(_object, image, depth)
            self._depth_slider.Enable()
        else :
            self._depth_slider.Disable()
            
        self._save_object_state(_object)
        
        
        self._viewer_3d.render_window_interactor.Render()
        
        
    def _on_medipy_images_list_changed(self, event):
        pass
        #in order to .Layout() the display panel ...
        
        
        
        
    #####################
    # Private interface #
    #####################
    
    
    def _clean_memento(self, liste):
        new_memento = {}
        for _object in liste :
            new_memento[_object] = self._memento[_object]
        self._memento = new_memento
        
        
    def _create_image_layers(self, state, _object) :
        image = _object.image
        colormap_name = image.metadata["colormap"] if "colormap" in image.metadata else "gray"
        colormap = Colormap(get_colormap_from_name(colormap_name), None)
        for axis, matrice in self.slicing_matrices.items() :
            image_layer = ImageLayer(matrice, image, colormap=colormap, display_coordinates=self._display_coordinates)
            image_layer.position = numpy.divide(_object.image.shape, 2.)
            matrix = vtkMatrix4x4()
            matrix.Identity()
            direction = image_layer._slice_to_world[::-1,::-1]
            for i in range(3):
                for j in range(3):
                    matrix.SetElement(i,j,direction[i,j])
            image_layer.actor.SetUserMatrix(matrix)
            image_layer.zero_transparency = True


            image_layer.actor.SetVisibility(False)
            
            center = list(_object.image.shape)
            center.reverse()
            center = numpy.divide(center, 2.)            

            image_layer.actor.SetPosition(center)
            
            self._viewer_3d.renderer.AddActor(image_layer.actor)
            
            state[axis + "_image_layer"] = image_layer
            state[axis + "_image_layer_center"] = center
        
        
    def _dispatch_interaction(self,dataset = vtkPolyData()):
        points = vtkPoints()
        dataset.SetPoints(points)
        object = Object3D(dataset, "New drawing")
        object.show_vertices()
        object.actor.GetProperty().SetPointSize(5)
        self._objects_3d.append(object)
        if self._button is not None :
            self._tools[self._button].dispatch_interaction()
            
    
    def _display_image_layers(self, axis, _object) :
        state = self._memento[_object]
        if not state.has_key(axis+"_image_layer") :
            self._create_image_layers(state, _object)
                        
        choice = getattr(self, "_"+axis+"_choice")
        choice = choice.GetStringSelection()
        
        if choice == "Inactive" :
            state[axis+"_image_layer"].actor.SetVisibility(False)
        else :
            state[axis+"_image_layer"].actor.SetVisibility(True)
            

    def _enable(self):
        # Disable the 3D Cine menu, because no object is loaded at the opening
        # or Enable it when objects (3d) are opened
        # Do the same for the object_editor_panel
        if len(self.objects_3d) > 0 :
            self._object_editor_panel.Enable(True)
            self._animation_menu_item.Enable(True)
            self._save_animation_menu_item.Enable(True)
        else :
            self._object_editor_panel.Enable(False)
            self._animation_menu_item.Enable(False)
            self._save_animation_menu_item.Enable(False)
        
        
    def _synchronized_cut(self, widget, axis):
        if self._synchronize_checkbox.IsChecked() : #ToCheck
            for _object in self._objects_3d :
                self._extended_cut(_object)
                if _object.image :
                    if widget == "choice" :
                        self._display_image_layers(axis, _object)
                        self._move_image_layers(axis, _object)
                    elif widget == "slider" :
                        self._move_image_layers(axis, _object)
        else :
            index = self._objects_check_list_box.GetSelection()
            _object = self._objects_3d[index]
            self._extended_cut(_object, axis)
            if _object.image :

                #_object.actor.SetPosition(_object.image.origin[::-1])
                #_object.actor.SetScale(_object.image.spacing[::-1])

                if widget == "choice" :
                    self._display_image_layers(axis, _object)
                    self._extended_cut(_object, axis)
                    self._update_image_layers(axis, _object)
                elif widget == "slider" :
                    self._update_image_layers(axis, _object)    
    
    def _extended_cut(self, _object, axis):
        if _object.image :
            shape_inf = [0,0,0] #_object.image.origin
            shape_sup = numpy.asarray(_object.image.shape)-1

            if self._display_coordinates in ["physical","nearest_axis_aligned"] :
                shape_sup = _object.image.index_to_physical(shape_sup)
                shape_inf = _object.image.index_to_physical(shape_inf)

            zmax, ymax, xmax = shape_sup
            zmin, ymin, xmin = shape_inf

        else :
            raise medipy.base.Exception("Cannot cut volume in physical coordinate system without an image attached to the Oject3D.")
            #xmin, xmax, ymin, ymax, zmin, zmax = _object.dataset.GetBounds() 
        
       
        bounds = [xmin, xmax, ymin, ymax, zmin, zmax]
        self._manage_clipping(bounds, _object)
        
        self._internal_cut(_object, bounds)
        
        self._save_object_state(_object)
        
        
    def _internal_cut(self, _object, bounds):
        if _object.clipping_functions :

            state = self._memento[_object]
            activated_acs = []
            for axis in ["axial","coronal","sagittal"] :
                if state[axis+"_image_layer"].actor.GetVisibility() :
                    activated_acs.append( axis )
            
            point = numpy.asarray(_object.image.shape)-1
            cut_options = []
            for axis in activated_acs :
                if axis == "axial" :
                    cut_options.append( self._axial_choice.GetStringSelection() )
                    axi_val = self._axial_slider.GetValue()
                    point[1] = axi_val
                elif axis == "coronal" :
                    cut_options.append( self._coronal_choice.GetStringSelection() )
                    cor_val = self._coronal_slider.GetValue() 
                    point[2] = cor_val
                elif axis == "sagittal" :
                    cut_options.append( self._sagittal_choice.GetStringSelection() )
                    sag_val = self._sagittal_slider.GetValue()
                    point[0] = sag_val
                else :
                    raise medipy.base.Exception("Uknown axis : %s"%(axis,))

            if self._display_coordinates in ["physical","nearest_axis_aligned"] :
                if _object.image :
                    point = _object.image.index_to_physical(point)
                else :
                    raise medipy.base.Exception("Cannot cut volume in physical coordinate system without an image attached to the Oject3D.")  
            point = point[::-1]
  
            b_inf = bounds[::2] 
            b_sup = bounds[1::2]
            key_inf = list(point != b_inf)
            key_sup = list(point != b_sup)
            key = key_inf and key_sup     

            bounds_sorted = [] # Space Problem ?
            axis = ["sagittal", "coronal", "axial"]
            new_axis = [] # signed ?
            for i in range(3) :
                if b_inf[i]<b_sup[i] :
                    bounds_sorted.append(b_inf[i])
                    bounds_sorted.append(b_sup[i])
                else :
                    bounds_sorted.append(b_sup[i])
                    bounds_sorted.append(b_inf[i])
 
                perm = numpy.where(_object.image._index_to_physical_matrix[i]!=0)[0][0]
                if axis[perm] in activated_acs :
                    new_axis.append(axis[perm])              
            bounds = numpy.asarray(bounds_sorted)
                
            while True in key :
                index = key.index(True)
                cut_option = cut_options[activated_acs.index(axis[index])]
                if cut_option in ["Keep superior", "Keep posterior", "Keep right"] : 
                    bounds[index*2] = point[index]
                else :
                    bounds[index*2+1] = point[index]
                key[index] = False
            
            _object.clipping_functions[0].SetBounds(bounds)
         
    
    def _manage_clipping(self, bounds, object):   
        if not (object.clipping_functions) :
            p = vtkPlanes()
            p.SetBounds(bounds)
            object.clipping_functions = [p]
            
        self._manage_clipping_sliders(object)
            

    def _manage_clipping_sliders(self, object) :
        axi = self._axial_choice.GetStringSelection()
        cor = self._coronal_choice.GetStringSelection()
        sag = self._sagittal_choice.GetStringSelection()
        
        axi_cond = (axi == "Inactive")
        cor_cond = (cor == "Inactive")
        sag_cond = (sag == "Inactive")
        
        if axi_cond :
            self._axial_slider.Enable(False)
        else :
            self._axial_slider.Enable(True)
        
        if cor_cond :
            self._coronal_slider.Enable(False)
        else :
            self._coronal_slider.Enable(True)
        
        if sag_cond :
            self._sagittal_slider.Enable(False)
        else :
            self._sagittal_slider.Enable(True)    
            
        if axi_cond and cor_cond and sag_cond :
            object.clipping_functions = []
            self._inside_out_checkbox.Enable(False)
            self._viewer_3d.render_window_interactor.Render()
        else :
            self._inside_out_checkbox.Enable(True)
            self._viewer_3d.render_window_interactor.Render()
            
    
    def _manage_images(self, state, _object):
        self._clipping_image_chooser.remove_observer("any", self._on_clipping_image_changed)
        if _object.image :
            self._clipping_image_chooser.may_be_empty_checked = False
        else :
            self._clipping_image_chooser.may_be_empty_checked = True
        self._clipping_image_chooser.value = _object.image
        self._clipping_image_chooser.add_observer("any", self._on_clipping_image_changed)
        
    
    def _manage_texture_from_depth(self, _object, image, depth):
        (min, max) = image.data.min(), image.data.max() 
        lut = medipy.vtk.build_vtk_colormap(medipy.gui.colormaps["gray"])
        _object.color_by_scalars(lut, min, max)
        texture_from_depth(_object, image, depth)
            
    def _update_image_layers(self, axis, _object):
        indexes  = {"axial" : 1, "coronal" : 2, "sagittal" : 0}
        i = indexes[axis]
        
        state = self._memento[_object]
        
        value = getattr(self, "_"+axis+"_slider")
        value = value.GetValue()          

        state[axis + "_image_layer"].position[i] = value
        new_position = state[axis + "_image_layer"].position
            
        if self._display_coordinates == "index" :
            state[axis + "_image_layer"]._set_index_position(new_position)
            direction = state[axis + "_image_layer"]._world_to_slice[::-1,::-1]
            state[axis + "_image_layer"].actor.SetPosition( numpy.dot(direction,new_position[::-1]) )
        else :
            state[axis + "_image_layer"]._set_physical_position(_object.image.index_to_physical(new_position))
            direction = state[axis + "_image_layer"]._world_to_slice[::-1,::-1]
            state[axis + "_image_layer"].actor.SetPosition( numpy.dot(direction,_object.image.index_to_physical(new_position)[::-1]) )
        state[axis + "_image_layer"]._update_change_information()
            

        
        
    def _setup_tool(self, tool) :
        self._viewer_3d.tools["LeftButton"].deselect()
        self._viewer_3d.tools["LeftButton"] = tool
        self._viewer_3d.tools["LeftButton"].select()
    
    
    def _on_active_object(self, event):
        object = event.object.active_object
        index = self._objects_3d.index(object)
        self.select_object(index)
        
    
    def _on_objects_3d_modified(self, event):
        self._enable()
        
        if event.event == "append" :
            self._update_memento(event.object)
        
        if event.event == "delete_item" :
            if event.old_value.image :
                state = self._memento[event.old_value]
                for axis in ("axial", "coronal", "sagittal") :
                    actor = state[axis + "_image_layer"].actor
                    self._viewer_3d.renderer.RemoveActor(actor)
            self._clean_memento(event.object)
            
        self._update_objects_check_list_box()
        
        was_empty = ((event.event == "append") or (event.event == "insert")) and (len(event.object) == 1) 
        if (was_empty) :
            self._objects_check_list_box.SetSelection(0)
            self.update_object_editor()
            Pmin, Pmax = self._viewer_3d.get_bounds()
            camera_position = ((Pmin[0]+Pmax[0])/2., 
                               Pmax[1]*250,
                               (Pmin[2]+Pmax[2])/2.,)  
            camera = self._viewer_3d.renderer.GetActiveCamera()
            camera.SetPosition(camera_position)
            camera.SetViewUp(0, 0, 1)
            self._viewer_3d.view_all()
        self._viewer_3d.render_window_interactor.Render()
        
    
    
    def _save_object_default_state(self, state, _object):
        if _object.image :
            zmin, ymin, xmin = [0,0,0] #_object.image.origin
            zmax, ymax, xmax = _object.image.shape
        else :
            _object.dataset.Update()
            xmin, xmax, ymin, ymax, zmin, zmax = _object.dataset.GetBounds()
            
        Pmin = (xmin, ymin, zmin)
        Pmax = (xmax, ymax, zmax)
        
        center = numpy.add(Pmin, Pmax)/2.
            
        state["check_box"] = False
        state["axial_slider_min"] = zmin
        state["axial_slider_max"] = ymax-1
        state["axial_slider_value"] = center[2]
        state["axial_choice"] = 0 #Inactive
        state["coronal_slider_min"] = ymin    
        state["coronal_slider_max"] = xmax-1
        state["coronal_slider_value"] = center[1]
        state["coronal_choice"] = 0 #Inactive
        state["sagittal_slider_min"] = xmin    
        state["sagittal_slider_max"] = zmax-1
        state["sagittal_slider_value"] = center[0]
        state["sagittal_choice"] = 0 #Inactive
        state["static_visibility"] = _object.visibility
        state["dynamic_visibility"] = False
        state["coloration"] = 0 #Color
        state["texture_image"] = None
        state["depth"] = 90
        
        if _object.image :
            self._create_image_layers(state, _object)
            
            
    
    def _save_object_state(self, _object, default=False):
        if (default) :
            state = {}
            self._save_object_default_state(state, _object)  
        else :
            state = self._memento[_object]
            flag = self._inside_out_checkbox.GetValue()
            state["check_box"] = flag
            state["axial_slider_min"] = self._axial_slider.GetMin()        
            state["axial_slider_max"] = self._axial_slider.GetMax()
            state["axial_slider_value"] = self._axial_slider.GetValue()
            state["axial_choice"] = self._axial_choice.GetSelection()
            state["coronal_slider_min"] = self._coronal_slider.GetMin()        
            state["coronal_slider_max"] = self._coronal_slider.GetMax()
            state["coronal_slider_value"] = self._coronal_slider.GetValue()
            state["coronal_choice"] = self._coronal_choice.GetSelection()
            state["sagittal_slider_min"] = self._sagittal_slider.GetMin()        
            state["sagittal_slider_max"] = self._sagittal_slider.GetMax()
            state["sagittal_slider_value"] = self._sagittal_slider.GetValue()
            state["sagittal_choice"] = self._sagittal_choice.GetSelection()
            state["coloration"] = self._color_type_choice.GetSelection()
            state["texture_image"] = self._texture_image_chooser.value
            state["depth"] = self._depth_slider.GetValue()
            
            _object.image = self._clipping_image_chooser.value 
            
        self._memento[_object] = state
        
            
    def _start_interaction(self, object, event):
        button = event[:-len("PressEvent")]
        self._button = button
        self._tools[button].start_interaction()
        
    
    def _stop_interaction(self, object, event):
        button = event[:-len("ReleaseEvent")]
        self._button = None
        self._tools[button].stop_interaction() 
    
    
    def _update_objects_check_list_box(self):
        self._objects_check_list_box.Clear()
        
        names = [object.name for object in self._objects_3d]
        self._objects_check_list_box.InsertItems(names, 0)
        for i in range(len(self._objects_3d)) :
            state = self._memento[self._objects_3d[i]]
            if state["static_visibility"] :
                self._objects_check_list_box.Check(i)
                
    
    def _update_color(self):
        """ Change one the color of the currently
            selected object
        """
        
        button = getattr(self, "_color_button")
        
        color_data = wx.ColourData()
        color_data.SetColour(button.GetBackgroundColour())
        
        dialog = wx.ColourDialog(self, color_data)
        if dialog.ShowModal() == wx.ID_OK :
            new_color = dialog.GetColourData().GetColour()
            button.SetBackgroundColour(new_color)
            index = self._objects_check_list_box.GetSelection()
            object = self._objects_3d[index]
            setattr(object, "color", [c/255. for c in new_color])
            self._viewer_3d.render_window_interactor.Render()
            
                
    def _update_memento(self, liste):
        for _object in liste :
            if _object not in self._memento :
                self._save_object_state(_object, default=True)
            
    
          
            
    
    
    
    
    



if __name__ == "__main__" :
    from vtk import (vtkDataReader, vtkPolyDataReader, 
                     vtkRectilinearGridReader, vtkStructuredGridReader, 
                     vtkUnstructuredGridReader, vtkVRMLImporter)
     
    app = wx.App()
    app.objects_3d = ObservableList()
    
    def load_object_3d(self, path, viewer_3d):
        generic_reader = vtkDataReader()
        generic_reader.SetFileName(path)
        
        if generic_reader.OpenVTKFile() and generic_reader.ReadHeader() :
            if generic_reader.IsFileStructuredPoints() :
                raise medipy.base.Exception("Cannot read VTK structured points")
            elif generic_reader.IsFilePolyData() :
                reader = vtkPolyDataReader()
                reader.SetFileName(path)
                object = Object3D(reader.GetOutput(), path)
            elif generic_reader.IsFileStructuredGrid() :
                reader = vtkStructuredGridReader()
                reader.SetFileName(path)
                object = Object3D(reader.GetOutput(), path)
            elif generic_reader.IsFileUnstructuredGrid() :
                reader = vtkUnstructuredGridReader()
                reader.SetFileName(path)
                object = Object3D(reader.GetOutput(), path)
            elif generic_reader.IsFileRectilinearGrid() :
                reader = vtkRectilinearGridReader()
                reader.SetFileName(path)
                object = Object3D(reader.GetOutput(), path)
            else : 
                raise medipy.base.Exception("Cannot read VTK file containing type %i"%generic_reader.GetFileType())
            self.objects_3d.append(object)
        else :
            importer = vtkVRMLImporter()
            importer.SetFileName(path)
            importer.Update()
            actors = importer.GetRenderer().GetActors()
            number_of_actors = actors.GetNumberOfItems()
            actors.InitTraversal()
            for i in range(number_of_actors) :
                actor = actors.GetNextItem()
                object = Object3D(path + ", %i"%i)
                object.actor.SetProperty(actor.GetProperty())
                object.actor.SetMapper(actor.GetMapper())
                self.objects_3d.append(object)
#        else :
#            raise Exception("Cannot load file %s of unknown type"%path)
    
    wx.App.load_object_3d = load_object_3d
    
    frame = Viewer3DFrame(None, app.objects_3d, size=(1200,800))
    frame.Show()
    app.SetTopWindow(frame)
    
    app.MainLoop()
