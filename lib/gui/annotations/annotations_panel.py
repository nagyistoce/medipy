##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import colorsys
import os
import random

import numpy
import wx
import wx.xrc

from medipy.base import find_resource
from medipy.base import ImageAnnotation
import medipy.gui.base
import medipy.gui.xrc_wrapper

from medipy.gui.control import Float 

class AnnotationsPanel(medipy.gui.base.Panel):
    
    class UI(medipy.gui.base.UI):
        def __init__(self):
            self.position = None
            self.label = None
            self.shape = None
            self.size = None
            self.color = None
            self.filled = None
            self.comment = None
            self.annotations = None
            
            self.controls = ["position", "label", "shape", "size", "color",
                             "filled", "comment", "annotations"]
    
    def __init__(self, parent=None, *args, **kwargs):
        
        # Image whose annotations will be changed
        self._image = None
        # Mapping from shape ID to name
        self._annotation_id_to_shape = dict( 
            [ (getattr(ImageAnnotation.Shape, name), name)
                for name in dir(ImageAnnotation.Shape) if name[:2] != "__"
            ]
        )
        # User interface
        self.ui = AnnotationsPanel.UI()
        
        xrc_file = find_resource(
            os.path.join("resources", "gui", "annotations_panel.xrc"))
        wrappers = medipy.gui.xrc_wrapper.FloatXMLHandler()
        medipy.gui.base.Panel.__init__(self, xrc_file, "annotations_panel", 
            [wrappers], self.ui, self.ui.controls, parent, *args, **kwargs)
        
        for shape in sorted(self._annotation_id_to_shape.values()) : 
            self.ui.shape.Append(shape)
        
        # Events
        self.ui.position.Bind(wx.EVT_TEXT, self.OnPosition)
        self.ui.label.Bind(wx.EVT_TEXT, self.OnLabel)
        self.ui.shape.Bind(wx.EVT_CHOICE, self.OnShape)
        self.ui.size.add_observer("value", self._on_size)
        self.ui.color.Bind(wx.EVT_BUTTON, self.OnColor)
        self.ui.filled.Bind(wx.EVT_CHECKBOX, self.OnFilled)
        self.ui.comment.Bind(wx.EVT_TEXT, self.OnComment)
        self.ui.annotations.Bind(wx.EVT_LISTBOX, self.OnAnnotation)
    
    def add_annotation(self):
        """ Add an annotation with a random color at the current cursor position
            of the image.
        """
        
        annotation = ImageAnnotation(label = "New annotation")
        annotation.position = (self._image.cursor_physical_position 
                               if self._image.display_coordinates == "physical"
                               else self._image.cursor_index_position)
        annotation.size = 5
        annotation.filled = True
        annotation.color = colorsys.hsv_to_rgb(random.random(), 1, 1)
        self._image.annotations.append(annotation)
        self._image.render()
        self.ui.annotations.Append(annotation.label)
        
        self.select_annotation(-1)
    
    def delete_selected_annotations(self):
        """ Delete the selected annotations.
        """
        
        items = list(self.ui.annotations.GetSelections())
        # Sort in reverse order so that last items are deleted first, and the
        # indices stay correct
        items.sort(reverse=True)
        for item in items :
            del self._image.annotations[item]
            self.ui.annotations.Delete(item)
        
        if self._image.annotations :
            self.select_annotation(min(items[-1], self.ui.annotations.GetCount()-1))
        else :
            self.ui.position.Clear()
            self.ui.label.Clear()
            self.ui.size.value = 1
            self.ui.color.SetBackgroundColour(wx.NullColor)
        
        self._image.render()
    
    def select_annotation(self, index):
        if index < 0 :
            index = len(self._image.annotations)+index
        
        annotation = self._image.annotations[index]
        
        position = "("+", ".join([str(x) for x in reversed(annotation.position)])+")"
        
        self.ui.position.ChangeValue(position)
        self.ui.label.ChangeValue(annotation.label)
        self.ui.shape.SetStringSelection(self._annotation_id_to_shape[annotation.shape])
        self.ui.size.value = annotation.size
        self.ui.color.SetBackgroundColour([int(255.*c) for c in annotation.color])
        self.ui.filled.SetValue(annotation.filled)
        self.ui.comment.ChangeValue(annotation.comment)
        self.ui.annotations.SetSelection(index)
        
        self._image.cursor_physical_position = annotation.position
    
    ##############
    # Properties #
    ##############
    
    def _get_image(self):
        return self._image
    
    def _set_image(self, image):
        self._image = image
        
        self.ui.annotations.Clear()
        for annotation in image.annotations :
            self.ui.annotations.Append(annotation.label)
    
    image = property(_get_image, _set_image)
    
    ##################
    # Event handlers #
    ##################
    
    def OnPosition(self, event):
        annotation = self._get_selected_annotation()
        if annotation is None :
            return
        
        string = self.ui.position.GetValue()
        
        position = None
        
        if string[0]+string[-1] != "()" and string[0]+string[-1] != "[]" :
            pass
        else :
            position = []
            elements = string[1:-1].split(",")
            for element in elements :
                try :
                    x = float(element)
                except ValueError :
                    position = None
                    break
                else :
                    position.insert(0, x)
        
        if position is not None :
            self.ui.position.SetBackgroundColour(wx.NullColour)
            annotation.position = position
            self._image.render()
        else :
            self.ui.position.SetBackgroundColour(wx.RED)
    
    def OnLabel(self, event):
        annotation = self._get_selected_annotation()
        if annotation is None :
            return
        
        self.ui.annotations.SetString(self.ui.annotations.GetSelection(), 
                                      self.ui.label.GetValue())
        
        annotation.label = self.ui.label.GetValue()
        self._image.render()
    
    def OnShape(self, event):
        annotation = self._get_selected_annotation()
        if annotation is None :
            return
        
        shape_string = self.ui.shape.GetStringSelection()
        shape = getattr(ImageAnnotation.Shape, shape_string)
        annotation.shape = shape
        self._image.render()
    
    def _on_size(self, event):
        
        annotation = self._get_selected_annotation()
        if annotation is None :
            return
        
        annotation.size = self.ui.size.value
        self._image.render()
    
    def OnColor(self, event):
        
        annotation = self._get_selected_annotation()
        if annotation is None :
            return
        
        color_data = wx.ColourData()
        color_data.SetColour(self.ui.color.GetBackgroundColour())
        
        dialog = wx.ColourDialog(self, color_data)
        if dialog.ShowModal() != wx.ID_OK :
            return
        
        new_color = dialog.GetColourData().GetColour()
        self.ui.color.SetBackgroundColour(new_color)
        annotation.color = [c/255. for c in new_color]
        self._image.render()
    
    def OnFilled(self, event):
        annotation = self._get_selected_annotation()
        if annotation is None :
            return
        
        annotation.filled = self.ui.filled.GetValue()
        self._image.render()
    
    def OnComment(self, event):
        annotation = self._get_selected_annotation()
        if annotation is None :
            return
        
        annotation.comment = self.ui.comment.GetValue()
    
    def OnAnnotation(self, event):
        index = self.ui.annotations.GetSelection()
        if index == -1 :
            return
        
        self.select_annotation(index)

    #####################
    # Private interface #
    #####################
    
    def _get_selected_annotation(self):
        if self.ui.annotations.GetSelection() == -1 :
            return None
        
        return self._image.annotations[self.ui.annotations.GetSelection()]
        