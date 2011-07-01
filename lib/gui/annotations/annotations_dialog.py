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
import xml.dom.minidom as minidom

import numpy
import wx
import wx.xrc

from medipy.base import find_resource
from medipy.base import ImageAnnotation

from medipy.gui.control import Float 

class AnnotationsDialog(wx.Dialog):
    
    def __init__(self, parent, wx_id=wx.ID_ANY, pos=wx.DefaultPosition, size=wx.DefaultSize,
                 style=wx.DEFAULT_DIALOG_STYLE, name="annotations_dialog"):
        
        self._image = None
        self._current_annotation = None
        
        # Find all possible shapes, build a mapping from ID to name
        self._annotation_id_to_shape = dict( 
            [ (getattr(ImageAnnotation.Shape, name), name)
                for name in dir(ImageAnnotation.Shape) if name[:2] != "__"
            ]
        )
        
        super(AnnotationsDialog, self).__init__(parent, wx_id, "", pos, size, style, name)
        
        # Widgets
        self._position = wx.TextCtrl(self)
        self._label = wx.TextCtrl(self)
        self._shape = wx.Choice(self)
        self._size = Float(self, 1, (1, 100))
        self._color = wx.Button(self, size=(20,20))
        self._filled = wx.CheckBox(self, label = "Filled")
        self._comment = wx.TextCtrl(self, style=wx.TE_MULTILINE)
        self._annotations = wx.ListBox(self, style=wx.LB_SINGLE)
        self._add = wx.Button(self, label="Add")
        self._delete = wx.Button(self, label = "Delete")
        
        for shape in self._annotation_id_to_shape.values() : 
            self._shape.Append(shape)
        
        # Layout
        controls_sizer = wx.FlexGridSizer(0, 2)
        controls_sizer.AddGrowableCol(1)
        controls_sizer.Add(wx.StaticText(self, label="Position : "), 0, wx.ALIGN_CENTER)
        controls_sizer.Add(self._position, flag=wx.EXPAND)
        controls_sizer.Add(wx.StaticText(self, label="Label : "), 0, wx.ALIGN_CENTER)
        controls_sizer.Add(self._label, flag=wx.EXPAND)
        controls_sizer.Add(wx.StaticText(self, label="Shape : "), 0, wx.ALIGN_CENTER)
        controls_sizer.Add(self._shape)
        controls_sizer.Add(wx.StaticText(self, label="Size : "), 0, wx.ALIGN_CENTER)
        controls_sizer.Add(self._size, flag=wx.EXPAND)
        
        controls_sizer.Add(wx.StaticText(self, label="Color : "), 0, wx.ALIGN_CENTER)
        color_and_filled_sizer = wx.BoxSizer(wx.HORIZONTAL)
        color_and_filled_sizer.Add(self._color, flag=wx.ALIGN_CENTER)
        color_and_filled_sizer.AddSpacer(30)
        color_and_filled_sizer.Add(self._filled, flag=wx.ALIGN_CENTER)
        controls_sizer.Add(color_and_filled_sizer)
        
        controls_sizer.Add(wx.StaticText(self, label="Comment : "), 0, wx.ALIGN_TOP)
        controls_sizer.Add(self._comment, 1, wx.EXPAND)
        
        controls_sizer.Layout()
        
        buttons_sizer = wx.BoxSizer(wx.HORIZONTAL)
        buttons_sizer.Add(self._add)
        buttons_sizer.Add(self._delete)
        
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(controls_sizer, 0, wx.EXPAND)
        main_sizer.Add(self._annotations, 0, wx.EXPAND)
        main_sizer.Add(buttons_sizer, 0, wx.EXPAND)
        
        self.SetSizer(main_sizer)
        main_sizer.SetSizeHints(self)
        
        # Events
        self._position.Bind(wx.EVT_TEXT, self.OnPosition)
        self._label.Bind(wx.EVT_TEXT, self.OnLabel)
        self._shape.Bind(wx.EVT_CHOICE, self.OnShape)
        self._size.add_observer("value", self._on_size)
        self._color.Bind(wx.EVT_BUTTON, self.OnColor)
        self._filled.Bind(wx.EVT_CHECKBOX, self.OnFilled)
        self._comment.Bind(wx.EVT_TEXT, self.OnComment)
        self._annotations.Bind(wx.EVT_LISTBOX, self.OnAnnotation)
        self._add.Bind(wx.EVT_BUTTON, self.OnAdd)
        self._delete.Bind(wx.EVT_BUTTON, self.OnDelete)
    
    def select_annotation(self, index):
        if index < 0 :
            index = len(self._image.annotations)+index
        
        annotation = self._image.annotations[index]
        
        position = "("+", ".join([str(x) for x in reversed(annotation.position)])+")"
        
        self._position.SetValue(position)
        self._label.ChangeValue(annotation.label)
        self._shape.SetStringSelection(self._annotation_id_to_shape[annotation.shape])
        self._size.value = annotation.size
        self._color.SetBackgroundColour([int(255.*c) for c in annotation.color])
        self._filled.SetValue(annotation.filled)
        self._comment.SetValue(annotation.comment)
        self._annotations.SetSelection(index)
        
        self._image.cursor_physical_position = annotation.position
    
    ##############
    # Properties #
    ##############
    
    def _get_image(self):
        return self._image
    
    def _set_image(self, image):
        self._image = image
        
        self._annotations.Clear()
        for annotation in image.annotations :
            self._annotations.Append(annotation.label)
    
    image = property(_get_image, _set_image)
    
    ##################
    # Event handlers #
    ##################
    
    def OnPosition(self, event):
        annotation = self._get_selected_annotation()
        if annotation is None :
            return
        
        string = self._position.GetValue()
        
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
            self._position.SetBackgroundColour(wx.NullColour)
            annotation.position = position
            self._image.render()
        else :
            self._position.SetBackgroundColour(wx.RED)
    
    def OnLabel(self, event):
        annotation = self._get_selected_annotation()
        if annotation is None :
            return
        
        self._annotations.SetString(self._annotations.GetSelection(), 
                                    self._label.GetValue())
        
        annotation.label = self._label.GetValue()
        self._image.render()
    
    def OnShape(self, event):
        annotation = self._get_selected_annotation()
        if annotation is None :
            return
        
        shape_string = self._shape.GetStringSelection()
        shape = getattr(ImageAnnotation.Shape, shape_string)
        annotation.shape = shape
        self._image.render()
    
    def _on_size(self, event):
        
        annotation = self._get_selected_annotation()
        if annotation is None :
            return
        
        annotation.size = self._size.value
        self._image.render()
    
    def OnColor(self, event):
        
        annotation = self._get_selected_annotation()
        if annotation is None :
            return
        
        color_data = wx.ColourData()
        color_data.SetColour(self._color.GetBackgroundColour())
        
        dialog = wx.ColourDialog(self, color_data)
        if dialog.ShowModal() != wx.ID_OK :
            return
        
        new_color = dialog.GetColourData().GetColour()
        self._color.SetBackgroundColour(new_color)
        annotation.color = [c/255. for c in new_color]
        self._image.render()
    
    def OnFilled(self, event):
        annotation = self._get_selected_annotation()
        if annotation is None :
            return
        
        annotation.filled = self._filled.GetValue()
        self._image.render()
    
    def OnComment(self, event):
        annotation = self._get_selected_annotation()
        if annotation is None :
            return
        
        annotation.comment = self._comment.GetValue()
    
    def OnAnnotation(self, event):
        index = self._annotations.GetSelection()
        if index == -1 :
            return
        
        self.select_annotation(index)
    
    def OnAdd(self, event):
        
        annotation = ImageAnnotation(label = "New annotation")
        annotation.position = (self._image.cursor_physical_position 
                               if self._image.display_coordinates == "physical"
                               else self._image.cursor_index_position)
        annotation.size = 5
        annotation.filled = True
        annotation.color = colorsys.hsv_to_rgb(random.random(), 1, 1)
        self._image.annotations.append(annotation)
        self._image.render()
        self._annotations.Append(annotation.label)
        
        self.select_annotation(-1)
        
    def OnDelete(self, event):
        items = list(self._annotations.GetSelections())
        # Sort in reverse order so that last items are deleted first, and the
        # indices stay correct
        items.sort(reverse=True)
        for item in items :
            del self._image.annotations[item]
            self._annotations.Delete(item)
        
        if self._image.annotations :
            self.select_annotation(min(items[-1], self._annotations.GetCount()-1))
        else :
            self._position.Clear()
            self._label.Clear()
            self._size.value = 1
            self._color.SetBackgroundColour(wx.NullColor)
        
        self._image.render()

    #####################
    # Private interface #
    #####################
    
    def _get_selected_annotation(self):
        if self._annotations.GetSelection() == -1 :
            return None
        
        return self._image.annotations[self._annotations.GetSelection()]
        