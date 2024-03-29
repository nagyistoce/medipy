##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import colorsys
import os
import random
import re
import xml.etree.ElementTree

import numpy
import wx
import wx.xrc

import medipy.base
from medipy.base import ImageAnnotation
import medipy.io.image_annotation
import medipy.gui.base
import medipy.gui.xrc_wrapper

from medipy.gui.control import Float 

class AnnotationsPanel(medipy.gui.base.Panel):
    """ GUI showing a list of annotations.
    """
    
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
            
            self.add = None
            self.delete = None
            self.load = None
            self.save = None
            
            self.controls = ["position", "label", "shape", "size", "color",
                             "filled", "comment", "annotations",
                             "add", "delete", "load", "save"]
    
    def __init__(self, parent=None, *args, **kwargs):
        # History of commands
        self.history = None
        # Image whose annotations will be changed
        self._image = None
        # User interface
        self.ui = AnnotationsPanel.UI()
        
        xrc_file = medipy.base.find_resource(
            os.path.join("resources", "gui", "annotations_panel.xrc"))
        wrappers = [medipy.gui.xrc_wrapper.FloatXMLHandler()]
        medipy.gui.base.Panel.__init__(self, xrc_file, "annotations_panel", 
            wrappers, self.ui, self.ui.controls, parent, *args, **kwargs)
        
        for shape in medipy.base.ImageAnnotation.Shape.members : 
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
        self.ui.add.Bind(wx.EVT_BUTTON, self.OnAdd)
        self.ui.delete.Bind(wx.EVT_BUTTON, self.OnDelete)
        self.ui.load.Bind(wx.EVT_BUTTON, self.OnLoad)
        self.ui.save.Bind(wx.EVT_BUTTON, self.OnSave)
        
        self.ui.delete.Enable(False)
        self.ui.save.Enable(False)
    
    def add_annotation(self):
        """ Add an annotation with a random color at the current cursor position
            of the image.
        """
        
        labels = []
        for annotation in self._image.annotations :
            match = re.match(r"Annotation (\d+)", annotation.label)
            if match :
                try :
                    value = int(match.group(1))
                except :
                    pass
                else :
                    labels.append(value)
        if labels :
            label = 1+max(labels)
        else :
            label = 1
        annotation = medipy.base.ImageAnnotation(label = "Annotation {0}".format(label))
        annotation.position = self._image.cursor_physical_position
        
        shape_string = self.ui.shape.GetStringSelection()
        shape = getattr(medipy.base.ImageAnnotation.Shape, shape_string)
        annotation.shape = shape
        
        annotation.size = self.ui.size.value
        annotation.filled = self.ui.filled.IsChecked()
        annotation.color = colorsys.hsv_to_rgb(random.random(), 1, 1)
        
        command = AddAnnotation(self, annotation)
        if self.history :
            self.history.add(command)
        else :
            command.execute()
    
    def delete_selected_annotations(self):
        """ Delete the selected annotations.
        """
        
        items = list(self.ui.annotations.GetSelections())
        # Sort in reverse order so that last items are deleted first, and the
        # indices stay correct
        items.sort(reverse=True)
        for item in items :
            command = DeleteAnnotation(self, item)
            if self.history :
                self.history.add(command)
            else :
                command.execute()
    
    def select_annotation(self, index):
        if index < 0 :
            index = len(self._image.annotations)+index
        
        annotation = self._image.annotations[index]
        
        position = "("+", ".join([str(x) for x in reversed(annotation.position)])+")"
        
        self.ui.position.ChangeValue(position)
        self.ui.label.ChangeValue(annotation.label)
        self.ui.shape.SetStringSelection(medipy.base.ImageAnnotation.Shape[annotation.shape])
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
        
        if self._image is not None :
            for annotation in image.annotations :
                self.ui.annotations.Append(annotation.label)
            
            self.ui.delete.Enable(len(image.annotations)!=0)
            self.ui.save.Enable(len(image.annotations)!=0)
        else :
            self.ui.delete.Disable()
            self.ui.save.Disable()
    
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
        shape = getattr(medipy.base.ImageAnnotation.Shape, shape_string)
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
    
    def OnAdd(self, event):
        self.add_annotation()
        
    def OnDelete(self, event):
        self.delete_selected_annotations()
    
    def OnLoad(self, event):
        filename = wx.LoadFileSelector("Annotations", "*", "", self)
        if filename :
            try :
                tree = xml.etree.ElementTree.parse(filename)
                annotations = medipy.io.image_annotation.annotations_from_xml(tree.getroot())
            except Exception, e :
                wx.MessageBox("Cannot load file : {0}".format(e), 
                              "Cannot load annotations", wx.OK, self)
            else :
                self.image.annotations[:] = annotations
                self.image.render()
                self.ui.annotations.Clear()
                for annotation in self.image.annotations :
                    self.ui.annotations.Append(annotation.label)
                self.ui.delete.Enable(len(self.image.annotations)!=0)
                self.ui.save.Enable(len(self.image.annotations)!=0)
                
    
    def OnSave(self, event):
        filename = wx.SaveFileSelector("Annotations", "*", "", self)
        if filename :
            try :
                element = medipy.io.image_annotation.annotations_to_xml(self.image.annotations)
                with open(filename, "w") as f :
                    f.write(xml.etree.ElementTree.tostring(element))
            except medipy.base, e :
                wx.MessageBox("Cannot save file : {0}".format(e), 
                              "Cannot save annotations", wx.OK, self)

    #####################
    # Private interface #
    #####################
    
    def _get_selected_annotation(self):
        if self.ui.annotations.GetSelection() == -1 :
            return None
        
        return self._image.annotations[self.ui.annotations.GetSelection()]

class AddAnnotation(medipy.base.UndoableCommand) :
    """ Command for adding an annotation.
        
        Undo effect : remove the annotation.
    """
    
    def __init__(self, annotations_panel, new_annotation):
        medipy.base.UndoableCommand.__init__(self, "Add annotation")
        
        self.annotations_panel = annotations_panel
        self.new_annotation = new_annotation
    
    def execute(self) :
        self.annotations_panel.image.annotations.append(self.new_annotation)
        self.annotations_panel.image.render()
        self.annotations_panel.ui.annotations.Append(self.new_annotation.label)
        
        self.annotations_panel.ui.delete.Enable(True)
        self.annotations_panel.ui.save.Enable(True)
        
        self.annotations_panel.select_annotation(-1)
        
    def undo(self) :
        try :
            index = self.annotations_panel.image.annotations.index(self.new_annotation)
        except ValueError, e:
            # Annotation is not in panel. Weird, but nothing to do
            return
        delete_command = DeleteAnnotation(self.annotations_panel, index)
        delete_command.execute()

class DeleteAnnotation(medipy.base.UndoableCommand):
    """ Command for deleting an annotation.
        
        Undo effect : add the annotation.
    """
    
    def __init__(self, annotations_panel, annotation_index):
        medipy.base.UndoableCommand.__init__(self, "Delete annotation")
        
        self.annotations_panel = annotations_panel
        self.annotation_index = annotation_index
        self._annotation = None
    
    def execute(self):
        
        self._annotation = self.annotations_panel.image.annotations[self.annotation_index]
        
        del self.annotations_panel.image.annotations[self.annotation_index]
        self.annotations_panel.ui.annotations.Delete(self.annotation_index)
        
        if self.annotations_panel.image.annotations :
            self.annotations_panel.select_annotation(
                 min(self.annotation_index, 
                     self.annotations_panel.ui.annotations.GetCount()-1))
            self.annotations_panel.ui.delete.Enable(True)
            self.annotations_panel.ui.save.Enable(True)
        else :
            self.annotations_panel.ui.position.Clear()
            self.annotations_panel.ui.label.Clear()
            self.annotations_panel.ui.size.value = 1
            self.annotations_panel.ui.color.SetBackgroundColour(wx.NullColor)
            self.annotations_panel.ui.delete.Enable(False)
            self.annotations_panel.ui.save.Enable(False)
        
        self.annotations_panel.image.render()
    
    def undo(self):
        if self._annotation is None :
            raise medipy.base.Exception("Cannot undo")
        
        command = AddAnnotation(self.annotations_panel, self._annotation)
        command.execute()
    
    def _get_annotation(self):
        """ Annotation that was deleted
        """
        
        return self._annotation

    annotation = property(_get_annotation)
