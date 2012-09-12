##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import xml.etree.ElementTree

import wx

import medipy.base
import medipy.io.image_annotation

from annotations_panel import AnnotationsPanel 

class AnnotationsDialog(wx.Dialog):
    
    def __init__(self, parent, wx_id=wx.ID_ANY, title="", 
                 pos=wx.DefaultPosition, size=wx.DefaultSize,
                 style=wx.DEFAULT_DIALOG_STYLE, name="annotations_dialog"):
        
        wx.Dialog.__init__(self, parent, wx_id, title, pos, size, style, name)
        
        # Widgets
        self._annotations_panel = AnnotationsPanel(self)
        self._add = wx.BitmapButton(
            self, bitmap = wx.ArtProvider.GetBitmap(wx.ART_NORMAL_FILE))
        self._delete = wx.BitmapButton(
            self, bitmap = wx.ArtProvider.GetBitmap(wx.ART_CROSS_MARK))
        self._load = wx.BitmapButton(
            self, bitmap = wx.ArtProvider.GetBitmap(wx.ART_FILE_OPEN))
        self._save = wx.BitmapButton(
            self, bitmap = wx.ArtProvider.GetBitmap(wx.ART_FILE_SAVE))
        
        # Layout
        buttons_sizer = wx.BoxSizer(wx.HORIZONTAL)
        buttons_sizer.Add(self._add)
        buttons_sizer.Add(self._delete)
        buttons_sizer.Add(self._load)
        buttons_sizer.Add(self._save)
        
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(self._annotations_panel, 1, wx.EXPAND)
        main_sizer.Add(buttons_sizer, 0, wx.EXPAND)
        
        self.SetSizer(main_sizer)
        main_sizer.SetSizeHints(self)
        
        # Events
        self._add.Bind(wx.EVT_BUTTON, self.OnAdd)
        self._delete.Bind(wx.EVT_BUTTON, self.OnDelete)
        self._load.Bind(wx.EVT_BUTTON, self.OnLoad)
        self._save.Bind(wx.EVT_BUTTON, self.OnSave)
    
    ##############
    # Properties #
    ##############
    
    def _get_image(self):
        return self._annotations_panel.image
    
    def _set_image(self, image):
        self._annotations_panel.image = image
    
    image = property(_get_image, _set_image)
    
    ##################
    # Event handlers #
    ##################
    
    def OnAdd(self, event):
        self._annotations_panel.add_annotation()
        
    def OnDelete(self, event):
        self._annotations_panel.delete_selected_annotations()
    
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
