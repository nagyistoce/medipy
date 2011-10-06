##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import wx

from annotations_panel import AnnotationsPanel 

class AnnotationsDialog(wx.Dialog):
    
    def __init__(self, parent, wx_id=wx.ID_ANY, title="", 
                 pos=wx.DefaultPosition, size=wx.DefaultSize,
                 style=wx.DEFAULT_DIALOG_STYLE, name="annotations_dialog"):
        
        wx.Dialog.__init__(self, parent, wx_id, title, pos, size, style, name)
        
        # Widgets
        self._annotations_panel = AnnotationsPanel(self)
        self._add = wx.Button(self, label="Add")
        self._delete = wx.Button(self, label = "Delete")
        
        # Layout
        buttons_sizer = wx.BoxSizer(wx.HORIZONTAL)
        buttons_sizer.Add(self._add)
        buttons_sizer.Add(self._delete)
        
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(self._annotations_panel, 1, wx.EXPAND)
        main_sizer.Add(buttons_sizer, 0, wx.EXPAND)
        
        self.SetSizer(main_sizer)
        main_sizer.SetSizeHints(self)
        
        # Events
        self._add.Bind(wx.EVT_BUTTON, self.OnAdd)
        self._delete.Bind(wx.EVT_BUTTON, self.OnDelete)
    
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
