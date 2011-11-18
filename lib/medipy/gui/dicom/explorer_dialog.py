##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import wx

from medipy.gui.dicom import Explorer
from medipy.gui import wxVTKRenderWindowInteractor

class ExplorerDialog(wx.Dialog):
    def __init__(self, *args, **kwargs):
        wx.Dialog.__init__(self, *args, **kwargs)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(sizer)
        sizer.SetSizeHints(self)
        
        self._explorer = Explorer(self)
        sizer.Add(self._explorer, 1, wx.EXPAND)
        
        buttons_sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(buttons_sizer, 0, wx.EXPAND)
        
        self._import_button = wx.Button(self, wx.ID_OK, "&Import")
        buttons_sizer.Add(self._import_button, 0)
        
        self._cancel_button = wx.Button(self, wx.ID_CANCEL, "&Cancel")
        buttons_sizer.Add(self._cancel_button, 0)
        
        self._explorer.hierarchy_tree.Bind(wx.EVT_TREE_SEL_CHANGED, self.OnInformationEntitySelected)
        self._cancel_button.Bind(wx.EVT_BUTTON, self.OnClose)
        self.Bind(wx.EVT_CLOSE, self.OnClose)
    
    def set_datasets(self, datasets):
        self._explorer.hierarchy_tree.set_datasets(datasets)
        self.OnInformationEntitySelected(wx.TreeEvent())
    
    def get_selected_datasets(self):
        return self._explorer.hierarchy_tree.selected_datasets
    
    ##################
    # Event handlers #
    ##################
    
    def OnInformationEntitySelected(self, event):
        event.Skip()
        if self.get_selected_datasets() :
            self._import_button.Enable(True)
        else :
            self._import_button.Enable(False)
    
    def OnClose(self, event):
        """ Fix a destruction bug that happens when the RenderWindowInteractors
            are destroyed too late : the drawable is not valid and an error is
            thrown.
        """
        
        # Close all renderwindow interactors beneath us
        queue = [self._explorer]
        while queue :
            window = queue.pop(0)
            for child in window.GetChildren() :
                if isinstance(child, wxVTKRenderWindowInteractor) :
                    if isinstance(event, wx.CloseEvent) :
                        child.Close(not event.CanVeto())
                    else :
                        child.Close()
                queue.append(child)
        
        event.Skip()