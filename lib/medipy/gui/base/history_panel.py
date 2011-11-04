##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import wx

import medipy.base

class HistoryPanel(wx.Panel) :
    def __init__(self, *args, **kwargs) :
        wx.Panel.__init__(self, *args, **kwargs)
        
        self._listbox = wx.ListBox(self)
        sizer = wx.BoxSizer()
        sizer.Add(self._listbox, 1, wx.EXPAND)
        self.SetSizer(sizer)
        
        self._listbox.Bind(wx.EVT_LISTBOX, self.OnListBox)
        self._previous_selection = None
        
        self._history = medipy.base.History()
        self._update_listbox()
    
    def add(self, command) :
        self._history.add(command)
        self._update_listbox()
        self._listbox.SetSelection(self._listbox.GetCount()-1)
        self._previous_selection = self._listbox.GetSelection()
    
    def undo(self, count=1) : 
        self._history.undo(count)
        self._listbox.SetSelection(self._previous_selection-count)
        self._previous_selection = self._listbox.GetSelection()
    
    def redo(self, count=1) :
        self._history.redo(count)
        self._listbox.SetSelection(self._previous_selection+count)
        self._previous_selection = self._listbox.GetSelection()
    
    def _get_can_undo(self) :
        return self._history.can_undo
    
    def _get_can_redo(self) :
        return self._history.can_redo
    
    def _get_labels(self) :
        return self._history.labels
        
    can_undo = property(_get_can_undo)
    can_redo = property(_get_can_redo)
    labels = property(_get_labels)
    
    def OnListBox(self, event) :
        event.Skip()
        if self._listbox.GetSelection() == wx.NOT_FOUND :
            return
        offset = self._listbox.GetSelection()-self._previous_selection
        if offset >=0 :
            function = self.redo
        else :
            offset = -offset
            function = self.undo
        function(offset)
    
    def _update_listbox(self) :
        self._listbox.Clear()
        for label in self._history.labels[::-1] :
            self._listbox.Append(label)