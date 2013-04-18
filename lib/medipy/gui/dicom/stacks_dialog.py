##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import sys

import numpy
import wx

import medipy.io.dicom

class StacksDialog(wx.Dialog):
    """ Display a dialog for all stacks of a DICOM series (i.e. sets of images
        having a coherent orientation).
    """
    
    def __init__(self, parent, single_stack=False, *args, **kwargs) :
        self._stacks = []
        self._item_data = {}
        
        wx.Dialog.__init__(self, parent, *args, **kwargs)
        
        # Widgets
        style = (wx.LC_REPORT | wx.LC_SINGLE_SEL) if single_stack else wx.LC_REPORT
        self._listctrl = wx.ListCtrl(self, style=style)
        buttons_sizer = self.CreateButtonSizer(wx.OK|wx.CANCEL)
        self._ok_button = self.FindWindowById(self.GetAffirmativeId())
        
        # Layout
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self._listctrl, 1, wx.EXPAND)
        sizer.Add(buttons_sizer, 0, flag=wx.EXPAND)
        self.SetSizerAndFit(sizer)
        if "size" not in kwargs :
            self.SetSize((350, 300))
        
        # Events
        self._listctrl.Bind(wx.EVT_LIST_ITEM_SELECTED, self.OnSelectionChanged)
        self._listctrl.Bind(wx.EVT_LIST_ITEM_DESELECTED, self.OnSelectionChanged)
        
        self.OnSelectionChanged(wx.CommandEvent())
    
    def set_stacks(self, stacks):
        self._stacks = stacks
        self._update_listctrl()
    
    def get_selected_stacks(self):
        item = -1
        
        stacks = []
        
        # Find all selected items
        while True :
            item = self._listctrl.GetNextItem(item, state=wx.LIST_STATE_SELECTED)
            if item == -1 :
                break
            stacks.append(self._item_data[item])

        return stacks
    
    #===========================================================================
    # Event handlers
    #===========================================================================
    def OnSelectionChanged(self, event):
        if self._ok_button is not None :
            self._ok_button.Enable(self._listctrl.GetSelectedItemCount() > 0)
    
    #===========================================================================
    # Private interface
    #===========================================================================
    
    def _update_listctrl(self):
        """ Update the listctrl with the data contained in self._groups
        """
        
        self._listctrl.ClearAll()
        self._listctrl.InsertColumn(sys.maxint, "Orientation")
        self._listctrl.InsertColumn(sys.maxint, "Number of images")
        
        for stack in self._stacks :
            orientation = stack[0].image_orientation_patient.value
            value = self._format_image_orientation_patient(orientation)
            item = self._listctrl.InsertStringItem(sys.maxint, value)
            self._listctrl.SetStringItem(item, 1, "{0}".format(len(stack)))
            self._item_data[item] = stack
        
        self.OnSelectionChanged(wx.CommandEvent())
        
    def _format_image_orientation_patient(self, value):
        max_norm = 0.05
        
        def comparator(v1, v2):
            return numpy.linalg.norm(numpy.subtract(numpy.abs(v1), 
                                                    numpy.abs(v2)), numpy.inf)
        
        axial = [1, 0, 0, 0, 1, 0]
        coronal = [1, 0, 0, 0, 0, 1]
        sagittal = [0, 1, 0, 0, 0, 1]
        if comparator(axial, value) <= 0.05 :
            return "axial"
        elif comparator(coronal, value) <= 0.05 :
            return "coronal"
        elif comparator(sagittal, value) <= 0.05 :
            return "sagittal"
        else :
            return "(%.2f, %.2f, %.2f, %.2f, %.2f, %.2f)"%tuple(value)
