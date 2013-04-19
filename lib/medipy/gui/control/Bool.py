##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import wx
from medipy.base import Observable

class Bool(wx.Panel, Observable):
    def __init__(self, parent, value=False, *args, **kwargs):
        ##############
        # Properties #
        ##############
        self._value = None
        self._default_value = None
        
        ##################
        # Initialization #
        ##################
        wx.Panel.__init__(self, parent, *args, **kwargs)
        Observable.__init__(self, ["value"])
        
        # Widgets
        self._checkbox = wx.CheckBox(self)
        # Layout
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.SetSizer(sizer)
        sizer.Add(self._checkbox, 0)
        
        # Events
        self._checkbox.Bind(wx.EVT_CHECKBOX, self.OnCheckBox)
        
        self._set_value(value)
        self._set_default_value(value)
    
    
    def validate(self):
        return isinstance(self._value, bool)
    
    def reset(self):
        self._set_value(self._default_value)
    
    ##############
    # Properties #
    ##############
    
    def _set_value(self, value):
        self._value = value
        self._checkbox.SetValue(value)
        self.notify_observers("value")
    
    def _set_default_value(self, default_value):
        self._default_value = default_value
    
    value = property(lambda x:x._value, _set_value)
    default_value = property(lambda x:x._default_value, _set_default_value)
    
    ##################
    # Event handlers #
    ##################
    
    def OnCheckBox(self,event):
        self._set_value(self._checkbox.IsChecked())
