##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import os.path

import wx

from medipy.base import Observable
import medipy.gui.base

class File(wx.Panel, Observable):
    def __init__(self, parent, value="", output=False, may_be_empty = False, wildcard="(all files)|*", *args, **kwargs):
        
        self.may_be_empty = may_be_empty
        self._preferences = medipy.gui.base.Preferences(
            wx.GetApp().GetAppName(), wx.GetApp().GetVendorName())
        
        ##############
        # Properties #
        ##############
        self._value = None
        self._default_value = None
        self._output = False 
        self._wildcard = None
        
        ##################
        # Initialization #
        ##################
        kwargs["style"] = kwargs.get("style", 0) | wx.TAB_TRAVERSAL
        wx.Panel.__init__(self, parent, *args, **kwargs)
        Observable.__init__(self, ["value"])
        
        # Widgets
        self._text = wx.TextCtrl(self)
        self._button = wx.Button(self, label="...", size=(30, -1))
        # Layout
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.SetSizer(sizer)
        sizer.Add(self._text, 1, wx.ALIGN_CENTER)
        sizer.Add(self._button, flag=wx.ALIGN_CENTER)
        # Events
        self._text.Bind(wx.EVT_TEXT, self.OnText)
        self._button.Bind(wx.EVT_BUTTON, self.OnButton)
        
        self._set_value(value)
        self._set_default_value(value)
        self._set_output(output)
        self._set_wildcard(wildcard)
    
    ####################
    # Public interface #
    ####################
    def validate(self):
        if not self._output and os.path.isfile(self._value) :
            self._text.SetBackgroundColour(None)
            return True
        elif self._output :
            self._text.SetBackgroundColour(None)
            return True
        elif not self._value and self.may_be_empty :
            self._text.SetBackgroundColour(None)
            return True
        else : # self._output and not os.path.isfile(self._value)
            self._text.SetBackgroundColour(wx.RED)
            return False
    
    def reset(self):
        self._set_value(self._default_value)
    
    ##############
    # Properties #
    ##############
    def _set_value(self, value):
        self._value = value
        try :
            insertion_point = self._text.GetInsertionPoint()
            self._text.ChangeValue(str(self._value))
            self._text.SetInsertionPoint(insertion_point)
        except :
            pass
        self.validate()
        self.notify_observers("value")
    
    def _set_default_value(self, value):
        self._default_value = value
    
    def _set_output(self, value):
        self._output = value
    
    def _set_wildcard(self, value):
        self._wildcard = value
    
    value = property(lambda x:x._value, _set_value)
    default_value = property(lambda x:x._default_value, _set_default_value)
    output = property(lambda x:x._output, _set_output)
    wildcard = property(lambda x:x._wildcard, _set_wildcard)
    
    ##########
    # Events #
    ##########
    def OnText(self, event):
        try :
            value = self._text.GetValue()
            self._set_value(value) 
        except Exception, e :
            self._value = None
        self.validate()
    
    def OnButton(self, event):
        path = self._preferences.get(self._get_preferences_entry(), "")
        
        dialog = wx.FileDialog(self, "Choose a File", wildcard = self._wildcard, defaultFile=path)
        if isinstance(self._value, (str, unicode)) :
            dialog.SetPath(self._value)
        if dialog.ShowModal() == wx.ID_OK:
            self._preferences.set(self._get_preferences_entry(), "")
            self._set_value(dialog.GetPath())
            self.validate()
    
    #####################
    # Private interface #
    #####################
    def _get_preferences_entry(self) :
        if self.output :
                entry = "IO/save_path"
        else :
            entry = "IO/load_path"
        return entry
