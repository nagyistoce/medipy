##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import wx
from medipy.base import Observable

class String(wx.Panel, Observable):
    def __init__(self, parent, value="", may_be_empty = True, *args, **kwargs):
        
        ####################
        # Public interface #
        ####################
        self.may_be_empty = may_be_empty
        
        ##############
        # Properties #
        ##############
        self._value = None
        self._default_value = None 
        
        ##################
        # Initialization #
        ##################
        kwargs["style"] = kwargs.get("style", 0) | wx.TAB_TRAVERSAL
        wx.Panel.__init__(self, parent, *args, **kwargs)
        Observable.__init__(self, ["value"])
        
        # Widgets
        self._text = wx.TextCtrl(self)
        # Layout
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.SetSizer(sizer)
        sizer.Add(self._text, 1)
        # Events
        self._text.Bind(wx.EVT_TEXT, self.OnText)
        
        self._set_value(value)
        self._set_default_value(value)
    
    ####################
    # Public interface #
    ####################
    def validate(self):
        if isinstance(self._value, (str, unicode)):
            if not self.may_be_empty and not self._value :
                self._text.SetBackgroundColour(wx.RED)
                return False
            else :
                self._text.SetBackgroundColour(None)
                return True
        else : 
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
    
    value = property(lambda x:x._value, _set_value)
    default_value = property(lambda x:x._default_value, _set_default_value)
    
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

if __name__ == "__main__" :
    import sys
    
    app = wx.App()
    
    frame = wx.Frame(None)
    app.SetTopWindow(frame)
    
    ok_button = wx.Button(frame, id=wx.ID_OK, label="OK")
    reset_button = wx.Button(frame, id=wx.ID_RESET, label="Reset")
    
    ok_button.Bind(wx.EVT_BUTTON, lambda event : sys.stdout.write("Value is %s (%s)\n"%(control.value, "valid" if control.validate() else "invalid")))
    reset_button.Bind(wx.EVT_BUTTON, lambda event : control.reset())
    
    control = String(frame)
    control.value = "42"
    control.default_value = "foo"
    
    sizer = wx.BoxSizer(wx.VERTICAL)
    
    sizer.Add(control, flag=wx.EXPAND)
    
    button_sizer = wx.BoxSizer(wx.HORIZONTAL)
    button_sizer.Add(ok_button)
    button_sizer.Add(reset_button)
    sizer.Add(button_sizer)
    
    frame.SetSizer(sizer)
    sizer.SetSizeHints(frame)
    sizer.Layout()
    
    frame.Show()
    app.MainLoop()