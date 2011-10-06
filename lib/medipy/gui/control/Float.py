##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import wx
from medipy.base import Observable

class Float(wx.Panel, Observable):
    def __init__(self, parent, value=0, range=None, *args, **kwargs):
        
        ##############
        # Properties #
        ##############
        self._value = None
        self._default_value = None 
        self._range = None
        
        ##################
        # Initialization #
        ##################
        kwargs["style"] = kwargs.get("style", 0) | wx.TAB_TRAVERSAL
        wx.Panel.__init__(self, parent, *args, **kwargs)
        Observable.__init__(self, ["value"])
        
        # Widgets
        self._text = wx.TextCtrl(self)
        self._slider = wx.Slider(self)
        # Layout
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.SetSizer(sizer)
        sizer.Add(self._text, 1, wx.ALIGN_CENTER)
        sizer.Add(self._slider, 2, wx.ALIGN_CENTER)
        # Events
        self._text.Bind(wx.EVT_TEXT, self.OnText)
        self._slider.Bind(wx.EVT_SLIDER, self.OnSlider)
        
        self._set_value(value)
        self._set_default_value(value)
        self._set_range(range)
    
    ####################
    # Public interface #
    ####################
    def validate(self):
        in_range = True
        if self._range is not None :
            in_range = (self._range[0] <= self._value <= self._range[1])
        
        if isinstance(self._value, (float, int)) and in_range :
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
        try :
            self._slider.SetValue(self._value)
        except :
            pass
        self.validate()
        self.notify_observers("value")
    
    def _set_default_value(self, value):
        self._default_value = value
    
    def _set_range(self, range):
        self._range = range
        if range is not None :
            if range[0] != range[1] :
                self._slider.SetMin(range[0])
                self._slider.SetMax(range[1])
                self._slider.Enable()
            else :
                self._slider.Disable()
            if self.validate() :
                self._slider.SetValue(self._value)
            self._slider.Show()
        else : 
            self._slider.Hide()
        self.Fit()
        self.validate()
    
    value = property(lambda x:x._value, _set_value)
    default_value = property(lambda x:x._default_value, _set_default_value)
    range = property(lambda x:x._range, _set_range)
    
    ##########
    # Events #
    ##########
    def OnText(self, event):
        try :
            value = float(self._text.GetValue())
            self._set_value(value) 
        except Exception, e :
            self._value = None
            self.notify_observers("value")
        self.validate()
    
    def OnSlider(self, event):
        self._set_value(self._slider.GetValue())
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
    
    control = Float(frame)
    control.value = 42.35
    control.range=(0.25,255.18)
    control.default_value = 0.25
    
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