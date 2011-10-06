##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import os.path

import wx
from medipy.base import Observable

class Directory(wx.Panel, Observable):
    def __init__(self, parent, value="", output=False, *args, **kwargs):
        
        ##############
        # Properties #
        ##############
        self._value = None
        self._default_value = None
        self._output = False 
        
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
        sizer.Add(self._text, 1)
        sizer.Add(self._button)
        # Events
        self._text.Bind(wx.EVT_TEXT, self.OnText)
        self._button.Bind(wx.EVT_BUTTON, self.OnButton)
        
        self._set_value(value)
        self._set_default_value(value)
        self._set_output(output)
    
    ####################
    # Public interface #
    ####################
    def validate(self):
        if not self._output and os.path.isdir(self._value) :
            self._text.SetBackgroundColour(None)
            return True
        elif self._output :
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
    
    value = property(lambda x:x._value, _set_value)
    default_value = property(lambda x:x._default_value, _set_default_value)
    output = property(lambda x:x._output, _set_output)
    
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
        path = wx.ConfigBase_Get().Read("LoadImagePath")
        
        dialog = wx.DirDialog(self, "Choose a Directory", defaultPath=path)
        if isinstance(self._value, (str, unicode)) :
            dialog.SetPath(self._value)
        if dialog.ShowModal() == wx.ID_OK:
            wx.ConfigBase_Get().Write("LoadImagePath", dialog.GetPath())
            wx.ConfigBase_Get().Flush()
            self._set_value(dialog.GetPath())
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
    
    control = Directory(frame)
    #control.output = True
    control.wildcard = "BMP and GIF files (*.bmp;*.gif)|*.bmp;*.gif|PNG files (*.png)|*.png"
    
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
