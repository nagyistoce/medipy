##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import wx
from medipy.base import Observable

class Enum(wx.Panel, Observable):
    def __init__(self, parent, choices, value=None, *args, **kwargs):
        
        ##############
        # Properties #
        ##############
        self._value = None
        self._default_value = None
        self._choices = None
        
        ##############
        # Initialize #
        ##############
        if value is None :
            value = choices[0]
        
        wx.Panel.__init__(self, parent, *args, **kwargs)
        Observable.__init__(self, ["value"])
        
        # Widgets
        self._choices_widget = wx.Choice(self)
        # Layout
        sizer = wx.BoxSizer()
        sizer.Add(self._choices_widget)
        self.SetSizer(sizer)
        self.Layout()
        # Events : 
        self._choices_widget.Bind(wx.EVT_CHOICE, self.OnChoiceChanged)
        
        self._set_choices(choices)
        self._set_default_value(value)
        self._set_value(value)
    
    ####################
    # Public interface #
    ####################
    def validate(self):
        if self._value in self._choices:
            self._choices_widget.SetBackgroundColour(None)
            return True
        else : 
            self._choices_widget.SetBackgroundColour(wx.RED)
            return False
    
    def reset(self):
        self._set_value(self._default_value)
    
    ##############
    # Properties #
    ##############
    def _set_value(self, value):
        self._value = value
        try :
            index = self._choices.index(value)
            self._choices_widget.SetSelection(index)
        except Exception, e:
            pass
        self.validate()
        self.notify_observers("value")
    
    def _set_default_value(self, value):
        self._default_value = value
    
    def _set_choices(self, choices):
        self._choices = choices
        self._choices_widget.Clear()
        for choice in self._choices :
            self._choices_widget.Append(str(choice))
    
    value = property(lambda x:x._value, _set_value)
    default_value = property(lambda x:x._default_value, _set_default_value)
    choices = property(lambda x:x._choices, _set_choices)
    
    ##########
    # Events #
    ##########
    def OnChoiceChanged(self, event):
        index = self._choices_widget.GetCurrentSelection()
        self._set_value(self._choices[index])
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
    
    control = Enum(frame, ["tata", "titi", "toto"])
#    control.choices = ["toto", "tutu"]
#    control.value = "tutu"
#    control.default_value = "toto"
    
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