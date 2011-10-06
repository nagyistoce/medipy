##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import math
import weakref

import wx
from medipy.base import Observable

class Object3D(wx.Panel, Observable):
    def __init__(self, parent, choices, value=None, output_checked=False, *args, **kwargs):
        """ Control to select a 3D scene to which the Object3D will be appended
        """
        ##############
        # Properties #
        ##############
        self._value = None
        self._default_value = None
        self._choices = choices
        self._default_output_checked = None
        
        ##############
        # Initialize #
        ##############
        if value is None and len(choices)>0 :
            value = choices[0]
        
        wx.Panel.__init__(self, parent, *args, **kwargs)
        Observable.__init__(self, ["value"])
        
        # Widgets
        self._checkbox = wx.CheckBox(self, label="(new)")
        # Layout
        sizer = wx.BoxSizer(wx.VERTICAL)
        self._radiobuttons_sizer = wx.GridSizer()
        sizer.Add(self._radiobuttons_sizer)
        sizer.Add(self._checkbox)
        self.SetSizer(sizer)
        # Events
        self._checkbox.Bind(wx.EVT_CHECKBOX, self.OnCheckBox)
        
        self._set_default_value(value)
        self._set_output_checked(output_checked)
        self._set_default_output_checked(output_checked)
        self._set_choices(choices)
        self._set_value(value)
        
        self._choices.add_observer("any", self._on_choices_modified)
        
    ####################
    # Public interface #
    ####################
    def validate(self):
        if self.value in self._choices or self._checkbox.IsChecked() :
            self.SetBackgroundColour(None)
            return True
        else : 
            self.SetBackgroundColour(wx.RED)
            return False
    
    def reset(self):
        self._set_value(self.default_value)
        self._set_output_checked(self._default_output_checked)
        self.validate()
    
    ##############
    # Properties #
    ##############
    def _set_value(self, value):
        if value is not None :
            self._value = weakref.ref(value)
        else : 
            self._value = None
        try :
            index = self._choices.index(value)
            button = self._radiobuttons_sizer.GetChildren()[index].GetWindow()
            button.SetValue(True)
        except Exception, e:
            pass
        self.validate()
        self.notify_observers("value")
    
    def _set_default_value(self, value):
        if value is not None : 
            self._default_value = weakref.ref(value)
        else : 
            self._default_value = None
    
    def _set_choices(self, choices):
        self._choices[:] = choices
        self._update_gui()
    
    def _set_output(self, value):
        self._output = value
        
        if self._output :
            self._checkbox.Show()
        else : 
            self._checkbox.Hide()
        self.Fit()
    
    def _set_output_checked(self, value):
        self._checkbox.SetValue(value)
        self._update_gui()
        self.notify_observers("value")
    
    def _set_default_output_checked(self, value):
        self._default_output_checked = value
    
    value = property(lambda x:x._value() if x._value is not None else None, _set_value)
    default_value = property(lambda x:x._default_value() if x._default_value is not None else None, _set_default_value)
    choices = property(lambda x:x._choices, _set_choices)
    output = property(lambda x:x._output, _set_output)
    output_checked = property(lambda x:x._checkbox.IsChecked(), _set_output_checked)
    default_output_checked = property(lambda x:x._default_output_checked, _set_default_output_checked)
    
    ##########
    # Events #
    ##########
    def OnRadioButton(self, event):
        index = int(event.GetEventObject().GetLabel())-1
        self._set_value(self._choices[index])
        self.validate()
    
    def OnCheckBox(self, event):
        self._update_gui()
        self.validate()
        self.notify_observers("value")
    
    def _on_choices_modified(self, event):
        self._update_gui()
    
    #####################
    # Private interface #
    #####################
    def _update_gui(self):
        self._radiobuttons_sizer.Clear(True)
        
        # Re-shape the sizer
        nb_objects = len(self.choices)
        rows = max(int(math.ceil(math.sqrt(nb_objects))),1)
        self._radiobuttons_sizer.SetRows(rows)
        self._radiobuttons_sizer.SetCols(rows)
        
        if len(self.choices) == 0 :
            label = wx.StaticText(self, label=("(no image loaded)"))
            font = label.GetFont()
            font.SetStyle(wx.FONTSTYLE_ITALIC)
            label.SetFont(font)
            self._radiobuttons_sizer.Add(label, 1, wx.EXPAND)
        else :
            style=wx.RB_GROUP
            for i in range(0, len(self.choices)) :
                button = wx.RadioButton(self, -1, str(i+1), style=style)
                style=0
                button.Bind(wx.EVT_RADIOBUTTON, self.OnRadioButton)
                self._radiobuttons_sizer.Add(button, 0)
                enable_widget = not self._checkbox.IsChecked()
                button.Enable(enable_widget)
        
        self.Fit()

if __name__ == "__main__" :
    import sys
    
    from medipy.base import ObservableList
    from medipy.base import Image as ImageBase
    from medipy.gui.viewer_3d_frame import Viewer3DFrame
    
    app = wx.App()
    
    frame = wx.Frame(None)
    app.SetTopWindow(frame)
    
    ok_button = wx.Button(frame, id=wx.ID_OK, label="OK")
    reset_button = wx.Button(frame, id=wx.ID_RESET, label="Reset")
    
    ok_button.Bind(wx.EVT_BUTTON, lambda event : sys.stdout.write("Value is %s (%s)\n"%(control.value, "valid" if control.validate() else "invalid")))
    reset_button.Bind(wx.EVT_BUTTON, lambda event : control.reset())
    
    scenes = ObservableList()
    for i in range(0) :
        scenes.append(Viewer3DFrame(None, ObservableList()))
    
    control = Object3D(frame, 
                    scenes, 
#                    value="machin", 
                    output_checked=False)
    
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
