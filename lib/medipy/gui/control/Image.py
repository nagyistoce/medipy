##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import math
import weakref

import wx

import medipy.base
from medipy.base import Observable

class Image(wx.Panel, Observable):
    def __init__(self, parent, choices, value=None, output=False, 
                 output_checked=False, may_be_empty=False, 
                 may_be_empty_checked=False, *args, **kwargs):
        ##############
        # Initialize #
        ##############
        self._last_choosed = None
        
        self._choices = choices
        self._value = value
        self._output = output
        self._output_checked = output_checked
        self._may_be_empty = may_be_empty
        self._may_be_empty_checked = may_be_empty_checked
        
        if self._output and self._may_be_empty :
            raise medipy.base.Exception("may_be_empty and output are mutually exclusive")
        
        if not self._output and self._output_checked :
            raise medipy.base.Exception("output must exist to change its value")
        
        if not self._may_be_empty and self._may_be_empty_checked :
            raise medipy.base.Exception("may_be_empty must exist to change its value")
        
        if len(choices) > 0 : 
            if value :
                self._value = weakref.ref(value)
                self._last_choosed = weakref.ref(value)
            else :
                self._value = weakref.ref(choices[0])
                self._last_choosed = weakref.ref(choices[0])
        
            if self._output_checked or self._may_be_empty_checked :
                self._value = None
        else : 
            if self._output : 
                self._output_checked = True
        
        self._default_value = self._value
        
        self._default_output_checked = self._output_checked
        self._default_may_be_empty_checked = self._may_be_empty_checked
        
        wx.Panel.__init__(self, parent, *args, **kwargs)
        Observable.__init__(self, ["value"])
        
        # Widgets
        self._new_checkbox = wx.CheckBox(self, label="(new)")
        self._empty_checkbox = wx.CheckBox(self, label="None")
        self._new_checkbox.Hide()
        self._empty_checkbox.Hide()
        # Layout
        sizer = wx.BoxSizer(wx.VERTICAL)
        self._radiobuttons_sizer = wx.GridSizer()
        sizer.Add(self._empty_checkbox)
        sizer.Add(self._radiobuttons_sizer)
        sizer.Add(self._new_checkbox)
        self.SetSizer(sizer)
        # Events
        self._new_checkbox.Bind(wx.EVT_CHECKBOX, self.OnNewCheckBox)
        self._empty_checkbox.Bind(wx.EVT_CHECKBOX, self.OnEmptyCheckBox)
        self._choices.add_observer("any", self._on_choices_modified)
        
        self._update_gui()
        self.validate()
        
    ####################
    # Public interface #
    ####################
    def validate(self):
        if self._value in self._choices or self._new_checkbox or self._may_be_empty:
            self.SetBackgroundColour(None)
            return True
        else : 
            self.SetBackgroundColour(wx.RED)
            return False
        self.Fit()
        
    
    def reset(self):
        self._output_checked = self._default_output_checked
        self._may_be_empty_checked = self._default_may_be_empty_checked
        if self._default_value is not None :
            self._set_value(self._default_value())
        else :
            self._set_value(None)
        self._update_gui()
        self.validate()
    
    
    
    ##############
    # Properties #
    ##############
    def _set_choices(self, choices):
        self._choices.remove_observer("any", self._on_choices_modified)
        self._choices = choices
        if len(self._choices) > 0 : 
            if self._value :
                self._last_choosed = weakref.ref(self._value)
            else :
                self._last_choosed = weakref.ref(self._choices[0])
                if not self._output_checked and not self._may_be_empty_checked :
                    self._value = self._choices[0]
                
        else :
            if self._output :
                self.output_checked = True
        self._choices.add_observer("any", self._on_choices_modified)
        self._update_gui()
        self.validate()
        
        
    def _get_value(self):
        if self._value is not None :
            return self._value()
        else :
            return None
    
    
    def _set_value(self, value):
        if value is not None :
            self._value = weakref.ref(value)
            self._last_choosed = weakref.ref(value)
        else :
            self._value = None
        
        self.validate()
        self.notify_observers("value")
    
        
    def _get_default_value(self):
        return self._default_value
    
    
    def _set_default_value(self, value):
        if value is not None : 
            self._default_value = weakref.ref(value)
        else : 
            self._default_value = None
        self._update_gui()
        self.validate()
    
    
    def _set_output(self, value):
        self._output = value
        if self._output and self._may_be_empty :
            raise medipy.base.Exception("output and may_be_empty are mutually exclusive")
        self._update_gui()
        self.validate()
        
    
    def _set_output_checked(self, flag):
        if not self._output :
            raise medipy.base.Exception("output must exist to change its value")
        self._output_checked = flag
        self.notify_observers("value")
        self._update_gui()
        self.validate()
        
        
    def _set_may_be_empty(self, value):
        self._may_be_empty = value    
        if self._output and self._may_be_empty :
            raise medipy.base.Exception("may_be_empty and output are mutually exclusive")
        self._update_gui()
        self.validate()
        
    
    def _set_may_be_empty_checked(self, flag):
        if not self._may_be_empty :
            raise medipy.base.Exception("may_be_empty must exist to change its value")
        self._may_be_empty_checked = flag
        self.notify_observers("value")
        self._update_gui()
        self.validate()
    
    
    choices = property(lambda x:x._choices, _set_choices)
    value = property(_get_value, _set_value)
    default_value = property(_get_default_value, _set_default_value)
    output = property(lambda x:x._output, _set_output)
    output_checked = property(lambda x:x._output_checked, _set_output_checked)
    may_be_empty = property(lambda x:x._may_be_empty, _set_may_be_empty)
    may_be_empty_checked = property(lambda x:x._may_be_empty_checked, _set_may_be_empty_checked)

    
    
    
    ##########
    # Events #
    ##########
    def OnRadioButton(self, event):
        index = int(event.GetEventObject().GetLabel())-1
        self._set_value(self._choices[index])
    
    
    def OnNewCheckBox(self, event):
        if len(self._choices) > 0 :
            self._output_checked = not self._output_checked
            if self._output_checked :
                self._set_value(None)
            else :
                self._set_value(self._last_choosed())
        self._update_gui()
        self.notify_observers("value")
        
        
    def OnEmptyCheckBox(self, event):
        self._may_be_empty_checked = not self._may_be_empty_checked
        if self._may_be_empty_checked :
            self._set_value(None)
        else :
            self._set_value(self._last_choosed())
        
        self._update_gui()
        self.notify_observers("value")
        
    
    def _on_choices_modified(self, event):
        if self._last_choosed is None :
            if self._choices:
                self._last_choosed = weakref.ref(self._choices[0])
        else :
            if self._last_choosed() not in self._choices :
                if self._choices:
                    self._last_choosed = weakref.ref(self._choices[0])
                else :
                    self._last_choosed = None
        self._update_gui()
    
    #####################
    # Private interface #
    #####################
    def _update_gui(self):
        self._radiobuttons_sizer.Clear(True)
            
        self._new_checkbox.SetValue(self._output_checked)
        self._empty_checkbox.SetValue(self._may_be_empty_checked)
        
        # Re-shape the sizer
        nb_objects = len(self.choices)
        rows = max(int(math.ceil(math.sqrt(nb_objects))),1)
        self._radiobuttons_sizer.SetRows(rows)
        self._radiobuttons_sizer.SetCols(rows)
        if len(self.choices) == 0 :
            if self._output :
                self._new_checkbox.Show()
                self._new_checkbox.SetValue(self._output_checked)
            else : 
                self._new_checkbox.Hide()
            
            label = wx.StaticText(self, label=("(no image loaded)"))
            font = label.GetFont()
            font.SetStyle(wx.FONTSTYLE_ITALIC)
            label.SetFont(font)
            self._radiobuttons_sizer.Add(label, 1, wx.EXPAND)
        else :
            
            if self._output :
                self._new_checkbox.Show()
                self._new_checkbox.SetValue(self._output_checked)
            else : 
                self._new_checkbox.Hide()
            
            if self._may_be_empty :
                self._empty_checkbox.Show()
                self._empty_checkbox.SetValue(self._may_be_empty_checked)
            else :
                self._empty_checkbox.Hide()

            style=wx.RB_GROUP
            for i in range(0, len(self.choices)) :
                button = wx.RadioButton(self, -1, str(i+1), style=style)
                style=0
                button.Bind(wx.EVT_RADIOBUTTON, self.OnRadioButton)
                self._radiobuttons_sizer.Add(button, 0)
                flag1 = not (self._output and self._new_checkbox.IsChecked())
                flag2 = not (self._may_be_empty and self._empty_checkbox.IsChecked())
                enable_widget = flag1 and flag2
                button.Enable(enable_widget)
            
            if self._last_choosed is not None :
                index = self._choices.index(self._last_choosed())
            else : 
                index = 0
            choosed_button = self._radiobuttons_sizer.GetChildren()[index].GetWindow()
            choosed_button.SetValue(True)
        
        self.Fit()
