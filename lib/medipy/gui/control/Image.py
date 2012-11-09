##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import math
import wx
import medipy.base

class Image(wx.Panel, medipy.base.Observable):
    """ Control allowing the user to choose an item from a list of objects of 
        type medipy.base.Image. If output is True, then a "(new)" checkbox will
        be displayed ; when checked, the value of the control will be set to 
        None, and it will be the caller's responsibility to perform the 
        necessary  actions. If may_be_empty is True, then a "(none)" checkbox 
        will be displayed ; when checked, the value of the control will be set
        to None.
    """
    
    def __init__(self, parent, choices, value=None, output=False, 
                 output_checked=False, may_be_empty=False, 
                 may_be_empty_checked=False, *args, **kwargs):
        
        self._choices = None
        self._value = None
        self._default_value = None
        
        # Initialize
        wx.Panel.__init__(self, parent, *args, **kwargs)
        medipy.base.Observable.__init__(self, ["value"])
        
        # Widgets
        self._new_checkbox = wx.CheckBox(self, label="(new)")
        self._empty_checkbox = wx.CheckBox(self, label="(none)")
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
        
        self.value = value
        self.default_value = value
        
        self.choices = choices
        
        self.output = output
        self.output_checked = output_checked
        
        self.may_be_empty = may_be_empty
        self.may_be_empty_checked = may_be_empty_checked
        
        self.validate()
    
    def validate(self):
        valid = (self.output_checked or self.may_be_empty_checked or
                 (self._value is not None and self._value in [id(x) for x in self.choices]))
        
        if valid :
            self.SetBackgroundColour(None)
        else :
            self.SetBackgroundColour(wx.RED)
        
        return valid
    
    def reset(self):
        """ Reset the current choice to the default value, uncheck the "(new)"
            and "(none)" checkboxes.
        """
        
        self.output_checked = False
        self.may_be_empty_checked = False
        self.value = self.default_value
        self.validate()
    
    def update_gui(self):
        """ Update the GUI to reflect the current state of the control (value, 
            choices, output parameters and may_be_empty parameters).
        """
        
        if self._choices is None :
            return
        
        self._radiobuttons_sizer.Clear(True)
            
        # Re-shape the sizer to be as square as possible
        nb_objects = len(self.choices)
        rows = max(math.sqrt(len(self.choices)),1)
        rows = math.ceil(rows)
        self._radiobuttons_sizer.SetRows(rows)
        self._radiobuttons_sizer.SetCols(rows)
        
        if len(self.choices) == 0 :
            label = wx.StaticText(self, label="(no image loaded)")
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
                button.Enable(not self.output_checked and not self.may_be_empty_checked)
            
            ids = [id(x) for x in self.choices]
            try :
                index = ids.index(self._value)
            except ValueError :
                # No such value, just keep the current button checked
                pass
            else :
                button = self._radiobuttons_sizer.GetChildren()[index].GetWindow()
                button.SetValue(True)
        
        self.Layout()
    
    ##############
    # Properties #
    ##############
    
    def _get_value(self):
        """ The chosen Image from the list, or None if either the "(new)" 
            checkbox or the "(none)" checkbox is checked.
        """
        
        if self._value is None :
            return self._value
        else :
            ids = [id(x) for x in self.choices]
            index = ids.index(self._value)
            return self._choices[index]
    
    def _set_value(self, value):
        if value is not None :
            ids = [id(x) for x in self.choices]
            if id(value) in ids :
                self._value = id(value)
            else :
                self._value = None
        else :
            self._value = None
        self.update_gui()
        self.validate()
        
        self.notify_observers("value")
    
    def _get_default_value(self):
        """ The default value of the control, used when resetting.
        """
        
        return self._default_value
    
    def _set_default_value(self, default_value):
        if default_value is not None :
            ids = [id(x) for x in self.choices]
            if id(default_value) in ids :
                self._default_value = id(default_value)
            else :
                self._default_value = None
        else :
            self._default_value = None
    
    def _get_choices(self):
        """ The list of choices displayed by the control. If this is an 
            ObservableList, then the control will become an observer of the list
            and update when the list is modified. Otherwise, update_gui must
            be called explicitly.
        """
        
        return self._choices
    
    def _set_choices(self, choices):
        if isinstance(self._choices, medipy.base.ObservableList) :
            self._choices.remove_observer("any", self._on_choices_modified)
        self._choices = choices
        if isinstance(self._choices, medipy.base.ObservableList) :
            self._choices.add_observer("any", self._on_choices_modified)
        self.update_gui()
        self.validate()
    
    def _get_output(self):
        """ The visibility of the "(new)" checkbox.
        """
        
        return self._new_checkbox.IsShown()
    
    def _set_output(self, output):
        self._new_checkbox.Show(output)
        self.Layout()
        self.validate()
    
    def _get_output_checked(self):
        """ Is the "(new)" checkbox checked ?
        """
        return self._new_checkbox.IsShown() and self._new_checkbox.IsChecked()
    
    def _set_output_checked(self, output_checked):
        if output_checked :
            self.value = None
        else :
            if self.choices :
                self.value = self.choices[0]
            else :
                self.value = None
        self._new_checkbox.SetValue(output_checked)
        self.update_gui()
        self.validate()
    
    def _get_may_be_empty(self):
        """ The visibility of the "(none)" checkbox.
        """
        
        return self._empty_checkbox.IsShown()
    
    def _set_may_be_empty(self, may_be_empty):
        self._empty_checkbox.Show(may_be_empty)
        self.Layout()
        self.validate()
    
    def _get_may_be_empty_checked(self):
        """ Is the "(none)" checkbox checked ?
        """
        
        return self._empty_checkbox.IsShown() and self._empty_checkbox.IsChecked()
    
    def _set_may_be_empty_checked(self, may_be_empty_checked):
        if may_be_empty_checked :
            self.value = None
        else :
            if self.choices :
                self.value = self.choices[0]
            else :
                self.value = None
        self._empty_checkbox.SetValue(may_be_empty_checked)
        self.update_gui()
        self.validate()
    
    value = property(_get_value, _set_value)
    default_value = property(_get_default_value, _set_default_value)
    choices = property(_get_choices, _set_choices)
    output = property(_get_output, _set_output)
    output_checked = property(_get_output_checked, _set_output_checked)
    may_be_empty = property(_get_may_be_empty, _set_may_be_empty)
    may_be_empty_checked = property(_get_may_be_empty_checked, 
                                    _set_may_be_empty_checked)

    ##################
    # Event handlers #
    ##################
    
    def OnRadioButton(self, event):
        index = int(event.GetEventObject().GetLabel())-1
        self._set_value(self._choices[index])
    
    def OnNewCheckBox(self, event):
        if self.output_checked :
            self.value = None
        else :
            if self._choices :
                self.value = self._choices[0]
            else :
                self.value = None
        self.update_gui()
        self.validate()
        event.Skip()
    
    def OnEmptyCheckBox(self, event):
        if self.may_be_empty_checked :
            self.value = None
        else :
            if self._choices :
                self.value = self._choices[0]
            else :
                self.value = None
        self.update_gui()
        self.validate()
        event.Skip()
    
    def _on_choices_modified(self, event):
        if self.choices and not self.output_checked and not self.may_be_empty_checked :
            self.value = self.choices[0]
        self.update_gui()
        self.validate()
