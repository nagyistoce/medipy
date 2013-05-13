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
    
    class UI(object) :
        def __init__(self, parent) :
            # Create widgets
            self.empty = wx.CheckBox(parent, label="(none)")
            self.choices = wx.Panel(parent)
            self.layers = wx.ComboBox(parent)
            self.new = wx.CheckBox(parent, label="(new)")
            
            # Hide empty and new at startup
            self.empty.Hide()
            self.new.Hide()
            
            # Layout
            self.choices_sizer = wx.GridSizer()
            self.choices.SetSizer(self.choices_sizer) 
            
            sizer = wx.BoxSizer(wx.VERTICAL)
            sizer.Add(self.empty)
            sizer.Add(self.choices)
            sizer.Add(self.layers)
            sizer.Add(self.new)
            
            parent.SetSizer(sizer)
        
    def __init__(self, parent, choices, value=None,
                 output=False, output_checked=False,
                 may_be_empty=False, may_be_empty_checked=False,
                 show_layers=False,
                 *args, **kwargs) :
    
        import medipy.gui.image
    
        self._choices = None
        self._value = None
        self._default_value = None
        
        wx.Panel.__init__(self, parent, *args, **kwargs)
        medipy.base.Observable.__init__(self, ["value"])
        
        self.ui = Image.UI(self)
        self.ui.new.Bind(wx.EVT_CHECKBOX, self._set_output_checked)
        self.ui.empty.Bind(wx.EVT_CHECKBOX, self._set_may_be_empty_checked)
        
        self.choices = choices
        
        self.value = value
        self.default_value = value
        
        self.output = output
        self.output_checked = output_checked
        
        self.may_be_empty = may_be_empty
        self.may_be_empty_checked = may_be_empty_checked
        
        self.show_layers = show_layers
    
    def validate(self):
        ids = [id(x) for x in self.choices]
        valid = (self.output_checked or self.may_be_empty_checked or
                 (self._value is not None and self._value in ids))
        
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
        if isinstance(value, wx.Event) :
            index = int(value.GetEventObject().GetLabel())-1
            self._set_value(self._choices[index])
        else :
            if value is None :
                self._value = None
            ids = [id(x) for x in self.choices]
            if id(value) in ids :
                index = ids.index(id(value))
                if self.ui.choices_sizer.GetChildren()[index].GetWindow().IsEnabled() :
                    self._value = id(value)
                else :
                    self._value = None
            else :
                self._value = None
            
            self.ui.layers.Clear()
            if isinstance(self.value, medipy.gui.image.Image) :
                self.ui.layers.Enable()
                for index, layer in enumerate(self.value.layers) :
                    self.ui.layers.Insert(str(1+index), 0)
                    self.ui.layers.SetSelection(self.ui.layers.GetCount()-1)
            else :
                self.ui.layers.Disable()
            self.validate()
            self.notify_observers("value")
    
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
        
        self._update_choices_ui()
    
    def _get_output(self):
        """ The visibility of the "(new)" checkbox.
        """
        
        return self.ui.new.IsShown()
    
    def _set_output(self, value):
        if value and self.may_be_empty :
            raise medipy.base.Exception("output and may_be_empty are mutually exclusive")
        
        self.ui.new.Show(value)
        if value :
            self.ui.new.SetValue(self.output_checked)
        
        self.Layout()
        self.GetParent().Layout()
        self.validate()
    
    def _get_output_checked(self):
        """ Is the "(new)" checkbox checked ?
        """
        return self.ui.new.GetValue()
    
    def _set_output_checked(self, value):
        if isinstance(value, wx.Event) :
            self._set_output_checked(self.ui.new.GetValue())
            value.Skip()
            return
        else :
            self.value = None
            if not value :
                for index, child in enumerate(self.ui.choices_sizer.GetChildren()) :
                    widget = child.GetWindow()
                    if isinstance(widget, wx.RadioButton) and widget.GetValue() :
                        self.value = self.choices[index]
                        break
            self.ui.new.SetValue(value)
            self.ui.choices.Enable(not value)
            self.validate()
    
    def _get_may_be_empty(self):
        """ The visibility of the "(none)" checkbox.
        """
        
        return self.ui.empty.IsShown()
    
    def _set_may_be_empty(self, value):
        if value and self.output :
            raise medipy.base.Exception("output and may_be_empty are mutually exclusive")
        
        self.ui.empty.Show(value)
        if value :
            self.ui.empty.SetValue(self.may_be_empty_checked)
        
        self.Layout()
        self.GetParent().Layout()
        self.validate()
    
    def _get_may_be_empty_checked(self):
        """ Is the "(none)" checkbox checked ?
        """
        
        return self.ui.empty.GetValue()
    
    def _set_may_be_empty_checked(self, value):
        if isinstance(value, wx.Event) :
            self._set_may_be_empty_checked(self.ui.empty.GetValue())
            value.Skip()
            return
        else :
            self.value = None
            if not value :
                for index, child in enumerate(self.ui.choices_sizer.GetChildren()) :
                    widget = child.GetWindow()
                    if isinstance(widget, wx.RadioButton) and widget.GetValue() :
                        self.value = self.choices[index]
                        break
            self.ui.empty.SetValue(value)
            self.ui.choices.Enable(not value)
            self.validate()
    
    def _get_show_layers(self) :
        return self.ui.layers.IsShown()
    
    def _set_show_layers(self, value) :
        self.ui.layers.Show(value)
        self._update_choices_ui()
        self.Layout()
        self.GetParent().Layout()
    
    value = property(_get_value, _set_value)
    choices = property(_get_choices, _set_choices)
    output = property(_get_output, _set_output)
    output_checked = property(_get_output_checked, _set_output_checked)
    may_be_empty = property(_get_may_be_empty, _set_may_be_empty)
    may_be_empty_checked = property(_get_may_be_empty_checked, 
                                    _set_may_be_empty_checked)
    show_layers = property(_get_show_layers, _set_show_layers)
    
    #####################
    # Private interface #
    #####################
    
    def _update_choices_ui(self) :
        self.ui.choices_sizer.Clear(True)
        
        # Shape the sizer to be as square as possible
        rows = max(math.sqrt(len(self.choices)),1)
        rows = math.ceil(rows)
        self.ui.choices_sizer.SetRows(rows)
        self.ui.choices_sizer.SetCols(rows)
        
        if len(self.choices) == 0 :
            # No choice, just display an info message
            label = wx.StaticText(self.ui.choices, label="(no image loaded)")
            font = label.GetFont()
            font.SetStyle(wx.FONTSTYLE_ITALIC)
            label.SetFont(font)
            self.ui.choices_sizer.Add(label, 1, wx.EXPAND)
        else :
            # Display radio buttons
            style=wx.RB_GROUP
            for i in range(len(self.choices)) :
                button = wx.RadioButton(self.ui.choices, -1, str(i+1), style=style)
                style=0
                button.Bind(wx.EVT_RADIOBUTTON, self._set_value)
                self.ui.choices_sizer.Add(button, 0)
                
                if isinstance(self.choices[i], medipy.base.Image) :
                    button.Enable()
                elif self.show_layers and isinstance(self.choices[i], medipy.gui.image.Image) :
                    button.Enable()
                else :
                    button.Disable()
            
            if self._value in [id(x) for x in self.choices] :
                index = [id(x) for x in self.choices].index(self._value)
            else :
                index = 0
            button = self.ui.choices_sizer.GetChildren()[index].GetWindow()
            button.SetValue(True)
        
        self.Layout()
        self.GetParent().Layout()
    
    def _on_choices_modified(self, event):
        self._update_choices_ui()
        if (self.choices and 
                self._value not in [id(x) for x in self.choices] and 
                not self.output_checked and not self.may_be_empty_checked) :
            self.value = self.choices[0]
        elif not self.choices :
            self.value = None
            
        self.validate()
