##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import re

import wx

import medipy.base

class Array(wx.PyPanel, medipy.base.Observable) :
    """ Control to enter an array of numbers (floats or ints).
    """
    
    def __init__(self, parent, type, min_length=None, max_length=None, value=None, 
                 *args, **kwargs):
        
        self._type = type
        self._min_length = min_length
        self._max_length = max_length
        
        self._value = None
        
        wx.PyPanel.__init__(self, parent, *args, **kwargs)
        medipy.base.Observable.__init__(self, ["value"])
        
        self._text = wx.TextCtrl(self)
        
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(self._text, 1, wx.ALIGN_CENTER_VERTICAL)
        self.SetSizer(sizer)
        
        self._text.Bind(wx.EVT_TEXT, self.OnText)
        
        # Set default value
        if value is None :
            value = []
            if self.min_length is not None :
                value = self.min_length*[self.type(0)]
            elif self.max_length is not None :
                value = self.max_length*[self.type(0)]
        self.value = value
    
    def validate(self):
        valid = False
        if isinstance(self._value, (list, tuple)) :
            valid = all([isinstance(x, self.type) for x in self._value])
            if valid and self.min_length is not None : 
                valid = (len(self.value) >= self.min_length)
            if valid and self.max_length is not None :
                valid = (len(self.value) <= self.max_length)
        
        if valid :
            self._text.SetBackgroundColour(None)
        else :
            self._text.SetBackgroundColour(wx.RED)
        
        return valid
    
    ##############
    # Properties #
    ##############
    
    def _get_type(self):
        """ Type of the array elements (may be int or float).
        """
        
        return self._type
    
    def _set_type(self, type):
        self._type = type
    
    def _get_min_length(self):
        """ Minimum length of the array. Set to None to remove constraint.
        """
        
        return self._min_length
    
    def _set_min_length(self, min_length):
        self._min_length = min_length
    
    def _get_max_length(self):
        """ Maximum length of the array. Set to None to remove constraint.
        """
        
        return self._max_length
    
    def _set_max_length(self, max_length):
        self._max_length = max_length
    
    def _get_value(self):
        """ Coordinates stored in the control
        """
        
        return self._value
    
    def _set_value(self, value):
        self._value = value
        
        if self.type == int :
            format = "{0:d}"
        elif self.type == float : 
            format = "{0:.2f}"
        else :
            format = "{0}"
        
        text_value = " ".join([format.format(self.type(x)) for x in reversed(self._value)])
        self._text.SetValue(text_value)
        
        self.notify_observers("value")
        
    type = property(_get_type, _set_type)
    min_length = property(_get_min_length, _set_min_length)
    max_length = property(_get_max_length, _set_max_length)
    value = property(_get_value, _set_value)
    
    ##################
    # Event handlers #
    ##################
    
    def OnText(self, dummy):
        elements = re.split(r"\s+", self._text.GetValue().strip())
        
        try :
            value = [self.type(x) for x in reversed(elements)]
        except ValueError :
            value = None
        self._value = value
        self.validate()
    