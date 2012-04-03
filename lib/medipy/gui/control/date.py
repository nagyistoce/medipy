##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import datetime

import wx

import medipy.base

class Date(wx.PyPanel, medipy.base.Observable) :
    """ Control to enter a date.
    """
    
    def __init__(self, parent, value=None, *args, **kwargs):
        self._value = None
        
        wx.PyPanel.__init__(self, parent, *args, **kwargs)
        medipy.base.Observable.__init__(self, ["value"])
        
        self._date = wx.DatePickerCtrl(self, style=wx.DP_DEFAULT|wx.DP_ALLOWNONE)
        
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(self._date, 1, wx.EXPAND)
        
        self._date.Bind(wx.EVT_DATE_CHANGED, self.OnDate)
        
        #self.value = value
    
    def validate(self):
        return self._date.GetValue().IsValid()
    
    ##############
    # Properties #
    ##############
    
    def _get_value(self):
        """ Date stored in the control
        """
        
        return self._value
    
    def _set_value(self, value):
        self._value = value
        
        time_tuple = value.timetuple()
        
        wx_date = wx.DateTimeFromDMY(
            time_tuple[2], time_tuple[1]-1, time_tuple[0])
        
        self._date.SetValue(wx_date)
        
        self.notify_observers("value")
    
    value = property(_get_value, _set_value)
    
    ##################
    # Event handlers #
    ##################
    
    def OnDate(self, event) :
        wx_date = event.GetDate()
        
        if wx_date.IsValid() :
            ymd = map(int, wx_date.FormatISODate().split('-'))
            value = datetime.date(*ymd)
        else :
            value = None
        
        self.value = value
