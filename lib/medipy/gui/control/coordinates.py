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

class Coordinates(wx.PyPanel, medipy.base.Observable) :
    """ Control to enter the coordinates of a point in an image.
    """
    
    def __init__(self, parent, value=None, image=None, display_coordinates="index",
                 *args, **kwargs):
        self._value = None
        self._image = None
        self._display_coordinates = None
        
        wx.PyPanel.__init__(self, parent, *args, **kwargs)
        medipy.base.Observable.__init__(self, ["value"])
        
        self._text = wx.TextCtrl(self)
        self._show = wx.Button(self, label="Show")
        self._set = wx.Button(self, label="Set")
        
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(self._text, 1, wx.ALIGN_CENTER_VERTICAL)
        sizer.Add(self._show, flag=wx.ALIGN_CENTER_VERTICAL)
        sizer.Add(self._set, flag=wx.ALIGN_CENTER_VERTICAL)
        self.SetSizer(sizer)
        
        self._text.Bind(wx.EVT_TEXT, self.OnText)
        self._show.Bind(wx.EVT_BUTTON, self.OnShow)
        self._set.Bind(wx.EVT_BUTTON, self.OnSet)
        
        self.display_coordinates = display_coordinates
        self.value = value or (0,0,0)
        self.image = image
    
    def validate(self):
        if self._display_coordinates == "index" :
            type = int
        elif self._display_coordinates == "physical" :
            type = float
        
        if isinstance(self._value, (list, tuple)) :
            valid = all([isinstance(x, type) for x in self._value])
        else :
            valid = False
        
        if valid :
            self._text.SetBackgroundColour(None)
        else :
            self._text.SetBackgroundColour(wx.RED)
        
        return valid
    
    ##############
    # Properties #
    ##############
    
    def _get_value(self):
        """ Coordinates stored in the control
        """
        
        return self._value
    
    def _set_value(self, value):
        self._value = value
        
        if self.display_coordinates == "index" :
            format = "{0:d}"
            cast = int
        elif self.display_coordinates == "physical" : 
            format = "{0:.2f}"
            cast = float
        
        text_value = " ".join([format.format(cast(x)) for x in reversed(self._value)])
        self._text.SetValue(text_value)
        
        self.notify_observers("value")
        
    def _get_image(self):
        """ An instance of medipy.gui.image.Image used to pick the coordinates
        """
        
        return self._image
    
    def _set_image(self, image):
        self._image = image
        
        show_buttons = (image is not None)
        self._show.Show(show_buttons)
        self._set.Show(show_buttons)
        self.Layout()
    
    def _get_display_coordinates(self) :
        return self._display_coordinates
    
    def _set_display_coordinates(self, display_coordinates):
        if display_coordinates not in ["physical", "index"] :
            raise medipy.base.Exception("Unknown display coordinates : %s"%(display_coordinates,))
        
        self._display_coordinates = display_coordinates
    
    value = property(_get_value, _set_value)
    image = property(_get_image, _set_image)
    display_coordinates = property(_get_display_coordinates, _set_display_coordinates)
    
    ##################
    # Event handlers #
    ##################
    
    def OnText(self, dummy):
        if self._display_coordinates == "index" :
            cast = int
        elif self._display_coordinates == "physical" :
            cast = float
        
        elements = re.split(r"\s+", self._text.GetValue().strip())
        
        try :
            value = [cast(x) for x in reversed(elements)]
        except ValueError :
            value = None
        self._value = value
        self.validate()
    
    def OnShow(self, dummy):
        name = "cursor_{0}_position".format(self.display_coordinates)
        setattr(self.image, name, self.value)
    
    def OnSet(self, dummy):
        name = "cursor_{0}_position".format(self.display_coordinates)
        self.value = getattr(self.image, name)

if __name__ == "__main__" :
    import sys
    import medipy.io
    import medipy.gui.image
    
    image = medipy.io.load(sys.argv[1])
    
    app = wx.PySimpleApp()
    frame = wx.Frame(None)
    viewer = medipy.gui.image.Image(frame, layers=[{"image":image}])
    viewer.SetMinSize((400,400))
    coordinates = Coordinates(frame, image=viewer, display_coordinates="physical")
    
    sizer = wx.BoxSizer(wx.VERTICAL)
    sizer.Add(viewer, 1, wx.EXPAND)
    sizer.Add(coordinates, 1, wx.EXPAND)
    frame.SetSizer(sizer)
    sizer.SetSizeHints(frame)
    
    frame.Show()
    app.MainLoop()