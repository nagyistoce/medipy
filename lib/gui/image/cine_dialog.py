##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import os
import xml.dom.minidom

import wx
import wx.xrc

from medipy.base import find_resource
from medipy.gui.control import Int

class CineDialog(wx.Dialog):
    """ Display the dialog to control the cine mode of images
    """
    
    def __init__(self, parent, wx_id, pos=wx.DefaultPosition, size=wx.DefaultSize,
                 style=wx.DEFAULT_DIALOG_STYLE, name="cine_dialog"):
        
        self._timer = wx.Timer()
        self._delay = None
        self._image = None
        self._current_layer = 0
        self._status = "stop"
        
        super(CineDialog, self).__init__(parent, wx_id, "", pos, size, style, name)
        
        # Resources
        self._play_bitmap = wx.Bitmap(find_resource(
            os.path.join("resources", "gui", "media-playback-start.png")))
        self._pause_bitmap = wx.Bitmap(find_resource(
            os.path.join("resources", "gui", "media-playback-pause.png")))
        stop_bitmap = wx.Bitmap(find_resource(
            os.path.join("resources", "gui", "media-playback-stop.png")))
        
        # Widgets        
        self._layer = Int(self, 1, (1,1))
        self._speed = Int(self, 1, (1, 20))
        
        self._play_pause = wx.BitmapButton(self, wx.ID_ANY, self._play_bitmap)
        self._stop = wx.BitmapButton(self, wx.ID_ANY, stop_bitmap)
#        self.close = wx.Button(self, wx.ID_ANY, "Close")
        
        # Layout
        layer_speed_sizer = wx.FlexGridSizer(2, 2)
        layer_speed_sizer.AddGrowableCol(1)
        
        layer_speed_sizer.Add(wx.StaticText(self, wx.ID_ANY, "Layer : "), flag=wx.ALIGN_CENTER)
        layer_speed_sizer.Add(self._layer, 1)
        
        layer_speed_sizer.Add(wx.StaticText(self, wx.ID_ANY, "Speed : "), flag=wx.ALIGN_CENTER)
        layer_speed_sizer.Add(self._speed, 1)
        
        buttons_sizer = wx.BoxSizer(wx.HORIZONTAL)
        buttons_sizer.Add(self._play_pause)
        buttons_sizer.Add(self._stop)
#        buttons_sizer.AddStretchSpacer(1)
#        buttons_sizer.Add(self.close, flag=wx.ALIGN_CENTER)
        
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(layer_speed_sizer, 0, wx.EXPAND)
        main_sizer.Add(buttons_sizer, 0, wx.EXPAND)
        
        self.SetSizer(main_sizer)
        main_sizer.SetSizeHints(self)
        
        # Events
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        self._timer.Bind(wx.EVT_TIMER, self.OnTimer)
        
        self._layer.add_observer("value", self.on_layer)
        self._speed.add_observer("value", self.on_speed)
        self._play_pause.Bind(wx.EVT_BUTTON, self.OnPlayPauseClicked)
        self._stop.Bind(wx.EVT_BUTTON, self.OnStopClicked)
        
        self._set_fps(self._speed.value)
    
    def play(self):
        self._set_fps(self._speed.value)
        self._timer.Start(self._delay)
        self._status = "play"
        self._play_pause.SetBitmapLabel(self._pause_bitmap)
    
    def pause(self):
        self._timer.Stop()
        self._status = "pause"
        self._play_pause.SetBitmapLabel(self._play_bitmap)
    
    def stop(self):
        self._timer.Stop()
        if self._image is not None :
            self.select_layer(0)
            for i in range(self._image.number_of_layers) :
                self._image.set_layer_visibility(i, True)
            self._image.render()
        self._status = "stop"
        self._play_pause.SetBitmapLabel(self._play_bitmap)
    
    def select_layer(self, layer):
        self._current_layer = layer
        if self._image is not None :
            for i in range(len(self._image.layers)) :
                self._image.set_layer_visibility(i, False)
            self._image.set_layer_visibility(self._current_layer, True)
            self._image.render()
        self._layer.value = self._current_layer+1
    
    ##############
    # Properties #
    ##############
    
    def _get_image(self):
        return self._image
    
    def _set_image(self, image):
        # Restore previous image
        if self._image is not None :
            self.stop()
        
        self._image = image
        self._layer.range = (1, self._image.number_of_layers)
        self.select_layer(0)
    
    image = property(lambda x:x._image, _set_image)
    
    ##################
    # Event handlers #
    ##################
    
    def OnClose(self, event):
        if self._image is not None :
            self.stop()
        event.Skip()
    
    def OnTimer(self, event):
        layer = (self._current_layer+1)%len(self._image.layers)
        self.select_layer(layer) 
    
    def OnPlayPauseClicked(self, event):
        if self._status in ["pause", "stop"] :
            self.play()
        else :
            self.pause()
    
    def OnStopClicked(self, event):
        self.stop()
    
    def on_layer(self, event):
        self.select_layer(self._layer.value-1)
    
    def on_speed(self, event):
        self._set_fps(self._speed.value)
    
    #####################
    # Private interface #
    #####################
    
    def _set_fps(self, fps):
        self._delay = 1000. / float(fps)
        if self._timer.IsRunning() :
            self._timer.Start(self._delay)

class OldCineDialog(wx.Dialog):
    """ Display the dialog to control the cine mode of images
    """
    
    def __init__(self, parent, image, *args, **kwargs):
        
        self._timer = wx.Timer()
        self._delay = None
        self._image = None
        self._layer = 0
        
        wx.Dialog.__init__(self, parent, *args, **kwargs)
        
        dialog_file = find_resource(os.path.join("resources", "gui", "cine_dialog.xrc"))
        self._dialog_xml = xml.dom.minidom.parse(dialog_file)
        
        self._dialog_xrc = wx.xrc.XmlResource("")
        self._dialog_xrc.LoadFromString(str(self._dialog_xml.toxml()))
        
        # Widgets
        self._panel = self._dialog_xrc.LoadPanel(self, "main_panel")
        
        controls = ["layer_slider", "layer_text", "speed_slider", "speed_text",
             "play_pause_button", "stop_button", "close_button"
        ]
        for control in controls :
            setattr(self, "_"+control, wx.xrc.XRCCTRL(self, control))
        
        
        self._speed_slider.SetValue(5)
        self._speed_text.SetValue(str(self._speed_slider.GetValue()))
        
        self._play_bitmap = wx.Bitmap(find_resource(os.path.join("resources", "gui", "media-playback-start.png")))
        self._pause_bitmap = wx.Bitmap(find_resource(os.path.join("resources", "gui", "media-playback-pause.png")))
        stop_bitmap = wx.Bitmap(find_resource(os.path.join("resources", "gui", "media-playback-stop.png")))
        
        self._play_pause_button.SetBitmapLabel(self._play_bitmap)
        
        self._stop_button.SetBitmapLabel(stop_bitmap)
        

        # Layout
        global_sizer = wx.BoxSizer(wx.VERTICAL)
        global_sizer.Add(self._panel, 1, wx.EXPAND)

        self.SetTitle("Cine Image")
        self.SetSizer(global_sizer)
        global_sizer.Layout()
        self.Fit()
        
        # Events
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        self._timer.Bind(wx.EVT_TIMER, self.OnTimer)
        
        self._layer_slider.Bind(wx.EVT_SLIDER, self.OnLayerSlider)
        self._layer_text.Bind(wx.EVT_TEXT, self.OnLayerText)
        
        self._speed_slider.Bind(wx.EVT_SLIDER, self.OnSpeedSlider)
        self._speed_text.Bind(wx.EVT_TEXT, self.OnSpeedText)
        
        self._play_pause_button.Bind(wx.EVT_BUTTON, self.OnPlayPauseButtonClick)
        self._stop_button.Bind(wx.EVT_BUTTON, self.OnStopButtonClick)
        self._close_button.Bind(wx.EVT_BUTTON, lambda event : self.Close())
        
        
        self._set_image(image)
        self._set_fps(self._speed_slider.GetValue())
        self._set_layer(0)
    
    ####################
    # Public interface #
    ####################
    def play(self) :
        self._set_fps(self._speed_slider.GetValue())
        self._timer.Start(self._delay)
        
        self._play_pause_button.SetBitmapLabel(self._pause_bitmap)
    
    def pause(self):
        self._timer.Stop()
        
        self._play_pause_button.SetBitmapLabel(self._play_bitmap)
    
    def stop(self):
        self._timer.Stop()
        self._set_layer(0)
        for i in range(self._image.number_of_layers) :
            self._image.show_layer(i)
        self._image.render()
        
        self._play_pause_button.SetBitmapLabel(self._play_bitmap)
    
    ##############
    # Properties #
    ##############
    
    def _set_image(self, image):
        # Restore previous image
        if self._image is not None :
            self.stop()
        
        self._image = image
        if self._image.number_of_layers > 1 :
            self._layer_slider.SetMax(self._image.number_of_layers-1)
            self._layer_slider.Enable()
            self._set_layer(0)
        else :
            self._layer_slider.Disable() 
    
    image = property(lambda x:x._image, _set_image)
    
    #####################
    # Private interface #
    #####################
    
    def _set_fps(self, fps):
        self._delay = 1000. / float(fps)
        if self._timer.IsRunning() :
            self._timer.Start(self._delay)
    
    def _set_layer(self, layer):
        self._layer = layer
        for i in range(self._image.number_of_layers) :
            self._image.hide_layer(i)
        self._image.show_layer(self._layer)
        self._image.render()
        self._layer_slider.SetValue(self._layer)
        self._layer_text.ChangeValue(str(self._layer+1))
        self._layer_text.SetInsertionPointEnd()
    
    ##################
    # Event handlers #
    ##################
    
    def OnPlayPauseButtonClick(self, event):
        if self._play_pause_button.GetBitmapLabel() == self._play_bitmap :
            self.play()
        else :
            self.pause()
    
    def OnStopButtonClick(self, event):
        self.stop()
    
    def OnSpeedSlider(self, event):
        self._speed_text.SetValue(str(self._speed_slider.GetValue()))
        self._set_fps(self._speed_slider.GetValue())
    
    def OnSpeedText(self, event):
        try :
            value = int(self._speed_text.GetValue())
        except :
            self._speed_text.SetBackgroundColour(wx.RED)
        else :
            if self._speed_slider.GetMin() <= value <= self._speed_slider.GetMax() : 
                self._speed_slider.SetValue(value)
                self._speed_text.SetBackgroundColour(None)
                self._set_fps(self._speed_slider.GetValue())
            else :  
                self._speed_text.SetBackgroundColour(wx.RED)
    
    def OnClose(self, event):
        self.stop()
        event.Skip()
    
    def OnTimer(self, event):
        layer = (self._layer+1)%self._image.number_of_layers
        self._set_layer(layer) 
    
    def OnLayerSlider(self, event):
        self._set_layer(self._layer_slider.GetValue())
    
    def OnLayerText(self, event):
        try :
            value = int(self._layer_text.GetValue())-1
        except :
            self._layer_text.SetBackgroundColour(wx.RED)
        else :
            if self._layer_slider.GetMin() <= value <= self._layer_slider.GetMax() : 
                self._layer_text.SetBackgroundColour(None)
                self._set_layer(value)
            else :  
                self._layer_text.SetBackgroundColour(wx.RED)
