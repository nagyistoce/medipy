##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import os

import wx

import medipy.base
from medipy.gui import get_colormap_from_name, colormaps, stage_colormaps
import medipy.gui.base
import medipy.gui.io
import medipy.gui.xrc_wrapper

class LayersPanel(medipy.gui.base.Panel):
    """ Display controls that allow to modify the layers of an image. 
        
        The user may add or delete layers, and change the colors and the 
        window/level of the image.
    """
    
    class UI(medipy.gui.base.UI):
        def __init__(self):
            self.colormaps = None
            self.transparent_background = None
            self.opacity = None
            self.cut_low = None
            self.display_range = None
            self.cut_high = None
            self.layers = None
            self.load = None
            self.delete = None
            self.move_up = None
            self.move_down = None

            self.layer_value = None
            self.layer_slider = None
            self.speed_value = None
            self.speed_slider = None
            self.play = None
            self.stop = None
           
            self.controls = ["colormaps", "transparent_background", "opacity",
                             "cut_low", "display_range", "cut_high", "layers",
                             "load", "delete", "move_up", "move_down","layer_value",
                             "layer_slider","speed_value","speed_slider","play","stop"]
        
        def from_window(self, window, names):
            medipy.gui.base.UI.from_window(self, window, names)
            
            all_colormaps = sorted(colormaps.keys())
            all_colormaps.append(20*"-")
            all_colormaps.extend(sorted(stage_colormaps.keys()))
            self.colormaps.SetItems(all_colormaps)
            
            self.opacity.range = (0,100)
            self.display_range.orientation = wx.VERTICAL
            window.Layout()
    
    def __init__(self, parent=None, *args, **kwargs):
        
        self._image = None
        self._current_layer = 0
        self._timer = wx.Timer() #
        self._delay = None #
        self._status = "stop" #
        # User interface
        self.ui = LayersPanel.UI()
        
        xrc_file = medipy.base.find_resource(
            os.path.join("resources", "gui", "layers_panel.xrc"))
        wrappers = [medipy.gui.xrc_wrapper.CheckListBoxXMLHandler(),
                    medipy.gui.xrc_wrapper.FloatXMLHandler(),
                    medipy.gui.xrc_wrapper.FloatIntervalXMLHandler()]
        medipy.gui.base.Panel.__init__(self, xrc_file, "layers_panel", 
            wrappers, self.ui, self.ui.controls, parent, *args, **kwargs)
        self.ui.layer_slider.SetRange(0,1)
        self.ui.speed_slider.SetRange(1,15)

        # Resources
        self._play_bitmap = wx.Bitmap(medipy.base.find_resource(
            os.path.join("resources", "gui", "media-playback-start.png")))
        self._pause_bitmap = wx.Bitmap(medipy.base.find_resource(
            os.path.join("resources", "gui", "media-playback-pause.png")))
        stop_bitmap = wx.Bitmap(medipy.base.find_resource(
            os.path.join("resources", "gui", "media-playback-stop.png")))
        self.ui.play.SetBitmapLabel(self._play_bitmap)
        self.ui.stop.SetBitmapLabel(stop_bitmap)
        
        # Events
        self.ui.colormaps.Bind(wx.EVT_CHOICE, self.OnColormaps)
        self.ui.display_range.add_observer("value", self.on_display_range)
        self.ui.opacity.add_observer("value", self.on_opacity)
        self.ui.transparent_background.Bind(wx.EVT_CHECKBOX, self.OnTransparentBackground)
        self.ui.cut_low.Bind(wx.EVT_CHECKBOX, self.OnCutLow)
        self.ui.cut_high.Bind(wx.EVT_CHECKBOX, self.OnCutHigh)
        self.ui.layers.Bind(wx.EVT_LISTBOX, self.OnLayers)
        self.ui.layers.Bind(wx.EVT_CHECKLISTBOX, self.OnLayersCheck)
        self.ui.load.Bind(wx.EVT_BUTTON, self.OnLoad)
        self.ui.delete.Bind(wx.EVT_BUTTON, self.OnDelete)
        self.ui.move_up.Bind(wx.EVT_BUTTON, self.OnMoveUp)
        self.ui.move_down.Bind(wx.EVT_BUTTON, self.OnMoveDown)

        self._timer.Bind(wx.EVT_TIMER, self.OnTimer)
        self.ui.layer_value.Bind(wx.EVT_TEXT, self.on_layer_value)
        self.ui.layer_slider.Bind(wx.EVT_SLIDER, self.on_layer_slider)
        self.ui.speed_value.Bind(wx.EVT_TEXT, self.on_speed_value)
        self.ui.speed_slider.Bind(wx.EVT_SLIDER, self.on_speed_slider)
        self.ui.play.Bind(wx.EVT_BUTTON, self.OnPlayPauseClicked)
        self.ui.stop.Bind(wx.EVT_BUTTON, self.OnStopClicked)
        self._set_fps(self.ui.speed_value.GetValue())

    def play(self):
        if self._image is not None and len(self._image.layers)>1 :
            layer = (self._current_layer+1)%len(self._image.layers)
            self.select_layer(layer, True)       
            self._set_fps(self.ui.speed_value.GetValue())
            self._timer.Start(self._delay)
        self._status = "play"
        self.ui.play.SetBitmapLabel(self._pause_bitmap)
    
    def pause(self):
        self._timer.Stop()
        self._status = "pause"
        self.ui.play.SetBitmapLabel(self._play_bitmap)
    
    def stop(self):
        self._timer.Stop()
        if self._image is not None and len(self._image.layers)>1 :
            for i in range(self._image.number_of_layers) :
                self._image.set_layer_visibility(i, False)
            self.select_layer(0, True)
            self._image.render()
        self._status = "stop"
        self.ui.play.SetBitmapLabel(self._play_bitmap)
    
    def select_layer(self, index, select_single_layer=False):
        """ Select the current layer in the layers list, update the GUI.
        
            index must be positive
        """
        
        if self._current_layer is not None:
            image = self._image.get_layer_image(self._current_layer)
            image.remove_observer("modified", self._update_display_range)
        
        self._current_layer = index
        
        self.ui.layers.SetSelection(len(self._image.layers)-index-1)
        
        layer = self._image.layers[self._current_layer]
        colormap = self._image.get_layer_colormap(self._current_layer)
        
        colormap_name = ""
        for name, other_colormap in colormaps.items() :
            if colormap.data == other_colormap :
                colormap_name = name
                break
        if colormap_name == "" :
            for name, other_colormap in stage_colormaps.items() :
                if colormap.data == other_colormap :
                    colormap_name = name
                    break

        if len(self._image.layers)>1 :
            self.ui.layer_slider.SetRange(1,len(self._image.layers))
        if select_single_layer :
            for i in range(len(self._image.layers)) :
                self._image.set_layer_visibility(i, False)
            self._image.set_layer_visibility(self._current_layer, True)
            self._image.render()
        
        image = self._image.get_layer_image(self._current_layer)
        
        self.ui.colormaps.SetStringSelection(colormap_name)
        self.ui.display_range.range = (image.data.min(), image.data.max())
        
        self.ui.display_range.remove_observer("value", self.on_display_range)
        self.ui.display_range.value = colormap.display_range
        self.ui.display_range.add_observer("value", self.on_display_range)
        
        self.ui.opacity.remove_observer("value", self.on_opacity)
        self.ui.opacity.value = layer["opacity"]*100
        self.ui.opacity.add_observer("value", self.on_opacity)
        
        self.ui.transparent_background.SetValue(colormap.zero_transparency)
        self.ui.cut_low.SetValue(colormap.cut_low)
        self.ui.cut_high.SetValue(colormap.cut_high)
        
        self.ui.move_up.Enable(index != len(self._image.layers)-1)
        self.ui.move_down.Enable(index != 0)
        
        image = self._image.get_layer_image(self._current_layer)
        image.add_observer("modified", self._update_display_range)
    
    def insert_layer(self, index, image, colormap=None, opacity=1.0) :
        self._image.insert_layer(index, image, colormap, opacity)
        
        for event in ["data", "display_range", "cut_low", "cut_high", "zero_transparency"] :
            self._image.get_layer_colormap(index).add_observer(event, self.on_colormap_event)
        
        self.ui.delete.Enable(len(self._image.layers)>1)
        
        self._image.render()
    
    def delete_layer(self, index):
        self._image.delete_layer(index)
        
        self.ui.layers.Delete(self.ui.layers.GetCount()-index-1)
        
        self._relabel_layers()
        self._update_checked_layers()
        
        self.ui.delete.Enable(len(self._image.layers)>1)
        
        self._image.render()
    
    def move_layer(self, source_index, destination_index):
        image_1 = self._image.get_layer_image(source_index)
        colormap_1 = self._image.get_layer_colormap(source_index)
        opacity_1 = self._image.get_layer_opacity(source_index)
        display_range_1 = (image_1.data.min(), image_1.data.max())
        label_1 = self.ui.layers.GetString(self.ui.layers.GetCount()-source_index-1)
        
        image_2 = self._image.get_layer_image(destination_index)
        colormap_2 = self._image.get_layer_colormap(destination_index)
        opacity_2 = self._image.get_layer_opacity(destination_index)
        display_range_2 = (image_2.data.min(), image_2.data.max())
        label_2 = self.ui.layers.GetString(self.ui.layers.GetCount()-destination_index-1)
        
        self._image.set_layer_image(source_index, image_2)
        self._image.set_layer_colormap(source_index, colormap_2)
        self._image.set_layer_opacity(source_index, opacity_2)
        self.ui.layers.SetString(self.ui.layers.GetCount()-source_index-1, label_2)
        
        self._image.set_layer_image(destination_index, image_1)
        self._image.set_layer_colormap(destination_index, colormap_1)
        self._image.set_layer_opacity(destination_index, opacity_1)
        self.ui.layers.SetString(self.ui.layers.GetCount()-destination_index-1, label_1)
        
        self._image.render()
        self._update_checked_layers()
    
    ##############
    # Properties #
    ##############
    
    def _get_image(self):
        return self._image
    
    def _set_image(self, image):

        if self._image is not None :
            self._image.layers.remove_observer("any", self._update_layers)
            self._image.remove_observer("layer_visibility", self.on_layer_visibility)
            for i in range(len(self._image.layers)) :
                self._image.get_layer_colormap(i).remove_observer(
                  "display_range", self.on_colormap_event)
        
        self._image = image
        self._update_layers()
        
        if self._image is not None :
            self._image.layers.add_observer("any", self._update_layers)
    
    image = property(_get_image, _set_image)
    
    ##################
    # Event handlers #
    ##################
    
    def OnColormaps(self, dummy):
        """ Called when the colormaps choice changes.
        """
        
        colormap = self._image.get_layer_colormap(self._current_layer)
        colormap_name = self.ui.colormaps.GetStringSelection()
        if colormap_name in colormaps :
            data = colormaps[colormap_name]
        else :
            data = stage_colormaps[colormap_name]
        colormap.data = data
        self._image.render()
    
    def on_display_range(self, dummy):
        """ Called when the display range changes
        """
        
        colormap = self._image.get_layer_colormap(self._current_layer)
        colormap.display_range = self.ui.display_range.value
        self._image.render()
    
    def on_opacity(self, dummy):
        """ Called when the opacity changes
        """
        
        self._image.set_layer_opacity(self._current_layer, 
                                      self.ui.opacity.value/100.)
        self._image.render()
    
    def OnTransparentBackground(self, dummy):
        """ Called when the transparent background flag changes
        """
        
        colormap = self._image.get_layer_colormap(self._current_layer)
        colormap.zero_transparency = self.ui.transparent_background.GetValue()
        self._image.render()
    
    def OnCutLow(self, dummy):
        """ Called when the cut low checkbox is (un)checked.
        """
        
        colormap = self._image.get_layer_colormap(self._current_layer)
        colormap.cut_low = self.ui.cut_low.GetValue()
        self._image.render()
    
    def OnCutHigh(self, dummy):
        """ Called when the cut high checkbox is (un)checked.
        """
        
        colormap = self._image.get_layer_colormap(self._current_layer)
        colormap.cut_high = self.ui.cut_high.GetValue()
        self._image.render()
    
    def OnLayers(self, dummy):
        """ Called when the layers selection is changed
        """
        
        if self.ui.layers.GetSelection() == -1 :
            return
        
        self.select_layer(len(self._image.layers)-self.ui.layers.GetSelection()-1)
    
    def OnLayersCheck(self, event):
        """ Called when a layer is checked or unchecked
        """
        
        wx_index = event.GetInt()
        visibility = self.ui.layers.IsChecked(wx_index)
        index = len(self._image.layers)-wx_index-1
        self._image.set_layer_visibility(index, visibility)
        self._image.render()
    
    def OnLoad(self, dummy):
        """ Called when the load button is clicked.
        """
        
        images = medipy.gui.io.load(self.GetParent(), None)
        
        if images == [] :
            return
        
        self.insert_layer(self._current_layer+1, images[0])
        self.select_layer(self._current_layer+1)
    
    def OnDelete(self, dummy):
        """ Called when the delete button is clicked.
        """
        
        self.delete_layer(self._current_layer)
        self.select_layer(min(self._current_layer, len(self._image.layers)-1))
    
    def OnMoveUp(self, dummy):
        """ Called when the move up button is clicked.
        """
        
        self.move_layer(self._current_layer, self._current_layer+1)
        self.select_layer(self._current_layer+1)
    
    def OnMoveDown(self, dummy):
        """ Called when the move down button is clicked.
        """

        self.move_layer(self._current_layer, self._current_layer-1)
        self.select_layer(self._current_layer-1)
        
    def on_colormap_event(self, event):
        """ Called when a layer's colormap has been changed
        """
        
        if id(event.object) != id(self._image.get_layer_colormap(self._current_layer)) :
            return
        
        colormap = event.object
        if event.event == "display_range" :
            self.ui.display_range.value = colormap.display_range
        elif event.event == "cut_low" :
            self.ui.cut_low.SetValue(colormap.cut_low)
        elif event.event == "cut_high" :
            self.ui.cut_high.SetValue(colormap.cut_high)
        elif event.event == "zero_transparency" :
            self.ui.transparent_background.SetValue(colormap.zero_transparency)
    
    def on_layer_visibility(self, event):
        self.ui.layers.Check(self.ui.layers.GetCount()-event.index-1, 
                             self._image.get_layer_visibility(event.index))
    
    def _relabel_layers(self):
        # Re-label the layers
        for i in range(self.ui.layers.GetCount()) :
            self.ui.layers.SetString(i, "Layer {0}".format(self.ui.layers.GetCount()-i))
    
    def _update_checked_layers(self):
        checked_layers = [self.ui.layers.GetCount()-i-1 
                          for i, _ in enumerate(self._image.layers) 
                          if self._image.get_layer_visibility(i)] 
        self.ui.layers.SetChecked(checked_layers)
    
    def _update_display_range(self, *args):
        """ Update the display_range range of the current image. This function
            can also be used as an event handler.
        """
        
        image = self._image.get_layer_image(self._current_layer)
        self.ui.display_range.range = (image.data.min(), image.data.max())
    
    def _update_layers(self, *args):
        """ Update the layers list in the GUI. This function can also be used as
            an event handler.
        """
        
        self.ui.layers.Clear()
        
        self._current_layer = 0
        
        if self.image is not None :
            for i in range(len(self.image.layers)) :
                for event in ["data", "display_range", "cut_low", "cut_high", "zero_transparency"] :
                    self.image.get_layer_colormap(i).add_observer(event, self.on_colormap_event)
            
            for index, layer in enumerate(self.image.layers) :
                self.ui.layers.Append(
                   "Layer {0}".format(len(self.image.layers)-index))
        
            self._update_checked_layers()
            self.ui.delete.Enable(len(self.image.layers) > 1)
            self.select_layer(self._current_layer)
            
            self.image.add_observer("layer_visibility", self.on_layer_visibility)
        else :
            self.ui.delete.Disable()

    def OnTimer(self, event):
        layer = (self._current_layer+1)%len(self._image.layers)
        self.ui.layer_value.SetValue(str(layer+1))
        self.select_layer(layer, True) 
  
    def OnPlayPauseClicked(self, event):
        if self._status in ["pause", "stop"] :
            self.play()
        else :
            self.pause()
    
    def OnStopClicked(self, event):
        self.stop()

    def on_layer_value(self, event) :
        value = self.ui.layer_value.GetValue()
        if len(value)>0 :
            index = int(value)
            if index<=self.ui.layer_slider.GetMax() :
                self.ui.layer_slider.SetValue(index)
            else :
                index = self.ui.layer_slider.GetMax()
                self.ui.layer_value.SetValue(str(index))
            if self._image is not None :
                self.select_layer(index-1, True)

    def on_layer_slider(self, event):
        self.ui.layer_value.SetValue(str(self.ui.layer_slider.GetValue()))
        if self._image is not None :
            self.select_layer(self.ui.layer_slider.GetValue()-1, True)

    def on_speed_value(self, event) :
        value = self.ui.speed_value.GetValue()
        if len(value)>0 :
            index = int(value)
            if index<=self.ui.speed_slider.GetMax() :
                self.ui.speed_slider.SetValue(index)
            else :
                index = self.ui.speed_slider.GetMax()
                self.ui.speed_value.SetValue(str(index))
            self._set_fps(index)

    def on_speed_slider(self, event):
        self.ui.speed_value.SetValue(str(self.ui.speed_slider.GetValue()))
        self._set_fps(self.ui.speed_slider.GetValue())
    
    #####################
    # Private interface #
    #####################
    
    def _set_fps(self, fps):
        self._delay = 1000. / float(fps)
        if self._timer.IsRunning() :
            self._timer.Start(self._delay)
    
