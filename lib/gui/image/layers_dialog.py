##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import wx

from medipy.gui import get_colormap_from_name, colormaps, stage_colormaps
from medipy.gui.control import Float, FloatInterval
from medipy.gui.io import load

class LayersDialog(wx.Dialog):
    """ Display controls that allow to modify the layers of an image. 
        
        The user may add or delete layers, and change the colors and the 
        window/level of the image.
    """
    
    def __init__(self, parent, wx_id, pos=wx.DefaultPosition, size=wx.DefaultSize,
                 style=wx.DEFAULT_DIALOG_STYLE, name="layers_dialog"):
        
        super(LayersDialog, self).__init__(parent, wx_id, "", pos, size, style, name)
        
        self._image = None
        self._current_layer = 0
        self._display_ranges = []
        
        # Widgets
        self.layers = wx.CheckListBox(self, wx.ID_ANY)
        self.display_range = FloatInterval(self, wx.ID_ANY, vertical=False)
        self.opacity = Float(self, 0, (0, 100)) #wx.Slider(self, wx.ID_ANY, 0, 0, 100)
        
        self.colormaps = wx.Choice(self, wx.ID_ANY)
        all_colormaps = sorted(colormaps.keys())
        all_colormaps.append(20*"-")
        all_colormaps.extend(sorted(stage_colormaps.keys()))
        self.colormaps.SetItems(all_colormaps)
        
        self.transparent_background = wx.CheckBox(self, wx.ID_ANY, "Transparent background")
        self.cut_low = wx.CheckBox(self, wx.ID_ANY, "Cut low")
        self.cut_high = wx.CheckBox(self, wx.ID_ANY, "Cut high")
        self.load = wx.BitmapButton(self, wx.ID_ANY, 
                                    wx.ArtProvider.GetBitmap(wx.ART_FILE_OPEN))
        self.delete = wx.BitmapButton(self, wx.ID_ANY, 
                                      wx.ArtProvider.GetBitmap(wx.ART_CROSS_MARK))
        self.move_up = wx.BitmapButton(self, wx.ID_ANY, 
                                       wx.ArtProvider.GetBitmap(wx.ART_GO_UP))
        self.move_down = wx.BitmapButton(self, wx.ID_ANY, 
                                         wx.ArtProvider.GetBitmap(wx.ART_GO_DOWN))
        
        # Layout
        colormap_sizer = wx.BoxSizer(wx.HORIZONTAL)
        colormap_sizer.Add(wx.StaticText(self, label="Colormap : "), 
                           flag=wx.ALIGN_CENTER | wx.LEFT, border=5)
        colormap_sizer.Add(self.colormaps, flag=wx.RIGHT, border=5)
        colormap_sizer.Add(self.transparent_background, flag=wx.RIGHT | wx.ALIGN_CENTER, border=5)
        
        opacity_sizer = wx.BoxSizer(wx.HORIZONTAL)
        opacity_sizer.Add(wx.StaticText(self, label="Opacity : "), 
                          flag=wx.ALIGN_CENTER|wx.LEFT, border=5)
        opacity_sizer.Add(self.opacity, 1, wx.EXPAND | wx.RIGHT, border=5)
        
        display_range_sizer = wx.BoxSizer(wx.HORIZONTAL)
        display_range_sizer.Add(self.cut_low, 
                                flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT, border=5)
        display_range_sizer.Add(self.display_range, 1, wx.EXPAND | wx.RIGHT, border=5)
        display_range_sizer.Add(self.cut_high,
                                flag=wx.ALIGN_CENTER | wx.RIGHT, border=5)
        
        buttons_sizer = wx.BoxSizer(wx.HORIZONTAL)
        buttons_sizer.Add(self.load, flag=wx.LEFT|wx.RIGHT, border=5)
        buttons_sizer.Add(self.delete, flag=wx.RIGHT, border=5)
        buttons_sizer.Add(self.move_up, flag=wx.RIGHT, border=5)
        buttons_sizer.Add(self.move_down, flag=wx.RIGHT, border=5)
        
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(colormap_sizer, flag=wx.EXPAND)
        main_sizer.Add(opacity_sizer, flag=wx.EXPAND)
        main_sizer.Add(display_range_sizer, flag=wx.EXPAND)
        main_sizer.Add(self.layers, 1, wx.EXPAND)
        main_sizer.Add(buttons_sizer, flag=wx.EXPAND)
        
        self.SetSizer(main_sizer)
        main_sizer.SetSizeHints(self)
        
        # Events
        self.colormaps.Bind(wx.EVT_CHOICE, self.OnColormaps)
        self.display_range.add_observer("value", self.on_display_range)
        self.opacity.add_observer("value", self.on_opacity)
        self.transparent_background.Bind(wx.EVT_CHECKBOX, self.OnTransparentBackground)
        self.cut_low.Bind(wx.EVT_CHECKBOX, self.OnCutLow)
        self.cut_high.Bind(wx.EVT_CHECKBOX, self.OnCutHigh)
        self.layers.Bind(wx.EVT_LISTBOX, self.OnLayers)
        self.layers.Bind(wx.EVT_CHECKLISTBOX, self.OnLayersCheck)
        self.load.Bind(wx.EVT_BUTTON, self.OnLoad)
        self.delete.Bind(wx.EVT_BUTTON, self.OnDelete)
        self.move_up.Bind(wx.EVT_BUTTON, self.OnMoveUp)
        self.move_down.Bind(wx.EVT_BUTTON, self.OnMoveDown)
    
    def select_layer(self, index):
        """ Select the current layer in the layers list, update the GUI.
        
            index must be positive
        """
        
        self._current_layer = index
        
        self.layers.SetSelection(len(self._image.layers)-index-1)
        
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
        
        self.colormaps.SetStringSelection(colormap_name)
        self.display_range.range = self._display_ranges[self._current_layer]
        
        self.display_range.remove_observer("value", self.on_display_range)
        self.display_range.value = colormap.display_range
        self.display_range.add_observer("value", self.on_display_range)
        
        self.opacity.remove_observer("value", self.on_opacity)
        self.opacity.value = layer["opacity"]*100
        self.opacity.add_observer("value", self.on_opacity)
        
        self.transparent_background.SetValue(colormap.zero_transparency)
        self.cut_low.SetValue(colormap.cut_low)
        self.cut_high.SetValue(colormap.cut_high)
        
        self.move_up.Enable(index != len(self._image.layers)-1)
        self.move_down.Enable(index != 0)
    
    def insert_layer(self, index, image, colormap=None, opacity=1.0) :
        self._image.insert_layer(index, image, colormap, opacity)
        self._display_ranges.insert(index, (image.data.min(), image.data.max()))
        
        for event in ["data", "display_range", "cut_low", "cut_high", "zero_transparency"] :
            self._image.get_layer_colormap(index).add_observer(event, self.on_colormap_event)
        
        if self.layers.GetCount() == 1 or index == self.layers.GetCount():
            self.layers.Insert("", 0)
        else :
            self.layers.Insert("", self.layers.GetCount()-index-1)
        
        self._relabel_layers()
        self._update_checked_layers()
        
        self.delete.Enable(len(self._image.layers)>1)
        
        self._image.render()
    
    def delete_layer(self, index):
        self._image.delete_layer(index)
        del self._display_ranges[index]
        
        self.layers.Delete(self.layers.GetCount()-index-1)
        
        self._relabel_layers()
        self._update_checked_layers()
        
        self.delete.Enable(len(self._image.layers)>1)
        
        self._image.render()
    
    def move_layer(self, source_index, destination_index):
        image_1 = self._image.get_layer_image(source_index)
        colormap_1 = self._image.get_layer_colormap(source_index)
        opacity_1 = self._image.get_layer_opacity(source_index)
        display_range_1 = self._display_ranges[source_index]
        label_1 = self.layers.GetString(self.layers.GetCount()-source_index-1)
        
        image_2 = self._image.get_layer_image(destination_index)
        colormap_2 = self._image.get_layer_colormap(destination_index)
        opacity_2 = self._image.get_layer_opacity(destination_index)
        display_range_2 = self._display_ranges[destination_index]
        label_2 = self.layers.GetString(self.layers.GetCount()-destination_index-1)
        
        self._image.set_layer_image(source_index, image_2)
        self._image.set_layer_colormap(source_index, colormap_2)
        self._image.set_layer_opacity(source_index, opacity_2)
        self._display_ranges[source_index] = display_range_2
        self.layers.SetString(self.layers.GetCount()-source_index-1, label_2)
        
        self._image.set_layer_image(destination_index, image_1)
        self._image.set_layer_colormap(destination_index, colormap_1)
        self._image.set_layer_opacity(destination_index, opacity_1)
        self._display_ranges[destination_index] = display_range_1
        self.layers.SetString(self.layers.GetCount()-destination_index-1, label_1)
        
        self._image.render()
        self._update_checked_layers()
    
    ##############
    # Properties #
    ##############
    
    def _get_image(self):
        return self._image
    
    def _set_image(self, image):

        if self._image is not None :
            for i in range(len(self._image.layers)) :
                self._image.get_layer_colormap(i).remove_observer(
                  "display_range", self.on_colormap_event)
        self._image = image
        self._current_layer = 0
        
        self._display_ranges = []
        for i in range(len(self._image.layers)) :
            self._display_ranges.append((image.get_layer_image(i).data.min(), 
                                         image.get_layer_image(i).data.max()))
            for event in ["data", "display_range", "cut_low", "cut_high", "zero_transparency"] :
                self._image.get_layer_colormap(i).add_observer(event, self.on_colormap_event)
        
        self.layers.Clear()
        for index, layer in enumerate(self._image.layers) :
            self.layers.Append(
               "Layer {0}".format(len(self._image.layers)-index))
        
        self._update_checked_layers()
        self.delete.Enable(len(self._image.layers) > 1)
        self.select_layer(self._current_layer)
    
    image = property(_get_image, _set_image)
    
    ##################
    # Event handlers #
    ##################
    
    def OnColormaps(self, dummy):
        """ Called when the colormaps choice changes.
        """
        
        colormap = self._image.get_layer_colormap(self._current_layer)
        colormap_name = self.colormaps.GetStringSelection()
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
        colormap.display_range = self.display_range.value
        self._image.render()
    
    def on_opacity(self, dummy):
        """ Called when the opacity changes
        """
        
        self._image.set_layer_opacity(self._current_layer, 
                                      self.opacity.value/100.)
        self._image.render()
    
    def OnTransparentBackground(self, dummy):
        """ Called when the transparent background flag changes
        """
        
        colormap = self._image.get_layer_colormap(self._current_layer)
        colormap.zero_transparency = self.transparent_background.GetValue()
        self._image.render()
    
    def OnCutLow(self, dummy):
        """ Called when the cut low checkbox is (un)checked.
        """
        
        colormap = self._image.get_layer_colormap(self._current_layer)
        colormap.cut_low = self.cut_low.GetValue()
        self._image.render()
    
    def OnCutHigh(self, dummy):
        """ Called when the cut high checkbox is (un)checked.
        """
        
        colormap = self._image.get_layer_colormap(self._current_layer)
        colormap.cut_high = self.cut_high.GetValue()
        self._image.render()
    
    def OnLayers(self, dummy):
        """ Called when the layers selection is changed
        """
        
        if self.layers.GetSelection() == -1 :
            return
        
        self.select_layer(len(self._image.layers)-self.layers.GetSelection()-1)
    
    def OnLayersCheck(self, dummy):
        """ Called when a layer is checked or unchecked
        """
        
        for index in range(len(self._image.layers)) :
            visibility = self.layers.IsChecked(len(self._image.layers)-index-1)
            self._image.set_layer_visibility(index, visibility)
        
        self._image.render()
    
    def OnLoad(self, dummy):
        """ Called when the load button is clicked.
        """
        
        images = load(self.GetParent(), None)
        
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
            self.display_range.value = colormap.display_range
        elif event.event == "cut_low" :
            self.cut_low.SetValue(colormap.cut_low)
        elif event.event == "cut_high" :
            self.cut_high.SetValue(colormap.cut_high)
        elif event.event == "zero_transparency" :
            self.transparent_background.SetValue(colormap.zero_transparency)
    
    def _relabel_layers(self):
        # Re-label the layers
        for i in range(self.layers.GetCount()) :
            self.layers.SetString(i, "Layer {0}".format(self.layers.GetCount()-i))
    
    def _update_checked_layers(self):
        checked_layers = [self.layers.GetCount()-i-1 
                          for i, _ in enumerate(self._image.layers) 
                          if self._image.get_layer_visibility(i)] 
        self.layers.SetChecked(checked_layers)