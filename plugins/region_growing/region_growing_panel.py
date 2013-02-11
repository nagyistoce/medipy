##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import copy
import os

import wx

import medipy.base
import medipy.gui.base
import medipy.gui.xrc_wrapper

class RegionGrowingPanel(medipy.gui.base.Panel):
    class UI(medipy.gui.base.UI):
        def __init__(self):
            self.source_layer = None
            self.destination_layer = None
            self.new_destination_layer = None
            self.threshold = None
            self.threshold_type = None
            self.radius = None
            self.controls = [
                "source_layer", "destination_layer", "new_destination_layer", 
                "threshold", "threshold_type",
                "radius"]
            medipy.gui.base.UI.__init__(self)
        
        def from_window(self, window, names):
            medipy.gui.base.UI.from_window(self, window, names)
            self.radius.value = 5
            self.radius.range = [2, 100]
    
    def __init__(self, parent, *args, **kwargs):
        self._image = None
        
        # Build the panel
        self.ui = RegionGrowingPanel.UI()
        xrc_file = medipy.base.find_resource(
            os.path.join("region_growing", "region_growing.xrc"))
        handlers = [medipy.gui.xrc_wrapper.BoolXMLHandler(),
                    medipy.gui.xrc_wrapper.FloatXMLHandler()]
        medipy.gui.base.Panel.__init__(self, xrc_file, "root", 
            handlers, self.ui, self.ui.controls, parent, *args, **kwargs)
        
        # Event handlers
        self.ui.source_layer.Bind(wx.EVT_CHOICE, self.OnSourceLayer)
        self.ui.new_destination_layer.add_observer(
            "value", self._on_new_destination_layer)
    
    ##############
    # Properties #
    ##############
    
    def _get_source_layer(self) :
        """ Return the index of the source layer currently selected.
        """
        
        index = len(self.image.layers)-self.ui.source_layer.GetSelection()-1
        return index
    
    def _get_destination_layer(self) :
        """ Return the index of the source layer currently selected.
        """
        
        index = len(self.image.layers)-self.ui.destination_layer.GetSelection()-1
        return index
    
    def _get_image(self) :
        return self._image
    
    def _set_image(self, image) :
        self._image = image
        self._image.layers.add_observer("any", self._on_image_layers_modified)
        self._image.add_observer("cursor_position", self._on_image_cursor_position)
        
        self._update_layers_choices()
        # Select the first layer for the source, and the last for destination
        # If source and destination are the same, then check "new destination"
        self.ui.source_layer.Select(self.ui.source_layer.GetCount()-1)
        self.ui.destination_layer.Select(0)
        if len(self.image.layers) == 1 :
            self.ui.new_destination_layer.value = True
        
        self._update_threshold_range()
        self._update_threshold_value()
        
    source_layer = property(_get_source_layer)
    destination_layer = property(_get_destination_layer)
    image = property(_get_image, _set_image)
    
    ##################
    # Event handlers #
    ##################
    
    def OnSourceLayer(self, dummy) :
        self._update_threshold_range()
    
    def _on_new_destination_layer(self, dummy) :
        self.ui.destination_layer.Enable(not self.ui.new_destination_layer.value)
    
    def _on_image_layers_modified(self, dummy) :
        self._update_layers_choices()
    
    def _on_image_cursor_position(self, dummy) :
        self._update_threshold_value()
    
    #####################
    # Private interface #
    #####################
    
    def _update_layers_choices(self) :
        self.ui.source_layer.Clear()
        self.ui.destination_layer.Clear()
        for index, layer in enumerate(self.image.layers) :
            label = "Layer {0}".format(len(self.image.layers)-index)
            self.ui.source_layer.Append(label)
            self.ui.destination_layer.Append(label)
    
    def _update_threshold_range(self) :
        index = self.ui.source_layer.GetSelection()
        index = self.ui.source_layer.GetCount()-index-1
        image = self.image.get_layer_image(index)
        self.ui.threshold.range = image.data.min(), image.data.max()
    
    def _update_threshold_value(self) :
        index = self.ui.source_layer.GetSelection()
        index = self.ui.source_layer.GetCount()-index-1
        image = self.image.get_layer_image(index)
        
        position = self.image.cursor_index_position
        
        values = []
        for d in range(image.ndim) :
            for offset in [-1, 0, 1] :
                neighbor = copy.copy(position)
                neighbor[d] += offset
                values.append(image[tuple(neighbor)])
        self.ui.threshold.value = min(values)
