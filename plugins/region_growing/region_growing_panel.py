##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import copy
import os

import numpy
import wx

import medipy.base
import medipy.gui
import medipy.gui.base
import medipy.gui.xrc_wrapper

from connected_threshold import connected_threshold_with_radius

class RegionGrowingPanel(medipy.gui.base.Panel):
    class UI(medipy.gui.base.UI):
        def __init__(self):
            self.source_layer = None
            self.destination_layer = None
            self.new_destination_layer = None
            self.threshold = None
            self.threshold_type = None
            self.radius = None
            self.accept = None
            self.controls = [
                "source_layer", "destination_layer", "new_destination_layer", 
                "threshold", "threshold_type",
                "radius", 
                "accept"]
            medipy.gui.base.UI.__init__(self)
        
        def from_window(self, window, names):
            medipy.gui.base.UI.from_window(self, window, names)
            self.radius.value = 5
            self.radius.range = [2, 100]
    
    def __init__(self, parent, *args, **kwargs):
        self._image = None
        self._region_layer = None
        
        # Build the panel
        self.ui = RegionGrowingPanel.UI()
        xrc_file = medipy.base.find_resource(
            os.path.join("region_growing", "region_growing.xrc"))
        handlers = [medipy.gui.xrc_wrapper.BoolXMLHandler(),
                    medipy.gui.xrc_wrapper.FloatXMLHandler()]
        medipy.gui.base.Panel.__init__(self, xrc_file, "root", 
            handlers, self.ui, self.ui.controls, parent, *args, **kwargs)
        
        # Event handlers
        self.ui.source_layer.Bind(wx.EVT_CHOICE, self._set_source_layer)
        self.ui.new_destination_layer.add_observer(
            "value", self._on_new_destination_layer)
        self.ui.accept.Bind(wx.EVT_BUTTON, self.accept_region)
        
        self.ui.threshold.add_observer("value", self.update_region)
        self.ui.threshold_type.Bind(wx.EVT_CHOICE, self.update_region)
        self.ui.radius.add_observer("value", self.update_region)
    
    def update_region(self, *args):
        """ Update the segmented region. This function can also be used as an
            observer.
        """
        
        # Check that input parameters are available
        if None in [self.source_layer, self.destination_layer, self.image] :
            return
        if self.ui.threshold_type.GetSelection() == -1 :
            return
        if self._image.cursor_index_position is None :
            return
        
        source_image = self.image.get_layer_image(self.source_layer)
        
        lower = None
        upper = None
        if self.ui.threshold_type.GetStringSelection() == "Above" :
            lower = self.ui.threshold.value
            upper = source_image.data.max()+1
        elif self.ui.threshold_type.GetStringSelection() == "Below" :
            lower = source_image.data.min()-1
            upper = self.ui.threshold.value
        
        roi_image = connected_threshold_with_radius(
            source_image, lower, upper, self.ui.radius.value, 1, 
            [self._image.cursor_index_position])
        
        if self._region_layer is None :
            # Append a layer for the newly-segmented region at the top of other
            # layers
            self._image.layers.remove_observer("any", self._on_image_layers_modified)
            self._region_layer = len(self.image.layers)
            colormap = medipy.gui.Colormap(medipy.gui.colormaps["green"], (0,1),
                                           zero_transparency=True)
            self.image.insert_layer(
                self._region_layer, medipy.base.Image((1,1,1)), colormap)
            self._image.layers.add_observer("any", self._on_image_layers_modified)
        # Update the region image
        self.image.set_layer_image(self._region_layer, roi_image)
        self.image.render() 
    
    def Close(self, force=False):
        self.image = None
        medipy.gui.base.Panel.Close(self, force)
    
    def accept_region(self, *args):
        if self._region_layer is None :
            return
        
        region_image = self.image.get_layer_image(self._region_layer)
        if self.ui.new_destination_layer.value :
            # Delete region_layer to make sure it will be re-created
            self.image.delete_layer(self._region_layer)
            self._region_layer = None
            
            # Create a new layer just below region_layer 
            colormap = medipy.gui.Colormap(medipy.gui.colormaps["rainbow"], None, 
                                           zero_transparency=True)
            self.image.insert_layer(len(self.image.layers), region_image, colormap)
            
            # Un-check the "new" checkbox
            self.ui.new_destination_layer.value = False
        
            # Set the destination layer to the newly-created layer
            self.destination_layer = self.image.number_of_layers-1
            
        max = numpy.max(self._image.get_layer_image(self.destination_layer))
        
        # Update the image
        destination = self.image.get_layer_image(self.destination_layer)
        region = numpy.nonzero(region_image)
        destination[region] = 1+max
        destination.modified()
        self.source_layer = 0
        
        # Update the colormap display range to include the new region 
        self.image.get_layer_colormap(self.destination_layer).display_range = (
            self.image.get_layer_colormap(self.destination_layer).display_range[0],
            float(1+max))
        
        self.image.render()
    
    ##############
    # Properties #
    ##############
    
    def _get_source_layer(self) :
        """ Return the index of the source layer currently selected.
        """
        
        index = self.ui.source_layer.GetCount()-self.ui.source_layer.GetSelection()-1
        return index
    
    def _set_source_layer(self, *args):
        if isinstance(args[0], wx.Event) :
            event = args[0]
            index = len(self.image.layers)-event.GetInt()-1
            self._set_source_layer(index)
        else :
            index = args[0]
            
            # Select the layer in the wx.Choice
            if index < 0 :
                index = len(self.image.layers)+index
                if self._region_layer is not None :
                    # Do not count the region layer
                    index -= 1
            self.ui.source_layer.SetSelection(self.ui.source_layer.GetCount()-index-1)
            
            # Update the threshold widget
            self._update_threshold_range()
            self._update_threshold_value()
    
    def _get_destination_layer(self) :
        """ Return the index of the source layer currently selected.
        """
        
        index = self.ui.source_layer.GetCount()-self.ui.destination_layer.GetSelection()-1
        return index
    
    def _set_destination_layer(self, index):
        if index < 0 :
            index = len(self.image.layers)+index
            if self._region_layer is not None :
                # Do not count the region layer
                index -= 1
        self.ui.destination_layer.SetSelection(self.ui.destination_layer.GetCount()-index-1)
    
    def _get_image(self) :
        return self._image
    
    def _set_image(self, image) :
        # Clean-up the region layers and the observers
        if self._region_layer :
            self._image.delete_layer(self._region_layer)
            self._region_layer = None
            self._image.render()
        if self._image :
            self._image.layers.remove_observer("any", self._on_image_layers_modified)
            self._image.remove_observer("cursor_position", self._on_image_cursor_position)
        
        self._image = image
        if self._image :
            self._image.layers.add_observer("any", self._on_image_layers_modified)
            self._image.add_observer("cursor_position", self._on_image_cursor_position)
        
            self._update_layers_choices()
            # Select the first layer for the source, and the last for destination
            # If source and destination are the same, then check "new destination"
            self.source_layer = 0
            self.destination_layer = 0
            if self.ui.source_layer.GetCount() == 1 :
                self.ui.new_destination_layer.value = True
            else :
                self.destination_layer = -1
                self.ui.new_destination_layer.value = False
            self.Enable()
        else :
            self.Disable()
        
    source_layer = property(_get_source_layer, _set_source_layer)
    destination_layer = property(_get_destination_layer, _set_destination_layer)
    image = property(_get_image, _set_image)
    
    ##################
    # Event handlers #
    ##################
    
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
        
        count = len(self.image.layers)
        if self._region_layer :
            count -= 1
        for index in range(count) :
            label = "Layer {0}".format(count-index)
            self.ui.source_layer.Append(label)
            self.ui.destination_layer.Append(label)
    
    def _update_threshold_range(self) :
        image = self.image.get_layer_image(self.source_layer)
        self.ui.threshold.range = image.data.min(), image.data.max()
    
    def _update_threshold_value(self) :
        index = self.source_layer
        image = self.image.get_layer_image(index)
        
        position = self.image.cursor_index_position
        
        values = []
        for d in range(image.ndim) :
            for offset in [-1, 0, 1] :
                neighbor = copy.copy(position)
                neighbor[d] += offset
                values.append(image[tuple(neighbor)])
        self.ui.threshold.value = float(min(values))
