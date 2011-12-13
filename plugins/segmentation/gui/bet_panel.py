##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import os

import itk
import numpy
import wx

import medipy.base
import medipy.gui.base
import medipy.gui.colormaps
from medipy.gui.colormap import Colormap
import medipy.gui.xrc_wrapper
import medipy.itk
import medipy.segmentation

class BETPanel(medipy.gui.base.Panel, medipy.base.Observable) :
    """ Panel to allow tweaking of BET parameters.
    """
    
    class UI(medipy.gui.base.UI):
        def __init__(self):
            self.show_cog = None
            self.estimate_cog = None
            self.set_cog_to_cursor = None
            
            self.intensity_range = None
            self.display_intensity_range = None
            self.estimate_intensity_range = None
            
            self.local_threshold = None
            self.smoothness = None
            
            self.controls = [
                "show_cog", "estimate_cog", "set_cog_to_cursor",
                "intensity_range", "display_intensity_range", "estimate_intensity_range",
                "local_threshold", "smoothness"]
    
    def __init__(self, parent=None, *args, **kwargs):
        self.ui = BETPanel.UI()
        
        xrc_file = medipy.base.find_resource(
            os.path.join("segmentation", "resources", "bet_panel.xrc"))
        handlers = [medipy.gui.xrc_wrapper.FloatIntervalXMLHandler(),
                    medipy.gui.xrc_wrapper.FloatXMLHandler(),
                    medipy.gui.xrc_wrapper.IntXMLHandler()]
        medipy.gui.base.Panel.__init__(self, xrc_file, "bet_panel", 
            handlers, self.ui, self.ui.controls,
            parent, *args, **kwargs)
        medipy.base.Observable.__init__(
            self, ["center_of_gravity", "t_2", "t_98", "b_t", "smoothness"])

        self.ui.intensity_range.orientation = "horizontal"
        self.ui.local_threshold.range = (0,1)
        self.ui.local_threshold.value = 0.5
        self.ui.smoothness.range = (0,10)
        self.ui.smoothness.value = 0

        self.ui.show_cog.Bind(wx.EVT_BUTTON, self.OnShowCOG)
        self.ui.estimate_cog.Bind(wx.EVT_BUTTON, self.OnEstimateCOG)
        self.ui.set_cog_to_cursor.Bind(wx.EVT_BUTTON, self.OnSetCOGToCursor)
        
        self.ui.intensity_range.add_observer("value", self.on_intensity_range_value)
        self.ui.display_intensity_range.Bind(wx.EVT_BUTTON, self.OnDisplayIntensityRange)
        self.ui.estimate_intensity_range.Bind(wx.EVT_BUTTON, self.OnEstimateIntensityRange)
        
        self.ui.local_threshold.add_observer("value", self.on_local_threshold)
        self.ui.smoothness.add_observer("value", self.on_smoothness)
        
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        
        self.image = None
        self.layer = None
        self._bet_filter = None
        self._intensity_range_layer = None
    
    def estimate_thresholds(self):
        self._bet_filter.EstimateThresholds()
        self.t_2 = self._bet_filter.GetT2()
        self.t_98 = self._bet_filter.GetT98() 
    
    def show_thresholds(self, show):
        self._image.set_layer_visibility(self._intensity_range_layer, show)
        self._image.render()
        if show :
            self.ui.display_intensity_range.SetLabel("Hide")
        else :
            self.ui.display_intensity_range.SetLabel("Show")
    
    def estimate_center_of_gravity(self):
        self._bet_filter.EstimateCenterOfGravity()
        self.notify_observers("center_of_gravity")
    
    def show_center_of_gravity(self):
        cog = list(reversed(self._bet_filter.GetCenterOfGravity()))
        self._image.cursor_physical_position = cog
        self._image.render()
    
    ##############
    # Properties #
    ##############
    
    def _get_image(self):
        """ Instance of medipy.gui.image.Image to estimate and tweak BET parameters.
        """
        
        return self._image
    
    def _set_image(self, image):
        self._image = image
        self._setup_filter_and_gui()
    
    def _get_layer(self):
        """ Index of the layer to estimate and tweak BET parameters.
        """
        
        return self._layer
    
    def _set_layer(self, layer):
        if self.image :
            if not (0 <= layer <= self.image.number_of_layers-1) :
                raise Exception("Incorrect layer index : {0}".format(layer))
        self._layer = layer
        self._setup_filter_and_gui()
    
    def _get_t_2(self):
        return self._bet_filter.GetT2()
    
    def _set_t_2(self, value):
        self._bet_filter.SetT2(value)
        
        colormap = self._image.get_layer_colormap(self._intensity_range_layer)
        colormap.display_range = (value, colormap.display_range[1])
        
        self.ui.intensity_range.remove_observer("value", self.on_intensity_range_value)
        self.ui.intensity_range.value = (value, self.ui.intensity_range.value[1])
        self.ui.intensity_range.add_observer("value", self.on_intensity_range_value)
        
        self.notify_observers("t_2")
    
    def _get_t_98(self):
        return self._bet_filter.GetT98()
    
    def _set_t_98(self, value):
        self._bet_filter.SetT98(value)
        
        colormap = self._image.get_layer_colormap(self._intensity_range_layer)
        colormap.display_range = (colormap.display_range[0], value)
        
        self.ui.intensity_range.remove_observer("value", self.on_intensity_range_value)
        self.ui.intensity_range.value = (self.ui.intensity_range.value[0], value)
        self.ui.intensity_range.add_observer("value", self.on_intensity_range_value)
        
        self.notify_observers("t_98")
    
    def _get_center_of_gravity(self):
        cog = self._bet_filter.GetCenterOfGravity()
        return list(reversed(cog))
    
    def _set_center_of_gravity(self, value):
        cog = list(reversed(value))
        self._bet_filter.SetCenterOfGravity(cog)
        self.notify_observers("center_of_gravity")
    
    def _get_b_t(self):
        return self._bet_filter.GetBT()
    
    def _set_b_t(self, value):
        self._bet_filter.SetBT(value)
        self.notify_observers("b_t")
    
    def _get_smoothness(self):
        return self._bet_filter.GetSmoothnessFactor()
    
    def _set_smoothness(self, value):
        self._bet_filter.SetSmoothnessFactor(value)
        self.notify_observers("smoothness")
     
    image = property(_get_image, _set_image)
    layer = property(_get_layer, _set_layer)
    t_2 = property(_get_t_2, _set_t_2)
    t_98 = property(_get_t_98, _set_t_98)
    center_of_gravity = property(_get_center_of_gravity, _set_center_of_gravity)
    b_t = property(_get_b_t, _set_b_t)
    smoothness = property(_get_smoothness, _set_smoothness)
    
    ##################
    # Event handlers #
    ##################
    
    def OnShowCOG(self, event):
        event.Skip()
        self.show_center_of_gravity()
    
    def OnEstimateCOG(self, event):
        event.Skip()
        self.estimate_center_of_gravity()
        self.show_center_of_gravity()
    
    def OnSetCOGToCursor(self, event):
        event.Skip()
        self.center_of_gravity = self._image.cursor_physical_position
    
    def on_intensity_range_value(self, dummy):
        if not self.ui.intensity_range.validate() :
            return
        
        self.ui.intensity_range.remove_observer("value", self.on_intensity_range_value)
        t_2, t_98 = self.ui.intensity_range.value
        self.t_2 = t_2
        self.t_98 = t_98
        self.ui.intensity_range.add_observer("value", self.on_intensity_range_value)
    
    def OnDisplayIntensityRange(self, event):
        event.Skip()
        show = (self.ui.display_intensity_range.GetLabel() == "Show")
        self.show_thresholds(show)
        
    def OnEstimateIntensityRange(self, event):
        event.Skip()
        self.estimate_thresholds()
    
    def on_local_threshold(self, dummy):
        if not self.ui.local_threshold.validate() :
            return 
        
        self.b_t = self.ui.local_threshold.value
    
    def on_smoothness(self, dummy):
        if not self.ui.smoothness.validate() :
            return 
        
        self.smoothness = self.ui.smoothness.value
    
    def OnClose(self, event):
        event.Skip()
        if None not in [self.image, self.layer] :
            self.image.delete_layer(self._intensity_range_layer)
            self.image.render()
    
    #####################
    # Private interface #
    #####################
    
    def _setup_filter_and_gui(self):
        """ If both image and layer are specified, create the BET filter, 
            estimate the thresholds and center of gravity, and update the GUI.
        """
        
        enabled = (self.image is not None and self.layer is not None) 
        for control in self.ui.controls :
            getattr(self.ui, control).Enable(enabled)
        
        if enabled :
            image = self._image.get_layer_image(self._layer)
            
            colormap = medipy.gui.Colormap(
                medipy.gui.colormaps["red"], 
                (numpy.NINF, numpy.PINF),
                True, True, True)
            self._image.append_layer(image, colormap, 0.5)
            self._intensity_range_layer = self._image.number_of_layers-1
            self._image.render()
            
            self.ui.intensity_range.range = (image.data.min(), image.data.max())
            self._setup_bet_filter()
            
            self.estimate_thresholds()
            self.show_thresholds(True)
            
            self.estimate_center_of_gravity()
            self.show_center_of_gravity()
            
    def _setup_bet_filter(self):
        medipy_image = self.image.get_layer_image(self.layer)
        self._itk_image = medipy.itk.medipy_image_to_itk_image(medipy_image, False)
        self._bet_filter = itk.BETImageFilter[self._itk_image, self._itk_image].New(
            Input = self._itk_image)
