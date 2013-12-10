##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import medipy.base
from medipy.base import PropertySynchronized
import medipy.vtk
from medipy.vtk import vtkEnhancedLookupTable, vtkColorTransferFunctionWithAlpha

class Colormap(PropertySynchronized) :
    """ A colormap, mapping a scalar value to a color. Can either correspond to
        a vtkEnhancedLookupTable or a vtkColorTransferFunctionWithAlpha.
        
        The properties this class can be synchronized on are :
          * display_range
          * cut_low
          * cut_high
          * zero_transparency
        
        TODO : explain what parameters do to scalar -> color mapping
    """
    
    def __init__(self, data, display_range, 
                 cut_low=False, cut_high=False, zero_transparency=False) :
        
        self._data = None
        self._vtk_colormap = None
        self._display_range = None
        self._cut_low = None
        self._cut_high = None
        self._zero_transparency = None
        
        super(Colormap, self).__init__(
            ["data", "display_range", "cut_low", "cut_high", "zero_transparency"])
        # Event triggered when the id of the vtk colormap has changed, and that
        # the pipeline should be rebuilt
        self.add_allowed_event("vtk_colormap")
        
        self._set_data(data)
        self._set_display_range(display_range)
        self._set_cut_low(cut_low)
        self._set_cut_high(cut_high)
        self._set_zero_transparency(zero_transparency)
    
    ##############
    # Properties #
    ##############
    
    def _get_data(self) :
        "Data of the colormap"
        return self._data
    
    def _set_data(self, data) :
        self._data = data
        
        previous_id = id(self._vtk_colormap)
        self._vtk_colormap = medipy.vtk.build_vtk_colormap(self._data)
        if id(self._vtk_colormap) != previous_id :
            self.notify_observers("vtk_colormap")
        self._update_vtk_colormap()
        
        self.notify_observers("data")
    
    def _get_vtk_colormap(self) :
        "VTK object representing the colormap"
        return self._vtk_colormap
        
    def _get_display_range(self) :
        "Interval of displayed values"
        return self._display_range
        
    def _set_display_range(self, display_range) :
        self._display_range = display_range
        self._update_vtk_colormap()
        self.notify_observers("display_range")
        
    def _get_cut_low(self) :
        "Mapping of image values below the display range"
        return self._cut_low
        
    def _set_cut_low(self, cut_low) :
        self._cut_low = cut_low
        self._update_vtk_colormap()
        self.notify_observers("cut_low")
        
    def _get_cut_high(self) :
        "Mapping of image values above the display range"
        return self._cut_high
        
    def _set_cut_high(self, cut_high) :
        self._cut_high = cut_high
        self._update_vtk_colormap()
        self.notify_observers("cut_high")
        
    def _get_zero_transparency(self) :
        "Determines if the value 0 is mapped to a fully transparent color or not"
        return self._zero_transparency
        
    def _set_zero_transparency(self, zero_transparency) :
        self._zero_transparency = zero_transparency
        self._update_vtk_colormap()
        self.notify_observers("zero_transparency")
        
    data = property(_get_data, _set_data)
    vtk_colormap = property(_get_vtk_colormap)
    display_range = property(_get_display_range, _set_display_range)
    cut_low = property(_get_cut_low, _set_cut_low)
    cut_high = property(_get_cut_high, _set_cut_high)
    zero_transparency = property(_get_zero_transparency, _set_zero_transparency)

    #####################
    # Private interface #
    #####################
    def _update_vtk_colormap(self):
        if None in [self._data, self._display_range, self._cut_low, 
                    self._cut_high, self._zero_transparency] :
            return 
        
        if isinstance(self._vtk_colormap, vtkEnhancedLookupTable) :
            self._vtk_colormap.SetCutLow(self._cut_low)
            self._vtk_colormap.SetCutHigh(self._cut_high)
            self._vtk_colormap.SetTableRange(*self._display_range)
            self._vtk_colormap.SetZeroTransparency(self._zero_transparency)
        elif isinstance(self._vtk_colormap, vtkColorTransferFunctionWithAlpha) :
            
            # Determine the color of the min and max values of the display range
            min_color = self._data[0][1]
            max_color = self._data[-1][1]
            for (a,b), color in self._data :
                if a <= self._display_range[0] <= b :
                    min_color = color
                if a <= self._display_range[1] <= b :
                    max_color = color
            
            data = []
            
            for (a,b), color in self._data :
                is_transparent = ((self._cut_low and a < self._display_range[0]) 
                                  or
                                  (self._cut_high and b > self._display_range[1]))
                
                if is_transparent :
                    data.append(((a,b), (0,0,0,0)))
                else :
                    # Color is not transparent, clip to min or max if necessary
                    if a < self._display_range[0] :
                        data.append(((a,b), min_color))
                    elif b > self._display_range[1] :
                        data.append(((a,b), max_color))
                    else :
                        data.append(((a,b), color))
            
            self._vtk_colormap = medipy.vtk.build_vtk_colormap(data, self._vtk_colormap)
        else :
            raise medipy.base.Exception("Unknown colormap type : %s"%self._vtk_colormap.GetClassName())
