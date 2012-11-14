##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import numpy
from vtk import vtkActor, vtkContourFilter, vtkPolyDataMapper

from layer import Layer

class ContourLayer(Layer):
    "Layer showing its data as a set of isocontours."
    
    @staticmethod
    def can_create(image):
        return (image.data_type == "scalar" and
                image.image_type == "spectroscopy" and
                image.ndim <= 3)
    
    def __init__(self, world_to_slice, image, display_coordinates="physical",
                 colormap=None, opacity = 1.0, levels=None) :
        
        ############################
        # Property-related members #
        ############################
        
        # List of levels, or None to use default-generated values
        self._levels = None 
        
        self._actor = vtkActor()
        
        ###################
        # Private members #
        ###################
        self._contour_filter = vtkContourFilter()
        self._mapper = vtkPolyDataMapper()
        
        ##################
        # Initialization #
        ##################
        
        super(ContourLayer, self).__init__(
            world_to_slice, image, display_coordinates, colormap, opacity)
        
        self._contour_filter.SetInput(self._change_information.GetOutput())
        self._contour_filter.UseScalarTreeOn()
        self._mapper.SetInputConnection(self._contour_filter.GetOutputPort())
        self._mapper.ScalarVisibilityOn()
        self._mapper.SetLookupTable(self._colormap.vtk_colormap)
        self._mapper.UseLookupTableScalarRangeOn()
        self._actor.SetMapper(self._mapper)
        
        if levels is None :
            # Generate 20 evenly spaced contour values
            levels = numpy.linspace(self._image.data.min(), 
                                    self._image.data.max(), 20)
        self._set_levels(levels)
    
    ##############
    # Properties #
    ##############
    
    def _set_image(self, image):
        super(ContourLayer, self)._set_image(image)
        self._mapper.SetScalarRange(image.data.min(), image.data.max())
    
    def _set_colormap(self, colormap):
        super(ContourLayer, self)._set_colormap(colormap)
        self._mapper.SetLookupTable(self._colormap.vtk_colormap)
        self.colormap.add_observer("vtk_colormap", self._on_vtk_colormap)
    
    def _set_opacity(self, opacity):
        super(ContourLayer, self)._set_opacity(opacity)
        self._actor.GetProperty().SetOpacity(opacity)
    
    def _get_levels(self) :
        "Levels at which contours are computed."
        return self._levels
    
    def _set_levels(self, levels):
        self._levels = levels
        for i, level in enumerate(self._levels) :
            self._contour_filter.SetValue(i, level)
        
    def _get_actor(self):
        return self._actor

    levels = property(_get_levels, _set_levels)
    
    #####################
    # Private interface #
    #####################
    
    def _on_vtk_colormap(self, dummy) :
        """ Event handler called when the ID of colormap.vtk_colormap changes.
        """
        
        self._set_colormap(self._colormap)

Layer.derived_classes.append(ContourLayer)
