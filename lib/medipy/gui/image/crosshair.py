##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import numpy
from vtk import vtkActor, vtkAssembly, vtkLineSource, vtkPolyDataMapper
import medipy.base

class Crosshair(object):
    """ A crosshair, to be displayed on a :class:`~medipy.gui.image.Slice`.
    """
    
    def __init__(self, altitude=0, extent=None, position=None, mode="full",
                 color=None, hole_size=1):
        self._altitude = None
        self._extent = None
        self._position = None
        self._mode = None
        self._color = None
        self._hole_size = None
        
        # TODO : color
        
        self._sources = {}
        self._actor = vtkAssembly()
        
        self.altitude = altitude
        self.extent = extent or [[-1,-1], [1,1]]
        self.position = position or [0,0]
        self.mode = mode
        self.color = color or (1,0,0)
        self.hole_size = hole_size
    
    ##############
    # Properties #
    ##############
    
    def _get_altitude(self):
        """ Altitude of the crosshair.
        """
        
        return self._altitude
    
    def _set_altitude(self, altitude):
        self._altitude = altitude
        self._update_lines()
    
    def _get_extent(self):
        """ Extent of the crosshair, given in TODO coordinates. The extent
            will always be converted to a numpy array.
        """
        
        return self._extent
    
    def _set_extent(self, extent):
        self._extent = numpy.asarray(extent)
        self._update_lines()
    
    def _get_position(self):
        """ Position of the center of the crosshair, in VTK world coordinates,
            VTK order. The position will always be converted to a numpy array.
        """
        
        return self._position
    
    def _set_position(self, position):
        self._position = numpy.asarray(position)
        self._update_lines()
    
    def _get_mode(self):
        """ Display mode of the cursor
        
            * ``"full"`` : the full horizontal and vertical lines are displayed,
              based on the extent
            * ``"partial"`` : a hole is left at the center of the horizontal and
              vertical lines. The size of the hole is given by 
              :attr:`hole_size`.
            * ``"none"`` : no lines are displayed.
        """
        
        return self._mode
    
    def _set_mode(self, mode):
        if mode not in ["full", "partial", "none"] :
            raise medipy.base.Exception("Unknown crosshair mode: {0!r}".format(value))
        
        self._mode = mode
        
        # Clean-up the assembly
        parts = self._actor.GetParts()
        parts.InitTraversal()
        for _ in range(parts.GetNumberOfItems()) :
            part = parts.GetNextProp3D()
            self._actor.RemovePart(part)
        self._sources = {}
        
        # Create the sources
        if self.mode == "full" :
            self._sources = {
                "horizontal" : vtkLineSource(), 
                "vertical" : vtkLineSource()
            }
        elif self.mode == "partial" :
            self._sources = {
                "horizontal_before" : vtkLineSource(), 
                "horizontal_after" : vtkLineSource(),
                "vertical_before" : vtkLineSource(),
                "vertical_after" : vtkLineSource()
            }
        
        # Create the mappers and actors, add parts to the assembly
        mappers = []
        for source in self._sources.values() :
            mapper = vtkPolyDataMapper()
            mapper.SetInputConnection(source.GetOutputPort())
            actor = vtkActor()
            actor.SetMapper(mapper)
            
            if self.color is not None :
                actor.GetProperty().SetColor(self.color)
            
            self._actor.AddPart(actor)
        
        self._update_lines()
            
    def _get_color(self):
        """ Color of the crosshair in RGB, where each component lies between 0
            and 1. The color will always be converted to a numpy array.
        """
        
        return self._color
    
    def _set_color(self, color):
        self._color = numpy.asarray(color)
        
        parts = self._actor.GetParts()
        parts.InitTraversal()
        for _ in range(parts.GetNumberOfItems()) :
            part = parts.GetNextProp3D()
            part.GetProperty().SetColor(self.color)
    
    def _get_hole_size(self):
        """ Size of the hole at the center of the lines in ``"partial"`` mode.
            This size is in VTK world units.
        """
        
        return self._hole_size
    
    def _set_hole_size(self, hole_size):
        self._hole_size = hole_size
        
        if self._mode == "partial" :
            self._update_lines()
    
    def _get_actor(self):
        """ Actor for the crosshair
        """
        return self._actor
    
    altitude  = property(_get_altitude, _set_altitude)
    extent = property(_get_extent, _set_extent)
    position = property(_get_position, _set_position)
    mode = property(_get_mode, _set_mode)
    color = property(_get_color, _set_color)
    hole_size = property(_get_hole_size, _set_hole_size)
    actor = property(_get_actor)
    
    #####################
    # Private interface #
    #####################
    
    def _update_lines(self):
        """ Update the line sources based on altitude, extent, position, mode
            and hole_size (if mode == "partial").
        """
        
        if None in [self.altitude, self.extent, self.position, self.mode] :
            # Object is not fully initialized
            return
        if self.mode == "partial" and self.hole_size is None :
            # Object is not fully initialized
            return
        
        if self._mode == "full" :
            self._sources["horizontal"].SetPoint1(
                self.extent[0][0], self.position[1], self.altitude)
            self._sources["horizontal"].SetPoint2(
                self.extent[1][0], self.position[1], self.altitude)
            
            self._sources["vertical"].SetPoint1(
                self.position[0], self.extent[0][1], self.altitude)
            self._sources["vertical"].SetPoint2(
                self.position[0], self.extent[1][1], self.altitude)
        elif self._mode == "partial" :
            self._sources["horizontal_before"].SetPoint1(
                self.extent[0][0], self.position[1], self.altitude)
            self._sources["horizontal_before"].SetPoint2(
                self.position[0]-self.hole_size, self.position[1], self.altitude)
            
            self._sources["horizontal_after"].SetPoint1(
                self.position[0]+self.hole_size, self.position[1], self.altitude)
            self._sources["horizontal_after"].SetPoint2(
                self.extent[1][0], self.position[1], self.altitude)
            
            self._sources["vertical_before"].SetPoint1(
                self.position[0], self.extent[0][1], self.altitude)
            self._sources["vertical_before"].SetPoint2(
                self.position[0], self.position[1]-self.hole_size, self.altitude)
            
            self._sources["vertical_after"].SetPoint1(
                self.position[0], self.position[1]+self.hole_size, self.altitude)
            self._sources["vertical_after"].SetPoint2(
                self.position[0], self.extent[1][1], self.altitude)
        # Otherwise do nothing, as mode is "none"
