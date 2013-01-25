##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import numpy

from base import Base

class Sphere(Base):
    """ Sphere-shaped brush.
    
        This brush is a nD sphere brush.
    """
    
    def __init__(self, value, radius, dimension=3):
        self._radius = None
        self._dimension = None
        self._indices = None
        
        super(Sphere, self).__init__(value)
        self._set_radius(radius)
        self._set_dimension(dimension)
    
    def indices(self, position=None):
        """ See :func:`medipy.gui.brushes.Base.indices`.
        """
        
        if position is None :
            position = numpy.zeros(self._dimension)
        return numpy.add(self._indices, position)
    
    ##############
    # Properties #
    ##############
    
    def _get_radius(self):
        """ Radius in pixels.
        """
        
        return self._radius
    
    def _set_radius(self, radius):
        self._radius = radius
        if None not in [self._radius, self._dimension] :
            self._compute_indices()
    
    def _get_dimension(self):
        """ Number of dimensions.
        """
        
        return self._dimension
    
    def _set_dimension(self, dimension):
        self._dimension = dimension
        if None not in [self._radius, self._dimension] :
            self._compute_indices()
    
    radius = property(_get_radius, _set_radius)
    dimension = property(_get_dimension, _set_dimension)

    #####################
    # Private interface #
    #####################
    
    def _compute_indices(self):
        diameter = 2*self._radius+1
        
        first = numpy.asarray(self._dimension*[-self._radius], float)
        last = numpy.asarray(self._dimension*[self._radius], float)
        
        bounding_box = numpy.transpose(numpy.indices(last-first+1))+first
        bounding_box = bounding_box.reshape(diameter**self._dimension, self._dimension)
        
        self._indices = [x.astype(int) for x in bounding_box 
                         if numpy.linalg.norm(x)<=self._radius]
        