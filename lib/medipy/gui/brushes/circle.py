##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import numpy

from base import Base

class Circle(Base):
    """ Circle-shaped brush.
    
        This brush is a 2D circle brush, where the circle plane is defined
        by its normal direction.
    """
    
    def __init__(self, value, image, radius, normal):
        self._radius = None
        self._normal = None
        self._indices = None
        
        super(Circle, self).__init__(value, image)
        self._set_radius(radius)
        self._set_normal(normal)
    
    def indices(self, position=None):
        if position is None :
            position = numpy.zeros(self._image.ndim)
        return numpy.add(self._indices, position)
    
    ##############
    # Properties #
    ##############
    
    def _set_image(self, image):
        super(Circle, self)._set_image(image)
        if None not in [self._image, self._radius, self._normal] :
            self._compute_indices()
    
    def _get_radius(self):
        return self._radius
    
    def _set_radius(self, radius):
        self._radius = radius
        if None not in [self._image, self._radius, self._normal] :
            self._compute_indices()
    
    def _get_normal(self):
        """ Normal of the circle plane (must be unit-length) 
        """
        
        return self._normal
    
    def _set_normal(self, normal):
        self._normal = normal
        if None not in [self._image, self._radius, self._normal] :
            self._compute_indices()
    
    radius = property(_get_radius, _set_radius)
    normal = property(_get_normal, _set_normal)

    #####################
    # Private interface #
    #####################
    
    def _compute_indices(self):
        diameter = 2*self._radius+1
        bounding_box = numpy.transpose(numpy.indices(self._image.ndim*[diameter]))
        bounding_box -= self._image.ndim*[self._radius]
        bounding_box = bounding_box.reshape(diameter**self._image.ndim, self._image.ndim)
        bounding_box = bounding_box.astype(float)
        
        self._indices = []
        
        for index in bounding_box :
            distance_to_plane = numpy.abs(numpy.dot(index, self._normal))
            distance_to_origin = numpy.linalg.norm(index)
            if distance_to_plane < 0.5 and distance_to_origin <= self._radius :
                self._indices.append(index.astype(int))
        