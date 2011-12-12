##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import medipy.base

class Base(object):
    """ Base class for all brushes
    """
    
    def __init__(self, value, image):
        self._value = None
        self._image = None
        
        self._set_value(value)
        self._set_image(image)
    
    def indices(self, position=None):
        """ Return a list of indices that will be affected when painting at the
            given position. If position is not specified, it is assumed to be
            the origin.
        """
        
        raise NotImplementedError()
    
    def paint(self, position):
        """ Set the value of the voxels in self.indices(position) to self.value. 
        """
        
        self._image[self.indices(position)] = self._value
    
    ##############
    # Properties #
    ##############
    
    def _get_value(self):
        """ Value the affected voxels will be set to.
        """
        
        return self._value
    
    def _set_value(self, value):
        self._value = value
    
    def _get_image(self):
        """ Image to draw on.
        """
        
        return self._image
    
    def _set_image(self, image):
        self._image = image
    
    value = medipy.base.LateBindingProperty(_get_value, _set_value)
    image = medipy.base.LateBindingProperty(_get_image, _set_image)
