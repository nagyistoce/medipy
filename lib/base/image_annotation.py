##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import logging

from medipy.base import Observable

class ImageAnnotation(Observable):
    """ Annotation localized on an image. Position is in mm, expressed as 
        (z,y,x).
        Can fire the following events : 
            * position : old_value
            * label : old_value
            * shape : old_value
            * size : old_value
            * color : old_value
            * filled : old_value
            * depth : old_value
            * comment : old_value
    """
    
    class Shape(object):
        sphere = 0
        cube = 1
        cross = 2
        point = 3
    
    def __init__(self, position = None, label = None, shape = None, size = None,
        color = None, filled = None, depth = None, comment = None) :
        
        Observable.__init__(self, ["position", "label", "shape",
            "size", "color", "filled", "depth", "comment"])
        
        # Position in image space
        self.position = position if position is not None else [0., 0. ,0.]
        self.label = label or ""
        # Integer, value in ImageAnnotation.Shape
        self.shape = shape or ImageAnnotation.Shape.sphere
        self.size = size or 0.
        # RGB color, each component in 0,1 range
        self.color = color or [0., 0., 0.]
        self.filled = filled or False
        self.depth = depth or 0
        self.comment = comment or ""
    
    def __del__(self):
        logging.debug("medipy.base.ImageAnnotation.__del__(%x)"%id(self))
        
    def __setattr__(self, attr, value):
        if attr in ["position", "label", "shape", "size", "color", "filled", "depth", "comment"] :
            old_value = getattr(self, attr, None) 
            object.__setattr__(self, attr, value)
            self.notify_observers(attr, old_value = old_value)
        else :
            Observable.__setattr__(self, attr, value)