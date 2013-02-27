##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import logging
import xml.etree.ElementTree

from enum import enum
from exception import Exception
from observable import Observable

class ImageAnnotation(Observable):
    """ Annotation localized on an image. Position is in mm, expressed as 
        (z,y,x). Can fire the following events : 
            
            * position : old_value
            * label : old_value
            * shape : old_value
            * size : old_value
            * color : old_value
            * filled : old_value
            * comment : old_value
    """
    
    Shape = enum("Shape", "sphere", "cube", "cross", "point")
    
    def __init__(self, position = None, label = None, shape = None, size = None,
        color = None, filled = None, comment = None) :
        
        Observable.__init__(self, ["position", "label", "shape",
            "size", "color", "filled", "comment"])
        
        # Position in image space
        self.position = position if position is not None else [0., 0. ,0.]
        self.label = label or ""
        # Integer, value in ImageAnnotation.Shape
        self.shape = shape or ImageAnnotation.Shape.sphere
        self.size = size or 0.
        # RGB color, each component in 0,1 range
        self.color = color or [0., 0., 0.]
        self.filled = filled or False
        self.comment = comment or ""
    
    def __setattr__(self, attr, value):
        if attr in ["position", "label", "shape", "size", "color", "filled", "comment"] :
            old_value = getattr(self, attr, None) 
            object.__setattr__(self, attr, value)
            self.notify_observers(attr, old_value = old_value)
        else :
            Observable.__setattr__(self, attr, value)
