##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import logging
import xml.etree.ElementTree

from exception import Exception
from observable import Observable

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
            * comment : old_value
    """
    
    class Shape(object):
        sphere = 0
        cube = 1
        cross = 2
        point = 3
        
        @staticmethod
        def to_name(value):
            dictionary = dict([(getattr(ImageAnnotation.Shape, name), name) 
                               for name in dir(ImageAnnotation.Shape) 
                               if isinstance(getattr(ImageAnnotation.Shape, name), int)])
            return dictionary[value]
            
    
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
    
    def to_xml(self):
        root = xml.etree.ElementTree.Element("ImageAnnotation")
        
        position = xml.etree.ElementTree.SubElement(root, "position")
        position.text = " ".join([str(x) for x in self.position])
        
        label = xml.etree.ElementTree.SubElement(root, "label")
        label.text = self.label.encode("utf-8")
        
        shape = xml.etree.ElementTree.SubElement(root, "shape")
        shape.text = ImageAnnotation.Shape.to_name(self.shape)
        
        size = xml.etree.ElementTree.SubElement(root, "size")
        size.text = str(self.size)
        
        color = xml.etree.ElementTree.SubElement(root, "color")
        color.text = " ".join([str(x) for x in self.color])
        
        filled = xml.etree.ElementTree.SubElement(root, "filled")
        filled.text = str(self.filled)
        
        comment = xml.etree.ElementTree.SubElement(root, "comment")
        comment.text = self.comment.encode("utf-8")
        
        return root
    
    @staticmethod
    def from_xml(root):
        if root.tag != "ImageAnnotation" :
            raise Exception("Cannot build an ImageAnnotation from a {0}".format(repr(root.tag)))
        
        annotation = ImageAnnotation()
        
        child_tags = dict([(x.tag,x) for x in root.getchildren()])
        if "position" in child_tags :
            annotation.position = [
                float(x) for x in child_tags["position"].text.split(" ")]
        if "label" in child_tags :
            annotation.label = child_tags["label"].text.decode("utf-8")
        if "shape" in child_tags :
            annotation.shape = getattr(ImageAnnotation.Shape, child_tags["shape"].text)
        if "size" in child_tags :
            annotation.size = float(child_tags["size"].text)
        if "color" in child_tags :
            annotation.color = [
                float(x) for x in child_tags["color"].text.split(" ")]
        if "filled" in child_tags :
            annotation.filled = (child_tags["filled"].text.lower() == "true")
        if "comment" in child_tags :
            annotation.comment = child_tags["comment"].text.decode("utf-8")
        
        return annotation
        