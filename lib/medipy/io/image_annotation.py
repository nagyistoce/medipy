##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import xml.etree.ElementTree
import numpy
import medipy.base

def annotation_from_xml(root):
    """ Return an ImageAnnotation from an XML document.
    """
    
    if root.tag != "ImageAnnotation" :
        raise medipy.base.Exception("Cannot build an ImageAnnotation from a {0}".format(repr(root.tag)))
    
    annotation = medipy.base.ImageAnnotation()
    
    child_tags = dict([(x.tag, x) for x in root.getchildren() if x.text is not None])
    if "position" in child_tags :
        annotation.position = numpy.asarray([
            float(x) for x in child_tags["position"].text.split(" ")])
    if "label" in child_tags :
        annotation.label = child_tags["label"].text.decode("utf-8")
    if "shape" in child_tags :
        annotation.shape = getattr(medipy.base.ImageAnnotation.Shape, 
                                   child_tags["shape"].text)
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

def annotations_from_xml(root):
    """ Return a list of ImageAnnotation from an XML document.
    """
    
    if root.tag != "ImageAnnotations" :
        raise medipy.base.Exception("Cannot build a list of ImageAnnotations from a {0}".format(repr(root.tag)))
    
    annotations = []
    
    for element in root:
        annotations.append(annotation_from_xml(element))
    
    return annotations

def annotation_to_xml(annotation):
    """ Convert an ImageAnnotation to an XML element.
    """
    
    root = xml.etree.ElementTree.Element("ImageAnnotation")
        
    position = xml.etree.ElementTree.SubElement(root, "position")
    position.text = " ".join([str(x) for x in annotation.position])
    
    label = xml.etree.ElementTree.SubElement(root, "label")
    label.text = annotation.label.encode("utf-8")
    
    shape = xml.etree.ElementTree.SubElement(root, "shape")
    shape.text = medipy.base.ImageAnnotation.Shape.to_name(annotation.shape)
    
    size = xml.etree.ElementTree.SubElement(root, "size")
    size.text = str(annotation.size)
    
    color = xml.etree.ElementTree.SubElement(root, "color")
    color.text = " ".join([str(x) for x in annotation.color])
    
    filled = xml.etree.ElementTree.SubElement(root, "filled")
    filled.text = str(annotation.filled)
    
    comment = xml.etree.ElementTree.SubElement(root, "comment")
    comment.text = annotation.comment.encode("utf-8")
    
    return root

def annotations_to_xml(annotations):
    """ Convert a list of ImageAnnotation to an XML element.
    """
    
    root = xml.etree.ElementTree.Element("ImageAnnotations")
    for annotation in annotations :
        element = annotation_to_xml(annotation)
        root.append(element)
    
    return root