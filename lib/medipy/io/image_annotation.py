import xml.etree.ElementTree

import medipy.base

def from_xml(root):
    """ Return an ImageAnnotation from an XML document. The root element is
        assumed to be an annotation.
    """
    
    annotation = medipy.base.ImageAnnotation()
    annotation.filled = False
    
    for element in root :
        if element.tag == "position" :
            position = []
            for name in ["z", "y", "x"] :
                if name in element.attrib :
                    position.append(float(element.attrib[name]))
            annotation.position = position
        elif element.tag == "label" :
            annotation.label = element.attrib["value"]
        elif element.tag == "shape" :
            annotation.shape = getattr(medipy.base.ImageAnnotation.Shape, 
                                       element.attrib["value"])
        elif element.tag == "size" :
            annotation.size = float(element.attrib["value"])
        elif element.tag == "color" :
            annotation.color = [float(element.attrib[x]) 
                                for x in ["red", "green", "blue"]]
        elif element.tag == "filled" :
            annotation.filled = True
        elif element.tag == "comment" :
            annotation.comment = element.text
    
    return annotation

def to_xml(annotation):
    """ Convert an ImageAnnotation to an XML element.
    """
    
    root = xml.etree.ElementTree.Element("annotation")
    
    position = xml.etree.ElementTree.Element("position")
    names = ["x", "y", "z"]
    for index, value in enumerate(annotation.position) :
        position.attrib[names[len(annotation.position)-index-1]] = str(value)
    root.append(position)
    
    label = xml.etree.ElementTree.Element("label", value=annotation.label)
    root.append(label)
    
    shape = xml.etree.ElementTree.Element("shape")
    shapes = dict([(getattr(medipy.base.ImageAnnotation.Shape, x), x) 
                   for x in dir(medipy.base.ImageAnnotation.Shape) 
                   if not x.startswith("_")])
    shape.attrib["value"] = shapes[annotation.shape]
    root.append(shape)
    
    size = xml.etree.ElementTree.Element("size", value=str(annotation.size))
    root.append(size)
    
    color = xml.etree.ElementTree.Element("color", 
        red=str(annotation.color[0]), green=str(annotation.color[1]),
        blue=str(annotation.color[2]))
    root.append(color)
    
    if annotation.filled :
        root.append(xml.etree.ElementTree.Element("filled"))
    
    comment = xml.etree.ElementTree.Element("comment")
    comment.text = annotation.comment
    root.append(comment)
    
    return root
