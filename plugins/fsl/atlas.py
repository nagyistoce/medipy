##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import os
from xml.etree.ElementTree import XMLParser

import medipy.base

class Atlas(object):
    """ Atlas from FSL with the following attributes :
    
        * ``name``(e.g. ``"Juelich Histological Atlas"``)
        * ``type``, either ``label`` (each voxel has a definite class) or
          ``probabilistic`` (each voxel has a list of probabilities of 
          belonging to a class)
        * ``images`` : a list of pair of filenames. For ``label`` atlases, the 
          two elements are the same, and correspond to the label image. For
          probabilistic atlases, the first element is the 4D image containing
          the probabilities for each class, and the second element is the label
          image of the maximal probability class in each voxel.
        * ``labels`` : a mapping of labels to their names
        * ``centers`` : a mapping of labels to their centers in the image.
    """
    
    Type = medipy.base.enum("Type", "label", "probabilistic")
    
    def __init__(self) : 
        self.name = None
        self.type = None
        self.images = []
        self.labels = {}
        self.centers = {}
        
    @staticmethod
    def read(filename):
        """ Read an atlas from a XML file.
        """
        
        builder = TreeBuilder(filename)
        parser = XMLParser(target=builder)
        
        data = open(filename).read()
        parser.feed(data)
        return parser.close()
    
class TreeBuilder(object):
    """ XML tree builder for the FSL atlas format. 
    """
    
    def __init__(self, filename):
        self._atlas = Atlas()
        self._filename = filename
        self._state = None
        
        self._image = None
        self._summary_image = None
        
        self._index = None
        self._label = None
    
    def start(self, tag, attributes):
        self._state = tag
        
        if tag == "atlas" :
            if "version" not in attributes :
                raise medipy.base.Exception("No version specified")
            if attributes["version"] != "1.0" :
                raise medipy.base.Exception("Unknown version {0}".format(attributes["version"]))
        elif tag == "label" :
            if "index" not in attributes :
                raise medipy.base.Exception("Attribute \"index\" missing from \"label\" element")
            try :
                self._index = int(attributes["index"])
            except ValueError :
                raise medipy.base.Exception("Cannot parse \"index\" attribute with value {0}".format(repr(attributes["index"])))
            
            center = (int(attributes.get("z", 0)),
                      int(attributes.get("y", 0)),
                      int(attributes.get("x", 0)))
            self._atlas.centers[self._index] = center
    
    def end(self, tag):
        if tag == "images" :
            self._atlas.images.append((self._image, self._summary_image))
        elif tag == "label" :
            self._atlas.labels[self._index] = self._label
        self._state = None
    
    def data(self, data):
        if self._state == "name" :
            self._atlas.name = data
        elif self._state == "type" :
            # The typo in "Probabalistic" is a "feature" of FSL
            types = { "Label" : Atlas.Type.label, 
                      "Probabalistic" : Atlas.Type.probabilistic}
            if data not in types.keys() :
                raise medipy.base.Exception("Unknown type {0}".format(repr(data)))
            self._atlas.type = types[data]
        elif self._state == "imagefile" :
            if data.startswith("/") :
                data = data[1:]
            root = os.path.join(os.path.dirname(self._filename), data)
            
            candidates = ["{0}.nii".format(root), "{0}.nii.gz".format(root)]
            image = None
            for candidate in candidates :
                if os.path.isfile(candidate) :
                    image = candidate
                    break
            
            if image is None :
                raise medipy.base.Exception("Cannot find image {0}".format(repr(root)))
            self._image = image
        elif self._state == "summaryimagefile" :
            if data.startswith("/") :
                data = data[1:]
            root = os.path.join(os.path.dirname(self._filename), data)
            
            candidates = ["{0}.nii".format(root), "{0}.nii.gz".format(root)]
            image = None
            for candidate in candidates :
                if os.path.isfile(candidate) :
                    image = candidate
                    break
            
            if image is None :
                raise medipy.base.Exception("Cannot find summary image {0}".format(repr(root)))
            self._summary_image = image
        elif self._state == "label" :
            self._label = data
    
    def close(self):
        return self._atlas
