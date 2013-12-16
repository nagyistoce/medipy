##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import logging
import sys

import numpy

import medipy.base

from tag import Tag
from dictionary import data_dictionary, name_dictionary
from private_dictionaries import private_dictionaries
from vr import *
from medipy.io.dicom import dictionary

class DataSet(dict):
    """ An Data Set is a collection (dictionary) of Data Elements values.
        
        The values can be accessed or set using one of the following methods :
            * Numerical value of tag : ``dataset[0x0010, 0x0010]``
            * Name of tag : ``dataset.patients_name``
            * Tag object : ``dataset[medipy.io.dicom.Tag(0x0010, 0x0010)]``
        
        Member requests are checked against the dicom dictionary. If the name 
        matches a DicomDictionary descriptive string, the corresponding tag is 
        used to get or set the Data Element's value.
    """
    
    def __init__(self, include_header=True, **kwargs):
        dict.__init__(self, {})
        if include_header :
            self.header = DataSet(include_header=False)
        for key, value in kwargs.items() :
            self[key] = value
        self.normalized = False
    
    def tags(self):
        """ Return the named tags in the dataset.
        """
        
        result = []
        for key in self.keys() :
            if key.group%2 == 1 :
                result.append(key)
            else :
                result.append(data_dictionary[key][4])
        
        return result
    
    ##############
    # Properties #
    ##############
    
    def _get_pixel_array(self):
        """ Return a numpy array with the pixel data, or None.
        """
        
        if "pixel_data" not in self :
            return None
        
        if "_pixel_array" not in self.__dict__ :
            
            shape = (self.rows.value, self.columns.value)
            if "number_of_frames" in self :
                shape = (self.number_of_frames.value,)+shape
            
            dtype = "int%d"%(self.bits_allocated.value)
            default_pixel_representation = {
                "1.2.840.10008.5.1.4.1.1.2": 1, # CT Image Storage
            }.get(self.get("sop_class_uid", UI(None)).value, 0)
            if self.get("pixel_representation", US(default_pixel_representation)).value == 0 :
                dtype = "u"+dtype
            dtype = numpy.dtype(dtype)
            
            self._pixel_array = numpy.fromstring(self.pixel_data.value, dtype).reshape(shape)
            
            dataset_is_little_endian = (
                self.header.get("transfer_syntax_uid", "") != "1.2.840.10008.1.2.2")
            sys_is_little_endian = (sys.byteorder == 'little')
            
            if sys_is_little_endian != dataset_is_little_endian :
                self._pixel_array.byteswap(True)
            
        return self._pixel_array
        
    pixel_array = property(_get_pixel_array)
    
    #########################
    # Item access functions #
    #########################
    
    def __dir__(self):
        """ Return the list of public elements.
        """
        
        attributes = []
        for tag in self.keys() :
            name = None
            if tag.group%2==0 :
                name = data_dictionary[tag][4]
                
            if name :
                attributes.append(name)
        return attributes
    
    def __getitem__(self, item):
        """ Access to an item using a numerical or named tag.
        """
        
        tag = self._get_tag(item)
        return dict.__getitem__(self, tag)
    
    def __setitem__(self, item, value):
        """ Access to an item using a numerical or named tag.
        """
        
        tag = self._get_tag(item)
        if tag is None :
            raise AttributeError("No such element: {0!r}".format(item))
        
        if not isinstance(value, VR) :
            if tag not in data_dictionary :
                raise medipy.base.Exception("value is not a VR")
            vr = data_dictionary[tag][0]
            value = globals()[vr](value)
        
        dict.__setitem__(self, tag, value)
    
    def __delitem__(self, item) :
        """ Delete an item using a numerical or named tag.
        """
        
        tag = self._get_tag(item)
        dict.__delitem__(self, tag)
    
    ###########################
    # Member access functions #
    ###########################
    
    def __getattr__(self, name):
        """ Access to an item using a named tag. 
        """
        
        # __getattr__ is only called if instance cannot find name in self.__dict__
        
        tag = name_dictionary.get(name, None)
        
        if tag is None :
            raise AttributeError("DataSet does not have attribute '%s'." % name)
        elif tag not in self :
            raise AttributeError("DataSet does not have attribute '%s'." % name)
        else :
            return self[tag]
    
    def __setattr__(self, name, value):
        """ Access to an item using a named tag.
        """
        tag = name_dictionary.get(name, None)
        
        if tag : 
            if not isinstance(value, VR) :
                if tag not in data_dictionary :
                    raise medipy.base.Exception("value is not a VR")
                vr = data_dictionary[tag][0]
                value = globals()[vr](value)
            self[tag] = value
        else:  
            # name is not in DICOM dictionary : set an instance member
            # XXX : no error can be raised if the user misspells a DICOM element
            self.__dict__[name] = value
    
    def __delattr__(self, name):
        """ Delete an item using a named tag.
        """
        
        tag = name_dictionary.get(name, None)
        if tag and tag in self:
            del self[tag]
        elif name in self.__dict__ :
            # Instance member
            del  self.__dict__[name]
        else:
            raise AttributeError("DataSet does not have attribute '%s'." % name)
    
    ########################
    # Comparison functions #
    ########################
    
    def __eq__(self, other):
        if not isinstance(other, DataSet) :
            return False
        if self.keys() != other.keys() :
            return False
        for key, value in self.items() :
            if other[key] != value :
                return False
        return True
    
    def __ne__(self, other):
        return not self.__eq__(other)
    
    ######################
    # Sequence functions #
    ######################
    
    def __contains__(self, key):
        """ Test the presence of numerical or named tag in the Data Set.
        """
        tag = self._get_tag(key)
        return dict.__contains__(self, tag)
    
    ########################
    # Dictionary functions #
    ########################
    
    def get(self, key, default=None):
        """ Access to an item using a numerical or named tag.
        """
        tag = self._get_tag(key)
        return dict.get(self, tag, default)
    
    def setdefault(self, key, value=None):
        """ Access to an item using a numerical or named tag.
        """
        
        if key not in self:
            self[key] = value
        return self[key]
    
    ###################
    # Misc. functions #
    ###################
    
    def __unicode__(self):
        
        result = []
        for tag in sorted(self.keys()):
            value = self[tag]
            
            vr = value.__class__.__name__
            
            if tag.private :
                if tag.element < 0x0100 :
                    name = "Private Creator"
                else :
                    # The private block for private tag (gggg,xxyy) (where gggg
                    # is odd and yy is in [00,ff] is xx. The private creator is
                    # then given in (gggg,00xx)
                    block = tag.element/0x100
                    private_creator = self.get(Tag(tag.group, block), None)
                    
                    if private_creator not in private_dictionaries :
                        private_creator = None
                    
                    private_tag = "{0:04x}xx{1:02x}".format(tag.group, 
                                                            tag.element & 0x00ff)
                    if (private_creator is not None and 
                        private_creator in private_dictionaries and
                        private_tag in private_dictionaries[private_creator]) :
                        
                        name = private_dictionaries[private_creator][private_tag][2]
                    else :
                        name = tag
            else :
                if tag.group/0x100 in [0x50, 0x60] :
                    # Repeating group element, cf. PS 3.5-2011, 7.6
                    tag_in_dictionary = "{0:02x}xx{1:04x}".format(
                        tag.group/0x100, tag.element)
                else :
                    tag_in_dictionary = tag
                name = data_dictionary.get(tag_in_dictionary,
                    ("UN", "1", unicode(tag_in_dictionary)))[2]
            
            if vr == "SQ" :
                elements = [unicode(item) for item in value.value]
                value = []
                for index, element in enumerate(elements) :
                    value.append(4*" "+"Item %i"%index)
                    value.extend([8*" "+unicode(x) for x in element.splitlines()])
                value = "\n"+"\n".join(value)
            elif vr in ["OB", "OW", "OB/OW", "OF", "UN"] :
                value = "<array of %i bytes>"%(len(value.value),)
            elif vr == "UI" and value.value in dictionary.uid_dictionary :
                value = "{0} ({1})".format(dictionary.uid_dictionary[value.value][0],
                                           value.value)
            else :
                value = value.value

            if tag.private :
                try :
                    value = unicode(value)
                except UnicodeDecodeError :
                    value = "<array of %i bytes>"%(len(value),)
            
            result.append("%s %s [%s]: %s"%(name, tag, vr, value))
        return "\n".join(result)
    
    def __str__(self):
        return unicode(self).encode('utf-8')
    
    def __repr__(self):
        return object.__repr__(self)
    
    @staticmethod
    def _get_tag(value):
        """ Create a tag from either a string or a numerical value
        """
        
        if isinstance(value, basestring):
            tag = name_dictionary.get(value, None)
            if tag is not None :
                tag = Tag(tag)
        else :
            tag = Tag(value)
            
        return tag
