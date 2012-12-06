##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import logging
import sys

import numpy

from tag import Tag
from dictionary import data_dictionary, name_dictionary
from private_dictionaries import private_dictionaries
from medipy.io.dicom import dictionary

class DataSet(dict):
    """ An Data Set is a collection (dictionary) of Data Elements values.
        
        Example of two ways to retrieve or set values:
        
        1. dataset[0x10, 0x10] --> patient's name
        2. dataset.patients_name --> patient's name
        
        Example (2) is referred to as *Named tags* in this documentation.
        patients_name is not actually a member of the object, but unknown member
        requests are checked against the dicom dictionary. If the name matches a
        DicomDictionary descriptive string, the corresponding tag is used
        to look up or set the Data Element's value.
    """
    
    @staticmethod
    def from_dict(dictionary, tags = "all", private_tags = "all", process_header="first_run"):
        """ Create an Data Set from a dictionary. tags can be either a
            sequence of tags (numerical or named) from the dictionary to be
            included in the Data Set or the string "all". In the
            latter case, all tags from the dictionary will be included in the
            Data Set. private_tags has a similar effect for private
            tags.
        """
        
        dataset = DataSet()
        
        if process_header == "first_run" :
            header_dict = dict([(tag, element) 
                                for tag, element in dictionary.items() 
                                if tag[0]==0x0002])
            if header_dict :
                header = DataSet.from_dict(header_dict, tags, private_tags, "header_only")
                dataset.header = header
        
        for tag, element in dictionary.items() :
            if process_header != "header_only" and tag[0] == 0x0002 :
                continue
            if tag[1] == 0 :
                # group length
                continue
            
            # Find VR and whether to include the element
            if tag[0]%2==1 :
                include_tag = (private_tags == "all" or
                               (private_tags and tag in private_tags))
            else :
                include_tag = (tags == "all" or
                               tag in tags or tag_to_underscored_name(tag) in tags)
            
            if include_tag :
                if isinstance(element, list) and element and isinstance(element[0], dict) :
                    value = []
                    for item in element :
                        sub_dataset = DataSet.from_dict(item, tags, private_tags, "no_header")
                        value.append(sub_dataset)
                else :
                    value = element
                dict.__setitem__(dataset, Tag(tag), value)
        
        return dataset
    
    def __init__(self):
        dict.__init__(self, {})
        self.header = {}
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
            
            shape = (self.rows, self.columns)
            if "number_of_frames" in self :
                shape = (self.number_of_frames,)+shape
            
            dtype = "int%d"%(self.bits_allocated)
            if self.get("pixel_representation", 0) == 0 :
                dtype = "u"+dtype
            dtype = numpy.dtype(dtype)
            
            self._pixel_array = numpy.fromstring(self.pixel_data, dtype).reshape(shape)
            
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
            
        if tag.private :
            pass
#            # See PS 3.5-2008 section 7.8.1 (p. 44) for how blocks are reserved
#            logging.debug("Setting private tag %r" % tag)
#            private_block = tag.elem >> 8
#            private_creator_tag = Tag(tag.group, private_block)
#            if private_creator_tag in self and tag != private_creator_tag:
#                data_element.private_creator = self[private_creator_tag]
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
    
    ###################
    # Misc. functions #
    ###################
    
    def __unicode__(self):
        
        result = []
        for tag in sorted(self.keys()):
            value = self[tag]
            
            if tag.private :
                if tag.element == 0x0010 :
                    name = "Private Creator"
                    vr = None
                else :
                    private_creator = self.get(Tag(tag.group, 0x0010), None)
                    if private_creator not in private_dictionaries :
                        private_creator = None
                    
                    private_tag = "{0:04x}xx{1:02x}".format(tag.group, 
                                                            tag.element & 0x00ff)
                    if (private_creator is not None and 
                        private_creator in private_dictionaries and
                        private_tag in private_dictionaries[private_creator]) :
                        
                        name = private_dictionaries[private_creator][private_tag][2]
                        vr = private_dictionaries[private_creator][private_tag][0]
                    else :
                        name = tag
                        vr = None
            else :
                if tag.group/0x100 in [0x50, 0x60] :
                    # Repeating group element, cf. PS 3.5-2011, 7.6
                    tag_in_dictionary = "{0:02x}xx{1:04x}".format(
                        tag.group/0x100, tag.element)
                else :
                    tag_in_dictionary = tag
                name = data_dictionary.setdefault(
                    tag_in_dictionary, ("UN", "1", unicode(tag_in_dictionary)))[2]
                vr = data_dictionary[tag_in_dictionary][0]
            
            if vr == "SQ" :
                elements = [unicode(item) for item in value]
                value = []
                for index, element in enumerate(elements) :
                    value.append(4*" "+"Item %i"%index)
                    value.extend([8*" "+unicode(x) for x in element.splitlines()])
                value = "\n"+"\n".join(value)
            elif vr in ["OB", "OW", "OB/OW", "OF", "UN"] :
                value = "<array of %i bytes>"%(len(value),)
            elif vr == "UI" and value in dictionary.uid_dictionary :
                value = "{0} ({1})".format(dictionary.uid_dictionary[value][0],
                                           value)
            else :
                value = self[tag]

            if tag.private :
                try :
                    value = unicode(value)
                except UnicodeDecodeError :
                    value = "<array of %i bytes>"%(len(value),)
            
            result.append("%s %s: %s"%(name, tag, value))
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
