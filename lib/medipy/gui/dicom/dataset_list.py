##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import sys
import wx
from medipy.io.dicom import Tag, dictionary
from medipy.io.dicom.private_dictionaries import private_dictionaries

class DataSetList(wx.ListCtrl) :
    """ List of DICOM attributes from a Data Set.
    """
    
    def __init__(self, parent, ID=wx.ID_ANY, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=wx.LC_REPORT, 
                 validator=wx.DefaultValidator, name = "datasetlist") :
        
        wx.ListCtrl.__init__(self, parent, ID, pos, size, style, validator, name)
        
        self.InsertColumn(0, "Name")
        self.InsertColumn(1, "Tag")
        self.InsertColumn(2, "Value")
    
    def set_dataset(self, dataset, frame_number=None) :
        """ Set the information entities whose modules will appear in the list.
            If frame_number is not None, then the image is considered to be
            multi-frame, and the Per-frame Functional Groups of the given frame
            number will be displayed in a separate entry. If frame_number is
            None, then all the Per-frame Functional Groups will be displayed in
            the Image entry.
        """
        
        self._update_ui(dataset, frame_number)
    
    #####################
    # Private interface #
    #####################
    
    def _update_ui(self, dataset, frame_number) :
        """ Update the contents of the list.
        """
        
        self.ClearAll()
        
        self.InsertColumn(0, "Name")
        self.InsertColumn(1, "Tag")
        self.InsertColumn(2, "Value")
        
        self._add_dataset(dataset)
        
        if frame_number is not None :
            item = self.InsertStringItem(sys.maxint, "Per-frame elements")
            self.SetItemBackgroundColour(item, (201,237,255))
            if "perframe_functional_groups_sequence" in dataset :
                self._add_dataset(dataset.perframe_functional_groups_sequence[frame_number])
        
        for i in range(self.GetColumnCount()):
            self.SetColumnWidth(i, wx.LIST_AUTOSIZE)
            width = self.GetColumnWidth(i)
            self.SetColumnWidth(i, wx.LIST_AUTOSIZE_USEHEADER)
            width = max(width, self.GetColumnWidth(i))
            self.SetColumnWidth(i, width)
        
    def _add_dataset(self, dataset, indent=""):
        
        for tag in sorted(dataset.keys()):
            value = dataset[tag]
            
            if tag.private :
                if tag.element == 0x0010 :
                    name = "Private Creator"
                    vr = None
                else :
                    private_creator = dataset.get(Tag(tag.group, 0x0010), None)
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
                        name = unicode(tag)
                        vr = None
            else :
                if tag.group/0x100 in [0x50, 0x60] :
                    # Repeating group element, cf. PS 3.5-2011, 7.6
                    tag_in_dictionary = "{0:02x}xx{1:04x}".format(
                        tag.group/0x100, tag.element)
                else :
                    tag_in_dictionary = tag
                name = dictionary.data_dictionary.setdefault(
                    tag_in_dictionary, ("UN", "1", unicode(tag_in_dictionary)))[2]
                vr = dictionary.data_dictionary[tag_in_dictionary][0]
            
            if vr in ["OB", "OW", "OB/OW", "OF", "UN"] :
                value = "<array of %i bytes>"%(len(value),)
            else :
                value = dataset[tag]

            if tag.private :
                try :
                    value = unicode(value)
                except UnicodeDecodeError :
                    value = "<array of %i bytes>"%(len(value),)
            
            description = indent+name
            item = self.InsertStringItem(sys.maxint, description)
            self.SetStringItem(item, 1, "({0:04x},{1:04x})".format(tag.group, tag.element))
            
            if vr != "SQ" :
                self.SetStringItem(item, 2, unicode(value))
            else :
                self.SetStringItem(item, 2, "{0} item{1}".format(
                    len(value), "s" if len(value)>1 else ""))
                self._add_sequence(dataset[tag], "    "+indent)
        
    def _add_sequence(self, sequence, indent):
        for index, item in enumerate(sequence) :
            self._add_dataset(item, indent)
            if index != len(sequence)-1 :
                listctrl_item = self.InsertStringItem(sys.maxint, "")
                self.SetStringItem(listctrl_item, 2, 20*"-")
    
    @staticmethod
    def _string_representation(value, vr):
        """ Return the string representation of DICOM value knowing its VR.
        """
        
        if vr in ["OB", "OF", "OW", "OW/OB", "UN"] : #"US or SS"]
            result = "<binary data, {0} bytes>".format(len(value))
        elif isinstance(value, (int, float, list, tuple)) :
            try :
                result = str(value)
            except ValueError :
                result = "<no string representation available>"
        else :
            # Should be string-like
            try :
                result = unicode(value)
            except UnicodeDecodeError :
                result = "<no string representation available>"
        
        return result
