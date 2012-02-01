##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import sys
import wx
from medipy.io.dicom import dictionary

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
            self._add_dataset(dataset.perframe_functional_groups_sequence[frame_number])
        
        for i in range(self.GetColumnCount()):
            self.SetColumnWidth(i, wx.LIST_AUTOSIZE)
            width = self.GetColumnWidth(i)
            self.SetColumnWidth(i, wx.LIST_AUTOSIZE_USEHEADER)
            width = max(width, self.GetColumnWidth(i))
            self.SetColumnWidth(i, width)
        
    def _add_dataset(self, dataset, indent=""):
        for tag in sorted(dataset) :
            if tag.private :
                continue
            if tag == 0x7fe00010 : # Pixel Data
                continue
            if tag == 0x52009230 : # Per-frame Functional Groups Sequence
                continue
            
            entry = dictionary.data_dictionary.get(tag, ("UN", "1", str(tag), False, str(tag)))
            
            description = indent+entry[2]
            item = self.InsertStringItem(sys.maxint, description)
            self.SetStringItem(item, 1, "({0:04x},{1:04x})".format(tag.group, tag.element))
            
            value = dataset[tag]
            vr = indent+entry[0]
            
            if vr != "SQ" :
                value = self._string_representation(value, vr)
                self.SetStringItem(item, 2, value)
            else :
                self.SetStringItem(item, 2, "{0} items".format(len(value)))
                self._add_sequence(value, "    "+indent)
        
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
