##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import sys

# Determine int class to use : DICOM tags higher than (7fff,ffff) cannot be
# represented by a 32-bits int
IntergerClass = int if sys.maxint >= 0xffffffffL else long 

class Tag(IntergerClass):
    """ Tag of a DICOM element. The tag can be either created by two 16-bits
        integers or by one 32-bits integer : ::
        
            # Both are tags for Patient's Name (0010,0010)
            tag1 = medipy.io.dicom.Tag(0x0010,0x0010)
            tag2 = medipy.io.dicom.Tag(0x00100010)
    """
    
    def __new__(cls, arg, arg2=None):
        if arg2 is None and isinstance(arg, (int, long)) :
            value = arg
        elif arg2 is None and isinstance(arg, (list, tuple)) :
            value = (arg[0]<<16) + arg[1]
        else :
            value = (arg<<16) + arg2
        
        return IntergerClass.__new__(cls, value)
    
    def __str__(self):
        return "(%04x,%04x)"%(self.group, self.element)
    
    def __repr__(self):
        return "(0x%04x,0x%04x)"%(self.group, self.element)
    
    def __eq__(self, other):
        if not isinstance(other, Tag) :
            other = Tag(other)
        return (int(self) == int(other))
    
    def __ne__(self, other):
        return not self.__eq__(other)
    
    def _get_group(self):
        """ Group of the tag, i.e. the first two bytes.
        
            >>> tag = medipy.io.dicom.Tag(0x0010, 0x0020)
            >>> tag.group
            0x0010
        """
        
        return self >> 16
    
    def _get_element(self):
        """ Element of the tag, i.e. the last two bytes.
        
            >>> tag = medipy.io.dicom.Tag(0x0010, 0x0020)
            >>> tag.element
            0x0020
        """
        return self & 0xffff
    
    def _get_private(self):
        """ Whether the tag is private (i.e. odd group).
        
            >>> tag1 = medipy.io.dicom.Tag(0x0010, 0x0010)
            >>> tag1.private
            False
            >>> tag2 = medipy.io.dicom.Tag(0x0029, 0x1047)
            >>> tag2.private
            True
        """
        return (self.group%2 == 1)
    
    group = property(_get_group)
    element = property(_get_element)
    private = property(_get_private)
