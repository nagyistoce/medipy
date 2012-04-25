##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

class Tag(int):
    
    def __new__(cls, arg, arg2=None):
        if arg2 is None and isinstance(arg, (int, long)) :
            value = arg
        elif arg2 is None and isinstance(arg, (list, tuple)) :
            value = (arg[0]<<16) + arg[1]
        else :
            value = (arg<<16) + arg2
        
        return int.__new__(cls, value)
    
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
        return self >> 16
    
    def _get_element(self):
        return self & 0xffff
    
    def _get_private(self):
        return (self.group%2 == 1)
    
    group = property(_get_group)
    element = property(_get_element)
    private = property(_get_private)
