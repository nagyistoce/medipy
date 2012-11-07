##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

from observable import Observable

class ObservableList(list, Observable) :
    """List object with observer support.
    
    This class mimics the behavior of the Python builtin list, and adds support
    for observers when the list is modified.
    
    The events are fired after the list is notified, and have the following
    arguments. The arguments are sufficient to have all information about the
    modification.
        
        * set_item : self, index, old_value
        * delete_item : self, index, old_value
        * append : self
        * pop : self, index, popped_value
        * extend : self, added_value
        * insert : self, index, value
        * remove : self, index_of_value, value
        * reverse : self
        * sort : self, cmpfunc, key, reverse, old_list
    """ 
    
    def __init__(self, sequence=[]) :
        list.__init__(self, sequence)
        Observable.__init__(self, ["set_item", "delete_item", "append", "pop",
            "extend", "insert", "remove", "reverse", "sort"])
    
    def __setitem__(self, index, value) :
        old_value = self.__getitem__(index)
        list.__setitem__(self, index, value)
        self.notify_observers("set_item", index=index, old_value=old_value)
    
    def __delitem__(self, index) :
        old_value = list.__getitem__(self, index)
        list.__delitem__(self, index)
        self.notify_observers("delete_item", index=index, old_value=old_value)
    
    def __setslice__(self, i, j, value) :
        self.__setitem__(slice(i,j), value)
    
    def __delslice__(self, i, j) :
        self.__delitem__(slice(i,j))

    def append(self, value) :
        list.append(self, value)
        self.notify_observers("append")
    
    def pop(self, *args) :
        popped_value = list.pop(self, *args)
        if len(args) > 0 :
            index = args[0]
        else :
            # Do not subtract one, as we just popped
            index = len(self)
        self.notify_observers("pop", index=index, popped_value=popped_value)
        return popped_value
    
    def extend(self, value) :
        list.extend(self, value)
        self.notify_observers("extend", value=value)
    
    def insert(self, index, value) :
        list.insert(self, index, value)
        self.notify_observers("insert", index=index, value=value)
            
    def remove(self, value) :
        index = list.index(self, value)
        list.remove(self, value)
        self.notify_observers("remove", index=index, value=value)
    
    def reverse(self) :
        list.reverse(self)
        self.notify_observers("reverse")
    
    def sort(self, cmpfunc=None, key=None, reverse=False) :
        old_list = self[:]
        list.sort(self,cmpfunc, key, reverse)
        self.notify_observers("sort", cmpfunc=cmpfunc, key=key, reverse=reverse, old_list=old_list)
