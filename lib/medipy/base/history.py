##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

from observable import Observable

class History(Observable) :
    """ History of Commands with multiple undo/redo levels.
    """
    
    def __init__(self, maximum_steps=None) :
        Observable.__init__(self, ["cursor"])
        self.maximum_steps = maximum_steps
        self._steps = []
        self._cursor = None
    
    def add(self, command) :
        """ Execute the command, and add it to the history.
        """
        
        previous = self._cursor
        
        command.execute()
        if self._cursor :
            del self._steps[0:self._cursor]
        self._cursor = 0
        if self.maximum_steps and len(self._steps) == self.maximum_steps :
            del self._steps[-1]
        self._steps.insert(0, command)
        self.notify_observers("cursor", previous=previous)
    
    def undo(self, count=1) :
        """ Undo the current command.
        """
        
        for _ in range(count) :
            self.cursor+=1
    
    def redo(self, count=1) :
        """ Redo the current command.
        """
        
        for _ in range(count) :
            self.cursor-=1

    def get_command(self, index):
        """ Return the command at the given index. The command MUST NOT be 
            executed or undone.
        """
        
        return self._steps[index]

    def _get_cursor(self):
        """ Current position in the history. Latest commands have lower indices
            than earliest commands : the earliest command has index steps_count,
            the latest has index 0.
        """
        
        return self._cursor
    
    def _set_cursor(self, cursor):
        previous = self._cursor
        
        offset = cursor - self._cursor
        if offset <= 0 :
            function = "execute"
        else :
            function = "undo"
        i = self._cursor
        while i != cursor :
            if function == "execute" :
                i -= 1
            step = self._steps[i]
            getattr(step, function)()
            if function == "undo" :
                i += 1
        self._cursor = cursor
        self.notify_observers("cursor", previous=previous)
    
    def _get_steps_count(self) :
        """ Number of steps in the history.
        """
        
        return len(self._steps)
    
    def _get_empty(self) :
        """ Test if the history is empty.
        """
        
        return not self._steps
    
    def _get_can_undo(self) :
        """ Test if an undo operation can be performed.
        """
        
        return (not self.empty) and (self._cursor < len(self._steps))
    
    def _get_can_redo(self) :
        """ Test if an redo operation can be performed.
        """
        
        return (not self.empty) and (self._cursor>0)

    def _get_labels(self) :
        """ Return a list of commands labels, ordered from latest to earliest.
        """
        
        return [command.label for command in self._steps]

    cursor = property(_get_cursor, _set_cursor)
    steps_count = property(_get_steps_count)
    empty = property(_get_empty)
    can_undo = property(_get_can_undo)
    can_redo = property(_get_can_redo)
    labels = property(_get_labels)
