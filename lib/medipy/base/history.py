##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

class History(object) :
    """ History of Commands with multiple undo/redo levels
    """
    
    def __init__(self, maximum_steps=None) :
        self.maximum_steps = maximum_steps
        self._steps = []
        self._cursor = None
    
    def add(self, command) :
        """ Execute the command, and add it to the history.
        """
        
        command.execute()
        if self._cursor :
            del self._steps[0:self._cursor]
        self._cursor = 0
        if self.maximum_steps and len(self._steps) == self.maximum_steps :
            del self._steps[-1]
        self._steps.insert(0, command)
    
    def undo(self, count=1) :
        """ Undo the current command.
        """
        
        for _ in range(count) :
            self._steps[self._cursor].undo()
            self._cursor+=1
    
    def redo(self, count=1) :
        """ Redo the current command.
        """
        
        for _ in range(count) :
            self._cursor-=1
            self._steps[self._cursor].execute()

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

    empty = property(_get_empty)
    can_undo = property(_get_can_undo)
    can_redo = property(_get_can_redo)
    labels = property(_get_labels)
