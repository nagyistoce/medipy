##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

""" Command classes following the Command design pattern.
"""

class Command(object) :
    """ Abstract base class for all commands.
    """
    
    def __init__(self, label="") : 
        self.label = label
    
    def execute(self) :
        """ Execute the command. This function *must* be reimplemented in 
            derived classes.
        """
        
        raise NotImplementedError()

class UndoableCommand(Command) :
    """ Abstract base class for all undo-able commands.
    """
    
    def __init__(self, label="") :
        Command.__init__(self, label)
        
    def undo(self) :
        """ Execute the command. This function *must* be reimplemented in 
            derived classes.
        """
        
        raise NotImplementedError()
