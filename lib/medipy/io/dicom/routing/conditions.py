##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

class Condition(object):
    """ Base class for all routing conditions. All concrete derived classes must
        implement the __call__ method.
    """
    
    def __call__(self, dataset):
        raise NotImplementedError()

class AlwaysTrue(Condition):
    """ Condition which is always True
    """
    
    def __call__(self, dataset):
        return True

class AlwaysFalse(Condition):
    """ Condition which is always False
    """
    
    def __call__(self, dataset):
        return False

class And(Condition):
    """ And-combination of several conditions. 
    """
    
    def __init__(self, *args):
        self._children = args
    
    def __call__(self, dataset):
        return all([child(dataset) for child in self._children])

class Or(Condition):
    """ And-combination of several conditions. 
    """
    
    def __init__(self, *args):
        self._children = args
    
    def __call__(self, dataset):
        return any([child(dataset) for child in self._children])

class Not(Condition):
    """ Inverse of a condition. 
    """
    
    def __init__(self, child):
        self._child = child
    
    def __call__(self, dataset):
        return not self._child(dataset)

class ElementMatch(Condition):
    """ Test if an element in a dataset matches a value. 
    """
    
    def __init__(self, tag, value):
        self.tag = tag
        self.value = value
    
    def __call__(self, dataset):
        # TODO : matching in the DICOM sense (cf. query)
        if self.tag not in dataset :
            return False
        else :
            return dataset[self.tag].value == self.value
