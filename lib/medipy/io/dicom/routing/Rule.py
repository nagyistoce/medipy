##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

class Rule(object):
    """ Apply actions on a dataset if it matches a condition
    """
    
    def __init__(self, condition, actions, enabled=True):
        self.condition = condition
        self.actions = actions
        self.enabled = enabled
    
    def __call__(self, dataset):
        if not self.enabled or not self.condition(dataset) :
            return
        
        for action in self.actions :
            action(dataset)
