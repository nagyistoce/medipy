##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import numpy

from base import Tool

class ContrastManager(Tool):
    """ Contrast manager to set up T and F contrasts
    """
    
    name = "stats.con"
    
    class TContrast(Tool.Config):
        name = "tcon"
        def __init__(self, name, vector, replicate="TODO"):
            self.name = name
            self.vector = vector
            self.replicate = replicate
        
        def _get_parameters(self):
            result = {
                "name" : self.name,
                "convec" : numpy.asarray(self.vector),
                "sessrep" : "TODO",
            }
            return result
    
    class FContrast(Tool.Config):
        pass
    
    class ConditionTContrast(Tool.Config):
        pass
    
    def __init__(self, root, design_matrix, contrasts=None, 
                 delete_existing_contrasts=False):
        Tool.__init__(self, root)
        self.design_matrix = design_matrix
        self.contrasts = contrasts or []
        self.delete_existing_contrasts = delete_existing_contrasts
    
    def _get_script(self):
        script = []
        
        script.extend(Tool._generate_script(
            self.name, {"spmmat" : [self.design_matrix]}))
        for index, contrast in enumerate(self.contrasts) :
            script.extend(Tool._generate_script(
                "{0}.{1}({2})".format(self.name, contrast.name, 1+index), 
                contrast.parameters))
        script.extend(Tool._generate_script(
            self.name, {"delete" : int(self.delete_existing_contrasts)}))
        
        return script