##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import numpy

from base import Tool

class ModelEstimation(Tool):
    """ Model estimation
    """
    
    name = "stats.fmri_est"
    
    class Method(Tool.Config):
        """ Estimation method, the method can be
              * classical
              * bayesian_1
              * bayesian_2
        """
        
        def __init__(self, method="classical"):
            self.method = method
        
        def _get_parameters(self):
            result = {}
            
            if self.method == "classical" :
                result["Classical"] = 1
            elif self.method == "bayesian_2" :
                result["Bayesian2"] = 1
            
            return result
    
    def __init__(self, root, design_matrix, method=None) :
        Tool.__init__(self, root)
        self.design_matrix = design_matrix
        self.method = method or ModelEstimation.Method()
    
    def _get_script(self):
        script = []
        script.extend(Tool._generate_script(
            self.name, {"spmmat" : [self.design_matrix]}))
        script.extend(Tool._generate_script(
            "{0}.method".format(self.name), self.method.parameters))
        
        return script