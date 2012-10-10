##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import numpy

from base import Tool

class ResultReport(Tool):
    """ Result report of statistical tests
    """
    
    name = "stats.results"
    
    class Contrast(Tool.Config):
        def __init__(self, name, contrasts, threshold_type="FWE", threshold=0.05,
                     extent=0, masking=None) :
            self.name = name
            self.contrasts = contrasts
            self.threshold_type = threshold_type
            self.threshold = threshold
            self.extent = extent
            self.masking = masking
        
        def _get_parameters(self):
            result = {
                "titlestr" : self.name,
                "contrasts" : numpy.asarray(self.contrasts),
                "threshdesc" : "none" if self.threshold_type is None else self.threshold_type,
                "thresh" : self.threshold,
                "extent" : self.extent,
                #"mask" : None
            }
            
            return result
    
    def __init__(self, root, design_matrix, contrasts, data_type=1, 
                 print_results=True) :
        Tool.__init__(self, root)
        self.design_matrix = design_matrix
        self.contrasts = contrasts
        self.data_type = data_type
        self.print_results = print_results
    
    def _get_script(self):
        script = []
        
        script.extend(Tool._generate_script(
            self.name, {"spmmat" : [self.design_matrix]}))
        for index, contrast in enumerate(self.contrasts) :
            script.extend(Tool._generate_script(
                "{0}.conspec({1})".format(self.name, 1+index), 
                contrast.parameters))
        
        return script