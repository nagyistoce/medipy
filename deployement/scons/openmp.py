##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import sys

def configuration_variables() :
    if sys.platform.startswith("win") :
        result = {
            "CXXFLAGS" : ["/openmp"],
			"CFLAGS" : ["/openmp"],			
        }
    else :
        result = {
            "CXXFLAGS" : ["-fopenmp"],
			"CFLAGS" : ["-fopenmp"],
            "LIBS" : ["gomp"]
        }
    
    return result
