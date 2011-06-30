##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import copy
import os
import sys

def split_environment_variable(variable):
    # List separator in environment variables
    if sys.platform == "win32" :
        separator = ";"
    else : 
        separator = ":"
    
    elements = os.environ[variable].split(separator) if variable in os.environ else []
    
    return elements

def merge_construction_variables(*args):
    result = copy.copy(args[0])
    
    for env in args[1:] :
        for key in env : 
            if key in result : 
                result[key].extend(env[key])
            else : 
                result[key] = copy.copy(env[key])
    
    return result