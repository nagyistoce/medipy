##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

""" Builders for Swig.
    
    The "Swig" dictionary will be added to the environment, with the following
    keys :
    * "SWIGFLAGS"
    * "SHLIBPREFIX"
    * "SHLIBSUFFIX" (conditional)
      
"""

import sys

def configuration_variables() :
    variables = {
        "SWIGFLAGS" : ["-python", "-c++"],
        "SHLIBPREFIX" : "",
    }
    
    if sys.platform == "win32" and sys.version >= '2.5':
        variables["SHLIBSUFFIX"] = ".pyd"
    
    return variables

def exists(env):
    return env.Detect("swig")

def generate(env):
    env["Swig"] = configuration_variables()
    