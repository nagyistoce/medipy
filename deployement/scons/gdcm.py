##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

""" Configuration for GDCM.
    
    The "GDCM" dictionary will be added to the environment, with the following
    keys :
    * "CPPPATH"
    * "LIBPATH"
    * "LIBS"
"""

import os.path
import re
import subprocess
import sys

from utils import split_environment_variable

def configuration_variables() :
    command = ["cmake", "-P", os.path.join(os.path.dirname(__file__), "gdcm.cmake")]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode != 0 :
        print "Could not find VTK, some modules will not be compiled."
        print stderr
        return None
    else :
        result = {}
        # Result of cmake -P is on stderr
        for line in stderr.split("\n") :
            if not line :
                continue
            key, value = line.split("=", 1)
            if key in ["CPPPATH", "LIBPATH", "LIBS" ] :
                value = list(set(value.split(";")))
            if key == "LIBS" :
                value = [x for x in value if x != "general"]
            result[key] = value
    # Remove absolute file paths from LIBS
    result["LIBS"] = [x for x in result["LIBS"] if os.path.sep not in x]
    return result

def exists(env):
    return (configuration_variables() is not None)

def generate(env):
    env["GDCM"] = configuration_variables()
    
