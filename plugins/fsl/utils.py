##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import os
import re
import subprocess

import medipy.base

def environment(fsl_sh = None) :
    """ Return a dictionary of the environment needed by FSL programs. fsl_sh 
        is the path to the fsl.sh script; if absent default locations are 
        searched.
    """
    
    fsl_sh = fsl_sh or "/etc/fsl/fsl.sh"
    
    if not os.path.isfile(fsl_sh) :
        raise medipy.base.Exception("No such file: {0}".format(fsl_sh))
    
    # Use sh commands and a string instead of a list since we're using shell=True
    # Pass empty environment to get only the FSL variables
    command = ". {0} ; /usr/bin/printenv".format(fsl_sh)
    process = subprocess.Popen(command, shell=True, env={}, 
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode != 0 :
        raise medipy.base.Exception("Could not parse fsl.sh: {0}".format(stderr))
    
    # Parse the output : each line should be of the form "VARIABLE_NAME=value"
    fsl_environment = {}
    for line in stdout.split(os.linesep) :
        match = re.match(r"^(\w+)=(\S*)$", line)
        if match :
            name, value = match.groups()
            if name != "PWD" :
                fsl_environment[name] = value
    
    return fsl_environment

def find(fsl_sh = None) :
    """ Return the root of FSL. fsl_sh, if supplied, is the path to the fsl.sh
        script.
        
        THIS FUNCTION IS DEPRECATED.
    """

    fsldir = None
    
    # Use environment if possible
    fsldir = os.environ.get("FSLDIR", None)
    
    # Try fsl.sh
    fsl_sh = fsl_sh or "/etc/fsl/fsl.sh"
    if fsldir is None :
        if os.path.isfile(fsl_sh) :
            lines = open(fsl_sh).readlines()
            for line in lines :
                match = re.match(r"^FSLDIR=(.*)", line)
                if match :
                    fsldir = match.group(1)
    
    return fsldir
