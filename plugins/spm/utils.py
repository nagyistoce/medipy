##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import subprocess

import medipy.base

def find(matlab=None, matlab_path=None):
    """ Return the root directory of SPM. matlab, if given, is the path to the
        MATLAB executable. matlab_path, if given, is a MATLAB expression fed to
        addpath.
    """
    
    script = ("if isempty(which('spm'));"
              "fprintf(1, '');"
              "exit();"
              "end;"
              "fprintf(1, '%s', spm('dir'));"
              "exit();")

    if matlab_path :
        script = "addpath({0});".format(matlab_path)+script

    command = [matlab or "matlab", 
               "-nodisplay", "-nosplash", "-nojvm", 
               "-r", script]
    
    try :
        process = subprocess.Popen(command, stdout=subprocess.PIPE)
    except OSError :
        # Matlab is not present
        raise medipy.base.Exception("Could not find SPM")
    stdout = process.communicate()[0]
    last_line = stdout.split("\n")[-1]
    # Weird data at the end of the line
    if '\x1b' in last_line :
        last_line = last_line[:last_line.index('\x1b')]
    if last_line == "" :
        raise medipy.base.Exception("Could not find SPM")
    
    return last_line
