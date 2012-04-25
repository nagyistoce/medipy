##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import os
import re

def find(fsl_sh = None) :
    """ Return the root of FSL. fsl_sh, if supplied, is the path to the fsl.sh
        script.
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