##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import operator

class VR(object):
    """ Base class for all VR classes.
    """
    
    def __init__(self, value):
        self.value = value

# Dynamically define all the VR classes
vrs = ["AE", "AS", "AT", "CS", "DA", "DS", "DT", "FL", "FD", "IS", 
       "LO", "LT", "OB", "OF", "OW", "PN", "SH", "SL", "SQ", "SS", 
       "ST", "TM", "UI", "UL" "UN", "US", "UT"]
for vr in vrs :
    # globals() is the dictionary of the current module.
    globals()[vr] = type(vr, (VR,), dict())
