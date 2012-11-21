##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import medipy.base

class IOBase(object):
    """ Base class for all IO classes.
    """
    
    def __init__(self, filename=None, report_progress=None):
        self._filename = filename
        self._report_progress = report_progress
    
    def _get_filename(self):
        return self._filename
    
    def _set_filename(self, filename):
        self._filename = filename
    
    filename = medipy.base.LateBindingProperty(_get_filename, _set_filename)