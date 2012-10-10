##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

_Exception = Exception

class Exception(_Exception):
    def __init__(self, *args, **kwargs):
        _Exception.__init__(self, *args, **kwargs)
