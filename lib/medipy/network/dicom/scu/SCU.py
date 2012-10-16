##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

class SCU(object) :
    """ Base class for all DICOM Service Class Users. 
    
        All derived classes must implement the __call__ function, which performs
        the actual operation.
    """
    
    def __init__(self, connection) :
        self.connection = connection
    
    def __call__(self, *args, **kwargs):
        raise NotImplementedError()
