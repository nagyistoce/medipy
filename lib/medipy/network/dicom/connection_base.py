##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

from medipy.base import LateBindingProperty

class ConnectionBase(object):
    """ Base class for DICOM connections from SCUs to SCPs.
    """

    def __init__(self, host, port, calling_ae_title, called_ae_title, connect=False) :
        
        #Properties private attributes
        self._host=None
        self._port=None
        self._calling_ae_title=None
        self._called_ae_title=None
        
        self._set_host(host)
        self._set_port(port)
        self._set_calling_ae_title(calling_ae_title)
        self._set_called_ae_title(called_ae_title)
        
        # add transfer syntaxes, timeout, ... ?
        # cf. PS 3.7-2011, 7.1 
        
        if connect :
            self.connect()
    
    def connect(self) :
        """ Open the connection.
        """
        
        raise NotImplementedError()
    
    def disconnect(self) :
        """ Close the connection.
        """
        
        raise NotImplementedError()
        
    ##############
    # Properties #
    ##############
    def _get_host(self):
        """ Host to connect to.
        """
        return self._host
        
    def _set_host(self,host):
        self._host = host

    def _get_port(self):
        """ TCP port to connect to.
        """
        return self._port
        
    def _set_port(self,port):
        self._port = port

    def _get_calling_ae_title(self):
        """ Calling (i.e. local) AE title.
        """
        return self._calling_ae_title
        
    def _set_calling_ae_title(self,calling_ae_title):
        self._calling_ae_title = calling_ae_title
        
    def _get_called_ae_title(self):
        """ Called (i.e. remote) AE title.
        """
        return self._called_ae_title
        
    def _set_called_ae_title(self,called_ae_title):
        self._called_ae_title = called_ae_title
        
    host = LateBindingProperty(_get_host, _set_host)
    port = LateBindingProperty(_get_port, _set_port)
    calling_ae_title = LateBindingProperty(_get_calling_ae_title, _set_calling_ae_title)
    called_ae_title = LateBindingProperty(_get_called_ae_title, _set_called_ae_title)
