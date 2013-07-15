##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

class ConnectionBase(object):
    """ Abstract Connection Class
        Use derived classes instead : Connection and SSHTunnelConnection
    """

    def __init__(self, host, port, calling_ae_title, called_ae_title) :
        self.host = host
        self.port = port
        self.calling_ae_title = calling_ae_title
        self.called_ae_title = called_ae_title
        # add transfer syntaxes, timeout, ... ?
        # cf. PS 3.7-2011, 7.1 
    
    def connect(self) :
        raise NotImplementedError()
    
    def disconnect(self) :
        raise NotImplementedError()
