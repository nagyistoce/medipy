##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import subprocess
import medipy.base

from SCU import SCU

class Echo(SCU) :
    """ Echo SCU, i.e. DICOM Ping.
    """
    
    def __init__(self, connection) :
        SCU.__init__(self, connection)
        # Nothing else for the SOP : DIMSE is empty, as is IOD
    
    def __call__(self) :
        command = ["echoscu"]
        
        command.extend(["--aetitle", self.connection.calling_ae_title])
        command.extend(["--call", self.connection.called_ae_title])
        
        command.append(self.connection.host)
        command.append(str(self.connection.port))
        
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        
        if process.returncode != 0 :
            raise medipy.base.Exception(stderr)
