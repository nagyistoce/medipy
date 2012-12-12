##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import os
import subprocess
import tempfile

import medipy.io.dicom

from SCU import SCU

class Store(SCU):
    """ Store SCU.
    """
    
    def __init__(self, connection, dataset) :
        SCU.__init__(self, connection)
        self.dataset = dataset
    
    def __call__(self) :
        
        f, filename = tempfile.mkstemp()
        os.close(f)
        medipy.io.dicom.write(self.dataset, filename)
        
        command = ["storescu"]
        
        command.extend(["--aetitle", self.connection.calling_ae_title])
        command.extend(["--call", self.connection.called_ae_title])
        
        command.append(self.connection.host)
        command.append(str(self.connection.port))
        
        command.append(filename)
        
        process = subprocess.Popen(command, 
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        
        os.unlink(filename)
        
        if process.returncode != 0 :
            raise medipy.base.Exception(stderr)
