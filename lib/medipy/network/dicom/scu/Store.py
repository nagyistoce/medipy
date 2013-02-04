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
    """ The Store SCU stores DataSets on a DICOM node. This SCU has only one
        parameter: :attr:`dataset`, the DataSet to be stored. The Store SCU does
        not return any value, but raises an exception if the data set could not
        be stored.::
        
                connection = medipy.io.dicom.Connection("pacs.example.com", 104, 
                    "MY_MACHINE", "REMOTE_PACS")
        
                dataset = medipy.io.dicom.DataSet()
                # Fill the dataset ...
                
                store = medipy.network.dicom.scu.Store(connection, dataset)
                store()
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
