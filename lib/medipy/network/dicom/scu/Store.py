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

# Configuration file template for DCMTK Store SCU.
# Parameters: transfer_syntax and sop_class_uid
config_template = r"""
[[TransferSyntaxes]]

[MediPyTS]
TransferSyntax1 = {transfer_syntax}

[[PresentationContexts]]

[MediPyContext]
PresentationContext1 = {sop_class_uid}\MediPyTS

[[Profiles]]

[MediPyProfile]

PresentationContexts = MediPyContext
"""

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
        
        f, dataset_filename = tempfile.mkstemp()
        os.close(f)
        medipy.io.dicom.write(self.dataset, dataset_filename)
        
        # Generate a custom configuration file to make sure that the presentation
        # context is correct, even for private SOP classes. The TS associated
        # with the presentation context defaults to Explicit VR Little Endian.
        transfer_syntax = getattr(self.dataset, "header", {}).get(
            "transfer_syntax_uid", medipy.io.dicom.UI("1.2.840.10008.1.2.1")).value
        f, config_filename = tempfile.mkstemp()
        config = config_template.format(
            transfer_syntax=transfer_syntax,
            sop_class_uid=self.dataset.sop_class_uid.value)
        os.write(f, config)
        os.close(f)
        
        command = ["storescu"]
        
        command.extend(["--config-file", config_filename, "MediPyProfile"])
        
        command.extend(["--aetitle", self.connection.calling_ae_title])
        command.extend(["--call", self.connection.called_ae_title])
        
        command.append(self.connection.host)
        command.append(str(self.connection.port))
        
        command.append(dataset_filename)
        
        process = subprocess.Popen(command, 
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        
        # Remove the temporary files
        os.unlink(dataset_filename)
        os.unlink(config_filename)
        
        if process.returncode != 0 :
            raise medipy.base.Exception(stderr)
