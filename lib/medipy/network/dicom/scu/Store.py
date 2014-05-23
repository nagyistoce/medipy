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
import shutil

import medipy.io.dicom

from scu import SCU

# Configuration file template for DCMTK Store SCU.
# Parameters: transfer_syntax and sop_class_uid
config_template = r"""
[[TransferSyntaxes]]

[MediPyTS]
{syntaxes}

[[PresentationContexts]]

[MediPyContext]
{classes}

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
    
    def __init__(self, connection, datasets) :
        SCU.__init__(self, connection)
        
        if isinstance(datasets,medipy.io.dicom.DataSet):
            datasets = [datasets]
        
        self.datasets = datasets
    
    def __call__(self) :
                
        # Generate a custom configuration file to make sure that the presentation
        # context is correct, even for private SOP classes. The TS associated
        # with the presentation context defaults to Explicit VR Little Endian.
        
        main_tmpdir = tempfile.mkdtemp()
        
        contexts={}
        for dataset in self.datasets:
            transfer_syntax = getattr(dataset, "header", {}).get(
            "transfer_syntax_uid", medipy.io.dicom.UI("1.2.840.10008.1.2.1")).value
            sop_class_uid = dataset.sop_class_uid.value
            contexts.setdefault((transfer_syntax,sop_class_uid),[]).append(dataset)
        
        keys = contexts.keys()
        packages = [keys[i:i+128] for i in range(0, len(keys), 128)]
        datasets = []
        for package in packages :
            datasets_tmpdir = tempfile.mkdtemp(dir=main_tmpdir)
            syntaxes = set()
            classes = set()
            for ts, class_ in package :
                syntaxes.add(ts)
                classes.add(class_)
            syntaxes = "\n".join(
                "TransferSyntax{0} = {1}".format(1+index, ts) for index, ts in enumerate(syntaxes))
            classes = "\n".join(
                "PresentationContext{0} = {1}\\MediPyTS".format(1+index, class_) for index, class_ in enumerate(classes))
            for item in package:
                datasets.extend(contexts[item])
        
            config = config_template.format(syntaxes=syntaxes, classes=classes)
            
            # Save config
            f, config_filename = tempfile.mkstemp(dir=main_tmpdir)
            os.write(f,config)
            os.close(f)
            
            # Save datasets
            for dataset in datasets:
                f, dataset_filename = tempfile.mkstemp(dir=datasets_tmpdir)
                os.close(f)
                medipy.io.dicom.write(dataset, dataset_filename)
            
            # Call storescu
            self._build_command(config_filename, datasets_tmpdir)

            # Remove the temporary files
            shutil.rmtree(datasets_tmpdir)
            os.unlink(config_filename)
            
        shutil.rmtree(main_tmpdir)
        
    def _build_command(self,config_filename,dataset_tmpdir):
        
        command = ["storescu"]
        
        command.extend(["--config-file", config_filename, "MediPyProfile"])
        command.append("--scan-directories")
        command.extend(["--aetitle", self.connection.calling_ae_title])
        command.extend(["--call", self.connection.called_ae_title])
        
        command.append(self.connection.host)
        command.append(str(self.connection.port))
        command.append(dataset_tmpdir)
        
        process = subprocess.Popen(command, 
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        
        if process.returncode != 0 :
            raise medipy.base.Exception(stderr)
