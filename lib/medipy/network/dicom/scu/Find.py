##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import glob
import os
import shutil
import subprocess
import tempfile

import medipy.base

from SCU import SCU

class Find(SCU) :
    """ Find SCU.
    """
    
    def __init__(self, connection, root="study", query_level="patient", query_parameters=None) :
        SCU.__init__(self, connection)
        
        # DIMSE
        self._root=None # patient, study. Patient-study not included, as it is retired
        
        # IOD
        self._query_level=None
        self._query_parameters=None
        # Patient ID (0010,0020) must be known for STUDY level
        # Study Instance UID (0020,000D) for SERIES level
        # Series Instance UID (0020,000E) for IMAGE level
        # TODO : define a list a allowed attributes for each level
        
        # Initialization
        self.root = root
        self.query_level = query_level
        self.query_parameters = query_parameters or medipy.io.dicom.DataSet()
    
    def __call__(self) :
        
        # Create a temporary directory for query files and results
        temporary_directory = tempfile.mkdtemp()
        print "Temporary directory: {0}".format(temporary_directory)
        
        # Build the command
        command = ["findscu"]
        
        command.extend(["--aetitle", self.connection.calling_ae_title])
        command.extend(["--call", self.connection.called_ae_title])
        
        root_option = { "patient" : "-P", "study" : "-S", "psonly" : "-O" }
        command.append(root_option[self.root])
        
        query_level_option = { "patient" : "PATIENT", "study" : "STUDY", 
                               "serie" : "SERIES", "image" : "IMAGE" }
        command.extend([
            "-k", 
            "(0008,0052)={0}".format(query_level_option[self.query_level])
        ])
        
        keys = [(tag, value) for tag, value in self.query_parameters.items()
                if value is not None]
        query = [tag for tag, value in self.query_parameters.items()
                 if value is None]

        for tag, value in keys :
            command.extend(["-k", "{0}={1}".format(str(tag), str(value))])
        
        command.append("-X")
        
        command.append(self.connection.host)
        command.append(str(self.connection.port))
        
        query_file = self._create_query_file(query, temporary_directory)
        command.append(query_file)
        
        #print " ".join(command), temporary_directory
        #return
        
        process = subprocess.Popen(command, 
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=temporary_directory)
        stdout, stderr = process.communicate()
        
        if process.returncode != 0 :
            raise medipy.base.Exception(stderr)
        
        # Process results
        datasets = []
        for filename in glob.glob(os.path.join(temporary_directory, "rsp*.dcm")) :
            dataset = medipy.io.dicom.parse(filename)
            datasets.append(dataset)
        
        # Remove temporary directory
        shutil.rmtree(temporary_directory)
        
        return datasets
    
    ##############
    # Properties #
    ##############
    
    def _get_root(self) :
        return self._root
    
    def _set_root(self, root) :
        if root not in ["patient", "study"] :
            raise medipy.base.Exception("Unknown Find root: {0}".format(root))
        self._root = root
    
    def _get_query_level(self) :
        return self._query_level
    
    def _set_query_level(self, query_level) :
        if query_level not in ["patient", "study", "serie", "image"] :
            raise medipy.base.Exception("Unknown query level: {0}".format(query_level))
        if self.query_parameters is not None :
            # Check mandatory attributes
            pass
        self._query_level = query_level
    
    def _get_query_parameters(self) :
        return self._query_parameters
    
    def _set_query_parameters(self, query_parameters) :
        if self.query_level is not None :
            # Check mandatory attributes
            pass
        # Check allowed attributes
        pass
        
        self._query_parameters = query_parameters
    
    root = property(_get_root, _set_root)
    query_level = property(_get_query_level, _set_query_level)
    query_parameters = property(_get_query_parameters, _set_query_parameters)
    
    #####################
    # Private interface #
    #####################
    
    def _create_query_file(self, query, temporary_directory) :
        f, query_file_txt = tempfile.mkstemp(dir=temporary_directory)
        for tag in query :
            vr = "LO"
            os.write(f, "{0} {1} []\n".format(tag, vr))
        os.close(f)
        
        f, query_file_dcm = tempfile.mkstemp(dir=temporary_directory)
        os.close(f)
        
        command = ["dump2dcm", query_file_txt, query_file_dcm]
        subprocess.call(command)
        # TODO : Popen
        
        return query_file_dcm
