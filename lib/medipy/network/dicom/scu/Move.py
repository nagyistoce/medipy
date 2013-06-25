##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
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
import medipy.io.dicom

from SCU import SCU

class Move(SCU):

    def __init__(self, connection, root="patient", query_level="patient",
        move_destination=None, query_parameters=None) :

        SCU.__init__(self, connection)
        
        self.root = root
        self.query_level = query_level
        self.move_destination = move_destination
        self.connection = connection
        self.query_parameters = query_parameters or medipy.io.dicom.DataSet()
        
    def __call__(self):
        
        # Create a temporary directory for query files and results
        temporary_directory = tempfile.mkdtemp()
        
        #Build the command    
        command = ["movescu"]
        command.append("-v")
        
        root_option = { "patient":"-P", "study":"-S" }
        command.append(root_option[self.root])
        
        command.extend(["--aetitle", self.connection.calling_ae_title])
        command.extend(["--port", str(self.connection.port)])
        command.extend(["--call", self.connection.called_ae_title])
        command.extend(["--move", self.move_destination])
        
        for key in self.query_parameters.keys():
            command.extend([
                "-k",
                str(key)[1:-1]+"={0}".format(self.query_parameters[key].value)
            ])
        
        query_level_option = { "patient" : "PATIENT", "study" : "STUDY", 
                               "series" : "SERIES", "image" : "IMAGE" }
        command.extend([
            "-k","0008,0052={0}".format(query_level_option[self.query_level])
        ])
  
        command.append(self.connection.host)
        command.append(str(self.connection.port))
        
        query_file = self._create_query_file(self.query_parameters, temporary_directory)
        command.append(str(query_file)) 
        
        #Subshell : exec cmd
        process = subprocess.Popen(command, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, cwd=temporary_directory)
        stdout, stderr = process.communicate()
        
        if process.returncode != 0:
            raise medipy.base.Exception(stderr)
        
        # Process results
        datasets = []
        for filename in glob.glob(os.path.join(temporary_directory, "*")):
            if not filename in glob.glob(os.path.join(temporary_directory, "tmp*")):
                dataset = medipy.io.dicom.read(filename)
                datasets.append(dataset)
        
        # Remove temporary directory
        shutil.rmtree(temporary_directory)
            
        return datasets


    def _create_query_file(self, query, temporary_directory) :
        f, query_file = tempfile.mkstemp(dir=temporary_directory)
        os.close(f)
        medipy.io.dicom.write(query, query_file)
        
        return query_file

if __name__ == "__main__":

    import medipy.io.dicom
    import medipy.network.dicom
    
    connection = medipy.network.dicom.Connection("aude.u-strasbg.fr", 11112,
        "forez", "PIIV-RECHERCHE")
        
    query = medipy.io.dicom.DataSet(patient_id = "010027LA",series_description="*3d")
    for key in ["patient_id", "study_instance_uid", "series_instance_uid",
            "sop_instance_uid"] :
        query.setdefault(key, None)
                    
    datasets =  medipy.network.dicom.query.relational(connection,"patient",
                "patient", query)
    
    print datasets.patient_id.value
                
    retrieve_query=[]
    for dataset in datasets:
        retrieve_query.append(medipy.io.dicom.DataSet(
                patient_id = dataset.patient_id.value,
                study_instance_uid = dataset.study_instance_uid.value,
                series_instance_uid = dataset.series_instance_uid.value,
                sop_instance_uid = dataset.sop_instance_uid.value))

    for query in retrieve_query:
        move = medipy.network.dicom.scu.Move(connection, "patient", "patient",
            "forez", query)
        results = move()

#    print results

