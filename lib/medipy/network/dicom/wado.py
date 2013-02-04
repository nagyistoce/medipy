##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

""" Other than Service Classes that can copy data 
    (:class:`~medipy.network.dicom.scu.Get` and 
    :class:`~medipy.network.dicom.scu.Move`), the DICOM standard contains a 
    method to access DICOM data through HTTP: *WADO*, for *Web Access to DICOM 
    Objects* (PS 3.18-2011). Through this service, the user can access DICOM
    data with and HTTP GET request that specifies the key of a given DICOM object.
"""

import os
import tempfile
import urllib
import urlparse

import medipy.io.dicom

def get(url, dataset, filename=None):
    """ Perform a WADO request using a root URL and a dataset to provide the
        necessary information (Study Instance UID, Series Instance UID and
        SOP Instance UID).
        
        If filename is None, then the resulting dataset will be returned, otherwise
        it will be stored in the given filename. ::
        
            wado_url = "http://pacs.example.com/wado"
    
            query = medipy.io.dicom.DataSet()
            query.patient_id = "12345"
            query.study_instance_uid = "..."
            query.series_instance_uid = "..."
            query.sop_instance_uid = "..."
            
            # Store the result in a file
            medipy.network.dicom.wado.get(wado_url, query, "result.dcm")
            
            # Return the resulting dataset.
            dataset = medipy.network.dicom.wado.get(wado_url, query)
    """
    
    query = {
        "requestType" : "WADO", # Required, cf. 3.18-2011, 8.1.1,
        "studyUID" : dataset.study_instance_uid.value,
        "seriesUID" : dataset.series_instance_uid.value,
        "objectUID" : dataset.sop_instance_uid.value,
        "contentType" : "application/dicom"
    }
    
    query = "&".join(["{0}={1}".format(*item) for item in query.items()])
    
    dicom_fd = urllib.urlopen("{0}?{1}".format(url, query))
    if dicom_fd.getcode() != 200 :
        return None
    
    if filename is None :
        fd, temp_name = tempfile.mkstemp()
        fd = os.fdopen(fd, "wb")
    else :
        fd = open(filename, "wb")
    
    fd.write(dicom_fd.read())
    
    dicom_fd.close()
    fd.close()
    
    if filename is None :
        dataset = medipy.io.dicom.read(temp_name)
        os.remove(temp_name)
    
        return dataset
    else :
        return None
