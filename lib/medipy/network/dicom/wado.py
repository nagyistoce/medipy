##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import os
import tempfile
import urllib
import urlparse

import medipy.io.dicom

def get(url, dataset):
    """ Perform a WADO request using a root URL and a dataset to provide the
        necessary information (Study Instance UID, Series Instance UID and
        SOP Instance UID).
    """
    
    query = {
        "requestType" : "WADO", # Required, cf. 3.18-2011, 8.1.1,
        "studyUID" : dataset.study_instance_uid,
        "seriesUID" : dataset.series_instance_uid,
        "objectUID" : dataset.sop_instance_uid,
        "contentType" : "application/dicom"
    }
    
    query = "&".join(["{0}={1}".format(*item) for item in query.items()])
    
    dicom_fd = urllib.urlopen("{0}?{1}".format(url, query))
    if dicom_fd.getcode() != 200 :
        return None
    
    temp_fd, temp_name = tempfile.mkstemp()
    os.write(temp_fd, dicom_fd.read())
    
    dicom_fd.close()
    os.close(temp_fd)
    
    dataset = medipy.io.dicom.parse(temp_name)
    os.remove(temp_name)
    
    return dataset