##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
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

def get(url, dataset, filename=None):
    """ Perform a WADO request using a root URL and a dataset to provide the
        necessary information (Study Instance UID, Series Instance UID and
        SOP Instance UID).
        
        If filename is None, then the resulting dataset will be returned, otherwise
        it will be stored in the given filename.
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
    
    if filename is None :
        fd, temp_name = tempfile.mkstemp()
        fd = os.fdopen(fd, "wb")
    else :
        fd = open(filename, "wb")
    
    fd.write(dicom_fd.read())
    
    dicom_fd.close()
    fd.close()
    
    if filename is None :
        dataset = medipy.io.dicom.parse(temp_name)
        os.remove(temp_name)
    
        return dataset
    else :
        return None
