##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

""" Load images from a dicom_series file.
    
    The path is interpreted as a local filesystem path to a dicom_series file.
    
    The fragment is of the form tag "=" value where tag can be one of :
      * uid
      * description
      * custom_name
"""

import csv
import os

import medipy.base
import medipy.io

def load(path, fragment) :
    """ Load an image.
    """
    
    # Check arguments
    
    if not os.path.isfile(path) :
        raise medipy.base.Exception("No such file : {0}".format(repr(path)))
    
    tag, value = fragment.split("=")
    
    if tag not in ["uid", "description", "custom_name"] :
        raise medipy.base.Exception("Unknown tag : {0}".format(repr(tag)))
    
    # Read root and series informations
    
    fd = open(path)
    root = fd.readline().strip()
    
    reader = csv.reader(fd)
    series = [x for x in reader]
    
    fd.close()
    
    # Look for a matching URL
    
    url = None
    for serie in series :
        uid = serie[0]
        
        if tag == "uid" and uid == value :
            url = "{0}#series_instance_uid={1}".format(root, uid)
            break
        elif tag == "description" and len(serie)>1 and serie[1] == value :
            url = "{0}#series_instance_uid={1}".format(root, uid)
        elif tag == "custom_name" and len(serie)>2 and serie[2] == value :
            url = "{0}#series_instance_uid={1}".format(root, uid)
    
    if url is None :
        raise medipy.base.Exception("No serie matching {0} in {1}".format(repr(fragment), repr(path)))
    
    return medipy.io.load(url)
    
def number_of_images(path, fragment) :
    """ Return the number of images.
    """
    
    pass

if __name__ == "__main__" :
    load("/base_image/sep/choline/01-001/Original/Exam01/dicom_series", "uid=1.3.12.2.1107.5.2.30.25842.30000010043011390948400000768")
    load("/base_image/sep/choline/01-001/Original/Exam01/dicom_series", "description=T1 MPR3D TRA GADO")
    load("/base_image/sep/choline/01-001/Original/Exam01/dicom_series", "custom_name=T1 3D gado")