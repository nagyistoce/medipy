##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

""" Load images from a dicom_series file.
    
    The path is interpreted as a local filesystem path to a dicom_series file.
    
    The fragment is of the form ``tag "=" value`` where tag can be one of :
      
      * ``uid``
      * ``description``
      * ``custom_name``
    
    The following examples are correct fragments : ::
    
        uid_fragment = "uid=2.25.169876028567313845957086020100076112905"
        description_fragment = "description=T1 3D SPGR"
        custom_name_fragment = "custom_name=my specific description"
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
        return 0
    else :
        return 1
