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
    """ Load images.
    """
    
    scheme_path, scheme_fragment = _get_real_path_and_fragment(path, fragment)
    return scheme.load(scheme_path, scheme_fragment)

def number_of_images(path, fragment) :
    """ Return the number of images.
    """
    
    scheme_path, scheme_fragment = _get_real_path_and_fragment(path, fragment)
    return scheme.number_of_images(url)

def _get_real_path_and_fragment(path, fragment):
    if not os.path.isfile(path) :
        raise medipy.base.Exception("No such file : {0}".format(repr(path)))
    
    tag, value = fragment.split("=")
    
    if tag not in ["uid", "description", "custom_name"] :
        raise medipy.base.Exception("Unknown tag : {0}".format(repr(tag)))
    
    # Read root and series informations
    
    fd = open(path)
    root = fd.readline().strip()
    scheme, scheme_path = root.split(":", 1)
    scheme = getattr(medipy.io.schemes, scheme)
    
    reader = csv.reader(fd)
    series = [x for x in reader]
    
    fd.close()
    
    # Look for a matching fragment
    scheme_fragment = None
    for serie in series :
        uid = serie[0]
        
        if tag == "uid" and uid == value :
            scheme_fragment = "series_instance_uid={0}".format(uid)
            break
        elif tag == "description" and len(serie)>1 and serie[1] == value :
            scheme_fragment = "series_instance_uid={0}".format(uid)
        elif tag == "custom_name" and len(serie)>2 and serie[2] == value :
            scheme_fragment = "series_instance_uid={0}".format(uid)
    
    if scheme_fragment is None :
        raise medipy.base.Exception("No serie matching {0} in {1}".format(repr(fragment), repr(path)))
    
    return scheme_path, scheme_fragment
