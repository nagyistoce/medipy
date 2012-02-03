##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import os
import re
import urlparse
import sys

import numpy

import schemes

def load(url, dtype=numpy.single) :
    """ Load an image.
        
        url : url to load from.
        dtype : type to which the data will be cast. Passing None will not cast.
        
        url uses the usual syntax of [scheme] "://" [authority] path [ "#" fragment]
        scheme can be one of :
          * "file" : the default value if no scheme is specified. Load the image
            from the filesystem.
          * "dicomdir" : load an image using a DICOMDIR
          * "dicom" : load an image using a local filesystem directory 
            containing DICOM files
        
        >>> import medipy.io
        >>> medipy.io.load("/some/where/image.nii") # Uses "file" scheme
        >>> medipy.io.load("dicomdir:/some/where/DICOMDIR#series_instance_uid=1.2.3.4")
        >>> medipy.io.load("dicom:/some/where/#series_instance_uid=1.2.3.4")
        
        Refer to the different schemes for the details of the URL syntax
    """
    
    scheme, path, fragment = _split(url)
    
    try :
        loader = getattr(scheme, "load")
    except AttributeError :
        raise Exception("Scheme \"{0}\" cannot load files".format(scheme))

    image = loader(path, fragment)
    
    if dtype :
        image.data = image.data.astype(dtype)
    
    image.metadata.setdefault("loader", {})["url"] = url
    
    return image    

def save(image, url) :
    """ Save an image.
        
        image : image to save.
        url : url to save to.
    """
    
    scheme, path, _ = _split(url)
    
    try :
        saver = getattr(scheme, "save")
    except AttributeError :
        raise Exception("Scheme \"{0}\" cannot save files".format(scheme))

    saver(image, path)

def number_of_images(url):
    """ Return the number of images contained at given URL
    """
    
    scheme, path, fragment = _split(url)
    
    try :
        function = getattr(scheme, "save")
    except AttributeError :
        raise Exception("Scheme \"{0}\" cannot compute number of images".format(scheme))

    function(path, fragment)

def _split(url):
    """ Return the scheme (as a Python module), the path and the fragment from
        given url.
    """

    slash_added = False    
    if sys.platform == "win32" and re.match(r"[a-zA-Z]:", url) :
        slash_added = True
        url = "/"+url
    
    # Parse the URL : first get the scheme, then parse the URL without scheme
    # to deal with fragment. Doing this avoids modifying urlparse.uses_fragment
    scheme, netloc, path, query, fragment = urlparse.urlsplit(url)
    if sys.platform == "win32" and re.match(r"[a-zA-Z]:", path) :
        path = "/"+path
        slash_added = True
    url_without_scheme = (scheme, netloc, path, query, fragment)
    _, _, path, _, fragment = urlparse.urlsplit(urlparse.urlunsplit(url_without_scheme))
    
    scheme = scheme or "file"
    try :
        scheme = getattr(schemes, scheme)
    except AttributeError :
        raise Exception("Unknown scheme : \"{0}\"".format(scheme))
    
    if slash_added :
        path = path[1:]
    
    return scheme, path, fragment
