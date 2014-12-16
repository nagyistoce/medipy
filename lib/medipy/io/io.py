##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import logging
import os
import re
import urlparse
import sys

import numpy

import medipy.base

import schemes

def load_serie(url, dtype=numpy.single):
    """ This function is deprecated. medipy.io.load should be used instead.
    """
    
    logging.warning("\"medipy.io.load_serie\" is a deprecated function")
    return load(url, dtype)
    
def save_serie(images, url) :
    """ Save a serie of images.
        
          * image : image to save.
          * url : url to save to.
          
        See :func:`load` for URL details.
    """
    
    scheme, path, _ = _split(url)
    
    try :
        saver = getattr(scheme, "save_serie")
    except AttributeError :
        raise medipy.base.Exception("Scheme \"{0}\" cannot save files".format(scheme))

    saver(images, path)    

def load(url, dtype=numpy.single) :
    """ Load an image.
        
        * ``url`` : url to load from, uses the usual syntax of 
          ``[scheme "://"] [authority] path [ "#" fragment]``
        * ``dtype`` : type to which the data will be cast. Passing ``None`` will not cast.
        
        The URL ``scheme`` can be one of :
          
        * :mod:`~medipy.io.schemes.file` : load the image from the filesystem.
        * :mod:`~medipy.io.schemes.dicomdir` : load an image using a DICOMDIR
        * :mod:`~medipy.io.schemes.dicom` : load an image using a local 
          filesystem directory containing DICOM files
        * :mod:`dicom-series <medipy.io.schemes.dicom_series>` : load an image
          using a :class:`dicom_series <medipy.io.dicom.DicomSeries>` file
        
        Refer to the documentation of the different schemes for the details of 
        the URL syntax.
        
        If no scheme is specified, the URL is assumed to be a filesystem path. 
    """
    
    scheme, path, fragment = _split(url)
    
    try :
        loader = getattr(scheme, "load")
    except AttributeError :
        raise medipy.base.Exception("Scheme \"{0}\" cannot load files".format(scheme))

    images = loader(path, fragment)
    
    if not images:
        raise medipy.base.Exception("No such file or directory: {0}".format(url))
    
    for image in images:
        if dtype :
            image.data = image.data.astype(dtype)
        image.metadata.setdefault("loader", {})["url"] = url
    
    if len(images)>1:
        return images
    else:
        return images[0]
    
def save(image, url) :
    """ Save an image.
        
          * image : image to save.
          * url : url to save to.
        
        See :func:`load` for URL details.
    """
    
    scheme, path, _ = _split(url)
    
    try :
        saver = getattr(scheme, "save")
    except AttributeError :
        raise medipy.base.Exception("Scheme \"{0}\" cannot save files".format(scheme))

    saver(image, path)

def number_of_images(url):
    """ Return the number of images contained at given URL
    """
    
    scheme, path, fragment = _split(url)
    
    try :
        function = getattr(scheme, "number_of_images")
    except AttributeError :
        raise medipy.base.Exception("Scheme \"{0}\" cannot compute number of images".format(scheme))

    return function(path, fragment)

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
    # Note that "_" is not allowed in schemes (cf. http://tools.ietf.org/html/rfc3986#section-3.1)
    # and that "-" is not allowed in Python names : any "-" in the scheme is
    # transformed to "_" in the Python name
    scheme, netloc, path, query, fragment = urlparse.urlsplit(url)
    scheme = scheme.replace("-", "_")
    
    if sys.platform == "win32" and re.match(r"[a-zA-Z]:", path) :
        path = "/"+path
        slash_added = True
    url_without_scheme = ("", netloc, path, query, fragment)
    _, _, path, _, fragment = urlparse.urlsplit(urlparse.urlunsplit(url_without_scheme))
    
    scheme = scheme or "file"
    try :
        scheme = getattr(schemes, scheme)
    except AttributeError :
        raise medipy.base.Exception("Unknown scheme : \"{0}\"".format(scheme))
    
    if slash_added :
        path = path[1:]
    
    return scheme, path, fragment
