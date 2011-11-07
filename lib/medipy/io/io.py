##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import numpy

from medipy.base import Image

from medipy.io.ipb import IPB
from medipy.io.itk_io import ITK
from medipy.io.nifti_io import Nifti
from medipy.io.nmr2D import Nmr2D
from medipy.io.wx_image import WXImage
# Nifti is quite verbose when testing if an image can be loaded, so let's test
# it last
io_classes = [ITK, IPB, WXImage, Nmr2D, Nifti]

def get_loader(filename, report_progress=None) :
    """Search for a loader in io_classes"""
    for loader_class in io_classes : 
        loader = loader_class(filename, report_progress=report_progress)
        if loader.can_load() :
            return loader
    
    # If we get here, no loader was found
    raise Exception("No loader available for %s"%filename)

def get_saver(image, filename) :
    """Search for a saver in io_classes"""
    for saver_class in io_classes : 
        saver = saver_class(filename)
        if saver.can_save(image) :
            return saver

def number_of_images(filename, loader_class=None, loader=None) :
    """Return the number of images contained in the given filename"""
    if loader_class is not None and loader is not None :
        raise Exception("Cannot specify both loader_class and loader")
    
    if loader_class is None and loader is None:
        loader = get_loader(filename)
    elif loader_class is not None and loader is None:
        loader = loader_class(filename)
    # else loader_class is None and loader is not None : do nothing
    
    return loader.number_of_images()
    

def load(filename, index=0, dtype=numpy.single, rescale_data=True, 
         loader_class=None, loader=None, report_progress=None) :
    """ Load an image from a file.
        
        filename : the name of the file to load from.
        index : in the case of formats where multiple images can be stored in a
            file, this is the index of the image to load.
        type : the type to which the data will be cast. Pass None to keep the 
            original type.
        rescale_data : if the file contains a slope and/or a shift, this flag
            indicates whether or not to apply the transformation.
        loader_class : specific loader class to use. If None, it is automatically
            determined.
        report_progress : a unary function taking a float between 0 and 1 to
            report the loading progress
    """
    
    if loader_class is not None and loader is not None :
        raise Exception("Cannot specify both loader_class and loader")
    
    if loader_class is None and loader is None:
        loader = get_loader(filename, report_progress=report_progress)
    elif loader_class is not None and loader is None:
        loader = loader_class(filename, report_progress=report_progress)
    # else loader_class is None and loader is not None : do nothing
    
    max_index = number_of_images(filename, loader=loader)
    if index >= max_index :
        raise Exception("Cannot access image %i. File \"%s\" contains %i images"%(index, filename, max_index))
    
    data = loader.load_data(index)
    metadata = loader.load_metadata(index)
    metadata["loader"] = {
        "filename" : filename,
        "index" : index,
        "dtype" : dtype,
        "rescale_data" : rescale_data,
        "loader" : loader
    }
    
    # Convert the buffer to float
    original_type = data.dtype
    
    if rescale_data :
        if data.dtype != numpy.single :
            data = data.astype(numpy.single)
        # Scale if necessary
        if metadata.has_key("slope") :
            data *= metadata["slope"]
        
        # Shift if necessary
        if metadata.has_key("shift") :
            data += metadata["shift"]
    
    if dtype is not None : 
        if dtype != data.dtype:
            data = data.astype(dtype)
    elif original_type != data.dtype : 
        data = data.astype(original_type)
    
    args = {}
    for name in ["direction", "origin", "spacing", "data_type", "image_type", "annotations"] :
        if name in metadata :
            args[name] = metadata[name]
            del metadata[name]
    args["metadata"] = metadata
    
    return Image(data=data, **args)

def save(image, filename, saver_class=None) :
    """ Save an image to a file
    """
    if saver_class is None :
        saver = get_saver(image, filename)
    else :
        saver = saver_class(filename)
    saver.save(image)
