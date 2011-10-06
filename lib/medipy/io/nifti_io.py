##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import fnmatch
import logging

import nifti
import numpy

# The following two tests ease the translation from earlier versions of PyNifti
# (i.e. with modules nifti.niftiimage and nifti.niftiformat) to newer 
# versions (i.e. with modules nifti.image and nifti.format)

if "image" in dir(nifti) :
    from nifti import NiftiImage
else :
    from nifti.niftiimage import NiftiImage

if "format" in dir(nifti) :
    from nifti.format import NiftiFormat
else : 
    from nifti.niftiformat import NiftiFormat

from medipy.base import ObservableList
from io_base import IOBase

class Nifti(IOBase) :
    """I/O class for Nifti format.
    
    Uses pynifti library."""
    
    filenames = ["*.nii", "*.hdr", "*.nii.gz"]
    
    def __init__(self, *args, **kwargs):
        IOBase.__init__(self, *args, **kwargs)
        self._filename = str(self.filename)
    
    def can_load(self) :
        loadable=False
        try :
            NiftiFormat(self._filename)
            loadable = True
        except :
            pass
        return loadable
    
    def number_of_images(self) :
        return 1
    
    def load_data(self, index=0) :
        # pynifti does not appear to have a hook in nifti.clib.nifti_image_read,
        # nor is it possible to feed it a buffer filled with the content of the
        # file
        if self._report_progress is not None :
            self._report_progress(0.)
            
        image = NiftiImage(self._filename)
        data = image.data
        
        if self._report_progress is not None :
            self._report_progress(1.)
        
        while data.shape[0] == 1 and len(data.shape)>3 :
            data = data.reshape(data.shape[1:])
        return data
    
    def load_metadata(self, index=0) :
        format = NiftiFormat(self._filename)
        
        metadata = {}
        
        # Get the number of dimensions in the image from the extent
        nb_dimensions = len(format.extent)
        for i in range(max(len(format.extent)-3, 0)) :
            if format.extent[-1] == 1 : 
                nb_dimensions -= 1
            else : 
                break
        # Get spacing and origin for those dimensions. Reverse to keep the order
        # of the numpy array
        metadata["spacing"] = format.asDict()["pixdim"][1:nb_dimensions+1]
        metadata["spacing"].reverse()
        metadata["spacing"] = numpy.asarray(metadata["spacing"])
        
        if format.asDict()["scl_slope"] != 0 :
            metadata["slope"] = format.asDict()["scl_slope"] 
            metadata["shift"] = format.asDict()["scl_inter"]
        
        metadata["header"] = format.asDict()
        metadata["annotations"] = ObservableList()
        
        return metadata
    
    def can_save(self, image) :
        found = False
        for pattern in self.filenames :
            if fnmatch.fnmatch(self.filename, pattern):
                found = True
                break
        return found
    
    def save(self, image) :
        nifti_image = NiftiImage(image.data)
        spacing = image.spacing.tolist()
        spacing.reverse()
        nifti_image.pixdim = spacing 
#        logging.warning("Image direction and origin will not be saved")
        nifti_image.save(self._filename)

if __name__ == "__main__" :
    import sys
    
    loader = Nifti(sys.argv[1])
    metadata = loader.load_metadata()
    #del metadata["header"]
    print metadata