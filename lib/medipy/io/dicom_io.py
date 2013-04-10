##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import medipy.base
import medipy.io.dicom
import medipy.io.dicom.normalize
from medipy.io.io_base import IOBase

class Dicom(IOBase) :
    
    filenames = ["*.DCM", "*.dcm"]
    
    def __init__(self, filename=None, report_progress=None) :
        IOBase.__init__(self, filename, report_progress)
        self._stacks = None
        self._image = None
        self._index = None
        self._set_filename(filename)
    
    def can_load(self):
        if self._filename is not None and self._stacks is None :
            try :
                dataset = medipy.io.dicom.read(self._filename)
            except medipy.base.Exception :
                return False
            else :
                datasets = medipy.io.dicom.normalize.normalize(dataset)
                if isinstance(datasets, medipy.io.dicom.DataSet) :
                    datasets = [datasets]
                self._stacks = medipy.io.dicom.stacks(datasets)
                return True
        else :
            return False
    
    def number_of_images(self) :
        return len(self._stacks)
    
    def load_data(self, index=0):
        if self._index != index :
            self._image = medipy.io.dicom.image(self._stacks[index])
            self._index = index
        return self._image.data
    
    def load_metadata(self, index=0):
        if self._index != index :
            self._image = medipy.io.dicom.image(self._stacks[index])
            self._index = index
        metadata = self._image.metadata
        metadata["direction"] = self._image.direction
        metadata["origin"] = self._image.origin
        metadata["spacing"] = self._image.spacing
        return metadata
    
    def can_save(self,image):
        return False
    
    ##############
    # Properties #
    ##############
    
    def _set_filename(self, filename):
        self._filename = filename
        self._stacks = None
        self._image = None
        self._index = None
