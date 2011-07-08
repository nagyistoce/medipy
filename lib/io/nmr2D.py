##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import os

from medipy.io.io_base import IOBase
import rbnmr

class Nmr2D(IOBase) :
    """I/O class for nmr2D format."""
    
    filenames = ["2rr", "2ir", "2ri", "2ii"]
    
    def __init__(self, filename=None, report_progress = None) :
        IOBase.__init__(self, filename, report_progress)
        self.__struct = None
    
    def can_load(self) :
        loadable=False

        if os.path.isfile(self.filename):
            loadable = True
        else :
            loadable = False
            
        if loadable :
            curdir = os.path.dirname(self.filename)
            if os.path.isfile(os.path.join(curdir, "proc2s")) and os.path.isfile(os.path.join(curdir, "procs")):
                pass
            else : 
                loadable = False
                
        return loadable
    
    def number_of_images(self) :
        return 1
    
    def load_data(self, index=0) :
        return self._struct["Data"]
    
    def load_metadata(self, index=0) :
        metadata = {}
        
        metadata["image_type"] = "spectroscopy"
        for key in self._struct:
            if key != "Data":
                metadata[key] = self._struct[key]
         
        return metadata
    
    def can_save(self):
        return False
    
    def _get__struct(self, index=0):
        if self.__struct is None :
            self.__struct = rbnmr.rbnmr(self._filename, self._report_progress)
        return self.__struct
    _struct = property(_get__struct)