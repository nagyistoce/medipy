##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import wx
from wx import GetTranslation as _

import medipy.io
import medipy.io.schemes.file

class ImageFileDialog(wx.FileDialog):
    """ Custom version of wx.FileDialog tuned for images. It loads a default
        wildcard corresponding to all known images and a default path from
    """
    
    _last_path = None
    _wildcard = None
    
    def __init__(self, *args, **kwargs):
        wx.FileDialog.__init__(self, *args, **kwargs)
        self._config = wx.Config("MediPy")
        
        if self._last_path is None :
            self._setup_default_directory()
    
        if self._wildcard is None :
            self._setup_wildcard()
        
        self.SetWildcard(self._wildcard)
        self.SetDirectory(self._last_path)
    
    def ShowModal(self):
        """ Shows the dialog, returning wxID_OK if the user pressed OK, and 
            wxID_CANCEL otherwise. If the return code was wx.ID_OK, then the
            directory will be saved
        """
        
        return_code = wx.FileDialog.ShowModal(self)
        
        if return_code == wx.ID_OK :
            self._last_path = self.GetDirectory()
            self._config.Write("ImageFileDialog/DefaultPath", 
                               self.GetDirectory())
            self._config.Flush()
        
        return return_code
    
    def _setup_default_directory(self):
        self._last_path = self._config.Read("ImageFileDialog/DefaultPath")
    
    def _setup_wildcard(self):
        
        self._wildcard = []
        
        all_filenames = set()
        for io_class in medipy.io.schemes.file.io_classes :
            for filename in io_class.filenames :
                all_filenames.add(filename)
        all_filenames.add("DICOMDIR")
        all_filenames.add("dicomdir")
        
        # All known images
        self._wildcard += [_("All known images"), ";".join(all_filenames)]
        # One entry for each filename
        for filename in all_filenames :
            self._wildcard += [filename, filename]
        # DICOMDIR
        self._wildcard += ["DICOMDIR", "DICOMDIR;dicomdir"]
        # Everything else
        self._wildcard += [_("All"), "*"]
        
        self._wildcard = "|".join(self._wildcard)
