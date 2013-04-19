##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import sys
import wx

class Importer(object):
    def __init__(self, with_dialog = True):
        self._with_dialog = with_dialog
        
        if self._with_dialog:
            self._dialog = wx.ProgressDialog("Loading Medipy",
                "Loading an.imported.module.name.CouldBeAsLongAsThis ...",
                style = wx.PD_APP_MODAL|wx.PD_AUTO_HIDE|wx.PD_SMOOTH)
            self._dialog.Update(0, 'Loading ...')
            self._dialog.Show()
            self.find_module = lambda name, path : self._find_module(name, path, True)
        else : 
            self.find_module = lambda name, path : self._find_module(name, path, False)
        self._count = 0
        
    def _find_module(self, fullname, path=None, update_dialog=False):
        self._count += 1
        if update_dialog :
            self._dialog.Pulse('Loading %s ...' % fullname)
    
    def end(self):
        if self._with_dialog:
            self._dialog.Hide()
            self._dialog.Destroy()
            
    count = property(lambda self : self._count)

def preimport():
    im = Importer()
    sys.meta_path.append(im)
    import imports
    sys.meta_path.remove(im)
    im.end()
    return im.count
