# -*- coding: latin-1 -*-

##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import os

import wx

import medipy.base
import medipy.gui.base
import medipy.gui.image
import medipy.gui.io

class MainFrame(medipy.gui.base.Frame):
    class UI(medipy.gui.base.UI):
        def __init__(self):
            self.panel = None
            self.image = None
            
            self.controls = ["panel"]
            
            medipy.gui.base.UI.__init__(self)
        
        def from_window(self, window, names):
            medipy.gui.base.UI.from_window(self, window, names)
            image = medipy.base.Image((256,256,256), value=0)
            self.image = medipy.gui.image.Image(
                self.panel, layers = [{"image" : image}], interpolation = True)
            sizer = wx.BoxSizer(wx.HORIZONTAL)
            sizer.Add(self.image, 1, wx.EXPAND)
            self.panel.SetSizer(sizer)
            self.panel.Layout()
    
    def __init__(self, parent=None, *args, **kwargs):
        self.ui = MainFrame.UI()
        self._save_location = None
        
        xrc_file = medipy.base.find_resource(
            os.path.join("resources", "viewer.xrc"))
        medipy.gui.base.Frame.__init__(self, xrc_file, "main_frame", 
            [], self.ui, self.ui.controls, parent, *args, **kwargs)
    
    def OnOpen(self, dummy):
        images = medipy.gui.io.load(self)
        if images :
            self.ui.image.delete_layer(0)
            self.ui.image.append_layer(images[0])
            self.ui.image.reset_view()
    
    def OnSave(self, event):
        if self._save_location is None :
            self.OnSaveAs(event)
        else :
            medipy.io.save(self.ui.image.get_layer_image(0), self._save_location)
    
    def OnSaveAs(self, dummy):
        result = medipy.gui.io.save(self.ui.image.get_layer_image(0), self)
        if result is not None and self._save_location is None :
            self._save_location = result
    
    def OnQuit(self, dummy):
        self.Close()
    
    def OnAbout(self, dummy):
        info = wx.AboutDialogInfo()
        info.AddDeveloper("Julien Lamy")
        info.AddDeveloper("Based on the MediPy framework (http://code.google.com/p/medipy/)")
        info.SetCopyright(u"Copyright (C) Université de Strasbourg, 2011")
        info.SetDescription("MediPy/iewer is a simple medical image viewer")
        info.SetLicense("Distributed under the terms of the CeCILL-B license, " 
                        "as published by the CEA-CNRS-INRIA.\n"
                        "Refer to the LICENSE file or to "
                        "http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html "
                        "for details.")
        info.SetName("MediPy/Viewer")
        info.SetWebSite("http://code.google.com/p/medipy/")
        
        wx.AboutBox(info)