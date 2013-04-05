# -*- coding: latin-1 -*-

##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
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
from medipy.gui.image.layers_panel import LayersPanel
import medipy.gui.io

class MainFrame(medipy.gui.base.Frame):
    class UI(medipy.gui.base.UI):
        def __init__(self):
            self.panel = None
            self.image = None
            self.layers_panel_container = None
            self.layers_panel = None
            
            self.controls = ["panel", "layers_panel_container"]
            
            medipy.gui.base.UI.__init__(self)
        
        def from_window(self, window, names):
            medipy.gui.base.UI.from_window(self, window, names)
            image = medipy.base.Image((256,256,256), value=0)
            self.image = medipy.gui.image.Image(
                self.panel, layers = [{"image" : image}], interpolation = False)
            sizer = wx.BoxSizer(wx.HORIZONTAL)
            sizer.Add(self.image, 1, wx.EXPAND)
            self.panel.SetSizer(sizer)
            
            self.layers_panel = LayersPanel(self.layers_panel_container)
            self.layers_panel.image = self.image
            
            self.layers_panel_container.GetSizer().Add(self.layers_panel, 1, wx.EXPAND)
            
            window.Layout()
    
    def __init__(self, parent=None, *args, **kwargs):
        self._title = "MediPy/Viewer"
        self.ui = MainFrame.UI()
        self._save_location = None
        
        xrc_file = medipy.base.find_resource(
            os.path.join("resources", "viewer.xrc"))
        medipy.gui.base.Frame.__init__(self, xrc_file, "main_frame", 
            [], self.ui, self.ui.controls, parent, *args, **kwargs)
    
    def OnOpen(self, dummy):
        images = medipy.gui.io.load(self)
        if images :
            image = images[0]
            
            self.ui.image.delete_layer(0)
            self.ui.image.insert_layer(0, image)
            self.ui.image.reset_view()
            
            if "loader" in image.metadata :
                url = image.metadata["loader"]["url"]
            else :
                url = None
            
            if url :
                self.SetTitle("{0} ({1})".format(self._title, url))
            else :
                self.SetTitle(self._title)
    
    def OnSave(self, event):
        if self._save_location is None :
            self.OnSaveAs(event)
        else :
            medipy.io.save(self.ui.image.get_layer_image(0), self._save_location)
    
    def OnSaveAs(self, dummy):
        result = medipy.gui.io.save(self.ui.image.get_layer_image(0), self)
        if result is not None and self._save_location is None :
            self._save_location = result
            self.SetTitle("{0} ({1})".format(self._title, result))
    
    def OnQuit(self, dummy):
        self.Close()
    
    def OnResetView(self, dummy):
        self.ui.image.reset_view()
        self.ui.image.render()
    
    def OnMultiplanar(self, dummy):
        self.ui.image.slice_mode = "multiplanar"
        self.ui.image.render()
    
    def OnAxial(self, dummy):
        self.ui.image.slice_mode = "axial"
        self.ui.image.render()
    
    def OnCoronal(self, dummy):
        self.ui.image.slice_mode = "coronal"
        self.ui.image.render()
    
    def OnSagittal(self, dummy):
        self.ui.image.slice_mode = "sagittal"
        self.ui.image.render()
    
    def OnInterpolation(self, event):
        self.ui.image.interpolation = event.IsChecked()
        self.ui.image.render()
    
    def OnAbout(self, dummy):
        info = wx.AboutDialogInfo()
        info.AddDeveloper("LINC-IPB")
        info.AddDeveloper("Based on the MediPy framework (http://code.google.com/p/medipy/)")
        info.SetCopyright(u"Copyright (C) Université de Strasbourg, 2011-2012")
        info.SetDescription("MediPy/Viewer is a simple medical image viewer")
        info.SetLicense("Distributed under the terms of the CeCILL-B license, " 
                        "as published by the CEA-CNRS-INRIA.\n"
                        "Refer to the LICENSE file or to "
                        "http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html "
                        "for details.")
        info.SetName("MediPy/Viewer")
        info.SetWebSite("http://code.google.com/p/medipy/")
        
        wx.AboutBox(info)