##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import wx
import wx.xrc

class DICOMExplorerXMLHandler(wx.xrc.XmlResourceHandler):
    def __init__(self):
        wx.xrc.XmlResourceHandler.__init__(self)
        # No style for this window
    
    def CanHandle(self, node):
        return self.IsOfClass(node, "medipy.gui.dicom.Explorer")
    
    def DoCreateResource(self):
        from medipy.gui.dicom import Explorer
        assert self.GetInstance() is None
        control = Explorer(
            parent = self.GetParentAsWindow(),
            ID = self.GetID(),
            pos = self.GetPosition(),
            size = self.GetSize(),
            style = self.GetStyle(),
            name = self.GetName()
        )
        self.SetupWindow(control)
        self.CreateChildren(control)
        
        return control
