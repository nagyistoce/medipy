##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import wx
import wx.xrc

from medipy.gui.viewer_3d import Viewer3D

class Viewer3DXMLHandler(wx.xrc.XmlResourceHandler):
    def __init__(self):
        wx.xrc.XmlResourceHandler.__init__(self)
        # No style for this window
    
    def CanHandle(self, node):
        return self.IsOfClass(node, "medipy.gui.viewer_3d.Viewer3D")
    
    def DoCreateResource(self):
        assert self.GetInstance() is None
        viewer_3d = Viewer3D(
            parent = self.GetParentAsWindow(),
            id = self.GetID(),
            pos = self.GetPosition(),
            size = self.GetSize(),
            style = self.GetStyle(),
            name = self.GetName()
        )
        
        self.SetupWindow(viewer_3d)
        self.CreateChildren(viewer_3d)

        return viewer_3d
