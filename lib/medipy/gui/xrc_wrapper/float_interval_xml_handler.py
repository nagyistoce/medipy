##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import wx
import wx.xrc

class FloatIntervalXMLHandler(wx.xrc.XmlResourceHandler):
    def __init__(self):
        wx.xrc.XmlResourceHandler.__init__(self)
        # No style for this window
    
    def CanHandle(self, node):
        return self.IsOfClass(node, "medipy.gui.control.FloatInterval")
    
    def DoCreateResource(self):
        # Late import to avoid circular dependencies
        from medipy.gui.control import FloatInterval
        
        assert self.GetInstance() is None
        widget = FloatInterval(
            parent = self.GetParentAsWindow(),
            id = self.GetID(),
            pos = self.GetPosition(),
            size = self.GetSize(),
            style = self.GetStyle(),
            name = self.GetName()
        )
        self.SetupWindow(widget)
        self.CreateChildren(widget)
        
        return widget