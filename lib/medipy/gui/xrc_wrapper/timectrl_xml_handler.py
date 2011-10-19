import wx
import wx.xrc

from wx.lib.masked import TimeCtrl

class TimeCtrlXMLHandler(wx.xrc.XmlResourceHandler):
    def __init__(self):
        wx.xrc.XmlResourceHandler.__init__(self)
        # No style for this window
    
    def CanHandle(self, node):
        return self.IsOfClass(node, "wx.lib.masked.TimeCtrl")
    
    def DoCreateResource(self):
        assert self.GetInstance() is None
        control = TimeCtrl(
            parent = self.GetParentAsWindow(),
            id = self.GetID(),
            pos = self.GetPosition(),
            size = self.GetSize(),
            style = self.GetStyle(),
            name = self.GetName()
        )
        self.SetupWindow(control)
        self.CreateChildren(control)
        
        return control