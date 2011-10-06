##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import wx
class EvtHandler(wx.EvtHandler):
    """ This class takes over control of an already instantiated object (such as
        those created by wx.xrc.XmlResource(<filename>).LoadFrame(<framename>)).
        
        This code is taken from http://wiki.wxpython.org/UsingXmlResources
    """ 
    def __init__(self, other):
        wx.EvtHandler.__init__(self)
        self.this = other.this
        self.thisown = 1
        del other.this
        self._setOORInfo(self)
