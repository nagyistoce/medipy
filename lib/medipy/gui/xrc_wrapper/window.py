##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import wx
from evt_handler import EvtHandler

class Window(wx.Window, EvtHandler):
    """ XRC Wrapper around wx.Window. Position, size, and name can specified, 
        overriding values defined in the XRC.
    """
    
    def __init__(self, other, id=wx.ID_ANY, pos=None, size=None, name=None) :
        EvtHandler.__init__(self, other)
        
        self.SetId(id)
        if pos : 
            self.SetPos(pos)
        if size : 
            self.SetSize(size)
        if name : 
            self.SetName(name)
