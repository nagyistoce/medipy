##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import wx
from top_level_window import TopLevelWindow

class Dialog(wx.Dialog, TopLevelWindow):
    def __init__(self, other, 
                 id=wx.ID_ANY, title=None, pos=None, size=None, name=None):
        TopLevelWindow.__init__(self, other, id, pos, size, name)
        if title : 
            self.SetTitle(title)
