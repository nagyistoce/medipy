##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import wx
from window import Window

class Panel(wx.Panel, Window):
    def __init__(self, other, 
                 id=wx.ID_ANY, pos=None, size=None, name=None):
        Window.__init__(self, other, id, pos, size, name)
