##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import wx
from layers_panel import LayersPanel

class LayersDialog(wx.Dialog):
    """ Display controls that allow to modify the layers of an image. 
        
        The user may add or delete layers, and change the colors and the 
        window/level of the image.
    """
    
    def __init__(self, parent, wx_id=wx.ID_ANY, title="",
                 pos=wx.DefaultPosition, size=wx.DefaultSize,
                 style=wx.DEFAULT_DIALOG_STYLE, name="layers_dialog"):
        
        super(LayersDialog, self).__init__(parent, wx_id, title, pos, size, style, name)
        
        self._layers_panel = LayersPanel(self)
        sizer = wx.BoxSizer()
        sizer.Add(self._layers_panel, 1, wx.EXPAND | wx.ALL, 5)
        self.SetSizer(sizer)
        sizer.SetSizeHints(self)
    
    ##############
    # Properties #
    ##############
    
    def _get_image(self):
        return self._layers_panel.image
    
    def _set_image(self, image):
        self._layers_panel.image = image
    
    image = property(_get_image, _set_image)