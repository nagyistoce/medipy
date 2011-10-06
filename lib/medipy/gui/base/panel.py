##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import medipy.gui.base.base
import medipy.gui.utilities
import medipy.gui.xrc_wrapper


class Panel(medipy.gui.xrc_wrapper.Panel):
    """ Basic panel loaded from XRC resource
    """
    
    def __init__(self, xrc_file, panel_name, handlers, ui, controls,
                 parent=None, *args, **kwargs):
        
        self.ui = ui
        
        # Load XRC
        self._xml_resource, frame = medipy.gui.base.base.load_xrc(
            xrc_file, handlers, "Panel", parent, panel_name)
        
        # Build dialog
        medipy.gui.xrc_wrapper.Panel.__init__(self, frame, *args, **kwargs)
        
        # Build UI object
        self.ui.from_window(self, controls)