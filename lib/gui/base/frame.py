##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import logging

import wx
import wx.xrc

import medipy.gui.base.base
import medipy.gui.utilities
import medipy.gui.xrc_wrapper


class Frame(medipy.gui.xrc_wrapper.Frame):
    """ Basic frame loaded from XRC resource
    """
    
    def __init__(self, xrc_file, frame_name, handlers, ui, controls,
                 parent=None, *args, **kwargs):
        
        self.ui = ui
        
        # Load XRC
        self._xml_resource, frame = medipy.gui.base.base.load_xrc(
            xrc_file, handlers, "Frame", parent, frame_name)
        
        # Build frame
        medipy.gui.xrc_wrapper.Frame.__init__(self, frame, *args, **kwargs)
        
        # Build UI object
        self.ui.from_window(self, controls)
        # Bind menu items events to event handlers
        self.bind_menu_items()
    
    def bind_menu_items(self):
        """ Bind the EVT_MENU event of all menu items to member functions with
            the (almost) same name. Menu items must use 
            underscore_separated_name, while function must be called 
            OnCamelCaseName.
        """
        
        for element in self._xml_resource.getElementsByTagName("object") :
            if element.getAttribute("class") != "wxMenuItem" :
                continue
            
            name = element.getAttribute("name")
                
            if name == "" :
                logging.warning("Menu item has invalid name \"%s\"", name)
                continue
                
            function_name = medipy.gui.utilities.underscore_to_camel_case(name)
            function_name = "On" + function_name
            
            item_id = wx.xrc.XRCID(name)
            
            # Bind to handler function
            if hasattr(self, function_name) :
                handler = getattr(self, function_name)
                self.Bind(wx.EVT_MENU, handler, id=item_id)
            else :
                logging.warning("No function called %s", function_name)