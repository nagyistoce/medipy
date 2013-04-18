##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import logging
import xml.dom.minidom

import wx.xrc

class UI(object):
    """ Simple object holding wx controls
    """
    
    def from_window(self, window, names):
        """ Set members from controls taken from wx window
        """
        
        for name in names :
            if not hasattr(self, name) :
                logging.warning("%s has no attribute \"%s\"", 
                                type(self), name)
            value = wx.xrc.XRCCTRL(window, name)
            if value is None :
                logging.warning("%s has no control \"%s\"", 
                                type(window), name)
            else :
                setattr(self, name, value)

def load_xrc(xrc_file, handlers, class_name, parent, window_name):
    """ Load an XRC file with given handlers, and load a window from that file
        Return the xml document and the window
    """
    
    xml_document = xml.dom.minidom.parse(xrc_file)
    resource = wx.xrc.EmptyXmlResource()
    for handler in handlers :
        resource.InsertHandler(handler)
    # Don't use LoadFromString so that the internal wxFileSystem is set to
    # the directory holding xrc_file
    resource.Load(xrc_file)
    
    # Build frame
    loader = getattr(resource, "Load%s"%class_name)
    if class_name in ["Bitmap", "Icon", "Menu"] :
        window = loader(window_name)
    else :
        window = loader(parent, window_name)
    
    return (xml_document, window) 
