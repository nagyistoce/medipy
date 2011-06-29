##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import logging
import optparse
import sys
import warnings

import wx

class Application(wx.PySimpleApp):
    def __init__(self):
        
        levels = {"debug": logging.DEBUG,
                  "info": logging.INFO,
                  "warning": logging.WARNING,
                  "error": logging.ERROR,
                  "critical": logging.CRITICAL}
        
        parser = optparse.OptionParser()
        parser.add_option("-d", "--debug", dest="debug",
                          help="Set the logging level "
                               "(debug, info, warning, error, or critical", 
                          metavar="LEVEL")
        parser.add_option("-r", "--redirect-to-console", dest="redirect", 
                          action="store_false", default=True,
                          help="Redirect all messages to the console")
        options = parser.parse_args()[0]
        
        if options.debug is not None :
            level = levels.get(options.debug, None)
            if level is None :
                raise Exception("Warning : unknown logging level %s" % options.debug)
            
            logging_format = ("[%(asctime)s] {%(pathname)s:%(lineno)d} "
                              "%(levelname)s - %(message)s")
            date_format = "%Y-%m-%d %H:%M:%S"
            
            logging.basicConfig(level=level, format=logging_format,
                                datefmt=date_format)
            
            if level != logging.DEBUG :
                # Disable deprecation warnings
                warnings.simplefilter("ignore", DeprecationWarning)
        else :
            logging.basicConfig(level=logging.ERROR)
            # Disable deprecation warnings
            warnings.simplefilter("ignore", DeprecationWarning)
        
        wx.PySimpleApp.__init__(self, redirect=options.redirect)
    
        if sys.platform.startswith("win") :
            import vtk
            vtk.vtkObject.GlobalWarningDisplayOff()