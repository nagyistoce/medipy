##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import ConfigParser
import imp
import new
import os
import sys

class Importer(object):
    
    def __init__(self):
        # The plugins search path is composed of
        #   * the environment variable MEDIPY_PLUGINS_PATH
        #   * on Linux, the plugins/path entry of /etc/medipy.cfg
        # When a plugin is requested, the different elements of the path are
        # searched in order. The first element where the plugin is found is
        # used.
        self._plugins_path = ""
        if "MEDIPY_PLUGINS_PATH" in os.environ :
            self._plugins_path = [
                x for x in os.environ["MEDIPY_PLUGINS_PATH"].split(os.pathsep)
                if os.path.isdir(x)
            ]
        global_path_getter = getattr(
            self, "_get_global_path_{0}".format(sys.platform), None)
        if global_path_getter :
            self._plugins_path.extend(global_path_getter())
    
    def find_module(self, fullname, path=None):
        if not fullname.startswith("medipy.") :
            # Not a medipy module, use default load
            return
        plugin = fullname.split(".")[1]
        for directory in self._plugins_path :
            if plugin in os.listdir(directory) :
                return self
    
    def load_module(self, fullname):
        try :
            path = self._plugins_path
            for level, package in enumerate(fullname.split(".")[1:]) :
                file, pathname, description = imp.find_module(package, path)
                path = [pathname]
                module = imp.load_module(".".join(fullname.split(".")[1:level+2]), file, pathname, description)
        except :
            raise
        else :
            return module
    
    #####################
    # Private interface #
    #####################
    
    def _get_global_path_linux2(self):
        """ The the global plugins path on Linux platform, based on the contents
            of /etc/medipy.cfg
        """
        
        global_plugins_path = []
        config_file = "/etc/medipy.cfg"
        if os.path.isfile(config_file) :
            config = ConfigParser.ConfigParser()
            config.read(config_file)
            if config.has_section("plugins") and config.has_option("plugins", "path") :
                global_plugins_path = [
                    x for x in config.get("plugins", "path").split(os.pathsep)
                    if os.path.isdir(x)
                ]
        
        return global_plugins_path
        
sys.meta_path.append(Importer())
