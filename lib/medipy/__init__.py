##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# If uuid module is loaded too late, TLS might be uninitialized and calls to
# uuid functions will cause a segfault on some Linux versions. See
# http://www.sourceware.org/bugzilla/show_bug.cgi?id=12453
import uuid

# If expat is loaded too late, WrapITK causes crashes when an exception is 
# raised. See https://issues.itk.org/jira/browse/HISTITK-834
import xml.parsers.expat

import ConfigParser
import imp
import marshal
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
        self.plugins_path = []
        global_path_getter = getattr(
            self, "_get_global_path_{0}".format(sys.platform), None)
        if global_path_getter :
            self.plugins_path = global_path_getter()
        if "MEDIPY_PLUGINS_PATH" in os.environ :
            self.plugins_path = [
                x for x in os.environ["MEDIPY_PLUGINS_PATH"].split(os.pathsep)
                if os.path.isdir(x)
            ]
    
    def find_module(self, fullname, path=None):
        # Find the items from ``path`` which are rooted in self.plugins_path
        plugins_path = []
        for p in (path or []) :
            ancestors = [os.path.commonprefix([p, x])
                         for x in self.plugins_path]
            plugins_path.extend([x for x in ancestors 
                                 if x in self.plugins_path])
        
        # If no such path exist and fullname is not a top-level medipy plugin,
        # it is not our job to load it
        if len(plugins_path)==0 and not fullname.startswith("medipy.") :
            return None
        
        # Check if we can find fullname in our plugins search paths.
        plugin_path = self._get_search_path(fullname, 
                                            plugins_path or self.plugins_path)
        # We manage fullname iff we cound find such a path
        return self if plugin_path is not None else None 
    
    def load_module(self, fullname):
        module = sys.modules.setdefault(fullname, imp.new_module(fullname))
        
        path = os.path.join(self._get_search_path(fullname), 
                            *fullname.split(".")[1:])
        if os.path.isfile(path+".so") :
            module = imp.load_dynamic(fullname, path+".so")
        else :
            code, filename = self._get_code(fullname)
            
            if self._is_package(fullname) :
                module.__path__ = []
                module.__file__ = filename
                module.__package__ = fullname
            else :
                module.__file__ = filename
                module.__package__ = fullname.rpartition(".")[0]
            
            exec(code, module.__dict__)
        
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
    
    def _get_search_path(self, fullname, plugins_path=None):
        """ Return the first path from ``plugins_path`` (or ``self.plugins_path``
            if ``plugins_path`` is None) that contains the plugin ``fullname``.
            If no such path exists, return None.
        """
        
        if plugins_path is None :
            plugins_path = self.plugins_path
        
        path = None
        for directory in plugins_path :
            plugin_path = os.path.join(directory, *fullname.split(".")[1:])
            candidates = []
            if os.path.isdir(plugin_path) :
                candidates = [os.path.join(plugin_path, "__init__"+suffix)
                              for suffix in [".py", ".pyc", ".so"]]
            else :
                candidates = [plugin_path+suffix for suffix in [".py", ".pyc", ".so"]]
            if any(os.path.isfile(x) for x in candidates) :
                path = directory
                break
        return path
    
    def _get_code(self, fullname):
        """ Return the code object associated with the plugin ``fullname``.
        """
        
        path = os.path.join(self._get_search_path(fullname), 
                            *fullname.split(".")[1:])
        
        if os.path.isdir(path) :
            path = os.path.join(path, "__init__")
        
        if os.path.isfile(path+".py") :
            path += ".py"
            source = True
        elif os.path.isfile(path+".pyc") :
            path += ".pyc"
            source = False
        
        code = None
        
        if source :
            with open(path, "rU") as f :
                source = f.read()
                code = compile(source, path, "exec")
        else :
            with open(path, "rb") as f :
                # Check the magic number
                magic = f.read(4)
                if magic != imp.get_magic() :
                    return None
                # Skip timestamp
                f.read(4)
                # Load the code
                code = marshal.load(f)

        return code, path
    
    def _is_package(self, fullname):
        """ Test whether the plugin ``fullname`` is a package.
        """
        
        path = os.path.join(self._get_search_path(fullname), 
                            *fullname.split(".")[1:])
        return os.path.isdir(path)

# Only include importer if we are not in a frozen app (e.g. py2exe), otherwise
# this messes up the regular modules import
if not hasattr(sys, "frozen") :
    sys.meta_path.append(Importer())
