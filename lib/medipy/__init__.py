import imp
import new
import os
import sys

class Importer(object):
    
    def __init__(self):
        self._plugins_path = ""
        if "MEDIPY_PLUGINS_PATH" in os.environ :
            self._plugins_path = os.environ["MEDIPY_PLUGINS_PATH"].split(os.pathsep)
    
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
        
sys.meta_path.append(Importer())
