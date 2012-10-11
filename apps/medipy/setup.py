import os
import deployement

if __name__ == "__main__" :
    includes = ["itk"]
    plugins_root = os.environ["MEDIPY_PLUGINS_PATH"]
    medipy_plugins = [entry 
                      for entry in os.listdir(plugins_root)
                      if os.path.isdir(os.path.join(plugins_root, entry)) and
                         os.path.isfile(os.path.join(plugins_root, entry, "api.py")) and
                         entry != "medimax"
                     ]
    includes.extend(["medipy.{0}.api".format(entry) for entry in medipy_plugins])
    if "hsqc" in medipy_plugins :
        # TODO : move this to the plugin directory
        includes.extend(["openopt", "scipy.interpolate",])
    deployement.setup("MediPy", "medipy", includes, medipy_plugins)