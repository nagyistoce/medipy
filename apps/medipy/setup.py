import os
import deployement

if __name__ == "__main__" :
    includes = ["itk"]
    plugins_root = os.environ["MEDIPY_PLUGINS_PATH"]
    medipy_plugins = [entry 
                      for entry in os.listdir(plugins_root)
                      if os.path.isdir(os.path.join(plugins_root, entry)) and
                         entry != "medimax"
                     ]
    if "hsqc" in medipy_plugins :
        # TODO : move this to the plugin directory
        includes.extend([
            "decimal", "openopt", 
            "scipy.interpolate", "scipy.signal", 
            "sklearn", "sklearn.cluster", "sklearn.metrics", "sklearn.svm"
        ])
    deployement.setup("MediPy", "medipy", includes, medipy_plugins)