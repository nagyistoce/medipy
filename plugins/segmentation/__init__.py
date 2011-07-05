from api import *

def load_wrapitk_module(path, name):
    """ Load a WrapITK module, and add its symbols to the itk module """
    
    import sys
    
    import itkBase
    import itk
    
    sys.path.append(path)
    
    data = {}
    
    import BETImageFilterConfig
    data = BETImageFilterConfig.__dict__
    
    itkBase.module_data[name] = data

    namespace = {}
    itkBase.LoadModule(name, namespace)

    for k, v in namespace.items():
        if not hasattr(itk, k) :
            setattr(itk, k, v)

import os.path
load_wrapitk_module(os.path.dirname(__file__), "BETImageFilter")
