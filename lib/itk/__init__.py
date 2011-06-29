##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

from numpy_bridge import (array_to_itk_image, array_to_itk_vector_image,
                          itk_image_to_array, itk_vector_image_to_array, 
                          medipy_image_to_itk_image, itk_image_to_medipy_image, 
                          itk_to_dtype, dtype_to_itk)
                         

def load_wrapitk_module(path, name):
    """ Load a WrapITK module, and add its symbols to the itk module """
    import sys
    
    import itkBase
    import itk
    
    sys.path.append(path)
    
    data = {}
    
    import MediPyBridgeConfig
    data = MediPyBridgeConfig.__dict__
    
    itkBase.module_data[name] = data

    namespace = {}
    itkBase.LoadModule(name, namespace)

    for k, v in namespace.items():
        if not hasattr(itk, k) :
            setattr(itk, k, v)

import os.path
load_wrapitk_module(os.path.dirname(__file__), "MediPyBridge")
