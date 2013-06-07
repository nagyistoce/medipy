##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import os
import sys

import itk
import itkBase

def load_wrapitk_module(path, name):
    """ Load a WrapITK module, and add its symbols to the itk module """

    old_path = sys.path
    
    sys.path = sys.path+[path]
    
    config_module_name = "{0}Config".format(name)
    config_module = __import__(config_module_name)
    
    path = os.path.dirname(os.path.abspath(config_module.__file__))
    if path not in sys.path :
        sys.path = sys.path + [path]
    
    data = config_module.__dict__
    itkBase.module_data[name] = data

    namespace = {}
    itkBase.LoadModule(name, namespace)

    for k, v in namespace.items():
        if not hasattr(itk, k) :
            setattr(itk, k, v)
    
    sys.path = old_path
