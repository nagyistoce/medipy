##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import imp
import sys

import itk
import itkBase

def load_wrapitk_module(path, name):
    """ Load a WrapITK module, and add its symbols to the itk module """
    
    sys.path.append(path)
    
    config_module_name = "{0}Config".format(name)
    file, pathname, description = imp.find_module(config_module_name, [path])
    config_module = imp.load_module(config_module_name, file, pathname, description)
    
    data = config_module.__dict__
    itkBase.module_data[name] = data

    namespace = {}
    itkBase.LoadModule(name, namespace)

    for k, v in namespace.items():
        if not hasattr(itk, k) :
            setattr(itk, k, v)