##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

from tag import Tag
from dataset import DataSet
from parse import parse
from split import series, stacks
from reconstruction import image
from misc import load_dicomdir_records, uid_and_description

def load_wrapitk_module(path, name):
    """ Load a WrapITK module, and add its symbols to the itk module """
    
    import sys
    
    import itkBase
    import itk
    
    sys.path.append(path)
    
    data = {}
    
    import AssembleTilesImageFilterConfig
    data = AssembleTilesImageFilterConfig.__dict__
    
    itkBase.module_data[name] = data

    namespace = {}
    itkBase.LoadModule(name, namespace)

    for k, v in namespace.items():
        if not hasattr(itk, k) :
            setattr(itk, k, v)

import os.path
load_wrapitk_module(os.path.dirname(__file__), "AssembleTilesImageFilter")
