##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import logging
import os
import sys
import traceback

def title_for(title):
    """ Create a title from a underscore-separated string.
    """
    
    return title.replace("_", " ").capitalize()

def build_menu(path, root) :
    """ Build a menu from given path : path is recursively scanned for files
        name api.py, which are parsed. Any function declared inside will be 
        added to the menu.
    """
    
    result = []
    
    names = os.listdir(path)
    
    for name in names :
        if os.path.isdir(os.path.join(path, name)) :
            sub_result = build_menu(os.path.join(path, name), root)
            if len(sub_result)>0 :
                label = title_for(name)
                result.append((label, None, sub_result))
    result.sort(lambda x,y : cmp(x[0], y[0]))
    
    if not any([x in names for x in ["api", "api.py", "api.pyc"]]) :
        # No api 
        return result
    
    module_name = ["medipy"] + path[len(os.path.normpath(root))+1:].split(os.path.sep) + ["api"]
    # Remove empty elements
    module_name = ".".join([x for x in module_name if x])
    
    try :
        __import__(module_name)
    except ImportError, e :
        # An api exists, but it cannot be imported
        logging.debug("Could not import {0}: {1}".format(module_name, e))
        return result
    
    module = sys.modules[module_name]
    
    items = []
    for name in dir(module) :
        if name.startswith("_") :
            continue
        
        label = title_for(name)
        items.append((label, getattr(module, name), []))
    items.sort(lambda x,y : cmp(x[0], y[0]))
    result.extend(items)
    
    return result
