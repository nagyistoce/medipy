##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import copy
import logging
import os
import sys

def title_for(title):
    """ Create a title from a underscore-separated string.
    """
    
    return title.replace("_", " ").capitalize()

def build_menu(path) :
    """ Build a menu from given path : path is recursively scanned for files
        name api.py, which are parsed. Any function declared inside will be 
        added to the menu.
    """
    
    result = []
    
    names = os.listdir(path)
    
    for name in names :
        if os.path.isdir(os.path.join(path, name)) :
            sub_result = build_menu(os.path.join(path, name))
            if len(sub_result)>0 :
                label = title_for(name)
                result.append((label, None, sub_result))
    result.sort(lambda x,y : cmp(x[0], y[0]))
            
    if "api.py" in names :
        api_file = os.path.join(path, "api.py")
        api_globals = {}
        api_locals = {}
        
        # Backup sys.path to avoid pollution
        old_path = copy.copy(sys.path)
        sys.path.insert(0, path)
        try :
            execfile(api_file, api_globals, api_locals)
        except Exception, e :
            logging.warn("Could not load %s/api.py : %s"%(path, e))
        else :
            functions = []
            for function_name in api_locals :
                label = title_for(function_name)
                functions.append((label, api_locals[function_name], []))
            functions.sort(lambda x,y : cmp(x[0], y[0]))
            result.extend(functions)
        sys.path = old_path
    
    return result
