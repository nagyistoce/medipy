##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import logging
import re

def get_level(line, indent) :
    """ Return the level of indentation of a given line
    """
    
    level = 0
    while line.startswith(indent) :
        line = line[len(indent):]
        level += 1
    if re.match(r"^(\s+)", line) :
        raise Exception("Syntax Error")
    return level

def build_hierarchy(lines, root, top_level, first_children, first_siblings, path=[]) :
    """ Build a hierarchy of lines rooted at given line. Each item is of the
        form (name, display_name, sub_items)
    """
    
    sub_items = []
    
    line = lines[root].strip()
    leaf, display_name = line.split(" ", 1) if " " in line else (line, None)
    
    if root in first_children :
        child = first_children[root]
        sub_item = build_hierarchy(lines, child, top_level, 
                                   first_children, first_siblings,
                                   path = path+[leaf])
        sub_items.append(sub_item)
        
        while child in first_siblings :
            child = first_siblings[child]
            sub_item = build_hierarchy(lines, child, top_level, 
                                       first_children, first_siblings,
                                       path = path+[leaf])
            sub_items.append(sub_item)
    
    if not display_name :
        display_name = display_name or leaf.split(".")[-1]
        display_name = display_name.split("_")
        display_name[0] = display_name[0].capitalize()
        display_name = " ".join(display_name)
    
    if not sub_items :
        module_name = (".".join(["medipy"]+path)
                       if not leaf.startswith("medipy.") 
                       else ".".join(leaf.split(".")[:-1]))
        function_name = (leaf if not leaf.startswith("medipy.") 
                         else leaf.split(".")[-1])
        
        namespace = {}
        try :
            exec "from %s import %s"%(module_name, function_name) in namespace
        except ImportError :
            logging.error(
                "No function named %s in module %s", function_name, module_name)
            function = None
        else :
            function = namespace[function_name]
    else : 
        function = None
    
    return (display_name, function, sub_items)

def build_menu(filename):
    """ Return a hierarchical representation of a menu described in the file.
        Such a file will look like the following :  
        arithmetic
            addition
            soustraction
            division
        medimax3
            recalage Registration
                recalage_lineaire_multi_IM
                recalage_Bspline_topo
                ApplyTransfo3d_GUI Apply 3D Transform
        segmentation
            bet BET
    
        Each line has the following syntax : <indentation><name> [displayed name]
    """
    
    # Read all lines
    try :
        file = open(filename)
    except IOError, e :
        logging.error("Cannot open file %s : %s", filename, e)
        return []
    try :
        lines = file.readlines()
        # Remove comments
        lines = [x.split("#")[0] for x in lines]
        # Remove trailing spaces and end of lines
        lines = [x.rstrip() for x in lines]
        # Remove blank lines
        lines = [x for x in lines if x]
    except IOError :
        file.close()
        return []
    file.close()
    
    # Determine indentation used in the file : use the indentation of the first
    # line
    indent = ""
    for line in lines :
        match = re.match(r"^(\s+)", line)
        if match :
            indent = match.group(0)
            break
    
    # Build tree hierarchy : top-level items, first child and first sibling of
    # each item
    top_level = []
    first_children = {}
    first_siblings = {}
    
    for index in range(len(lines)-1) :
        
        line = lines[index]
        
        level = get_level(line, indent)
        
        if level == 0 :
            top_level.append(index)
        
        next_level = get_level(lines[index+1], indent)
        if next_level == level+1 :
            first_children[index] = index+1 
        
        # Look for siblings
        for other_line_index in range(index+1, len(lines)) :
            other_line = lines[other_line_index]
            
            other_level = get_level(other_line, indent)
            if level == other_level :
                # Found first sibling
                first_siblings[index] = other_line_index
                break
            elif level > other_level :
                # Went too far, in another menu
                break
    
    return [build_hierarchy(lines, x, top_level, first_children, first_siblings)
            for x in top_level]

if __name__ == "__main__" :
    import wx
    from medipy.base import find_resource
    from medipy.gui.menu_builder import fill_treectrl
    
    menu = build_menu(find_resource("gui/full.menu"))
    
    app = wx.PySimpleApp()
    frame = wx.Frame(None)
    sizer = wx.BoxSizer()
    frame.SetSizer(sizer)
    
    treectrl = wx.TreeCtrl(frame, 
        style=wx.TR_DEFAULT_STYLE|wx.TR_HAS_BUTTONS|wx.TR_SINGLE|wx.TR_HIDE_ROOT)
    fill_treectrl(menu, treectrl)
    sizer.Add(treectrl, 1, wx.EXPAND)
    
    frame.Show()
    app.MainLoop()