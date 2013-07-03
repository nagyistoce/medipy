##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

def fill_treectrl(menu, treectrl, match):
    """ Fill a treectrl with the menu. Each menu item is a triplet 
        (label, function, children). If no function is associated with this item
        it is set to None.
    """
    root = treectrl.AddRoot("Functions")
    add_tree_nodes(root, menu, treectrl,match)

def add_tree_nodes(parentItem, items, treectrl,match):
    for label,function,children in items:
        if match in label.lower() or search_in_menu(children,match) or match=="" :
            tree_item = treectrl.AppendItem(parentItem, label)
            treectrl.SetPyData(tree_item, function)
            add_tree_nodes(tree_item, children, treectrl,match)

def search_in_menu(items,match) :
    key = False
    for label,function,children in items:
        key = key or ( match in label.lower() or any( search_in_menu([child],match) for child in children ) )
    return key
