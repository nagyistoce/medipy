##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

def fill_treectrl(menu, treectrl):
    """ Fill a treectrl with the menu. Each menu item is a triplet 
        (label, function, children). If no function is associated with this item
        it is set to None.
    """
    root = treectrl.AddRoot("Functions")
    add_tree_nodes(root, menu, treectrl)

def add_tree_nodes(parentItem, items, treectrl):
    for label, function, children in items:
        tree_item = treectrl.AppendItem(parentItem, label)
        treectrl.SetPyData(tree_item, function)
        add_tree_nodes(tree_item, children, treectrl)