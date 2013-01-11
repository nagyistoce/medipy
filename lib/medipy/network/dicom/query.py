##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import itertools

import medipy.io.dicom
import medipy.network.dicom

def relational(connection, root, level, query) :
    """ Perform a relation DICOM query, as described in PS 3.4-2011, C.4.1.2.2.1
        
        The resulting datasets will contain the highest level below the given
        level where all arguments from the query are matched. The keys for all
        those levels will be present.
    """
    
    # Add the current level key if needed
    keys = medipy.network.dicom.scu.Find.keys[root][level]
    for key in keys :
        if key not in query :
            query[key] = ""
    
    # Perform the query at the current level
    find = medipy.network.dicom.scu.Find(connection, root, level, query)
    result = find()
    
    # Check what attributes have been matched
    matched = {}
    all_matched = True
    for attribute in query :
        everywhere = True
        for dataset in result :
            value = dataset.get(attribute, medipy.io.dicom.UN(None))
            if isinstance(value, medipy.io.dicom.SQ) :
                query_item = query[attribute].value[0]
                for item in value :
                    for query_item_attribute in query_item :
                        if not item.get(query_item_attribute, None).value :
                            everywhere = False
                            break
            elif not value.value :
                everywhere = False
                break
        matched[attribute] = everywhere
        if not everywhere :
            all_matched = False
    
    # If everything has been matched, we are done
    if all_matched :
        return result
    
    # Otherwise, if we are already at the image level, we have no match
    if level == "image" :
        return []
    
    # Build a query at the sub-level
    sub_query = medipy.io.dicom.DataSet()
    # Set the key for the current level using the current level result
    keys = medipy.network.dicom.scu.Find.keys[root][level]
    for key in keys :
        value = "\\".join(x[key].value for x in result)
        sub_query[key] = value
    # Only keep the non-matched, non-key attributes
    for attribute, value in query.items() : 
        if not matched[attribute] :
            sub_query[attribute] = value
        else :
            keys = medipy.network.dicom.scu.Find.keys[root].values()
            if attribute in itertools.chain(*keys) and attribute not in sub_query:
                sub_query[attribute] = "\\".join(x[attribute].value for x in result)
    # Add the key for the lower-level if not already specified
    sub_level = {
        "patient" : "study",
        "study" : "series",
        "series" : "image"
    }[level]
    keys = medipy.network.dicom.scu.Find.keys[root][sub_level]
    for key in keys :
        if key not in sub_query :
            sub_query[key] = ""
    
    sub_result = relational(connection, root, sub_level, sub_query)
    
    # Index the result of the current level with the key of the current level
    result_index = {}
    keys = medipy.network.dicom.scu.Find.keys[root][level]
    for dataset in result :
        result_index[tuple([dataset[key].value for key in keys])] = dataset
    # Update the sub_result with the matched attributes of the current level
    for sub_dataset in sub_result :
        dataset = result_index[tuple([sub_dataset[key].value for key in keys])]
        for attribute, attribute_matched in matched.items() :
            if attribute_matched :
                sub_dataset[attribute] = dataset[attribute]
    
    return sub_result
