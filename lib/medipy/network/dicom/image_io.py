##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import copy
import multiprocessing

import medipy.base
import medipy.io.dicom
import medipy.network.dicom.scu
import medipy.network.dicom.query

def load(connection, query, retrieve_method, retrieve_data=None) :
    """ Load an image from a DICOM node using the given query. The query must
        return a single series. If the series contains several stacks, a list of 
        images is returned; otherwise an image is returned.
        
        The value of ``retrieve_method`` and ``retrieve_data``  can be one of 
        the following:

        * ``"GET"``: the datasets corresponding to the query results are
          retrieved using C-GET, no extra argument will be provided.
        * ``("MOVE", destination)``: the datasets corresponding to the query 
          results are retrieved using C-MOVE. The extra argument will contain 
          the destination AE title.
        * ``("WADO", url)``: the datasets corresponding to the query results are 
          retrieved using WADO, using the given URL.

        Usage example using WADO retrieval: ::

            query = medipy.io.dicom.DataSet()
            query.patients_name = "Doe^John"
            query.series_description = "T1 3D 1mm"
            
            image = medipy.network.dicom.image_io.load(connection, query, "WADO", "http://www.example.com/wado")
    """

    functions = {
        "GET": _get,
        "MOVE": _move,
        "WADO": _wado
    }
    
    try:
        function = functions[retrieve_method]
    except KeyError, e:
        raise medipy.base.Exception(
            "Unknown retrieve method: {0}".format(retrieve_method))

    datasets = function(connection, query, retrieve_data)

    # Check that only a single series matches the query
    series_uids =  set([x.get("series_instance_uid", medipy.io.dicom.UI(None)).value for x in datasets])
    if len(series_uids) > 1 :
        raise medipy.base.Exception("Only one series must match query")

    stacks = medipy.io.dicom.stacks(datasets)
    
    if len(stacks) == 0 :
        raise medipy.base.Exception("No stack matching query")
    
    pool = multiprocessing.Pool()
    async_results = [pool.apply_async(medipy.io.dicom.image, (stack,))
                     for stack in stacks]

    pool.close()
    pool.join()

    if len(stacks) == 1 :
        return async_results[0].get()
    else :
        return [x.get() for x in async_results]

def _get(connection, query, dummy):
    query_results = medipy.network.dicom.query.relational(
        connection, "patient", "patient", query)
    datasets = []
    for item in query_results:
        get = medipy.network.dicom.scu.Get(connection, "patient", "series", item)
        results = get()
        datasets.extend(results)
    
    return datasets

def _move(connection, query, destination):
    query_results = medipy.network.dicom.query.relational(
        connection, "patient", "patient", query)
    datasets = []
    for item in query_results:
        move = medipy.network.dicom.scu.Move(
            connection, "patient", "series", destination, item)
        results = move()
        datasets.extend(results)
    
    return datasets

def _wado(connection, query, url):
    query = copy.deepcopy(query)
    query.setdefault("sop_instance_uid", None)
    query_results = medipy.network.dicom.query.relational(
        connection, "patient", "patient", query)
    
    pool = multiprocessing.Pool()
    async_results = [pool.apply_async(medipy.network.dicom.wado.get, (url, item)) 
                     for item in query_results]
    pool.close()
    pool.join()
    datasets = [x.get() for x in async_results]

    return datasets
