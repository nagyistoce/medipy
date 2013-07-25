##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import multiprocessing

import numpy

import medipy.base
import medipy.io.dicom
import query as query_module

def load(connection, query, retrieve) :
    """ Load an image from a DICOM node using the given query. The query must
        return a single serie. If the serie contains several stacks, a list of 
        images is returned ; otherwise an image is returned.
        
        The value of ``retrieve`` can be one of the following :

        * ``"GET"`` : the datasets corresponding to the query results are
          retrieved using C-GET, no extra argument will be provided.
        * ``("MOVE", destination)`` : the datasets corresponding to the query 
          results are retrieved using C-MOVE. The extra argument will contain 
          the destination AE title.
        * ``("WADO", url)`` : the datasets corresponding to the query results are 
          retrieved using WADO, using the given URL.

        Usage example using WADO retrieval : ::

            query = medipy.io.dicom.DataSet()
            query.patients_name = "Doe^John"
            query.series_description = "T1 3D 1mm"
            
            retrieve = ("WADO", "http://www.example.com/wado")

            image = medipy.network.dicom.image_io.load(connection, query, retrieve)
    """

    # Make sure we retrieve the SOP Instance UID
    if "sop_instance_uid" not in query :
        new_query = medipy.io.dicom.DataSet()
        new_query.update(query)
        new_query.sop_instance_uid = None
        
        query = new_query

    query_results = query_module.relational(connection, "patient", "patient", query)
    if not query_results :
        raise medipy.base.Exception("No image matching query")

    # Check that only a single series matches the query
    series_uids =  set([x.get("series_instance_uid", None) 
                        for x in query_results])
    if len(series_uids) > 1 :
        raise medipy.base.Exception("Only one series must match query")

    if retrieve[0] == "GET" :
        raise NotImplementedError()
        
    elif retrieve[0] == "MOVE" :
        move_query = medipy.io.dicom.DataSet(sop_instance_uid='')
        for item in query_results:
            sop_uid = str(item.sop_instance_uid.value)
            mv_sop = str(move_query.sop_instance_uid.value) + '\\' + sop_uid
            move_query.__setattr__('sop_instance_uid',mv_sop)
            
        move = medipy.network.dicom.scu.Move(connection,"patient","image",
                retrieve[1],move_query)

        results = move()
        datasets = medipy.io.dicom.split.images(results)

    elif retrieve[0] == "WADO" :
        pool = multiprocessing.Pool()
        async_results = [pool.apply_async(medipy.network.dicom.wado.get, 
                                          (retrieve[1], wado_query)) 
                         for wado_query in query_results]
        pool.close()
        pool.join()
        datasets = [x.get() for x in async_results]
        
    else :
        raise medipy.base.Exception("Unknown retrieve mode {0!r}".format(retrieve))

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
