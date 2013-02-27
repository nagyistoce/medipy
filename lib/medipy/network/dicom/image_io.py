##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import multiprocessing.dummy

import medipy.base
import medipy.io.dicom
import medipy.io.dicom.split
import query as query_module

def load(query, connection, retrieve, *args) :
    """ Load an image from a DICOM node using the given query.
        The value of ``retrieve`` can be one of the following :

        * ``"GET"`` : the datasets corresponding to the query results are
          retrieved using C-GET, no extra argument will be provided.
        * ``"MOVE"`` : the datasets corresponding to the query results are
          retrieved using C-MOVE. The extra argument will contain the 
          destination AE title.
        * an URL : the datasets corresponding to the query results are 
          retrieved using WADO. The extra argument will contain the WADO URL.

        Usage example using WADO retrieval : ::

            query = medipy.io.dicom.DataSet()
            query.patients_name = "Doe^John"
            query.series_description = "T1 3D 1mm"

            image = load_image(query, connection, "WADO", "http://www.example.com/wado")
    """

    # Make sure we retrieve the SOP Instance UID
    if "sop_instance_uid" not in query :
        query.sop_instance_uid = None

    query_results = query_module.relational(connection, "patient", "patient", query)

    if not query_results :
        raise medipy.base.Exception("No image matching query")

    # Check that only a single series matches the query
    series_uids =  set([x.get("series_instance_uid", None) 
                        for x in query_results])
    if len(series_uids) > 1 :
        raise medipy.base.Exception("Only one series must match query")

    if retrieve == "GET" :
        raise NotImplementedError()
    elif retrieve == "MOVE" :
        raise NotImplementedError()
    elif retrieve == "WADO" :
        wado_url = args[0]
        def worker(dataset) :
            return medipy.network.dicom.wado.get(wado_url, dataset)
        pool = multiprocessing.dummy.Pool(2*multiprocessing.cpu_count())
        async_results = [pool.apply_async(worker, (result,)) 
                         for result in query_results]
        pool.close()
        pool.join()
        datasets = [x.get() for x in async_results]
    else :
        raise medipy.base.Exception("Unknown retrieve mode {0!r}".format(retrieve))

    stacks = medipy.io.dicom.stacks(datasets)
    if len(stacks) == 0 :
        raise medipy.base.Exception("No stack matching query")
    elif len(stacks) > 1 :
        raise medipy.base.Exception("Only one stack must match query")

    return medipy.io.dicom.image(stacks[0])

