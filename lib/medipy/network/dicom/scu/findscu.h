/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _65e3324c_d558_4860_a1cd_4892155c1c25
#define _65e3324c_d558_4860_a1cd_4892155c1c25

#include "scu.h"

#include <Python.h>

#include <dcmtk/dcmdata/dcdatset.h>

class Find
:public SCU
{
public :

    Find(   PyObject* connection,
            PyObject* root_lvl,
            PyObject* query_lvl,
            PyObject* dataset) throw(SCUException);
    virtual ~Find() {}
    
    //Python __call__ equivalent
    virtual PyObject* operator()();

private :
    //Find specific attributes
    PyObject* _root_level;
    PyObject* _query_level;
    PyObject* _query_dataset;
};

#endif // _65e3324c_d558_4860_a1cd_4892155c1c25
