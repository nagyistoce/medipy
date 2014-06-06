/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _2b6e0502_6350_4124_9e72_24d739a638b5
#define _2b6e0502_6350_4124_9e72_24d739a638b5

#include "scu.h"

#include <Python.h>

class Get
: public SCU
{
public:
    Get(PyObject* connection,
        PyObject* root_lvl,
        PyObject* query_lvl,
        PyObject* dataset) throw(SCUException);
    virtual ~Get() {}

    //python __call__ equivalent
    virtual PyObject* operator()();

private:
    PyObject* _root_level;
    PyObject* _query_level;
    PyObject* _query_dataset;
};

#endif // _2b6e0502_6350_4124_9e72_24d739a638b5
