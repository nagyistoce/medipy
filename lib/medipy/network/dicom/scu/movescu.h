/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _7f96cbf1_8e15_4b55_9815_9ecb3bedc61f
#define _7f96cbf1_8e15_4b55_9815_9ecb3bedc61f

#include "scu.h"

#include <Python.h>

#include <dcmtk/dcmdata/dcdatset.h>
#include <dcmtk/dcmnet/dimse.h>
#include <dcmtk/dcmdata/dcfilefo.h>

//Dcmtk need it
struct QuerySyntax {
    const char *findSyntax;
    const char *moveSyntax;
};

//Our own Move definition (useful stuff down here)
class Move
:public SCU
{
public :
    Move( PyObject* connection,
          PyObject* root_lvl,
          PyObject* query_lvl,
          PyObject* destination = Py_None,
          PyObject* dataset = Py_None);
    virtual ~Move() {}
    
    //Python __call__ equivalent
    virtual PyObject* operator()();

private :
    //Move specific attributes
    PyObject* _root_level;
    PyObject* _query_level;
    PyObject* _query_dataset;
    PyObject* _destination;
    
    //TCP Option
    QuerySyntax _querySyntax;
};

#endif // _7f96cbf1_8e15_4b55_9815_9ecb3bedc61f
