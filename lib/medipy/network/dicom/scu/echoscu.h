/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _9472f354_77b6_436d_ba13_5bbb36b8b8a5
#define _9472f354_77b6_436d_ba13_5bbb36b8b8a5

#include <Python.h>
#include "scu.h"

class Echo
:public SCU
{
public :
    Echo(PyObject* connection);
    virtual ~Echo() {}
    
    //Python __call__ equivalent
    virtual PyObject* operator()();
};

#endif //_9472f354_77b6_436d_ba13_5bbb36b8b8a5
