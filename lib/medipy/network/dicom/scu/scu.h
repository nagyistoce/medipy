/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _1139136d_802f_41ae_9d88_12e44050f186
#define _1139136d_802f_41ae_9d88_12e44050f186

#include <Python.h>
#include <string>

#include <dcmtk/config/osconfig.h>
#include "MedipyDcmSCU.h"
#include "scuexcept.h"

class SCU
{
public :
    //Constructor
    SCU();
    SCU(PyObject* connection)throw(SCUException);
    
    //Destructor
    virtual ~SCU();
    
    PyObject* GetConnection() const;
    PyObject* SetConnection(PyObject* connection);
    
    //Python __call__ equivalent
    virtual PyObject* operator()();
    
protected :
    //medipy.network.dicom.Connection object
    PyObject* _connection;
    
    PyObject* _medipy_network_dicom_connection;
    PyObject* _medipy_base_exception;
    
    //Dcmtk specific class : handle association
    MedipyDcmSCU _scu;
    
    const char* _abstractSyntax;
    
    void CloseConnection();
    
    void GetPythonClasses();
    
    std::string GetAttr(PyObject* connection,
                        const char* attribute) const throw(SCUException);

};

#endif // _1139136d_802f_41ae_9d88_12e44050f186
