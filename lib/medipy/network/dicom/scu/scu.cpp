/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "scu.h"
#include "scuexcept.h"
 
#include <sstream>
#include <stdexcept>
#include <typeinfo>

#include <Python.h>
#include "DcmSCU.h"

//Constructor
SCU
::SCU()
: _connection(Py_None)
{
    this->GetPythonClasses();
}

SCU
::SCU(PyObject* connection) throw(SCUException)
: _connection(Py_None)
{
    this->GetPythonClasses();
    this->SetConnection(connection);
    if(PyErr_Occurred())
    {
        PyObject* type; PyObject* message; PyObject* trace;
        PyErr_Fetch(&type, &message, &trace);
        char* string_err = PyString_AsString(message);
        throw SCUException(string_err);
    }
}

//Destructor
SCU
::~SCU()
{
    Py_DECREF(this->_medipy_network_dicom_connection);
    Py_DECREF(this->_medipy_base_exception);
    this->CloseConnection();
}


PyObject* 
SCU
::GetConnection() const
{
    return this->_connection;
}

std::string
SCU
::GetAttr(PyObject* connection, const char* attribute) const  throw(SCUException)
{
    std::string value;
    PyObject* data = NULL;
    
    PyObject* object = PyObject_GetAttrString(connection,attribute);
    if(object == NULL || !(PyString_Check(object)||PyUnicode_Check(object)))
    {
        std::stringstream message;
        message << "Cannot access " << attribute << " as string.";
        throw SCUException(message.str().c_str());
    }
    
    if(PyUnicode_Check(object))
    {
        data = PyUnicode_AsEncodedString(object, "utf-8", "strict");
    }
    else
    {
        data = object;
    }
    
    value = std::string(PyString_AsString(data));
    
    if(PyUnicode_Check(object))
    {
        Py_DECREF(data);
    }
    return value;
}

PyObject*
SCU
::SetConnection(PyObject* connection)
{
    //Reset Params and Drop Network
    if(this->_connection != Py_None)
    {
        try
        {
            this->CloseConnection();
        }
        catch(std::exception const & e)
        {
            PyErr_SetString(this->_medipy_base_exception, e.what());
            return NULL;
        }
    }
    
    if(!PyObject_IsInstance(connection, this->_medipy_network_dicom_connection))
    {
        PyErr_SetString(this->_medipy_base_exception, 
            "Connection is not an instance of medipy.network.dicom.ConnectionBase");
        return NULL;
    }
    std::string const host = this->GetAttr(connection,"host");
    std::string const called = this->GetAttr(connection,"called_ae_title");
    std::string const calling = this->GetAttr(connection,"calling_ae_title");
    
    long port;
    {
        PyObject* object = PyObject_GetAttrString(connection,"port");
        if(object == NULL || !PyInt_Check(object))
        {
            PyErr_SetString(this->_medipy_base_exception, 
                "Cannot access port as long");
            return NULL;
        }
        port = PyInt_AsLong(object);
    }
    
    //Set DcmSCU association parameters
    this->_scu.setAETitle(calling.c_str());
    this->_scu.setPeerHostName(host.c_str());
    this->_scu.setPeerPort(port);
    this->_scu.setPeerAETitle(called.c_str());
    
    this->_connection = connection;
    Py_INCREF(this->_connection);
    
    return Py_None;
}

void
SCU
::CloseConnection()
{

    Py_DECREF(this->_connection);
    this->_connection = NULL;
}

void
SCU
::GetPythonClasses()
{
    {
        PyObject* module = PyImport_ImportModule("medipy.network.dicom"); //New reference
        PyObject* moduleDict = PyModule_GetDict(module);
        this->_medipy_network_dicom_connection = PyDict_GetItemString(moduleDict, "ConnectionBase");
        Py_INCREF(this->_medipy_network_dicom_connection);
        Py_DECREF(module);
    }
    
    {
        PyObject* module = PyImport_ImportModule("medipy.base"); //New reference
        PyObject* moduleDict = PyModule_GetDict(module);
        this->_medipy_base_exception = PyDict_GetItemString(moduleDict, "Exception");
        Py_INCREF(this->_medipy_base_exception);
        Py_DECREF(module);
    }
}

PyObject*
SCU
::operator()()
{
    OFCondition cond = this->_scu.initNetwork();
    if (cond.bad())
    {
        PyErr_SetString(this->_medipy_base_exception,cond.text());
        return NULL;
    }

    cond = this->_scu.negotiateAssociation();
    if (cond.bad())
    {
        PyErr_SetString(this->_medipy_base_exception,cond.text());
        return NULL;
    }
    
    return Py_None;
}
