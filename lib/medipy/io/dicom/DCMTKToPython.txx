#ifndef _2aae91c6_2f12_4888_9d6f_7d4c61fd0bff
#define _2aae91c6_2f12_4888_9d6f_7d4c61fd0bff

#include "DCMTKToPython.h"

#include <Python.h>

#include <dcmtk/config/osconfig.h>
#include <dcmtk/dcmdata/dctk.h>

template<typename TValue>
PyObject *
DCMTKToPython
::_to_python_number(TValue const & value, DcmEVR const & valid_vr) const
{
    if(valid_vr == EVR_FD || valid_vr == EVR_FL)
    {
        return PyFloat_FromDouble(value);
    }
    else
    {
        return PyInt_FromLong(value);
    }
}

template<typename TValue>
PyObject * 
DCMTKToPython
::_to_python_number(DcmElement * element, 
                    OFCondition (DcmElement::*getter)(TValue &, unsigned long)) const
{
    PyObject * python_value = NULL;
    
    DcmEVR const vr(DcmVR(element->getVR()).getValidEVR());
    
    unsigned long count = element->getVM();
    if(count > 1)
    {
        python_value = PyList_New(count);
        for(unsigned long i=0; i<count; ++i)
        {
            TValue value;
            (element->*getter)(value, i);
            
            PyList_SetItem(python_value, i, this->_to_python_number<TValue>(value, vr));
        }
    }
    else
    {
        TValue value;
        (element->*getter)(value, 0);
        
        python_value = this->_to_python_number<TValue>(value, vr);
    }
    
    return python_value;
}

#endif // _2aae91c6_2f12_4888_9d6f_7d4c61fd0bff
