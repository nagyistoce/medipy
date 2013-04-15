#ifndef _2aae91c6_2f12_4888_9d6f_7d4c61fd0bff
#define _2aae91c6_2f12_4888_9d6f_7d4c61fd0bff

#include "DCMTKToPython.h"

#include <Python.h>

#include <dcmtk/config/osconfig.h>
#include <dcmtk/dcmdata/dctk.h>

template<typename TValue>
PyObject * 
DCMTKToPython
::_to_python_number(DcmElement * element, 
                    OFCondition (DcmElement::*getter)(TValue &, unsigned long)) const
{
    PyObject * python_value = NULL;
    
    DcmEVR const vr(element->getVR());
    
    unsigned long count = element->getVM();
    if(count > 1)
    {
        python_value = PyList_New(count);
        for(unsigned long i=0; i<count; ++i)
        {
            TValue value;
            (element->*getter)(value, i);
            
            if(vr == EVR_FD || vr == EVR_FL)
            {
                PyList_SetItem(python_value, i, PyFloat_FromDouble(value));
            }
            else
            {
                PyList_SetItem(python_value, i, PyInt_FromLong(value));
            }
        }
    }
    else
    {
        TValue value;
        (element->*getter)(value, 0);
        
        if(vr == EVR_FD || vr == EVR_FL)
        {
            python_value = PyFloat_FromDouble(value);
        }
        else
        {
            python_value = PyInt_FromLong(value);
        }
    }
    
    return python_value;
}

#endif // _2aae91c6_2f12_4888_9d6f_7d4c61fd0bff
