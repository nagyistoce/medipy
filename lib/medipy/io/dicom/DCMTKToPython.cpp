/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "DCMTKToPython.h"

#include <Python.h>

#include <stdexcept>
#include <string>

#include <dcmtk/config/osconfig.h>
#include <dcmtk/dcmdata/dctk.h>

PyObject *
DCMTKToPython
::operator()(DcmObject * dataset)
{
    PyObject* args = PyTuple_New(0);
    PyObject* python_dataset = PyObject_CallObject(this->_medipy_io_dicom_DataSet, args);
    Py_DECREF(args);
    
    DcmObject * it = NULL;
    while(NULL != (it = dataset->nextInContainer(it)))
    {
        if(it->getTag() == DCM_SpecificCharacterSet)
        {
            // Specific Character Set: setup internal iconv converter
        }
        
        if(it->getETag() == 0)
        {
            // Group length, do nothing
            continue;
        }
        else
        {
            this->_add_element(it, python_dataset);
        }
    }
    
    return python_dataset;
}

/*******************************************************************************
 * Specializations of DCMTKToPython::_to_python for the different VRs.
 ******************************************************************************/
 
template<>
PyObject *
DCMTKToPython
::_to_python<EVR_AE>(DcmObject * element) const
{
    return this->_to_python_text(dynamic_cast<DcmByteString*>(element), false);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_AS>(DcmObject * element) const
{
    return this->_to_python_text(dynamic_cast<DcmByteString*>(element), false);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_AT>(DcmObject * element) const
{
    return this->_to_python_number(dynamic_cast<DcmElement*>(element),
                                   &DcmElement::getUint32);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_CS>(DcmObject * element) const
{
    return this->_to_python_text(dynamic_cast<DcmByteString*>(element), false);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_DA>(DcmObject * element) const
{
    return this->_to_python_text(dynamic_cast<DcmByteString*>(element), false);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_DT>(DcmObject * element) const
{
    return this->_to_python_text(dynamic_cast<DcmByteString*>(element), false);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_DS>(DcmObject * element) const
{
    PyObject * python_value = NULL;
    
    DcmDecimalString * ds = dynamic_cast<DcmDecimalString*>(element);
    unsigned long count = ds->getVM();

    if(count > 1)
    {
        python_value = PyList_New(count);
        
        for(unsigned long i=0; i<count; ++i)
        {
            Float64 value;
            ds->getFloat64(value, i);
            PyList_SetItem(python_value, i, PyFloat_FromDouble(value));
        }
    }
    else
    {
        Float64 value;
        ds->getFloat64(value, 0);
        python_value = PyFloat_FromDouble(value);
    }
    
    return python_value;
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_FD>(DcmObject * element) const
{
    return this->_to_python_number(dynamic_cast<DcmElement*>(element),
                                   &DcmElement::getFloat64);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_FL>(DcmObject * element) const
{
    return this->_to_python_number(dynamic_cast<DcmElement*>(element),
                                   &DcmElement::getFloat32);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_IS>(DcmObject * element) const
{
    PyObject * python_value = NULL;
    
    DcmIntegerString * is = dynamic_cast<DcmIntegerString*>(element);
    unsigned long count = is->getVM();

    if(count > 1)
    {
        python_value = PyList_New(count);
        
        for(unsigned long i=0; i<count; ++i)
        {
            Sint32 value;
            is->getSint32(value, i);
            PyList_SetItem(python_value, i, PyInt_FromLong(value));
        }
    }
    else
    {
        Sint32 value;
        is->getSint32(value, 0);
        python_value = PyInt_FromLong(value);
    }
    
    return python_value;
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_LO>(DcmObject * element) const
{
    return this->_to_python_text(dynamic_cast<DcmByteString*>(element), true);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_LT>(DcmObject * element) const
{
    return this->_to_python_text(dynamic_cast<DcmByteString*>(element), true);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_OB>(DcmObject * element) const
{
    return this->_to_python_binary(dynamic_cast<DcmElement*>(element));
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_OF>(DcmObject * element) const
{
    return this->_to_python_binary(dynamic_cast<DcmElement*>(element));
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_OW>(DcmObject * element) const
{
    return this->_to_python_binary(dynamic_cast<DcmElement*>(element));
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_PN>(DcmObject * element) const
{
    return this->_to_python_text(dynamic_cast<DcmByteString*>(element), true);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_SH>(DcmObject * element) const
{
    return this->_to_python_text(dynamic_cast<DcmByteString*>(element), true);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_SL>(DcmObject * element) const
{
    return this->_to_python_number(dynamic_cast<DcmElement*>(element),
                                   &DcmElement::getSint32);
}

// SQ is not processed here

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_SS>(DcmObject * element) const
{
    return this->_to_python_number(dynamic_cast<DcmElement*>(element),
                                   &DcmElement::getSint16);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_ST>(DcmObject * element) const
{
    return this->_to_python_text(dynamic_cast<DcmByteString*>(element), true);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_TM>(DcmObject * element) const
{
    return this->_to_python_text(dynamic_cast<DcmByteString*>(element), false);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_UI>(DcmObject * element) const
{
    return this->_to_python_text(dynamic_cast<DcmByteString*>(element), false);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_UL>(DcmObject * element) const
{
    return this->_to_python_number(dynamic_cast<DcmElement*>(element),
                                   &DcmElement::getUint32);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_UN>(DcmObject * element) const
{
    return this->_to_python_binary(dynamic_cast<DcmElement*>(element));
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_US>(DcmObject * element) const
{
    return this->_to_python_number(dynamic_cast<DcmElement*>(element),
                                   &DcmElement::getUint16);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_UT>(DcmObject * element) const
{
    return this->_to_python_text(dynamic_cast<DcmByteString*>(element), true);
}

/*******************************************************************************
 * End of specializations of DCMTKToPython::_to_python for the different VRs.
 ******************************************************************************/

PyObject *
DCMTKToPython
::_to_python_text(DcmByteString * element, bool use_utf8) const
{
    PyObject * python_value = NULL;
    
    return python_value;
}

PyObject *
DCMTKToPython
::_to_python_binary(DcmElement * element) const
{
    PyObject * python_value = NULL;
    
    DcmOtherByteOtherWord* byte_string = dynamic_cast<DcmOtherByteOtherWord*>(element);
    if(byte_string != NULL)
    {
        char const * begin(NULL);
        if(element->getVR() != EVR_OW)
        {
            Uint8 * data;
            byte_string->getUint8Array(data);
            begin = reinterpret_cast<char const *>(data);
        }
        else
        {
            Uint16 * data;
            byte_string->getUint16Array(data);
            begin = reinterpret_cast<char const *>(data);
        }
        
        python_value = PyString_FromStringAndSize(begin, element->getLength());
    }
    else
    {
        throw std::runtime_error(
            std::string("Cannot handle conversion of ")+
                DcmVR(element->getVR()).getValidVRName());
    }
    
    return python_value;
}

void
DCMTKToPython
::_add_element(DcmObject * element, PyObject * python_dataset) const
{
    DcmEVR const dcmtk_vr(element->getVR());
    char const * const vr_name = DcmVR(dcmtk_vr).getValidVRName();

    PyObject * python_value = NULL;
    
    if(dcmtk_vr == EVR_SQ)
    {
        python_value = PyList_New(0);
        
        DcmSequenceOfItems * sequence = dynamic_cast<DcmSequenceOfItems*>(element);
        DcmObject * sequence_it = NULL;
        while(NULL != (sequence_it = sequence->nextInContainer(sequence_it)))
        {
            DCMTKToPython converter;
            PyObject * item_value = converter(sequence_it);
            PyList_Append(python_value, item_value);
        }
    }
    else if(element->getLength() == 0)
    {
        python_value = Py_None;
    }
    else
    {
        if(dcmtk_vr == EVR_AE) { python_value = this->_to_python<EVR_AE>(element); }
        if(dcmtk_vr == EVR_AS) { python_value = this->_to_python<EVR_AS>(element); }
        if(dcmtk_vr == EVR_AT) { python_value = this->_to_python<EVR_AT>(element); }
        if(dcmtk_vr == EVR_CS) { python_value = this->_to_python<EVR_CS>(element); }
        if(dcmtk_vr == EVR_DA) { python_value = this->_to_python<EVR_DA>(element); }
        if(dcmtk_vr == EVR_DT) { python_value = this->_to_python<EVR_DT>(element); }
        if(dcmtk_vr == EVR_DS) { python_value = this->_to_python<EVR_DS>(element); }
        if(dcmtk_vr == EVR_FD) { python_value = this->_to_python<EVR_FD>(element); }
        if(dcmtk_vr == EVR_FL) { python_value = this->_to_python<EVR_FL>(element); }
        if(dcmtk_vr == EVR_IS) { python_value = this->_to_python<EVR_IS>(element); }
        if(dcmtk_vr == EVR_LO) { python_value = this->_to_python<EVR_LO>(element); }
        if(dcmtk_vr == EVR_LT) { python_value = this->_to_python<EVR_LT>(element); }
        if(dcmtk_vr == EVR_OB) { python_value = this->_to_python<EVR_OB>(element); }
        if(dcmtk_vr == EVR_OF) { python_value = this->_to_python<EVR_OF>(element); }
        if(dcmtk_vr == EVR_OW) { python_value = this->_to_python<EVR_OW>(element); }
        if(dcmtk_vr == EVR_LO) { python_value = this->_to_python<EVR_LO>(element); }
        if(dcmtk_vr == EVR_LT) { python_value = this->_to_python<EVR_LT>(element); }
        if(dcmtk_vr == EVR_OB) { python_value = this->_to_python<EVR_OB>(element); }
        if(dcmtk_vr == EVR_OF) { python_value = this->_to_python<EVR_OF>(element); }
        if(dcmtk_vr == EVR_OW) { python_value = this->_to_python<EVR_OW>(element); }
        if(dcmtk_vr == EVR_PN) { python_value = this->_to_python<EVR_PN>(element); }
        if(dcmtk_vr == EVR_SH) { python_value = this->_to_python<EVR_SH>(element); }
        // SQ is not processed here
        if(dcmtk_vr == EVR_SL) { python_value = this->_to_python<EVR_SL>(element); }
        if(dcmtk_vr == EVR_SS) { python_value = this->_to_python<EVR_SS>(element); }
        if(dcmtk_vr == EVR_ST) { python_value = this->_to_python<EVR_ST>(element); }
        if(dcmtk_vr == EVR_TM) { python_value = this->_to_python<EVR_TM>(element); }
        if(dcmtk_vr == EVR_UI) { python_value = this->_to_python<EVR_UI>(element); }
        if(dcmtk_vr == EVR_UL) { python_value = this->_to_python<EVR_UL>(element); }
        if(dcmtk_vr == EVR_UN) { python_value = this->_to_python<EVR_UN>(element); }
        if(dcmtk_vr == EVR_US) { python_value = this->_to_python<EVR_US>(element); }
        if(dcmtk_vr == EVR_UT) { python_value = this->_to_python<EVR_UT>(element); }
    }
    
    if(python_value == NULL)
    {
        throw std::runtime_error(std::string("Unhandled VR:") + vr_name);
    }

    // Build the tag
    PyObject* key = Py_BuildValue("(II)", element->getGTag(), element->getETag());
    PyObject* tag_args = Py_BuildValue("(O)", key);
    PyObject* tag = PyObject_CallObject(this->_medipy_io_dicom_Tag, tag_args);
    Py_DECREF(tag_args);
    Py_DECREF(key);
    
    // Build the value
    std::map<DcmEVR, PyObject *>::const_iterator const vr_it = 
        this->_medipy_io_dicom_vr.find(dcmtk_vr);
    if(vr_it==this->_medipy_io_dicom_vr.end())
    {
        throw std::runtime_error(std::string("Unknown MediPy VR:") + vr_name);
    }
    PyObject* MediPyVR = vr_it->second;

    PyObject* args = Py_BuildValue("(O)", python_value);
    PyObject* typed_value = PyObject_CallObject(MediPyVR, args);
    
    Py_DECREF(args);
    Py_DECREF(MediPyVR);
    Py_DECREF(python_value);

    // Update the dictionary
    PyDict_SetItem(python_dataset, tag, typed_value);
    Py_DECREF(tag);
    Py_DECREF(typed_value);
}
