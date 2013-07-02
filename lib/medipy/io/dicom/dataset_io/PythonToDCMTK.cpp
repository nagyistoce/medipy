/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "PythonToDCMTK.h"

#include <Python.h>

#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>

#include <dcmtk/config/osconfig.h>
#include <dcmtk/dcmdata/dctk.h>

PythonToDCMTK
::PythonToDCMTK()
{
    // Nothing else.
}

PythonToDCMTK
::PythonToDCMTK(PythonToDCMTK const & )
{
    // Nothing else.
}

PythonToDCMTK
::~PythonToDCMTK()
{
    // Nothing to do.
}

PythonToDCMTK &
PythonToDCMTK
::operator=(PythonToDCMTK const & )
{
    // Nothing else
    
    return *this;
}

DcmDataset
PythonToDCMTK
::operator()(PyObject * python_dataset)
{
    DcmDataset dataset;
    
    PyObject * tags = PyDict_Keys(python_dataset); // New reference
    PyList_Sort(tags);
    for(unsigned int index=0; index<PyList_Size(tags); ++index)
    {
        PyObject * tag = PyList_GetItem(tags, index);
        PyObject * value = PyDict_GetItem(python_dataset, tag);
        this->_add_element(tag, value, dataset);
    }
    
    Py_DECREF(tags);
    
    return dataset;
}

/*******************************************************************************
 * Specializations of PythonToDCMTK::_to_dcmtk for the different VRs.
 ******************************************************************************/

template<DcmEVR VVR> struct ElementTrait;
template<> struct ElementTrait<EVR_AE> { typedef DcmApplicationEntity Type; };
template<> struct ElementTrait<EVR_AS> { typedef DcmAgeString Type; };
template<> struct ElementTrait<EVR_AT> { typedef DcmAttributeTag Type; };
template<> struct ElementTrait<EVR_CS> { typedef DcmCodeString Type; };
template<> struct ElementTrait<EVR_DA> { typedef DcmDate Type; };
template<> struct ElementTrait<EVR_DS> { typedef DcmDecimalString Type; };
template<> struct ElementTrait<EVR_DT> { typedef DcmDateTime Type; };
template<> struct ElementTrait<EVR_FL> { typedef DcmFloatingPointSingle Type; };
template<> struct ElementTrait<EVR_FD> { typedef DcmFloatingPointDouble Type; };
template<> struct ElementTrait<EVR_IS> { typedef DcmIntegerString Type; };
template<> struct ElementTrait<EVR_LO> { typedef DcmLongString Type; };
template<> struct ElementTrait<EVR_LT> { typedef DcmLongText Type; };
template<> struct ElementTrait<EVR_OB> { typedef DcmOtherByteOtherWord Type; };
template<> struct ElementTrait<EVR_OF> { typedef DcmOtherByteOtherWord Type; };
template<> struct ElementTrait<EVR_OW> { typedef DcmOtherByteOtherWord Type; };
template<> struct ElementTrait<EVR_PN> { typedef DcmPersonName Type; };
template<> struct ElementTrait<EVR_SH> { typedef DcmShortString Type; };
template<> struct ElementTrait<EVR_SL> { typedef DcmSignedLong Type; };
template<> struct ElementTrait<EVR_SQ> { typedef DcmSequenceOfItems Type; };
template<> struct ElementTrait<EVR_SS> { typedef DcmSignedShort Type; };
template<> struct ElementTrait<EVR_ST> { typedef DcmShortText Type; };
template<> struct ElementTrait<EVR_TM> { typedef DcmTime Type; };
template<> struct ElementTrait<EVR_UI> { typedef DcmUniqueIdentifier Type; };
template<> struct ElementTrait<EVR_UL> { typedef DcmUnsignedLong Type; };
template<> struct ElementTrait<EVR_UN> { typedef DcmOtherByteOtherWord Type; };
template<> struct ElementTrait<EVR_US> { typedef DcmUnsignedShort Type; };
template<> struct ElementTrait<EVR_UT> { typedef DcmUnlimitedText Type; };

template<DcmEVR VVR>
DcmElement * 
createElement(DcmTag const & tag)
{
    // Set correct VR (useful for private tags)
    DcmTag newTag(tag);
    newTag.setVR(VVR);
    
    return new typename ElementTrait<VVR>::Type(newTag);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_AE>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_text(python_element, false, ' ', dataset, createElement<EVR_AE>(tag));
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_AS>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_text(python_element, false, ' ', dataset, createElement<EVR_AS>(tag));
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_AT>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    DcmAttributeTag * dcmtk_element = 
        dynamic_cast<DcmAttributeTag*>(createElement<EVR_AT>(tag));
    if(PyList_Check(python_element))
    {
        for(Py_ssize_t index=0; index<PyList_Size(python_element); ++index)
        {
            PyObject * python_tag = PyList_GetItem(python_element, index);
            unsigned long const numeric_tag = PyLong_AsLong(python_tag);
            DcmTagKey const tag(numeric_tag>>16, numeric_tag&0xffff);
            dcmtk_element->putTagVal(tag, index);
        }
    }
    else
    {
        PyObject * python_tag = python_element;
        unsigned long const numeric_tag = PyLong_AsLong(python_tag);
        DcmTagKey const tag(numeric_tag>>16, numeric_tag&0xffff);
        dcmtk_element->putTagVal(tag, 0);
    }
    dataset.insert(dcmtk_element);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_CS>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_text(python_element, false, ' ', dataset, createElement<EVR_CS>(tag));
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_DA>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_text(python_element, false, ' ', dataset, createElement<EVR_DA>(tag));
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_DS>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_number_string(python_element, dataset, createElement<EVR_DS>(tag));
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_DT>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_text(python_element, false, ' ', dataset, createElement<EVR_DT>(tag));
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_FD>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_binary(python_element, PyFloat_AsDouble, dataset, 
                     createElement<EVR_FD>(tag), &DcmElement::putFloat64);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_FL>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_binary(python_element, PyFloat_AsDouble, dataset, 
                     createElement<EVR_FL>(tag), &DcmElement::putFloat32);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_IS>(PyObject * python_element, DcmDataset & dataset,
                   DcmTag const & tag) const
{
    this->_to_number_string(python_element, dataset, createElement<EVR_IS>(tag));
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_LO>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_text(python_element, true, ' ', dataset, createElement<EVR_LO>(tag));
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_LT>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_text(python_element, true, ' ', dataset, createElement<EVR_LT>(tag));
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_OB>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_raw_8(python_element, dataset, createElement<EVR_OB>(tag));
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_OF>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_raw_8(python_element, dataset, createElement<EVR_OF>(tag));
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_OW>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_raw_16(python_element, dataset, createElement<EVR_OW>(tag));
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_PN>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_text(python_element, true, ' ', dataset, createElement<EVR_PN>(tag));
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_SH>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_text(python_element, true, ' ', dataset, createElement<EVR_SH>(tag));
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_SL>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_binary(python_element, PyInt_AsLong, dataset,
                     createElement<EVR_SL>(tag), &DcmElement::putSint32);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_SQ>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_sequence(python_element, dataset, 
        dynamic_cast<DcmSequenceOfItems*>(createElement<EVR_SQ>(tag)));
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_SS>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_binary(python_element, PyInt_AsLong, dataset,
                     createElement<EVR_SS>(tag), &DcmElement::putSint16);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_ST>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_text(python_element, true, ' ', dataset, createElement<EVR_ST>(tag));
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_TM>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_text(python_element, false, ' ', dataset, createElement<EVR_TM>(tag));
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_UI>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_text(python_element, false, '\0', dataset, createElement<EVR_UI>(tag));
}


template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_UL>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_binary(python_element, PyInt_AsLong, dataset,
                     createElement<EVR_UL>(tag), &DcmElement::putUint32);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_UN>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_raw_8(python_element, dataset, createElement<EVR_UN>(tag));
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_US>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_binary(python_element, PyInt_AsLong, dataset,
                     createElement<EVR_US>(tag), &DcmElement::putUint16);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_UT>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_text(python_element, true, ' ', dataset, createElement<EVR_UT>(tag));
}

/*******************************************************************************
 * End of specializations of BSONToDataSet::_to_dcmtk for the different VRs.
 ******************************************************************************/

void
PythonToDCMTK
::_add_element(PyObject * python_tag, PyObject * python_value, DcmDataset & dataset)
{
    // Get the tag
    unsigned long const numeric_tag = PyLong_AsLong(python_tag);
    DcmTag const tag(numeric_tag>>16, numeric_tag&0xffff);
    
    if(tag == DCM_SpecificCharacterSet)
    {
        // We encode everything to UTF-8, so we overwrite the Specific Character Set
        dataset.putAndInsertOFStringArray(tag, "ISO_IR 192");
        return;
    }
    
    // Get the VR
    PyObject* value_type = PyObject_GetAttrString(python_value, "__class__"); // New reference
    PyObject* python_vr = PyObject_GetAttrString(value_type, "__name__"); // New reference
    DcmVR const vr(PyString_AsString(python_vr));
    DcmEVR const evr(vr.getValidEVR());
    Py_DECREF(python_vr);
    Py_DECREF(value_type);

    // Get the value
    PyObject* nested_value = PyObject_GetAttrString(python_value, "value"); // New reference
    if(nested_value == NULL)
    {
        throw std::runtime_error("Cannot access element value");
    }

    if(nested_value != Py_None)
    {
        if(evr == EVR_AE) this->_to_dcmtk<EVR_AE>(nested_value, dataset, tag);
        else if(evr == EVR_AS) this->_to_dcmtk<EVR_AS>(nested_value, dataset, tag);
        else if(evr == EVR_AT) this->_to_dcmtk<EVR_AT>(nested_value, dataset, tag);
        else if(evr == EVR_CS) this->_to_dcmtk<EVR_CS>(nested_value, dataset, tag);
        else if(evr == EVR_DA) this->_to_dcmtk<EVR_DA>(nested_value, dataset, tag);
        else if(evr == EVR_DS) this->_to_dcmtk<EVR_DS>(nested_value, dataset, tag);
        else if(evr == EVR_DT) this->_to_dcmtk<EVR_DT>(nested_value, dataset, tag);
        else if(evr == EVR_FD) this->_to_dcmtk<EVR_FD>(nested_value, dataset, tag);
        else if(evr == EVR_FL) this->_to_dcmtk<EVR_FL>(nested_value, dataset, tag);
        else if(evr == EVR_IS) this->_to_dcmtk<EVR_IS>(nested_value, dataset, tag);
        else if(evr == EVR_LO) this->_to_dcmtk<EVR_LO>(nested_value, dataset, tag);
        else if(evr == EVR_LT) this->_to_dcmtk<EVR_LT>(nested_value, dataset, tag);
        else if(evr == EVR_OB) this->_to_dcmtk<EVR_OB>(nested_value, dataset, tag);
        else if(evr == EVR_OF) this->_to_dcmtk<EVR_OF>(nested_value, dataset, tag);
        else if(evr == EVR_OW) this->_to_dcmtk<EVR_OW>(nested_value, dataset, tag);
        else if(evr == EVR_PN) this->_to_dcmtk<EVR_PN>(nested_value, dataset, tag);
        else if(evr == EVR_SH) this->_to_dcmtk<EVR_SH>(nested_value, dataset, tag);
        else if(evr == EVR_SQ) this->_to_dcmtk<EVR_SQ>(nested_value, dataset, tag);
        else if(evr == EVR_SL) this->_to_dcmtk<EVR_SL>(nested_value, dataset, tag);
        else if(evr == EVR_SS) this->_to_dcmtk<EVR_SS>(nested_value, dataset, tag);
        else if(evr == EVR_ST) this->_to_dcmtk<EVR_ST>(nested_value, dataset, tag);
        else if(evr == EVR_TM) this->_to_dcmtk<EVR_TM>(nested_value, dataset, tag);
        else if(evr == EVR_UI) this->_to_dcmtk<EVR_UI>(nested_value, dataset, tag);
        else if(evr == EVR_UL) this->_to_dcmtk<EVR_UL>(nested_value, dataset, tag);
        else if(evr == EVR_UN) this->_to_dcmtk<EVR_UN>(nested_value, dataset, tag);
        else if(evr == EVR_US) this->_to_dcmtk<EVR_US>(nested_value, dataset, tag);
        else if(evr == EVR_UT) this->_to_dcmtk<EVR_UT>(nested_value, dataset, tag);
        else
        {
            throw std::runtime_error(std::string("Unhandled VR: ") + vr.getValidVRName());
        }
    }
    else
    {
        dataset.insertEmptyElement(tag);
    }

    Py_DECREF(nested_value);
}

/** 
 * @brief Return a Python string containing the UTF-8 encoded version of the 
 * string or unicode passed to the function.
 * @return New reference
 */
PyObject* _to_utf8(PyObject * data)
{
    PyObject * result = NULL;
    
    if(PyString_Check(data))
    {
        // Add a reference, to make sure that a new reference is always returned
        Py_INCREF(data);
        result = data;
    }
    else if(PyUnicode_Check(data))
    {
        result = PyUnicode_AsEncodedString(data, "utf-8", "strict");
    }
    
    if(result == NULL)
    {
        throw std::runtime_error("Could not encode object to UTF-8");
    }
    
    return result;
}

void
PythonToDCMTK
::_to_text(PyObject * python_value, bool use_utf8, char padding,
           DcmDataset & dataset, DcmElement * element) const
{
    std::ostringstream stream;
    stream.imbue(std::locale("C"));
    
    if(PyList_Check(python_value))
    {
        Py_ssize_t const size = PyList_Size(python_value);
        for(Py_ssize_t index=0; index<size; ++index)
        {
            PyObject * python_item = PyList_GetItem(python_value, index);
            if(use_utf8)
            {
                PyObject * decoded_string = _to_utf8(python_item);
                stream << PyString_AsString(decoded_string);
                Py_DECREF(decoded_string);
            }
            else
            {
                stream << PyString_AsString(python_item);
            }
            
            if(index != size-1)
            {
                stream << "\\";
            }
        }
    }
    else
    {
        if(use_utf8)
        {
            PyObject * decoded_string = _to_utf8(python_value);
            stream << PyString_AsString(decoded_string);
            Py_DECREF(decoded_string);
        }
        else
        {
            stream << PyString_AsString(python_value);
        }
    }
    
    OFString value(stream.str().c_str());

    if(value.size()%2!=0)
    {
        value += padding;
    }
    
    OFCondition const condition = element->putString(value.c_str());
    if(condition.bad())
    {
        throw std::runtime_error(
            std::string("Cannot convert text element: ")+condition.text());
    }
    
    dataset.insert(element);
}

template<typename TGetter, typename TInserter>
void
PythonToDCMTK
::_to_binary(PyObject * python_value, TGetter getter,
             DcmDataset & dataset, DcmElement * element, TInserter inserter) const
{
    if(PyList_Check(python_value))
    {
        for(Py_ssize_t index=0; index<PyList_Size(python_value); ++index)
        {
            PyObject * item = PyList_GetItem(python_value, index);
            
            OFCondition const condition = (element->*inserter)(getter(item), index);
            if(condition.bad())
            {
                std::ostringstream message;
                message << "Cannot convert item " << index << " of numeric element "
                        << element->getTag() 
                        << " (" << element->getTag().getVR().getValidVRName() << "): "
                        << condition.text();
                throw std::runtime_error(message.str());
            }
        }
    }
    else
    {
        OFCondition const condition = (element->*inserter)(getter(python_value), 0);
        if(condition.bad())
        {
            std::ostringstream message;
            message << "Cannot convert numeric element "
                    << element->getTag() 
                    << " (" << element->getTag().getVR().getValidVRName() << "): "
                    << condition.text();
            throw std::runtime_error(message.str());
        }
    }
    
    dataset.insert(element);
}

void
PythonToDCMTK
::_to_raw_8(PyObject * python_value, DcmDataset & dataset,
            DcmElement * element) const
{
    unsigned long size = PyString_Size(python_value)/sizeof(Uint8);
    Uint8 const * data = reinterpret_cast<Uint8 const *>(PyString_AsString(python_value));
    
    element->putUint8Array(data, size);
    dataset.insert(element);
}

void
PythonToDCMTK
::_to_raw_16(PyObject * python_value, DcmDataset & dataset,
             DcmElement * element) const
{
    unsigned long size = PyString_Size(python_value)/sizeof(Uint16);
    Uint16 const * data = reinterpret_cast<Uint16 const *>(PyString_AsString(python_value));
    
    element->putUint16Array(data, size);
    dataset.insert(element);
}

void
PythonToDCMTK
::_to_number_string(PyObject * python_value, DcmDataset & dataset,
                    DcmElement * element) const
{
    std::ostringstream stream;
    stream.imbue(std::locale("C"));
    
    if(PyList_Check(python_value))
    {
        Py_ssize_t const size = PyList_Size(python_value);
        for(Py_ssize_t index=0; index<size; ++index)
        {
            PyObject * python_item = PyList_GetItem(python_value, index);
            PyObject * python_string = PyObject_Str(python_item); // New reference
            stream << PyString_AsString(python_string);
            Py_DECREF(python_string);
            
            if(index != size-1)
            {
                stream << "\\";
            }
        }
    }
    else
    {
        PyObject * python_string = PyObject_Str(python_value); // New reference
        stream << PyString_AsString(python_string);
        Py_DECREF(python_string);
    }
    
    OFString value(stream.str().c_str());

    if(value.size()%2!=0)
    {
        value += ' ';
    }
    
    element->putOFStringArray(value);
    dataset.insert(element);
}

void 
PythonToDCMTK
::_to_sequence(PyObject * python_value, DcmDataset & dataset, DcmSequenceOfItems * element) const
{
    if(!PyList_Check(python_value))
    {
        throw std::runtime_error("Cannot convert non-list SQ");
    }
    for(Py_ssize_t index=0; index<PyList_Size(python_value); ++index)
    {
        PythonToDCMTK converter(*this);
        // DcmSequenceOfItems::append takes ownership of item, hence the new.
        // Caution : do NOT append a DcmDataset to a DcmSequenceOfItems, but a DcmItem
        DcmItem * item = new DcmItem(converter(PyList_GetItem(python_value, index)));
        element->append(item);
    }
    
    dataset.insert(element);
}
