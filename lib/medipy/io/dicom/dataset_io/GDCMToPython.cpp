/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "GDCMToPython.h"

#include <Python.h>

#include <algorithm>
#include <map>
#include <stdexcept>
#include <string>

#include <gdcmGlobal.h>

#include <gdcmDataElement.h>
#include <gdcmDataSet.h>
#include <gdcmDictEntry.h>
#include <gdcmDicts.h>
#include <gdcmSequenceOfItems.h>
#include <gdcmTag.h>
#include <gdcmVR.h>

template<typename TIterator1, typename TIterator2>
TIterator1 find_first_not_of(TIterator1 first1, TIterator1 last1,
                            TIterator2 first2, TIterator2 last2)
{
    TIterator1 it=first1;
    while(it != last1)
    {
        if(std::find(first2, last2, *it) == last2)
        {
            break;
        }
        else
        {
            ++it;
        }
    }

    return it;
}

template<typename TIterator1, typename TIterator2>
TIterator1 find_last_not_of(TIterator1 first1, TIterator1 last1,
                            TIterator2 first2, TIterator2 last2)
{
    if(first1 == last1)
    {
        // Empty sequence
        return first1;
    }

    // If nothing is found,
    TIterator1 result=last1;
    // Start at the last element
    TIterator1 it=last1;
    --it;

    while(it != first1)
    {
        if(std::find(first2, last2, *it) == last2)
        {
            result = it;
            break;
        }
        --it;
    }

    // First element of the sequence
    if(it == first1)
    {
        if(std::find(first2, last2, *it) == last2)
        {
            result = it;
        }
    }

    return result;
}

GDCMToPython
::GDCMToPython()
: _medipy_io_dicom_vr(), _medipy_io_dicom_DataSet(NULL), 
  _medipy_io_dicom_Tag(NULL), _specific_character_set(""), _python_encoding(""),
  _pixel_representation(0)
{
    PyObject* modules = PyImport_GetModuleDict(); // Borrowed reference
    PyObject* medipy = PyMapping_GetItemString(modules, const_cast<char*>("medipy")); // New reference
    PyObject* medipy_io = PyObject_GetAttrString(medipy, "io"); // New reference
    PyObject* medipy_io_dicom = PyObject_GetAttrString(medipy_io, "dicom"); // New reference
    
    this->_medipy_io_dicom_DataSet = PyObject_GetAttrString(medipy_io_dicom, "DataSet"); // New reference
    this->_medipy_io_dicom_Tag = PyObject_GetAttrString(medipy_io_dicom, "Tag"); // New reference
    
    // New references
    this->_medipy_io_dicom_vr[gdcm::VR::AE] = PyObject_GetAttrString(medipy_io_dicom, "AE");
    this->_medipy_io_dicom_vr[gdcm::VR::AS] = PyObject_GetAttrString(medipy_io_dicom, "AS");
    this->_medipy_io_dicom_vr[gdcm::VR::AT] = PyObject_GetAttrString(medipy_io_dicom, "AT");
    this->_medipy_io_dicom_vr[gdcm::VR::CS] = PyObject_GetAttrString(medipy_io_dicom, "CS");
    this->_medipy_io_dicom_vr[gdcm::VR::DA] = PyObject_GetAttrString(medipy_io_dicom, "DA");
    this->_medipy_io_dicom_vr[gdcm::VR::DS] = PyObject_GetAttrString(medipy_io_dicom, "DS");
    this->_medipy_io_dicom_vr[gdcm::VR::DT] = PyObject_GetAttrString(medipy_io_dicom, "DT");
    this->_medipy_io_dicom_vr[gdcm::VR::FD] = PyObject_GetAttrString(medipy_io_dicom, "FD");
    this->_medipy_io_dicom_vr[gdcm::VR::FL] = PyObject_GetAttrString(medipy_io_dicom, "FL");
    this->_medipy_io_dicom_vr[gdcm::VR::IS] = PyObject_GetAttrString(medipy_io_dicom, "IS");
    this->_medipy_io_dicom_vr[gdcm::VR::LO] = PyObject_GetAttrString(medipy_io_dicom, "LO");
    this->_medipy_io_dicom_vr[gdcm::VR::LT] = PyObject_GetAttrString(medipy_io_dicom, "LT");
    this->_medipy_io_dicom_vr[gdcm::VR::OB] = PyObject_GetAttrString(medipy_io_dicom, "OB");
    this->_medipy_io_dicom_vr[gdcm::VR::OF] = PyObject_GetAttrString(medipy_io_dicom, "OF");
    this->_medipy_io_dicom_vr[gdcm::VR::OW] = PyObject_GetAttrString(medipy_io_dicom, "OW");
    this->_medipy_io_dicom_vr[gdcm::VR::PN] = PyObject_GetAttrString(medipy_io_dicom, "PN");
    this->_medipy_io_dicom_vr[gdcm::VR::SH] = PyObject_GetAttrString(medipy_io_dicom, "SH");
    this->_medipy_io_dicom_vr[gdcm::VR::SQ] = PyObject_GetAttrString(medipy_io_dicom, "SQ");
    this->_medipy_io_dicom_vr[gdcm::VR::SL] = PyObject_GetAttrString(medipy_io_dicom, "SL");
    this->_medipy_io_dicom_vr[gdcm::VR::SS] = PyObject_GetAttrString(medipy_io_dicom, "SS");
    this->_medipy_io_dicom_vr[gdcm::VR::ST] = PyObject_GetAttrString(medipy_io_dicom, "ST");
    this->_medipy_io_dicom_vr[gdcm::VR::TM] = PyObject_GetAttrString(medipy_io_dicom, "TM");
    this->_medipy_io_dicom_vr[gdcm::VR::UI] = PyObject_GetAttrString(medipy_io_dicom, "UI");
    this->_medipy_io_dicom_vr[gdcm::VR::UL] = PyObject_GetAttrString(medipy_io_dicom, "UL");
    this->_medipy_io_dicom_vr[gdcm::VR::UN] = PyObject_GetAttrString(medipy_io_dicom, "UN");
    this->_medipy_io_dicom_vr[gdcm::VR::US] = PyObject_GetAttrString(medipy_io_dicom, "US");
    this->_medipy_io_dicom_vr[gdcm::VR::UT] = PyObject_GetAttrString(medipy_io_dicom, "UT");
    
    Py_DECREF(medipy_io_dicom);
    Py_DECREF(medipy_io);
    Py_DECREF(medipy);
}

GDCMToPython
::GDCMToPython(GDCMToPython const & other)
: _medipy_io_dicom_vr(other._medipy_io_dicom_vr), 
  _medipy_io_dicom_DataSet(other._medipy_io_dicom_DataSet), 
  _medipy_io_dicom_Tag(other._medipy_io_dicom_Tag),
  _specific_character_set(""), _python_encoding(""),
  _pixel_representation(other._pixel_representation)
{
    for(VRMap::iterator it=this->_medipy_io_dicom_vr.begin();
        it!=this->_medipy_io_dicom_vr.end(); ++it)
    {
        Py_XINCREF(it->second);
    }
    Py_XINCREF(this->_medipy_io_dicom_DataSet);
    Py_XINCREF(this->_medipy_io_dicom_Tag);
    
    this->set_specific_character_set(other.get_specific_character_set());
}

GDCMToPython
::~GDCMToPython()
{
    for(VRMap::iterator it=this->_medipy_io_dicom_vr.begin();
        it!=this->_medipy_io_dicom_vr.end(); ++it)
    {
        Py_XDECREF(it->second);
    }
    Py_XDECREF(this->_medipy_io_dicom_DataSet);
    Py_XDECREF(this->_medipy_io_dicom_Tag);
}

GDCMToPython &
GDCMToPython
::operator=(GDCMToPython const & other)
{
    if(this != &other)
    {
        for(VRMap::iterator it=this->_medipy_io_dicom_vr.begin();
            it!=this->_medipy_io_dicom_vr.end(); ++it)
        {
            Py_XDECREF(it->second);
        }
        Py_XDECREF(this->_medipy_io_dicom_DataSet);
        Py_XDECREF(this->_medipy_io_dicom_Tag);
        
        this->_medipy_io_dicom_vr = other._medipy_io_dicom_vr;
        this->_medipy_io_dicom_DataSet = other._medipy_io_dicom_DataSet;
        this->_medipy_io_dicom_Tag = other._medipy_io_dicom_Tag;
        
        for(VRMap::iterator it=this->_medipy_io_dicom_vr.begin();
            it!=this->_medipy_io_dicom_vr.end(); ++it)
        {
            Py_XDECREF(it->second);
        }
        Py_XINCREF(this->_medipy_io_dicom_DataSet);
        Py_XINCREF(this->_medipy_io_dicom_Tag);
        
        this->set_specific_character_set(other.get_specific_character_set());
        this->_pixel_representation = other._pixel_representation;
    }
    return *this;
}

std::string const &
GDCMToPython
::get_specific_character_set() const
{
    return this->_specific_character_set;
}

void
GDCMToPython
::set_specific_character_set(std::string const & charset)
{
    if(charset != this->get_specific_character_set())
    {
        this->_specific_character_set = charset;
        
        if(charset == "") this->_python_encoding = "ascii";
        else if(charset == "ISO_IR 100") this->_python_encoding = "latin_1";
        else if(charset == "ISO_IR 101") this->_python_encoding = "iso8859_2";
        else if(charset == "ISO_IR 109") this->_python_encoding = "iso8859_3";
        else if(charset == "ISO_IR 110") this->_python_encoding = "iso8859_4";
        else if(charset == "ISO_IR 144") this->_python_encoding = "iso8859_5";
        else if(charset == "ISO_IR 127") this->_python_encoding = "iso8859_6";
        else if(charset == "ISO_IR 126") this->_python_encoding = "iso8859_7";
        else if(charset == "ISO_IR 138") this->_python_encoding = "iso8859_8";
        else if(charset == "ISO_IR 148") this->_python_encoding = "iso8859_9";
        else if(charset == "ISO_IR 13") this->_python_encoding = "iso2022_jp";
        // CP874 seems to be a superset of TIS-620/ISO-IR-166 (e.g. presence of 
        // the euro sign in the CP874 at an unassigned place of TIS-620), but we
        // should get away with it.
        else if(charset == "ISO_IR 166") this->_python_encoding = "cp874";
        // Single-byte character sets with code extensions (Table C.12-3)
        // ISO 2022 IR 6
        // ISO 2022 IR 100
        // ISO 2022 IR 101
        // ISO 2022 IR 109
        // ISO 2022 IR 110
        // ISO 2022 IR 144
        // ISO 2022 IR 127
        // ISO 2022 IR 126
        // ISO 2022 IR 138
        // ISO 2022 IR 148
        // ISO 2022 IR 113
        // ISO 2022 IR 166
        // Multi-byte character sets without code extensions (Table C.12-4)
        // ISO 2022 IR 87 // Kanji
        // ISO 2022 IR 159 // Supplementary Kanji
        // ISO 2022 IR 149 // Hangul, Hanja
        // Multi-byte character sets without code extensions (Table C.12-5)
        else if(charset == "ISO_IR 192") this->_python_encoding = "utf_8";
        else if(charset == "GB18030") this->_python_encoding = "gb18030";
    }
}

PyObject *
GDCMToPython
::operator()(gdcm::DataSet const & dataset)
{
    PyObject* python_dataset = PyObject_CallFunction(
        this->_medipy_io_dicom_DataSet, "()");
    
    for(gdcm::DataSet::ConstIterator it=dataset.Begin(); it!=dataset.End(); ++it)
    {
        gdcm::DataElement const & gdcm_element = *it;
        gdcm::Tag const & gdcm_tag=gdcm_element.GetTag();
        if(gdcm_tag.GetElement() == 0)
        {
            // Group length, do nothing
            continue;
        }
        else if(gdcm_tag == 0x00080005)
        {
            // Specific Character Set: setup internal iconv converter
            if(gdcm_element.GetByteValue() != NULL && 
               gdcm_element.GetByteValue()->GetPointer() != NULL)
            {
                this->set_specific_character_set(std::string(
                    gdcm_element.GetByteValue()->GetPointer(), 
                    gdcm_element.GetVL()));
            }
            else
            {
                this->set_specific_character_set("");
            }
        }
        else if(gdcm_tag == 0x00280103)
        {
            // Pixel representation: used to determine VR for implicit transfer
            // syntaxes
            if(gdcm_element.GetByteValue() != NULL && 
               gdcm_element.GetByteValue()->GetPointer() != NULL)
            {
                this->_pixel_representation = 
                    *reinterpret_cast<unsigned short const *>(
                        gdcm_element.GetByteValue()->GetPointer());
            }
        }
        else
        {
            PyObject * python_tag = this->operator()(gdcm_tag);
            PyObject * python_element = this->operator()(gdcm_element, dataset);
            // Update the dictionary
            PyDict_SetItem(python_dataset, python_tag, python_element);
            Py_DECREF(python_tag);
            Py_DECREF(python_element);
        }
    }
    
    return python_dataset;
}

PyObject * 
GDCMToPython
::operator()(gdcm::DataElement const & gdcm_element,
             gdcm::DataSet const & dataset) const
{
    PyObject* value=NULL;

    gdcm::VR gdcm_vr = gdcm_element.GetVR();
    
    if(gdcm_vr == gdcm::VR::INVALID)
    {
        gdcm::Global const & global = gdcm::GlobalInstance;
        gdcm::Dicts const & dicts = global.GetDicts();
        
        std::string strowner;
        const char *owner = 0;
        gdcm::Tag const & tag = gdcm_element.GetTag();
        if(tag.IsPrivate() && !tag.IsPrivateCreator())
        {
            strowner = dataset.GetPrivateCreator(tag);
            owner = strowner.c_str();
        }
        gdcm::DictEntry const & entry = dicts.GetDictEntry(tag, owner);
        gdcm_vr = entry.GetVR();
        
        if(gdcm_vr == gdcm::VR::US_SS)
        {
            if(this->_pixel_representation == 0)
            {
                gdcm_vr = gdcm::VR::US;
            }
            else if(this->_pixel_representation == 1)
            {
                gdcm_vr = gdcm::VR::SS;
            }
            else
            {
                throw std::runtime_error("Unknown Pixel Representation");
            }
        }
        else if(gdcm_vr == gdcm::VR::OB_OW)
        {
            // FIXME : this should only be correct for Images, and should use
            // Bits Allocated
            gdcm_vr = gdcm::VR::OW;
        }
        else if(gdcm_vr == gdcm::VR::INVALID)
        {
            // Could not get VR from dictionary, assume an unknown private tag
            gdcm_vr = gdcm::VR::OB;
        }
    }
    
    if(gdcm_vr != gdcm::VR::INVALID && gdcm_vr != gdcm::VR::SQ && 
       (gdcm_element.GetByteValue() == NULL ||
        gdcm_element.GetByteValue()->GetPointer() == NULL))
    {
        value = Py_None;
        Py_INCREF(value);
    }
    else if(gdcm_vr & (gdcm::VR::OB | gdcm::VR::OF | gdcm::VR::OW | gdcm::VR::UN))
    {
        // Return str, to be used as sequence of bytes
        value = PyString_FromStringAndSize(
            gdcm_element.GetByteValue()->GetPointer(), gdcm_element.GetVL());
    }
    else if(gdcm_vr == gdcm::VR::SQ)
    {
        gdcm::SmartPointer<gdcm::SequenceOfItems> sequence = 
            gdcm_element.GetValueAsSQ();

        if(sequence == NULL)
        {
            value = PyList_New(0);
        }
        else
        {
            value = PyList_New(sequence->GetNumberOfItems());

            unsigned int index=0;
            for(gdcm::SequenceOfItems::ConstIterator sequence_it=sequence->Begin();
                sequence_it!=sequence->End(); ++sequence_it, ++index)
            {
                GDCMToPython converter(*this);
                PyObject * python_item = converter(sequence_it->GetNestedDataSet());
                if(python_item == NULL)
                {
                    Py_DECREF(value);
                    value = NULL;
                    break;
                }
                PyList_SetItem(value, index, python_item);
                // Don't Py_DECREF(python_item) : according to Python doc :
                // This function “steals” a reference to item and discards a
                // reference to an item already in the list at the affected position.
            }
        }
    }
    else if(gdcm::VR::IsBinary(gdcm_vr) || gdcm::VR::IsASCII(gdcm_vr))
    {
        gdcm::ByteValue const * byte_value = gdcm_element.GetByteValue();

        unsigned int count = 0;
        if(byte_value == NULL)
        {
            count = 0;
        }
        else 
        {
            if(gdcm::VR::IsBinary(gdcm_vr))
            {
                count = byte_value->GetLength()/gdcm_vr.GetSize();
            }
            else if(gdcm_vr & (gdcm::VR::LT | gdcm::VR::ST | gdcm::VR::UT))
            {
                // LT, ST and UT may not be multi-valued and 
                // gdcm::VM::GetNumberOfElementsFromArray does not use the VR !
                count = 1;
            }
            else 
            {
                count = gdcm::VM::GetNumberOfElementsFromArray(
                    byte_value->GetPointer(), byte_value->GetLength());
            }
        }
        
        if(count == 0)
        {
            value = Py_None;
            Py_INCREF(value);
        }
        else if(count == 1)
        {
            char const * begin = byte_value->GetPointer();
            char const * end = begin+byte_value->GetLength();
            value = this->_to_python_value(begin, end, gdcm_vr);
        }
        else // count > 1
        {
            value = PyList_New(count);
            char const * begin = byte_value->GetPointer();
            char const * end = begin+byte_value->GetLength();
            for(unsigned int i=0; i<count; ++i)
            {
                char const * item_end = begin+this->_get_length(begin, end, gdcm_vr);
                PyObject* item = this->_to_python_value(begin, item_end, gdcm_vr);
                PyList_SetItem(value, i, item);
                // Don't Py_DECREF(python_item) : according to Python doc :
                // This function “steals” a reference to item and discards a
                // reference to an item already in the list at the affected position.
                begin = item_end;
                if(gdcm::VR::IsASCII(gdcm_vr))
                {
                    // Skip the backslash
                    ++begin;
                }
            }
        }
    }
    
    if(value == NULL)
    {
        std::ostringstream message;
        message << "Could not convert data to VR " << gdcm_vr 
                << " for tag " << gdcm_element.GetTag();
        throw std::runtime_error(message.str());
    }
    
    // Build the value
    VRMap::const_iterator const vr_it = this->_medipy_io_dicom_vr.find(gdcm_vr);
    if(vr_it==this->_medipy_io_dicom_vr.end())
    {
        std::ostringstream message;
        message << "Unknown MediPy VR on tag " << gdcm_element.GetTag() 
                << " : " << gdcm_vr;
        throw std::runtime_error(message.str());
    }
    PyObject* MediPyVR = vr_it->second; // Borrowed reference

    PyObject* python_element = PyObject_CallFunction(MediPyVR, "(O)", value);
    Py_DECREF(value);
    
    if(python_element == NULL)
    {
        std::ostringstream message;
        message << "Could not build Python element of VR " << gdcm_vr 
                << "for tag " << gdcm_element.GetTag();
        throw std::runtime_error(message.str());
    }

    return python_element;
}

PyObject *
GDCMToPython
::operator()(gdcm::Tag const & gdcm_tag) const
{
    PyObject* python_tag = PyObject_CallFunction(
        this->_medipy_io_dicom_Tag, "(II)", 
        gdcm_tag.GetGroup(), gdcm_tag.GetElement());
    if(python_tag == NULL)
    {
        std::ostringstream message;
        message << "Could not create Python tag from " << gdcm_tag;
        throw std::runtime_error(message.str());
    }
    return python_tag;
}

PyObject*
GDCMToPython
::_to_python_value(char const * begin, char const * end, gdcm::VR const & vr) const
{
    static std::string const whitespace(" \0", 2);
    static std::string const whitespace_and_backslash(" \0\\", 3);

    if(vr & (gdcm::VR::AE | gdcm::VR::AS | gdcm::VR::CS | gdcm::VR::DA |
             gdcm::VR::DT | gdcm::VR::TM | gdcm::VR::UI))
    {
        // ASCII text content, leading and trailing whitespaces are not significant
        char const * first = find_first_not_of(
            begin, end,
            whitespace_and_backslash.begin(), whitespace_and_backslash.end());

        unsigned int size=0;
        if(first != end)
        {
            char const * last = find_last_not_of(
                begin, end,
                whitespace_and_backslash.begin(), whitespace_and_backslash.end());
            size = last-first+1;
        }
        return PyString_FromStringAndSize(first, size);
    }
    else if(vr & (gdcm::VR::LO | gdcm::VR::LT | gdcm::VR::PN | gdcm::VR::SH |
                  gdcm::VR::ST | gdcm::VR::UT))
    {
        // Non-ASCII text content, translate to Python Unicode
        char const * first;
        if(vr & (gdcm::VR::LT | gdcm::VR::ST | gdcm::VR::UT))
        {
            // Leading spaces are significant for LT, ST, and UT
            first = find_first_not_of(
                begin, end,
                whitespace.begin(), whitespace.end());
        }
        else
        {
            first = find_first_not_of(
                begin, end,
                whitespace_and_backslash.begin(), whitespace_and_backslash.end());
        }

        unsigned int size=0;
        if(first != end)
        {
            char const * last;
            if(vr & (gdcm::VR::LT | gdcm::VR::ST | gdcm::VR::UT))
            {
                // LT, ST and UT may not be multi-valued
                last = end-1;
            }
            else
            {
                last = find_last_not_of(
                    begin, end,
                    whitespace_and_backslash.begin(), whitespace_and_backslash.end());
            }
            size = last-first+1;
        }

#define TRY_DECODE(o, e) \
    if(object==NULL) \
    { \
        /* Avoid propagating the error. */ \
        PyErr_Clear(); \
        object = PyUnicode_Decode(first, size, e, "strict"); \
    }
        PyObject* object = NULL;
        TRY_DECODE(object, this->_python_encoding.c_str());
        TRY_DECODE(object, "ascii");
        TRY_DECODE(object, "latin_1");
        TRY_DECODE(object, "iso8859_2");
        TRY_DECODE(object, "iso8859_3");
        TRY_DECODE(object, "iso8859_4");
        TRY_DECODE(object, "iso8859_5");
        TRY_DECODE(object, "iso8859_6");
        TRY_DECODE(object, "iso8859_7");
        TRY_DECODE(object, "iso8859_8");
        TRY_DECODE(object, "iso8859_9");
        TRY_DECODE(object, "iso2022_jp");
        TRY_DECODE(object, "utf_8");
        TRY_DECODE(object, "gb18030");
#undef TRY_DECODE
        if(object == NULL)
        {
            throw std::runtime_error("Could not decode string");
        }

        return object;
    }
    else if(vr == gdcm::VR::AT)
    {
        uint16_t e0 = *reinterpret_cast<uint16_t const *>(begin);
        uint16_t e1 = *(1+reinterpret_cast<uint16_t const *>(begin));
        return PyObject_CallFunction(this->_medipy_io_dicom_Tag, "(II)", e0, e1);
    }
    else if(vr == gdcm::VR::DS)
    {
        char * old_numeric = setlocale(LC_NUMERIC, NULL);
        setlocale(LC_NUMERIC, "C");
        char* endptr;
        double const d = std::strtod(begin, &endptr);
        setlocale(LC_NUMERIC, old_numeric);
        if(endptr == begin)
        {
            PyErr_SetString(PyExc_Exception, "Cannot parse DS");
            return NULL;
        }

        return PyFloat_FromDouble(d);
    }
    else if(vr == gdcm::VR::FD)
    {
        // TODO : correct 64 bits type
        return PyFloat_FromDouble(*reinterpret_cast<double const *>(begin));
    }
    else if(vr == gdcm::VR::FL)
    {
        // TODO : correct 32 bits type
        return PyFloat_FromDouble(*reinterpret_cast<float const *>(begin));
    }
    else if(vr == gdcm::VR::IS)
    {
        char * old_numeric = setlocale(LC_NUMERIC, NULL);
        setlocale(LC_NUMERIC, "C");
        char* endptr;
        long const d = std::strtol(begin, &endptr, 10);
        setlocale(LC_NUMERIC, old_numeric);
        if(endptr == begin)
        {
            PyErr_SetString(PyExc_Exception, "Cannot parse IS");
            return NULL;
        }
        return PyInt_FromLong(d);
    }
    else if(vr == gdcm::VR::SL)
    {
        return PyInt_FromLong(*reinterpret_cast<int32_t const *>(begin));
    }
    else if(vr == gdcm::VR::SS)
    {
        return PyInt_FromLong(*reinterpret_cast<int16_t const *>(begin));
    }
    else if(vr == gdcm::VR::UL)
    {
        return PyInt_FromLong(*reinterpret_cast<uint32_t const *>(begin));
    }
    else if(vr == gdcm::VR::US)
    {
        return PyInt_FromLong(*reinterpret_cast<uint16_t const *>(begin));
    }
    else
    {
        std::ostringstream message;
        message << "Cannot parse VR " << vr;
        PyErr_SetString(PyExc_Exception, message.str().c_str());
        return NULL;
    }
}

unsigned long
GDCMToPython
::_get_length(char const * begin, char const * end, gdcm::VR const & vr) const
{
    if(vr & (gdcm::VR::AE | gdcm::VR::AS | gdcm::VR::CS | gdcm::VR::DA |
             gdcm::VR::DS | gdcm::VR::DT | gdcm::VR::IS | gdcm::VR::LO | 
             gdcm::VR::PN | gdcm::VR::SH | gdcm::VR::TM | gdcm::VR::UI))
    {
        char const * value_end = std::find(begin, end, '\\');
        return (value_end-begin);
    }
    else if(vr & (gdcm::VR::LT | gdcm::VR::ST | gdcm::VR::UT))
    {
        // LT, ST and UT may not be multi-valued
        return end-begin;
    }
    else if(vr == gdcm::VR::AT)
    {
        return 4;
    }
    else if(vr == gdcm::VR::FD)
    {
        return 8;
    }
    else if(vr == gdcm::VR::FL)
    {
        return 4;
    }
    else if(vr == gdcm::VR::SL)
    {
        return 4;
    }
    else if(vr == gdcm::VR::SS)
    {
        return 2;
    }
    else if(vr == gdcm::VR::UL)
    {
        return 4;
    }
    else if(vr == gdcm::VR::US)
    {
        return 2;
    }
    else
    {
        std::ostringstream message;
        message << "Cannot get length for VR " << vr;
        PyErr_SetString(PyExc_Exception, message.str().c_str());
        return 0;
    }

}
