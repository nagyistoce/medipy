/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2012
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "DataSetBridge.h"

#include <Python.h>

#include <algorithm>
#include <cstdlib>
#include <sstream>
#include <stdint.h>

#include <gdcmAttribute.h>
#include <gdcmByteValue.h>
#include <gdcmDataElement.h>
#include <gdcmDataSet.h>
#include <gdcmSequenceOfItems.h>
#include <gdcmVR.h>

#include <numpy/arrayobject.h>

// Initialize PyArray_API : since the symbol is declared as static unless
// PY_ARRAY_UNIQUE_SYMBOL is defined to a unique value and NO_IMPORT is defined
// (cf. __multiarray_api.h), it is easier to call import_array in each
// translation unit using the C API of NumPy.
namespace
{

bool initialize()
{
    // Use import_array1 instead of import_array so that we have a return value
    // even if the import fails (cf. __multiarray_api.h).
    import_array1(false);
    return true;
}

static const bool initialized=initialize();

}


DataSetBridge
::DataSetBridge()
: _python_dataset(NULL)
{
    // Nothing else.
}

DataSetBridge
::DataSetBridge(gdcm::DataSet const & dataset)
: _python_dataset(NULL)
{
    this->set_gdcm(dataset);
}

DataSetBridge
::DataSetBridge(PyObject* dataset)
{
    this->set_python(dataset);
}

gdcm::DataSet const &
DataSetBridge
::get_gdcm() const
{
    return this->_gdcm_dataset;
}

void
DataSetBridge
::set_gdcm(gdcm::DataSet const & dataset)
{
    this->_gdcm_dataset = dataset;
}

PyObject*
DataSetBridge
::get_python() const
{
    return this->_python_dataset;
}

void
DataSetBridge
::set_python(PyObject* dataset)
{
    if(this->_python_dataset != NULL)
    {
        Py_DECREF(this->_python_dataset);
    }
    this->_python_dataset = dataset;
    Py_INCREF(this->_python_dataset);
}

std::string const &
DataSetBridge
::get_encoding() const
{
    return this->_encoding;
}

void
DataSetBridge
::set_encoding(std::string const & encoding)
{
    this->_encoding = encoding;
}

PyObject*
DataSetBridge
::to_python()
{
    PyObject* result = PyDict_New();

    for(gdcm::DataSet::ConstIterator dataset_it=this->_gdcm_dataset.Begin();
        dataset_it!=this->_gdcm_dataset.End(); ++dataset_it)
    {
        uint16_t const tag_group = dataset_it->GetTag().GetGroup();
        uint16_t const tag_element = dataset_it->GetTag().GetElement();

        if(tag_element == 0)
        {
            // Group length, do nothing
            continue;
        }
        else
        {
            PyObject* key = Py_BuildValue("(II)", tag_group, tag_element);

            if(tag_group == 0x0008 && tag_element == 0x0005)
            {
                // Specific Character Set: map to Python encoding
                gdcm::Attribute<0x0008,0x0005> attribute;
                attribute.SetFromDataElement(*dataset_it);
                // TODO : multi-valued Specific Character Set
                std::string const specific_character_set=attribute.GetValue();
                std::string dataset_encoding;
                if(specific_character_set == "") dataset_encoding = "ascii";
                // Single-byte character sets with code extensions (PS 3.3, Table C.12-2)
                else if(specific_character_set == "ISO_IR 100") dataset_encoding = "latin_1";
                else if(specific_character_set == "ISO_IR 101") dataset_encoding = "iso8859_2";
                else if(specific_character_set == "ISO_IR 109") dataset_encoding = "iso8859_3";
                else if(specific_character_set == "ISO_IR 110") dataset_encoding = "iso8859_4";
                else if(specific_character_set == "ISO_IR 144") dataset_encoding = "iso8859_5";
                else if(specific_character_set == "ISO_IR 127") dataset_encoding = "iso8859_6";
                else if(specific_character_set == "ISO_IR 126") dataset_encoding = "iso8859_7";
                else if(specific_character_set == "ISO_IR 138") dataset_encoding = "iso8859_8";
                else if(specific_character_set == "ISO_IR 148") dataset_encoding = "iso8859_9";
                else if(specific_character_set == "ISO_IR 13") dataset_encoding = "iso2022_jp";
                // CP874 seems to be a superset of TIS-620/ISO-IR-166 (e.g.
                // presence of the euro sign in the CP874 at an unassigned place
                // of TIS-620), but we should get away with it.
                else if(specific_character_set == "ISO_IR 166") dataset_encoding = "cp874";
                // Single-byte character sets with code extensions (PS 3.3, Table C.12-3)
//                ISO 2022 IR 6
//                ISO 2022 IR 100
//                ISO 2022 IR 101
//                ISO 2022 IR 109
//                ISO 2022 IR 110
//                ISO 2022 IR 144
//                ISO 2022 IR 127
//                ISO 2022 IR 126
//                ISO 2022 IR 138
//                ISO 2022 IR 148
//                ISO 2022 IR 113
//                ISO 2022 IR 166
                // Multi-byte character sets without code extensions (PS 3.3, Table C.12-4)
//                ISO 2022 IR 87
//                ISO 2022 IR 159
//                ISO 2022 IR 149
                // Multi-byte character sets without code extensions (PS 3.3, Table C.12-5)
                else if(specific_character_set == "ISO_IR 192") dataset_encoding = "utf_8";
                else if(specific_character_set == "GB18030") dataset_encoding = "gb18030";

                this->set_encoding(dataset_encoding);
            }

            PyObject* value = this->to_python(*dataset_it);

            if(value == NULL)
            {
                // An error occurred
                Py_DECREF(result);
                result = NULL;
                break;
            }

            PyDict_SetItem(result, key, value);
            Py_DECREF(key);
            Py_DECREF(value);
        }
    }

    this->set_python(result);

    return this->get_python();
}

PyObject*
DataSetBridge
::to_python(gdcm::DataElement const & data_element) const
{
    PyObject* value=NULL;

    gdcm::VR const & vr = data_element.GetVR();

    if(vr & (gdcm::VR::OB | gdcm::VR::OF | gdcm::VR::OW))
    {
        // Return a NumPy array

        gdcm::ByteValue const * byte_value = data_element.GetByteValue();

        // Make a copy, since data_element keeps ownership of its data
        char* data = new char[byte_value->GetLength()];
        char const * begin = byte_value->GetPointer();
        char const * end = byte_value->GetPointer()+byte_value->GetLength();
        std::copy(begin, end, data);

        // Determine type and dimension
        int item_type;
        npy_intp* dimensions = new npy_intp[1];
        if(vr == gdcm::VR::OB)
        {
            item_type = NPY_UBYTE;
            dimensions[0] = data_element.GetVL();
        }
        else if(vr == gdcm::VR::OF)
        {
            item_type = NPY_FLOAT32;
            dimensions[0] = data_element.GetVL()/4;
        }
        else // vr == gdcm::VR::OW
        {
            item_type = NPY_UINT16;
            dimensions[0] = data_element.GetVL()/2;
        }

        // Create array, set data ownership (since we made a copy earlier)
        PyObject* array = PyArray_SimpleNewFromData(1, dimensions, item_type, data);
        reinterpret_cast<PyArrayObject*>(array)->flags |= NPY_OWNDATA;

        delete[] dimensions;

        return array;
    }
    else if (vr == gdcm::VR::UN)
    {
        // Return str, to be used as sequence of bytes
        value = PyString_FromStringAndSize(
            data_element.GetByteValue()->GetPointer(), data_element.GetVL());
    }
    else if(vr == gdcm::VR::SQ)
    {
        gdcm::SmartPointer<gdcm::SequenceOfItems> sequence = data_element.GetValueAsSQ();

        value = PyList_New(sequence->GetNumberOfItems());

        unsigned int index=0;
        for(gdcm::SequenceOfItems::ConstIterator sequence_it=sequence->Begin();
            sequence_it!=sequence->End(); ++sequence_it, ++index)
        {
            DataSetBridge bridge(sequence_it->GetNestedDataSet());
            bridge.set_encoding(this->get_encoding());
            PyObject* python_item = bridge.to_python();
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
    else if(gdcm::VR::IsBinary(vr) || gdcm::VR::IsASCII(vr))
    {
        gdcm::ByteValue const * byte_value = data_element.GetByteValue();

        unsigned int count = 0;
        if(gdcm::VR::IsBinary(vr))
        {
            if(byte_value != NULL)
            {
                count = byte_value->GetLength()/vr.GetSize();
            }
            else
            {
                count = 0;
            }
        }
        else
        {
            if(byte_value != NULL)
            {
                count = gdcm::VM::GetNumberOfElementsFromArray(
                    byte_value->GetPointer(), byte_value->GetLength());
            }
            else
            {
                count = 0;
            }
        }

        if(count == 0)
        {
            value = Py_None;
        }
        else if(count == 1)
        {
            char const * begin = byte_value->GetPointer();
            char const * end = begin+byte_value->GetLength();
            value = this->_to_python(begin, end, vr);
        }
        else // count > 1
        {
            value = PyList_New(count);
            char const * begin = byte_value->GetPointer();
            char const * end = begin+byte_value->GetLength();
            for(unsigned int i=0; i<count; ++i)
            {
                char const * item_end = begin+this->_get_length(begin, end, vr)+1;
                PyObject* item = this->_to_python(begin, item_end, vr);
                PyList_SetItem(value, i, item);
                // Don't Py_DECREF(python_item) : according to Python doc :
                // This function “steals” a reference to item and discards a
                // reference to an item already in the list at the affected position.
                begin = item_end;
            }
        }
    }
    else
    {
        // We don't know how to handle this VR.
        value = NULL;
        PyErr_SetString(PyExc_Exception, "Unknown VR");
    }

    return value;
}

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

PyObject*
DataSetBridge
::_to_python(char const * begin, char const * end, gdcm::VR const & vr) const
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
        TRY_DECODE(object, this->_encoding.c_str());
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

        return object;
    }
    else if(vr == gdcm::VR::AT)
    {
        uint16_t e0 = *reinterpret_cast<uint16_t const *>(begin);
        uint16_t e1 = *(1+reinterpret_cast<uint16_t const *>(begin));
        return Py_BuildValue("(II)", int(e0), int(e1));
    }
    else if(vr == gdcm::VR::DS)
    {
        char* endptr;
        double const d = std::strtod(begin, &endptr);
        if(endptr == begin)
        {
            PyErr_SetString(PyExc_Exception, "Cannot parse DS");
            return NULL;
        }

        return PyFloat_FromDouble(d);
    }
    else if(vr == gdcm::VR::FD)
    {
        // TODO : correct 32 bits type
        return PyFloat_FromDouble(*reinterpret_cast<float const *>(begin));
    }
    else if(vr == gdcm::VR::FL)
    {
        // TODO : correct 64 bits type
        return PyFloat_FromDouble(*reinterpret_cast<double const *>(begin));
    }
    else if(vr == gdcm::VR::IS)
    {
        char* endptr;
        long const d = std::strtol(begin, &endptr, 10);
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
DataSetBridge
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
        return 1;
    }
    else if(vr == gdcm::VR::AT)
    {
        return 4;
    }
    else if(vr == gdcm::VR::FD)
    {
        return 4;
    }
    else if(vr == gdcm::VR::FL)
    {
        return 8;
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
