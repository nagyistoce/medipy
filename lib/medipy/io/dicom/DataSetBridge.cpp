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

#include <gdcmBinEntry.h>
#include <gdcmDocEntry.h>
#include <gdcmSeqEntry.h>
#include <gdcmSQItem.h>
#include <gdcmValEntry.h>

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
::DataSetBridge(gdcm::Document * dataset)
: _python_dataset(NULL)
{
    this->set_gdcm(dataset);
}

DataSetBridge
::DataSetBridge(PyObject* dataset)
{
    this->set_python(dataset);
}

gdcm::Document *
DataSetBridge
::get_gdcm() const
{
    return this->_gdcm_dataset;
}

void
DataSetBridge
::set_gdcm(gdcm::Document * dataset)
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
    return this->to_python(this->_gdcm_dataset);
}

PyObject*
DataSetBridge
::to_python(gdcm::DocEntrySet * doc_entry_set) const
{
    PyObject* result = PyDict_New();

    gdcm::DocEntry* entry = doc_entry_set->GetFirstEntry();
    while(entry != NULL)
    {
        if(entry->GetElement() == 0)
        {
            // Group length, do nothing
        }
        else if(entry->GetGroup() == 0xFFFE && (entry->GetElement() == 0xE000 ||
                entry->GetElement() == 0xE00D || entry->GetElement() == 0xE0DD))
         {
             // Item, Item Delimitation Item and Sequence Delimitation Item
             // Skip them
         }
        else
        {
            if(entry->GetGroup() == 0x0008 && entry->GetElement() == 0x0005)
            {
                gdcm::ValEntry* charset(dynamic_cast<gdcm::ValEntry*>(entry));
                // Single-byte character sets without code extensions (Table C.12-2)
                if(charset->GetValue() == "") this->_encoding = "ascii";
                else if(charset->GetValue() == "ISO_IR 100") this->_encoding = "latin_1";
                else if(charset->GetValue() == "ISO_IR 101") this->_encoding = "iso8859_2";
                else if(charset->GetValue() == "ISO_IR 109") this->_encoding = "iso8859_3";
                else if(charset->GetValue() == "ISO_IR 110") this->_encoding = "iso8859_4";
                else if(charset->GetValue() == "ISO_IR 144") this->_encoding = "iso8859_5";
                else if(charset->GetValue() == "ISO_IR 127") this->_encoding = "iso8859_6";
                else if(charset->GetValue() == "ISO_IR 126") this->_encoding = "iso8859_7";
                else if(charset->GetValue() == "ISO_IR 138") this->_encoding = "iso8859_8";
                else if(charset->GetValue() == "ISO_IR 148") this->_encoding = "iso8859_9";
                else if(charset->GetValue() == "ISO_IR 13") this->_encoding = "iso2022_jp";
                // ISO_IR 166
                // Single-byte character sets with code extensions (Table C.12-3)
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
                // Multi-byte character sets without code extensions (Table C.12-4)
//                ISO 2022 IR 87
//                ISO 2022 IR 159
//                ISO 2022 IR 149
            // Multi-byte character sets without code extensions (Table C.12-5)
                else if(charset->GetValue() == "ISO_IR 192") this->_encoding = "utf_8";
                else if(charset->GetValue() == "GB18030") this->_encoding = "gb18030";
            }

            PyObject* key = Py_BuildValue("(II)", entry->GetGroup(), entry->GetElement());

            PyObject* value = NULL;
            if(dynamic_cast<gdcm::BinEntry*>(entry))
            {
                value = this->to_python(dynamic_cast<gdcm::BinEntry*>(entry));
            }
            else if(dynamic_cast<gdcm::ValEntry*>(entry))
            {
                value = this->to_python(dynamic_cast<gdcm::ValEntry*>(entry));
            }
            else if(dynamic_cast<gdcm::SeqEntry*>(entry))
            {
                gdcm::SeqEntry & sequence = *dynamic_cast<gdcm::SeqEntry*>(entry);
                gdcm::SQItem* item = sequence.GetFirstSQItem();

                value = PyList_New(0);
                while(item != NULL)
                {
                    PyObject* python_item = this->to_python(item);
                    if(python_item == NULL)
                    {
                        return NULL;
                    }
                    PyList_Append(value, python_item);
                    Py_DECREF(python_item);
                    item = sequence.GetNextSQItem();
                }
            }
            else
            {
                PyErr_SetString(PyExc_Exception, "Unknown entry type");
                return NULL;
            }

            if(value == NULL)
            {
                // An error occurred
                return NULL;
            }

            PyDict_SetItem(result, key, value);
            Py_DECREF(key);
            Py_DECREF(value);
        }

        entry = doc_entry_set->GetNextEntry();
    }

    return result;
}

PyObject*
DataSetBridge
::to_python(gdcm::BinEntry * data_element) const
{
    char* data = reinterpret_cast<char*>(data_element->GetBinArea());
    if(data_element->GetVR() == "AT")
    {
        unsigned int const vm = data_element->GetLength()/4;
        PyObject* at = NULL;
        if(vm == 1)
        {
            uint16_t e0 = *reinterpret_cast<uint16_t*>(data_element->GetBinArea());
            uint16_t e1 = *reinterpret_cast<uint16_t*>(data_element->GetBinArea()+2);
            at = Py_BuildValue("(II)", int(e0), int(e1));
        }
        else
        {
            at = PyList_New(0);
            for(unsigned int i=0; i<vm; ++i)
            {
                uint16_t e0 = *reinterpret_cast<uint16_t*>(data_element->GetBinArea()+4*i);
                uint16_t e1 = *reinterpret_cast<uint16_t*>(data_element->GetBinArea()+4*i+2);
                PyObject* item = Py_BuildValue("(II)", int(e0), int(e1));
                PyList_Append(at, item);
                Py_DECREF(item);
            }
        }
        return at;
    }
    else
    {
        // Return a string of bytes
        return PyString_FromStringAndSize(data, data_element->GetLength());
    }
}

PyObject*
DataSetBridge
::to_python(gdcm::ValEntry * data_element) const
{
    std::string const & value = data_element->GetValue();
    std::string const & vr = data_element->GetVR();

    if(vr != "LT" && vr != "ST" && vr != "UT" &&
       value.find('\\') != std::string::npos)
    {
        return this->_parse_multi_valued(value, vr);
    }
    else
    {
        // LT, ST and UT may not be multi-valued
        return this->_parse_single_valued(value, vr);
    }
}

PyObject*
DataSetBridge
::_parse_single_valued(std::string const & value, std::string const & vr) const
{
    // Explicit size is necessary due to the presence of NUL.
    // NUL is used to pad UI
    static std::string const whitespace(" \0", 2);

    if(vr == "AE" || vr == "AS" || vr == "CS" || vr == "DA" || vr == "DT" ||
       vr == "LO" || vr == "LT" || vr == "PN" || vr == "SH" || vr == "ST" ||
       vr == "TM" || vr == "UI" || vr == "UT")
    {
        std::string::size_type first=0;
        std::string::size_type size=0;

        std::string::size_type const last = value.find_last_not_of(whitespace);
        if(last != std::string::npos)
        {
            if(vr != "LT" && vr != "ST" && vr != "UT")
            {
                // Leading spaces are significant for LT, ST, and UT
                first = value.find_first_not_of(whitespace);
                if(first == std::string::npos)
                {
                    first = 0;
                }
            }
            size = (last-first)+1;
        }

#define TRY_DECODE(o, e) \
    if(object==NULL) \
    { \
        /* Avoid propagating the error. */ \
        PyErr_Clear(); \
        /*std::cout << "Trying " << e << " for " << value.c_str()+first << " ";*/ \
        object = PyUnicode_Decode(value.c_str()+first, size, e, "strict"); \
        /*std::cout << object << std::endl;*/ \
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

        return object;
    }
    else if(vr == "DS" || vr == "FD" || vr == "FL")
    {
        std::istringstream in(value);
        in.imbue(std::locale("C"));
        double d;
        in >> d;
        return PyFloat_FromDouble(d);
    }
    else if(vr == "IS" || vr == "SL" || vr == "SS" || vr == "UL" || vr == "US")
    {
        std::istringstream in(value);
        in.imbue(std::locale("C"));
        long l;
        in >> l;
        return PyInt_FromLong(l);
    }
    else
    {
        PyErr_SetString(PyExc_Exception, ("Cannot parse VR "+vr).c_str());
        return NULL;
    }
}

PyObject*
DataSetBridge
::_parse_multi_valued(std::string const & value, std::string const & vr) const
{
    PyObject* list = PyList_New(0);

    std::string::size_type start = 0;
    while(start != std::string::npos)
    {
        // Split on '\\', parse each element
        std::string::size_type stop = value.find('\\', start);
        std::string::size_type size = (stop==std::string::npos)?std::string::npos:(stop-start);
        PyObject* item = this->_parse_single_valued(value.substr(start, size), vr);
        if(item == NULL)
        {
            // An error occured
            return NULL;
        }
        PyList_Append(list, item);
        Py_DECREF(item);
        if(stop != std::string::npos)
        {
            start = stop+1;
        }
        else
        {
            start = stop;
        }
    }

    return list;
}

