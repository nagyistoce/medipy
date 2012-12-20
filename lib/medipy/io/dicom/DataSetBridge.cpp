/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "DataSetBridge.h"

#include <Python.h>

#include <algorithm>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <stdint.h>
#include <vector>

#include <gdcmBinEntry.h>
#include <gdcmDictSet.h>
#include <gdcmDocEntry.h>
#include <gdcmGlobal.h>
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
: _gdcm_dataset(NULL), _python_dataset(NULL)
{
    // Nothing else.
}

DataSetBridge
::DataSetBridge(gdcm::Document * dataset)
: _gdcm_dataset(NULL), _python_dataset(NULL)
{
    this->set_gdcm(dataset);
}

DataSetBridge
::DataSetBridge(PyObject* dataset)
: _gdcm_dataset(NULL), _python_dataset(NULL)
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
    this->_python_dataset = dataset;
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
    PyObject* modules = PyImport_GetModuleDict(); // Borrowed reference
    PyObject* medipy = PyMapping_GetItemString(modules, "medipy"); // New reference
    PyObject* medipy_io = PyObject_GetAttrString(medipy, "io"); // New reference
    PyObject* medipy_io_dicom = PyObject_GetAttrString(medipy_io, "dicom"); // New reference

    PyObject* Tag = PyObject_GetAttrString(medipy_io_dicom, "Tag"); // New reference

    if(Tag == NULL)
    {
        return NULL;
    }

    PyObject* DataSet = PyObject_GetAttrString(medipy_io_dicom, "DataSet"); // New reference

    if(DataSet == NULL)
    {
        return NULL;
    }

    PyObject* args = PyTuple_New(0);
    PyObject* result = PyObject_CallObject(DataSet, args);
    Py_DECREF(args);

    if(result == NULL)
    {
        return NULL;
    }

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

            // Build the tag
            PyObject* key = Py_BuildValue("(II)", entry->GetGroup(), entry->GetElement());
            PyObject* tag_args = Py_BuildValue("(O)", key);
            PyObject* tag = PyObject_CallObject(Tag, tag_args);
            Py_DECREF(tag_args);
            Py_DECREF(key);

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

            std::string vr;
            if(entry->IsImplicitVR())
            {
                gdcm::Dict* dict = gdcm::Global::GetDicts()->GetDefaultPubDict();
                gdcm::DictEntry * dictEntry = dict->GetEntry(entry->GetGroup(), entry->GetElement());
                if(dictEntry == NULL)
                {
                    vr = "UN";
                }
                else
                {
                    vr = dictEntry->GetVR();
                }

            }
            else
            {
                vr = entry->GetVR();
            }

            PyObject* MediPyVR = PyObject_GetAttrString(medipy_io_dicom, &(vr[0]));
            if(MediPyVR==NULL)
            {
                return NULL;
            }

            PyObject* args = Py_BuildValue("(O)", value);
            PyObject* typed_value = PyObject_CallObject(MediPyVR, args);

            Py_DECREF(args);
            Py_DECREF(MediPyVR);
            Py_DECREF(value);

            PyDict_SetItem(result, tag, typed_value);
            Py_DECREF(tag);
            Py_DECREF(typed_value);
        }

        entry = doc_entry_set->GetNextEntry();
    }

    Py_DECREF(Tag);
    Py_DECREF(DataSet);
    Py_DECREF(medipy_io_dicom);
    Py_DECREF(medipy_io);
    Py_DECREF(medipy);

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

void
DataSetBridge
::to_gdcm(gdcm::DocEntrySet & document)
{
    return this->_to_gdcm(this->_python_dataset, document);
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

std::vector<uint8_t> stringVR(PyObject* object, bool encodeAsUTF8, char padding)
{
    std::vector<uint8_t> result;

    if(PyList_Check(object))
    {
        int const itemsCount = PyList_Size(object);
        for(int itemIndex=0; itemIndex<itemsCount; ++itemIndex)
        {
            PyObject* pythonItem = PyList_GetItem(object, itemIndex);
            std::vector<uint8_t> const itemData =
                stringVR(pythonItem, encodeAsUTF8, padding);
            std::copy(itemData.begin(), itemData.end(),
                      std::back_inserter(result));
            if(itemIndex!=itemsCount-1)
            {
                result.push_back('\\');
            }
        }
    }
    else if(object == Py_None)
    {
        // Nothing to do : we need an empty result.
    }
    else
    {
        bool object_is_unicode = PyUnicode_Check(object);

        PyObject* string;
        if(object_is_unicode)
        {
            if(encodeAsUTF8)
            {
                string = PyUnicode_AsUTF8String(object);
            }
            else
            {
                string = PyUnicode_AsASCIIString(object);
            }
        }
        else if(PyString_Check(object))
        {
            string = object;
        }
        else
        {
            throw std::runtime_error("Object is not an string");
        }

        char* buffer = PyString_AsString(string);
        result.resize(PyString_Size(string));
        std::copy(buffer, buffer+result.size(), result.begin());

        if(object_is_unicode)
        {
            Py_DECREF(string);
        }
    }

    if(result.size()%2==1)
    {
        result.push_back(padding);
    }

    return result;
}

template<typename T>
std::vector<uint8_t> intVR(PyObject* object)
{
    std::vector<uint8_t> result;

    if(PyList_Check(object))
    {
        int const itemsCount = PyList_Size(object);
        for(int itemIndex=0; itemIndex<itemsCount; ++itemIndex)
        {
            PyObject* pythonItem = PyList_GetItem(object, itemIndex);
            std::vector<uint8_t> const itemData = intVR<T>(pythonItem);
            std::copy(itemData.begin(), itemData.end(),
                      std::back_inserter(result));
        }
    }
    else if(PyInt_Check(object))
    {
        T const value = T(PyInt_AsLong(object));
        uint8_t const * buffer = reinterpret_cast<uint8_t const *>(&value);
        result.resize(sizeof(T));
        std::copy(buffer, buffer+result.size(), result.begin());
    }
    else
    {
        throw std::runtime_error("Object is not an int");
    }

    return result;
}

template<typename T>
std::vector<uint8_t> floatVR(PyObject* object)
{
    std::vector<uint8_t> result;

    if(PyList_Check(object))
    {
        int const itemsCount = PyList_Size(object);
        for(int itemIndex=0; itemIndex<itemsCount; ++itemIndex)
        {
            PyObject* pythonItem = PyList_GetItem(object, itemIndex);
            std::vector<uint8_t> const itemData = floatVR<T>(pythonItem);
            std::copy(itemData.begin(), itemData.end(),
                      std::back_inserter(result));
        }
    }
    else if(PyFloat_Check(object))
    {
        T const value = T(PyFloat_AsDouble(object));
        uint8_t const * buffer = reinterpret_cast<uint8_t const *>(&value);
        result.resize(sizeof(T));
        std::copy(buffer, buffer+result.size(), result.begin());
    }
    else
    {
        throw std::runtime_error("Object is not an float");
    }

    return result;
}

std::vector<uint8_t> bufferVR(PyObject* object)
{
    if(!PyString_Check(object))
    {
        throw std::runtime_error("Object is not an buffer string");
    }
    std::vector<uint8_t> result(PyString_Size(object));
    std::copy(PyString_AsString(object), PyString_AsString(object)+result.size(),
              result.begin());

    if(result.size()%2==1)
    {
        result.push_back('\0');
    }

    return result;
}

std::vector<uint8_t> AE(PyObject* object) { return stringVR(object, false, ' '); }
std::vector<uint8_t> AS(PyObject* object) { return stringVR(object, false, ' '); }

std::vector<uint8_t> AT(PyObject* object)
{
    std::vector<uint8_t> result;

    if(PyList_Check(object))
    {
        int const itemsCount = PyList_Size(object);
        for(int itemIndex=0; itemIndex<itemsCount; ++itemIndex)
        {
            PyObject* pythonItem = PyList_GetItem(object, itemIndex);
            std::vector<uint8_t> itemData = AT(pythonItem);
            std::copy(itemData.begin(), itemData.end(),
                      std::back_inserter(result));
        }
    }
    else
    {
        PyObject* first = PyTuple_GetItem(object, 0);
        uint16_t first_value = uint16_t(PyInt_AsLong(first));
        PyObject* second = PyTuple_GetItem(object, 1);
        uint16_t second_value = uint16_t(PyInt_AsLong(second));

        result.resize(2*sizeof(uint16_t));
        std::copy(reinterpret_cast<uint8_t*>(&first_value),
            reinterpret_cast<uint8_t*>(&first_value)+result.size(),
            result.begin());
        std::copy(reinterpret_cast<uint8_t*>(&second_value),
            reinterpret_cast<uint8_t*>(&second_value)+result.size(),
            result.begin()+sizeof(uint16_t));
    }

    return result;
}

std::vector<uint8_t> CS(PyObject* object) { return stringVR(object, false, ' '); }
std::vector<uint8_t> DA(PyObject* object) { return stringVR(object, false, ' '); }

std::vector<uint8_t> DS(PyObject* object)
{
    std::vector<uint8_t> result;

    if(PyList_Check(object))
    {
        int const itemsCount = PyList_Size(object);
        for(int itemIndex=0; itemIndex<itemsCount; ++itemIndex)
        {
            PyObject* pythonItem = PyList_GetItem(object, itemIndex);
            std::vector<uint8_t> itemData = DS(pythonItem);
            std::copy(itemData.begin(), itemData.end(),
                      std::back_inserter(result));
            if(itemIndex!=itemsCount-1)
            {
                result.push_back('\\');
            }
        }
    }
    else
    {
        if(object == Py_None)
        {
            // Nothing to do : we need an empty result.
        }
        else if(PyUnicode_Check(object))
        {
            if(PyUnicode_GetSize(object)!=0)
            {
                throw std::runtime_error("Cannot convert to DS: "
                                         "unicode objects must be empty");
            }
            // Otherwise do nothing : we need an empty result
        }
        else if(PyString_Check(object))
        {
            if(PyString_Size(object)!=0)
            {
                throw std::runtime_error("Cannot convert to DS: "
                                         "str objects must be empty");
            }
            // Otherwise do nothing : we need an empty result
        }
        else if(PyFloat_Check(object))
        {
            double value = PyFloat_AsDouble(object);
            std::ostringstream stream;
            stream.imbue(std::locale("C"));
            stream.precision(std::numeric_limits<double>::digits10);
            stream << value;
            std::string const string = stream.str();
            result.resize(string.size());
            std::copy(string.begin(), string.end(), result.begin());
        }
        else
        {
            throw std::runtime_error("Cannot convert to DS: "
                                     "Python object is neither a float "
                                     "nor an empty string/unicode");
        }
    }

    if(result.size()%2==1)
    {
        result.push_back(' ');
    }

    return result;
}

std::vector<uint8_t> DT(PyObject* object) { return stringVR(object, false, ' '); }
// FIXME : check that double is 64 bits
std::vector<uint8_t> FD(PyObject* object) { return floatVR<double>(object); }
// FIXME : check that float is 32 bits
std::vector<uint8_t> FL(PyObject* object) { return floatVR<float>(object); }

std::vector<uint8_t> IS(PyObject* object)
{
    std::vector<uint8_t> result;

    if(PyList_Check(object))
    {
        int const itemsCount = PyList_Size(object);
        for(int itemIndex=0; itemIndex<itemsCount; ++itemIndex)
        {
            PyObject* pythonItem = PyList_GetItem(object, itemIndex);
            std::vector<uint8_t> itemData = IS(pythonItem);
            std::copy(itemData.begin(), itemData.end(),
                      std::back_inserter(result));
            if(itemIndex!=itemsCount-1)
            {
                result.push_back('\\');
            }
        }
    }
    else
    {
        if(object == Py_None)
        {
            // Nothing to do : we need an empty result.
        }
        else if(PyUnicode_Check(object))
        {
            if(PyUnicode_GetSize(object)!=0)
            {
                throw std::runtime_error("Cannot convert to IS: "
                                         "unicode objects must be empty");
            }
            // Otherwise do nothing : we need an empty result
        }
        else if(PyString_Check(object))
        {
            if(PyString_Size(object)!=0)
            {
                throw std::runtime_error("Cannot convert to IS: "
                                         "str objects must be empty");
            }
            // Otherwise do nothing : we need an empty result
        }
        else if(PyInt_Check(object))
        {
            long value = PyInt_AsLong(object);
            std::ostringstream stream;
            stream.imbue(std::locale("C"));
            stream << value;
            std::string const string = stream.str();
            result.resize(string.size());
            std::copy(string.begin(), string.end(), result.begin());
        }
        else
        {
            throw std::runtime_error("Cannot convert to IS: "
                                     "Python object is neither an int "
                                     "nor an empty string/unicode");
        }
    }

    if(result.size()%2==1)
    {
        result.push_back(' ');
    }

    return result;
}

std::vector<uint8_t> LO(PyObject* object) { return stringVR(object, true, ' '); }
// LT should not have VM > 1, but better safe than sorry.
std::vector<uint8_t> LT(PyObject* object) { return stringVR(object, true, ' '); }
std::vector<uint8_t> OB(PyObject* object) { return bufferVR(object); }
std::vector<uint8_t> OF(PyObject* object) { return bufferVR(object); }
std::vector<uint8_t> OW(PyObject* object) { return bufferVR(object); }
std::vector<uint8_t> PN(PyObject* object) { return stringVR(object, true, ' '); }
std::vector<uint8_t> SH(PyObject* object) { return stringVR(object, true, ' '); }
std::vector<uint8_t> SL(PyObject* object) { return intVR<int32_t>(object); }
std::vector<uint8_t> SS(PyObject* object) { return intVR<int16_t>(object); }
// ST should not have VM > 1, but better safe than sorry.
std::vector<uint8_t> ST(PyObject* object) { return stringVR(object, true, ' '); }
std::vector<uint8_t> TM(PyObject* object) { return stringVR(object, false, ' '); }
std::vector<uint8_t> UI(PyObject* object) { return stringVR(object, false, '\0'); }
std::vector<uint8_t> UL(PyObject* object) { return intVR<uint32_t>(object); }
std::vector<uint8_t> UN(PyObject* object) { return bufferVR(object); }
std::vector<uint8_t> US(PyObject* object) { return intVR<uint16_t>(object); }
// UT should not have VM > 1, but better safe than sorry.
std::vector<uint8_t> UT(PyObject* object) { return stringVR(object, true, ' '); }

void
DataSetBridge
::_to_gdcm(PyObject* dictionary, gdcm::DocEntrySet & document)
{
    PyObject *key, *value;
    Py_ssize_t pos = 0;

    while(PyDict_Next(dictionary, &pos, &key, &value))
    {
        long const tag = PyLong_AsLong(key);
        long const tagGroup = tag>>16;
        long const tagElement = tag&0xffff;

        PyObject* value_type = PyObject_GetAttrString(value, "__class__"); // New reference
        PyObject* vr_python = PyObject_GetAttrString(value_type, "__name__"); // New reference
        std::string const vr = PyString_AsString(vr_python);
        Py_DECREF(vr_python);
        Py_DECREF(value_type);

        if(tag == 0x00080005)
        {
            // We encode everything to UTF-8, so we overwrite the
            // Specific Character Set
            std::string const value = "ISO_IR 192";
            document.InsertValEntry((value.size()%2==0)?value:(value+" "),
                tagGroup, tagElement, vr);
            continue;
        }

        PyObject* nested_value = PyObject_GetAttrString(value, "value"); // New reference

        if(vr == "SQ")
        {
            gdcm::SeqEntry * seqEntry = document.InsertSeqEntry(
                tagGroup, tagElement);

            int const itemsCount = PyList_Size(nested_value);
            for(int itemIndex=0; itemIndex<itemsCount; ++itemIndex)
            {
                PyObject* pythonItem = PyList_GetItem(nested_value, itemIndex);
                // Depth level is only for printing, so we set it to a default
                // value
                gdcm::SQItem * gdcmItem = new gdcm::SQItem(0);
                this->_to_gdcm(pythonItem, *gdcmItem);
                seqEntry->AddSQItem(gdcmItem, itemIndex);
            }
        }
        else
        {
            this->_to_gdcm(nested_value, tagGroup, tagElement, vr, document);
        }

        Py_DECREF(nested_value);
    }
}

void
DataSetBridge
::_to_gdcm(PyObject* object,
           unsigned long group, unsigned long element, std::string vr,
           gdcm::DocEntrySet & document)
{
#define MEDIPY_VR_ACTION(vr_value) \
    if(vr == #vr_value) \
    { \
        std::vector<uint8_t> data = vr_value(object); \
        document.InsertBinEntry(&data[0], data.size(), group, element, vr); \
    }

    MEDIPY_VR_ACTION(AE)
    MEDIPY_VR_ACTION(AS)
    MEDIPY_VR_ACTION(AT)
    MEDIPY_VR_ACTION(CS)
    MEDIPY_VR_ACTION(DA)
    MEDIPY_VR_ACTION(DT)
    MEDIPY_VR_ACTION(DS)
    MEDIPY_VR_ACTION(FD)
    MEDIPY_VR_ACTION(FL)
    MEDIPY_VR_ACTION(IS)
    MEDIPY_VR_ACTION(LO)
    MEDIPY_VR_ACTION(LT)
    MEDIPY_VR_ACTION(OB)
    MEDIPY_VR_ACTION(OF)
    MEDIPY_VR_ACTION(OW)
    MEDIPY_VR_ACTION(PN)
    MEDIPY_VR_ACTION(SH)
    MEDIPY_VR_ACTION(SL)
    // SQ is not processed here !
    MEDIPY_VR_ACTION(SS)
    MEDIPY_VR_ACTION(ST)
    MEDIPY_VR_ACTION(TM)
    MEDIPY_VR_ACTION(UI)
    MEDIPY_VR_ACTION(UL)
    MEDIPY_VR_ACTION(UN)
    MEDIPY_VR_ACTION(US)
    MEDIPY_VR_ACTION(UT)

#undef MEDIPY_VR_ACTION
}
