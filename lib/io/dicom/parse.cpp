/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "parse.h"
#include <locale>
#include <sstream>
#include <string>

#include <gdcmBinEntry.h>
#include <gdcmContentEntry.h>
#include <gdcmDocEntry.h>
#include <gdcmDocEntrySet.h>
#include <gdcmFile.h>
#include <gdcmSeqEntry.h>
#include <gdcmSQItem.h>
#include <gdcmValEntry.h>
#include <Python.h>

/**
 * @brief Parse a DICOM element with VM=1.
 * @return A new reference
 */
PyObject* parse_single_valued(std::string const & value, std::string const & vr,
                              std::string const & encoding = "ascii",
                              std::string const & errors = "strict")
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
        object = PyUnicode_Decode(value.c_str()+first, size, e, errors.c_str()); \
        /*std::cout << object << std::endl;*/ \
    }

        PyObject* object = NULL;
        TRY_DECODE(object, encoding.c_str());
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

/**
 * @brief Parse a DICOM element with VM>1
 * @return A new reference
 */
PyObject* parse_multi_valued(std::string const & value, std::string const & vr,
                             std::string const & encoding = "ascii",
                             std::string const & errors = "strict")
{
    PyObject* list = PyList_New(0);

    std::string::size_type start = 0;
    while(start != std::string::npos)
    {
        // Split on '\\', parse each element
        std::string::size_type stop = value.find('\\', start);
        std::string::size_type size = (stop==std::string::npos)?std::string::npos:(stop-start);
        PyObject* item = parse_single_valued(value.substr(start, size), vr,
                                             encoding, errors);
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

/**
 * @brief Parse a DICOM element with a string representation (i.e. gdcm::ValEntry).
 * @return A new reference
 *
 * According to gdcm::VR::IsVROfStringRepresentable, these are (21 VRs) :
 *   * AE
 *   * AS
 *   * CS
 *   * DA
 *   * DS
 *   * DT
 *   * FD
 *   * FL
 *   * IS
 *   * LO
 *   * LT
 *   * PN
 *   * SH
 *   * SL
 *   * SS
 *   * ST
 *   * TM
 *   * UI
 *   * UL
 *   * US
 *   * UT
 *
 * The 6 remaining VRs are 5 VRs without a string representation and SQ.
 */
PyObject* parse_val_entry(gdcm::ValEntry & entry,
                          std::string const & encoding = "ascii",
                          std::string const & errors = "strict")
{
    std::string const & value = entry.GetValue();
    std::string const & vr = entry.GetVR();

    if(vr != "LT" && vr != "ST" && vr != "UT" &&
       value.find('\\') != std::string::npos)
    {
        return parse_multi_valued(value, vr, encoding, errors);
    }
    else
    {
        // LT, ST and UT may not be multi-valued
        return parse_single_valued(value, vr, encoding, errors);
    }
}

/**
 * @brief Parse a DICOM element without a string representation (i.e. gdcm::BinEntry).
 * @return A new reference
 *
 * According to gdcm::VR::IsVROfStringRepresentable, these are (5 VRs):
 *   * AT
 *   * OB
 *   * OF
 *   * OW
 *   * UN
 *
 * The 22 remaining VRs are 21 VRs with a string representation and SQ.
 */
PyObject* parse_bin_entry(gdcm::BinEntry & entry)
{
    char* data = reinterpret_cast<char*>(entry.GetBinArea());
    if(entry.GetVR() == "AT")
    {
        unsigned int const vm = entry.GetLength()/4;
        PyObject* at = NULL;
        if(vm == 1)
        {
            uint16_t e0 = *reinterpret_cast<uint16_t*>(entry.GetBinArea());
            uint16_t e1 = *reinterpret_cast<uint16_t*>(entry.GetBinArea()+2);
            at = Py_BuildValue("(II)", int(e0), int(e1));
        }
        else
        {
            at = PyList_New(0);
            for(unsigned int i=0; i<vm; ++i)
            {
                uint16_t e0 = *reinterpret_cast<uint16_t*>(entry.GetBinArea()+4*i);
                uint16_t e1 = *reinterpret_cast<uint16_t*>(entry.GetBinArea()+4*i+2);
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
        return PyString_FromStringAndSize(data, entry.GetLength());
    }
}

/**
 * @brief Parse a container of DICOM elements.
 * @return A new reference
 */
PyObject* parse_doc_entry_set(gdcm::DocEntrySet & doc_entry_set,
                              std::string const & encoding = "ascii",
                              std::string const & errors = "strict")
{
    PyObject* result = PyDict_New();

    std::string dataset_encoding = encoding;

    gdcm::DocEntry* entry = doc_entry_set.GetFirstEntry();
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
                if(charset->GetValue() == "") dataset_encoding = "ascii";
                else if(charset->GetValue() == "ISO_IR 100") dataset_encoding = "latin_1";
                else if(charset->GetValue() == "ISO_IR 101") dataset_encoding = "iso8859_2";
                else if(charset->GetValue() == "ISO_IR 109") dataset_encoding = "iso8859_3";
                else if(charset->GetValue() == "ISO_IR 110") dataset_encoding = "iso8859_4";
                else if(charset->GetValue() == "ISO_IR 144") dataset_encoding = "iso8859_5";
                else if(charset->GetValue() == "ISO_IR 127") dataset_encoding = "iso8859_6";
                else if(charset->GetValue() == "ISO_IR 126") dataset_encoding = "iso8859_7";
                else if(charset->GetValue() == "ISO_IR 138") dataset_encoding = "iso8859_8";
                else if(charset->GetValue() == "ISO_IR 148") dataset_encoding = "iso8859_9";
                else if(charset->GetValue() == "ISO_IR 13") dataset_encoding = "iso2022_jp";
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
                else if(charset->GetValue() == "ISO_IR 192") dataset_encoding = "utf_8";
                else if(charset->GetValue() == "GB18030") dataset_encoding = "gb18030";
            }

            PyObject* key = Py_BuildValue("(II)", entry->GetGroup(), entry->GetElement());

            PyObject* value = NULL;
            if(dynamic_cast<gdcm::BinEntry*>(entry))
            {
                value = parse_bin_entry(*dynamic_cast<gdcm::BinEntry*>(entry));
            }
            else if(dynamic_cast<gdcm::ValEntry*>(entry))
            {
                value = parse_val_entry(*dynamic_cast<gdcm::ValEntry*>(entry),
                                        dataset_encoding, errors);
            }
            else if(dynamic_cast<gdcm::SeqEntry*>(entry))
            {
                gdcm::SeqEntry & sequence = *dynamic_cast<gdcm::SeqEntry*>(entry);
                gdcm::SQItem* item = sequence.GetFirstSQItem();

                value = PyList_New(0);
                while(item != NULL)
                {
                    PyObject* python_item = parse_doc_entry_set(*item, dataset_encoding, errors);
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

        entry = doc_entry_set.GetNextEntry();
    }

    return result;
}

PyObject* parse_file(std::string const & filename)
{
    gdcm::Document document;
    document.SetFileName(filename);
    document.Load();

    if(!document.IsReadable())
    {
        PyErr_SetString(PyExc_Exception, ("Cannot parse "+filename).c_str());
        return NULL;
    }

    // Load /all/ binary entries
    gdcm::DocEntry* entry = document.GetFirstEntry();
    while(entry != NULL)
    {
        gdcm::BinEntry* bin_entry = dynamic_cast<gdcm::BinEntry*>(entry);
        if(bin_entry != NULL)
        {
            document.LoadEntryBinArea(bin_entry);
        }
        entry = document.GetNextEntry();
    }

    return parse_doc_entry_set(document);
}
