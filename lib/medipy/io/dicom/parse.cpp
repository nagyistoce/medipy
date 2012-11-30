/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "parse.h"

#include <string>

#include <gdcmBinEntry.h>
#include <gdcmDocEntry.h>
#include <gdcmDocument.h>
#include <Python.h>

#include "DataSetBridge.h"

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

    DataSetBridge data_set_bridge(&document);
    return data_set_bridge.to_python();
}
