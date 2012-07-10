/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "parse.h"

#include <Python.h>

#include <gdcmFile.h>
#include <gdcmFileMetaInformation.h>
#include <gdcmReader.h>

#include "DataSetBridge.h"

PyObject* parse_file(std::string const & filename)
{
    gdcm::Reader reader;
    reader.SetFileName(filename.c_str());
    reader.Read();

    gdcm::File const & file  = reader.GetFile();

    gdcm::FileMetaInformation const & header = file.GetHeader();
    DataSetBridge header_bridge(header);
    PyObject* header_dict = header_bridge.to_python();

    gdcm::DataSet const & data_set = file.GetDataSet();
    DataSetBridge data_set_bridge(data_set);
    PyObject* data_set_dict = data_set_bridge.to_python();

    int const merged = PyDict_Update(data_set_dict, header_dict);

    if(merged == 0)
    {
        Py_DECREF(header_dict);
        return data_set_dict;
    }
    else
    {
        Py_DECREF(header_dict);
        Py_DECREF(data_set_dict);
        return NULL;
    }
}
