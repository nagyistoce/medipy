/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "dataset_io.h"

#include <sstream>
#include <stdexcept>
#include <string>

#include <dcmtk/config/osconfig.h>
#include <dcmtk/dcmdata/dctk.h>
#include <dcmtk/ofstd/ofcond.h>

#include <Python.h>

#include "DCMTKToPython.h"
#include "PythonToDCMTK.h"

PyObject* read(std::string const & filename)
{
    OFCondition condition;
    
    DcmFileFormat reader;
    condition = reader.loadFile(filename.c_str());
    if(condition.bad())
    {
        std::ostringstream message;
        message << "Cannot read '" << filename << "': " << condition.text();
        throw std::runtime_error(message.str());
    }
    
    DCMTKToPython converter;
    return converter(reader.getDataset());
}

void write(PyObject* dataset, std::string const & filename)
{
    PythonToDCMTK converter;
    DcmDataset dcmtk_dataset = converter(dataset);
    {
        DcmFileFormat writer(&dcmtk_dataset);
        OFCondition const condition = writer.saveFile(filename.c_str());
        if(condition.bad())
        {
            std::ostringstream message;
            message << "Cannot write '" << filename << "': " << condition.text();
            throw std::runtime_error(message.str());
        }
    }
}
