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

bool can_read(std::string const & filename)
{
    bool result = true;
    
    FILE* file = fopen(filename.c_str(), "rb");
    if(file==NULL)
    {
        result = false;
    }
    else
    {
        char signature[4];
        if(fseek(file, DCM_PreambleLen, SEEK_SET) < 0)
        {
            result = false;
        }
        else if(fread(signature, 1, DCM_MagicLen, file) != DCM_MagicLen ||
                strncmp(signature, DCM_Magic, DCM_MagicLen) != 0)
        {
            result = false;
        }
        fclose(file);
    }
    
    return result;
}

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
    PyObject * dataset = converter(reader.getDataset());
    PyObject_SetAttrString(dataset, "header", converter(reader.getMetaInfo()));
    
    return dataset;
}

void write(PyObject* dataset, std::string const & filename)
{
    PythonToDCMTK converter;
    DcmDataset dcmtk_dataset = converter(dataset);
    DcmFileFormat writer(&dcmtk_dataset);
    
    // Update the writer meta-info with the dataset header, if present
    if(PyObject_HasAttrString(dataset, "header"))
    {
        // New reference
        PyObject * python_header = PyObject_GetAttrString(dataset, "header"); 
        DcmDataset dcmtk_header = converter(python_header);
        Py_XDECREF(python_header);
        
        DcmObject * it = NULL;
        while(NULL != (it = dcmtk_header.nextInContainer(it)))
        {
            // The DcmMetaInfo takes ownership of the inserted element, hence
            // the hoop-jumping with clone/dynamic_cast
            DcmElement * element = dynamic_cast<DcmElement*>(it->clone());
            if(element == NULL)
            {
                throw std::runtime_error("Cannot insert element");
            }
            writer.getMetaInfo()->insert(element);
        }
    }

    // If the transfer syntax is missing, default to Little Endian Explicit VR
    if(!writer.getMetaInfo()->tagExists(DCM_TransferSyntaxUID))
    {
        writer.getMetaInfo()->putAndInsertOFStringArray(DCM_TransferSyntaxUID, 
            DcmXfer(EXS_LittleEndianExplicit).getXferID());
    }
    
    OFString value;
    writer.getMetaInfo()->findAndGetOFStringArray(DCM_TransferSyntaxUID, value);
    DcmXfer const transfer_syntax(value.c_str());
    
    // Transfer syntax /must/ be specified if it is not present in the dataset.
    OFCondition const condition = writer.saveFile(
        filename.c_str(), transfer_syntax.getXfer());
    if(condition.bad())
    {
        std::ostringstream message;
        message << "Cannot write '" << filename.c_str() << "': " << condition.text();
        throw std::runtime_error(message.str());
    }
}
