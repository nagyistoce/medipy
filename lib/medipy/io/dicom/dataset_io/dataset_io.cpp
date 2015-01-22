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
#include <dcmtk/dcmdata/dcistrmb.h>
#include <dcmtk/ofstd/ofcond.h>

#include <gdcmAttribute.h>
#include <gdcmUIDGenerator.h>
#include <gdcmReader.h>
#include <gdcmWriter.h>

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

PyObject* get_medipy_base_exception()
{
    PyObject * base = PyImport_ImportModule("medipy.base"); // New reference
    PyObject * exception = PyObject_GetAttrString(base, "Exception"); // New reference
    Py_DECREF(base);
    return exception;
}

PyObject* read(std::string const & filename)
{
    // GDCM is more robust when reading weird files. Use it to "convert" the
    // file so that DCMTK can read it.
    
    static PyObject * medipy_base_exception = get_medipy_base_exception();
    
    gdcm::Trace::DebugOn();
    gdcm::Trace::WarningOn();
    gdcm::Trace::ErrorOn();
    
    gdcm::Reader reader;
    reader.SetFileName(filename.c_str());
    bool const gdcm_read_ok = reader.Read();
    if(!gdcm_read_ok)
    {
        PyErr_SetString(medipy_base_exception, "Could not read file");
        return NULL;
    }
    
    gdcm::DataSet & gdcm_dataset = reader.GetFile().GetDataSet();
    
    // If SOP Instance UID is not in the dataset, then the Writer fails
    // silently
    if(!gdcm_dataset.FindDataElement(gdcm::Tag(0x0008,0x0018)))
    {
        gdcm::UIDGenerator generator;
        gdcm::Attribute<0x0008,0x0018> attribute = {generator.Generate()};
        gdcm_dataset.Insert(attribute.GetAsDataElement());
    }
    
    gdcm::Writer writer;
    std::ostringstream gdcm_stream;
    writer.SetFile(reader.GetFile());
    writer.SetStream(gdcm_stream);
    bool const write_ok = writer.Write();
    if(!write_ok)
    {
        PyErr_SetString(medipy_base_exception, "Could not write to stream");
        return NULL;
    }
    
    std::string const data = gdcm_stream.str();
    
    // Read the data with DCMTK
    DcmInputBufferStream dcmtk_stream;
    dcmtk_stream.setBuffer(&data[0], data.size());
    dcmtk_stream.setEos();
    
    DcmFileFormat format;
    format.transferInit();
    OFCondition const dcmtk_read_ok = format.read(dcmtk_stream);
    format.transferEnd();
    dcmtk_stream.releaseBuffer();
    
    if(!dcmtk_read_ok.good())
    {
        PyErr_SetString(medipy_base_exception, "Could not read from stream");
        return NULL;
    }
    
    format.loadAllDataIntoMemory();
    
    // Convert header and dataset
    DCMTKToPython converter;
    PyObject * header = converter(format.getMetaInfo());
    PyObject * dataset = converter(format.getDataset());
    
    PyObject_SetAttrString(dataset, "header", header);
    Py_DECREF(header);
    
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
