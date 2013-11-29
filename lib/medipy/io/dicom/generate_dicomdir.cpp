/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "generate_dicomdir.h"

#include <string>
#include <vector>

#include <Python.h>

#include <dcmtk/config/osconfig.h>
#include <dcmtk/dcmdata/dctk.h>
#include <dcmtk/oflog/oflog.h>

#include "DicomDirGenerator.h"

void 
generate_dicomdir_cpp(PyObject * files, PyObject * root, 
    PyObject * dicomdir, PyObject * patient_extra_attributes, 
    PyObject * study_extra_attributes, PyObject * series_extra_attributes)
{
    /*********************************************************
     * Caution : no error checking is made in this function. *
     *********************************************************/
    
    /***********************
     * Argument conversion *
     ***********************/
     
    // files
    std::vector<std::string> files_cpp;
    for(Py_ssize_t i=0; i<PySequence_Length(files); ++i)
    {
        PyObject* it = PySequence_GetItem(files, i);
        files_cpp.push_back(PyString_AsString(it));
        Py_DECREF(it);
    }
    
    // root
    std::string const root_cpp(PyString_AsString(root));
    
    // dicomdir
    std::string const dicomdir_cpp(PyString_AsString(dicomdir));
    
    // patient_extra_attributes
    std::vector<DcmTagKey> const patient_extra_attributes_cpp(
        _convert_attributes_list(patient_extra_attributes));
    
    // study_extra_attributes
    std::vector<DcmTagKey> const study_extra_attributes_cpp(
        _convert_attributes_list(study_extra_attributes));
    
    // series_extra_attributes
    std::vector<DcmTagKey> const series_extra_attributes_cpp(
        _convert_attributes_list(series_extra_attributes));
    
    /************************************************
     * Create and configure the DICOMDIR generator. *
     ************************************************/
    
    OFLog::configure(OFLogger::WARN_LOG_LEVEL);
    DicomDirGenerator generator;
    generator.enableMapFilenamesMode();
    generator.createNewDicomDir(DicomDirGenerator::AP_GeneralPurpose, 
        PyString_AsString(dicomdir));
    generator.setPatientExtraAttributes(patient_extra_attributes_cpp);
    generator.setStudyExtraAttributes(study_extra_attributes_cpp);
    generator.setSeriesExtraAttributes(series_extra_attributes_cpp);
    
    /**********************************
     *  Add each file to the DICOMDIR *
     **********************************/
    for(std::vector<std::string>::const_iterator it=files_cpp.begin(); 
        it!=files_cpp.end(); ++it)
    {
        generator.addDicomFile(it->c_str());
    }
    
    // Write the DICOMDIR file.
    generator.writeDicomDir();
}

std::vector<DcmTagKey> _convert_attributes_list(PyObject * list)
{
    std::vector<DcmTagKey> attributes;
    for(Py_ssize_t i=0; i<PySequence_Length(list); ++i)
    {
        PyObject* it = PySequence_GetItem(list, i);
        long const tag = PyInt_AsLong(it);
        attributes.push_back(DcmTagKey(tag>>16, tag&0xffff));
        Py_DECREF(it);
    }
    
    return attributes;
}
