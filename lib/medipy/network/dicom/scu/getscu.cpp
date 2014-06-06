/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "getscu.h"

#include "../../../io/dicom/dataset_io/DCMTKToPython.h"
#include "../../../io/dicom/dataset_io/PythonToDCMTK.h"
#include "DcmSCU.h"

#include <stdexcept>
#include <vector>

#include <dcmtk/dcmdata/dcvr.h>
#include <dcmtk/dcmdata/dctag.h>
#include <dcmtk/dcmdata/dcvrcs.h>
#include <dcmtk/dcmdata/dcdatset.h>


//Constructor
Get
:: Get(PyObject* connection, PyObject* root_lvl, PyObject* query_lvl, PyObject* dataset) throw(SCUException)
:SCU(connection)
{
    this->_root_level = root_lvl;
    this->_query_dataset = dataset;
    
    std::string const level = PyString_AsString(query_lvl);
    if(level == "patient") this->_query_level = PyString_FromString("PATIENT");
    else if(level == "study") this->_query_level = PyString_FromString("STUDY");
    else if(level == "series") this->_query_level = PyString_FromString("SERIES");
    else if(level == "image") this->_query_level = PyString_FromString("IMAGE");
    else
    {
        throw SCUException(
            "Specified query level is none of the following : patient, study, series or image");
    }
    
    this->_scu.SetMemoryFlag(true);
    
    OFList<OFString> syntaxes;
    syntaxes.push_back(UID_LittleEndianImplicitTransferSyntax);
    
    //Set abstract syntax associated to root_level
    std::string const root = PyString_AsString(this->_root_level);
    if(root == "patient")
    {
        this->_abstractSyntax = UID_GETPatientRootQueryRetrieveInformationModel;
    }
    else if(root == "study")
    {
        this->_abstractSyntax = UID_GETStudyRootQueryRetrieveInformationModel;
    }
    else
    {   
        throw SCUException(
            "Only 'patient' and 'study' root levels are available");
    }
    
    this->_scu.addPresentationContext(this->_abstractSyntax,syntaxes);
    
    // Add storage presentation contexts (long list of storage SOP classes, uncompressed)
    for (Uint16 j = 0; j < numberOfDcmLongSCUStorageSOPClassUIDs; j++)
    {
        this->_scu.addPresentationContext(dcmLongSCUStorageSOPClassUIDs[j], syntaxes, ASC_SC_ROLE_SCP);
    }
}

//Real work : python __call__ equivalent
PyObject*
Get
:: operator ()()
{
    PyObject* parent_result = this->SCU::operator()();
    if(parent_result==NULL)
    {
        return NULL;
    }

    T_ASC_PresentationContextID presId = this->_scu.findPresentationContextID(this->_abstractSyntax, "");
    if (presId == 0)
    {
        PyErr_SetString(this->_medipy_base_exception, "No adequate Presentation Contexts");
        return NULL;

    }

    //Translate Python DataSet into DCMTK's
    PythonToDCMTK converter;
    DcmDataset dcmtk_dataset = converter(this->_query_dataset);
    
    //Set query level option (-k # overridekeys)
    const char* value = PyString_AsString(this->_query_level);
    
    //Create DcmElement
    DcmTag tag(0x0008, 0x0052);
    tag.setVR(EVR_CS);
    DcmCodeString* element = new DcmCodeString(tag);
    element->putString(value);
    
    //Put the override key into dset replacing existing tags
    dcmtk_dataset.insert(element, OFTrue);
        
    OFList<RetrieveResponse*> responses;
    OFCondition cond = this->_scu.sendCGETRequest(presId, &dcmtk_dataset, &responses);
    if(cond.bad())
    {
        PyErr_SetString(this->_medipy_base_exception, cond.text());
        return NULL;
    }
    
    //Translate DCMTK DataSet into Python's
    std::vector<DcmDataset*> received_datasets = this->_scu.GetDatasets();
    
    DCMTKToPython py_converter;
    PyObject* results = PyList_New(received_datasets.size());//New reference
    
    for(unsigned int index=0; index<received_datasets.size(); index++)
    {
        DcmDataset* dataset = received_datasets[index];
        PyList_SetItem(results, index, py_converter(dataset));
    }
    
    return results;
}
