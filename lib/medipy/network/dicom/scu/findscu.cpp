/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "findscu.h"
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
Find
::Find(PyObject* connection, PyObject* root_lvl, PyObject* query_lvl, PyObject* dataset) throw(SCUException)
:SCU(connection)
{
    std::string const root_level = PyString_AsString(root_lvl);
    this->_query_dataset = dataset;
    
    std::string const level = PyString_AsString(query_lvl);
    if(level == "patient") this->_query_level = "PATIENT";
    else if(level == "study") this->_query_level = "STUDY";
    else if(level == "series") this->_query_level = "SERIES";
    else if(level == "image") this->_query_level = "IMAGE";
    else
    {
        throw SCUException(
            "Specified query level is none of the following : patient, study, series or image");
    }
    
    OFList<OFString> syntaxes;
    syntaxes.push_back(UID_LittleEndianImplicitTransferSyntax);
    
    //Set abstract syntax associated to root_level
    if(root_level == "patient")
    {
        this->_abstractSyntax = UID_FINDPatientRootQueryRetrieveInformationModel;
    }
    else if(root_level == "study")
    {
        this->_abstractSyntax = UID_FINDStudyRootQueryRetrieveInformationModel;
    }
    else
    {   
        throw SCUException(
            "Only 'patient' and 'study' root levels are available");
    }
    
    this->_scu.addPresentationContext(this->_abstractSyntax,syntaxes);
}

//Real work : python __call__ equivalent
PyObject*
Find
::operator()()
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
    DcmCodeString* element = new DcmCodeString(DCM_QueryRetrieveLevel);
    element->putString(this->_query_level.c_str());
    
    //Put the override key into dset replacing existing tags
    dcmtk_dataset.insert(element, OFTrue);
        
    OFList<QRResponse*> responses;
    OFCondition cond = this->_scu.sendFINDRequest(presId, &dcmtk_dataset, &responses);
    if(cond.bad())
    {
        PyErr_SetString(this->_medipy_base_exception, cond.text());
        return NULL;
    }
    
    //Translate DCMTK DataSet into Python's
    std::vector<DcmDataset*> dcmtk_results;
    for(OFIterator<QRResponse*> iter=responses.begin(); iter!=responses.end(); iter++)
    {
        if((*iter)->m_dataset!=NULL)    dcmtk_results.push_back((*iter)->m_dataset);
    }
    DCMTKToPython py_converter;
    PyObject* results = PyList_New(dcmtk_results.size());//New reference
    
    for(unsigned int index=0; index<dcmtk_results.size(); index++)
    {
        DcmDataset* dataset = dcmtk_results[index];
        PyList_SetItem(results, index, py_converter(dataset));
    }
    
    // Clean-up the responses
    for(OFIterator<QRResponse*> iter=responses.begin(); iter!=responses.end(); iter++)
    {
        delete *iter;
    }
    
    return results;
}
