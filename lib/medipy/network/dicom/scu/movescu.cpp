/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

//Check move RSP and beyond to set call return type

#include "movescu.h"
#include "../../../io/dicom/dataset_io/DCMTKToPython.h"
#include "../../../io/dicom/dataset_io/PythonToDCMTK.h"

#include <Python.h>
#include <stdexcept>
#include <string>

#include <dcmtk/dcmnet/assoc.h>
#include <dcmtk/dcmdata/dcvr.h>
#include <dcmtk/dcmdata/dctag.h>
#include <dcmtk/dcmdata/dcvrcs.h>
#include <dcmtk/dcmdata/dcdatset.h>
#include <dcmtk/ofstd/ofcond.h>
#include <dcmtk/dcmdata/dcfilefo.h>
#include <dcmtk/dcmnet/diutil.h>

#include <dcmtk/oflog/oflog.h>


/* ------------------ Move definition (real work) ----------------------- */

//Constructor
Move
::Move( PyObject* connection,
        PyObject* root_lvl,
        PyObject* query_lvl,
        PyObject* destination,
        PyObject* dataset)
:SCU(connection)
{
    this->_root_level = root_lvl;
    this->_query_dataset = dataset;
    this->_destination = destination;

    std::string const level = PyString_AsString(query_lvl);
    if(level == "patient") this->_query_level = PyString_FromString("PATIENT");
    else if(level == "study") this->_query_level = PyString_FromString("STUDY");
    else if(level == "series") this->_query_level = PyString_FromString("SERIES");
    else if(level == "image") this->_query_level = PyString_FromString("IMAGE");
    else
    {
        throw std::runtime_error(
            "Specified query level is none of the following : patient, study, series or image");
    }
    
    OFList<OFString> syntaxes;
    syntaxes.push_back(UID_LittleEndianImplicitTransferSyntax);
    
    //Set Presentation Context
    std::string const root = PyString_AsString(this->_root_level);
    if(root == "patient")
    {
        this->_querySyntax = { UID_FINDPatientRootQueryRetrieveInformationModel,
                        UID_MOVEPatientRootQueryRetrieveInformationModel };
    }
    else if(root == "study")
    {
        this->_querySyntax = { UID_FINDStudyRootQueryRetrieveInformationModel,
                        UID_MOVEStudyRootQueryRetrieveInformationModel };
    }
    else
    {   
        throw std::runtime_error(
            "Only 'patient' and 'study' root levels are available");
    }
    this->_scu.addPresentationContext(this->_querySyntax.findSyntax, syntaxes);
    this->_scu.addPresentationContext(this->_querySyntax.moveSyntax, syntaxes);
}

//Real work : python __call__ equivalent
PyObject*
Move
:: operator()()
{
    PyObject* parent_result = this->SCU::operator()();
    if(parent_result==NULL)
    {
        return NULL;
    }
    
    const char* sopClass = this->_querySyntax.moveSyntax;
    T_ASC_PresentationContextID presId = this->_scu.findPresentationContextID(sopClass, "");
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
    OFCondition cond = this->_scu.sendMOVERequest(
                presId,
                PyString_AsString(this->_destination),
                &dcmtk_dataset,
                &responses);
    if(cond.bad())
    {
        PyErr_SetString(this->_medipy_base_exception,cond.text());
        return NULL;
    }
    
    //Translate DCMTK DataSet into Python's
    std::vector<DcmDataset*> dcmtk_results;
    for(OFIterator<RetrieveResponse*> iter=responses.begin(); iter!=responses.end(); iter++)
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
    
    return results;
}
