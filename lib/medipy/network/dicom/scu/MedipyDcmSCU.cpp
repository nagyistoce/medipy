/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/
 
#include "MedipyDcmSCU.h"

MedipyDcmSCU
::MedipyDcmSCU()
:DcmSCU()
{
    this->_ReceiveInMemory = false;
}

MedipyDcmSCU
::~MedipyDcmSCU()
{
    this->Clear();
}

void
MedipyDcmSCU
::Clear()
{
    for(unsigned int i=0; i<this->_datasets.size(); i++)
    {
        delete this->_datasets[i];
    }
}

OFCondition
MedipyDcmSCU
::handleSTORERequest(
        const T_ASC_PresentationContextID presID,
        DcmDataset *incomingObject,
        OFBool &continueCGETSession,
        Uint16 &cStoreReturnStatus)
{
    OFCondition result = EC_Normal;
    if(this->_ReceiveInMemory)
    {
        DcmDataset* new_dataset = reinterpret_cast<DcmDataset*>(incomingObject->clone());
        this->_datasets.push_back(new_dataset);
    }
    else
    {
        result = DcmSCU::handleSTORERequest(presID, incomingObject, continueCGETSession, cStoreReturnStatus);
    }
    return result;
}

bool
MedipyDcmSCU
::GetMemoryFlag()
{
    return this->_ReceiveInMemory;
}

void
MedipyDcmSCU
::SetMemoryFlag(bool value)
{
    this->_ReceiveInMemory = value;
}
    
std::vector<DcmDataset*>
MedipyDcmSCU
::GetDatasets()
{
    return this->_datasets;
}
