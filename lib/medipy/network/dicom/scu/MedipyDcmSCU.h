/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _8622808a_288a_469d_af20_077478e7214e
#define _8622808a_288a_469d_af20_077478e7214e

#include "DcmSCU.h"
#include <vector>

class MedipyDcmSCU : public DcmSCU
{
public:
    MedipyDcmSCU();
    virtual ~MedipyDcmSCU();
    
    void Clear();
    virtual OFCondition handleSTORERequest(
        const T_ASC_PresentationContextID presID,
        DcmDataset *incomingObject,
        OFBool &continueCGETSession,
        Uint16 &cStoreReturnStatus);

    bool GetMemoryFlag();
    void SetMemoryFlag(bool);
    
    std::vector<DcmDataset*> GetDatasets();

private:
    bool _ReceiveInMemory;
    std::vector<DcmDataset*> _datasets;

};

#endif //_8622808a_288a_469d_af20_077478e7214e
