/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "echoscu.h"
#include "DcmSCU.h"

#include <Python.h>

#include <dcmtk/config/osconfig.h>
#include <dcmtk/dcmdata/dcdatset.h>
#include <dcmtk/dcmnet/assoc.h>
#include <dcmtk/dcmnet/dimse.h>
#include <dcmtk/dcmnet/diutil.h>

//Constructor
Echo
::Echo(PyObject* connection)
:SCU(connection)
{
    //Set Specific Presentation Context
    OFList<OFString> syntaxes;
    syntaxes.push_back(UID_LittleEndianImplicitTransferSyntax);
    
    this->_abstractSyntax = UID_VerificationSOPClass;
    
    this->_scu.addPresentationContext(this->_abstractSyntax,syntaxes);
}

PyObject*
Echo
::operator()()
{
    PyObject* parent_result = this->SCU::operator()();
    if(parent_result==NULL)
    {
        return NULL;
    }

    T_ASC_PresentationContextID pcid = this->_scu.findPresentationContextID(this->_abstractSyntax, "");
    if (pcid == 0)
    {
        PyErr_SetString(this->_medipy_base_exception, "No adequate Presentation Contexts");
        return NULL;

    }

    OFCondition cond=this->_scu.sendECHORequest(pcid);
    if (cond.bad())
    {
        PyErr_SetString(this->_medipy_base_exception,cond.text());
        return NULL;
    }

    return Py_None;
}
