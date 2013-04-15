/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _96c8515b_6fb2_41fc_94d1_a1b3df5cadea
#define _96c8515b_6fb2_41fc_94d1_a1b3df5cadea

#include <Python.h>

#include <map>

#include <dcmtk/config/osconfig.h>
#include <dcmtk/dcmdata/dctk.h>

/**
 * @brief Convert a DCMTK DataSet to a Python object.
 */
class DCMTKToPython
{
public :
    DCMTKToPython();
    DCMTKToPython(DCMTKToPython const & other);
    ~DCMTKToPython();
    DCMTKToPython & operator=(DCMTKToPython const & other);
    
    PyObject * operator()(DcmObject * dataset);

private :
    template<typename TValue>
    struct ElementValueGetter
    {
        typedef OFCondition (DcmElement::*Type)(TValue &, unsigned long);
    };
    
    std::map<DcmEVR, PyObject *> _medipy_io_dicom_vr;
    PyObject * _medipy_io_dicom_DataSet;
    PyObject * _medipy_io_dicom_Tag;
    
    /// @brief Convert data from a DICOM element to Python.
    template<DcmEVR VVR>
    PyObject * _to_python(DcmObject * element) const;
    
    /**
     * @brief Convert data from a text DICOM element to Python.
     *
     * This is used for AE, AS, CS, DA, DT, LO, LT, PN, SH, ST, TM, UI, UT
     */
    PyObject * _to_python_text(DcmByteString * element, bool use_utf8) const;
    
    /**
     * @brief Convert data from a binary DICOM element to Python.
     *
     * This is used for OB, OF, OW, UN
     */
    PyObject * _to_python_binary(DcmElement * element) const;

    /**
     * @brief Convert data from a numeric DICOM element to Python.
     *
     * This is used for FD, FL, SL, SS, UL, US
     */
    template<typename TValue>
    PyObject * _to_python_number(DcmElement * element, 
                                 OFCondition (DcmElement::*getter)(TValue &, unsigned long)) const;
    
    /**
     * @brief 
     * Since _to_bson is specialized and instantiated in _add_element,
     * this function must be declared after the the specializations.
     */
    void _add_element(DcmObject * element, PyObject * python_dataset) const;
};

#include "DCMTKToPython.txx"

#endif // _96c8515b_6fb2_41fc_94d1_a1b3df5cadea
