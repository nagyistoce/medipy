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
#include <string>

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
    
    /// @brief Return the Specific Character Set.
    std::string const & get_specific_character_set() const;
    
    /// @brief Set the Specific Character Set.
    void set_specific_character_set(std::string const & charset);
    
    /// @brief Convert the DCMTK dataset to a MediPy DataSet.
    PyObject * operator()(DcmObject * dataset);

private :
    /// @brief Map from DCMTK VR to MediPy VR class.
    std::map<DcmEVR, PyObject *> _medipy_io_dicom_vr;
    
    /// @brief MediPy DataSet class.
    PyObject * _medipy_io_dicom_DataSet;
    
    /// @brief MediPy Tag class.
    PyObject * _medipy_io_dicom_Tag;
    
    /// @brief Specific Character Set.
    std::string _specific_character_set;
    
    /// @brief Python encoding equivalent to the Specific Character Set.
    std::string _python_encoding;
    
    /// @brief Add the DCMTK element to the Python DataSet.
    void _add_element(DcmObject * element, PyObject * python_dataset) const;
    
    /// @brief Convert a DCMTK tag to a Python tag.
    PyObject * _to_python_tag(DcmTagKey const & dcmtk_tag) const;
    
    /// @brief Convert data from a DICOM element to Python.
    template<DcmEVR VVR>
    PyObject * _to_python(DcmObject * element) const;
    
    /**
     * @brief Convert data from a text DICOM element to Python.
     *
     * This is used for AE, AS, CS, DA, DT, LO, LT, PN, SH, ST, TM, UI and UT.
     */
    PyObject * _to_python_text(DcmByteString * element, bool use_utf8) const;
    
    /// @brief Convert data from a text DICOM element to Python.
    PyObject * _to_python_text(OFString const & value, DcmEVR const & vr, bool use_utf8) const;
    
    /**
     * @brief Convert data from a binary DICOM element to Python.
     *
     * This is used for OB, OF, OW and UN.
     */
    PyObject * _to_python_binary(DcmElement * element) const;

    /**
     * @brief Convert data from a numeric DICOM element to Python.
     *
     * This is used for FD, FL, SL, SS, UL and US.
     */
    template<typename TValue>
    PyObject * _to_python_number(DcmElement * element, 
                                 OFCondition (DcmElement::*getter)(TValue &, unsigned long)) const;
    
    /// @brief Convert data from a numeric DICOM element to Python.
    template<typename TValue>
    PyObject * _to_python_number(TValue const & value, DcmEVR const & valid_vr) const;
};

#include "DCMTKToPython.txx"

#endif // _96c8515b_6fb2_41fc_94d1_a1b3df5cadea
