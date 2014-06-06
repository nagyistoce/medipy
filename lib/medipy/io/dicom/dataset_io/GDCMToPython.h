/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _2038797a_854e_474c_a392_42ab11ba8bb1
#define _2038797a_854e_474c_a392_42ab11ba8bb1

#include <Python.h>

#include <map>
#include <string>

#include <gdcmDataElement.h>
#include <gdcmDataSet.h>
#include <gdcmTag.h>
#include <gdcmVR.h>

/**
 * @brief Convert a GDCM DataSet to a Python object.
 */
class GDCMToPython
{
public :
    GDCMToPython();
    GDCMToPython(GDCMToPython const & other);
    ~GDCMToPython();
    GDCMToPython & operator=(GDCMToPython const & other);
    
    /// @brief Return the Specific Character Set.
    std::string const & get_specific_character_set() const;
    
    /// @brief Set the Specific Character Set.
    void set_specific_character_set(std::string const & charset);
    
    /// @brief Convert a GDCM dataset to a MediPy DataSet.
    PyObject * operator()(gdcm::DataSet const & dataset);
    
    /// @brief Convert a GDCM element to a MediPy DICOM value.
    PyObject * operator()(gdcm::DataElement const & element, gdcm::DataSet const & dataset) const;
    
    /// @brief Convert a GDCM tag to a MediPy tag.
    PyObject * operator()(gdcm::Tag const & tag) const;

private :
    typedef std::map<gdcm::VR::VRType, PyObject *> VRMap;
    
    /// @brief Map from GDCM VR to MediPy VR class.
    VRMap _medipy_io_dicom_vr;
    
    /// @brief MediPy DataSet class.
    PyObject * _medipy_io_dicom_DataSet;
    
    /// @brief MediPy Tag class.
    PyObject * _medipy_io_dicom_Tag;
    
    /// @brief Specific Character Set.
    std::string _specific_character_set;
    
    /// @brief Python encoding equivalent to the Specific Character Set.
    std::string _python_encoding;
    
    /// @brief Pixel Representation
    unsigned short _pixel_representation;
    
    /// @brief Convert data from a GDCM element to Python.
    PyObject * _to_python_value(
        char const * begin, char const * end, gdcm::VR const & vr) const;
    
    /// @brief Return the length of the first item from DICOM element.
    unsigned long _get_length(
        char const * begin, char const * end, gdcm::VR const & vr) const;
};

//#include "GDCMToPython.txx"

#endif // _96c8515b_6fb2_41fc_94d1_a1b3df5cadea
