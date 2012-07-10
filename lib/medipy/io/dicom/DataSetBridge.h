/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2012
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _eab6ea72_375e_4eb8_afbc_ef22e0eddbfa
#define _eab6ea72_375e_4eb8_afbc_ef22e0eddbfa

#include <Python.h>

#include <string>

#include <gdcmDataElement.h>
#include <gdcmDataSet.h>
#include <gdcmVR.h>

/**
 * @brief Convert between GDCM and Python DICOM Data Set.
 *
 * The GDCM Data Sets are instances of gdcm::DataSet, while Python Data Sets
 * are dictionaries of (tag, value). In Python Data Sets, the tags are the
 * numeric (group, element), and according to the VR, the values are :
 *   * AE, AS, CS, DA, DT, TM, UI : str objects, with leading and trailing
 *     whitespaces removed.
 *   * LO, LT, PN, SH, ST, UT : unicode objects, with leading and trailing
 *     whitespaces removed (at the exception of leading whitespaces of LT, ST
 *     and UT which are considered significant).
 *   * DS, FD, FL : float objects.
 *   * AT : tuples of two int objects.
 *   * IS, SL, SL, UL, US : int objects.
 *   * OB : NumPy arrays of unsigned 8-bits integers.
 *   * OF : NumPy arrays of 32-bits floats.
 *   * OW : NumPy arrays of unsigned 16-bits integers.
 *   * UN : str objects, to be used as sequence of bytes
 *   * SQ : lists of dictionaries, each item being a Python Data Set.
 */
class DataSetBridge
{
public :
    /// @brief Build a bridge with empty GDCM and Python Data Sets.
    DataSetBridge();
    /// @brief Build a bridge with given GDCM Data Set.
    DataSetBridge(gdcm::DataSet const & dataset);
    /// @brief Build a bridge with given Python Data Set.
    DataSetBridge(PyObject* dataset);

    /// @brief Return the GDCM Data Set.
    gdcm::DataSet const & get_gdcm() const;
    /// @brief Set the GDCM Data Set.
    void set_gdcm(gdcm::DataSet const & dataset);

    /// @brief Return the Python Data Set.
    PyObject* get_python() const;
    /// @brief Set the Python Data Set.
    void set_python(PyObject* dataset);

    /// @brief Return the Data Set encoding.
    std::string const & get_encoding() const;
    /**
     * @brief Set the Data Set encoding (must be a Python Standard Encoding),
     * cf. http://docs.python.org/library/codecs.html#standard-encodings
     */
    void set_encoding(std::string const & encoding);

    /// @brief Convert the GDCM Data Set to a Python Data Set.
    PyObject* to_python();

    /// @brief Convert a GDCM Data Element to a Python Data Element.
    PyObject* to_python(gdcm::DataElement const & data_element) const;

    // TODO : gdcm::DataSet const & to_gdcm();
private :
    /// @brief GDCM version of the Data Set.
    gdcm::DataSet _gdcm_dataset;

    /// @brief Python version of the Data Set : a dictionary of (tag, value).
    PyObject* _python_dataset;

    /// @brief Encoding of the Data Set.
    std::string _encoding;

    /// @brief Convert a DICOM Value Field of known VR to a Python object.
    PyObject* _to_python(char const * begin, char const * end, gdcm::VR const & vr) const;

    /// @brief Return the length of an object from a (subset of) DICOM Value Field of known VR.
    unsigned long _get_length(char const * begin, char const * end, gdcm::VR const & vr) const;
};

#endif // _eab6ea72_375e_4eb8_afbc_ef22e0eddbfa
