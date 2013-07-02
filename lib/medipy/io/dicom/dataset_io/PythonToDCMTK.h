/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _30f65a36_44ce_495d_9270_d16009ee2bcf
#define _30f65a36_44ce_495d_9270_d16009ee2bcf

#include <Python.h>

#include <string>

#include <dcmtk/config/osconfig.h>
#include <dcmtk/dcmdata/dctk.h>

/**
 * @brief Convert a Python DataSet to a DCMTK DataSet.
 */
class PythonToDCMTK
{
public :

    PythonToDCMTK();
    PythonToDCMTK(PythonToDCMTK const & other);
    ~PythonToDCMTK();
    PythonToDCMTK & operator=(PythonToDCMTK const & other);

    DcmDataset operator()(PyObject * python_dataset);
private :
    /// @brief Add the Python element to the DCMTK DataSet.
    void _add_element(PyObject * python_tag, PyObject * python_value, DcmDataset & dataset);

    /// @brief Convert data from a Python element to DCMTK.
    template<DcmEVR VVR>
    void _to_dcmtk(PyObject * python_element, 
                   DcmDataset & dataset, DcmTag const & tag) const;

    void _to_text(PyObject * python_value, bool use_utf8, char padding,
                  DcmDataset & dataset, DcmElement * element) const;

    template<typename TGetter, typename TInserter>
    void _to_binary(PyObject * python_value, TGetter getter,
                    DcmDataset & dataset, DcmElement * element,
                    TInserter inserter) const;

    void _to_raw_8(PyObject * python_value, DcmDataset & dataset,
                   DcmElement * element) const;
    
    void _to_raw_16(PyObject * python_value, DcmDataset & dataset,
                    DcmElement * element) const;

    void _to_number_string(PyObject * python_value, 
                           DcmDataset & dataset, DcmElement * element) const;
    
    void _to_sequence(PyObject * python_value,
                      DcmDataset & dataset, DcmSequenceOfItems * element) const;
};

#endif // _30f65a36_44ce_495d_9270_d16009ee2bcf
