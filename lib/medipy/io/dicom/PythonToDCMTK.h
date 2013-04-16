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
                  DcmDataset & dataset, DcmTag const & tag) const;

    template<typename TInserter, typename TGetter>
    void _to_binary(PyObject * python_value, TGetter getter,
                    DcmDataset & dataset, DcmTag const & tag,
                    TInserter inserter) const;

    void _to_raw(PyObject * python_value, DcmDataset & dataset,
                 DcmTag const & tag) const;

    void _to_number_string(PyObject * python_value, 
                           DcmDataset & dataset, DcmTag const & tag) const;
};

#endif // _30f65a36_44ce_495d_9270_d16009ee2bcf
