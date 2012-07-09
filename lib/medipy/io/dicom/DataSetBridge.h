#ifndef _eab6ea72_375e_4eb8_afbc_ef22e0eddbfa
#define _eab6ea72_375e_4eb8_afbc_ef22e0eddbfa

#include <Python.h>

#include <string>

#include <gdcmDataElement.h>
#include <gdcmDataSet.h>
#include <gdcmVR.h>

class DataSetBridge
{
public :
    DataSetBridge(gdcm::DataSet const & dataset);
    DataSetBridge(PyObject* dataset);

    gdcm::DataSet const & get_gdcm() const;
    void set_gdcm(gdcm::DataSet const & dataset);

    std::string const & get_encoding() const;
    void set_encoding(std::string const & encoding);

    PyObject* get_python() const;
    void set_python(PyObject* dataset);

    PyObject* to_python();
    PyObject* to_python(gdcm::DataElement const & data_element) const;

    // gdcm::DataSet const & to_gdcm();
private :
    gdcm::DataSet _gdcm_dataset;
    PyObject* _python_dataset;

    std::string _encoding;

    PyObject* _to_python(char const * begin, char const * end, gdcm::VR const & vr) const;
    unsigned long _get_length(char const * begin, char const * end, gdcm::VR const & vr) const;
};

#endif // _eab6ea72_375e_4eb8_afbc_ef22e0eddbfa
