#ifndef _512a4521_0310_449e_9911_5829b39b589c
#define _512a4521_0310_449e_9911_5829b39b589c

#include <vector>

#include <Python.h>

#include <dcmtk/config/osconfig.h>
#include <dcmtk/dcmdata/dctk.h>

void generate_dicomdir_cpp(PyObject * files, PyObject * root, 
    PyObject * dicomdir, PyObject * patient_extra_attributes, 
    PyObject * study_extra_attributes, PyObject * series_extra_attributes);

std::vector<DcmTagKey> _convert_attributes_list(PyObject * list);

#endif // _512a4521_0310_449e_9911_5829b39b589c
