/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "conversion.hpp"

#include <cmath>
#include <sstream>

#include <Python.h>
#include <boost/python.hpp>

#include <numpy/arrayobject.h>

#include <noyau/imx_3d.h>
#include <noyau/imx_types.h>


namespace MedimaxWrapper
{

void inner_initialize()
{
    // import_array is a macro expanding to the body of a void function
    import_array();
}

bool initialize()
{
    inner_initialize();
    return true;
}

grphic3d* allocateGrphic3d(PyArrayObject* array)
{
    unsigned int const nd = array->nd;

    unsigned int const width = array->dimensions[nd-1];
    unsigned int const height = (nd>=2)?array->dimensions[nd-2]:1;
    unsigned int const depth = (nd>=3)?array->dimensions[nd-3]:1;
    grphic3d* result = cr_grphic3d_modif(width, height, depth, 0, 1, 0);

    return result;
}

grphic3d* imageToGrphic3D(PyObject* object)
{
    grphic3d* result;

    if (object == Py_None) 
        {
        result=NULL;
        return result;
        }
    
    PyObject* data = PyObject_GetAttrString(object, "data");
    if(!PyArray_ISCONTIGUOUS(data))
    {
        throw std::runtime_error("Cannot convert a non C-contiguous array");
    }

    // Copy data to grphic3d
    PyArrayObject* array = (PyArrayObject*)(
        PyArray_FromAny(data, NULL, 1, 3, NPY_CARRAY, NULL));

    ArrayToGrphic3DFunctor arrayToGrphic3D;
    result = executeArrayFunctor(arrayToGrphic3D, array);

    Py_DECREF(data); // From PyArray_FromAny
    Py_DECREF(data); // From PyObject_GetAttrString

    // Copy meta-data to grphic3d
    // TODO: patient information, DTI, ...
    PyObject* spacing = PyObject_GetAttrString(object, "spacing");
    array = (PyArrayObject*)(PyArray_FromAny(spacing, NULL, 1, 1, NPY_CARRAY, NULL));
    MedimaxWrapper::CopySpacingToGrphic3dFunctor copySpacingFunctor;
    copySpacingFunctor.dest = result;
    MedimaxWrapper::executeArrayFunctor(copySpacingFunctor, array);
    Py_DECREF(spacing); // From PyArray_FromAny
    Py_DECREF(spacing); // From PyObject_GetAttrString

    return result;
}

void grphic3DToImage(grphic3d* grphicImage, PyObject* medipyImage)
{
    
    if (grphicImage==NULL)
    {
        medipyImage=Py_None;
    }
    else 
    {
    // Copy data from grphic3d
    PyObject* data = PyObject_GetAttrString(medipyImage, "data");
    PyArrayObject* array = (PyArrayObject*)(
        PyArray_FromAny(data, NULL, 1, 3, NPY_CARRAY, NULL));

    npy_intp newshape[] = { grphicImage->depth, grphicImage->height, grphicImage->width };
    PyArray_Dims newdims = {newshape, 3};
    PyArray_Resize(array, &newdims, 0, PyArray_CORDER);

    MedimaxWrapper::CopyDataFromGrphic3dFunctor copyDataFunctor;
    copyDataFunctor.source = grphicImage;
    MedimaxWrapper::executeArrayFunctor(copyDataFunctor, array);

    Py_DECREF(data); // From PyArray_FromAny
    Py_DECREF(data); // From PyObject_GetAttrString

    // Copy meta-data from grphic3d
    // TODO: patient information, DTI, ...
    npy_intp dimension = 3;
    PyObject* spacing = PyArray_SimpleNew(1, &dimension, NPY_FLOAT);
    MedimaxWrapper::CopySpacingFromGrphic3dFunctor copySpacingFunctor;
    copySpacingFunctor.source = grphicImage;
    MedimaxWrapper::executeArrayFunctor(copySpacingFunctor, (PyArrayObject*)spacing);
    PyObject_SetAttrString(medipyImage, "spacing", spacing);
    }
}


//t_dti dtiFromPython(PyObject* object)
//{
//  t_dti dti;
//
//  PyObject* type = PyObject_GetAttrString(object, "type");
//  if(type != Py_None)
//  {
//      dti.is_dti = PyInt_AsLong(type);
//  }
//
//  PyObject* gradient = PyObject_GetAttrString(object, "gradient");
//  if(gradient != Py_None)
//  {
//      PyArrayObject* gradientArray = (PyArrayObject*)(
//              PyArray_FromAny(gradient, NULL, 1, 1, NPY_CARRAY, NULL));
//      if(gradientArray != NULL)
//      {
//          DTIGradientFromPythonFunctor functor(dti);
//          executeArrayFunctor(functor, gradientArray);
//      }
//  }
//
//  PyObject* b = PyObject_GetAttrString(object, "b");
//  if(b != Py_None)
//  {
//      dti.b = PyFloat_AsDouble(b);
//  }
//
//  PyObject* number_of_directions = PyObject_GetAttrString(object, "number_of_directions");
//  if(type != Py_None)
//  {
//      dti.nb_directions = PyInt_AsLong(number_of_directions);
//  }
//
//
//  return dti;
//}
//
//
//PyObject* dtiToPython(t_dti const & dti)
//{
//  PyRun_SimpleString("from medipy.base import DTI ; dti = DTI()");
//  PyObject* module = PyImport_AddModule("__main__");
//  PyObject* dict = PyModule_GetDict(module);
//  PyObject* pyDTI= PyDict_GetItemString(dict, "patient");
//
//  PyObject_SetAttrString(pyDTI, "type", PyInt_FromLong(dti.is_dti));
//
//  npy_intp dimension = 3;
//  PyObject* pyGradient = PyArray_SimpleNew(1, &dimension, NPY_FLOAT);
//  DTIGradientToPythonFunctor functor(dti);
//  executeArrayFunctor(functor, (PyArrayObject*)pyGradient);
//  PyObject_SetAttrString(pyDTI, "gradient", pyGradient);
//
//  PyObject_SetAttrString(pyDTI, "b", PyFloat_FromDouble(dti.b));
//  PyObject_SetAttrString(pyDTI, "number_of_directions", PyInt_FromLong(dti.nb_directions));
//
//  return pyDTI;
//}
//
//
//Patient patientFromPython(PyObject* object)
//{
//  Patient patient;
//
//  namespace bp = boost::python;
//
//  bp::object o(bp::handle<>(bp::borrowed(object)));
//  bp::dict metadata = bp::extract<bp::dict>(o);
//
//  if(metadata.has_key("patient"))
//  {
//      bp::dict metadataPatient = bp::extract<bp::dict>(metadata.get("patient"));
//      if(metadataPatient.has_key("patients_name"))
//        {
//          strncpy(patient.name,
//                  bp::extract<char*>(metadataPatient.get("patients_name")),
//                  256);
//        }
//      if(metadataPatient.has_key("patients_birth_date"))
//        {
//            strncpy(patient.d_birth,
//                    bp::extract<char*>(metadataPatient.get("patients_birth_date")),
//                    35);
//        }
//  }
//
//  if(metadata.has_key("study"))
//    {
//        bp::dict metadataStudy = bp::extract<bp::dict>(metadata.get("study"));
//        if(metadataStudy.has_key("referring_physicians_name"))
//        {
//            strncpy(patient.medecin,
//                    bp::extract<char*>(metadataStudy.get("referring_physicians_name")),
//                    35);
//        }
//    }
//
//  if(metadata.has_key("series"))
//    {
//        bp::dict metadataSeries = bp::extract<bp::dict>(metadata.get("series"));
//        if(metadataSeries.has_key("series_date"))
//        {
//            strncpy(patient.d_examen,
//                    bp::extract<char*>(metadataSeries.get("series_date")),
//                    35);
//        }
//    }
//
//  // TODO: process d_stim and percent_pem.
//
//  return patient;
//
//}
//
//
//PyObject* patientToPython(Patient const & patient)
//{
//    namespace bp = boost::python;
//
//    bp::object medipy = bp::import("medipy");
//    bp::object medipy_base = bp::import("medipy.base");
//
////    bp::dict globals;
////    bp::object pythonPatient = bp::exec("from medipy.base import Metadata\n"
////                                        "patient = Metadata()", globals);
//
//  PyRun_SimpleString("from medipy.base import Patient ; patient = Patient()");
//  PyObject* module = PyImport_AddModule("__main__");
//  PyObject* dict = PyModule_GetDict(module);
//  PyObject* pyPatient = PyDict_GetItemString(dict, "patient");
//
//  PyObject_SetAttrString(pyPatient, "name", PyString_FromString(patient.name));
//  PyObject_SetAttrString(pyPatient, "birth_date", PyString_FromString(patient.d_birth));
//  PyObject_SetAttrString(pyPatient, "physician_name", PyString_FromString(patient.medecin));
//  PyObject_SetAttrString(pyPatient, "exam_date", PyString_FromString(patient.d_examen));
//  PyObject_SetAttrString(pyPatient, "stimulation_date", PyString_FromString(patient.d_stim));
//  PyObject_SetAttrString(pyPatient, "stimulation_intensity", PyFloat_FromDouble(patient.percent_pem));
//
//  return pyPatient;
//}


float computeIcomp(float min, float max)
{
    float const epsilon=1e-5;

    float const max_abs = std::max(std::abs(min), max);

    // Find the optimal icomp and rcoeff
    float rcoeff, icomp;
    if(max_abs == 0.0)
    {
        rcoeff = 1.;
        icomp = 0.;
    }
    else
    {
        rcoeff = max_abs/MAXMRI3D;
        icomp = std::log(rcoeff)/std::log(2.f);
    }

    // Get the closest integer icomp
    if(icomp<=0.0)
    {
        // If fractional part is larger than epsilon, use the integer
        // immediately above, else round
        if(std::abs(icomp-std::floor(icomp+0.5))>epsilon)
        {
            icomp = std::ceil(icomp);
        }
        else
        {
            icomp = std::floor(icomp+0.5);
        }
    }
    else
    {
        // If fractional part is larger than epsilon, use the integer
        // immediately above, else round
        if(std::abs(icomp-std::floor(icomp+0.5))>epsilon)
        {
            icomp = std::ceil(icomp+1.0);
        }
        else
        {
            icomp = std::floor(icomp+1.0+0.5);
        }
    }

    return icomp;
}

}
