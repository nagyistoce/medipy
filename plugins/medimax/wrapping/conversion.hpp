/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef medimax3_conversion_hpp
#define medimax3_conversion_hpp

#include <config.h>

#include <string>
#include <typeinfo>
#include <utility>

#include <Python.h>
#include <numpy/arrayobject.h>

#include <noyau/imx_3d.h>
#include <noyau/imx_types.h>

#define IS_SAME_TYPE_MACRO(t1, t2) \
	(typeid(t1) == typeid(t2))

#define TYPE_MATCH_MACRO(c_type, reference_c_type, array_type, reference_array_type) \
	(IS_SAME_TYPE_MACRO(c_type, reference_c_type) && (array_type==reference_array_type))

namespace MedimaxWrapper
{

/**
 * Execute a functor on the array, passing the correct element type to the
 * functor as template parameter.
 */
template<typename TFunctor>
typename TFunctor::return_type
executeArrayFunctor(TFunctor & functor, PyArrayObject* array);

/**
 * Return a grphic3d* holding the data from given array.
 * icomp and rcoeff are correctly positioned in this image
 */
struct ArrayToGrphic3DFunctor
{
    typedef grphic3d* return_type;

    // Min and max values in the array, set by call
    std::pair<float, float> minAndMax;

    template<typename T>
    return_type call(PyArrayObject * array);
};

/**
 * Return a grphic3d* holding the data and metadata from given image.
 */
DllExport grphic3d* imageToGrphic3D(PyObject* object);

/**
 * Modify given medipyImage to be a copy of grphic3d
 */
DllExport void grphic3DToImage(grphic3d* grphicImage, PyObject* medipyImage);

/**
 * Allocate a grphic3d structure according to the given array.
 */
DllExport grphic3d* allocateGrphic3d(PyArrayObject* array);

/**
 * Create a struct s_dti from a base.DTI
 */
//DllExport t_dti dtiFromPython(PyObject* object);

/**
 * Create a base.DTI from a struct s_dti
 */
//DllExport PyObject* dtiToPython(t_dti const & patient);

/**
 * Create a struct s_patient from a base.Patient
 */
//DllExport Patient patientFromPython(PyObject* object);

/**
 * Create a base.Patient from a struct s_patient
 */
//DllExport PyObject* patientToPython(Patient const & patient);


/**
 * Find the icomp and rcoeff(=2^icomp) values that allow to fit the values of
 * given range in the range of values of the type TYPEMRI3D.
 */
DllExport float computeIcomp(float min, float max);


////----- Change Detection algorithm (Marcel Bosc & Christophe Felder--------
//DllExport AlgorithmParameters changedetectionparametersFromPython(PyObject* object);



/**
 * Return the value held by the numpy array, cast as a float.
 */
struct GetValueFunctor
{
	typedef float return_type;

	unsigned int x,y,z;

	template<typename T>
	return_type call(PyArrayObject* array) const;
};

/**
 * Return the min and max of the array, as floats.
 *
 * The first element of the pair is the min, the second is the max.
 */
struct MinMaxFunctor
{
    typedef std::pair<float, float> return_type;

    template<typename T>
    return_type call(PyArrayObject * array) const;
};

/**
 * Copy the data from the Python array to the grphic3d.
 */
struct CopyDataToGrphic3dFunctor
{
    typedef void return_type;

    grphic3d* dest;

    template<typename T>
    return_type call(PyArrayObject * array) const;
};

/**
 * Copy the data from the grphic3d to the Python array.
 */
struct CopyDataFromGrphic3dFunctor
{
    typedef void return_type;

    grphic3d* source;

    template<typename T>
    return_type call(PyArrayObject * array) const;
};

/**
 * Copy the spacing from the Python array to the grphic3d.
 */
struct CopySpacingToGrphic3dFunctor
{
	typedef void return_type;

	grphic3d* dest;

	template<typename T>
	return_type call(PyArrayObject * array) const;
};

/**
 * Copy the spacing from the grphic3d to the Python array.
 */
struct CopySpacingFromGrphic3dFunctor
{
	typedef void return_type;

	grphic3d* source;

	template<typename T>
	return_type call(PyArrayObject * array) const;
};

struct DTIGradientFromPythonFunctor
{
	typedef void return_type;
	DTIGradientFromPythonFunctor(t_dti & dti);
	t_dti & dti;

	template<typename T>
	return_type call(PyArrayObject * array);
};

struct DTIGradientToPythonFunctor
{
	typedef void return_type;
	DTIGradientToPythonFunctor(t_dti const & dti);
	t_dti const & dti;

	template<typename T>
	return_type call(PyArrayObject * array) const;
};




} //end namespace

#include "conversion.txx"

#endif // medimax3_conversion_hpp
