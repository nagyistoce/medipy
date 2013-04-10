/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef medimax3_conversion_txx
#define medimax3_conversion_txx

#include "conversion.hpp"

#include <algorithm>
#include <iostream>
#include <cmath>
#include <numeric>
#include <utility>
#include <sstream>
#include <stdexcept>

#include <numpy/arrayobject.h>

namespace MedimaxWrapper
{

template<typename TFunctor>
typename TFunctor::return_type
executeArrayFunctor(TFunctor & functor, PyArrayObject* array)
{
	int const array_type = PyArray_TYPE(array);

	if(array_type == NPY_BYTE)
	{
		return functor.template call<signed char>(array);
	}
	else if(array_type == NPY_UBYTE)
	{
		return functor.template call<unsigned char>(array);
	}
	else if(array_type == NPY_SHORT)
	{
		return functor.template call<signed short>(array);
	}
	else if(array_type == NPY_USHORT)
	{
		return functor.template call<unsigned short>(array);
	}
	else if(array_type == NPY_INT)
	{
		return functor.template call<signed int>(array);
	}
	else if(array_type == NPY_UINT)
	{
		return functor.template call<unsigned int>(array);
	}
	else if(array_type == NPY_LONG)
	{
		return functor.template call<signed long>(array);
	}
	else if(array_type == NPY_ULONG)
	{
		return functor.template call<unsigned long>(array);
	}
	else if(array_type == NPY_LONGLONG)
	{
		return functor.template call<signed long long>(array);
	}
	else if(array_type == NPY_ULONGLONG)
	{
		return functor.template call<unsigned long long>(array);
	}
	else if(array_type == NPY_FLOAT)
	{
		return functor.template call<float>(array);
	}
	else
	{
		std::stringstream message;
		message << "Cannot operate on PyArray of type " << array_type;
		throw std::runtime_error(message.str());
	}
}

template<typename T>
ArrayToGrphic3DFunctor::return_type
ArrayToGrphic3DFunctor
::call(PyArrayObject * array)
{
    grphic3d* result = allocateGrphic3d(array);

    // Find min and max. Will be used for icomp, rcoeff and meta-data
    MinMaxFunctor minMaxFunctor;
    this->minAndMax = executeArrayFunctor(minMaxFunctor, array);

    // icomp : set to 0 (i.e. rcoeff=1) if TYPEMRI3D and array_type match,
    // compute from min and max otherwise
    int const arrayType = PyArray_TYPE(array);
    if(TYPE_MATCH_MACRO(TYPEMRI3D, signed char, arrayType, NPY_BYTE) ||
       TYPE_MATCH_MACRO(TYPEMRI3D, unsigned char, arrayType, NPY_UBYTE) ||
       TYPE_MATCH_MACRO(TYPEMRI3D, signed short, arrayType, NPY_SHORT) ||
       TYPE_MATCH_MACRO(TYPEMRI3D, unsigned short, arrayType, NPY_USHORT) ||
       TYPE_MATCH_MACRO(TYPEMRI3D, signed int, arrayType, NPY_INT) ||
       TYPE_MATCH_MACRO(TYPEMRI3D, unsigned int, arrayType, NPY_INT) ||
       TYPE_MATCH_MACRO(TYPEMRI3D, signed long, arrayType, NPY_LONG) ||
       TYPE_MATCH_MACRO(TYPEMRI3D, unsigned long, arrayType, NPY_LONG) ||
       TYPE_MATCH_MACRO(TYPEMRI3D, signed long long, arrayType, NPY_LONGLONG) ||
       TYPE_MATCH_MACRO(TYPEMRI3D, unsigned long long, arrayType, NPY_ULONGLONG) ||
       TYPE_MATCH_MACRO(TYPEMRI3D, float, arrayType, NPY_FLOAT)
    )
    {
        result->icomp = 0;
    }
    else
    {
        result->icomp = MedimaxWrapper::computeIcomp(this->minAndMax.first, this->minAndMax.second);
    }
    result->rcoeff = std::pow(2.f, result->icomp);

    result->min_pixel = this->minAndMax.first / result->rcoeff;
    result->max_pixel = this->minAndMax.second / result->rcoeff;

    // Data

    MedimaxWrapper::CopyDataToGrphic3dFunctor copyDataFunctor;
    copyDataFunctor.dest = result;
    MedimaxWrapper::executeArrayFunctor(copyDataFunctor, array);

    return result;
}

template<typename T>
MinMaxFunctor::return_type
MinMaxFunctor
::call(PyArrayObject * array) const
{
	unsigned long nb_elements = std::accumulate(
			array->dimensions, array->dimensions+array->nd,
			1, std::multiplies<int>());
	T* begin = reinterpret_cast<T*>(array->data);
	T* end = begin+nb_elements;
	float min = static_cast<float>(*std::min_element(begin, end));
	float max = static_cast<float>(*std::max_element(begin, end));
	return std::make_pair(min, max);
}

template<typename T>
GetValueFunctor::return_type
GetValueFunctor
::call(PyArrayObject * array) const
{
	T* p = reinterpret_cast<T*>(PyArray_GETPTR3(array, this->z, this->y, this->x));
	return float(*p);
}

template<typename T>
CopyDataToGrphic3dFunctor::return_type
CopyDataToGrphic3dFunctor
::call(PyArrayObject * array) const
{
	T* data = reinterpret_cast<T*>(array->data);
	for(unsigned int z=0; z<this->dest->depth; ++z)
	{
		for(unsigned int y=0; y<this->dest->height; ++y)
		{
			for(unsigned int x=0; x<this->dest->width; ++x, ++data)
			{
				float const value =  float(*data);
				TYPEMRI3D const dest_value = (TYPEMRI3D)(value/this->dest->rcoeff+0.5);
				this->dest->mri[x][y][z] = dest_value;
			}
		}
	}
}

template<typename T>
CopyDataFromGrphic3dFunctor::return_type
CopyDataFromGrphic3dFunctor
::call(PyArrayObject * array) const
{
	// In Medimax, the z coordinates varies the most rapidly
	// In Numpy C-type arrays (and ImLib and ITK), the x coordinate varies
	// the most rapidly.
	// A quick test indicates that it is faster to use Medimax order and generate
	// cache misses in the Numpy array than the opposite.
	T* data = reinterpret_cast<T*>(array->data);
	for(unsigned int x=0; x<this->source->width; ++x)
	{
		for(unsigned int y=0; y<this->source->height; ++y)
		{
			for(unsigned int z=0; z<this->source->depth; ++z, ++data)
			{
				TYPEMRI3D const source_value = this->source->mri[x][y][z];
				float const float_value =  float(source_value)*this->source->rcoeff;
				T const dest_value = T(float_value);
				T* dest_ptr = reinterpret_cast<T*>(PyArray_GETPTR3(array, z,y,x));
				*dest_ptr = dest_value;
			}
		}
	}
}

template<typename T>
CopySpacingToGrphic3dFunctor::return_type
CopySpacingToGrphic3dFunctor
::call(PyArrayObject * array) const
{
	int const spacing_array_size = array->dimensions[0];
	T* data = (T*) array->data;
	this->dest->dx = data[spacing_array_size-1];
	this->dest->dy = (spacing_array_size>=2)?data[spacing_array_size-2]:1;
	this->dest->dz = (spacing_array_size>=3)?data[spacing_array_size-3]:1;
}

template<typename T>
CopySpacingFromGrphic3dFunctor::return_type
CopySpacingFromGrphic3dFunctor
::call(PyArrayObject * array) const
{
	T* data = (T*) array->data;
	data[0] = this->source->dz;
	data[1] = this->source->dy;
	data[2] = this->source->dx;
}

DTIGradientFromPythonFunctor
::DTIGradientFromPythonFunctor(t_dti & dti)
: dti(dti)
{
	// Nothing else
}

template<typename T>
DTIGradientFromPythonFunctor::return_type
DTIGradientFromPythonFunctor
::call(PyArrayObject * array)
{
	T* data = reinterpret_cast<T*>(array->data);
	this->dti.gx = data[0];
	this->dti.gy = data[1];
	this->dti.gz = data[2];
}


DTIGradientToPythonFunctor
::DTIGradientToPythonFunctor(t_dti const & dti)
: dti(dti)
{
	// Nothing else
}

template<typename T>
DTIGradientToPythonFunctor::return_type
DTIGradientToPythonFunctor
::call(PyArrayObject * array) const
{
	T* data = reinterpret_cast<T*>(array->data);
	data[0] = this->dti.gx;
	data[1] = this->dti.gy;
	data[2] = this->dti.gz;
}

}

#endif // medimax3_conversion_txx
