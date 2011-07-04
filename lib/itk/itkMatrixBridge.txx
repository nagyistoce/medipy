/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef itk_itkmatrixbridge_txx
#define itk_itkmatrixbridge_txx

#include "itkMatrixBridge.h"

#include <typeinfo>
#include <itkPixelTraits.h>

namespace itk

{

template<typename TMatrix>
PyObject*
MatrixBridge<TMatrix>
::GetArrayFromMatrix(MatrixType const & matrix)
{
    npy_intp dimensions[] = { matrix.RowDimensions, matrix.ColumnDimensions };
    int item_type = GetPyType();
    PyObject * obj = PyArray_New(&PyArray_Type, 2, dimensions, item_type, NULL, NULL, 0,
                                 NPY_CARRAY, NULL);

    for(int r=0; r<matrix.RowDimensions; ++r)
    {
        for(int c=0; c<matrix.ColumnDimensions; ++c)
        {
            *reinterpret_cast<typename TMatrix::ValueType*>(PyArray_GETPTR2(obj, r, c)) = matrix(r,c);
        }
    }

    return obj;
}


template<typename TMatrix>
typename MatrixBridge<TMatrix>::MatrixType
MatrixBridge<TMatrix>
::GetMatrixFromArray(PyObject* obj)
{
    MatrixType matrix;

    // Get the array object
    int element_type = GetPyType();
    PyArrayObject * parray = (PyArrayObject *) PyArray_ContiguousFromAny(obj, element_type, 2, 2);

    if(parray == NULL)
    {
        throw std::runtime_error("Contiguous array couldn't be created from input python object");
    }

    for(int r=0; r<matrix.RowDimensions; ++r)
    {
        for(int c=0; c<matrix.ColumnDimensions; ++c)
        {
            matrix(r,c) = *reinterpret_cast<typename TMatrix::ValueType*>(PyArray_GETPTR2(obj, r, c));
        }
    }

    return matrix;
}

template<typename TMatrix>
typename MatrixBridge<TMatrix>::PyArrayType
MatrixBridge<TMatrix>
::GetPyType(void)
{

  PyArrayType item_type;
  typedef typename PixelTraits<typename TMatrix::ValueType>::ValueType ScalarType;
  if(typeid(ScalarType) == typeid(double))
    {
    item_type = PyArray_DOUBLE;
    }
  else if(typeid(ScalarType) == typeid(float))
    {
    item_type = PyArray_FLOAT;
    }
  else if(typeid(ScalarType) == typeid(long))
    {
    item_type = PyArray_LONG;
    }
  else if(typeid(ScalarType) == typeid(unsigned long))
    {
#ifdef NDARRAY_VERSION
    item_type = PyArray_ULONG;
#else
    throw std::runtime_error("Type currently not supported");
#endif
    }
  else if(typeid(ScalarType) == typeid(int))
    {
    item_type = PyArray_INT;
    }
  else if(typeid(ScalarType) == typeid(unsigned int))
    {
    item_type = PyArray_UINT;
    }
  else if(typeid(ScalarType) == typeid(short))
    {
    item_type = PyArray_SHORT;
    }
  else if(typeid(ScalarType) == typeid(unsigned short))
    {
    item_type = PyArray_USHORT;
    }
  else if(typeid(ScalarType) == typeid(signed char))
    {
    item_type = PyArray_BYTE;
    }
  else if(typeid(ScalarType) == typeid(unsigned char))
    {
    item_type = PyArray_UBYTE;
    }
  else
    {
    item_type = PyArray_NOTYPE;
    throw std::runtime_error("Type currently not supported");
    }
  return item_type;
}

}

#endif // itk_itkmatrixbridge_txx
