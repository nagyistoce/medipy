/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _07ff966d_f8c1_4b24_b2ae_e7f068e6ed91
#define _07ff966d_f8c1_4b24_b2ae_e7f068e6ed91

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
    int item_type = Self::GetPyType();
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
    int element_type = Self::GetPyType();
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

#define MATRIX_BRIDGE_SET_ITEM_TYPE(image_type, c_type, item_type, python_type) \
    if(typeid(image_type) == typeid(c_type)) { item_type = python_type; }

template<typename TMatrix>
typename MatrixBridge<TMatrix>::PyArrayType
MatrixBridge<TMatrix>
::GetPyType(void)
{
    PyArrayType item_type=PyArray_NOTYPE;
    typedef typename PixelTraits<typename TMatrix::ValueType>::ValueType ScalarType;

    MATRIX_BRIDGE_SET_ITEM_TYPE(ScalarType, double, item_type, PyArray_DOUBLE)
    MATRIX_BRIDGE_SET_ITEM_TYPE(ScalarType, float, item_type, PyArray_FLOAT)
    MATRIX_BRIDGE_SET_ITEM_TYPE(ScalarType, unsigned char, item_type, PyArray_UBYTE)
    MATRIX_BRIDGE_SET_ITEM_TYPE(ScalarType, unsigned short, item_type, PyArray_USHORT)
    MATRIX_BRIDGE_SET_ITEM_TYPE(ScalarType, unsigned int, item_type, PyArray_UINT)
    MATRIX_BRIDGE_SET_ITEM_TYPE(ScalarType, unsigned long, item_type, PyArray_ULONG)
    MATRIX_BRIDGE_SET_ITEM_TYPE(ScalarType, char, item_type, PyArray_BYTE)
    MATRIX_BRIDGE_SET_ITEM_TYPE(ScalarType, short, item_type, PyArray_SHORT)
    MATRIX_BRIDGE_SET_ITEM_TYPE(ScalarType, int, item_type, PyArray_INT)
    MATRIX_BRIDGE_SET_ITEM_TYPE(ScalarType, long, item_type, PyArray_LONG)

    if(item_type == PyArray_NOTYPE)
    {
        throw std::runtime_error("Type currently not supported");
    }
    return item_type;
}

#undef MATRIX_BRIDGE_SET_ITEM_TYPE

}

#endif // _07ff966d_f8c1_4b24_b2ae_e7f068e6ed91
