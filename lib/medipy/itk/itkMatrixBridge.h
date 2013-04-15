/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _f9007e4f_e53e_45ab_9997_e289ee31ae8f
#define _f9007e4f_e53e_45ab_9997_e289ee31ae8f

#include <itkMatrix.h>

// The python header defines _POSIX_C_SOURCE without a preceding #undef
#undef _POSIX_C_SOURCE
#include <Python.h>
#include <numpy/arrayobject.h>

namespace itk
{

template<typename TMatrix>
class MatrixBridge
{
public:
    /** Standard class typedefs */
    typedef MatrixBridge Self;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<const Self> ConstPointer;

    /// Type of the matrix from where the buffer will be converted
    typedef TMatrix MatrixType;

    /**
     * Get an Array with the content of the matrix
     */
    static PyObject* GetArrayFromMatrix(MatrixType const & matrix);

    /**
     * Get an ITK matrix from a Python array
     */
    static MatrixType GetMatrixFromArray(PyObject* obj);

protected:
    typedef enum PyArray_TYPES PyArrayType;
    static PyArrayType GetPyType(void);

    MatrixBridge(const Self&); // Not implemented.
    void operator=(const Self&); // Not implemented.
};

}

#include "itkMatrixBridge.txx"

#endif // _f9007e4f_e53e_45ab_9997_e289ee31ae8f
