/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _bcd8dc5f_fe20_4542_99ff_db406e1c0329
#define _bcd8dc5f_fe20_4542_99ff_db406e1c0329

#include <itkImage.h>
#include <itkImportImageFilter.h>

// The python header defines _POSIX_C_SOURCE without a preceding #undef
#undef _POSIX_C_SOURCE
#include <Python.h>
#include <numpy/arrayobject.h>

namespace itk
{

template<typename TImage>
class NumpyBridge
{
public:
    /** Standard class typedefs */
    typedef NumpyBridge Self;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<const Self> ConstPointer;

    /// Type of the image from where the buffer will be converted
    typedef TImage                              ImageType;
    typedef typename ImageType::PixelType       PixelType;
    typedef typename ImageType::IOPixelType     IOPixelType;
    typedef typename ImageType::SizeType        SizeType;
    typedef typename ImageType::IndexType       IndexType;
    typedef typename ImageType::RegionType      RegionType;
    typedef typename ImageType::PointType       PointType;
    typedef typename ImageType::SpacingType     SpacingType;
    typedef typename ImageType::Pointer         ImagePointer;

     /** Image dimension. */
    itkStaticConstMacro(ImageDimension, unsigned int, ImageType::ImageDimension);

    /// Type of the import image filter
    typedef ImportImageFilter<PixelType, ImageDimension> ImporterType;

    typedef typename ImporterType::Pointer ImporterPointer;

    /**
     * Get an Array with the content of the image buffer
     */
    static PyObject* GetArrayFromImage(ImageType* image, bool transferOwnership);

    /**
     * Get an ITK image from a Python array
     */
    static const ImagePointer GetImageFromArray(PyObject* obj, bool transferOwnership);

    /**
     * Test if an ITK image and a Python array share the same buffer
     */
    static bool IsBufferShared(PyObject* obj, ImageType* image);

protected:
    typedef enum PyArray_TYPES PyArrayType;
    static PyArrayType GetPyType(void);
    static void set_array_ownership(PyArrayObject* array, bool value);
    static PyArrayObject* GetPyArrayObject(PyObject* obj);

    NumpyBridge(const Self&); // Not implemented.
    void operator=(const Self&); // Not implemented.
};

}

#include "itkNumpyBridge.txx"

#endif // _bcd8dc5f_fe20_4542_99ff_db406e1c0329
