/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _6be648b4_a39e_41bb_9d88_a83e6873e3a1
#define _6be648b4_a39e_41bb_9d88_a83e6873e3a1

#include <cstring>
#include "itkNumpyBridge.h"
#include <itkPixelTraits.h>

// Deal with slight incompatibilites between NumPy (the future, hopefully),
// Numeric (old version) and Numarray's Numeric compatibility module (also old).
#ifndef NDARRAY_VERSION
// NDARRAY_VERSION is only defined by NumPy's arrayobject.h
// In non NumPy arrayobject.h files, PyArray_SBYTE is used instead of BYTE.
#define PyArray_BYTE PyArray_SBYTE
#endif

namespace itk
{

template<typename TImage>
PyObject *
NumpyBridge<TImage>
::GetArrayFromImage(ImageType* image, bool transferOwnership)
{
    if(!image)
    {
        throw std::runtime_error("Input image is null");
    }

    image->Update();

    // Determine number of dimensions of the numpy.ndarray
    unsigned int dataDimension=-1;
    unsigned int* dataShape=NULL;
    if(image->GetNameOfClass() == std::string("Image"))
    {
        dataDimension = 0;
        dataShape = new unsigned int[dataDimension];
    }
    else if(image->GetNameOfClass() == std::string("VectorImage"))
    {
        dataDimension = 1;
        dataShape = new unsigned int[dataDimension];
        dataShape[0] = image->GetNumberOfComponentsPerPixel();
    }
    else if(image->GetNameOfClass() == std::string("MatrixImage"))
    {
        dataDimension = 2;
        //dataShape = new unsigned int[dataDimension];
    }
    else
    {
        throw std::runtime_error(std::string("Cannot convert object of type") + image->GetNameOfClass());
    }

    unsigned int const arrayDimension = Self::ImageDimension+dataDimension;

    // Fill array dimensions
    npy_intp* dimensions = new npy_intp[arrayDimension];

    // First the image dimensions
    SizeType size = image->GetBufferedRegion().GetSize();
    for (unsigned int d = 0; d < ImageDimension; d++)
    {
        dimensions[ImageDimension - d - 1] = size[d];
    }
    // Then the (non-scalar) data dimensions
    for(unsigned int d=0; d<dataDimension; ++d)
    {
        dimensions[arrayDimension-d-1] = dataShape[d];
    }

    void* data = reinterpret_cast<void*>(image->GetBufferPointer());

    // Create the Python object, obj will not own the memory by default
    int item_type = Self::GetPyType();
    PyObject * obj = PyArray_SimpleNewFromData(arrayDimension, dimensions, item_type, data);
    delete[] dimensions;

    // Transfer ownership only if ITK image manages memory
    if(transferOwnership && image->GetPixelContainer()->GetContainerManageMemory())
    {
        Self::set_array_ownership(reinterpret_cast<PyArrayObject*>(obj), true);
        image->GetPixelContainer()->ContainerManageMemoryOff();
    }

    return obj;
}

template<typename TImage>
const typename NumpyBridge<TImage>::ImagePointer
NumpyBridge<TImage>
::GetImageFromArray(PyObject* obj, bool transferOwnership)
{
    PyArrayObject* array = Self::GetPyArrayObject(obj);

    // Fill image size, leave out the (non-scalar) data dimensions
    SizeType size;
    for(unsigned int d = 0; d<TImage::GetImageDimension(); d++)
    {
        size[TImage::GetImageDimension() - d - 1] = array->dimensions[d];
    }

    IndexType start;
    start.Fill(0);

    RegionType region;
    region.SetIndex(start);
    region.SetSize(size);

    PointType origin;
    origin.Fill(0.0);

    SpacingType spacing;
    spacing.Fill(1.0);

    ImagePointer result = TImage::New();

    unsigned int arrayDimension=0;
    if(result->GetNameOfClass() == std::string("Image"))
    {
        arrayDimension = TImage::GetImageDimension();
    }
    else if(result->GetNameOfClass() == std::string("VectorImage"))
    {
        arrayDimension = TImage::GetImageDimension()+1;
    }
    else if(result->GetNameOfClass() == std::string("MatrixImage"))
    {
        arrayDimension = TImage::GetImageDimension()+2;
    }
    else
    {
        throw std::runtime_error(std::string("Cannot convert object of type") + result->GetNameOfClass());
    }

    unsigned long numberOfPixels = 1;
    for(unsigned int d=0; d<arrayDimension; ++d)
    {
        numberOfPixels *= array->dimensions[d];
    }

    typedef itk::ImportImageContainer<unsigned long, IOPixelType> ContainerType;
    typename ContainerType::Pointer container = ContainerType::New();
    IOPixelType* data = (IOPixelType*) array->data;
    container->SetImportPointer(data, numberOfPixels, transferOwnership);

    result->SetRegions(region);
    result->SetOrigin(origin);
    result->SetSpacing(spacing);
    result->SetPixelContainer(container);
    if(result->GetNameOfClass() == std::string("VectorImage"))
    {
        result->SetNumberOfComponentsPerPixel(array->dimensions[arrayDimension-1]);
    }

    if(transferOwnership)
    {
        Self::set_array_ownership(array, false);
    }

    return result;
}


template<typename TImage>
bool
NumpyBridge<TImage>
::IsBufferShared(PyObject* obj, ImageType* image)
{
    PyArrayObject* array = NULL;
    try
    {
        PyArrayObject* array = Self::GetPyArrayObject(obj);
    }
    catch(...)
    {
        // Could not convert : buffer is not shared
        PyErr_Clear();
        return false;
    }

    if(array == NULL)
    {
        return false;
    }

    return (array->data==reinterpret_cast<char*>(image->GetBufferPointer()));
}

#define NUMPY_BRIDGE_SET_ITEM_TYPE(image_type, c_type, item_type, python_type) \
    if(typeid(image_type) == typeid(c_type)) { item_type = python_type; }

template<typename TImage>
typename NumpyBridge<TImage>::PyArrayType
NumpyBridge<TImage>
::GetPyType(void)
{
    PyArrayType item_type=PyArray_NOTYPE;
    typedef typename PixelTraits<IOPixelType>::ValueType ScalarType;

    NUMPY_BRIDGE_SET_ITEM_TYPE(ScalarType, double, item_type, PyArray_DOUBLE)
    NUMPY_BRIDGE_SET_ITEM_TYPE(ScalarType, float, item_type, PyArray_FLOAT)
    NUMPY_BRIDGE_SET_ITEM_TYPE(ScalarType, unsigned char, item_type, PyArray_UBYTE)
    NUMPY_BRIDGE_SET_ITEM_TYPE(ScalarType, unsigned short, item_type, PyArray_USHORT)
    NUMPY_BRIDGE_SET_ITEM_TYPE(ScalarType, unsigned int, item_type, PyArray_UINT)
    NUMPY_BRIDGE_SET_ITEM_TYPE(ScalarType, unsigned long, item_type, PyArray_ULONG)
    NUMPY_BRIDGE_SET_ITEM_TYPE(ScalarType, char, item_type, PyArray_BYTE)
    NUMPY_BRIDGE_SET_ITEM_TYPE(ScalarType, short, item_type, PyArray_SHORT)
    NUMPY_BRIDGE_SET_ITEM_TYPE(ScalarType, int, item_type, PyArray_INT)
    NUMPY_BRIDGE_SET_ITEM_TYPE(ScalarType, long, item_type, PyArray_LONG)

    if(item_type == PyArray_NOTYPE)
    {
        throw std::runtime_error("Type currently not supported");
    }
    return item_type;
}

#undef NUMPY_BRIDGE_SET_ITEM_TYPE

template<typename TImage>
void
NumpyBridge<TImage>
::set_array_ownership(PyArrayObject* array, bool value)
{
    if(value)
    {
        array->flags |= NPY_OWNDATA;
    }
    else
    {
        array->flags &= ~NPY_OWNDATA;
    }
}

template<typename TImage>
PyArrayObject*
NumpyBridge<TImage>
::GetPyArrayObject(PyObject* obj)
{
    if(strcmp(obj->ob_type->tp_name, "numpy.ndarray") != 0)
    {
        throw std::runtime_error(
            "Cannot convert a " + std::string(obj->ob_type->tp_name));
    }

    PyArrayObject* array = reinterpret_cast<PyArrayObject*>(obj);
    // Check that we have a C array
    if(!PyArray_ISCONTIGUOUS(array))
    {
        throw std::runtime_error("Cannot convert a non C-contiguous array");
    }
    // Check that array and image types match
    PyArrayType const py_type = Self::GetPyType();
    if(array->descr->type_num != Self::GetPyType())
    {
        std::cout << "array: " << array->descr->type_num << " "
                  << "image: " << Self::GetPyType() << std::endl;
        throw std::runtime_error("Array and image types do not match");
    }
    // Check that array and image dimensions match
    ImagePointer image = TImage::New();
    unsigned int arrayDimension=0;
    if(image->GetNameOfClass() == std::string("Image"))
    {
        arrayDimension = TImage::GetImageDimension();
    }
    else if(image->GetNameOfClass() == std::string("VectorImage"))
    {
        arrayDimension = TImage::GetImageDimension()+1;
    }
    else if(image->GetNameOfClass() == std::string("MatrixImage"))
    {
        arrayDimension = TImage::GetImageDimension()+2;
    }
    else
    {
        throw std::runtime_error(
            std::string("Cannot convert object of type")+image->GetNameOfClass());
    }
    if(arrayDimension != array->nd)
    {
        throw std::runtime_error("Array and image dimensions do not match");
    }

    return array;
}

} // namespace itk

#endif // _6be648b4_a39e_41bb_9d88_a83e6873e3a1
