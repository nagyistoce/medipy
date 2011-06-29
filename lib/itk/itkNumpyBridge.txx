/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef itk_itknumpybridge_txx
#define itk_itknumpybridge_txx

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

    unsigned int const arrayDimension = ImageDimension+dataDimension;

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
    int item_type = GetPyType();
    PyObject * obj = PyArray_SimpleNewFromData(arrayDimension, dimensions, item_type, data);
    delete[] dimensions;

    // Transfer ownership only if ITK image manages memory
    if(transferOwnership && image->GetPixelContainer()->GetContainerManageMemory())
    {
        Self::set_array_ownership(obj, true);
        image->GetPixelContainer()->ContainerManageMemoryOff();
    }

    return obj;
}


template<typename TImage>
const typename NumpyBridge<TImage>::ImagePointer
NumpyBridge<TImage>
::GetImageFromArray(PyObject* obj, bool transferOwnership)
{
    if(!PyArray_ISCONTIGUOUS(obj))
    {
        throw std::runtime_error("Cannot convert a non C-contiguous array");
    }

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

    // Get the array object
    int element_type = GetPyType();
    PyArrayObject * parray = (PyArrayObject *) PyArray_ContiguousFromAny(obj, element_type, arrayDimension, arrayDimension);

    if(parray == NULL)
    {
        throw std::runtime_error("Contiguous array couldn't be created from input python object");
    }


    // Fill image size, leave out the (non-scalar) data dimensions
    SizeType size;
    for(unsigned int d = 0; d<TImage::GetImageDimension(); d++)
    {
        size[TImage::GetImageDimension() - d - 1] = parray->dimensions[d];
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

    unsigned long numberOfPixels = 1;
    for(unsigned int d=0; d<arrayDimension; ++d)
    {
        numberOfPixels *= parray->dimensions[d];
    }

    typedef itk::ImportImageContainer<unsigned long, IOPixelType> ContainerType;
    typename ContainerType::Pointer container = ContainerType::New();
    IOPixelType* data = (IOPixelType*) parray->data;
    container->SetImportPointer(data, numberOfPixels, transferOwnership);

    result->SetRegions(region);
    result->SetOrigin(origin);
    result->SetSpacing(spacing);
    result->SetPixelContainer(container);
    if(result->GetNameOfClass() == std::string("VectorImage"))
    {
        result->SetNumberOfComponentsPerPixel(parray->dimensions[arrayDimension-1]);
    }

    if(transferOwnership)
    {
        Self::set_array_ownership(obj, false);
    }

    // has been Py_INCREF'ed in PyArray_FromArray, called by PyArray_FromAny,
    // called by PyArray_ContiguousFromObject
    Py_XDECREF(obj);

    return result;
}


template<typename TImage>
bool
NumpyBridge<TImage>
::IsBufferShared(PyObject* obj, ImageType* image)
{
    if(!PyArray_ISCONTIGUOUS(obj))
   {
       throw std::runtime_error("Cannot convert a non C-contiguous array");
   }

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
        throw std::runtime_error(std::string("Cannot convert object of type") + image->GetNameOfClass());
    }

    // Get the array object
    int element_type = GetPyType();
    PyArrayObject * array = (PyArrayObject *) PyArray_ContiguousFromAny(obj, element_type, arrayDimension, arrayDimension);

    // has been Py_INCREF'ed in PyArray_FromArray, called by PyArray_FromAny,
    // called by PyArray_ContiguousFromObject
    Py_XDECREF(obj);

    return (array->data==reinterpret_cast<char*>(image->GetBufferPointer()));
}


template<typename TImage>
typename NumpyBridge<TImage>::PyArrayType
NumpyBridge<TImage>
::GetPyType(void)
{

  PyArrayType item_type;
  typedef typename PixelTraits< IOPixelType >::ValueType    ScalarType;
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

template<typename TImage>
void
NumpyBridge<TImage>
::set_array_ownership(PyObject* array, bool value)
{
    PyArrayObject* pyArray = reinterpret_cast<PyArrayObject*>(array);

    if(value)
    {
        pyArray->flags |= NPY_OWNDATA;
    }
    else
    {
        pyArray->flags &= ~NPY_OWNDATA;
    }
}

} // namespace itk

#endif // itk_itknumpybridge_txx
