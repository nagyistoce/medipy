/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2012
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _f600a6dd_4a7f_4d0e_8f76_80e092ee1118
#define _f600a6dd_4a7f_4d0e_8f76_80e092ee1118

#include "itkPyArrayFileReader.h"

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageIOFactory.h>
#include <itkImageRegionConstIterator.h>
#include <itkImportImageContainer.h>
#include <itkRGBPixel.h>

namespace itk
{

/**
 * @brief Component type of a pixel.
 *
 * Generic class, will be specialized for non-scalar pixels.
 */
template<typename TPixel>
struct ComponentTrait { typedef TPixel Type; };

/**
 * @brief Component type of a pixel.
 *
 * Specialization for itk::RGBPixel
 */
template<typename TComponentType>
struct ComponentTrait<itk::RGBPixel<TComponentType> > { typedef TComponentType Type; };

/**
 * @brief Functor filling an array from an itk container.
 *
 * Generic class, will be specialized for non-scalar pixels.
 */
template<typename TElementIdentifier, typename TElement>
struct FillArray
{
    typedef typename ComponentTrait<TElement>::Type ComponentType;
    typedef itk::ImportImageContainer<TElementIdentifier, TElement> ContainerType;
    
    static void
    action(ContainerType * container, PyObject * array)
    {
        TElement * source = container->GetBufferPointer();
        TElement * end = source + container->Size();
        ComponentType * destination = 
            reinterpret_cast<ComponentType*>(PyArray_DATA(array));
        std::copy(source, end, destination);
    }
};

/**
 * @brief Functor filling an array from an itk container.
 *
 * Specialization for itk::RGBPixel
 */
template<typename TElementIdentifier, typename TComponent>
struct FillArray<TElementIdentifier, itk::RGBPixel<TComponent> >
{
    typedef itk::RGBPixel<TComponent> TElement;
    
    typedef typename ComponentTrait<TElement>::Type ComponentType;
    typedef itk::ImportImageContainer<TElementIdentifier, TElement> ContainerType;
    
    static void
    action(ContainerType * container, PyObject * array)
    {
        TElement * source = container->GetImportPointer();
        TElement * end = source + container->Capacity();
        ComponentType * destination = 
            reinterpret_cast<ComponentType*>(PyArray_DATA(array));
        
        for(/* No initialization */; source != end; ++source)
        {
            destination = std::copy(source->Begin(), source->End(), destination);
            ++source;
        }
    }
};


template<typename TPixel, unsigned int VImageDimension>
void
PyArrayFileReader<TPixel, VImageDimension>
::SetImageIO(ImageIOBase* imageIO)
{
	itkDebugMacro("setting ImageIO to " << imageIO );
	if(this->m_ImageIO != imageIO )
	{
		this->m_ImageIO = imageIO;
		this->Modified();
	}
	m_UserSpecifiedImageIO = true;
}

template<typename TPixel, unsigned int VImageDimension>
void
PyArrayFileReader<TPixel, VImageDimension>
::Update()
{
    if(!this->m_UserSpecifiedImageIO)
    {
        this->m_ImageIO = itk::ImageIOFactory::CreateImageIO(
            this->GetFileName(), itk::ImageIOFactory::ReadMode);
    }

    this->m_ImageIO->SetFileName(this->GetFileName());
    this->m_ImageIO->ReadImageInformation();

    if(this->m_ImageIO->GetPixelType() == ImageIOBase::SCALAR)
    {
        typedef Image<TPixel, VImageDimension> ImageType;
        this->template GenerateData<ImageType>();
    }
    else if(this->m_ImageIO->GetPixelType() == ImageIOBase::RGB)
    {
        typedef Image<RGBPixel<TPixel>, VImageDimension> ImageType;
        this->template GenerateData<ImageType>();
    }
    else
    {
        throw std::logic_error("Cannot process images with PixelType=="+
            this->m_ImageIO->GetComponentTypeAsString(this->m_ImageIO->GetComponentType()));
    }
}

template<typename TPixel, unsigned int VImageDimension>
PyObject*
PyArrayFileReader<TPixel, VImageDimension>
::GetArray() const
{
    Py_XINCREF(this->m_Array);
    return this->m_Array;
}

template<typename TPixel, unsigned int VImageDimension>
PyObject*
PyArrayFileReader<TPixel, VImageDimension>
::GetOrigin() const
{
    Py_XINCREF(this->m_Origin);
    return this->m_Origin;
}

template<typename TPixel, unsigned int VImageDimension>
PyObject*
PyArrayFileReader<TPixel, VImageDimension>
::GetSpacing() const
{
    Py_XINCREF(this->m_Spacing);
    return this->m_Spacing;
}

template<typename TPixel, unsigned int VImageDimension>
PyObject*
PyArrayFileReader<TPixel, VImageDimension>
::GetDirection() const
{
    Py_XINCREF(this->m_Direction);
    return this->m_Direction;
}

template<typename TPixel, unsigned int VImageDimension>
PyArrayFileReader<TPixel, VImageDimension>
::PyArrayFileReader()
: m_FileName(""), m_ImageIO(NULL), m_UserSpecifiedImageIO(false),
  m_Array(NULL), m_Origin(NULL), m_Spacing(NULL), m_Direction(NULL)
{
	// Nothing else
}

template<typename TPixel, unsigned int VImageDimension>
PyArrayFileReader<TPixel, VImageDimension>
::~PyArrayFileReader()
{
    Py_XDECREF(this->m_Array);
    Py_XDECREF(this->m_Origin);
    Py_XDECREF(this->m_Spacing);
    Py_XDECREF(this->m_Direction);
}

template<typename TPixel, unsigned int VImageDimension>
void
PyArrayFileReader<TPixel, VImageDimension>
::PrintSelf(std::ostream& os, Indent indent) const
{
	Superclass::PrintSelf(os, indent);

	if(m_ImageIO)
	{
		os << indent << "ImageIO: \n";
		m_ImageIO->Print(os, indent.GetNextIndent());
	}
	else
	{
		os << indent << "ImageIO: (null)" << "\n";
	}

	os << indent << "UserSpecifiedImageIO flag: " << m_UserSpecifiedImageIO << "\n";
	os << indent << "FileName: " << m_FileName << "\n";
}

template<typename TPixel, unsigned int VImageDimension>
template<typename T>
typename PyArrayFileReader<TPixel, VImageDimension>::PyArrayType
PyArrayFileReader<TPixel, VImageDimension>
::GetPyType(void)
{

	PyArrayType item_type;
	typedef typename PixelTraits<T>::ValueType ScalarType;
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

template<typename TPixel, unsigned int VImageDimension>
template<typename TImage>
void
PyArrayFileReader<TPixel, VImageDimension>
::GenerateData()
{
    typedef ImageFileReader<TImage> ReaderType;

    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetImageIO(this->m_ImageIO);
    reader->SetFileName(this->m_ImageIO->GetFileName());
    reader->Update();

    typename TImage::Pointer image = reader->GetOutput();

    // Get geometric informations from image
    this->template SetOrigin<TImage>(image.GetPointer());
    this->template SetSpacing<TImage>(image.GetPointer());
    this->template SetDirection<TImage>(image.GetPointer());

    this->template SetArray<TImage>(image.GetPointer());
}

template<typename TPixel, unsigned int VImageDimension>
template<typename TImage>
void
PyArrayFileReader<TPixel, VImageDimension>
::SetArray(TImage * image)
{
    Py_XDECREF(this->m_Array);
    
    // Get the number of dimensions
    unsigned int ndim;
    if(this->m_ImageIO->GetPixelType() == ImageIOBase::SCALAR)
    {
        ndim = VImageDimension;
    }
    else if(this->m_ImageIO->GetPixelType() == ImageIOBase::RGB)
    {
        ndim = 1+VImageDimension;
    }
    else
    {
        throw std::logic_error("Cannot assign array with PixelType=="+
            this->m_ImageIO->GetComponentTypeAsString(this->m_ImageIO->GetComponentType()));
    }

    // Get the array size in NumPy order
    npy_intp * array_size = new npy_intp[ndim];
    typename TImage::SizeValueType const * source =
        image->GetRequestedRegion().GetSize().m_Size;
    std::reverse_copy(source, source+VImageDimension, array_size);
    if(this->m_ImageIO->GetPixelType() == ImageIOBase::RGB)
    {
        array_size[VImageDimension]=3;
    }

	// Create the array
    typedef typename TImage::PixelType PixelType;
	typedef typename ComponentTrait<PixelType>::Type ComponentType;
	PyArrayType const pyArrayType = Self::template GetPyType<ComponentType>();

    this->m_Array = PyArray_SimpleNew(ndim, array_size, pyArrayType);
    delete[] array_size;
    //reinterpret_cast<PyArrayObject*>(this->m_Array)->flags |= NPY_OWNDATA;
    
    // Get the data
    typedef typename TImage::PixelContainer ContainerType;
    ContainerType * container = image->GetPixelContainer();
    typedef FillArray<
        typename ContainerType::ElementIdentifier, 
        typename TImage::PixelType> FillArrayType;
    FillArrayType::action(container, this->m_Array);
}

template<typename TPixel, unsigned int VImageDimension>
template<typename TImage>
void
PyArrayFileReader<TPixel, VImageDimension>
::SetOrigin(TImage * image)
{
    Py_XDECREF(this->m_Origin);
    
	typedef typename TImage::PointType::ValueType ValueType;

	// Create the array
	npy_intp size = VImageDimension;
	this->m_Origin = PyArray_SimpleNew(1, &size, Self::template GetPyType<ValueType>());

	// Copy from ITK object, reversing the order to get the NumPy order
	std::reverse_copy(image->GetOrigin().GetDataPointer(),
		image->GetOrigin().GetDataPointer()+VImageDimension,
		reinterpret_cast<ValueType*>(PyArray_DATA(this->m_Origin)));
}

template<typename TPixel, unsigned int VImageDimension>
template<typename TImage>
void
PyArrayFileReader<TPixel, VImageDimension>
::SetSpacing(TImage * image)
{
    Py_XDECREF(this->m_Spacing);
	typedef typename TImage::SpacingType::ValueType ValueType;

	// Create the array
	npy_intp size = VImageDimension;
	this->m_Spacing = PyArray_SimpleNew(1, &size, Self::template GetPyType<ValueType>());

	// Copy from ITK object, reversing the order to get the NumPy order
	std::reverse_copy(image->GetSpacing().GetDataPointer(),
		image->GetSpacing().GetDataPointer()+VImageDimension,
		reinterpret_cast<ValueType*>(PyArray_DATA(this->m_Spacing)));
}

template<typename TPixel, unsigned int VImageDimension>
template<typename TImage>
void
PyArrayFileReader<TPixel, VImageDimension>
::SetDirection(TImage * image)
{
    Py_XDECREF(this->m_Direction);
    
	typedef typename TImage::DirectionType::ValueType ValueType;
	// Create the array
	npy_intp size[] = {VImageDimension, VImageDimension};
	this->m_Direction= PyArray_SimpleNew(2, size, Self::template GetPyType<ValueType>());

	// Copy from ITK object, reversing the order to get the NumPy order
	std::reverse_copy(&image->GetDirection()(0,0),
			&image->GetDirection()(0,0)+VImageDimension*VImageDimension,
		reinterpret_cast<ValueType*>(PyArray_DATA(this->m_Direction)));
}

}

#endif // _f600a6dd_4a7f_4d0e_8f76_80e092ee1118
