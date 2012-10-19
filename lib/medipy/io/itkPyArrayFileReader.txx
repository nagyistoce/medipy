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

namespace itk
{

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
	typedef Image<TPixel, VImageDimension> ImageType;
	typedef ImageFileReader<ImageType> ReaderType;

	typename ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(this->GetFileName());
	if(this->m_UserSpecifiedImageIO)
	{
		reader->SetImageIO(this->GetImageIO());
	}

	reader->Update();

	this->m_ImageIO = reader->GetImageIO();

	typename ImageType::Pointer image = reader->GetOutput();

	// Get geometric informations from image
	this->template SetOrigin(image.GetPointer());
	this->template SetSpacing(image.GetPointer());
	this->template SetDirection(image.GetPointer());

	// Get data from image, transfer data ownership
	this->template SetArray(image.GetPointer());
	reinterpret_cast<PyArrayObject*>(this->m_Array)->flags |= NPY_OWNDATA;
	image->GetPixelContainer()->ContainerManageMemoryOff();
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
	// Nothing else
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
::SetArray(TImage * image)
{
	// Get the array size in NumPy order
	typename TImage::SizeType const region_size =
		image->GetRequestedRegion().GetSize();
	npy_intp array_size[VImageDimension];
	std::reverse_copy(region_size.m_Size, region_size.m_Size+VImageDimension,
					  array_size);

	// Create the array
	this->m_Array = PyArray_SimpleNewFromData(VImageDimension, array_size,
		Self::template GetPyType<typename TImage::PixelType>(),
		reinterpret_cast<void*>(image->GetBufferPointer()));
}

template<typename TPixel, unsigned int VImageDimension>
template<typename TImage>
void
PyArrayFileReader<TPixel, VImageDimension>
::SetOrigin(TImage * image)
{
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
