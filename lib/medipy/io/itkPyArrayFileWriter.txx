/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2012
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _77446535_7eef_424b_8156_9c2df6a5a28f
#define _77446535_7eef_424b_8156_9c2df6a5a28f

#include "itkPyArrayFileWriter.h"

#include <itkImage.h>
#include <itkImageFileWriter.h>

namespace itk
{

template<typename TPixel, unsigned int VImageDimension>
void
PyArrayFileWriter<TPixel, VImageDimension>
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
PyArrayFileWriter<TPixel, VImageDimension>
::Update()
{
	this->m_Image->Print(std::cout);

	typedef ImageFileWriter<ImageType> WriterType;

	typename WriterType::Pointer writer = WriterType::New();
	writer->SetInput(this->m_Image);
	writer->SetFileName(this->GetFileName());
	if(this->m_UserSpecifiedImageIO)
	{
		writer->SetImageIO(this->GetImageIO());
	}

	writer->Update();

	this->m_ImageIO = writer->GetImageIO();
}

template<typename TPixel, unsigned int VImageDimension>
template<typename TOutputIterator>
void
PyArrayFileWriter<TPixel, VImageDimension>
::copyPyArray(PyObject* array,
              TOutputIterator const & destination,
              bool reverse)
{
    int const array_type = PyArray_TYPE(array);

    if(array_type == NPY_BYTE)
    {
        return copyPyArrayInner<signed char, TOutputIterator>(array, destination, reverse);
    }
    else if(array_type == NPY_UBYTE)
    {
        return copyPyArrayInner<unsigned char, TOutputIterator>(array, destination, reverse);
    }
    else if(array_type == NPY_SHORT)
    {
        return copyPyArrayInner<signed short, TOutputIterator>(array, destination, reverse);
    }
    else if(array_type == NPY_USHORT)
    {
        return copyPyArrayInner<unsigned short, TOutputIterator>(array, destination, reverse);
    }
    else if(array_type == NPY_INT)
    {
        return copyPyArrayInner<signed int, TOutputIterator>(array, destination, reverse);
    }
    else if(array_type == NPY_UINT)
    {
        return copyPyArrayInner<unsigned int, TOutputIterator>(array, destination, reverse);
    }
    else if(array_type == NPY_LONG)
    {
        return copyPyArrayInner<signed long, TOutputIterator>(array, destination, reverse);
    }
    else if(array_type == NPY_ULONG)
    {
        return copyPyArrayInner<unsigned long, TOutputIterator>(array, destination, reverse);
    }
    else if(array_type == NPY_LONGLONG)
    {
        return copyPyArrayInner<signed long long, TOutputIterator>(array, destination, reverse);
    }
    else if(array_type == NPY_ULONGLONG)
    {
        return copyPyArrayInner<unsigned long long, TOutputIterator>(array, destination, reverse);
    }
    else if(array_type == NPY_FLOAT)
    {
        return copyPyArrayInner<float, TOutputIterator>(array, destination, reverse);
    }
    else
    {
        std::stringstream message;
        message << "Cannot operate on PyArray of type " << array_type;
        throw std::runtime_error(message.str());
    }
}

template<typename TPixel, unsigned int VImageDimension>
template<typename TArrayValue, typename TOutputIterator>
void
PyArrayFileWriter<TPixel, VImageDimension>
::copyPyArrayInner(PyObject* array,
              TOutputIterator const & destination,
              bool reverse)
{
    TArrayValue* begin = reinterpret_cast<TArrayValue*>(PyArray_DATA(array));
    TArrayValue* end = begin+PyArray_Size(array);
    if(reverse)
    {
        std::reverse_copy(begin, end, destination);
    }
    else
    {
        std::copy(begin, end, destination);
    }
}

template<typename TPixel, unsigned int VImageDimension>
void
PyArrayFileWriter<TPixel, VImageDimension>
::SetArray(PyObject* value)
{
	if(value != this->m_Array)
	{
		this->m_Array = value;

		typename ImageType::IndexType index; index.Fill(0);
		typename ImageType::SizeType size;
		std::reverse_copy(PyArray_DIMS(value),
						  PyArray_DIMS(value)+VImageDimension,
						  size.m_Size);
		typename ImageType::RegionType region(index, size);
		this->m_Image->SetRegions(region);

		typedef typename ImageType::PixelContainer ContainerType;
		typename ContainerType::Pointer container = ContainerType::New();
		container->SetImportPointer(
			reinterpret_cast<TPixel*>(PyArray_DATA(value)), PyArray_Size(value), false);

		this->m_Image->SetPixelContainer(container);


		this->Modified();
	}
}

template<typename TPixel, unsigned int VImageDimension>
void
PyArrayFileWriter<TPixel, VImageDimension>
::SetOrigin(PyObject* value)
{
	if(value != this->m_Array)
	{
		this->m_Origin = value;

		typename ImageType::PointType origin;
		copyPyArray<typename ImageType::PointType::Iterator>(value, origin.Begin(), true);
		this->m_Image->SetOrigin(origin);

		this->Modified();
	}

}

template<typename TPixel, unsigned int VImageDimension>
void
PyArrayFileWriter<TPixel, VImageDimension>
::SetSpacing(PyObject* value)
{
	if(value != this->m_Array)
	{
		this->m_Spacing = value;

		typename ImageType::SpacingType spacing;
		copyPyArray<typename ImageType::SpacingType::Iterator>(value, spacing.Begin(), true);
		this->m_Image->SetSpacing(spacing);

		this->Modified();
	}
}

template<typename TPixel, unsigned int VImageDimension>
void
PyArrayFileWriter<TPixel, VImageDimension>
::SetDirection(PyObject* value)
{
	if(value != this->m_Array)
	{
		this->m_Direction = value;

		typename ImageType::DirectionType direction;
		copyPyArray<typename ImageType::DirectionType::ValueType*>(value, &direction(0,0), true);
		this->m_Image->SetDirection(direction);

		this->Modified();
	}

}

template<typename TPixel, unsigned int VImageDimension>
PyArrayFileWriter<TPixel, VImageDimension>
::PyArrayFileWriter()
: m_FileName(""), m_ImageIO(NULL), m_UserSpecifiedImageIO(false),
  m_Array(NULL), m_Origin(NULL), m_Spacing(NULL), m_Direction(NULL),
  m_Image(ImageType::New())
{
	// Nothing else
}

template<typename TPixel, unsigned int VImageDimension>
PyArrayFileWriter<TPixel, VImageDimension>
::~PyArrayFileWriter()
{
	// Nothing else
}

template<typename TPixel, unsigned int VImageDimension>
void
PyArrayFileWriter<TPixel, VImageDimension>
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

}

#endif // _77446535_7eef_424b_8156_9c2df6a5a28f
