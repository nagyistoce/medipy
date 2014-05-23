/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2012
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _7e43e11a_30b4_4457_89c8_549a2ebde293
#define _7e43e11a_30b4_4457_89c8_549a2ebde293

#include <Python.h>

#include <itkImage.h>
#include <itkImageIOBase.h>
#include <itkProcessObject.h>

#include <numpy/arrayobject.h>

namespace itk
{

/** \brief Data source that writes image data to a single file, using image
 *  information (data, origin, spacing and direction) stored in NumPy arrays.
 *
 *  \warning This is NOT a ProcessObject
 */
template<typename TPixel, unsigned int VImageDimension>
class PyArrayFileWriter : public Object
{
public :
	/** Standard class typedefs. */
	typedef PyArrayFileWriter Self;
	typedef Object Superclass;
	typedef SmartPointer<Self> Pointer;

	/** Method for creation through the object factory. */
	itkNewMacro(Self);

	/** Run-time type information (and related methods). */
	itkTypeMacro(PyArrayFileWriter, Object);

	/** Specify the file to read. This is forwarded to the IO instance. */
	itkSetStringMacro(FileName);
	itkGetStringMacro(FileName);

	/** Set/Get the ImageIO helper class. Often this is created via the object
	 * factory mechanism that determines whether a particular ImageIO can
	 * read a certain file. This method provides a way to get the ImageIO
	 * instance that is created. Or you can directly specify the ImageIO
	 * to use to read a particular file in case the factory mechanism will
	 * not work properly (e.g., unknown or unusual extension). */
	void SetImageIO(ImageIOBase* imageIO);
	itkGetObjectMacro(ImageIO, ImageIOBase);

	/** Read the image data. */
	void Update();

	/** Return the image data. */
	itkGetConstMacro(Array, PyObject*);
	void SetArray(PyObject* value);

	/** Return the image origin. */
	itkGetConstMacro(Origin, PyObject*);
	void SetOrigin(PyObject* value);

	/** Return the image spacing. */
	itkGetConstMacro(Spacing, PyObject*);
	void SetSpacing(PyObject* value);

	/** Return the image direction. */
	itkGetConstMacro(Direction, PyObject*);
	void SetDirection(PyObject* value);

protected:
	typedef enum PyArray_TYPES PyArrayType;
	typedef Image<TPixel, VImageDimension> ImageType;

	std::string m_FileName;

	ImageIOBase::Pointer m_ImageIO;
	bool m_UserSpecifiedImageIO;

	PyObject* m_Array;
	PyObject* m_Origin;
	PyObject* m_Spacing;
	PyObject* m_Direction;

	typename ImageType::Pointer m_Image;

	PyArrayFileWriter();
	~PyArrayFileWriter();

	void PrintSelf(std::ostream& os, Indent indent) const;

	template<typename TOutputIterator>
	static void copyPyArray(PyObject* array,
							TOutputIterator const & destination,
							bool reverse=false);

	template<typename TArrayValue, typename TOutputIterator>
	static void copyPyArrayInner(PyObject* array,
							     TOutputIterator const & destination,
							     bool reverse=false);

	/** Set the image data from an ITK image. */
	template<typename TImage>
	void SetArray(PyObject* array);

	/** Set the origin from an ITK image. */
	template<typename TImage>
	void SetOrigin(PyObject* origin);

	/** Set the spacing from an ITK image. */
	template<typename TImage>
	void SetSpacing(PyObject* origin);

	/** Set the direction from an ITK image. */
	template<typename TImage>
	void SetDirection(PyObject* origin);

private:
	PyArrayFileWriter(const Self&); //purposely not implemented
	void operator=(const Self&); //purposely not implemented
};

}

#include "itkPyArrayFileWriter.txx"

#endif // _7e43e11a_30b4_4457_89c8_549a2ebde293
