/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _3bb18cbd_6c72_44be_953e_9c713b891872
#define _3bb18cbd_6c72_44be_953e_9c713b891872

#include <Python.h>

#include <itkImageIOBase.h>
#include <itkProcessObject.h>

#include <numpy/arrayobject.h>

namespace itk
{

/** \brief Data source that reads image data from a single file, and stores
 *  image information (data, origin, spacing and direction) in NumPy arrays.
 *
 *  \warning This is NOT a ProcessObject
 */
template<typename TPixel, unsigned int VImageDimension>
class PyArrayFileReader : public Object
{
public :
	/** Standard class typedefs. */
	typedef PyArrayFileReader Self;
	typedef Object Superclass;
	typedef SmartPointer<Self> Pointer;

	/** Method for creation through the object factory. */
	itkNewMacro(Self);

	/** Run-time type information (and related methods). */
	itkTypeMacro(PyArrayFileReader, Object);

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
    PyObject* GetArray() const;

	/** Return the image origin. */
	PyObject* GetOrigin() const;

	/** Return the image spacing. */
	PyObject* GetSpacing() const;

	/** Return the image direction. */
	PyObject* GetDirection() const;

protected:
	typedef enum PyArray_TYPES PyArrayType;

	std::string m_FileName;

	ImageIOBase::Pointer m_ImageIO;
	bool m_UserSpecifiedImageIO;

	PyObject* m_Array;
	PyObject* m_Origin;
	PyObject* m_Spacing;
	PyObject* m_Direction;

	PyArrayFileReader();
	~PyArrayFileReader();

	void PrintSelf(std::ostream& os, Indent indent) const;

	/** Return a NumPy array type from a C type. */
	template<typename T>
	static PyArrayType GetPyType(void);

	/** Generate the array, origin, spacing, and direction. */
	template<typename TImage>
	void GenerateData();

	/** Set the image data from an ITK image. */
	template<typename TImage>
	void SetArray(TImage * image);

	/** Set the origin from an ITK image. */
	template<typename TImage>
	void SetOrigin(TImage * image);

	/** Set the spacing from an ITK image. */
	template<typename TImage>
	void SetSpacing(TImage * image);

	/** Set the direction from an ITK image. */
	template<typename TImage>
	void SetDirection(TImage * image);

private:
	PyArrayFileReader(const Self&); //purposely not implemented
	void operator=(const Self&); //purposely not implemented
};

}

#include "itkPyArrayFileReader.txx"

#endif // _3bb18cbd_6c72_44be_953e_9c713b891872
