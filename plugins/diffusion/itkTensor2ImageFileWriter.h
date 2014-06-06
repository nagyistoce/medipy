/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _69a7dc45_975d_4e55_bd73_7c699d23c56a
#define _69a7dc45_975d_4e55_bd73_7c699d23c56a

#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkImageIOBase.h>
#include <itkProcessObject.h>
#include <itkSymmetricSecondRankTensor.h>

namespace itk
{

template<typename TInputImage>
class Tensor2ImageFileWriter : public ProcessObject
{
public:
    /** Standard class typedefs. */
    typedef Tensor2ImageFileWriter Self;
    typedef ProcessObject Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<const Self> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(Tensor2ImageFileWriter, ProcessObject);

    /** Some convenient typedefs. */
    typedef TInputImage InputImageType;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename InputImageType::RegionType InputImageRegionType;
    typedef typename InputImageType::PixelType InputImagePixelType;

    /** Set/Get the image input of this writer.  */
    void SetInput(InputImageType const * input);
    InputImageType const * GetInput();
    InputImageType const * GetInput(unsigned int idx);

    /** Specify the name of the output file to write. */
    virtual void SetFileName(char const* filename);
    virtual void SetFileName(std::string const & filename);
    virtual const char* GetFileName() const;

    /** Set/Get the ImageIO helper class. Usually this is created via the object
     * factory mechanism that determines whether a particular ImageIO can
     * write a certain file. This method provides a way to get the ImageIO
     * instance that is created, or one can be manually set where the
     * IO factory mechanism may not work (for example, raw image files or
     * image files with non-standard filename suffix's.
     * If the user specifies the ImageIO, we assume she makes the
     * correct choice and will allow a file to be created regardless of
     * the file extension. If the factory has set the ImageIO, the
     * extension must be supported by the specified ImageIO. */
    void SetImageIO(ImageIOBase* imageIO);
    virtual ImageIOBase* GetImageIO();

    /** A special version of the Update() method for writers.  It
     * invokes start and end events and handles releasing data. It
     * eventually calls GenerateData() which does the actual writing.
     * Note: the write method will write data specified by the
     * IORegion. If not set, then then the whole image is written.  Note
     * that the region will be cropped to fit the input image's
     * LargestPossibleRegion. */
    virtual void Write();

    /** Specify the region to write. If left NULL, then the whole image
     * is written. */
    void SetIORegion(ImageIORegion const & region);
    ImageIORegion const & GetIORegion() const;

    /** Aliased to the Write() method to be consistent with the rest of the
     * pipeline. */
    virtual void Update();

protected:
    Tensor2ImageFileWriter();
    ~Tensor2ImageFileWriter();

    void PrintSelf(std::ostream& os, Indent indent) const;

private :
    typedef typename TInputImage::ConstPointer InputImageConstPointer;
    typedef itk::SymmetricSecondRankTensor<typename InputImagePixelType::ValueType, 3> TensorType;

    typedef itk::Image<TensorType, TInputImage::ImageDimension> TensorImageType;
    typedef typename TensorImageType::Pointer TensorImagePointer;

    typedef ImageFileWriter<TensorImageType> TensorWriter;

    typename TensorWriter::Pointer m_TensorWriter;

    Tensor2ImageFileWriter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    void upperDiagonalToTensor(InputImageConstPointer input,
                               TensorImagePointer tensors) const;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTensor2ImageFileWriter.txx"
#endif

#endif // _69a7dc45_975d_4e55_bd73_7c699d23c56a

