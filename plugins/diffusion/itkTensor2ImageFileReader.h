/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _e33eb69a_0d02_4e28_885f_614b5dbbff57
#define _e33eb69a_0d02_4e28_885f_614b5dbbff57

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageSource.h>
#include <itkSymmetricSecondRankTensor.h>

namespace itk
{

template<typename TOutputImage>
class Tensor2ImageFileReader : public ImageSource<TOutputImage>
{
public :
    /** Standard class typedefs. */
    typedef Tensor2ImageFileReader Self;
    typedef ImageSource<TOutputImage> Superclass;
    typedef SmartPointer<Self> Pointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(Tensor2ImageFileReader, ImageSource);

    /** The size of the output image. */
    typedef typename TOutputImage::SizeType SizeType;

    /** The size of the output image. */
    typedef typename TOutputImage::IndexType IndexType;

    /** The region of the output image. */
    typedef typename TOutputImage::RegionType ImageRegionType;

    /** The pixel type of the output image. */
    typedef typename TOutputImage::PixelType OutputImagePixelType;

    /** Specify the file to read. This is forwarded to the IO instance. */
    virtual void SetFileName(std::string const & filename);
    virtual std::string const & GetFileName() const;

    /** Set/Get the ImageIO helper class. Often this is created via the object
     * factory mechanism that determines whether a particular ImageIO can
     * read a certain file. This method provides a way to get the ImageIO
     * instance that is created. Or you can directly specify the ImageIO
     * to use to read a particular file in case the factory mechanism will
     * not work properly (e.g., unknown or unusual extension). */
    void SetImageIO(ImageIOBase* imageIO);
    virtual ImageIOBase* GetImageIO();

    /** Prepare the allocation of the output image during the first back
     * propagation of the pipeline. */
    virtual void GenerateOutputInformation();

    /** Give the reader a chance to indicate that it will produce more
     * output than it was requested to produce. ImageFileReader cannot
     * currently read a portion of an image (since the ImageIO objects
     * cannot read a portion of an image), so the ImageFileReader must
     * enlarge the RequestedRegion to the size of the image on disk. */
    virtual void EnlargeOutputRequestedRegion(DataObject* output);

protected :
    Tensor2ImageFileReader();
    ~Tensor2ImageFileReader();
    void PrintSelf(std::ostream& os, Indent indent) const;

    /** Does the real work. */
    virtual void GenerateData();

private :
    typedef TOutputImage OutputImageType;
    typedef typename OutputImageType::Pointer OutputImagePointer;

    typedef typename TOutputImage::PixelType PixelType;
    typedef itk::SymmetricSecondRankTensor<typename PixelType::ValueType, 3> TensorType;

    typedef itk::Image<TensorType, TOutputImage::ImageDimension> TensorImageType;
    typedef typename TensorImageType::Pointer TensorImagePointer;

    typedef ImageFileReader<TensorImageType> TensorReader;

    typename TensorReader::Pointer m_TensorReader;

    Tensor2ImageFileReader(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    void tensorToUpperDiagonal(TensorImagePointer tensors,
                               OutputImagePointer output) const;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTensor2ImageFileReader.txx"
#endif

#endif // _e33eb69a_0d02_4e28_885f_614b5dbbff57
