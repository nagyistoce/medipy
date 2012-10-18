/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/


#ifndef _itkFractionalAnisotropyImageFilter_h
#define _itkFractionalAnisotropyImageFilter_h


#include <itkImageToImageFilter.h>
#include "itkSymmetricSpectralAnalysisImageFilter.h"
#include <itkSmartPointer.h>


namespace itk
{

/**
 * \class FractionalAnisotropyImageFilter
 * \brief Compute the fractional anisotropy of a second order diffusion tensor image
 * 
 */

template<typename TInputImage, typename TOutputImage>
class FractionalAnisotropyImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public :
    /** Standard class typedefs. */
    typedef FractionalAnisotropyImageFilter Self;
    typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<Self const> ConstPointer;

    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(FractionalAnisotropyImageFilter, ImageToImageFilter);

    /** Useful typedefs */
    typedef typename Superclass::InputImageType     InputImageType;
    typedef typename Superclass::OutputImageType    OutputImageType;
    typedef typename TOutputImage::PixelType        OutputPixelType;
    typedef typename TInputImage::PixelType         InputPixelType;
    typedef typename InputPixelType::ValueType      InputValueType;

    typedef typename Superclass::InputImageType EigenValueImageType;
    typedef typename Superclass::InputImageType EigenVectorImageType;

    /** Accessors */
    void SetEigenValue(EigenValueImageType *val);
    void SetEigenVector(EigenVectorImageType *vec);
    EigenValueImageType* GetEigenValue();
    EigenVectorImageType* GetEigenVector();
    itkStaticConstMacro(OutputImageDimension, unsigned int, TOutputImage::ImageDimension);

protected :
    FractionalAnisotropyImageFilter();
    ~FractionalAnisotropyImageFilter() {}
    void PrintSelf(std::ostream& os, Indent indent) const;
    void BeforeThreadedGenerateData();
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, int);

private :
    bool m_EigenSystem;
    typename EigenValueImageType::Pointer m_EigVal;
    typename EigenVectorImageType::Pointer m_EigVec;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFractionalAnisotropyImageFilter.txx"
#endif

#endif 

