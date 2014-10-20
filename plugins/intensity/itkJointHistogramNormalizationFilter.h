/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _4bdca163_86fb_4f3f_a6bb_9d02a071a78e
#define _4bdca163_86fb_4f3f_a6bb_9d02a071a78e

#include "itkImageToImageFilter.h"
#include "itkJointHistogramCalculator.h"

namespace itk
{

template<typename TInputImage, typename TMask, typename TOutputImage>
class JointHistogramNormalizationFilter: 
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
    /** Standard typedefs */
    typedef JointHistogramNormalizationFilter Self;
    typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<Self const> ConstPointer;

    /** standard New() method support */
    itkNewMacro(Self);
    
    /** Run-time type information (and related methods). */
    itkTypeMacro(JointHistogramCalculator, ImageToImageFilter);
    
    typedef typename Superclass::InputImageType InputImageType;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename InputImageType::ConstPointer InputImageConstPointer;
    typedef typename Superclass::OutputImageType OutputImageType;
    
    typedef TMask MaskType;
    typedef typename MaskType::Pointer MaskPointer;
    typedef typename MaskType::PixelType MaskPixelType;
    
    typedef JointHistogramCalculator<TInputImage, TMask> 
        JointHistogramCalculatorType;
    typedef typename JointHistogramCalculatorType::Method Method;
    
    InputImageConstPointer GetFixedImage() { return this->GetInput(0); }
    void SetFixedImage(InputImageType const * image) { this->SetInput(0, image); }
    
    InputImageConstPointer GetMovingImage() { return this->GetInput(1); }
    void SetMovingImage(InputImageType const * image) { this->SetInput(1, image); }
    
    itkGetConstMacro(BinsCountFixed, unsigned int);
    itkSetMacro(BinsCountFixed, unsigned int);
    
    itkGetConstMacro(BinsCountMoving, unsigned int);
    itkSetMacro(BinsCountMoving, unsigned int);
    
    itkGetConstObjectMacro(Mask, MaskType);
    itkSetObjectMacro(Mask, MaskType);
    
    itkGetConstMacro(MaskValue, MaskPixelType);
    itkSetMacro(MaskValue, MaskPixelType);
    
    itkGetEnumMacro(Method, typename Method::Type);
    itkSetEnumMacro(Method, typename Method::Type);
    void SetMethodToNearestNeighbor() { this->SetMethod(Method::NEAREST_NEIGHBOR); }
    void SetMethodToLinearInterpolation() { this->SetMethod(Method::LINEAR_INTERPOLATION); }
    
protected:
    unsigned int m_BinsCountFixed;
    unsigned int m_BinsCountMoving;
    
    MaskPointer m_Mask;
    MaskPixelType m_MaskValue;
    typename Method::Type m_Method;
    
    JointHistogramNormalizationFilter();
    void PrintSelf(std::ostream & os, Indent indent);
    void GenerateData();
    
private:
    JointHistogramNormalizationFilter(Self const &); // purposely not implemented
    Self const & operator=(Self const &); // purposely not implemented
};

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkJointHistogramNormalizationFilter.txx"
#endif

#endif // _4bdca163_86fb_4f3f_a6bb_9d02a071a78e
