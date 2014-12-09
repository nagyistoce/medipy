/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef dca70582_dd05_40fe_aed4_0d84150ba81a
#define dca70582_dd05_40fe_aed4_0d84150ba81a

#include <itkImageToImageFilter.h>

namespace itk
{

/**
 * @brief Mono-modal statistical change detection.
 *
 * Implementation of the method described in "Automatic Chance Detection in
 * Multi-Modal Serial {MRI}: Application to Multiple Sclerosis Lesion Follow-up"
 * by Bosc & al. (NeuroImage, 2(20), pp. 643-656, Oct. 2003).
 */
template<typename TInputImage, typename TMaskImage, typename TOutputImage>
class MonoModalStatisticalChangeDetectionImageFilter :
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public :
    /** Standard class typedefs. */
    typedef MonoModalStatisticalChangeDetectionImageFilter Self;
    typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<Self const> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(MonoModalStatisticalChangeDetectionImageFilter, ImageToImageFilter);

    /** Type of the input image. */
    typedef TInputImage InputImageType;

    /** Type of the output image. */
    typedef TOutputImage OutputImageType;

    /** Type of the mask image. */
    typedef TMaskImage MaskImageType;

    void SetInput1(InputImageType const * image);
    void SetInput2(InputImageType const * image);

    /** Mask image. */
    itkGetConstObjectMacro(Mask, MaskImageType);
    itkSetObjectMacro(Mask, MaskImageType);

    /** Background value of the mask, default to 0. */
    itkGetMacro(MaskBackgroundValue, typename MaskImageType::PixelType);
    itkSetMacro(MaskBackgroundValue, typename MaskImageType::PixelType);

protected :
    MonoModalStatisticalChangeDetectionImageFilter();
    ~MonoModalStatisticalChangeDetectionImageFilter();
    void PrintSelf(std::ostream& os, itk::Indent indent) const;

    void BeforeThreadedGenerateData();
    void ThreadedGenerateData(typename InputImageType::RegionType const & region, ThreadIdType);

private :

    typedef typename NumericTraits<typename InputImageType::PixelType>::RealType
        InputRealType;
    typedef itk::Image<InputRealType, InputImageType::ImageDimension>
        DifferenceImageType;
    typedef itk::Neighborhood<typename DifferenceImageType::PixelType,
        DifferenceImageType::ImageDimension> DifferenceNeighborhoodType;

    typename MaskImageType::Pointer m_Mask;
    typename MaskImageType::PixelType m_MaskBackgroundValue;

    InputRealType variance_;
    typename DifferenceImageType::Pointer difference_image_;

    MonoModalStatisticalChangeDetectionImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented
};

}

#include "itkMonoModalStatisticalChangeDetectionImageFilter.txx"

#endif // dca70582_dd05_40fe_aed4_0d84150ba81a
