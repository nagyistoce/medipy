/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef dbb086f8_040b_43f5_9ce7_b015deb978fe
#define dbb086f8_040b_43f5_9ce7_b015deb978fe

#include "itkMonoModalStatisticalChangeDetectionImageFilter.h"

#include <cmath>

#include <itkSubtractImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkConstNeighborhoodIterator.h>
#include <itkNumericTraits.h>
#include <itkStatisticsImageFilter.h>

namespace itk
{

template<typename TInputImage, typename TMaskImage, typename TOutputImage>
void
MonoModalStatisticalChangeDetectionImageFilter<TInputImage, TMaskImage, TOutputImage>
::SetInput1(TInputImage const * image)
{
  // Process object is not const-correct so the const casting is required.
  this->SetNthInput(0, const_cast<TInputImage*>(image));
}

template<typename TInputImage, typename TMaskImage, typename TOutputImage>
void
MonoModalStatisticalChangeDetectionImageFilter<TInputImage, TMaskImage, TOutputImage>
::SetInput2(TInputImage const * image)
{
  // Process object is not const-correct so the const casting is required.
  this->SetNthInput(1, const_cast<TInputImage*>(image));
}

template<typename TInputImage, typename TMaskImage, typename TOutputImage>
MonoModalStatisticalChangeDetectionImageFilter<TInputImage, TMaskImage, TOutputImage>
::MonoModalStatisticalChangeDetectionImageFilter()
: m_Mask(0), m_MaskBackgroundValue(0), variance_(0), difference_image_(0)
{
}

template<typename TInputImage, typename TMaskImage, typename TOutputImage>
MonoModalStatisticalChangeDetectionImageFilter<TInputImage, TMaskImage, TOutputImage>
::~MonoModalStatisticalChangeDetectionImageFilter()
{
    // Nothing to do
}

template<typename TInputImage, typename TMaskImage, typename TOutputImage>
void
MonoModalStatisticalChangeDetectionImageFilter<TInputImage, TMaskImage, TOutputImage>
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
    this->Superclass::PrintSelf(os, indent);
    os << indent << "Mask : " << this->m_Mask << "\n";
    os << indent << "Mask background value : " << this->m_MaskBackgroundValue << "\n";
}

template<typename TInputImage, typename TMaskImage, typename TOutputImage>
void
MonoModalStatisticalChangeDetectionImageFilter<TInputImage, TMaskImage, TOutputImage>
::BeforeThreadedGenerateData()
{
    InputImageType const * input_1 = this->GetInput(0);
    InputImageType const * input_2 = this->GetInput(1);

    // Compute difference image
    typedef SubtractImageFilter<InputImageType, InputImageType, DifferenceImageType>
        SubtractFilterType;
    typename SubtractFilterType::Pointer subtract_filter = SubtractFilterType::New();
    subtract_filter->SetInput1(input_1);
    subtract_filter->SetInput2(input_2);
    subtract_filter->Update();
    this->difference_image_ = subtract_filter->GetOutput();

    // Compute variance on masked difference image, or on whole image if mask
    // has not been defined.
    if(this->m_Mask.IsNull())
    {
        itkWarningMacro(<<"No mask supplied, computing variance on the whole image");
        typedef itk::StatisticsImageFilter<DifferenceImageType> StatisticsFilterType;
        typename StatisticsFilterType::Pointer statistics = StatisticsFilterType::New();
        statistics->SetInput(this->difference_image_);
        statistics->Update();
        this->variance_ = statistics->GetVariance();
    }
    else
    {
        typedef itk::ImageRegionConstIterator<DifferenceImageType> InputIteratorType;
        typedef itk::ImageRegionConstIterator<MaskImageType> MaskIteratorType;

        // Make sure mask is up-to-date
        this->m_Mask->Update();

        InputIteratorType input_it(this->difference_image_, input_1->GetRequestedRegion());
        MaskIteratorType mask_it(this->m_Mask, input_1->GetRequestedRegion());

        unsigned long count=0;
        InputRealType sum = NumericTraits<InputRealType>::Zero;
        InputRealType sum_of_squares = NumericTraits<InputRealType>::Zero;

        while(!input_it.IsAtEnd())
        {
            if(mask_it.Get() != this->m_MaskBackgroundValue)
            {
                InputRealType const value(input_it.Get());
                sum += value;
                sum_of_squares += (value*value);
                ++count;
            }

            ++input_it;
            ++mask_it;
        }

        if(count != 0)
        {
            InputRealType const real_count(count);
            // mean = sum/real_count;
            this->variance_ = (sum_of_squares-sum*sum/real_count)/real_count;
        }
        else
        {
            itkWarningMacro(<<"No pixel included in the mask : setting variance to 0");
            this->variance_ = 0;
        }
    }

    itkDebugMacro(<< "Variance : " << this->variance_);

    this->GetOutput()->FillBuffer(0);
}

template<typename TInputImage, typename TMaskImage, typename TOutputImage>
void
MonoModalStatisticalChangeDetectionImageFilter<TInputImage, TMaskImage, TOutputImage>
::ThreadedGenerateData(typename InputImageType::RegionType const & region, int thread_id)
{
    typedef ConstNeighborhoodIterator<DifferenceImageType> SourceIterator;
    typedef ImageRegionConstIterator<MaskImageType> MaskIterator;
    typedef ImageRegionIterator<OutputImageType> DestinationIterator;

    typedef typename DifferenceImageType::IndexType IndexType;
    typedef typename DifferenceImageType::PixelType PixelType;

    typename SourceIterator::RadiusType radius;
    radius.Fill(1);
    SourceIterator source_it(radius, this->difference_image_, region);

    MaskIterator mask_it;
    if(!this->m_Mask.IsNull())
    {
        mask_it = MaskIterator(this->m_Mask, region);
    }

    DestinationIterator destination_it(this->GetOutput(), region);

    while(!source_it.IsAtEnd())
    {
        if(this->m_Mask.IsNull() || mask_it.Get() != this->m_MaskBackgroundValue)
        {
            InputRealType mean = 0;
            for(unsigned int i=0; i<source_it.Size(); ++i)
            {
                PixelType value;

                IndexType const index = source_it.GetIndex(i);
                if(this->difference_image_->GetRequestedRegion().IsInside(index))
                {
                    value = this->difference_image_->GetPixel(index);
                }
                else
                {
                    value = source_it.GetCenterPixel();
                }

                mean += value;
            }

            mean /= source_it.Size();

            // Likelihood of change
            double likelihood=std::pow(mean, 2.)*source_it.Size()*0.5/this->variance_;
            destination_it.Set(likelihood);
        }

        ++source_it;
        if(!this->m_Mask.IsNull())
        {
            ++mask_it;
        }
        ++destination_it;
    }
}

}

#endif // dbb086f8_040b_43f5_9ce7_b015deb978fe
