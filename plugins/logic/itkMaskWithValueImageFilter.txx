/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef b5f15137_0725_48ab_a822_92f27b5897d0
#define b5f15137_0725_48ab_a822_92f27b5897d0

#include "itkMaskWithValueImageFilter.h"

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

namespace itk
{

template<typename TInputImage, typename TMaskImage, typename TOutputImage>
TMaskImage const *
MaskWithValueImageFilter<TInputImage, TMaskImage, TOutputImage>
::GetMask()
{
    return this->GetInput(1);
}

template<typename TInputImage, typename TMaskImage, typename TOutputImage>
void
MaskWithValueImageFilter<TInputImage, TMaskImage, TOutputImage>
::SetMask(TMaskImage const * image)
{
    this->SetNthInput(1, const_cast<TMaskImage*>(image));
}

template<typename TInputImage, typename TMaskImage, typename TOutputImage>
void
MaskWithValueImageFilter<TInputImage, TMaskImage, TOutputImage>
::PrintSelf(std::ostream & os, Indent indent) const
{
    this->Superclass::PrintSelf(os, indent);
    os << indent << "BackgroundValue: "  << this->GetBackgroundValue() << std::endl;
    os << indent << "OutsideValue: "  << this->GetOutsideValue() << std::endl;
}

template<typename TInputImage, typename TMaskImage, typename TOutputImage>
void
MaskWithValueImageFilter<TInputImage, TMaskImage, TOutputImage>
::ThreadedGenerateData(typename InputImageType::RegionType const & region, ThreadIdType)
{
    typedef itk::ImageRegionConstIterator<InputImageType> SourceIteratorType;
    typedef itk::ImageRegionConstIterator<MaskImageType> MaskIteratorType;
    typedef itk::ImageRegionIterator<OutputImageType> DestinationIteratorType;

    SourceIteratorType source_it(this->GetInput(), region);
    MaskIteratorType mask_it(this->GetMask(), region);
    DestinationIteratorType destination_it(this->GetOutput(), region);

    while(!destination_it.IsAtEnd())
    {
        if(mask_it.Get() != this->m_BackgroundValue)
        {
            destination_it.Set(static_cast<OutputImagePixelType>(source_it.Get()));
        }
        else
        {
            destination_it.Set(this->m_OutsideValue);
        }
        ++source_it;
        ++mask_it;
        ++destination_it;
    }
}

template<typename TInputImage, typename TMaskImage, typename TOutputImage>
MaskWithValueImageFilter<TInputImage, TMaskImage, TOutputImage>
::MaskWithValueImageFilter()
: m_BackgroundValue(NumericTraits<MaskImagePixelType>::Zero),
  m_OutsideValue(NumericTraits<OutputImagePixelType>::Zero)
{
    this->SetNumberOfRequiredInputs(2);
}

template<typename TInputImage, typename TMaskImage, typename TOutputImage>
MaskWithValueImageFilter<TInputImage, TMaskImage, TOutputImage>
::~MaskWithValueImageFilter()
{
    // Nothing to do
}

}

#endif // b5f15137_0725_48ab_a822_92f27b5897d0
