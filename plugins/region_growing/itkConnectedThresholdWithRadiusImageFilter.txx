/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _5fd06903_5e55_421a_90aa_8c8a42877661
#define _5fd06903_5e55_421a_90aa_8c8a42877661

#include "itkConnectedThresholdWithRadiusImageFilter.h"

#include <itkNumericTraits.h>
#include <itkFloodFilledImageFunctionConditionalIterator.h>
#include <itkProgressReporter.h>

#include "itkBinaryThresholdWithRadiusImageFunction.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
typename ConnectedThresholdWithRadiusImageFilter<TInputImage, TOutputImage>::InputImagePixelType
ConnectedThresholdWithRadiusImageFilter<TInputImage, TOutputImage>
::GetLower() const
{
    typename InputPixelObjectType::Pointer lower = const_cast<Self*>(this)->GetLowerInput();
    return lower->Get();
}

template <class TInputImage, class TOutputImage>
void
ConnectedThresholdWithRadiusImageFilter<TInputImage, TOutputImage>
::SetLower(InputImagePixelType threshold)
{
    // first check to see if anything changed
    typename InputPixelObjectType::Pointer lower=this->GetLowerInput();
    if(lower && lower->Get() == threshold)
    {
        return;
    }

    // create a data object to use as the input and to store this
    // threshold. we always create a new data object to use as the input
    // since we do not want to change the value in any current input
    // (the current input could be the output of another filter or the
    // current input could be used as an input to several filters)
    lower = InputPixelObjectType::New();
    this->ProcessObject::SetNthInput(1, lower);

    lower->Set(threshold);
    this->Modified();
}

template <class TInputImage, class TOutputImage>
typename ConnectedThresholdWithRadiusImageFilter<TInputImage, TOutputImage>::InputImagePixelType
ConnectedThresholdWithRadiusImageFilter<TInputImage, TOutputImage>
::GetUpper() const
{
    typename InputPixelObjectType::Pointer upper = const_cast<Self*>(this)->GetUpperInput();
    return upper->Get();
}

template <class TInputImage, class TOutputImage>
void
ConnectedThresholdWithRadiusImageFilter<TInputImage, TOutputImage>
::SetUpper(InputImagePixelType threshold)
{
    // first check to see if anything changed
    typename InputPixelObjectType::Pointer upper=this->GetUpperInput();
    if(upper && upper->Get() == threshold)
    {
        return;
    }

    // create a data object to use as the input and to store this
    // threshold. we always create a new data object to use as the input
    // since we do not want to change the value in any current input
    // (the current input could be the output of another filter or the
    // current input could be used as an input to several filters)
    upper = InputPixelObjectType::New();
    this->ProcessObject::SetNthInput(2, upper);

    upper->Set(threshold);
    this->Modified();
}

template <class TInputImage, class TOutputImage>
typename ConnectedThresholdWithRadiusImageFilter<TInputImage, TOutputImage>::InputPixelObjectType *
ConnectedThresholdWithRadiusImageFilter<TInputImage, TOutputImage>
::GetLowerInput()
{
    typename InputPixelObjectType::Pointer lower =
        static_cast<InputPixelObjectType *>(this->ProcessObject::GetInput(1));
    if(!lower)
    {
        // no input object available, create a new one and set it to the
        // default threshold
        lower = InputPixelObjectType::New();
        lower->Set(NumericTraits<InputImagePixelType>::NonpositiveMin());
        this->ProcessObject::SetNthInput(1, lower);
    }

    return lower;
}

template <class TInputImage, class TOutputImage>
void
ConnectedThresholdWithRadiusImageFilter<TInputImage, TOutputImage>
::SetLowerInput(InputPixelObjectType const * input)
{
    if(input != this->GetLowerInput())
    {
        this->ProcessObject::SetNthInput(1, const_cast<InputPixelObjectType*> (input));
        this->Modified();
    }
}

template <class TInputImage, class TOutputImage>
typename ConnectedThresholdWithRadiusImageFilter<TInputImage, TOutputImage>::InputPixelObjectType *
ConnectedThresholdWithRadiusImageFilter<TInputImage, TOutputImage>
::GetUpperInput()
{
    typename InputPixelObjectType::Pointer upper =
        static_cast<InputPixelObjectType *>(this->ProcessObject::GetInput(2));
    if(!upper)
    {
        // no input object available, create a new one and set it to the
        // default threshold
        upper = InputPixelObjectType::New();
        upper->Set(NumericTraits<InputImagePixelType>::max());
        this->ProcessObject::SetNthInput(2, upper);
    }

    return upper;
}

template <class TInputImage, class TOutputImage>
void
ConnectedThresholdWithRadiusImageFilter<TInputImage, TOutputImage>
::SetUpperInput(InputPixelObjectType const * input)
{
    if(input != this->GetUpperInput())
    {
        this->ProcessObject::SetNthInput(2, const_cast<InputPixelObjectType*> (input));
        this->Modified();
    }
}

template <class TInputImage, class TOutputImage>
void
ConnectedThresholdWithRadiusImageFilter<TInputImage, TOutputImage>
::SetSeed(InputImageIndexType const & seed)
{
    this->ClearSeeds();
    this->AddSeed(seed);
}

template <class TInputImage, class TOutputImage>
void
ConnectedThresholdWithRadiusImageFilter<TInputImage, TOutputImage>
::AddSeed(InputImageIndexType const & seed)
{
    this->m_Seeds.push_back(seed);
    this->Modified();
}

template <class TInputImage, class TOutputImage>
void
ConnectedThresholdWithRadiusImageFilter<TInputImage, TOutputImage>
::ClearSeeds()
{
    if(!this->m_Seeds.empty())
    {
        this->m_Seeds.clear();
        this->Modified();
    }
}

template <class TInputImage, class TOutputImage>
void
ConnectedThresholdWithRadiusImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream & os, Indent indent) const
{
    this->Superclass::PrintSelf(os, indent);
    os << indent << "Radius: " << this->m_Radius << "\n";
}

template <class TInputImage, class TOutputImage>
void 
ConnectedThresholdWithRadiusImageFilter<TInputImage,TOutputImage>
::GenerateData()
{
    typedef typename Superclass::OutputImageType OutputImageType;
    typedef typename Superclass::InputImageConstPointer InputImageConstPointer;
    typedef typename Superclass::OutputImagePointer OutputImagePointer;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    InputImageConstPointer inputImage = this->GetInput();
    OutputImagePointer outputImage = this->GetOutput();

    typename InputPixelObjectType::Pointer lowerThreshold=this->GetLowerInput();
    typename InputPixelObjectType::Pointer upperThreshold=this->GetUpperInput();

    this->m_Lower = lowerThreshold->Get();
    this->m_Upper = upperThreshold->Get();

    // Zero the output
    OutputImageRegionType region = outputImage->GetRequestedRegion();
    outputImage->SetBufferedRegion(region);
    outputImage->Allocate();
    outputImage->FillBuffer(NumericTraits<OutputImagePixelType>::Zero);
  
    typedef BinaryThresholdWithRadiusImageFunction<InputImageType, double> FunctionType;

    typename FunctionType::Pointer function = FunctionType::New();
    function->SetInputImage(inputImage);
    function->ThresholdBetween(this->m_Lower, this->m_Upper);
    function->SetRadius(this->m_Radius);
    for(typename std::vector<InputImageIndexType>::const_iterator seeds_it = this->m_Seeds.begin();
        seeds_it != this->m_Seeds.end(); ++seeds_it)
    {
        function->AddSeed(*seeds_it);
    }

    ProgressReporter progress(this, 0, region.GetNumberOfPixels());

    typedef FloodFilledImageFunctionConditionalIterator<OutputImageType, FunctionType> IteratorType;
    IteratorType it(outputImage, function, this->m_Seeds);
    it.GoToBegin();

    while(!it.IsAtEnd())
    {
        it.Set(this->m_ReplaceValue);
        ++it;
        progress.CompletedPixel();
    }
}

template <class TInputImage, class TOutputImage>
ConnectedThresholdWithRadiusImageFilter<TInputImage, TOutputImage>
::ConnectedThresholdWithRadiusImageFilter()
: m_Lower(NumericTraits<InputImagePixelType>::NonpositiveMin()),
  m_Upper(NumericTraits<InputImagePixelType>::max()),
  m_Radius(1.f),
  m_ReplaceValue(NumericTraits<OutputImagePixelType>::One)
{
    typename InputPixelObjectType::Pointer lower = InputPixelObjectType::New();
    lower->Set(NumericTraits<InputImagePixelType>::NonpositiveMin());
    this->ProcessObject::SetNthInput(1, lower);

    typename InputPixelObjectType::Pointer upper = InputPixelObjectType::New();
    upper->Set(NumericTraits<InputImagePixelType>::max());
    this->ProcessObject::SetNthInput(2, upper);
}

template <class TInputImage, class TOutputImage>
ConnectedThresholdWithRadiusImageFilter<TInputImage, TOutputImage>
::~ConnectedThresholdWithRadiusImageFilter()
{
    // Nothing to do
}

} // namespace itk

#endif // _5fd06903_5e55_421a_90aa_8c8a42877661
