#ifndef _8fb5cf79_e4f0_4c47_ba23_39d6d96b6a7a
#define _8fb5cf79_e4f0_4c47_ba23_39d6d96b6a7a

#include "itkBinaryThresholdWithRadiusImageFunction.h"

#include <itkNumericTraits.h>

namespace itk
{

template<typename TInputImage, typename TCoordRep>
bool
BinaryThresholdWithRadiusImageFunction<TInputImage, TCoordRep>
::Evaluate(PointType const & point) const
{
    IndexType index;
    this->ConvertPointToNearestIndex(point, index);
    return this->EvaluateAtIndex(index);
}

template<typename TInputImage, typename TCoordRep>
bool
BinaryThresholdWithRadiusImageFunction<TInputImage, TCoordRep>
::EvaluateAtContinuousIndex(ContinuousIndexType const & continuous_index) const
{
    IndexType index;
    this->ConvertContinuousIndexToNearestIndex(continuous_index, index);
    return this->EvaluateAtIndex(index);
}

template<typename TInputImage, typename TCoordRep>
bool
BinaryThresholdWithRadiusImageFunction<TInputImage, TCoordRep>
::EvaluateAtIndex(IndexType const & index) const
{
    // Check for intensity
    InputPixelType const value = this->GetInputImage()->GetPixel(index);
    bool const threshold_ok = (this->m_Lower <= value && value <= this->m_Upper);

    if(threshold_ok)
    {
        // Check for distance
        bool distance_ok = false;
        for(typename std::vector<IndexType>::const_iterator seeds_it = this->m_Seeds.begin();
            seeds_it != this->m_Seeds.end(); ++seeds_it)
        {
            typename InputImageType::OffsetType const offset = index-(*seeds_it);
            float distance = 0;
            for(unsigned int d=0; d<InputImageType::ImageDimension; ++d)
            {
                distance += offset[d]*offset[d];
            }

            if(distance < this->m_RadiusSquared)
            {
                distance_ok = true;
                break;
            }
        }

        return distance_ok;
    }
    else
    {
        return false;
    }
}

template<typename TInputImage, typename TCoordRep>
void
BinaryThresholdWithRadiusImageFunction<TInputImage, TCoordRep>
::SetRadius(float radius)
{
    this->m_Radius = radius;
    this->m_RadiusSquared = radius*radius;
}

template<typename TInputImage, typename TCoordRep>
void
BinaryThresholdWithRadiusImageFunction<TInputImage, TCoordRep>
::SetSeed(IndexType const & seed)
{
    this->ClearSeeds();
    this->AddSeed(seed);
}

template<typename TInputImage, typename TCoordRep>
void
BinaryThresholdWithRadiusImageFunction<TInputImage, TCoordRep>
::AddSeed(IndexType const & seed)
{
    this->m_Seeds.push_back(seed);
    this->Modified();
}

template<typename TInputImage, typename TCoordRep>
void
BinaryThresholdWithRadiusImageFunction<TInputImage, TCoordRep>
::ClearSeeds()
{
    if(!m_Seeds.empty())
    {
        this->m_Seeds.clear();
        this->Modified();
    }
}

template<typename TInputImage, typename TCoordRep>
void
BinaryThresholdWithRadiusImageFunction<TInputImage, TCoordRep>
::ThresholdAbove(InputPixelType threshold)
{
    if(this->m_Lower != threshold ||
       this->m_Upper != NumericTraits<InputPixelType>::max())
    {
        this->m_Lower = threshold;
        this->m_Upper = NumericTraits<InputPixelType>::max();
        this->Modified();
    }
}

template<typename TInputImage, typename TCoordRep>
void
BinaryThresholdWithRadiusImageFunction<TInputImage, TCoordRep>
::ThresholdBelow(InputPixelType threshold)
{
    if(this->m_Lower != NumericTraits<InputPixelType>::NonpositiveMin() ||
       this->m_Upper != threshold)
    {
        this->m_Lower = NumericTraits<InputPixelType>::NonpositiveMin();
        this->m_Upper = threshold;
        this->Modified();
    }
}

template<typename TInputImage, typename TCoordRep>
void
BinaryThresholdWithRadiusImageFunction<TInputImage, TCoordRep>
::ThresholdBetween(InputPixelType lower, InputPixelType upper)
{
    if(this->m_Lower != lower || this->m_Upper != upper)
    {
        this->m_Lower = lower;
        this->m_Upper = upper;
        this->Modified();
    }
}

template<typename TInputImage, typename TCoordRep>
BinaryThresholdWithRadiusImageFunction<TInputImage, TCoordRep>
::BinaryThresholdWithRadiusImageFunction()
: m_Lower(NumericTraits<InputPixelType>::NonpositiveMin()),
  m_Upper(NumericTraits<InputPixelType>::max()),
  m_Radius(1.f), m_RadiusSquared(1.f)
{
}

template<typename TInputImage, typename TCoordRep>
BinaryThresholdWithRadiusImageFunction<TInputImage, TCoordRep>
::~BinaryThresholdWithRadiusImageFunction()
{
    // Nothing to do
}

template<typename TInputImage, typename TCoordRep>
void
BinaryThresholdWithRadiusImageFunction<TInputImage, TCoordRep>
::PrintSelf(std::ostream & os, Indent indent) const
{
    Superclass::PrintSelf(os, indent);

    os << indent << "Lower: " << this->m_Lower << std::endl;
    os << indent << "Upper: " << this->m_Upper << std::endl;
    os << indent << "Radius: " << this->m_Radius << std::endl;
    os << indent << "Seeds: " << this->m_Seeds.size() << std::endl;
}

} // namespace itk

#endif // _8fb5cf79_e4f0_4c47_ba23_39d6d96b6a7a
