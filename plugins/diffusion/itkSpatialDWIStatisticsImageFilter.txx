/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _fe52d872_9465_40b4_a6b8_741d2c51e3e6
#define _fe52d872_9465_40b4_a6b8_741d2c51e3e6

#include "itkSpatialDWIStatisticsImageFilter.h"

#include <ostream>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkNeighborhood.h>

namespace itk
{

template<typename TInputImage, typename TMeanImage, 
         typename TStandardDeviationImage, typename TMaskImage>
SpatialDWIStatisticsImageFilter<TInputImage, TMeanImage, 
                                TStandardDeviationImage, TMaskImage>
::SpatialDWIStatisticsImageFilter()
{
    this->SetSizePlane(3);
    this->SetSizeDepth(3);
}

template<typename TInputImage, typename TMeanImage, 
         typename TStandardDeviationImage, typename TMaskImage>
void
SpatialDWIStatisticsImageFilter<TInputImage, TMeanImage, 
                                TStandardDeviationImage, TMaskImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
    this->Superclass::PrintSelf(os,indent);
    os << indent << "Size in plane: " << this->GetSizePlane() << "\n";
    os << indent << "Size in depth: " << this->GetSizeDepth() << "\n";
}

template<typename TInputImage, typename TMeanImage, 
         typename TStandardDeviationImage, typename TMaskImage>
void 
SpatialDWIStatisticsImageFilter<TInputImage, TMeanImage, 
                                TStandardDeviationImage, TMaskImage>
::ThreadedGenerateData(OutputImageRegionType const & outputRegionForThread, int)
{
    typedef typename TInputImage::PixelType InputTensorType;
    typedef typename TMeanImage::PixelType OutputTensorType;
    
    typename TInputImage::ConstPointer tensor_image = this->GetInput();
    typename TInputImage::RegionType const tensor_region = 
        tensor_image->GetLargestPossibleRegion();
        
    typename TMaskImage::ConstPointer mask_image = this->GetMaskImage();
        
    typename TMeanImage::Pointer mean_image = 
        dynamic_cast<TMeanImage*>(this->ProcessObject::GetOutput(0));
    mean_image->FillBuffer(0);
    typename TStandardDeviationImage::Pointer standard_deviation_image = 
        dynamic_cast<TStandardDeviationImage*>(this->ProcessObject::GetOutput(1));
    standard_deviation_image->FillBuffer(0);

    // Build the neighborhood
    typedef itk::Neighborhood<typename TInputImage::PixelType, 
                              TInputImage::ImageDimension> NeighborhoodType;
    NeighborhoodType neighborhood;
    
    typename NeighborhoodType::SizeType neighborhood_size;
    neighborhood_size[0]=(this->m_SizePlane-1)/2; 
    neighborhood_size[1]=(this->m_SizePlane-1)/2; 
    neighborhood_size[2]=(this->m_SizeDepth-1)/2;
    neighborhood.SetRadius(neighborhood_size);

    typedef ImageRegionConstIteratorWithIndex<TMeanImage> IteratorType;
    for(IteratorType it(mean_image, outputRegionForThread); !it.IsAtEnd(); ++it)
    {
        typename TMeanImage::IndexType const & idx = it.GetIndex();

        // Skip pixels that are outside the mask or whose neighborhood is not 
        // fully inside the region
        bool skip=false;
        if(!mask_image.IsNull()) 
        {
            typename TInputImage::PointType point; 
            mean_image->TransformIndexToPhysicalPoint(idx, point);
            
            typename TMaskImage::IndexType mask_index;
            mask_image->TransformPhysicalPointToIndex(point, mask_index);
            
            if(!mask_image->GetLargestPossibleRegion().IsInside(mask_index) || 
               mask_image->GetPixel(mask_index) == 0)
            {
                skip=true;
            }
        }
        for(unsigned int i=0; i<neighborhood.Size() && !skip; ++i)
        {
            if(!tensor_region.IsInside(idx+neighborhood.GetOffset(i)))
            {
                skip=true;
            }
        }
        if(skip)
        {
            continue;
        }
        
        // first compute the mean
        OutputTensorType mean = mean_image->GetPixel(idx);
        mean.Fill(0.);
        for(unsigned int i=0; i<neighborhood.Size(); ++i)
        {
            InputTensorType const & tensor = tensor_image->GetPixel(
                idx+neighborhood.GetOffset(i));
            mean += tensor/neighborhood.Size();
        }

        // then compute the variance
        typename TStandardDeviationImage::PixelType variance = 0.;
        for(unsigned int i=0; i<neighborhood.Size(); ++i)
        {
            InputTensorType const & tensor = tensor_image->GetPixel(
                idx+neighborhood.GetOffset(i));
            variance += (
                (tensor[0] - mean[0])*(tensor[0] - mean[0]) +
                (tensor[3] - mean[3])*(tensor[3] - mean[3]) +
                (tensor[5] - mean[5])*(tensor[5] - mean[5]) +
                2.0*(tensor[1] - mean[1])*(tensor[1] - mean[1]) +
                2.0*(tensor[2] - mean[2])*(tensor[2] - mean[2]) +
                2.0*(tensor[4] - mean[4])*(tensor[4] - mean[4])
            )/neighborhood.Size();
        }
        // Clamp variance and compute standard deviation
        variance = (variance>=0)?variance:0;
        standard_deviation_image->SetPixel(idx, variance/6.0);
    }
}

}

#endif // _fe52d872_9465_40b4_a6b8_741d2c51e3e6
