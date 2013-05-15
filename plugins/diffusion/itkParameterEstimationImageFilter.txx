/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _73aa1937_258f_46fd_9080_4783e5971991
#define _73aa1937_258f_46fd_9080_4783e5971991

#include "itkParameterEstimationImageFilter.h"

#include <itkImageRegionIteratorWithIndex.h>
#include <itkNeighborhood.h>

namespace itk
{

template<typename TTensorImage, typename TMeanImage, typename TVarianceImage>
TMeanImage const *
ParameterEstimationImageFilter<TTensorImage, TMeanImage, TVarianceImage>
::GetMeanImage() const
{
    return dynamic_cast<TMeanImage const *>(this->ProcessObject::GetOutput(0));
}

template<typename TTensorImage, typename TMeanImage, typename TVarianceImage>
TVarianceImage const *
ParameterEstimationImageFilter<TTensorImage, TMeanImage, TVarianceImage>
::GetVarianceImage() const
{
    return dynamic_cast<TVarianceImage const *>(this->ProcessObject::GetOutput(1));
}

template<typename TTensorImage, typename TMeanImage, typename TVarianceImage>
ParameterEstimationImageFilter<TTensorImage, TMeanImage, TVarianceImage>
::ParameterEstimationImageFilter()
{
    this->SetSizePlane(3);
    this->SetSizeDepth(3);
    this->SetMask(NULL);

    this->SetNumberOfRequiredOutputs(2);
    this->SetNthOutput(0,this->MakeOutput(0));
    this->SetNthOutput(1,this->MakeOutput(1));
}

template<typename TTensorImage, typename TMeanImage, typename TVarianceImage>
void
ParameterEstimationImageFilter<TTensorImage, TMeanImage, TVarianceImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
    Superclass::PrintSelf(os,indent);
    os << indent << "Size in plane: " << this->GetSizePlane() << "\n";
    os << indent << "Size in depth: " << this->GetSizeDepth() << "\n";
    os << indent << "Mask: \n";
    if(this->m_Mask.IsNull())
    {
        os << indent.GetNextIndent() << "None\n";
    }
    else
    {
        this->m_Mask->Print(os, indent.GetNextIndent());
    }
}

template<typename TTensorImage, typename TMeanImage, typename TVarianceImage>
DataObject::Pointer 
ParameterEstimationImageFilter<TTensorImage, TMeanImage, TVarianceImage>
::MakeOutput(unsigned int index)
{
    DataObject::Pointer output;
    
    if(index == 0)
    {
        output = TMeanImage::New().GetPointer();
    }
    else if(index == 1)
    {
        output = TVarianceImage::New().GetPointer();
    }
    else
    {
        std::cerr << "No output " << index << std::endl;
        output = NULL;
    }
    
    return output.GetPointer();
}

template<typename TTensorImage, typename TMeanImage, typename TVarianceImage>
void
ParameterEstimationImageFilter<TTensorImage, TMeanImage, TVarianceImage>
::AllocateOutputs()
{
    typename TMeanImage::Pointer mean_image = dynamic_cast<TMeanImage*>(
        this->ProcessObject::GetOutput(0));
    if(mean_image) 
    {
        mean_image->SetBufferedRegion(mean_image->GetRequestedRegion());
        mean_image->SetVectorLength(6);
        mean_image->Allocate();
    }
    
    typename TVarianceImage::Pointer variance_image = dynamic_cast<TVarianceImage*>(
        this->ProcessObject::GetOutput(1));
    if(variance_image)
    {
        variance_image->SetBufferedRegion(variance_image->GetRequestedRegion());
        variance_image->Allocate();
    }
}

template<typename TTensorImage, typename TMeanImage, typename TVarianceImage>
void 
ParameterEstimationImageFilter<TTensorImage, TMeanImage, TVarianceImage>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int )
{
    typedef typename TTensorImage::PixelType InputTensorType;
    typedef typename TMeanImage::PixelType OutputTensorType;
    
    typename TTensorImage::ConstPointer tensor_image = this->GetInput();
    typename TTensorImage::RegionType const tensor_region = 
        tensor_image->GetLargestPossibleRegion();
        
    typename TMeanImage::Pointer mean_image = dynamic_cast<TMeanImage*>(
        this->ProcessObject::GetOutput(0));
    mean_image->FillBuffer(0);
    typename TVarianceImage::Pointer variance_image = dynamic_cast<TVarianceImage*>(
        this->ProcessObject::GetOutput(1));
    variance_image->FillBuffer(0);

    // Build the neighborhood
    typedef itk::Neighborhood<typename TTensorImage::PixelType, 
        TTensorImage::ImageDimension> NeighborhoodType;
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
        if(!this->m_Mask.IsNull()) 
        {
            typename TTensorImage::PointType point; 
            mean_image->TransformIndexToPhysicalPoint(idx, point);
            
            typename MaskType::IndexType mask_index;
            this->m_Mask->TransformPhysicalPointToIndex(point, mask_index);
            
            if(!this->m_Mask->GetLargestPossibleRegion().IsInside(mask_index) || 
               this->m_Mask->GetPixel(mask_index) == 0)
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
        typename TVarianceImage::PixelType variance = 0.;
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
        variance_image->SetPixel(idx, variance);
    }
}

}

#endif // _73aa1937_258f_46fd_9080_4783e5971991
