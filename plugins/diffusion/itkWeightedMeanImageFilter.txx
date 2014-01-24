/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _itkWeightedMeanImageFilter_txx
#define _itkWeightedMeanImageFilter_txx

#include "itkWeightedMeanImageFilter.h"

#include <ostream>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkNeighborhood.h>

namespace itk
{

template<typename TInputImage, typename TOutputImage, typename TTensorImage>
WeightedMeanImageFilter<TInputImage, TOutputImage, TTensorImage>
::WeightedMeanImageFilter()
{
    this->SetSizePlaneX(3);
    this->SetSizePlaneY(3);
    this->SetSizeDepth(3);
    this->SetTensorImage(NULL);
}

template<typename TInputImage, typename TOutputImage, typename TTensorImage>
void WeightedMeanImageFilter<TInputImage, TOutputImage, TTensorImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
    Superclass::PrintSelf(os,indent);
    os << indent << "Size X in plane: " << this->GetSizePlaneX() << "\n";
    os << indent << "Size Y in plane: " << this->GetSizePlaneY() << "\n";
    os << indent << "Size in depth: " << this->GetSizeDepth() << "\n";
    
    // Print Tensor Image
    os << indent << "Tensor image: \n";
    if(this->m_TensorImage.IsNull())
    {
        os << indent.GetNextIndent() << "None\n";
    }
    else
    {
        this->m_TensorImage->Print(os, indent.GetNextIndent());
    }
}

template<typename TInputImage, typename TOutputImage, typename TTensorImage>
void WeightedMeanImageFilter<TInputImage, TOutputImage, TTensorImage>
::AllocateOutputs()
{
    OutputImagePointer output = dynamic_cast< OutputImageType *>(this->ProcessObject::GetOutput(0));
    if(output) 
    {
        output->SetBufferedRegion(output->GetRequestedRegion());
        output->SetVectorLength(6);
        output->Allocate();
    }
}

template<typename TInputImage, typename TOutputImage, typename TTensorImage>
void WeightedMeanImageFilter<TInputImage, TOutputImage, TTensorImage>
::BeforeThreadedGenerateData()
{
    OutputImagePointer output = this->GetOutput();
    OutputImagePixelType zero(6);
    zero.Fill(0.);
    output->FillBuffer(zero);
}

template<typename TInputImage, typename TOutputImage, typename TTensorImage>
void WeightedMeanImageFilter<TInputImage, TOutputImage, TTensorImage>
::ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, int)
{
    typedef typename TensorImageType::PixelType TensorImagePixelType;
    typedef typename TTensorImage::ConstPointer TensorImageConstPointer;
    typedef typename TTensorImage::RegionType const TensorImageRegionType;
    
    InputImageConstPointer weighted_image = this->GetInput();
    
    TensorImageConstPointer tensor_image = this->GetTensorImage();
    TensorImageRegionType tensor_region = tensor_image->GetRequestedRegion();
    
    OutputImagePointer output = this->GetOutput();
    
    // Build the neighborhood
    typedef itk::Neighborhood<typename TTensorImage::PixelType, 
                              TTensorImage::ImageDimension> NeighborhoodType;
    NeighborhoodType neighborhood;
    
    typename NeighborhoodType::SizeType neighborhood_size;
    neighborhood_size[0]=(this->m_SizePlaneX-1)/2; 
    neighborhood_size[1]=(this->m_SizePlaneY-1)/2; 
    neighborhood_size[2]=(this->m_SizeDepth-1)/2;
    neighborhood.SetRadius(neighborhood_size);
    
    typedef ImageRegionConstIteratorWithIndex<TOutputImage> IteratorType;
    for(IteratorType it(output, outputRegionForThread); !it.IsAtEnd(); ++it)
    {
        typename TOutputImage::IndexType const & idx = it.GetIndex();
        
        // Skip pixels that are outside the mask or whose neighborhood is not 
        // fully inside the region
        bool process=true;
        for(unsigned int i=0; i<neighborhood.Size() && process; ++i)
        {
            if(!tensor_region.IsInside(idx+neighborhood.GetOffset(i)))
            {
                process=false;
            }
        }
        
        typename TOutputImage::PointType point; 
        output->TransformIndexToPhysicalPoint(idx, point);
        
        typename TInputImage::IndexType weighted_index;
        weighted_image->TransformPhysicalPointToIndex(point, weighted_index);
        
        if(!weighted_image->GetLargestPossibleRegion().IsInside(weighted_index) || weighted_image->GetPixel(weighted_index) == 0)
        {
            process=false;
        }
        
        if(process)
        {
            //compute the mean
            InputImagePixelType mean_proba = 0.0;
            OutputImagePixelType mean(6);
            mean.Fill(0.0);
            for(unsigned int i=0; i<neighborhood.Size(); ++i)
            {
                TensorImagePixelType const & tensor = tensor_image->GetPixel(idx+neighborhood.GetOffset(i));
                InputImagePixelType proba = weighted_image->GetPixel(weighted_index+neighborhood.GetOffset(i));
                mean += tensor * proba;
                mean_proba += proba;
            }
            mean = mean / mean_proba;
            output->SetPixel(idx, mean);
        }
     }
    
}

}

#endif
