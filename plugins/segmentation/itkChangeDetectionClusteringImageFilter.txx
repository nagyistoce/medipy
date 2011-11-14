/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef e3abd139_d767_4fa9_97c5_d6ce573e6099
#define e3abd139_d767_4fa9_97c5_d6ce573e6099

#include "itkChangeDetectionClusteringImageFilter.h"

#include <iterator>
#include <limits>
#include <set>

#include <itkConnectedComponentImageFilter.h>
#include <itkConstNeighborhoodIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkHistogram.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkThresholdImageFilter.h>

namespace itk
{

template<typename TInputImage, typename TMaskImage, typename TOutputImage>
ChangeDetectionClusteringImageFilter<TInputImage, TMaskImage, TOutputImage>
::ChangeDetectionClusteringImageFilter()
: m_Mask(0), m_MaskBackgroundValue(0), m_MinimumClusterSize(0),
  m_MaximumNumberOfClusters(std::numeric_limits<unsigned long>::max())
{
}

template<typename TInputImage, typename TMaskImage, typename TOutputImage>
ChangeDetectionClusteringImageFilter<TInputImage, TMaskImage, TOutputImage>
::~ChangeDetectionClusteringImageFilter()
{
    // Nothing to do.
}

template<typename TInputImage, typename TMaskImage, typename TOutputImage>
void
ChangeDetectionClusteringImageFilter<TInputImage, TMaskImage, TOutputImage>
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
    this->Superclass::PrintSelf(os, indent);
    os << indent << "Mask : " << this->m_Mask << "\n";
    os << indent << "Mask background value : " << this->m_MaskBackgroundValue << "\n";
}

template<typename TInputImage, typename TMaskImage, typename TOutputImage>
void
ChangeDetectionClusteringImageFilter<TInputImage, TMaskImage, TOutputImage>
::GenerateData()
{
    InputImagePixelType const threshold = this->compute_threshold();

    // Threshold input
    typedef itk::ThresholdImageFilter<InputImageType> InputThresholdFilterType;
    typename InputThresholdFilterType::Pointer input_threshold_filter = InputThresholdFilterType::New();
    input_threshold_filter->SetInput(this->GetInput());
    input_threshold_filter->SetOutsideValue(0);
    input_threshold_filter->ThresholdBelow(threshold);

    // Label the connected components, use mask if available
    // Use a specific image type to avoid problems if OutputImageType is float:
    // with VC2008, static_cast<unsigned long int>(std::numeric_limits<float>::max())
    // is 0. This test is used in ConnectedComponentImageFilter::ThreadedGenerateData
    typedef itk::Image<unsigned long, InputImageType::ImageDimension>
    	ConnectedComponentsImageType;
    typedef itk::ConnectedComponentImageFilter<InputImageType, ConnectedComponentsImageType, MaskImageType>
        ConnectedComponentFilterType;
    typename ConnectedComponentFilterType::Pointer connected_component_filter =
        ConnectedComponentFilterType::New();
    connected_component_filter->SetInput(input_threshold_filter->GetOutput());
    if(!this->m_Mask.IsNull())
    {
        connected_component_filter->SetMaskImage(this->m_Mask);
    }
    connected_component_filter->SetBackgroundValue(0);

    // Remove small objects
    typedef itk::RelabelComponentImageFilter<ConnectedComponentsImageType, OutputImageType>
        RelabelFilterType;

    typename RelabelFilterType::Pointer relabel_filter_1 = RelabelFilterType::New();
    relabel_filter_1->SetInput(connected_component_filter->GetOutput());
    relabel_filter_1->SetMinimumObjectSize(this->m_MinimumClusterSize);
    relabel_filter_1->Update();

    this->discard_clusters_close_to_mask_boundary(relabel_filter_1->GetOutput());

    // Sort the label by size
    typedef itk::RelabelComponentImageFilter<OutputImageType, OutputImageType>
		RelabelFilterOutputImageType;
    typename RelabelFilterOutputImageType::Pointer relabel_filter_2 = RelabelFilterOutputImageType::New();
    relabel_filter_2->SetInput(relabel_filter_1->GetOutput());
    // Keep only n labels
    typedef itk::ThresholdImageFilter<OutputImageType> OutputThresholdFilterType;
    typename OutputThresholdFilterType::Pointer output_threshold_filter = OutputThresholdFilterType::New();
    output_threshold_filter->SetInput(relabel_filter_2->GetOutput());
    output_threshold_filter->SetOutsideValue(0);
    output_threshold_filter->ThresholdAbove(this->m_MaximumNumberOfClusters);

    output_threshold_filter->GraftOutput(this->GetOutput());
    output_threshold_filter->Update();
    this->GraftOutput(output_threshold_filter->GetOutput());
}

template<typename TInputImage, typename TMaskImage, typename TOutputImage>
typename TInputImage::PixelType
ChangeDetectionClusteringImageFilter<TInputImage, TMaskImage, TOutputImage>
::compute_threshold()
{
    typedef typename itk::ImageRegionConstIterator<InputImageType>
        InputConstIteratorType;
    typedef typename itk::ImageRegionConstIterator<MaskImageType>
        MaskConstIteratorType;
    // Only keep the top 2% values, based on cumulative histogram

    typedef itk::MinimumMaximumImageCalculator<InputImageType> MinimumMaximumImageCalculator;
    typename MinimumMaximumImageCalculator::Pointer calculator = MinimumMaximumImageCalculator::New();
    calculator->SetImage(this->GetInput());
    calculator->Compute();

    // Compute histogram
    typedef itk::Statistics::Histogram<> Histogram;
    Histogram::Pointer histogram = Histogram::New();
    {
        Histogram::SizeType size; size.Fill(1000);
        Histogram::MeasurementVectorType min; min.Fill(calculator->GetMinimum());
        Histogram::MeasurementVectorType max; max.Fill(calculator->GetMaximum());
        histogram->Initialize(size, min, max);
    }

    if(this->m_Mask.IsNull())
    {
        itkWarningMacro(<<"No mask supplied, using the whole image");

        // Exclude background
        InputConstIteratorType input_it(this->GetInput(),
                                        this->GetInput()->GetRequestedRegion());
        while(!input_it.IsAtEnd())
        {
            if(input_it.Get() != 0)
            {
                histogram->IncreaseFrequency(input_it.Get(), 1);
            }
            ++input_it;
        }
    }
    else
    {
        // Make sure mask is up-to-date
        this->m_Mask->Update();

        InputConstIteratorType input_it(this->GetInput(),
                                        this->GetInput()->GetRequestedRegion());
        MaskConstIteratorType mask_it(this->m_Mask,
                                      this->GetInput()->GetRequestedRegion());

        while(!input_it.IsAtEnd())
        {
            if(mask_it.Get() != this->m_MaskBackgroundValue)
            {
                histogram->IncreaseFrequency(input_it.Get(), 1);
            }
            ++input_it;
            ++mask_it;
        }
    }

    // Find 0.98 threshold
    Histogram::InstanceIdentifier n=0;
    float cdf_value = 0;
    do
    {
        cdf_value += histogram->GetFrequency(n)/histogram->GetTotalFrequency();
        ++n;
    }
    while(cdf_value < 0.98);
    InputImagePixelType const threshold = histogram->GetBinMax(0, n);

    return threshold;
}

template<typename TInputImage, typename TMaskImage, typename TOutputImage>
void
ChangeDetectionClusteringImageFilter<TInputImage, TMaskImage, TOutputImage>
::discard_clusters_close_to_mask_boundary(TOutputImage * image)
{

    std::set<OutputImagePixelType> discarded_labels;
    MaskImageType * mask = this->m_Mask;

    // Find labels connected to the mask boundary
    typedef itk::ConstNeighborhoodIterator<OutputImageType> NeighborhoodIteratorType;
    typename NeighborhoodIteratorType::RadiusType radius;
    radius.Fill(1);

    for(NeighborhoodIteratorType nit(radius, image, image->GetRequestedRegion());
        !nit.IsAtEnd(); ++nit)
    {
        typename OutputImageType::PixelType const center_pixel = nit.GetCenterPixel();

        if(center_pixel == 0)
        {
            // No cluster at this point, skip the neighborhood scan.
            continue;
        }

        for(unsigned int i=0; i<nit.Size(); ++i)
        {
            // TODO : skip center
            typename NeighborhoodIteratorType::IndexType const neighbor = nit.GetIndex(i);
            if(mask->GetRequestedRegion().IsInside(neighbor) &&
               mask->GetPixel(neighbor) == this->m_MaskBackgroundValue)
            {
                discarded_labels.insert(center_pixel);
                break;
            }
        }
    }

    // Discard the labels that are connected to the mask boundary
    for(itk::ImageRegionIterator<OutputImageType> it(image, image->GetRequestedRegion());
        !it.IsAtEnd(); ++it)
    {
        if(discarded_labels.count(it.Get()) != 0)
        {
            it.Set(0);
        }
    }
}

}

#endif // e3abd139_d767_4fa9_97c5_d6ce573e6099
