#ifndef _7c033d14_c586_4b12_ae5b_4147868c33a9
#define _7c033d14_c586_4b12_ae5b_4147868c33a9

#include "itkClustersToAnnotationsCalculator.h"

#include <map>
#include <vector>

#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>

namespace itk
{

template<typename TImage>
void
ClustersToAnnotationsCalculator<TImage>
::Compute()
{
    this->m_Annotations.clear();

    typedef std::map<PixelType, std::vector<IndexType> > RegionsMapType;
    RegionsMapType regions;

    for(ImageRegionConstIteratorWithIndex<ImageType> it(
            this->m_Image, this->m_Image->GetRequestedRegion());
        !it.IsAtEnd(); ++it)
    {
        PixelType const value = it.Get();
        if(value != 0)
        {
            regions[value].push_back(it.GetIndex());
        }
    }

    // Compute a distance map for each region, use maximum as location of annotation
    for(typename RegionsMapType::const_iterator regions_it=regions.begin();
        regions_it!=regions.end(); ++regions_it)
    {
        typename RegionsMapType::mapped_type const & region = regions_it->second;

        // Find bounding box
        IndexType first = *(region.begin());
        IndexType last = *(region.begin());
        for(typename RegionsMapType::mapped_type::const_iterator it=region.begin();
            it != region.end(); ++it)
        {
            IndexType const index = *it;
            for(unsigned int d=0; d<ImageType::ImageDimension; ++d)
            {
                first[d] = std::min(first[d], index[d]);
                last[d] = std::max(last[d], index[d]);
            }
        }

        // Create sub-image
        typedef Image<signed int, ImageType::ImageDimension> SubImageType;
        typename SubImageType::Pointer sub_image = SubImageType::New();
        typename SubImageType::RegionType::SizeType size;
        for(unsigned int d=0; d<SubImageType::ImageDimension; ++d)
        {
            size[d] = last[d]-first[d]+1;
        }
        sub_image->SetRegions(typename SubImageType::RegionType(first, size));
        sub_image->Allocate();
        sub_image->FillBuffer(0);
        for(typename RegionsMapType::mapped_type::const_iterator it=region.begin();
            it != region.end(); ++it)
        {
            sub_image->SetPixel(*it, 1);
        }

        // Compute distance map
        typedef SignedMaurerDistanceMapImageFilter<SubImageType, SubImageType>
            DistanceMapFilterType;
        typename DistanceMapFilterType::Pointer distance_map_filter =
            DistanceMapFilterType::New();
        distance_map_filter->SetInput(sub_image);
        distance_map_filter->Update();

        // Find minimum, since "the inside is considered as having negative
        // distances" (cf. doc of SignedMaurerDistanceMapImageFilter)
        typedef MinimumMaximumImageCalculator<SubImageType> MinimumCalculatorType;
        typename MinimumCalculatorType::Pointer minimum_calculator = MinimumCalculatorType::New();
        minimum_calculator->SetImage(distance_map_filter->GetOutput());
        minimum_calculator->ComputeMinimum();

        this->m_Annotations[regions_it->first] = minimum_calculator->GetIndexOfMinimum();
    }
}

template<typename TImage>
std::vector<typename ClustersToAnnotationsCalculator<TImage>::MapType::key_type>
ClustersToAnnotationsCalculator<TImage>
::GetAnnotationsLabels() const
{
    std::vector<typename MapType::key_type> labels;
    labels.reserve(this->m_Annotations.size());

    for(typename MapType::const_iterator it=this->m_Annotations.begin();
        it!=this->m_Annotations.end(); ++it)
    {
        labels.push_back(it->first);
    }

    return labels;
}

template<typename TImage>
typename ClustersToAnnotationsCalculator<TImage>::IndexType
ClustersToAnnotationsCalculator<TImage>
::GetAnnotation(typename ClustersToAnnotationsCalculator<TImage>::PixelType label) const
{
    return this->m_Annotations.find(label)->second;
}

template<typename TImage>
ClustersToAnnotationsCalculator<TImage>
::ClustersToAnnotationsCalculator()
{
    // Nothing to do.
}

template<typename TImage>
ClustersToAnnotationsCalculator<TImage>
::~ClustersToAnnotationsCalculator()
{
    // Nothing to do.
}

template<typename TImage>
typename ClustersToAnnotationsCalculator<TImage>::MapType const &
ClustersToAnnotationsCalculator<TImage>
::GetAnnotations() const
{
    return this->m_Annotations;
}

}

#endif // _7c033d14_c586_4b12_ae5b_4147868c33a9
