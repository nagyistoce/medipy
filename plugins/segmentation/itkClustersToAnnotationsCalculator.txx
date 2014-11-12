#ifndef _7c033d14_c586_4b12_ae5b_4147868c33a9
#define _7c033d14_c586_4b12_ae5b_4147868c33a9

#include "itkClustersToAnnotationsCalculator.h"

#include <algorithm>
#include <functional>
#include <map>
#include <stdexcept>
#include <vector>

#include <itkFloodFilledImageFunctionConditionalConstIterator.h>
#include <itkImageFunction.h>
#include <itkImageRegionConstIteratorWithIndex.h>

namespace itk
{

/**
 * @brief Image function always returning true.
 *
 * This can be used to get an "infinite" propagation with flood-filled iterators.
 */
template<typename TInputImage, typename TCoordRep=float>
class AlwaysTrueImageFunction : public ImageFunction<TInputImage, bool, TCoordRep>
{
public:
    /** Standard class typedefs. */
    typedef AlwaysTrueImageFunction Self;
    typedef ImageFunction<TInputImage, bool, TCoordRep> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<const Self> ConstPointer;

    /** Run-time type information (and related methods). */
    itkTypeMacro(AlwaysTrueImageFunction, ImageFunction);

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** InputImageType typedef support. */
    typedef typename Superclass::InputImageType InputImageType;

    /** Typedef to describe the type of pixel. */
    typedef typename TInputImage::PixelType PixelType;

    /** Dimension underlying input image. */
    itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

    /** Point typedef support. */
    typedef typename Superclass::PointType PointType;

    /** Index typedef support. */
    typedef typename Superclass::IndexType IndexType;

    /** ContinuousIndex typedef support. */
    typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

    virtual bool Evaluate(PointType const & point) const
    {
        return true;
    }

    virtual bool EvaluateAtContinuousIndex(ContinuousIndexType const & index) const
    {
        return true;
    }

    virtual bool EvaluateAtIndex(IndexType const & index) const
    {
        return true;
    }

protected:
    AlwaysTrueImageFunction() {}
    ~AlwaysTrueImageFunction() {};

private:
    AlwaysTrueImageFunction(Self const &); //purposely not implemented
    void operator=(Self const &); //purposely not implemented
};

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

    // Compute the center of mass of each region
    std::map<PixelType, IndexType> centers;
    for(typename RegionsMapType::const_iterator regions_it=regions.begin();
        regions_it!=regions.end(); ++regions_it)
    {
        itk::ContinuousIndex<float, TImage::ImageDimension> center;
        center.Fill(0);

        typename RegionsMapType::mapped_type const & region = regions_it->second;
        for(typename RegionsMapType::mapped_type::const_iterator it=region.begin();
            it != region.end(); ++it)
        {
            IndexType const & index = *it;
            std::transform(center.Begin(), center.End(), index.m_Index,
                           center.Begin(), std::plus<float>());
        }

        std::transform(center.Begin(), center.End(),
            center.Begin(), std::bind2nd(std::divides<float>(), region.size()));

        // Convert to index : round the continuous index
        IndexType center_index;
        for(unsigned int i=0; i<TImage::ImageDimension; ++i)
        {
            center_index[i] = int(center[i]+0.5);
        }
        centers[regions_it->first] = center_index;
    }

    // Adjust the center of mass to the closest voxel in the region
    for(typename RegionsMapType::const_iterator regions_it=regions.begin();
        regions_it!=regions.end(); ++regions_it)
    {
        PixelType const & label = regions_it->first;

        IndexType & center = centers[label];
        if(this->m_Image->GetPixel(center) != label)
        {
            // Center of mass is not in region, find the closest voxel.

            // FIXME : function should always return true
            typedef AlwaysTrueImageFunction<ImageType> FunctionType;
            typename FunctionType::Pointer function = FunctionType::New();

            FloodFilledImageFunctionConditionalConstIterator<
                ImageType, FunctionType> image_it(this->m_Image, function, center);
            while(!image_it.IsAtEnd())
            {
                if(image_it.Get() == label)
                {
                    break;
                }
                ++image_it;
            }

            // Make sure we could find a matching voxel
            if(image_it.IsAtEnd())
            {
                throw std::runtime_error("Cannot find matching label");
            }

            // Adjust the center
            center = image_it.GetIndex();
        }
    }

    // Create the annotations
    for(typename RegionsMapType::const_iterator regions_it=regions.begin();
        regions_it!=regions.end(); ++regions_it)
    {
        PixelType const & label = regions_it->first;

        Annotation annotation;
        annotation.position = centers[label];
        // FIXME : compute annotation size
        annotation.size = 1;

        this->m_Annotations[label] = annotation;
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
::GetAnnotationPosition(typename ClustersToAnnotationsCalculator<TImage>::PixelType label) const
{
    return this->m_Annotations.find(label)->second.position;
}

template<typename TImage>
float
ClustersToAnnotationsCalculator<TImage>
::GetAnnotationSize(typename ClustersToAnnotationsCalculator<TImage>::PixelType label) const
{
    return this->m_Annotations.find(label)->second.size;
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
