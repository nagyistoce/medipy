/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef d2922d4c_30b9_4fea_8e30_0e5e76af06e4
#define d2922d4c_30b9_4fea_8e30_0e5e76af06e4

#include <itkImageToImageFilter.h>

namespace itk
{

/**
 * @brief Clustering algorithm for change detection.
 */
template<typename TInputImage, typename TMaskImage, typename TOutputImage>
class ChangeDetectionClusteringImageFilter :
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public :
    /** Standard class typedefs. */
    typedef ChangeDetectionClusteringImageFilter Self;
    typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<Self const> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(ChangeDetectionClusteringImageFilter, ImageToImageFilter);

    /** Superclass typedefs. */
    typedef typename Superclass::InputImageType InputImageType;
    typedef typename Superclass::InputImagePointer InputImagePointer;
    typedef typename Superclass::InputImageConstPointer InputImageConstPointer;
    typedef typename Superclass::InputImageRegionType InputImageRegionType;
    typedef typename Superclass::InputImagePixelType InputImagePixelType;

    typedef typename Superclass::OutputImageType OutputImageType;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
    typedef typename Superclass::OutputImagePixelType  OutputImagePixelType;

    /** Type of the mask image. */
    typedef TMaskImage MaskImageType;
    typedef typename MaskImageType::PixelType MaskImagePixelType;

    /** Mask image. */
    itkGetConstObjectMacro(Mask, MaskImageType);
    itkSetObjectMacro(Mask, MaskImageType);

    /** Background value of the mask, default to 0. */
    itkGetMacro(MaskBackgroundValue, typename MaskImageType::PixelType);
    itkSetMacro(MaskBackgroundValue, typename MaskImageType::PixelType);

    /** Quantile for threshold selection. Defaults to 0.99. */
    itkGetMacro(Quantile, float);
    itkSetMacro(Quantile, float);

    /** Minimum cluster size, default to 0. */
    itkGetMacro(MinimumClusterSize, unsigned long);
    itkSetMacro(MinimumClusterSize, unsigned long);

    /** Minimum cluster size, default to std::numeric_limits<unsigned long>::max(). */
    itkGetMacro(MaximumNumberOfClusters, unsigned long);
    itkSetMacro(MaximumNumberOfClusters, unsigned long);

protected :
    ChangeDetectionClusteringImageFilter();
    ~ChangeDetectionClusteringImageFilter();
    void PrintSelf(std::ostream& os, Indent indent) const;

    void GenerateData();

private :
    typename MaskImageType::Pointer m_Mask;
    typename MaskImageType::PixelType m_MaskBackgroundValue;

    float m_Quantile;
    unsigned long m_MinimumClusterSize;
    unsigned long m_MaximumNumberOfClusters;

    ChangeDetectionClusteringImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    InputImagePixelType compute_threshold();
    
    template<typename TImage>
    void discard_clusters_close_to_mask_boundary(TImage * image);
    
    template<typename TImage>
    void discard_small_clusters(TImage * image);
};

}

#include "itkChangeDetectionClusteringImageFilter.txx"

#endif // d2922d4c_30b9_4fea_8e30_0e5e76af06e4
