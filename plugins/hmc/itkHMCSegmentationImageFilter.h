/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _7e14fdcd_5a32_4318_9136_2231262e46aa
#define _7e14fdcd_5a32_4318_9136_2231262e46aa

#include <vector>

#include <itkImage.h>
#include <itkImageToImageFilter.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include "HilbertCurveChainGenerator.h"

namespace itk
{

class HMCSegmentationImageFilter : 
    public ImageToImageFilter<Image<float, 3>, Image<float, 3> >
{
public:

    /** Standard class typedefs. */
    typedef HMCSegmentationImageFilter Self;
    typedef ImageToImageFilter<Image<float, 3>, Image<float, 3> > Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<Self const> ConstPointer;
    
    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(HMCSegmentationImageFilter, ImageToImageFilter);

    /** Type of the input image. */
    typedef Image<float, 3> InputImageType;

    typedef Image<float, 3> MaskImageType;
    typedef typename MaskImageType::Pointer MaskImagePointer;
    typedef typename MaskImageType::ConstPointer MaskImageConstPointer;

    /** Type of the output image. */
    typedef Image<float, 3> OutputImageType;
    typedef typename OutputImageType::IndexType OutputImageIndexType;
    
    /// @brief Return the mask. 
    itkGetConstObjectMacro(MaskImage, MaskImageType);
    
    /// @brief Set the mask, default to NULL (no mask is used). 
    itkSetObjectMacro(MaskImage, MaskImageType);
    
    /// @brief Return the index of the Flair image.
    itkGetConstMacro(FlairImage, int);
    
    /// @brief Set the index of the Flair image, default to -1 (no Flair image).
    itkSetMacro(FlairImage, int);
    
    /// @brief Return the maximum number of iterations.
    itkGetConstMacro(Iterations, unsigned int);
    
    /// @brief Set the maximum number of iterations, default to 5.
    itkSetMacro(Iterations, unsigned int);
    
    /// @brief Return the number of modalities.
    itkGetConstMacro(Modalities, unsigned int);
    
    /// @brief Set the number of modalities, default to 0.
    itkSetMacro(Modalities, unsigned int);
    
    itkGetConstMacro(DisplayOutliers, unsigned int);
    itkSetMacro(DisplayOutliers, unsigned int);
    itkBooleanMacro(DisplayOutliers)
    
    /// @brief Return the outliers criterion.
    itkGetConstMacro(OutliersCriterion, int);
    
    /** 
     * @brief Set the outliers criterion, default to 0 (keep a percentage of the
     * residuals), can also be 1 (use a probability threshold).
     */
    itkSetMacro(OutliersCriterion, int);
    
    /// @brief Return the threshold value for outliers.
    itkGetConstMacro(Threshold, float);
    
    /// @brief Set the threshold value for outliers, default to 0.
    itkSetMacro(Threshold, float);
    
    /// @brief Return the segmentation image.
    OutputImageType::Pointer GetSegmentationImage();
    
    /// @brief Return the outliers image.
    OutputImageType::Pointer GetOutliersImage();

protected:
    HMCSegmentationImageFilter();
    ~HMCSegmentationImageFilter();
    
    void GenerateData();
    DataObject::Pointer MakeOutput(unsigned int index);
    
private:
    MaskImageConstPointer m_MaskImage;    
    int m_FlairImage;
    unsigned int m_Iterations;
    unsigned int m_Modalities;
    bool m_DisplayOutliers;
    int m_OutliersCriterion;
    float m_Threshold;
    
    static void _chain_to_image(vnl_vector<int> const & chain, 
                                vnl_vector<int> const & chain_mask,
                                HilbertCurveChainGenerator::ScanConstPointer const & scan, 
                                OutputImageType::Pointer image);
};

}

#endif // _7e14fdcd_5a32_4318_9136_2231262e46aa
