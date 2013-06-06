/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _2dd8b74d_fbaf_4ee0_ac0b_e4e910c4f229
#define _2dd8b74d_fbaf_4ee0_ac0b_e4e910c4f229

#include <ostream>

#include <itkDataObject.h>
#include <itkImageToImageFilter.h>
#include <itkSmartPointer.h>

namespace itk
{

/**
 * @brief Base class to computate the mean and standard deviation of diffusion data.
 * 
 * The mean image is a VectorImage containing an estimation of the mean 
 * second-order tensor at each voxel of the input image; the standard deviation
 * image is a scalar image.
 * 
 * An optional mask image can be specified to restrict the estimation: only 
 * non-zero voxels in the mask image will be processed. The input image and the
 * mask image must be in the same physical space.
 */
template<typename TInputImage, typename TMeanImage, 
         typename TStandardDeviationImage, typename TMaskImage>
class DWIStatisticsImageFilter: 
    public ImageToImageFilter<TInputImage, TStandardDeviationImage>
{
public:
    /// @brief Standard class typedefs.
    typedef DWIStatisticsImageFilter Self;
    typedef ImageToImageFilter<TInputImage, TStandardDeviationImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<Self const> ConstPointer;
    
    typedef TMeanImage MeanImageType;
    typedef typename MeanImageType::Pointer MeanImagePointer;
    
    typedef TStandardDeviationImage StandardDeviationImageType;
    typedef typename StandardDeviationImageType::Pointer StandardDeviationImagePointer;
    
    typedef TMaskImage MaskImageType;
    typedef typename MaskImageType::Pointer MaskImagePointer;
    typedef typename MaskImageType::ConstPointer MaskImageConstPointer;
    
    /// @brief Run-time type information (and related methods).
    itkTypeMacro(DWIStatisticsImageFilter, ImageToImageFilter);

    /// @brief Return the mask. 
    itkGetConstObjectMacro(MaskImage, MaskImageType);
    
    /// @brief Set the mask, default to NULL (no mask is used). 
    itkSetObjectMacro(MaskImage, MaskImageType);

    /// @brief Return the mean image.
    MeanImageType const * GetMeanImage() const;
    
    /// @brief Return the standard deviation image.
    StandardDeviationImageType const * GetStandardDeviationImage() const;

protected:
    DWIStatisticsImageFilter();
    ~DWIStatisticsImageFilter() {}
    void PrintSelf(std::ostream& os, Indent indent) const;
    DataObject::Pointer MakeOutput(unsigned int index);
    void AllocateOutputs();

private:
    MaskImagePointer m_MaskImage;
    
    DWIStatisticsImageFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDWIStatisticsImageFilter.txx"
#endif

#endif // _2dd8b74d_fbaf_4ee0_ac0b_e4e910c4f229
