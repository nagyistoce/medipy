/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _itkTensorReconstructionImageFilter_h
#define _itkTensorReconstructionImageFilter_h

#include <ostream>

#include <itkDataObject.h>
#include <itkImageToImageFilter.h>
#include <itkSmartPointer.h>

namespace itk
{

/**
 * @brief Base class to estimate the tensors image from diffusion data.
 * 
 * The tensors image is a VectorImage containing an estimation of the mean 
 * second-order tensor at each voxel of the input image; the baseline
 * image is a scalar image.
 * 
 * An optional mask image can be specified to restrict the estimation: only 
 * non-zero voxels in the mask image will be processed. The input image and the
 * mask image must be in the same physical space.
 */
template<typename TInputImage, typename TTensorsImage,
        typename TBaselineImage, typename TMaskImage>
class TensorReconstructionImageFilter: 
    public ImageToImageFilter<TInputImage, TTensorsImage>
{
public:
    /// @brief Standard class typedefs.
    typedef TensorReconstructionImageFilter Self;
    typedef ImageToImageFilter<TInputImage, TTensorsImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<Self const> ConstPointer;
    
    typedef TTensorsImage TensorsImageType;
    typedef typename TensorsImageType::Pointer TensorsImagePointer;
    
    typedef TBaselineImage BaselineImageType;
    typedef typename BaselineImageType::Pointer BaselineImagePointer;
    
    typedef TMaskImage MaskImageType;
    typedef typename MaskImageType::Pointer MaskImagePointer;
    typedef typename MaskImageType::ConstPointer MaskImageConstPointer;
    
    /// @brief Run-time type information (and related methods).
    itkTypeMacro(TensorReconstructionImageFilter, ImageToImageFilter);

    /// @brief Return the mask. 
    itkGetConstObjectMacro(MaskImage, MaskImageType);
    
    /// @brief Set the mask, default to NULL (no mask is used). 
    itkSetObjectMacro(MaskImage, MaskImageType);

    /// @brief Return the tensors image.
    TensorsImageType const * GetTensorsImage() const;
    
    /// @brief Return the baseline image.
    BaselineImageType const * GetBaselineImage() const;

protected:
    TensorReconstructionImageFilter();
    ~TensorReconstructionImageFilter() {}
    void PrintSelf(std::ostream& os, Indent indent) const;
    DataObject::Pointer MakeOutput(unsigned int index);
    void AllocateOutputs();
    
    /// @brief Return non-const pointer to the tensors image.
    TensorsImageType * GetTensorsImage();
    
    /// @brief Return non-const pointer to the baseline image.
    BaselineImageType * GetBaselineImage();

private:
    MaskImagePointer m_MaskImage;
    
    TensorReconstructionImageFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTensorReconstructionImageFilter.txx"
#endif

#endif // _itkTensorReconstructionImageFilter_h
