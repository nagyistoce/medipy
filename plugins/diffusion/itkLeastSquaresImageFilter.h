/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/


#ifndef _itkLeastSquaresImageFilter_h
#define _itkLeastSquaresImageFilter_h

#include <ostream>
#include <vector>

#include <itkTensorReconstructionImageFilter.h>
#include <itkPoint.h>

#include <vnl/vnl_matrix.h>

namespace itk
{

/**
 * \class LeastSquaresImageFilter
 * \brief Least Square Second Order Symmetric Tensor Reconstruction Filter
 * 
 * This filter must have as many input images as it has gradient directions.
 * All inputs are supposed to have the same Region.
 */

template<typename TInputImage, typename TTensorsImage,
        typename TBaselineImage, typename TMaskImage>
class LeastSquaresImageFilter :
    public TensorReconstructionImageFilter<TInputImage,
        TTensorsImage, TBaselineImage, TMaskImage>
{
public :
    /** Standard class typedefs. */
    typedef LeastSquaresImageFilter Self;
    typedef TensorReconstructionImageFilter<
        TInputImage, TTensorsImage, TBaselineImage,TMaskImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<Self const> ConstPointer;
    
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
    
    typedef TMaskImage MaskImageType;
    typedef typename MaskImageType::Pointer MaskImagePointer;
    typedef typename MaskImageType::ConstPointer MaskImageConstPointer;
    /// @brief Return the mask. 
    itkGetConstObjectMacro(MaskImage, MaskImageType);
    /// @brief Set the mask, default to NULL (no mask is used). 
    itkSetObjectMacro(MaskImage, MaskImageType);
    
    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(LeastSquaresImageFilter, TensorReconstructionImageFilter);

    /** Useful typedefs */
    typedef typename Superclass::InputImageType InputImageType;
    typedef typename Superclass::InputImagePixelType InputImagePixelType;

    typedef typename TTensorsImage::PixelType TensorsPixelType;
    typedef TTensorsImage TensorsImageType;
    typedef typename TensorsImageType::Pointer TensorsImagePointer;
    
    typedef TBaselineImage BaselineImageType;
    typedef typename TBaselineImage::PixelType BaselinePixelType;
    typedef typename BaselineImageType::Pointer BaselineImagePointer;

    /** Intern types */
    typedef Point<float,3> DirectionType;
    typedef float BValueType;
    typedef std::pair<BValueType, DirectionType> MetaDiffusionType;
    typedef vnl_matrix<float> BMatrixType;

    /** Return the number of gradient directions. */
    unsigned int GetNumberOfGradientDirections() const;
    /** Set the i^th gradient direction. */
    void SetGradientDirection(unsigned int i, DirectionType bvec);
    /** Return the i^th gradient direction. */
    DirectionType const & GetGradientDirection(unsigned int i) const;

    void SetBvalueAndGradientDirection(unsigned int i, BValueType BVal, DirectionType bvec);

protected :
    LeastSquaresImageFilter();
    ~LeastSquaresImageFilter() {}
    void PrintSelf(std::ostream& os, Indent indent) const;
    void BeforeThreadedGenerateData();
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, int);

private :
    std::vector<DirectionType> directions;
    std::vector<MetaDiffusionType> metadata_diffusion;
    
    BMatrixType bmatrix;
    BMatrixType invbmatrix;
    
    MaskImagePointer m_MaskImage;
    
    LeastSquaresImageFilter(Self const &); // purposely not implemented
    Self const & operator=(Self const &); // purposely not implemented

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLeastSquaresImageFilter.txx"
#endif

#endif 

