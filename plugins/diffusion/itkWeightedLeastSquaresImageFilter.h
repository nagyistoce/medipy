/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/


#ifndef _itkWeightedLeastSquaresImageFilter_h
#define _itkWeightedLeastSquaresImageFilter_h

#include <ostream>
#include <vector>

#include <itkTensorReconstructionImageFilter.h>
#include <itkVector.h>

#include <vnl/vnl_matrix.h>

namespace itk
{

/**
 * \class WeightedLeastSquaresImageFilter
 * \brief Weighted Least Squares Second Order Symmetric Tensor Reconstruction Filter
 * 
 * This filter must have as many input images as it has gradient directions.
 * All inputs are supposed to have the same Region.
 */

template<typename TInputImage, typename TTensorsImage, typename TBaselineImage>
class WeightedLeastSquaresImageFilter :
    public TensorReconstructionImageFilter<TInputImage, TTensorsImage, TBaselineImage>
{
public :
    /** Standard class typedefs. */
    typedef WeightedLeastSquaresImageFilter Self;
    typedef TensorReconstructionImageFilter<
        TInputImage, TTensorsImage, TBaselineImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<Self const> ConstPointer;
    
    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(WeightedLeastSquaresImageFilter, TensorReconstructionImageFilter);

    /* Typedefs from superclass */
    typedef typename Superclass::InputImageType InputImageType;
    typedef typename Superclass::TensorsImageType TensorsImageType;
    typedef typename Superclass::TensorsImagePointer TensorsImagePointer;
    typedef typename Superclass::BaselineImageType BaselineImageType;
    typedef typename Superclass::BaselineImagePointer BaselineImagePointer;

    /** Gradient direction type */
    typedef Vector<float,3> DirectionType;
    /** B-value type */
    typedef float BValueType;
    
    /** Return the nb of iter for WLS estimation. */
    itkGetConstMacro(IterationCount, unsigned int);
    /** Set the nb of iter for WLS estimation. */
    itkSetMacro(IterationCount, unsigned int);

    /** Return the number of gradient directions. */
    unsigned int GetNumberOfGradientDirections() const;
    /** Return the i^th gradient direction. */
    DirectionType const & GetGradientDirection(unsigned int i) const;
    /** Return the i^th b-value. */
    BValueType GetBvalue(unsigned int i) const;
    /** Set the i^th b-value and gradient direction. */
    void SetBvalueAndGradientDirection(unsigned int i, BValueType BVal, DirectionType bvec);

protected :
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
    
    WeightedLeastSquaresImageFilter() {}
    ~WeightedLeastSquaresImageFilter() {}
    void PrintSelf(std::ostream& os, Indent indent) const;
    void BeforeThreadedGenerateData();
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, int);

private :
    typedef vnl_matrix<float> BMatrixType;
    typedef std::pair<BValueType, DirectionType> MetaDiffusionType;
    
    std::vector<MetaDiffusionType> metadata_diffusion;
    
    unsigned int m_IterationCount;
    
    BMatrixType bmatrix;
    BMatrixType invbmatrix;

    WeightedLeastSquaresImageFilter(Self const &); // purposely not implemented
    Self const & operator=(Self const &); // purposely not implemented

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkWeightedLeastSquaresImageFilter.txx"
#endif

#endif 
