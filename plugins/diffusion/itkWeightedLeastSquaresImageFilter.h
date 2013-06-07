/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/


#ifndef _itkWeightedLeastSquaresImageFilter_h
#define _itkWeightedLeastSquaresImageFilter_h

#include <vector>

#include <itkImageToImageFilter.h>
#include <itkPoint.h>
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

template<typename TInputImage, typename TOutputImage>
class WeightedLeastSquaresImageFilter :
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public :
    /** Standard class typedefs. */
    typedef WeightedLeastSquaresImageFilter Self;
    typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<Self const> ConstPointer;

     typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(WeightedLeastSquaresImageFilter, ImageToImageFilter);

    /** Useful typedefs */
    typedef typename Superclass::InputImageType InputImageType;
    typedef typename Superclass::OutputImageType OutputImageType;
    typedef typename Superclass::InputImagePixelType InputImagePixelType;

    /** Intern types */
    typedef typename TOutputImage::PixelType OutputPixelType;
    typedef Point<float,3> DirectionType;
    typedef vnl_matrix<float> BMatrixType;

    /** Return the b-bvalue. */
    itkGetConstMacro(BVal, float);
    /** Set the b-value. */
    itkSetMacro(BVal, float);
    
    /** Return the nb of iter for WLS estimation. */
    itkGetConstMacro(IterationCount, unsigned int);
    /** Set the nb of iter for WLS estimation. */
    itkSetMacro(IterationCount, unsigned int);

    /** Return the number of gradient directions. */
    unsigned int GetNumberOfGradientDirections() const;
    /** Set the i^th gradient direction. */
    void SetGradientDirection(unsigned int i, DirectionType bvec);
    /** Return the i^th gradient direction. */
    DirectionType const & GetGradientDirection(unsigned int i) const;

protected :
    WeightedLeastSquaresImageFilter() {}
    ~WeightedLeastSquaresImageFilter() {}
    void PrintSelf(std::ostream& os, Indent indent) const;
    void AllocateOutputs();
    void BeforeThreadedGenerateData();
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, int);

private :
    std::vector<DirectionType> directions;
    float m_BVal;
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

