/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/


#ifndef _itkSymmetricSpectralAnalysisImageFilter_h
#define _itkSymmetricSpectralAnalysisImageFilter_h

#include <itkImageToImageFilter.h>
#include <itkSmartPointer.h>
#include <itkSymmetricEigenAnalysis.h>

namespace itk
{

/**
 * \class SymmetricSpectralAnalysisImageFilter
 * \brief Compute the eigenvalues and correspondind eigenvectors of a symmetric matrix
 * 
 */

template<typename TInputImage, typename TOutputImage>
class SymmetricSpectralAnalysisImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public :
    /** Standard class typedefs. */
    typedef SymmetricSpectralAnalysisImageFilter Self;
    typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<Self const> ConstPointer;

    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(SymmetricSpectralAnalysisImageFilter, ImageToImageFilter);

    /** Useful typedefs */
    typedef typename Superclass::InputImageType     InputImageType;
    typedef typename Superclass::OutputImageType    OutputImageType;
    typedef typename TOutputImage::PixelType        OutputPixelType;
    typedef typename TInputImage::PixelType         InputPixelType;
    typedef typename InputPixelType::ValueType      InputValueType;
    typedef typename OutputPixelType::ValueType     OutputValueType;

    typedef vnl_matrix<double> InputMatrixType;
    typedef FixedArray<double, 3> EigenValuesArrayType;
    typedef Matrix<double, 3, 3> EigenVectorMatrixType;
    typedef SymmetricEigenAnalysis<InputMatrixType, EigenValuesArrayType, EigenVectorMatrixType> CalculatorType;

    /** Intern types */
    enum EigenValueOrderType {
        OrderByValue,
        OrderByMagnitude,
        DoNotOrder
    };

    /** Accessors */
    itkGetMacro(SortOrder, EigenValueOrderType);
    itkSetMacro(SortOrder, EigenValueOrderType);
    itkStaticConstMacro(OutputImageDimension, unsigned int, TOutputImage::ImageDimension);
    
    /// @brief Return the eigenvalues image.
    TOutputImage * GetEigenValuesImage();
    /// @brief Return the eigenvectors image.
    TOutputImage * GetEigenVectorsImage();

protected :
    SymmetricSpectralAnalysisImageFilter();
    ~SymmetricSpectralAnalysisImageFilter() {}
    void PrintSelf(std::ostream& os, Indent indent) const;
    void AllocateOutputs();
    void BeforeThreadedGenerateData();
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, int);
    

private :
    EigenValueOrderType m_SortOrder;
    CalculatorType m_Calculator;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSymmetricSpectralAnalysisImageFilter.txx"
#endif

#endif 

