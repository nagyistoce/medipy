/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/


#ifndef _itkSecondOrderSymmetricTensorReconstructionFilter_h
#define _itkSecondOrderSymmetricTensorReconstructionFilter_h

#include <vector>

#include <itkImageToImageFilter.h>
#include <itkImageSource.h>
#include <itkSmartPointer.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_vector_fixed_ref.h>

namespace itk
{

/**
 * \class SecondOrderSymmetricTensorReconstructionFilter
 * \brief Least Square Second Order Symmetric Tensor Reconstruction Filter
 * 
 */

template<typename TInputImage, typename TOutputImage>
class SecondOrderSymmetricTensorReconstructionFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public :
    /** Standard class typedefs. */
    typedef SecondOrderSymmetricTensorReconstructionFilter Self;
    typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<Self const> ConstPointer;

     typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(SecondOrderSymmetricTensorReconstructionFilter, ImageToImageFilter);

    /** Useful typedefs */
    typedef typename Superclass::InputImageType InputImageType;
    typedef typename Superclass::OutputImageType OutputImageType;
    typedef typename Superclass::InputImagePixelType InputImagePixelType;

    /** Intern types */
    typedef typename TOutputImage::PixelType OutputPixelType;
    typedef Point<float,3> DirectionType;
    typedef vnl_matrix<float> BMatrixType;

    /** Accessors */
    itkGetMacro(BVal, float);
    itkSetMacro(BVal, float);
    itkStaticConstMacro(OutputImageDimension, unsigned int, TOutputImage::ImageDimension);

    void SetGradientDirection(unsigned int i, DirectionType bvec);
    DirectionType GetGradientDirection(unsigned int i);

protected :
    SecondOrderSymmetricTensorReconstructionFilter() {}
    ~SecondOrderSymmetricTensorReconstructionFilter() {}
    void PrintSelf(std::ostream& os, Indent indent) const;
    void AllocateOutputs();
    void BeforeThreadedGenerateData();
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, int);

private :
    std::vector<DirectionType> directions;
    float m_BVal;
    BMatrixType bmatrix;
    BMatrixType invbmatrix;

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSecondOrderSymmetricTensorReconstructionFilter.txx"
#endif

#endif 

