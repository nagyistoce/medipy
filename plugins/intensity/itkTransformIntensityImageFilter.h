/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _8091ee4c_bc6b_44b4_af97_1431f37978ea
#define _8091ee4c_bc6b_44b4_af97_1431f37978ea

#include <itkImageToImageFilter.h>
#include <vtkPiecewiseFunction.h>
#include <vtkSmartPointer.h>

namespace itk
{

template<typename TInputImage, typename TOutputImage>
class TransformIntensityImageFilter: 
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
    /** Standard typedefs */
    typedef TransformIntensityImageFilter Self;
    typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<Self const> ConstPointer;

    /** standard New() method support */
    itkNewMacro(Self);
    
    /** Run-time type information (and related methods). */
    itkTypeMacro(TransformIntensityImageFilter, ImageToImageFilter);
    
    typedef typename Superclass::InputImageType InputImageType;
    typedef typename Superclass::InputImagePointer InputImagePointer;
    typedef typename Superclass::InputImageConstPointer InputImageConstPointer;
    typedef typename Superclass::InputImageRegionType InputImageRegionType;
    typedef typename Superclass::InputImagePixelType InputImagePixelType;
    
    typedef typename Superclass::OutputImageType OutputImageType;
    typedef typename Superclass::OutputImagePointer OutputImagePointer;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
    typedef typename Superclass::OutputImagePixelType OutputImagePixelType;
    
    typedef vtkPiecewiseFunction TransferFunctionType;
    
    TransferFunctionType const * GetTransferFunction() const;
    void SetTransferFunction(TransferFunctionType const * function);

protected:
    typedef vtkSmartPointer<TransferFunctionType> TransferFunctionPointer;
    
    TransformIntensityImageFilter();
    void PrintSelf(std::ostream & os, Indent indent);
    void ThreadedGenerateData(
        OutputImageRegionType const & outputRegionForThread, int);
    
    TransferFunctionPointer m_TransferFunction;

private:
    TransformIntensityImageFilter(Self const &); // purposely not implemented
    Self const & operator=(Self const &); // purposely not implemented
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTransformIntensityImageFilter.txx"
#endif

#endif // _8091ee4c_bc6b_44b4_af97_1431f37978ea
