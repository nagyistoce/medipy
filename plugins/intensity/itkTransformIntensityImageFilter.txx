/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _2d43251d_1b2b_40da_bbb0_c6201ebd3b12
#define _2d43251d_1b2b_40da_bbb0_c6201ebd3b12

#include "itkTransformIntensityImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

namespace itk
{

template<typename TInputImage, typename TOutputImage>
typename TransformIntensityImageFilter<TInputImage, TOutputImage>::TransferFunctionType const *
TransformIntensityImageFilter<TInputImage, TOutputImage>
::GetTransferFunction() const
{
    return this->m_TransferFunction;
}

template<typename TInputImage, typename TOutputImage>
void
TransformIntensityImageFilter<TInputImage, TOutputImage>
::SetTransferFunction(TransferFunctionType const * function)
{
    // vtkPiecewiseFunction::GetValue is not const...
    this->m_TransferFunction = const_cast<TransferFunctionType *>(function);
    this->Modified();
}

template<typename TInputImage, typename TOutputImage>
TransformIntensityImageFilter<TInputImage, TOutputImage>
::TransformIntensityImageFilter()
: m_TransferFunction(NULL)
{
    this->SetNumberOfRequiredInputs(1);
}

template<typename TInputImage, typename TOutputImage>
void
TransformIntensityImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream & os, Indent indent)
{
    this->Superclass::PrintSelf(os, indent);
    os << indent << "Transfer function:";
    this->m_TransferFunction->PrintSelf(os, vtkIndent());
}

template<typename TInputImage, typename TOutputImage>
void
TransformIntensityImageFilter<TInputImage, TOutputImage>
::ThreadedGenerateData(
    OutputImageRegionType const & outputRegionForThread, ThreadIdType)
{
    typedef itk::ImageRegionConstIterator<InputImageType> InputIteratorType;
    InputIteratorType input_iterator(this->GetInput(), outputRegionForThread);
    
    typedef itk::ImageRegionIterator<OutputImageType> OutputIteratorType;
    OutputIteratorType output_iterator(this->GetOutput(), outputRegionForThread);
    
    while(!output_iterator.IsAtEnd())
    {
        InputImagePixelType const input = input_iterator.Get();
        OutputImagePixelType const output = 
            this->m_TransferFunction->GetValue(input);
        
        output_iterator.Set(output);
        
        ++input_iterator;
        ++output_iterator;
    }
}

}

#endif // _2d43251d_1b2b_40da_bbb0_c6201ebd3b12
