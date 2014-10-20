/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _7b40aeac_4aca_4a7f_851d_e4d7399ca04d
#define _7b40aeac_4aca_4a7f_851d_e4d7399ca04d

#include "itkJointHistogramNormalizationFilter.h"

#include <fstream>

#include "itkJointHistogramTransferFunctionCalculator.h"

namespace itk
{

template<typename TInputImage, typename TMask, typename TOutputImage>
JointHistogramNormalizationFilter<TInputImage, TMask, TOutputImage>
::JointHistogramNormalizationFilter()
: m_BinsCountFixed(100), m_BinsCountMoving(100),
  m_Mask(0), m_MaskValue(1), 
  m_Method(Self::Method::NEAREST_NEIGHBOR)
{
    this->SetNumberOfRequiredInputs(2);
}

template<typename TInputImage, typename TMask, typename TOutputImage>
void
JointHistogramNormalizationFilter<TInputImage, TMask, TOutputImage>
::PrintSelf(std::ostream & os, Indent indent)
{
    this->Superclass::PrintSelf(os, indent);
    os << "Bins Count: " 
        << this->m_BinsCountFixed << ", " << this->m_BinsCountMoving << "\n";
    
    os << "Mask:";
    if(!this->m_Mask.IsNull())
    {
        os << "\n";
        this->m_Mask->Print(os, indent.GetNextIndent());
    }
    else
    {
        os << " (none)\n";
    }
    
    os << "Mask Value: " << this->m_MaskValue << "\n";
    os << "Method: " << this->m_Method << "\n";
}

template<typename TInputImage, typename TMask, typename TOutputImage>
void 
JointHistogramNormalizationFilter<TInputImage, TMask, TOutputImage>
::GenerateData()
{
    InputImageConstPointer fixed = this->GetFixedImage();
    InputImageConstPointer moving = this->GetMovingImage();
    
    typename JointHistogramCalculatorType::Pointer jh_calculator = 
        JointHistogramCalculatorType::New();
    jh_calculator->SetImage1(moving);
    jh_calculator->SetImage2(fixed);
    jh_calculator->SetMask(this->m_Mask);
    jh_calculator->SetMaskValue(this->m_MaskValue);
    jh_calculator->SetBinsCount1(this->m_BinsCountMoving);
    jh_calculator->SetBinsCount2(this->m_BinsCountFixed);
    jh_calculator->SetMethod(this->m_Method);

    jh_calculator->Compute();

    typedef typename JointHistogramCalculatorType::HistogramType HistogramType;
    
    typedef itk::JointHistogramTransferFunctionCalculator<HistogramType>
        JointHistogramTransferFunctionCalculatorType;
    typedef typename JointHistogramTransferFunctionCalculatorType::Pointer
        JointHistogramTransferFunctionCalculatorPointer;
        
    JointHistogramTransferFunctionCalculatorPointer tf_calculator = 
        JointHistogramTransferFunctionCalculatorType::New();
    tf_calculator->SetHistogram(jh_calculator->GetHistogram());
    tf_calculator->SetResampleFactor(10);
    
    tf_calculator->Compute();
    
    typedef typename JointHistogramTransferFunctionCalculatorType::TransferFunctionType
        TransferFunctionType;
    
    // vtkPiecewiseFunction::GetValue is not const...
    TransferFunctionType * transfer_function = 
        const_cast<TransferFunctionType *>(tf_calculator->GetTransferFunction());
    
    this->AllocateOutputs();
    typename OutputImageType::Pointer output = this->GetOutput();
    
    // Transfer function is from moving intensity to fixed intensity
    typedef itk::ImageRegionConstIterator<InputImageType> InputIteratorType;
    typedef itk::ImageRegionIterator<OutputImageType> OutputIteratorType;
    
    InputIteratorType input_it(moving, output->GetRequestedRegion());
    OutputIteratorType output_it(output, output->GetRequestedRegion());
    while(!output_it.IsAtEnd())
    {
        output_it.Set(transfer_function->GetValue(input_it.Get()));
        
        ++input_it;
        ++output_it;
    }
}

} // namespace itk

#endif // _7b40aeac_4aca_4a7f_851d_e4d7399ca04d
