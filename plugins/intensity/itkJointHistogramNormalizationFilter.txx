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
: m_BinsCount1(100), m_BinsCount2(100),
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
    os << "Bins Count: " << this->m_BinsCount1 << ", " << this->m_BinsCount2 << "\n";
    
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
    typename InputImageType::Pointer image1 = 
        const_cast<InputImageType *>(this->GetInput(0));
    typename InputImageType::Pointer image2 = 
        const_cast<InputImageType *>(this->GetInput(1));
    
    typename JointHistogramCalculatorType::Pointer jh_calculator = 
        JointHistogramCalculatorType::New();
    jh_calculator->SetImage1(image1);
    jh_calculator->SetImage2(image2);
    jh_calculator->SetMask(this->m_Mask);
    jh_calculator->SetMaskValue(this->m_MaskValue);
    jh_calculator->SetBinsCount1(this->m_BinsCount1);
    jh_calculator->SetBinsCount2(this->m_BinsCount2);
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
    typedef typename TransferFunctionType::Pointer TransferFunctionPointer;
    
    TransferFunctionPointer transfer_function = tf_calculator->GetTransferFunction();
    /*
    {
        typedef typename JointHistogramTransferFunctionCalculatorType::TransferFunctionType
            TransferFunctionType;
        std::ofstream stream("data.txt");
        typedef itk::ImageRegionConstIteratorWithIndex<TransferFunctionType> InputIteratorType;
        for(InputIteratorType it(transfer_function, transfer_function->GetRequestedRegion());
            !it.IsAtEnd(); ++it)
        {
            typename TransferFunctionType::PointType point;
            transfer_function->TransformIndexToPhysicalPoint(it.GetIndex(), point);
            stream << point[0] << " " << it.Get() << "\n";
        }
        stream << "\n";
    }
    */
    
    this->AllocateOutputs();
    typename OutputImageType::Pointer output = this->GetOutput();
    
    // Transfer function is from image1 intensity to image2 intensity
    typedef itk::ImageRegionConstIterator<InputImageType> InputIteratorType;
    typedef itk::ImageRegionIterator<OutputImageType> OutputIteratorType;
    
    InputIteratorType input_it(image1, output->GetRequestedRegion());
    OutputIteratorType output_it(output, output->GetRequestedRegion());
    while(!output_it.IsAtEnd())
    {
        typename TransferFunctionType::PointType point;
		point.Fill(input_it.Get());
		
		typename TransferFunctionType::IndexType index;
		transfer_function->TransformPhysicalPointToIndex(point, index);
        // Clamp inside of transfer function
        index[0] = std::min<typename TransferFunctionType::IndexValueType>(
            std::max<typename TransferFunctionType::IndexValueType>(index[0], 0),
            transfer_function->GetRequestedRegion().GetSize()[0]-1);
        
        output_it.Set(transfer_function->GetPixel(index));
        
        ++input_it;
        ++output_it;
    }
}

} // namespace itk

#endif // _7b40aeac_4aca_4a7f_851d_e4d7399ca04d
