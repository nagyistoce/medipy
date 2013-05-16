/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _881946a8_b5ce_4de4_98e8_2ffc0cc32092
#define _881946a8_b5ce_4de4_98e8_2ffc0cc32092

#include "itkDWIStatisticsImageFilter.h"

#include <ostream>

namespace itk
{

template<typename TInputImage, typename TMeanImage, typename TStandardDeviationImage, typename TMaskImage>
TMeanImage const * 
DWIStatisticsImageFilter<TInputImage, TMeanImage, TStandardDeviationImage, TMaskImage>
::GetMeanImage() const
{
    return dynamic_cast<TMeanImage const *>(this->ProcessObject::GetOutput(0));
}

template<typename TInputImage, typename TMeanImage, typename TStandardDeviationImage, typename TMaskImage>
TStandardDeviationImage const * 
DWIStatisticsImageFilter<TInputImage, TMeanImage, TStandardDeviationImage, TMaskImage>
::GetStandardDeviationImage() const
{
    return dynamic_cast<TStandardDeviationImage const *>(
        this->ProcessObject::GetOutput(1));
}

template<typename TInputImage, typename TMeanImage, typename TStandardDeviationImage, typename TMaskImage>
DWIStatisticsImageFilter<TInputImage, TMeanImage, TStandardDeviationImage, TMaskImage>
::DWIStatisticsImageFilter()
{
    this->SetNumberOfRequiredOutputs(2);
    this->SetNthOutput(0,this->MakeOutput(0));
    this->SetNthOutput(1,this->MakeOutput(1));
    
    this->SetMaskImage(NULL);
}

template<typename TInputImage, typename TMeanImage, typename TStandardDeviationImage, typename TMaskImage>
void 
DWIStatisticsImageFilter<TInputImage, TMeanImage, TStandardDeviationImage, TMaskImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
    this->Superclass::PrintSelf(os, indent);
    os << indent << "Mask image: \n";
    if(this->m_MaskImage.IsNull())
    {
        os << indent.GetNextIndent() << "None\n";
    }
    else
    {
        this->m_MaskImage->Print(os, indent.GetNextIndent());
    }
}

template<typename TInputImage, typename TMeanImage, typename TStandardDeviationImage, typename TMaskImage>
DataObject::Pointer 
DWIStatisticsImageFilter<TInputImage, TMeanImage, TStandardDeviationImage, TMaskImage>
::MakeOutput(unsigned int index)
{
    DataObject::Pointer output;
    
    if(index == 0)
    {
        output = TMeanImage::New().GetPointer();
    }
    else if(index == 1)
    {
        output = TStandardDeviationImage::New().GetPointer();
    }
    else
    {
        std::cerr << "No output " << index << std::endl;
        output = NULL;
    }
    
    return output.GetPointer();
}

template<typename TInputImage, typename TMeanImage, typename TStandardDeviationImage, typename TMaskImage>
void 
DWIStatisticsImageFilter<TInputImage, TMeanImage, TStandardDeviationImage, TMaskImage>
::AllocateOutputs()
{
    typename TMeanImage::Pointer mean_image = 
        dynamic_cast<TMeanImage*>(this->ProcessObject::GetOutput(0));
    if(mean_image) 
    {
        mean_image->SetBufferedRegion(mean_image->GetRequestedRegion());
        mean_image->SetVectorLength(6);
        mean_image->Allocate();
    }
    
    typename TStandardDeviationImage::Pointer standard_deviation_image = 
        dynamic_cast<TStandardDeviationImage*>(this->ProcessObject::GetOutput(1));
    if(standard_deviation_image)
    {
        standard_deviation_image->SetBufferedRegion(standard_deviation_image->GetRequestedRegion());
        standard_deviation_image->Allocate();
    }
}

}

#endif // _881946a8_b5ce_4de4_98e8_2ffc0cc32092
