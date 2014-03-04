/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _itkTensorReconstructionImageFilter_txx
#define _itkTensorReconstructionImageFilter_txx

#include "itkTensorReconstructionImageFilter.h"

#include <ostream>

namespace itk
{

template<typename TInputImage, typename TTensorsImage, typename TBaselineImage>
TTensorsImage const * 
TensorReconstructionImageFilter<TInputImage, TTensorsImage, TBaselineImage>
::GetTensorsImage() const
{
    return dynamic_cast<TTensorsImage const *>(this->ProcessObject::GetOutput(0));
}

template<typename TInputImage, typename TTensorsImage, typename TBaselineImage>
TBaselineImage const * 
TensorReconstructionImageFilter<TInputImage, TTensorsImage, TBaselineImage>
::GetBaselineImage() const
{
    return dynamic_cast<TBaselineImage const *>(this->ProcessObject::GetOutput(1));
}

template<typename TInputImage, typename TTensorsImage, typename TBaselineImage>
TensorReconstructionImageFilter<TInputImage, TTensorsImage, TBaselineImage>
::TensorReconstructionImageFilter()
{
    this->SetNumberOfRequiredOutputs(2);
    this->SetNthOutput(0,this->MakeOutput(0));
    this->SetNthOutput(1,this->MakeOutput(1));
}

template<typename TInputImage, typename TTensorsImage, typename TBaselineImage>
void 
TensorReconstructionImageFilter<TInputImage, TTensorsImage, TBaselineImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
    this->Superclass::PrintSelf(os, indent);
}

template<typename TInputImage, typename TTensorsImage, typename TBaselineImage>
DataObject::Pointer 
TensorReconstructionImageFilter<TInputImage, TTensorsImage, TBaselineImage>
::MakeOutput(unsigned int index)
{
    DataObject::Pointer output;
    
    if(index == 0)
    {
        output = TTensorsImage::New().GetPointer();
    }
    else if(index == 1)
    {
        output = TBaselineImage::New().GetPointer();
    }
    else
    {
        std::cerr << "No output " << index << std::endl;
        output = NULL;
    }
    
    return output.GetPointer();
}

template<typename TInputImage, typename TTensorsImage, typename TBaselineImage>
void 
TensorReconstructionImageFilter<TInputImage, TTensorsImage, TBaselineImage>
::AllocateOutputs()
{
    TensorsImagePointer TensorsPtr =
        dynamic_cast< TensorsImageType *>( this->ProcessObject::GetOutput(0) );
    if ( TensorsPtr )
    {
        TensorsPtr->SetBufferedRegion( TensorsPtr->GetRequestedRegion() );
        TensorsPtr->SetVectorLength(6);
        TensorsPtr->Allocate();
    }
    
    BaselineImagePointer BaselinePtr = 
        dynamic_cast< TBaselineImage *>( this->ProcessObject::GetOutput(1) );
    if(BaselinePtr)
    {
        BaselinePtr->SetBufferedRegion(BaselinePtr->GetRequestedRegion());
        BaselinePtr->Allocate();
    }
}

template<typename TInputImage, typename TTensorsImage, typename TBaselineImage>
TTensorsImage * 
TensorReconstructionImageFilter<TInputImage, TTensorsImage, TBaselineImage>
::GetTensorsImage() 
{
    return dynamic_cast<TTensorsImage *>(this->ProcessObject::GetOutput(0));
}

template<typename TInputImage, typename TTensorsImage, typename TBaselineImage>
TBaselineImage * 
TensorReconstructionImageFilter<TInputImage, TTensorsImage, TBaselineImage>
::GetBaselineImage()
{
    return dynamic_cast<TBaselineImage *>(this->ProcessObject::GetOutput(1));
}

}

#endif // _itkTensorReconstructionImageFilter_txx
