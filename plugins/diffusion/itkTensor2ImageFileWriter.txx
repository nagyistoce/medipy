/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _25aa7409_0157_424b_9c22_7e2e20f2712e
#define _25aa7409_0157_424b_9c22_7e2e20f2712e

#include "itkTensor2ImageFileWriter.h"

#include <algorithm>

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

namespace itk
{

template<typename TInputImage>
void
Tensor2ImageFileWriter<TInputImage>
::SetInput(InputImageType const * input)
{
    // ProcessObject is not const_correct so this cast is required here.
    this->ProcessObject::SetNthInput(0, const_cast<TInputImage *>(input));
}

template<typename TInputImage>
typename Tensor2ImageFileWriter<TInputImage>::InputImageType const *
Tensor2ImageFileWriter<TInputImage>
::GetInput()
{
    if(this->GetNumberOfInputs() < 1)
    {
        return 0;
    }

    return static_cast<TInputImage*>(this->ProcessObject::GetInput(0));
}

template<typename TInputImage>
typename Tensor2ImageFileWriter<TInputImage>::InputImageType const *
Tensor2ImageFileWriter<TInputImage>
::GetInput(unsigned int idx)
{
    return static_cast<TInputImage*>(this->ProcessObject::GetInput(idx));
}

template<typename TInputImage>
void
Tensor2ImageFileWriter<TInputImage>
::SetFileName(char const* filename)
{
    if(this->m_TensorWriter->GetFileName() != filename)
    {
        this->m_TensorWriter->SetFileName(filename);
        this->Modified();
    }
}

template<typename TInputImage>
void
Tensor2ImageFileWriter<TInputImage>
::SetFileName(std::string const & filename)
{
    this->SetFileName(filename.c_str());
}

template<typename TInputImage>
const char*
Tensor2ImageFileWriter<TInputImage>
::GetFileName() const
{
    return this->m_TensorWriter->GetFileName();
}

template<typename TInputImage>
void
Tensor2ImageFileWriter<TInputImage>
::SetImageIO(ImageIOBase* imageIO)
{
    if(this->m_TensorWriter->GetImageIO() != imageIO)
    {
        this->m_TensorWriter->SetImageIO(imageIO);
        this->Modified();
    }
}

template<typename TInputImage>
ImageIOBase*
Tensor2ImageFileWriter<TInputImage>
::GetImageIO()
{
    return this->m_TensorWriter->GetImageIO();
}

template<typename TInputImage>
void
Tensor2ImageFileWriter<TInputImage>
::Write()
{
    InputImageType const * input = this->GetInput();

    TensorImagePointer tensors = TensorImageType::New();
    tensors->CopyInformation(input);
    tensors->SetRegions(input->GetRequestedRegion());
    tensors->Allocate();

    this->upperDiagonalToTensor(input, tensors);
    this->m_TensorWriter->SetInput(tensors);
    this->m_TensorWriter->Write();
}

template<typename TInputImage>
void
Tensor2ImageFileWriter<TInputImage>
::SetIORegion(ImageIORegion const & region)
{
    if(this->m_TensorWriter->GetIORegion() != region)
    {
        this->m_TensorWriter->SetIORegion(region);
        this->Modified();
    }
}

template<typename TInputImage>
ImageIORegion const &
Tensor2ImageFileWriter<TInputImage>
::GetIORegion() const
{
    return this->m_TensorWriter->GetIORegion();
}

template<typename TInputImage>
void
Tensor2ImageFileWriter<TInputImage>
::Update()
{
    this->Write();
}

template<typename TInputImage>
Tensor2ImageFileWriter<TInputImage>
::Tensor2ImageFileWriter()
: m_TensorWriter(TensorWriter::New())
{
    // Nothing else
}

template<typename TInputImage>
Tensor2ImageFileWriter<TInputImage>
::~Tensor2ImageFileWriter()
{
    // Nothing to do
}

template<typename TInputImage>
void
Tensor2ImageFileWriter<TInputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{

}

template<typename TInputImage>
void
Tensor2ImageFileWriter<TInputImage>
::upperDiagonalToTensor(InputImageConstPointer input,
                        TensorImagePointer tensors) const
{
    typedef ImageRegionConstIterator<InputImageType> InputImageIterator;
    typedef ImageRegionIterator<TensorImageType> TensorImageIterator;

    InputImageIterator inputIt(input, tensors->GetRequestedRegion());
    TensorImageIterator tensorsIt(tensors, tensors->GetRequestedRegion());

    while(!tensorsIt.IsAtEnd())
    {
        InputImagePixelType vector = inputIt.Get();
        TensorType tensor;

        // itk::SymmetricSecondRankTensor stores the upper diagonal
        // (cf. operator()(unsigned int row, unsigned int col))
        std::copy(&vector[0], &vector[6], tensor.Begin());

        tensorsIt.Set(tensor);

        ++inputIt;
        ++tensorsIt;
    }
}

}

#endif // _25aa7409_0157_424b_9c22_7e2e20f2712e

