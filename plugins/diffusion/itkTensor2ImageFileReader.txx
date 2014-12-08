/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _11f8ba94_d03b_45f8_845e_19b8f9635285
#define _11f8ba94_d03b_45f8_845e_19b8f9635285

#include "itkTensor2ImageFileReader.h"

#include <algorithm>

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

namespace itk
{

template<typename TOutputImage>
void
Tensor2ImageFileReader<TOutputImage>
::SetFileName(std::string const & filename)
{
    if(this->m_TensorReader->GetFileName() != filename)
    {
        this->m_TensorReader->SetFileName(filename);
        this->Modified();
    }
}

template<typename TOutputImage>
std::string const &
Tensor2ImageFileReader<TOutputImage>
::GetFileName() const
{
    return this->m_TensorReader->GetFileName();
}

template<typename TOutputImage>
void
Tensor2ImageFileReader<TOutputImage>
::SetImageIO(ImageIOBase* imageIO)
{
    if(this->m_TensorReader->GetImageIO() != imageIO)
    {
        this->m_TensorReader->SetImageIO(imageIO);
        this->Modified();
    }
}

template<typename TOutputImage>
ImageIOBase*
Tensor2ImageFileReader<TOutputImage>
::GetImageIO()
{
    return this->m_TensorReader->GetImageIO();
}

template<typename TOutputImage>
void
Tensor2ImageFileReader<TOutputImage>
::GenerateOutputInformation()
{
    this->m_TensorReader->Update();

    OutputImagePointer output = this->GetOutput();
    output->CopyInformation(this->m_TensorReader->GetOutput());
    output->SetVectorLength(6);
}

template<typename TOutputImage>
void
Tensor2ImageFileReader<TOutputImage>
::EnlargeOutputRequestedRegion(DataObject* output)
{
    OutputImagePointer outputImage = dynamic_cast<TOutputImage*>(output);
    outputImage->SetRequestedRegion(this->m_TensorReader->GetOutput()->GetRequestedRegion());
}

template<typename TOutputImage>
Tensor2ImageFileReader<TOutputImage>
::Tensor2ImageFileReader()
: m_TensorReader(TensorReader::New())
{
    // Nothing else
}

template<typename TOutputImage>
Tensor2ImageFileReader<TOutputImage>
::~Tensor2ImageFileReader()
{
    // Nothing to do
}

template<typename TOutputImage>
void
Tensor2ImageFileReader<TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{

}

template<typename TOutputImage>
void
Tensor2ImageFileReader<TOutputImage>
::GenerateData()
{
    typename TOutputImage::Pointer output = this->GetOutput();
    this->AllocateOutputs();

    this->m_TensorReader->Update();
    this->tensorToUpperDiagonal(this->m_TensorReader->GetOutput(), output);
}

template<typename TOutputImage>
void
Tensor2ImageFileReader<TOutputImage>
::tensorToUpperDiagonal(TensorImagePointer tensors,
                        OutputImagePointer output) const
{
    typedef ImageRegionConstIterator<TensorImageType> TensorImageIterator;
    typedef ImageRegionIterator<OutputImageType> OutputImageIterator;

    TensorImageIterator tensorsIt(tensors, output->GetRequestedRegion());
    OutputImageIterator outputIt(output, output->GetRequestedRegion());

    while(!outputIt.IsAtEnd())
    {
        TensorType const & tensor = tensorsIt.Get();
        OutputImagePixelType vector = outputIt.Get();

        // itk::SymmetricSecondRankTensor stores the upper diagonal
        // (cf. operator()(unsigned int row, unsigned int col))
        std::copy(tensor.Begin(), tensor.End(), &vector[0]);

        ++tensorsIt;
        ++outputIt;
    }
}

}

#endif // _11f8ba94_d03b_45f8_845e_19b8f9635285
