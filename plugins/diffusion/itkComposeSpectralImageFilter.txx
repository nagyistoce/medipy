/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _ComposeSpectralImageFilter_txx
#define _ComposeSpectralImageFilter_txx

#include "itkComposeSpectralImageFilter.h"

#include <itkImageRegionIterator.h>

namespace itk
{

template<typename TInputImage, typename TOutputImage>
ComposeSpectralImageFilter<TInputImage, TOutputImage>
::ComposeSpectralImageFilter()
{
    // Duplicate the code from itk::ImageSource since we do not derive
    // from it
    this->SetNumberOfRequiredInputs(2);
    
    typename TOutputImage::Pointer output =
        static_cast<TOutputImage*>(this->MakeOutput(0).GetPointer()); 
    this->ProcessObject::SetNumberOfRequiredOutputs(1);
    this->ProcessObject::SetNthOutput(0, output.GetPointer());
}

template<typename TInputImage, typename TOutputImage>
void
ComposeSpectralImageFilter<TInputImage, TOutputImage>
::SetInput(unsigned int idx, InputImageType const * input)
{
    // ProcessObject is not const_correct so this cast is required here.
    this->ProcessObject::SetNthInput(idx, const_cast<TInputImage *>(input));
}

template<typename TInputImage, typename TOutputImage>
void
ComposeSpectralImageFilter<TInputImage, TOutputImage>
::SetEigenValues(InputImageType const * input)
{
    // ProcessObject is not const_correct so this cast is required here.
    this->ProcessObject::SetNthInput(0, const_cast<TInputImage *>(input));
}

template<typename TInputImage, typename TOutputImage>
void
ComposeSpectralImageFilter<TInputImage, TOutputImage>
::SetEigenVectors(InputImageType const * input)
{
    // ProcessObject is not const_correct so this cast is required here.
    this->ProcessObject::SetNthInput(1, const_cast<TInputImage *>(input));
}

template<typename TInputImage, typename TOutputImage>
typename ComposeSpectralImageFilter<TInputImage, TOutputImage>::InputImageType const *
ComposeSpectralImageFilter<TInputImage, TOutputImage>
::GetInput()
{
    if(this->GetNumberOfInputs() < 1)
    {
        return 0;
    }
    return static_cast<TInputImage*>(this->ProcessObject::GetInput(0));
}

template<typename TInputImage, typename TOutputImage>
typename ComposeSpectralImageFilter<TInputImage, TOutputImage>::InputImageType const *
ComposeSpectralImageFilter<TInputImage, TOutputImage>
::GetInput(unsigned int idx)
{
    return static_cast<TInputImage*>(this->ProcessObject::GetInput(idx));
}

template<typename TInputImage, typename TOutputImage>
void
ComposeSpectralImageFilter<TInputImage, TOutputImage>
::AllocateOutputs()
{
    OutputImagePointer TensorsPtr = this->GetOutput();
    if( TensorsPtr )
    {
        // Explicitely copy information from input since we do not
        // derive from ImageToImageFilter
        TensorsPtr->SetRegions( this->GetInput(0)->GetRequestedRegion() );
        TensorsPtr->SetVectorLength(6);
        TensorsPtr->Allocate();
    }
}

template<typename TInputImage, typename TOutputImage>
TInputImage * 
ComposeSpectralImageFilter<TInputImage, TOutputImage>
::GetEigenValuesImage() 
{
    return dynamic_cast<TInputImage *>(this->ProcessObject::GetInput(0));
}

template<typename TInputImage, typename TOutputImage>
TInputImage * 
ComposeSpectralImageFilter<TInputImage, TOutputImage>
::GetEigenVectorsImage()
{
    return dynamic_cast<TInputImage *>(this->ProcessObject::GetInput(1));
}

template<typename TInputImage, typename TOutputImage >
void
ComposeSpectralImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
    Superclass::PrintSelf(os,indent);
    
    // Print Tensor Image
    os << indent << "Tensor image: \n";
    if(this->GetOutput() == NULL)
    {
        os << indent.GetNextIndent() << "None\n";
    }
    else
    {
        this->GetOutput()->Print(os, indent.GetNextIndent());
    }
}

template<typename TInputImage, typename TOutputImage>
typename ComposeSpectralImageFilter<TInputImage, TOutputImage>::DataObjectPointer
ComposeSpectralImageFilter<TInputImage, TOutputImage>
::MakeOutput(unsigned int)
{
    // Explicitely copy information from input since we do not
    // derive from ImageToImageFilter
    return static_cast<DataObject*>(TOutputImage::New().GetPointer());
}

template<typename TInputImage, typename TOutputImage>
typename ComposeSpectralImageFilter<TInputImage, TOutputImage>::OutputImageType * 
ComposeSpectralImageFilter<TInputImage, TOutputImage>
::GetOutput()
{
    // Explicitely copy information from input since we do not
    // derive from ImageToImageFilter
    return static_cast<TOutputImage*>(this->ProcessObject::GetOutput(0));
}

template<typename TInputImage, typename TOutputImage>
typename ComposeSpectralImageFilter<TInputImage, TOutputImage>::OutputImageType const *
ComposeSpectralImageFilter<TInputImage, TOutputImage>
::GetOutput() const
{
    // Explicitely copy information from input since we do not
    // derive from ImageToImageFilter
    return static_cast<TOutputImage const*>(this->ProcessObject::GetOutput(0));
}

template<typename TInputImage, typename TOutputImage>
void
ComposeSpectralImageFilter<TInputImage, TOutputImage>
::Update()
{
    // This class does not derive from itk::ImageSource, call
    // AllocateOutputs explicitely
    this->AllocateOutputs();
    
    // Input images
    InputImagePointer eigenvaluesImage = this->GetEigenValuesImage();
    InputImagePointer eigenvectorsImage = this->GetEigenVectorsImage();
    
    typedef ImageRegionIterator<InputImageType> EigenValueIterator;
    EigenValueIterator EigenValueIt(eigenvaluesImage, eigenvaluesImage->GetRequestedRegion());
    EigenValueIt.GoToBegin();
    
    typedef ImageRegionIterator<InputImageType> EigenVectorIterator;
    EigenVectorIterator EigenVectorIt(eigenvectorsImage, eigenvectorsImage->GetRequestedRegion());
    EigenVectorIt.GoToBegin();
    
    // Output image
    OutputImagePointer output = this->GetOutput();
    
    typedef ImageRegionIterator<OutputImageType> outputIterator;
    outputIterator outputIt(output, output->GetRequestedRegion());
    outputIt.GoToBegin();
    
    
    while(!outputIt.IsAtEnd())
    {
        InputImagePixelType eigenvalues = EigenValueIt.Get();
        InputImagePixelType eigenvectors = EigenVectorIt.Get();
        OutputImagePixelType tensor = outputIt.Get();
        
        // Reconstruction of the tensor
        OutputImagePixelType dt6_new(6);
        dt6_new.Fill(0);
        
        dt6_new[0] = eigenvectors[0]*eigenvalues[0]*eigenvectors[0]
                    + eigenvectors[3]*eigenvalues[1]*eigenvectors[3]
                    + eigenvectors[6]*eigenvalues[2]*eigenvectors[6];
        dt6_new[1] = eigenvectors[0]*eigenvalues[0]*eigenvectors[1]
                    + eigenvectors[3]*eigenvalues[1]*eigenvectors[4]
                    + eigenvectors[6]*eigenvalues[2]*eigenvectors[7];
        dt6_new[2] = eigenvectors[0]*eigenvalues[0]*eigenvectors[2]
                    + eigenvectors[3]*eigenvalues[1]*eigenvectors[5]
                    + eigenvectors[6]*eigenvalues[2]*eigenvectors[8];
        dt6_new[3] = eigenvectors[1]*eigenvalues[0]*eigenvectors[1]
                    + eigenvectors[4]*eigenvalues[1]*eigenvectors[4]
                    + eigenvectors[7]*eigenvalues[2]*eigenvectors[7];
        dt6_new[4] = eigenvectors[1]*eigenvalues[0]*eigenvectors[2]
                    + eigenvectors[4]*eigenvalues[1]*eigenvectors[5]
                    + eigenvectors[7]*eigenvalues[2]*eigenvectors[8];
        dt6_new[5] = eigenvectors[2]*eigenvalues[0]*eigenvectors[2]
                    + eigenvectors[5]*eigenvalues[1]*eigenvectors[5]
                    + eigenvectors[8]*eigenvalues[2]*eigenvectors[8];
        
        std::copy(&dt6_new[0], &dt6_new[dt6_new.Size()-1], &tensor[0]);

        outputIt.Set(tensor);

        ++EigenValueIt;
        ++EigenVectorIt;
        ++outputIt;
    }
}

template<typename TInputImage, typename TOutputImage>
TInputImage const * 
ComposeSpectralImageFilter<TInputImage, TOutputImage>
::GetEigenValuesImage() const
{
    return dynamic_cast<TInputImage const *>(this->ProcessObject::GetInput(0));
}

template<typename TInputImage, typename TOutputImage>
TInputImage const * 
ComposeSpectralImageFilter<TInputImage, TOutputImage>
::GetEigenVectorsImage() const
{
    return dynamic_cast<TInputImage const *>(this->ProcessObject::GetInput(1));
}

}

#endif // _ComposeSpectralImageFilter_txx
