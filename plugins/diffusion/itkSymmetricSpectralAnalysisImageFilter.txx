/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _SymmetricSpectralAnalysisImageFilter_txx
#define _SymmetricSpectralAnalysisImageFilter_txx

#include "itkSymmetricSpectralAnalysisImageFilter.h"

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkVectorImage.h>

namespace itk
{

template<typename TInputImage, typename TOutputImage>
TOutputImage * 
SymmetricSpectralAnalysisImageFilter<TInputImage, TOutputImage>
::GetEigenValuesImage() 
{
    return dynamic_cast<TOutputImage *>(this->ProcessObject::GetOutput(0));
}

template<typename TInputImage, typename TOutputImage>
TOutputImage * 
SymmetricSpectralAnalysisImageFilter<TInputImage, TOutputImage>
::GetEigenVectorsImage() 
{
    return dynamic_cast<TOutputImage *>(this->ProcessObject::GetOutput(1));
}

template<typename TInputImage, typename TOutputImage>
SymmetricSpectralAnalysisImageFilter<TInputImage, TOutputImage>
::SymmetricSpectralAnalysisImageFilter()
{
    this->SetNumberOfRequiredOutputs(2);
    this->SetNthOutput(0,this->MakeOutput(0));
    this->SetNthOutput(1,this->MakeOutput(1));

    this->m_Calculator.SetDimension(3);
    
    this->m_SortOrder = OrderByValue;
}

template<typename TInputImage, typename TOutputImage>
void
SymmetricSpectralAnalysisImageFilter<TInputImage, TOutputImage>
::AllocateOutputs()
{
    typename OutputImageType::Pointer outputPtr_val;
    typename OutputImageType::Pointer outputPtr_vec;

    outputPtr_val = this->GetEigenValuesImage();
    outputPtr_vec = this->GetEigenVectorsImage();

    if ( outputPtr_val ) {
        outputPtr_val->SetBufferedRegion( outputPtr_val->GetRequestedRegion() );
        outputPtr_val->SetVectorLength(3);
        outputPtr_val->Allocate();
    }
    if ( outputPtr_vec ) {
        outputPtr_vec->SetBufferedRegion( outputPtr_vec->GetRequestedRegion() );
        outputPtr_vec->SetVectorLength(9);
        outputPtr_vec->Allocate();
    }
}

template<typename TInputImage, typename TOutputImage>
void
SymmetricSpectralAnalysisImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
    Superclass::PrintSelf(os,indent);
}

template<typename TInputImage, typename TOutputImage>
void 
SymmetricSpectralAnalysisImageFilter<TInputImage, TOutputImage>
::BeforeThreadedGenerateData()
{
    if( this->m_SortOrder == OrderByMagnitude ) {
        this->m_Calculator.SetOrderEigenMagnitudes( true );
    }
    else if( this->m_SortOrder == DoNotOrder ) {
        this->m_Calculator.SetOrderEigenValues( false );
    }
}

template<typename TInputImage, typename TOutputImage>
void 
SymmetricSpectralAnalysisImageFilter<TInputImage, TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, ThreadIdType )
{
    //Const Iterator on Input 
    InputImageType const * tensors = this->GetInput(0);
    
    typedef ImageRegionConstIterator<InputImageType> InputIterator;
    InputIterator inputIt(tensors, outputRegionForThread);
    inputIt.GoToBegin();
    
    //Iterator on the two Output
    typename OutputImageType::Pointer output_val = this->GetEigenValuesImage();
    typename OutputImageType::Pointer output_vec = this->GetEigenVectorsImage();
    
    typedef ImageRegionIterator< OutputImageType > OutputIterator;
    OutputIterator outputValIt(output_val, outputRegionForThread);
    outputValIt.GoToBegin();
    OutputIterator outputVecIt(output_vec, outputRegionForThread);
    outputVecIt.GoToBegin();
    
    while( !inputIt.IsAtEnd() )
    {
        InputPixelType dt6 = inputIt.Get();
        InputMatrixType dt33(3,3);
        dt33[0][0] =  (double) dt6[0];
        dt33[1][1] =  (double) dt6[3];
        dt33[2][2] =  (double) dt6[5];
        dt33[0][1] =  (double) dt6[1];
        dt33[0][2] =  (double) dt6[2];
        dt33[1][2] =  (double) dt6[4];
        dt33[1][0] =  dt33[0][1];
        dt33[2][0] =  dt33[0][2];
        dt33[2][1] =  dt33[1][2];

        EigenValuesArrayType eigenvalues;
        EigenVectorMatrixType eigenvectors;

        this->m_Calculator.ComputeEigenValuesAndVectors(dt33,eigenvalues,eigenvectors);

        OutputPixelType val = outputValIt.Get();
        for (unsigned int i=0; i<3; i++)
        {
            val[i] = (OutputValueType) eigenvalues[i];
        }
        
        OutputPixelType vec = outputVecIt.Get();
        for (unsigned int i=0; i<3; i++)
        {
            for (unsigned int j=0; j<3; j++)
            {
                vec[3*i+j] = (OutputValueType) eigenvectors[i][j]; 
            }
        }
        
        ++inputIt;
        ++outputValIt;
        ++outputVecIt;
    }

}

}

#endif
