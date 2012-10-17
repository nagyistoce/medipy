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

#include "itkVectorImage.h"
#include "itkImageRegionIteratorWithIndex.h"

namespace itk
{


template<typename TInputImage, typename TOutputImage>
SymmetricSpectralAnalysisImageFilter<TInputImage, TOutputImage>
::SymmetricSpectralAnalysisImageFilter()
{
    this->SetNumberOfRequiredOutputs(2);
    this->SetNthOutput(0,this->MakeOutput(0));
    this->SetNthOutput(1,this->MakeOutput(1));

    this->m_Calculator.SetDimension(3);

    this->m_Order = OrderByValue;
}

template<typename TInputImage, typename TOutputImage>
void
SymmetricSpectralAnalysisImageFilter<TInputImage, TOutputImage>
::AllocateOutputs()
{
    typename OutputImageType::Pointer outputPtr_val;
    typename OutputImageType::Pointer outputPtr_vec;

    outputPtr_val = dynamic_cast< OutputImageType *>( this->ProcessObject::GetOutput(0) );
    outputPtr_vec = dynamic_cast< OutputImageType *>( this->ProcessObject::GetOutput(1) );

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
    if( m_Order == OrderByMagnitude ) {
        this->m_Calculator.SetOrderEigenMagnitudes( true );
    }
    else if( m_Order == DoNotOrder ) {
        this->m_Calculator.SetOrderEigenValues( false );
    }
}

template<typename TInputImage, typename TOutputImage>
void 
SymmetricSpectralAnalysisImageFilter<TInputImage, TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int )
{
    typename OutputImageType::Pointer output_val = static_cast< OutputImageType * >(this->ProcessObject::GetOutput(0));
    typename OutputImageType::Pointer output_vec = static_cast< OutputImageType * >(this->ProcessObject::GetOutput(1));

    ImageRegionConstIteratorWithIndex< OutputImageType > it(output_val, outputRegionForThread);
    it.GoToBegin();

    while( !it.IsAtEnd() ) {
        InputPixelType dt6 = this->GetInput(0)->GetPixel(it.GetIndex());
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

        OutputPixelType val = output_val->GetPixel(it.GetIndex());
        for (unsigned int i=0; i<3; i++) { val[i] = (OutputValueType) eigenvalues[i]; }
        OutputPixelType vec = output_vec->GetPixel(it.GetIndex());
        for (unsigned int i=0; i<3; i++) {
            for (unsigned int j=0; j<3; j++) { 
                vec[3*i+j] = (OutputValueType) eigenvectors[i][j]; 
            }
        }

        ++it;
    }

}

}

#endif

							
