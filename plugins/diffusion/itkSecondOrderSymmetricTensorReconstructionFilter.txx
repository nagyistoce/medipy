/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _itkSecondOrderSymmetricTensorReconstructionFilter_txx
#define _itkSecondOrderSymmetricTensorReconstructionFilter_txx

#include "itkSecondOrderSymmetricTensorReconstructionFilter.h"

#include "itkVectorImage.h"
#include "itkImageRegionIteratorWithIndex.h"

namespace itk
{

template<typename TInputImage, typename TOutputImage>
void
SecondOrderSymmetricTensorReconstructionFilter<TInputImage, TOutputImage>
::AllocateOutputs()
{
    typename OutputImageType::Pointer outputPtr;

    outputPtr = dynamic_cast< OutputImageType *>( this->ProcessObject::GetOutput() );

    if ( outputPtr ) {
        outputPtr->SetBufferedRegion( outputPtr->GetRequestedRegion() );
        outputPtr->SetVectorLength(6);
        outputPtr->Allocate();
    }
}

template<typename TInputImage, typename TOutputImage>
void
SecondOrderSymmetricTensorReconstructionFilter<TInputImage, TOutputImage>
::SetGradientDirection(unsigned int i, DirectionType bvec)
{
    if (i>=this->directions.size()) {
        this->directions.resize(i);
    }
    this->directions.insert(this->directions.begin()+i,bvec);
}

template<typename TInputImage, typename TOutputImage>
typename SecondOrderSymmetricTensorReconstructionFilter<TInputImage, TOutputImage>::DirectionType
SecondOrderSymmetricTensorReconstructionFilter<TInputImage, TOutputImage>
::GetGradientDirection(unsigned int i)
{
    if (i<this->directions.size()) {
        return this->directions[i];
    }
}

template<typename TInputImage, typename TOutputImage>
void
SecondOrderSymmetricTensorReconstructionFilter<TInputImage, TOutputImage>
::GenerateData()
{
    const PixelType min_signal = 5;
    const unsigned int VectorLength = 6;

    unsigned int nb_dir = this->directions.size();

    this->bmatrix.set_size(nb_dir-1,6);
    for (unsigned int i=1; i<nb_dir; i++) {					
        DirectionType bvec = this->directions[i];
        this->bmatrix(i-1,0) = (float) this->m_BVal*bvec[0]*bvec[0];        //Dxx
        this->bmatrix(i-1,1) = (float) this->m_BVal*2.0*bvec[0]*bvec[1];    //Dxy
        this->bmatrix(i-1,2) = (float) this->m_BVal*2.0*bvec[0]*bvec[2];    //Dxz
        this->bmatrix(i-1,3) = (float) this->m_BVal*bvec[1]*bvec[1];        //Dyy
        this->bmatrix(i-1,4) = (float) this->m_BVal*2.0*bvec[1]*bvec[2];    //Dyz
        this->bmatrix(i-1,5) = (float) this->m_BVal*bvec[2]*bvec[2];        //Dzz
    }
    this->invbmatrix.set_size(this->bmatrix.cols(),this->bmatrix.rows()); 
    BMatrixType b1 = this->bmatrix.transpose();
    BMatrixType b2 = vnl_matrix_inverse<float>(b1*this->bmatrix);
    this->invbmatrix = b2*b1;

    InputImageType const * im_ref = this->GetInput(0);
    typename InputImageType::SizeType size = im_ref->GetLargestPossibleRegion().GetSize();
    typename InputImageType::IndexType start = im_ref->GetLargestPossibleRegion().GetIndex();

    //VariableVectorType f(VectorLength);
    //for( unsigned int i=0; i<VectorLength; i++ ) { f[i] = 0; }

    OutputImageType * output = this->GetOutput();
    output->Print(std::cout);
    /*typename OutputImageType::RegionType region;
    region.SetSize( size );
    region.SetIndex( start );
    output->SetVectorLength(VectorLength);
    output->SetSpacing(im_ref->GetSpacing());
    output->SetOrigin(im_ref->GetOrigin());
    output->SetDirection(im_ref->GetDirection());
    output->SetRegions( region );
    output->Allocate();*/
    //output->FillBuffer(f);  

    ImageRegionConstIteratorWithIndex<InputImageType> it(im_ref, im_ref->GetRequestedRegion());
    //ImageRegionConstIteratorWithIndex<InputImageType> it(im_ref, regionThread);
    it.GoToBegin();
    while( !it.IsAtEnd() ) {

        vnl_matrix<float> S(nb_dir-1, 1, 0.0); 
        PixelType S0 = it.Get();
        if (S0<min_signal) { S0=min_signal; }
        for (unsigned int i=1; i<nb_dir; ++i) {
            float Si = this->GetInput(i)->GetPixel(it.GetIndex());
            if (Si<min_signal) { Si=min_signal; }
            if (S0>=Si) { S(i-1,0) = log(S0/Si); }
        }
        vnl_matrix<float> dt6(VectorLength, 1);
        dt6 = this->invbmatrix*S;

        OutputPixelType vec = output->GetPixel(it.GetIndex());
        for( unsigned int i=0; i<VectorLength; i++ ) { vec[i] = (typename OutputPixelType::ValueType) dt6(i,0); }

        ++it;
    }
}


}

#endif

							
