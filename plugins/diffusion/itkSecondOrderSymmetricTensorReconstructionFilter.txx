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

#include <vector>

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include <vnl/vnl_vector_fixed.h>

namespace itk
{

template<typename TInputImage, typename TOutputImage>
void
SecondOrderSymmetricTensorReconstructionFilter<TInputImage, TOutputImage>
::AllocateOutputs()
{
    typename OutputImageType::Pointer outputPtr;

    outputPtr = dynamic_cast< OutputImageType *>( this->ProcessObject::GetOutput(0) );

    if ( outputPtr ) {
        outputPtr->SetBufferedRegion( outputPtr->GetRequestedRegion() );
        outputPtr->SetVectorLength(6);
        outputPtr->Allocate();
    }
}

template<typename TInputImage, typename TOutputImage>
void
SecondOrderSymmetricTensorReconstructionFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
    std::locale C("C");
    std::locale originalLocale = os.getloc();
    os.imbue(C);

    Superclass::PrintSelf(os,indent);
 
    os << indent << "BValue: " << m_BVal << "\n";
    os << indent << "NumberOfGradientDirections: " << this->directions.size() << "\n";
    for(unsigned int i=0; i<this->directions.size(); ++i)
    {
        os << indent.GetNextIndent()
           << "Direction " << (i+1) << ": " << this->directions[i] << "\n";
    }

    os.imbue( originalLocale );
}

template<typename TInputImage, typename TOutputImage>
unsigned int
SecondOrderSymmetricTensorReconstructionFilter<TInputImage, TOutputImage>
::GetNumberOfGradientDirections() const
{
    return this->directions.size();
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
typename SecondOrderSymmetricTensorReconstructionFilter<TInputImage, TOutputImage>::DirectionType const &
SecondOrderSymmetricTensorReconstructionFilter<TInputImage, TOutputImage>
::GetGradientDirection(unsigned int i) const
{
    return this->directions[i];
}

template<typename TInputImage, typename TOutputImage>
void 
SecondOrderSymmetricTensorReconstructionFilter<TInputImage, TOutputImage>
::BeforeThreadedGenerateData()
{
    const unsigned int VectorLength = 6;
    unsigned int nb_dir = this->directions.size();

    this->bmatrix.set_size(nb_dir-1,VectorLength);
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
}

template<typename TInputImage, typename TOutputImage>
void 
SecondOrderSymmetricTensorReconstructionFilter<TInputImage, TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int )
{
    const InputImagePixelType min_signal = 5;
    const unsigned int VectorLength = 6;
    unsigned int const nb_dir = this->directions.size();

    typename OutputImageType::Pointer output = static_cast<OutputImageType *>(
        this->ProcessObject::GetOutput(0));

    // Create an iterator for each input
    typedef ImageRegionConstIterator<InputImageType> InputIterator;
    std::vector<InputIterator> inputIterators;
    for (unsigned int i=0; i<this->GetNumberOfInputs(); ++i)
    {
        InputIterator iterator(this->GetInput(i), outputRegionForThread);
        inputIterators.push_back(iterator);
    }

    vnl_vector<float> S(nb_dir-1, 0.0);

    typedef ImageRegionIterator<OutputImageType> OutputIterator;
    for(OutputIterator outputIt(output, outputRegionForThread);
        !outputIt.IsAtEnd(); ++outputIt)
    {
        // Set the signal vector to 0, to avoid using previous values if S0<Si
        S.fill(0.);

        InputImagePixelType S0 = inputIterators[0].Get();
        if (S0<min_signal)
        {
            S0=min_signal;
        }

        for (unsigned int i=1; i<nb_dir; ++i) {
            InputImagePixelType Si = inputIterators[i].Get();
            if (Si<min_signal)
            {
                Si=min_signal;
            }
            if (S0>=Si)
            {
                S(i-1) = log(S0/Si);
            }
        }

        vnl_vector_fixed<float, VectorLength> dt6;
        dt6 = this->invbmatrix*S;

        OutputPixelType vec = outputIt.Get();
        std::copy(dt6.begin(), dt6.end(), &vec[0]);

        for(typename std::vector<InputIterator>::iterator inputIteratorsIt=inputIterators.begin();
            inputIteratorsIt!=inputIterators.end(); ++inputIteratorsIt)
        {
            ++(*inputIteratorsIt);
        }
    }

}

}

#endif

							
