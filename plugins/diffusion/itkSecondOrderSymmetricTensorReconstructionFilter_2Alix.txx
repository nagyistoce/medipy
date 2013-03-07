

#ifndef _itkSecondOrderSymmetricTensorReconstructionFilter_2Alix_txx
#define _itkSecondOrderSymmetricTensorReconstructionFilter_2Alix_txx

#include "itkSecondOrderSymmetricTensorReconstructionFilter_2Alix.h"

#include <vector>

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include <vnl/vnl_vector.h>

namespace itk
{

template<typename TInputImage, typename TOutputImage>
void
SecondOrderSymmetricTensorReconstructionFilter_2Alix<TInputImage, TOutputImage>
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
SecondOrderSymmetricTensorReconstructionFilter_2Alix<TInputImage, TOutputImage>
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
SecondOrderSymmetricTensorReconstructionFilter_2Alix<TInputImage, TOutputImage>
::GetNumberOfGradientDirections() const
{
    return this->directions.size();
}

template<typename TInputImage, typename TOutputImage>
void
SecondOrderSymmetricTensorReconstructionFilter_2Alix<TInputImage, TOutputImage>
::SetGradientDirection(unsigned int i, DirectionType bvec)
{
    if (i>=this->directions.size()) {
        this->directions.resize(i);
    }
    this->directions.insert(this->directions.begin()+i,bvec);
}

template<typename TInputImage, typename TOutputImage>
typename SecondOrderSymmetricTensorReconstructionFilter_2Alix<TInputImage, TOutputImage>::DirectionType const &
SecondOrderSymmetricTensorReconstructionFilter_2Alix<TInputImage, TOutputImage>
::GetGradientDirection(unsigned int i) const
{
    return this->directions[i];
}

template<typename TInputImage, typename TOutputImage>
void 
SecondOrderSymmetricTensorReconstructionFilter_2Alix<TInputImage, TOutputImage>
::BeforeThreadedGenerateData()
{
    const unsigned int VectorLength = 6;
    unsigned int nb_dir = this->directions.size();

    this->bmatrix.set_size(nb_dir-1,VectorLength);
 
    for (unsigned int i=1; i<nb_dir; i++) 
    {					
        DirectionType bvec = this->directions[i];
        this->bmatrix(i-1,0) = (float) -this->m_BVal*bvec[0]*bvec[0];        //Dxx
        this->bmatrix(i-1,1) = (float) -this->m_BVal*2.0*bvec[0]*bvec[1];    //Dxy
        this->bmatrix(i-1,2) = (float) -this->m_BVal*2.0*bvec[0]*bvec[2];    //Dxz
        this->bmatrix(i-1,3) = (float) -this->m_BVal*bvec[1]*bvec[1];        //Dyy
        this->bmatrix(i-1,4) = (float) -this->m_BVal*2.0*bvec[1]*bvec[2];    //Dyz
        this->bmatrix(i-1,5) = (float) -this->m_BVal*bvec[2]*bvec[2];        //Dzz
    }
    this->invbmatrix.set_size(this->bmatrix.cols(),this->bmatrix.rows()); 
    BMatrixType b1 = this->bmatrix.transpose();
    BMatrixType b2 = vnl_matrix_inverse<float>(b1*this->bmatrix);
    this->invbmatrix = b2*b1;
}

template<typename TInputImage, typename TOutputImage>
void 
SecondOrderSymmetricTensorReconstructionFilter_2Alix<TInputImage, TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int)
{
    const InputImagePixelType min_signal = 5;
    const unsigned int VectorLength = 6;
    unsigned int const nb_dir = this->directions.size();
    unsigned int nb_iter = this->m_IterationCount;
    
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
        
        for (unsigned int i=1; i<nb_dir; ++i) 
        {
            InputImagePixelType Si = inputIterators[i].Get();
            if (Si<min_signal)
            {
                Si=min_signal;
            }
            if (S0>=Si)
            {
                S(i-1) = (float) log(Si/S0);
            }
        }

        vnl_vector<float> dt6 = this->invbmatrix*S;
        
        // changement WLS:
        //int nb_iter = 5;
        vnl_vector<float> W(nb_dir-1,0.0);
        BMatrixType tmp1;
        tmp1.set_size(VectorLength,VectorLength);
        vnl_vector<float> tmp2(VectorLength,0.0);
        
        for(unsigned int iter=0; iter<nb_iter; iter++)
        {
            W.fill(0.);
            for(unsigned int i=0; i<nb_dir-1; i++)
            {
                W(i) = (float) exp(2.0*inner_product(this->bmatrix.get_row(i), dt6));
            }
            BMatrixType b1 = this->bmatrix.transpose();
            
            tmp1.fill(0.);
            tmp2.fill(0.);
            
            for(unsigned int i=0; i<nb_dir-1; i++)
            {
                tmp1 = tmp1 + W(i) * outer_product(b1.get_column(i), this->bmatrix.get_row(i));
                tmp2 = tmp2 + W(i) * b1.get_column(i) * S(i);
            }
            
            BMatrixType itmp1;
            itmp1.set_size(tmp1.cols(),tmp1.rows()); 
            itmp1 = vnl_matrix_inverse<float>(tmp1);
            dt6 = itmp1*tmp2;
        }
        // fin changement WLS
        
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

							
