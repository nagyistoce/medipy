/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _itkLeastSquaresImageFilter_txx
#define _itkLeastSquaresImageFilter_txx

#include "itkLeastSquaresImageFilter.h"

#include <vector>

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <vnl/vnl_vector_fixed.h>

namespace itk
{
template<typename TInputImage, typename TTensorsImage, typename TBaselineImage, typename TMaskImage>
void
LeastSquaresImageFilter<TInputImage, TTensorsImage, TBaselineImage, TMaskImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
    std::locale C("C");
    std::locale originalLocale = os.getloc();
    os.imbue(C);

    Superclass::PrintSelf(os,indent);
    
    os << indent << "NumberOfGradientDirections: " << this->directions.size() << "\n";
    for(unsigned int i=0; i<this->directions.size(); ++i)
    {
        os << indent.GetNextIndent()
           << "Bvalue " << (i+1) << ": " << this->metadata_diffusion[i].first << "\n"
           << "Direction " << (i+1) << ": " << this->metadata_diffusion[i].second << "\n";
    }
    
    os.imbue( originalLocale );
}

template<typename TInputImage, typename TTensorsImage, typename TBaselineImage, typename TMaskImage>
unsigned int
LeastSquaresImageFilter<TInputImage, TTensorsImage, TBaselineImage, TMaskImage>
::GetNumberOfGradientDirections() const
{
    return this->directions.size();
}

template<typename TInputImage, typename TTensorsImage, typename TBaselineImage, typename TMaskImage>
void
LeastSquaresImageFilter<TInputImage, TTensorsImage, TBaselineImage, TMaskImage>
::SetGradientDirection(unsigned int i, DirectionType bvec)
{
    if (i>=this->directions.size())
    {
        this->directions.resize(i);
    }
    this->directions.insert(this->directions.begin()+i,bvec);
}

template<typename TInputImage, typename TTensorsImage, typename TBaselineImage, typename TMaskImage>
typename LeastSquaresImageFilter<TInputImage, TTensorsImage, TBaselineImage, TMaskImage>::DirectionType const &
LeastSquaresImageFilter<TInputImage, TTensorsImage, TBaselineImage, TMaskImage>
::GetGradientDirection(unsigned int i) const
{
    return this->directions[i];
}

template<typename TInputImage, typename TTensorsImage, typename TBaselineImage, typename TMaskImage>
void
LeastSquaresImageFilter<TInputImage, TTensorsImage, TBaselineImage, TMaskImage>
::SetBvalueAndGradientDirection(unsigned int i, BValueType BVal, DirectionType vec)
{
    if (i>=this->metadata_diffusion.size())
    {
        this->metadata_diffusion.resize(i);
    }
    MetaDiffusionType p = std::make_pair(BVal, vec);
    this->metadata_diffusion.push_back(p);
}

template<typename TInputImage, typename TTensorsImage, typename TBaselineImage, typename TMaskImage>
LeastSquaresImageFilter<TInputImage, TTensorsImage, TBaselineImage, TMaskImage>
::LeastSquaresImageFilter() 
{
    //this->SetMaskImage(NULL);
}

template<typename TInputImage, typename TTensorsImage, typename TBaselineImage, typename TMaskImage>
void 
LeastSquaresImageFilter<TInputImage, TTensorsImage, TBaselineImage, TMaskImage>
::BeforeThreadedGenerateData()
{
    const unsigned int VectorLength = 6;
    unsigned int nb_dir = this->directions.size();
    
    this->bmatrix.set_size(nb_dir,VectorLength+1);
    for (unsigned int i=0; i<nb_dir; i++)
    {
        MetaDiffusionType p = this->metadata_diffusion[i];
        BValueType BVal = p.first;
        DirectionType bvec = p.second;
        this->bmatrix(i,0) = (float) 1.0;                              //log(S0)
        this->bmatrix(i,1) = (float) -1.0*BVal*bvec[0]*bvec[0];        //Dxx
        this->bmatrix(i,2) = (float) -1.0*BVal*2.0*bvec[0]*bvec[1];    //Dxy
        this->bmatrix(i,3) = (float) -1.0*BVal*2.0*bvec[0]*bvec[2];    //Dxz
        this->bmatrix(i,4) = (float) -1.0*BVal*bvec[1]*bvec[1];        //Dyy
        this->bmatrix(i,5) = (float) -1.0*BVal*2.0*bvec[1]*bvec[2];    //Dyz
        this->bmatrix(i,6) = (float) -1.0*BVal*bvec[2]*bvec[2];        //Dzz
    }
    this->invbmatrix.set_size(this->bmatrix.cols(),this->bmatrix.rows()); 
    BMatrixType b1 = this->bmatrix.transpose();
    BMatrixType b2 = vnl_matrix_inverse<float>(b1*this->bmatrix);
    this->invbmatrix = b2*b1;
    
    TensorsImagePointer tensors =
        dynamic_cast<TTensorsImage*>(this->ProcessObject::GetOutput(0));
    TensorsPixelType zero(6);
    zero.Fill(0);
    tensors->FillBuffer(zero);
    
    BaselineImagePointer baseline = 
        dynamic_cast<TBaselineImage*>(this->ProcessObject::GetOutput(1));
    baseline->FillBuffer(0);
}

template<typename TInputImage, typename TTensorsImage, typename TBaselineImage, typename TMaskImage>
void 
LeastSquaresImageFilter<TInputImage, TTensorsImage, TBaselineImage, TMaskImage>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int )
{
    const InputImagePixelType min_signal = 5;
    const unsigned int VectorLength = 6;
    unsigned int const nb_dir = this->directions.size();

    TensorsImagePointer tensors =
        dynamic_cast<TensorsImageType*>(this->ProcessObject::GetOutput(0));
    
    BaselineImagePointer baseline = 
        dynamic_cast<TBaselineImage*>(this->ProcessObject::GetOutput(1));

    // Create an iterator for each input
    typedef ImageRegionConstIterator<InputImageType> InputIterator;
    std::vector<InputIterator> inputIterators;
    for (unsigned int i=0; i<this->GetNumberOfInputs(); ++i)
    {
        InputIterator iterator(this->GetInput(i), outputRegionForThread);
        inputIterators.push_back(iterator);
    }

    vnl_vector<float> S(nb_dir, 0.0);

    typedef ImageRegionIteratorWithIndex<TensorsImageType> OutputIterator;
    for(OutputIterator outputIt(tensors, outputRegionForThread);
        !outputIt.IsAtEnd(); ++outputIt)
    {
        typename TensorsImageType::IndexType const & idx = outputIt.GetIndex();
        
        // Skip pixels that are not in mask
        bool process=true;
        if(this->m_MaskImage.GetPointer() != 0)
        {
            typename TTensorsImage::PointType point; 
            tensors->TransformIndexToPhysicalPoint(idx, point);
            
            typename TMaskImage::IndexType mask_index;
            this->m_MaskImage->TransformPhysicalPointToIndex(point, mask_index);
            
            if(!this->m_MaskImage->GetLargestPossibleRegion().IsInside(mask_index) || 
               this->m_MaskImage->GetPixel(mask_index) == 0)
            {
                process=false;
            }
        }
        
        if(process)
        {
            // Set the signal vector to 0, to avoid using previous values if S0<Si
            S.fill(0.);
            
            for (unsigned int i=0; i<nb_dir; ++i)
            {
                InputImagePixelType Si = inputIterators[i].Get();
                if (Si<min_signal)
                {
                    Si=min_signal;
                }
                S(i) = (float) log(Si);
            }
            
            vnl_vector_fixed<float, VectorLength+1> X;
            X = this->invbmatrix*S;
            
            
            TensorsPixelType vec = tensors->GetPixel(idx);
            vec.Fill(0.);
            std::copy(X.begin()+1, X.end(), &vec[0]);
            
            baseline->SetPixel(idx, exp(X(0)));
        }

        for(typename std::vector<InputIterator>::iterator inputIteratorsIt=inputIterators.begin();
            inputIteratorsIt!=inputIterators.end(); ++inputIteratorsIt)
        {
            ++(*inputIteratorsIt);
        }
    }

}

}

#endif
