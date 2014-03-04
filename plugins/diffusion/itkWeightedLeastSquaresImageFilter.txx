

#ifndef _itkWeightedLeastSquaresImageFilter_txx
#define _itkWeightedLeastSquaresImageFilter_txx

#include "itkWeightedLeastSquaresImageFilter.h"

#include <vector>

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matlab_print.h>

namespace itk
{

template<typename TInputImage, typename TTensorsImage, typename TBaselineImage, typename TMaskImage>
unsigned int
WeightedLeastSquaresImageFilter<TInputImage, TTensorsImage, TBaselineImage, TMaskImage>
::GetNumberOfGradientDirections() const
{
    return this->metadata_diffusion.size();
}

template<typename TInputImage, typename TTensorsImage, typename TBaselineImage, typename TMaskImage>
typename WeightedLeastSquaresImageFilter<TInputImage, TTensorsImage, TBaselineImage, TMaskImage>::DirectionType const &
WeightedLeastSquaresImageFilter<TInputImage, TTensorsImage, TBaselineImage, TMaskImage>
::GetGradientDirection(unsigned int i) const
{
    return this->metadata_diffusion[i].second;
}

template<typename TInputImage, typename TTensorsImage, typename TBaselineImage, typename TMaskImage>
typename WeightedLeastSquaresImageFilter<TInputImage, TTensorsImage, TBaselineImage, TMaskImage>::BValueType
WeightedLeastSquaresImageFilter<TInputImage, TTensorsImage, TBaselineImage, TMaskImage>
::GetBvalue(unsigned int i) const
{
    return this->metadata_diffusion[i].first;
}

template<typename TInputImage, typename TTensorsImage, typename TBaselineImage, typename TMaskImage>
void
WeightedLeastSquaresImageFilter<TInputImage, TTensorsImage, TBaselineImage, TMaskImage>
::SetBvalueAndGradientDirection(unsigned int i, BValueType BVal, DirectionType vec)
{
    if (i>=this->metadata_diffusion.size())
    {
        this->metadata_diffusion.resize(i+1);
    }
    this->metadata_diffusion[i] = std::make_pair(BVal, vec);
}

template<typename TInputImage, typename TTensorsImage, typename TBaselineImage, typename TMaskImage>
void
WeightedLeastSquaresImageFilter<TInputImage, TTensorsImage, TBaselineImage, TMaskImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
    std::locale C("C");
    std::locale originalLocale = os.getloc();
    os.imbue(C);

    Superclass::PrintSelf(os,indent);
 
    os << indent << "NumberOfGradientDirections: " 
       << this->GetNumberOfGradientDirections() << "\n";
    for(unsigned int i=0; i<this->GetNumberOfGradientDirections(); ++i)
    {
        os << indent.GetNextIndent()
           << "Bvalue " << (i+1) << ": " << this->GetBvalue(i) << "\n"
           << "Direction " << (i+1) << ": " << this->GetGradientDirection(i) << "\n";
    }
    
    os.imbue( originalLocale );
}

template<typename TInputImage, typename TTensorsImage, typename TBaselineImage, typename TMaskImage>
void 
WeightedLeastSquaresImageFilter<TInputImage, TTensorsImage, TBaselineImage, TMaskImage>
::BeforeThreadedGenerateData()
{
    // Fill the b-matrix
    this->bmatrix = BMatrixType(
        this->GetNumberOfGradientDirections(), 
        this->GetTensorsImage()->GetVectorLength()+1, 0.);
    
    for (unsigned int i=0; i<this->GetNumberOfGradientDirections(); i++) 
    {
        BValueType const BVal = this->metadata_diffusion[i].first;
        DirectionType const bvec = this->metadata_diffusion[i].second;
        
        this->bmatrix(i,0) = (float) 1.0;                              //log(S0)
        this->bmatrix(i,1) = (float) -1.0*BVal*bvec[0]*bvec[0];        //Dxx
        this->bmatrix(i,2) = (float) -1.0*BVal*2.0*bvec[0]*bvec[1];    //Dxy
        this->bmatrix(i,3) = (float) -1.0*BVal*2.0*bvec[0]*bvec[2];    //Dxz
        this->bmatrix(i,4) = (float) -1.0*BVal*bvec[1]*bvec[1];        //Dyy
        this->bmatrix(i,5) = (float) -1.0*BVal*2.0*bvec[1]*bvec[2];    //Dyz
        this->bmatrix(i,6) = (float) -1.0*BVal*bvec[2]*bvec[2];        //Dzz
    }
    
    // Compute the pseudo-inverse of the b-matrix
    BMatrixType const b1 = this->bmatrix.transpose();
    BMatrixType const b2 = vnl_matrix_inverse<float>(b1*this->bmatrix);
    this->invbmatrix = b2*b1;
    
    // Zero-initialize the outputs of the filter
    TensorsImagePointer tensors =
        dynamic_cast<TTensorsImage*>(this->ProcessObject::GetOutput(0));
    typename TensorsImageType::PixelType zero(6);
    zero.Fill(0);
    tensors->FillBuffer(zero);
    
    BaselineImagePointer baseline = 
        dynamic_cast<TBaselineImage*>(this->ProcessObject::GetOutput(1));
    baseline->FillBuffer(0);
}

template<typename TInputImage, typename TTensorsImage, typename TBaselineImage, typename TMaskImage>
void 
WeightedLeastSquaresImageFilter<TInputImage, TTensorsImage, TBaselineImage, TMaskImage>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int)
{
    // Create an iterator for each input
    typedef ImageRegionConstIterator<InputImageType> InputIterator;
    std::vector<InputIterator> inputIterators;
    for (unsigned int i=0; i<this->GetNumberOfGradientDirections(); ++i)
    {
        InputIterator iterator(this->GetInput(i), outputRegionForThread);
        inputIterators.push_back(iterator);
    }

    typedef ImageRegionIteratorWithIndex<TensorsImageType> OutputIterator;
    for(OutputIterator outputIt(this->GetTensorsImage(), outputRegionForThread);
        !outputIt.IsAtEnd(); ++outputIt)
    {
        typename TensorsImageType::IndexType const & idx = outputIt.GetIndex();
        
        // Skip pixels that are not in mask
        bool process=true;
        if(this->GetMaskImage() != 0)
        {
            typename TensorsImageType::PointType point; 
            this->GetTensorsImage()->TransformIndexToPhysicalPoint(idx, point);
            
            typename MaskImageType::IndexType mask_index;
            this->GetMaskImage()->TransformPhysicalPointToIndex(point, mask_index);
            
            if(!this->GetMaskImage()->GetRequestedRegion().IsInside(mask_index) || 
               this->GetMaskImage()->GetPixel(mask_index) == 0)
            {
                process=false;
            }
        }

        if(process)
        {
            vnl_vector<float> S(this->GetNumberOfGradientDirections(), 0.0);
            
            for (unsigned int i=0; i<this->GetNumberOfGradientDirections(); ++i) 
            {
                typename InputImageType::PixelType const min_signal = 5;
                typename InputImageType::PixelType const Si = std::max(
                    min_signal, inputIterators[i].Get());
                S(i) = log(Si);
            }

            vnl_vector<float> X = this->invbmatrix*S;
            
            for(unsigned int iter=0; iter<this->m_IterationCount; iter++)
            {
                // Use double instead of float to avoid overflow in the 
                // exponential
                vnl_vector<double> W(this->GetNumberOfGradientDirections(), 0.0);
                for(unsigned int i=0; i<this->GetNumberOfGradientDirections(); i++)
                {
                    W(i) = exp(2.0*inner_product(this->bmatrix.get_row(i), X));
                }
                
                unsigned int const vector_length = 
                    this->GetTensorsImage()->GetVectorLength();
                vnl_matrix<double> tmp1(vector_length+1, vector_length+1, 0.);
                vnl_vector<double> tmp2(vector_length+1, 0.);
                
                for(unsigned int i=0; i<this->GetNumberOfGradientDirections(); i++)
                {
                    // Copy the row from float to double
                    vnl_vector<float> const row_float = this->bmatrix.get_row(i);
                    vnl_vector<double> row(row_float.size());
                    
                    std::copy(row_float.begin(), row_float.end(), row.begin());
                    
                    tmp1 += W(i) * outer_product(row, row);
                    tmp2 += W(i) * row * S(i);
                }
                
                // Copy X from double to float
                vnl_vector<double> const X_double = 
                    vnl_matrix_inverse<double>(tmp1)*tmp2;
                std::copy(X_double.begin(), X_double.end(), X.begin());
            }
            
            typename TensorsImageType::PixelType vec = 
                this->GetTensorsImage()->GetPixel(idx);
            std::copy(X.begin()+1, X.end(), &vec[0]);
            
            this->GetBaselineImage()->SetPixel(idx, exp(X(0)));
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
