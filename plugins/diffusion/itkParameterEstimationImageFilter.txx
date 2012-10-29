/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _ParameterEstimationImageFilter_txx
#define _ParameterEstimationImageFilter_txx

#include "itkParameterEstimationImageFilter.h"

#include "itkVectorImage.h"
#include "itkImageRegionIteratorWithIndex.h"

#include <vnl/vnl_vector_fixed.h>
#include <math.h>

namespace itk
{

template<typename TInputImage, typename TOutputImage>
ParameterEstimationImageFilter<TInputImage, TOutputImage>
::ParameterEstimationImageFilter()
{
    this->m_SizePlane = 3;
    this->m_SizeDepth = 3;

    this->SetNumberOfRequiredOutputs(2);
    this->SetNthOutput(0,this->MakeOutput(0));
    this->SetNthOutput(1,this->MakeOutput(1));

    this->masked = false;
}

template<typename TInputImage, typename TOutputImage>
void
ParameterEstimationImageFilter<TInputImage, TOutputImage>
::AllocateOutputs()
{
    typename OutputImageType::Pointer outputPtr_mean;
    typename OutputImageType::Pointer outputPtr_var;

    outputPtr_mean = dynamic_cast< OutputImageType *>( this->ProcessObject::GetOutput(0) );
    outputPtr_var = dynamic_cast< OutputImageType *>( this->ProcessObject::GetOutput(1) );

    if ( outputPtr_mean ) {
        outputPtr_mean->SetBufferedRegion( outputPtr_mean->GetRequestedRegion() );
        outputPtr_mean->SetVectorLength(6);
        outputPtr_mean->Allocate();
    }
    if ( outputPtr_var ) {
        outputPtr_var->SetBufferedRegion( outputPtr_var->GetRequestedRegion() );
        outputPtr_var->SetVectorLength(1);
        outputPtr_var->Allocate();
    }
}

template<typename TInputImage, typename TOutputImage>
void
ParameterEstimationImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
    Superclass::PrintSelf(os,indent);
}

template<typename TInputImage, typename TOutputImage>
void
ParameterEstimationImageFilter<TInputImage, TOutputImage>
::SetMask(MaskType *m)
{
    this->mask = m;
    this->masked = true;
}

template<typename TInputImage, typename TOutputImage>
typename ParameterEstimationImageFilter<TInputImage, TOutputImage>::MaskType*
ParameterEstimationImageFilter<TInputImage, TOutputImage>
::GetMask() const
{
    return this->mask;
}

template<typename TInputImage, typename TOutputImage>
void 
ParameterEstimationImageFilter<TInputImage, TOutputImage>
::BeforeThreadedGenerateData()
{
    this->size = this->GetInput(0)->GetLargestPossibleRegion().GetSize();
    this->shift_plane = (this->m_SizePlane-1)/2;
    this->shift_depth = (this->m_SizeDepth-1)/2;
    this->m_NumberOfElements = this->m_SizeDepth*this->m_SizePlane*this->m_SizePlane;
    if (!this->masked) {
        this->mask = MaskType::New();
        this->mask->SetRegions( this->GetInput(0)->GetLargestPossibleRegion() );
        this->mask->Allocate();
        this->mask->FillBuffer(1);
    }
}

template<typename TInputImage, typename TOutputImage>
void 
ParameterEstimationImageFilter<TInputImage, TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int )
{
    const unsigned int VectorLength = 6;

    typename OutputImageType::Pointer output_mean = static_cast< OutputImageType * >(this->ProcessObject::GetOutput(0));
    typename OutputImageType::Pointer output_var = static_cast< OutputImageType * >(this->ProcessObject::GetOutput(1));

    ImageRegionConstIteratorWithIndex< OutputImageType > it(output_mean, outputRegionForThread);
    it.GoToBegin();

    while( !it.IsAtEnd() ) {
        typename OutputImageType::IndexType idx = it.GetIndex();

        if (this->mask->GetPixel(idx)==1) {
            typename OutputImageType::IndexType idx_;

            // test if the neigborhood is correct
            if ( ((idx[2]-this->shift_depth)>=0) && ((idx[1]-this->shift_plane)>=0) && ((idx[0]-this->shift_plane)>=0) && 
                 ((idx[2]+this->shift_depth)<size[2]) && ((idx[1]+this->shift_plane)<size[1]) && ((idx[0]+this->shift_plane)<size[0]) ) {

                // init mean and variance
                vnl_vector_fixed<float, VectorLength> temp_mean;
                temp_mean.fill(0);
                float temp_var = 0;
	
                // first compute the mean
                for (int kk=idx[2]-this->shift_depth; kk<=idx[2]+this->shift_depth; kk++) {
                    for (int jj=idx[1]-this->shift_plane; jj<=idx[1]+this->shift_plane; jj++) {
                        for (int ii=idx[0]-this->shift_plane; ii<=idx[0]+this->shift_plane; ii++) {
                            idx_[0]=ii; idx_[1]=jj; idx_[2]=kk;
                            InputPixelType const tensor = this->GetInput(0)->GetPixel(idx_);
                            for (unsigned int l=0; l<VectorLength; l++) { temp_mean[l] += (float)tensor[l]/(float)this->m_NumberOfElements; }
                        }
                    }
                }

                // return the mean
                OutputPixelType mean = output_mean->GetPixel(idx);
                for (unsigned int l=0; l<VectorLength; l++) { mean[l] = (OutputValueType)temp_mean[l]; }

                // then compute the variance
                for (int kk=idx[2]-this->shift_depth; kk<=idx[2]+this->shift_depth; kk++) {
                    for (int jj=idx[1]-this->shift_plane; jj<=idx[1]+this->shift_plane; jj++) {
                        for (int ii=idx[0]-this->shift_plane; ii<=idx[0]+this->shift_plane; ii++) {
                            idx_[0]=ii; idx_[1]=jj; idx_[2]=kk;
                            InputPixelType const tensor = this->GetInput(0)->GetPixel(idx_);
                            temp_var += ((float)tensor[0] - temp_mean[0])*((float)tensor[0] - temp_mean[0])/(float)this->m_NumberOfElements
                                     + ((float)tensor[3] - temp_mean[3])*((float)tensor[3] - temp_mean[3])/(float)this->m_NumberOfElements
                                     + ((float)tensor[5] - temp_mean[5])*((float)tensor[5] - temp_mean[5])/(float)this->m_NumberOfElements 
                                     + 2.0*((float)tensor[1] - temp_mean[1])*((float)tensor[1] - temp_mean[1])/(float)this->m_NumberOfElements
                                     + 2.0*((float)tensor[2] - temp_mean[2])*((float)tensor[2] - temp_mean[2])/(float)this->m_NumberOfElements
                                     + 2.0*((float)tensor[4] - temp_mean[4])*((float)tensor[4] - temp_mean[4])/(float)this->m_NumberOfElements;
                        }
                    }
                }

                // return the variance
			    OutputPixelType var = output_var->GetPixel(idx);
                var[0] = (OutputValueType)temp_var;

            }
        }
        ++it;
    }
}


}

#endif

							
