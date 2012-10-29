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
    this->vector_size = 6;
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
::BeforeThreadedGenerateData()
{
      this->shift_plane = (w_size_plane-1)/2;
      this->shift_depth = (w_size_depth-1)/2;
      this->nb_elments = w_size_depth*w_size_plane*w_size_plane;
}

template<typename TInputImage, typename TOutputImage>
void 
ParameterEstimationImageFilter<TInputImage, TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int )
{
    InputImageType::SizeType size = this->GetInput(0)->GetLargestPossibleRegion().GetSize();

    typename OutputImageType::Pointer output_mean = static_cast< OutputImageType * >(this->ProcessObject::GetOutput(0));
    typename OutputImageType::Pointer output_var = static_cast< OutputImageType * >(this->ProcessObject::GetOutput(1));

    typename OutputImageType::Pointer output = static_cast< OutputImageType * >(this->ProcessObject::GetOutput(0));
    ImageRegionConstIteratorWithIndex< OutputImageType > it(output, outputRegionForThread);
    it.GoToBegin();

    while( !it.IsAtEnd() ) {
        OutputImageType::IndexType idx = it.GetIndex();
        OutputImageType::IndexType idx_;

        // test if the neigborhood is correct
        if ( ((idx[2]-shift_depth)>=0) && ((idx[1]-shift_plane)>=0) && ((idx[0]-shift_plane)>=0) && 
             ((idx[2]+shift_depth)<size[2]) && ((idx[1]+shift_plane)<size[1]) && ((idx[0]+shift_plane)<size[0]) ) {

            // init mean and variance
            vnl_vector_fixed<float, this->vector_size> temp_mean;
            temp_mean.fill(0);
            float temp_var = 0;
	
            // first compute the mean
            for (int kk=idx[2]-shift_depth; kk<=idx[2]+shift_depth; kk++) {
                for (int jj=idx[1]-shift_plane; jj<=idx[1]+shift_plane; jj++) {
                    for (int ii=idx[0]-shift_plane; ii<=idx[0]+shift_plane; ii++) {
                        idx_[0]=ii; idx_[1]=jj; idx_[2]=kk;
                        InputPixelType const tensor = this->GetInput(0)->GetPixel(idx_);
                        for (unsigned int l=0; l<this->vector_size; l++) { temp_mean[l] += (float)tensor[l]/(float)N; }
                    }
                }
            }

            // return the mean
            OutputPixelType mean = output_mean->GetPixel(idx);
            for (unsigned int l=0; l<this->vector_size; l++) { mean[l] = (OutputValueType)temp_mean[l]; }

            // then compute the variance
            for (int kk=idx[2]-shift_depth; kk<=idx[2]+shift_depth; kk++) {
                for (int jj=idx[1]-shift_plane; jj<=idx[1]+shift_plane; jj++) {
                    for (int ii=idx[0]-shift_plane; ii<=idx[0]+shift_plane; ii++) {
                        idx_[0]=ii; idx_[1]=jj; idx_[2]=kk;
                        InputPixelType const tensor = this->GetInput(0)->GetPixel(idx_);
                        temp_var += ((float)tensor[0] - temp_mean[0])*((float)tensor[0] - temp_mean[0])/(float)N 
                                 + ((float)tensor[3] - temp_mean[3])*((float)tensor[3] - temp_mean[3])/(float)N 
                                 + ((float)tensor[5] - temp_mean[5])*((float)tensor[5] - temp_mean[5])/(float)N 
                                 + 2.0*((float)tensor[1] - temp_mean[1])*((float)tensor[1] - temp_mean[1])/(float)N 
                                 + 2.0*((float)tensor[2] - temp_mean[2])*((float)tensor[2] - temp_mean[2])/(float)N
                                 + 2.0*((float)tensor[4] - temp_mean[4])*((float)tensor[4] - temp_mean[4])/(float)N;
                    }
                }
            }

            //return the variance
			OutputPixelType var = output_var->GetPixel(idx);
            var[0] = (OutputValueType)temp_var;

        }
    }
}


}

#endif

							
