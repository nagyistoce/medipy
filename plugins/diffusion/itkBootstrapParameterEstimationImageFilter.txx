/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _BootstrapParameterEstimationImageFilter_txx
#define _BootstrapParameterEstimationImageFilter_txx

#include "itkBootstrapParameterEstimationImageFilter.h"

#include "itkVectorImage.h"
#include "itkImageRegionIteratorWithIndex.h"

#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <math.h>

namespace itk
{

template<typename TInputImage, typename TOutputImage>
BootstrapParameterEstimationImageFilter<TInputImage, TOutputImage>
::BootstrapParameterEstimationImageFilter()
: m_SizePlane(3), m_SizeDepth(3), m_NumberOfBootstrap(200), m_Mask(NULL), m_UseSpatialBootstrap(true)
{
    this->SetNumberOfRequiredOutputs(2);
    this->SetNthOutput(0,this->MakeOutput(0));
    this->SetNthOutput(1,this->MakeOutput(1));
}

//----------------------------------------------------------------------------------------

template<typename TInputImage, typename TOutputImage>
void
BootstrapParameterEstimationImageFilter<TInputImage, TOutputImage>
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

//----------------------------------------------------------------------------------------

template<typename TInputImage, typename TOutputImage>
void
BootstrapParameterEstimationImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
    Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------------------

template<typename TInputImage, typename TOutputImage>
void
BootstrapParameterEstimationImageFilter<TInputImage, TOutputImage>
::SetGradientDirection(unsigned int i, DirectionType bvec)
{
    if (i>=this->directions.size()) {
        this->directions.resize(i);
    }
    this->directions.insert(this->directions.begin()+i,bvec);
}

//----------------------------------------------------------------------------------------

template<typename TInputImage, typename TOutputImage>
typename BootstrapParameterEstimationImageFilter<TInputImage, TOutputImage>::DirectionType
BootstrapParameterEstimationImageFilter<TInputImage, TOutputImage>
::GetGradientDirection(unsigned int i)
{
    if (i<this->directions.size()) {
        return this->directions[i];
    }
    else {
        throw "Wrong index to access gradients!";
    }
}

//----------------------------------------------------------------------------------------

template<typename TInputImage, typename TOutputImage>
void 
BootstrapParameterEstimationImageFilter<TInputImage, TOutputImage>
::BeforeThreadedGenerateData()
{
    this->size = this->GetInput(0)->GetLargestPossibleRegion().GetSize();
    this->shift_plane = (this->m_SizePlane-1)/2;
    this->shift_depth = (this->m_SizeDepth-1)/2;
    this->m_NumberOfElements = this->m_SizeDepth*this->m_SizePlane*this->m_SizePlane;
    if (this->m_Mask.IsNull()) {
        this->m_Mask = MaskType::New();
        this->m_Mask->SetRegions( this->GetInput(0)->GetLargestPossibleRegion() );
        this->m_Mask->Allocate();
        this->m_Mask->FillBuffer(1);
    }

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

//----------------------------------------------------------------------------------------

template<typename TInputImage, typename TOutputImage>
std::vector< vnl_matrix<float> > 
BootstrapParameterEstimationImageFilter<TInputImage, TOutputImage>
::SpatialBootstrapInit(typename OutputImageType::IndexType idx)
{
    unsigned int nb_dir = this->directions.size();
    std::vector< vnl_matrix<float> > signals;
    vnl_matrix<float> signal(nb_dir, 1); 

    if (this->m_Mask->GetPixel(idx)==1) {
        typename OutputImageType::IndexType idx_;

        // test if the neigborhood is correct
        if ( (((unsigned int)idx[2]-this->shift_depth)>=0) && (((unsigned int)idx[1]-this->shift_plane)>=0) && (((unsigned int)idx[0]-this->shift_plane)>=0) && 
             (((unsigned int)idx[2]+this->shift_depth)<size[2]) && (((unsigned int)idx[1]+this->shift_plane)<size[1]) && (((unsigned int)idx[0]+this->shift_plane)<size[0]) ) {

            // first get the neighborhood signals
            for (int kk=idx[2]-this->shift_depth; kk<=idx[2]+this->shift_depth; kk++) {
                for (int jj=idx[1]-this->shift_plane; jj<=idx[1]+this->shift_plane; jj++) {
                    for (int ii=idx[0]-this->shift_plane; ii<=idx[0]+this->shift_plane; ii++) {
                        idx_[0]=ii; idx_[1]=jj; idx_[2]=kk;
                        for (unsigned int l=0; l<nb_dir; ++l) { signal(l,0) = (float)this->GetInput(l)->GetPixel(idx_); }
                        signals.push_back( signal );
                    }
                }
            }

        }
    }
    return signals;
}

template<typename TInputImage, typename TOutputImage>
std::vector<unsigned int> 
BootstrapParameterEstimationImageFilter<TInputImage, TOutputImage>
::RandomSpatial(unsigned int N)
{
    static int first = 0;
    std::vector<unsigned int> index(N,0);

    for (unsigned int i=0; i<N; i++) {
        if (first==0){
            srand (time (NULL));
            first = 1;
        }
        index[i] = (unsigned int)( (rand()/(double)RAND_MAX)*(this->m_NumberOfElements-1.0));
    }
    return index;
}

template<typename TInputImage, typename TOutputImage>
std::vector< typename BootstrapParameterEstimationImageFilter<TInputImage, TOutputImage>::TensorType > 
BootstrapParameterEstimationImageFilter<TInputImage, TOutputImage>
::SpatialBootstrapGenerator(std::vector< vnl_matrix<float> > signals)
{
    const unsigned int VectorLength = 6;
    unsigned int nb_dir = this->directions.size();
    const float min_signal = 5;
    vnl_matrix<float> signal(nb_dir-1,1);
    TensorType dt6(VectorLength,1);
    const double epsi = 1e-5;
    std::vector< TensorType > tensors;

    // loop over the bootsrap sample
    for (unsigned int c=0; c<this->m_NumberOfBootstrap; ++c) {

        // draw with replacement
        std::vector<unsigned int> sampling = RandomSpatial(nb_dir);

        // compute signal at iteration c
        for (unsigned int l=1; l<nb_dir; ++l) {
            float S0 = signals[sampling[l]](0,0);
            float Si = signals[sampling[l]](l,0);
            if (S0<min_signal) { S0=min_signal; }
            if (Si<min_signal) { Si=min_signal; }
            if (S0>=Si) { signal(l-1,0) = log(S0/Si); }
            else { signal(l-1,0) = 0.0; }
        }

        // estimate tensor
        dt6 = this->invbmatrix*signal;

        // log eucliden space
	    vnl_matrix<double> m_L(3,3);

	    m_L[0][0] =  dt6(0,0);
      	m_L[1][1] =  dt6(3,0);
      	m_L[2][2] =  dt6(5,0);
      	m_L[0][1] =  dt6(1,0);
      	m_L[0][2] =  dt6(2,0);
      	m_L[1][2] =  dt6(4,0);
      	m_L[1][0] =  m_L[0][1];
      	m_L[2][0] =  m_L[0][2];
      	m_L[2][1] =  m_L[1][2];

        vnl_symmetric_eigensystem<double> eig(m_L);

        if ( eig.D(0,0)<epsi ) { eig.D(0,0) = epsi; }
        if ( eig.D(1,1)<epsi ) { eig.D(1,1) = epsi; }
        if ( eig.D(2,2)<epsi ) { eig.D(2,2) = epsi; }
        eig.D(0,0) = log(eig.D(0,0));
        eig.D(1,1) = log(eig.D(1,1));
        eig.D(2,2) = log(eig.D(2,2));

        m_L = eig.recompose();

        dt6(0,0) = (float)m_L[0][0];
        dt6(1,0) = (float)m_L[0][1];
        dt6(2,0) = (float)m_L[0][2];
        dt6(3,0) = (float)m_L[1][1];
        dt6(4,0) = (float)m_L[1][2];
        dt6(5,0) = (float)m_L[2][2];

        tensors.push_back(dt6);
    }

    return tensors;
}

//----------------------------------------------------------------------------------------

template<typename TInputImage, typename TOutputImage>
std::vector< vnl_matrix<float> >
BootstrapParameterEstimationImageFilter<TInputImage, TOutputImage>
::LocalBootstrapInit(typename OutputImageType::IndexType idx)
{
    const float min_signal = 5;
    unsigned int nb_dir = this->directions.size();
    std::vector< vnl_matrix<float> > signals;
    vnl_matrix<float> signal(nb_dir, 1); 
    vnl_matrix<float> adc(nb_dir-1,1,0.0);

    if (this->m_Mask->GetPixel(idx)==1) {
        for (unsigned int l=0; l<nb_dir; ++l) { signal(l,0) = (float)this->GetInput(l)->GetPixel(idx); }
        // compute signal
        for (unsigned int l=1; l<nb_dir; ++l) {
            float S0 = signal(0,0);
            float Si = signal(l,0);
            if (S0<min_signal) { S0=min_signal; }
            if (Si<min_signal) { Si=min_signal; }
            if (S0>=Si) { adc(l-1,0) = log(S0/Si); }
        }
        signals.push_back(adc);
    }
    return signals;
}

template<typename TInputImage, typename TOutputImage>
std::vector<unsigned int> 
BootstrapParameterEstimationImageFilter<TInputImage, TOutputImage>
::RandomLocal(unsigned int N)
{
    static int first = 0;
    std::vector<unsigned int> index(N,0);

    for (unsigned int i=0; i<N; i++) {
        if (first==0){
            srand (time (NULL));
            first = 1;
        }
        index[i] = (unsigned int)( (rand()/(double)RAND_MAX)*(N-1.0));
    }
    return index;
}

template<typename TInputImage, typename TOutputImage>
std::vector< typename BootstrapParameterEstimationImageFilter<TInputImage, TOutputImage>::TensorType > 
BootstrapParameterEstimationImageFilter<TInputImage, TOutputImage>
::LocalBootstrapGenerator(std::vector< vnl_matrix<float> > signals)
{
    const unsigned int VectorLength = 6;

    vnl_matrix<float> adc = signals[0];
    TensorType dt6(VectorLength,1);
    const double epsi = 1e-5;
    std::vector< TensorType > tensors;

    // loop over the bootsrap sample
    for (unsigned int c=0; c<this->m_NumberOfBootstrap; ++c) {

        // draw with replacement
        std::vector<unsigned int> sampling = RandomLocal(adc.rows());

        // precompute system + organized gradients at iteration c
        BMatrixType bmatrix_c;
        BMatrixType invbmatrix_c;
        bmatrix_c.set_size(sampling.size(),VectorLength);
        for (unsigned int i=0; i<sampling.size(); i++) {					
            DirectionType bvec = this->directions[sampling[i]+1];
            bmatrix_c(i,0) = (float) this->m_BVal*bvec[0]*bvec[0];        //Dxx
            bmatrix_c(i,1) = (float) this->m_BVal*2.0*bvec[0]*bvec[1];    //Dxy
            bmatrix_c(i,2) = (float) this->m_BVal*2.0*bvec[0]*bvec[2];    //Dxz
            bmatrix_c(i,3) = (float) this->m_BVal*bvec[1]*bvec[1];        //Dyy
            bmatrix_c(i,4) = (float) this->m_BVal*2.0*bvec[1]*bvec[2];    //Dyz
            bmatrix_c(i,5) = (float) this->m_BVal*bvec[2]*bvec[2];        //Dzz
        }
        invbmatrix_c.set_size(bmatrix_c.cols(),bmatrix_c.rows()); 
        BMatrixType b1 = bmatrix_c.transpose();
        BMatrixType b2 = vnl_matrix_inverse<float>(b1*bmatrix_c);
        invbmatrix_c = b2*b1;

	    // organized adc at iteration c
	    vnl_matrix<float> adc_c(adc.rows(),adc.cols());
	    for (unsigned int l=0; l<adc.rows(); l++) { adc_c(l,0) = adc(sampling[l],0); }

        // estimate tensor
        dt6 = invbmatrix_c*adc_c;

        // log eucliden space
	    vnl_matrix<double> m_L(3,3);

	    m_L[0][0] =  dt6(0,0);
      	m_L[1][1] =  dt6(3,0);
      	m_L[2][2] =  dt6(5,0);
      	m_L[0][1] =  dt6(1,0);
      	m_L[0][2] =  dt6(2,0);
      	m_L[1][2] =  dt6(4,0);
      	m_L[1][0] =  m_L[0][1];
      	m_L[2][0] =  m_L[0][2];
      	m_L[2][1] =  m_L[1][2];

        vnl_symmetric_eigensystem<double> eig(m_L);

        if ( eig.D(0,0)<epsi ) { eig.D(0,0) = epsi; }
        if ( eig.D(1,1)<epsi ) { eig.D(1,1) = epsi; }
        if ( eig.D(2,2)<epsi ) { eig.D(2,2) = epsi; }
        eig.D(0,0) = log(eig.D(0,0));
        eig.D(1,1) = log(eig.D(1,1));
        eig.D(2,2) = log(eig.D(2,2));

        m_L = eig.recompose();

        dt6(0,0) = (float)m_L[0][0];
        dt6(1,0) = (float)m_L[0][1];
        dt6(2,0) = (float)m_L[0][2];
        dt6(3,0) = (float)m_L[1][1];
        dt6(4,0) = (float)m_L[1][2];
        dt6(5,0) = (float)m_L[2][2];

        tensors.push_back(dt6);
    }

    return tensors;
}


//----------------------------------------------------------------------------------------

template<typename TInputImage, typename TOutputImage>
void 
BootstrapParameterEstimationImageFilter<TInputImage, TOutputImage>
::ComputeParameters( std::vector< TensorType > tensors, typename OutputImageType::IndexType idx, 
                     typename OutputImageType::Pointer &output_mean, typename OutputImageType::Pointer &output_var )
{
    const unsigned int VectorLength = 6;

    // init mean and variance
    vnl_vector_fixed<float, VectorLength> temp_mean;
    temp_mean.fill(0);
    float temp_var = 0;

    // first compute the mean
    for (unsigned int i=0; i<tensors.size(); ++i) {
        TensorType tensor = tensors[i];
        for (unsigned int l=0; l<VectorLength; l++) { temp_mean[l] += (float)tensor(l,0)/(float)tensors.size(); }
    }

    // return the mean
    OutputPixelType mean = output_mean->GetPixel(idx);
    for (unsigned int l=0; l<VectorLength; l++) { mean[l] = (OutputValueType)temp_mean[l]; }

    // then compute the variance
    for (unsigned int i=0; i<tensors.size(); ++i) {
        TensorType tensor = tensors[i];
        temp_var += ((float)tensor(0,0) - temp_mean[0])*((float)tensor(0,0) - temp_mean[0])/(float)tensors.size()
                 + ((float)tensor(3,0) - temp_mean[3])*((float)tensor(3,0) - temp_mean[3])/(float)tensors.size()
                 + ((float)tensor(5,0) - temp_mean[5])*((float)tensor(5,0) - temp_mean[5])/(float)tensors.size()
                 + 2.0*((float)tensor(1,0) - temp_mean[1])*((float)tensor(1,0) - temp_mean[1])/(float)tensors.size()
                 + 2.0*((float)tensor(2,0) - temp_mean[2])*((float)tensor(2,0) - temp_mean[2])/(float)tensors.size()
                 + 2.0*((float)tensor(4,0) - temp_mean[4])*((float)tensor(4,0) - temp_mean[4])/(float)tensors.size();
    }

    // return the variance
    OutputPixelType var = output_var->GetPixel(idx);
    var[0] = (OutputValueType)temp_var;
}

//----------------------------------------------------------------------------------------

template<typename TInputImage, typename TOutputImage>
void 
BootstrapParameterEstimationImageFilter<TInputImage, TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int )
{
    typename OutputImageType::Pointer output_mean = static_cast< OutputImageType * >(this->ProcessObject::GetOutput(0));
    typename OutputImageType::Pointer output_var = static_cast< OutputImageType * >(this->ProcessObject::GetOutput(1));

    ImageRegionConstIteratorWithIndex< OutputImageType > it(output_mean, outputRegionForThread);
    it.GoToBegin();

    while( !it.IsAtEnd() ) {
        typename OutputImageType::IndexType idx = it.GetIndex();

        if (m_UseSpatialBootstrap) {
            std::vector< vnl_matrix<float> > signals = SpatialBootstrapInit(idx);
            if (signals.size()>0) {
                std::vector< TensorType > tensors = SpatialBootstrapGenerator(signals);
                ComputeParameters(tensors,idx,output_mean,output_var);
            }
        }
        else {
            std::vector< vnl_matrix<float> > signals = LocalBootstrapInit(idx);
            if (signals.size()>0) {
                std::vector< TensorType > tensors = LocalBootstrapGenerator(signals);
                ComputeParameters(tensors,idx,output_mean,output_var);
            }
        }

        ++it;
    }
}


}

#endif

							
