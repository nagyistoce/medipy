#ifndef _66088227_effb_471e_98a9_bb891ea79ad9
#define _66088227_effb_471e_98a9_bb891ea79ad9

#include "itkBootstrapDWIStatisticsImageFilter.h"

#include <stdlib.h>
#include <time.h>
#include <vector>

#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkNeighborhood.h>

#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

namespace itk
{

template<typename TInputImage, typename TMeanImage,
         typename TStandardDeviationImage, typename TMaskImage>
typename BootstrapDWIStatisticsImageFilter<
    TInputImage, TMeanImage, TStandardDeviationImage, TMaskImage>::DirectionType 
BootstrapDWIStatisticsImageFilter<TInputImage, TMeanImage,
                                  TStandardDeviationImage, TMaskImage>
::GetGradientDirection(unsigned int i)
{
    if (i<this->directions.size()) {
        return this->directions[i];
    }
    else {
        throw "Wrong index to access gradients!";
    }
}

template<typename TInputImage, typename TMeanImage,
         typename TStandardDeviationImage, typename TMaskImage>
void 
BootstrapDWIStatisticsImageFilter<TInputImage, TMeanImage,
                                  TStandardDeviationImage, TMaskImage>
::SetGradientDirection(unsigned int i, DirectionType bvec)
{
    this->directions.resize(std::max<unsigned int>(this->directions.size(), i+1));
    this->directions[i] = bvec;
}

template<typename TInputImage, typename TMeanImage,
         typename TStandardDeviationImage, typename TMaskImage>
BootstrapDWIStatisticsImageFilter<TInputImage, TMeanImage,
                                  TStandardDeviationImage, TMaskImage>
::BootstrapDWIStatisticsImageFilter()
{
    this->SetSizePlane(3);
    this->SetSizeDepth(3);
    this->SetSamplesCount(200);
    this->UseSpatialBootstrapOn();
}

template<typename TInputImage, typename TMeanImage,
         typename TStandardDeviationImage, typename TMaskImage>
void 
BootstrapDWIStatisticsImageFilter<TInputImage, TMeanImage,
                                  TStandardDeviationImage, TMaskImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
    this->Superclass::PrintSelf(os,indent);
    os << indent << "Size in plane: " << this->GetSizePlane() << "\n";
    os << indent << "Size in depth: " << this->GetSizeDepth() << "\n";
    os << indent << "Samples count: " << this->GetSamplesCount() << "\n";
    os << indent << "Use spatial bootstrap: " 
                 << std::boolalpha << this->GetUseSpatialBootstrap() << "\n";
}

template<typename TInputImage, typename TMeanImage,
         typename TStandardDeviationImage, typename TMaskImage>
void
BootstrapDWIStatisticsImageFilter<TInputImage, TMeanImage,
                                  TStandardDeviationImage, TMaskImage>
::BeforeThreadedGenerateData()
{
    typename TMeanImage::Pointer mean_image = 
        dynamic_cast<TMeanImage*>(this->ProcessObject::GetOutput(0));
    typename TMeanImage::PixelType zero(6);
    zero.Fill(0);
    mean_image->FillBuffer(zero);
    
    typename TStandardDeviationImage::Pointer standard_deviation_image = 
        dynamic_cast<TStandardDeviationImage*>(this->ProcessObject::GetOutput(1));
    standard_deviation_image->FillBuffer(0);
    
    unsigned int const VectorLength = 6;
    unsigned int const nb_dir = this->directions.size();

    this->bmatrix.set_size(nb_dir-1,VectorLength);
    for (unsigned int i=1; i<nb_dir; i++) 
    {					
        DirectionType const & bvec = this->directions[i];
        this->bmatrix(i-1,0) = (float) this->m_BValue*bvec[0]*bvec[0];        //Dxx
        this->bmatrix(i-1,1) = (float) this->m_BValue*2.0*bvec[0]*bvec[1];    //Dxy
        this->bmatrix(i-1,2) = (float) this->m_BValue*2.0*bvec[0]*bvec[2];    //Dxz
        this->bmatrix(i-1,3) = (float) this->m_BValue*bvec[1]*bvec[1];        //Dyy
        this->bmatrix(i-1,4) = (float) this->m_BValue*2.0*bvec[1]*bvec[2];    //Dyz
        this->bmatrix(i-1,5) = (float) this->m_BValue*bvec[2]*bvec[2];        //Dzz
    }
    
    this->invbmatrix.set_size(this->bmatrix.cols(),this->bmatrix.rows()); 
    BMatrixType b1 = this->bmatrix.transpose();
    BMatrixType b2 = vnl_matrix_inverse<float>(b1*this->bmatrix);
    this->invbmatrix = b2*b1;
}

template<typename TInputImage, typename TMeanImage,
         typename TStandardDeviationImage, typename TMaskImage>
void 
BootstrapDWIStatisticsImageFilter<TInputImage, TMeanImage,
                                  TStandardDeviationImage, TMaskImage>
::ThreadedGenerateData(OutputImageRegionType const & outputRegionForThread, ThreadIdType)
{
    MaskImageConstPointer mask_image = this->GetMaskImage();
    MeanImagePointer mean_image = 
        dynamic_cast<TMeanImage*>(this->ProcessObject::GetOutput(0));
    StandardDeviationImagePointer standard_deviation_image = 
        dynamic_cast<TStandardDeviationImage*>(this->ProcessObject::GetOutput(1));

    BootstrapInitializer initializer = NULL;
    BootstrapGenerator generator = NULL;
    if(m_UseSpatialBootstrap)
    {
        initializer = &Self::SpatialBootstrapInit;
        generator = &Self::SpatialBootstrapGenerator;
    }
    else
    {
        initializer = &Self::LocalBootstrapInit;
        generator = &Self::LocalBootstrapGenerator;
    }

    typedef ImageRegionConstIteratorWithIndex<MeanImageType> IteratorType;
    for(IteratorType it(mean_image, outputRegionForThread); !it.IsAtEnd(); ++it)
    {
        typename MeanImageType::IndexType const & idx = it.GetIndex();
        
        bool process=true;
        if(!mask_image.IsNull()) 
        {
            typename TInputImage::PointType point; 
            mean_image->TransformIndexToPhysicalPoint(idx, point);
            
            typename TMaskImage::IndexType mask_index;
            mask_image->TransformPhysicalPointToIndex(point, mask_index);
            
            if(!mask_image->GetLargestPossibleRegion().IsInside(mask_index) || 
               mask_image->GetPixel(mask_index) == 0)
            {
                process=false;
            }
        }

        if(process)
        {
            std::vector<SignalType> const signals = (this->*initializer)(idx);
            if(signals.size()>0) 
            {
                std::vector<TensorType> const tensors = (this->*generator)(signals);
                this->ComputeParameters(tensors,idx);
            }
        }
    }
}

template<typename TInputImage, typename TMeanImage,
         typename TStandardDeviationImage, typename TMaskImage>
std::vector<unsigned int> 
BootstrapDWIStatisticsImageFilter<TInputImage, TMeanImage,
                                  TStandardDeviationImage, TMaskImage>
::Random(unsigned int size, unsigned int max_value) const
{
    static bool first = true;
    if(first)
    {
        srand(time(NULL));
        first = false;
    }
    
    std::vector<unsigned int> indices(size);
    for (unsigned int i=0; i<size; i++) 
    {
        indices[i] = float(rand())/float(RAND_MAX)*float(max_value);
    }
    return indices;
}

template<typename TInputImage, typename TMeanImage,
         typename TStandardDeviationImage, typename TMaskImage>
std::vector<typename BootstrapDWIStatisticsImageFilter<
    TInputImage, TMeanImage, TStandardDeviationImage, TMaskImage>::SignalType> 
BootstrapDWIStatisticsImageFilter<TInputImage, TMeanImage,
                                  TStandardDeviationImage, TMaskImage>
::SpatialBootstrapInit(MeanImageIndexType idx)
{
    typename MaskImageType::ConstPointer mask_image = this->GetMaskImage();
    
    typename TInputImage::RegionType const input_region = 
        this->GetInput(0)->GetLargestPossibleRegion();
    
    // Build the neighborhood
    typedef itk::Neighborhood<typename TInputImage::PixelType, 
                              TInputImage::ImageDimension> NeighborhoodType;
    NeighborhoodType neighborhood;
    
    typename NeighborhoodType::SizeType neighborhood_size;
    neighborhood_size[0]=(this->m_SizePlane-1)/2; 
    neighborhood_size[1]=(this->m_SizePlane-1)/2; 
    neighborhood_size[2]=(this->m_SizeDepth-1)/2;
    neighborhood.SetRadius(neighborhood_size);
    
    std::vector<SignalType> signals;
    SignalType signal(this->directions.size(), 1); 

    bool process_pixel = true;
    for(unsigned int i=0; i<neighborhood.Size() && process_pixel; ++i)
    {
        if(!input_region.IsInside(idx+neighborhood.GetOffset(i)))
        {
            process_pixel = false;
        }
    }

    if(process_pixel)
    {
        for(unsigned int i=0; i<neighborhood.Size(); ++i)
        {
            for(unsigned int l=0; l<this->directions.size(); ++l) 
            {
                signal(l, 0) = this->GetInput(l)->GetPixel(idx+neighborhood.GetOffset(i));
            }
            signals.push_back(signal);
        }
    }
    
    return signals;
}

template<typename TInputImage, typename TMeanImage,
         typename TStandardDeviationImage, typename TMaskImage>
std::vector<typename BootstrapDWIStatisticsImageFilter<
    TInputImage, TMeanImage, TStandardDeviationImage, TMaskImage>::TensorType> 
BootstrapDWIStatisticsImageFilter<TInputImage, TMeanImage,
                                  TStandardDeviationImage, TMaskImage>
::SpatialBootstrapGenerator(std::vector<SignalType> signals)
{
    unsigned int const VectorLength = 6;
    unsigned int const nb_dir = this->directions.size();
    float const min_signal = 5;
    double const epsi = 1e-5;
    
    SignalType signal(nb_dir-1,1);
    TensorType dt6(VectorLength,1);
    std::vector<TensorType> tensors;

    unsigned int const neighborhood_size = 
        this->m_SizePlane*this->m_SizePlane*this->m_SizeDepth;

    // loop over the bootsrap sample
    for(unsigned int c=0; c<this->m_SamplesCount; ++c) {

        // draw with replacement
        std::vector<unsigned int> const sampling = 
            this->Random(nb_dir, neighborhood_size-1);

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

template<typename TInputImage, typename TMeanImage,
         typename TStandardDeviationImage, typename TMaskImage>
std::vector<typename BootstrapDWIStatisticsImageFilter<
    TInputImage, TMeanImage, TStandardDeviationImage, TMaskImage>::SignalType> 
BootstrapDWIStatisticsImageFilter<TInputImage, TMeanImage,
                                  TStandardDeviationImage, TMaskImage>
::LocalBootstrapInit(MeanImageIndexType idx)
{
    float const min_signal = 5;
    unsigned int nb_dir = this->directions.size();
    std::vector<SignalType> signals;
    SignalType signal(nb_dir, 1); 
    vnl_matrix<float> adc(nb_dir-1,1,0.0);

    for (unsigned int l=0; l<nb_dir; ++l) 
    { 
        signal(l,0) = this->GetInput(l)->GetPixel(idx); 
    }
    // compute signal
    for (unsigned int l=1; l<nb_dir; ++l) 
    {
        float S0 = signal(0,0);
        float Si = signal(l,0);
        if (S0<min_signal) { S0=min_signal; }
        if (Si<min_signal) { Si=min_signal; }
        if (S0>=Si) { adc(l-1,0) = log(S0/Si); }
    }
    signals.push_back(adc);

    return signals;
}

template<typename TInputImage, typename TMeanImage,
         typename TStandardDeviationImage, typename TMaskImage>
std::vector<typename BootstrapDWIStatisticsImageFilter<
    TInputImage, TMeanImage, TStandardDeviationImage, TMaskImage>::TensorType> 
BootstrapDWIStatisticsImageFilter<TInputImage, TMeanImage,
                                  TStandardDeviationImage, TMaskImage>
::LocalBootstrapGenerator(std::vector<SignalType> signals)
{
    const unsigned int VectorLength = 6;

    vnl_matrix<float> adc = signals[0];
    TensorType dt6(VectorLength,1);
    const double epsi = 1e-5;
    std::vector< TensorType > tensors;

    // loop over the bootsrap sample
    for (unsigned int c=0; c<this->m_SamplesCount; ++c) {

        // draw with replacement
        std::vector<unsigned int> const sampling = 
            this->Random(adc.rows(), adc.rows()-1);

        // precompute system + organized gradients at iteration c
        BMatrixType bmatrix_c;
        BMatrixType invbmatrix_c;
        bmatrix_c.set_size(sampling.size(),VectorLength);
        for (unsigned int i=0; i<sampling.size(); i++) {					
            DirectionType bvec = this->directions[sampling[i]+1];
            bmatrix_c(i,0) = (float) this->m_BValue*bvec[0]*bvec[0];        //Dxx
            bmatrix_c(i,1) = (float) this->m_BValue*2.0*bvec[0]*bvec[1];    //Dxy
            bmatrix_c(i,2) = (float) this->m_BValue*2.0*bvec[0]*bvec[2];    //Dxz
            bmatrix_c(i,3) = (float) this->m_BValue*bvec[1]*bvec[1];        //Dyy
            bmatrix_c(i,4) = (float) this->m_BValue*2.0*bvec[1]*bvec[2];    //Dyz
            bmatrix_c(i,5) = (float) this->m_BValue*bvec[2]*bvec[2];        //Dzz
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

template<typename TInputImage, typename TMeanImage,
         typename TStandardDeviationImage, typename TMaskImage>
void 
BootstrapDWIStatisticsImageFilter<TInputImage, TMeanImage,
                                  TStandardDeviationImage, TMaskImage>
::ComputeParameters(std::vector<TensorType> const & tensors, MeanImageIndexType idx)
{
    typedef typename TMeanImage::PixelType OutputTensorType;
    
    typename TMeanImage::Pointer mean_image = 
        dynamic_cast<TMeanImage*>(this->ProcessObject::GetOutput(0));
    typename TStandardDeviationImage::Pointer standard_deviation_image = 
        dynamic_cast<TStandardDeviationImage*>(this->ProcessObject::GetOutput(1));

    // first compute the mean
    OutputTensorType mean = mean_image->GetPixel(idx);
    mean.Fill(0.);
    for(typename std::vector<TensorType>::const_iterator it=tensors.begin(); 
        it!=tensors.end(); ++it)
    {
        for(unsigned int i=0; i<mean.Size(); ++i)
        {
            mean[i] += (*it)(i,0)/tensors.size();
        }
    }
    
    // then compute the variance
    typename TStandardDeviationImage::PixelType variance = 0.;
    for(typename std::vector<TensorType>::const_iterator it=tensors.begin(); 
        it!=tensors.end(); ++it)
    {
        TensorType const & tensor = *it;
        variance += (
            (tensor(0,0) - mean[0])*(tensor(0,0) - mean[0]) +
            (tensor(3,0) - mean[3])*(tensor(3,0) - mean[3]) +
            (tensor(5,0) - mean[5])*(tensor(5,0) - mean[5]) +
            2.0*(tensor(1,0) - mean[1])*(tensor(1,0) - mean[1]) +
            2.0*(tensor(2,0) - mean[2])*(tensor(2,0) - mean[2]) +
            2.0*(tensor(4,0) - mean[4])*(tensor(4,0) - mean[4])
        )/tensors.size();
    }
    // Clamp variance and compute standard deviation
    variance = (variance>=0)?variance:0;
    standard_deviation_image->SetPixel(idx, variance/6.0);
}

}

#endif // _66088227_effb_471e_98a9_bb891ea79ad9
