/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _8aac2d03_ab63_42a4_8848_6c4955d19365
#define _8aac2d03_ab63_42a4_8848_6c4955d19365

#include "itkJointHistogramTransferFunctionCalculator.h"

#include <algorithm>

#include <itkIdentityTransform.h>
#include <itkImageDuplicator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkMeanImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkResampleImageFilter.h>

namespace itk
{

template<typename THistogram>
void 
JointHistogramTransferFunctionCalculator<THistogram>
::Compute()
{
    TransferFunctionPointer max_prob = this->MaximumProbabilityTransferFunction();
    this->m_TransferFunction = max_prob;

    std::vector<MeasurementType> const weights = this->ConfidenceWeights();
    
    float const smooth_sigma = this->m_Histogram->GetSize(0)/50.0;
    float const outliers_sigma = (
        this->m_Histogram->GetBinMax(0,this->m_Histogram->GetSize(0)-1)-
        this->m_Histogram->GetBinMin(0,0)
    )/25.0;
    
    typedef ImageDuplicator<TransferFunctionType> ImageDuplicatorType;
    typename ImageDuplicatorType::Pointer duplicator = ImageDuplicatorType::New();
    duplicator->SetInputImage(max_prob); 
    duplicator->Update(); 
    TransferFunctionPointer result = duplicator->GetOutput();
    
    // first compute reliable but imprecise curve
	// then refine curve, using  outlier downweighting 
	for(unsigned int iter=0; iter<2; ++iter)
	{
		double sigma;
		if(iter==0) 
        { 
            sigma = smooth_sigma*5;
        }
		else if(iter==1) 
        { 
            sigma = smooth_sigma;
        }
		
		// predicts values from neighboring tangents, using weights
		// and using gaussian smoothing
        unsigned int const size = result->GetRequestedRegion().GetSize(0);
		typename TransferFunctionType::IndexType x; 
        for(x[0]=0; x[0]<size; ++x[0])
		{
			double s=0;
			double totcoef=0;
            
            typename TransferFunctionType::IndexType x0; 
			for(x0[0]=0; x0[0]<size-1; ++x0[0])
			{
                typename TransferFunctionType::IndexType x1;
                x1[0] = x0[0]+1;
                
				// tangent
				double const dx = x[0] - 0.5*(x0[0]+x1[0]);
				double const dy = max_prob->GetPixel(x1)-max_prob->GetPixel(x0);
				
                /*
                if(dy<0)
                {
                    continue;
                }
                */
                
				// confidence weights
				double const w=std::min(weights[x0[0]], weights[x1[0]]);
                
				// gaussian smoothing
				double const coef=exp(-dx*dx/(2*sigma*sigma));
				double const predictedValue=dy*(float)(x[0]-x0[0])/(float)(x1[0]-x0[0]) + max_prob->GetPixel(x0);
                
				// outlier downweighting
				double reweight=1;
				if(iter>0)
				{
					double const error = predictedValue-result->GetPixel(x);
					reweight = exp(-error*error/(2*outliers_sigma*outliers_sigma));
				}
				double const finalCoef = w*reweight*coef;

				s += finalCoef*predictedValue;
				totcoef += finalCoef;
			}
			result->SetPixel(x, totcoef>0 ? s/totcoef : 0);
		}
	}
    
    /*
    // impose increasing function, starting at max confidence
    typename std::vector<MeasurementType>::const_iterator const max_it = 
        std::max_element(weights.begin(), weights.end());
    unsigned int const wxmax = std::distance(weights.begin(), max_it);

    typename TransferFunctionType::IndexType x; 
	for(x[0]=wxmax+1; x[0] < static_cast<int>(weights.size()); ++x[0])
    {
        typename TransferFunctionType::IndexType n;
        n[0] = x[0]-1;
        result->SetPixel(x, std::max(result->GetPixel(x), result->GetPixel(n)));
    }
    for(x[0]=wxmax-1; x[0]>=0; --x[0])
    {
        typename TransferFunctionType::IndexType n;
        n[0] = x[0]+1;
        result->SetPixel(x, std::min(result->GetPixel(x), result->GetPixel(n)));
    }
    */
    
	// smooth result
    result = this->Smooth(result, 2);
    result = this->Smooth(result, 1);
    
    this->m_TransferFunction = result;
}

template<typename THistogram>
JointHistogramTransferFunctionCalculator<THistogram>
::JointHistogramTransferFunctionCalculator()
{
    this->m_Histogram = HistogramPointer(NULL);
    this->m_ResampleFactor = 10;
    this->m_TransferFunction = TransferFunctionType::New();
}

template<typename THistogram>
void
JointHistogramTransferFunctionCalculator<THistogram>
::PrintSelf(std::ostream & os, Indent indent) const
{
    this->Superclass::PrintSelf(os, indent);
}

template<typename THistogram>
typename JointHistogramTransferFunctionCalculator<THistogram>::TransferFunctionPointer
JointHistogramTransferFunctionCalculator<THistogram>
::MaximumProbabilityTransferFunction() const
{
    TransferFunctionPointer result = TransferFunctionType::New();
    
    // Region
    typename TransferFunctionType::IndexType index; index.Fill(0);
    typename TransferFunctionType::SizeType size; size.Fill(this->m_Histogram->GetSize(0));
    result->SetRegions(typename TransferFunctionType::RegionType(index, size));
    result->Allocate();
    
    // Origin
    typename TransferFunctionType::PointType origin;
    origin.Fill(this->m_Histogram->GetBinMin(0,0));
    result->SetOrigin(origin);
    
    // Spacing
    typename TransferFunctionType::SpacingType spacing;
    spacing.Fill(this->m_Histogram->GetBinMax(0,0)-this->m_Histogram->GetBinMin(0,0));
    result->SetSpacing(spacing);
    
    typedef ImageRegionIteratorWithIndex<TransferFunctionType> IteratorType;
    for(IteratorType it(result, result->GetRequestedRegion()); !it.IsAtEnd(); ++it)
    {
        int const x = it.GetIndex()[0];
        
        TransferFunctionPointer column = this->Column(x);
        TransferFunctionPointer resampled = this->Resample(column);
        
        TransferFunctionPointer smoothed = this->Smooth(resampled, 2*this->m_ResampleFactor);
        smoothed = this->Smooth(smoothed, 2*this->m_ResampleFactor);
        smoothed = this->Smooth(smoothed, 2*this->m_ResampleFactor);
        
        typedef MinimumMaximumImageCalculator<TransferFunctionType>
            MinimumMaximumImageCalculatorType;
        typename MinimumMaximumImageCalculatorType::Pointer max_calculator = 
            MinimumMaximumImageCalculatorType::New();
        max_calculator->SetImage(smoothed);
        max_calculator->ComputeMaximum();
        
        typename TransferFunctionType::PointType point;
        smoothed->TransformIndexToPhysicalPoint(max_calculator->GetIndexOfMaximum(), point);
        it.Set(point[0]);
    }
    
    return result;
}

template<typename THistogram>
std::vector<typename JointHistogramTransferFunctionCalculator<THistogram>::MeasurementType>
JointHistogramTransferFunctionCalculator<THistogram>
::ConfidenceWeights() const
{
    std::vector<MeasurementType> result;
    result.reserve(this->m_Histogram->GetSize(0));
    
    typename HistogramType::IndexType index;
    for(index[0]=0; 
		index[0] < static_cast<int>(this->m_Histogram->GetSize(0)); ++index[0])
    {
        MeasurementType total_frequency=0;
        
        for(index[1]=0; 
			index[1] < static_cast<int>(this->m_Histogram->GetSize(1)); ++index[1])
        {
            total_frequency += this->m_Histogram->GetFrequency(index);
        }
        
        result.push_back(total_frequency);
    }
    
    return result;
}

template<typename THistogram>
typename JointHistogramTransferFunctionCalculator<THistogram>::TransferFunctionPointer
JointHistogramTransferFunctionCalculator<THistogram>
::Column(int const x) const
{
    TransferFunctionPointer column = TransferFunctionType::New();
        
    // Region
    typename TransferFunctionType::IndexType index; index.Fill(0);
    typename TransferFunctionType::SizeType size; size.Fill(this->m_Histogram->GetSize(1));
    column->SetRegions(typename TransferFunctionType::RegionType(index, size));
    column->Allocate();
    
    // Origin
    typename TransferFunctionType::PointType origin;
    origin.Fill(this->m_Histogram->GetBinMin(1, 0));
    column->SetOrigin(origin);
    
    // Spacing
    typename TransferFunctionType::SpacingType spacing;
    spacing.Fill(this->m_Histogram->GetBinMax(1, 0)-this->m_Histogram->GetBinMin(1, 0));
    column->SetSpacing(spacing);
    
    typedef ImageRegionIteratorWithIndex<TransferFunctionType> IteratorType;
    for(IteratorType column_it(column, column->GetRequestedRegion()); 
        !column_it.IsAtEnd(); ++column_it)
    {
        int const y = column_it.GetIndex()[0];
        typename HistogramType::IndexType index; index[0] = x; index[1] = y;
        column_it.Set(this->m_Histogram->GetFrequency(index));
    }
    
    return column;
}

template<typename THistogram>
typename JointHistogramTransferFunctionCalculator<THistogram>::TransferFunctionPointer
JointHistogramTransferFunctionCalculator<THistogram>
::Resample(TransferFunctionPointer function) const
{
    typedef ResampleImageFilter<TransferFunctionType, TransferFunctionType>
        ResampleFilterType;
    
    typename ResampleFilterType::Pointer resample = ResampleFilterType::New();
    
    resample->SetInput(function);
    
    // Since we only wish to supersample, the transform is set to identity, and
    // the origin is unchanged while the spacing and the size are changed
    resample->SetTransform(IdentityTransform<double, 1>::New());
    
    typename TransferFunctionType::SizeType size = 
        function->GetRequestedRegion().GetSize();
    size[0] *= this->m_ResampleFactor;
    resample->SetSize(size);
    
    resample->SetOutputOrigin(function->GetOrigin());
    resample->SetOutputDirection(function->GetDirection());
    
    typename TransferFunctionType::SpacingType spacing = function->GetSpacing();
    spacing /= this->m_ResampleFactor;
    resample->SetOutputSpacing(spacing);
    
    resample->SetDefaultPixelValue(-1);
    resample->Update();
    return resample->GetOutput();
}

template<typename THistogram>
typename JointHistogramTransferFunctionCalculator<THistogram>::TransferFunctionPointer
JointHistogramTransferFunctionCalculator<THistogram>
::Smooth(TransferFunctionPointer function, int factor) const

{
    typedef MeanImageFilter<TransferFunctionType, TransferFunctionType>
        MeanImageFilterType;
    
    typename MeanImageFilterType::InputSizeType radius;
    radius.Fill(factor);
    
    typename MeanImageFilterType::Pointer mean = MeanImageFilterType::New();
    mean->SetRadius(radius);
    mean->SetInput(function);
    mean->Update();
    
    return mean->GetOutput();
}

}

#endif // _8aac2d03_ab63_42a4_8848_6c4955d19365
