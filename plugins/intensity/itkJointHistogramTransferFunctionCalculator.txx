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
#include <iterator>

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
    std::vector<MeasurementType> const weights = this->ConfidenceWeights();
    
    float const smooth_sigma = this->m_Histogram->GetSize(0)/50.0;
    float const outliers_sigma = (
        this->m_Histogram->GetBinMax(0,this->m_Histogram->GetSize(0)-1)-
        this->m_Histogram->GetBinMin(0,0)
    )/25.0;
    
    TransferFunctionPointer result = TransferFunctionPointer::New();
    result->DeepCopy(max_prob);
    
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
        for(int x=0; x<result->GetSize(); ++x)
		{
			double s=0;
			double totcoef=0;
            
            for(int x0=0; x0<result->GetSize()-1; ++x0)
			{
                int const x1 = x0+1;
                
                double node_value[4];
                
                max_prob->GetNodeValue(x0, node_value);
                double const y0 = node_value[1];
                
                max_prob->GetNodeValue(x1, node_value);
                double const y1 = node_value[1];
                
				// tangent
				double const dx = x - 0.5*(x0+x1);
				double const dy = y1-y0;
				
                /*
                if(dy<0)
                {
                    continue;
                }
                */
                
				// confidence weights
				double const w=std::min(weights[x0], weights[x1]);
                
				// gaussian smoothing
				double const coef=exp(-dx*dx/(2*sigma*sigma));
                double const predictedValue=dy*(float)(x-x0)/(float)(x1-x0) + y0;
                
				// outlier downweighting
				double reweight=1;
				if(iter>0)
				{
                    result->GetNodeValue(x, node_value);
                    double const value = node_value[1];
                    
					double const error = predictedValue-value;
					reweight = exp(-error*error/(2*outliers_sigma*outliers_sigma));
				}
				double const finalCoef = w*reweight*coef;

				s += finalCoef*predictedValue;
				totcoef += finalCoef;
			}
            
            double node_value[4];                
            result->GetNodeValue(x, node_value);
            
            node_value[1] = totcoef>0 ? s/totcoef : 0;
            result->SetNodeValue(x, node_value);
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
typename JointHistogramTransferFunctionCalculator<THistogram>::TransferFunctionType const * 
JointHistogramTransferFunctionCalculator<THistogram>
::GetTransferFunction() const
{
    return this->m_TransferFunction;
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
    
    for(unsigned int i=0; i<this->m_Histogram->GetSize(0); ++i)
    {
        float const x = 0.5*(
            this->m_Histogram->GetBinMin(0, i)+
            this->m_Histogram->GetBinMax(0, i));
        
        std::vector<MeasurementType> const column = this->Column(i);
        std::vector<MeasurementType> const resampled = this->Resample(column);
        std::vector<MeasurementType> smoothed = 
            this->Smooth(resampled, 2*this->m_ResampleFactor);
        smoothed = this->Smooth(smoothed, 2*this->m_ResampleFactor);
        smoothed = this->Smooth(smoothed, 2*this->m_ResampleFactor);
        
        unsigned int const index_of_maximum = std::distance(smoothed.begin(),
            std::max_element(smoothed.begin(), smoothed.end()));
        
        float const index_of_maximum_original = 
            float(index_of_maximum)/float(this->m_ResampleFactor);
        
        float const alpha = index_of_maximum_original-int(index_of_maximum_original);
        
        float const y = 
            (1-alpha)*this->m_Histogram->GetBinMin(1, int(index_of_maximum_original))+
            alpha*this->m_Histogram->GetBinMax(1, int(index_of_maximum_original));
        
        result->AddPoint(x, y);
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
std::vector<typename JointHistogramTransferFunctionCalculator<THistogram>::MeasurementType>
JointHistogramTransferFunctionCalculator<THistogram>
::Column(int const x) const
{
    std::vector<MeasurementType> column(this->m_Histogram->GetSize(1));
    
    for(unsigned int y=0; y<this->m_Histogram->GetSize(1); ++y)
    {
        typename HistogramType::IndexType index; index[0] = x; index[1] = y;
        column[y] = this->m_Histogram->GetFrequency(index);
    }
    
    return column;
}

template<typename THistogram>
std::vector<typename JointHistogramTransferFunctionCalculator<THistogram>::MeasurementType>
JointHistogramTransferFunctionCalculator<THistogram>
::Resample(std::vector<MeasurementType> const & function) const
{
    std::vector<MeasurementType> resampled(this->m_ResampleFactor*function.size());
    
    for(unsigned int i=0; i<resampled.size(); ++i)
    {
        float p = float(i)/float(this->m_ResampleFactor);
        
        MeasurementType const before = function[
            std::min<int>(function.size()-1, int(p))];
        MeasurementType const after = function[
            std::min<int>(function.size()-1, 1+int(p))];
        
        float const alpha = p-int(p);
        float const value = (1-alpha)*before + alpha*after;
        resampled[i] = static_cast<MeasurementType>(value);
    }
    
    return resampled;
}

template<typename THistogram>
std::vector<typename JointHistogramTransferFunctionCalculator<THistogram>::MeasurementType>
JointHistogramTransferFunctionCalculator<THistogram>
::Smooth(std::vector<MeasurementType> const & function, int radius) const
{
    std::vector<MeasurementType> smoothed(function.size());
    for(unsigned int i=0; i<smoothed.size(); ++i)
    {
        float value=0;
        unsigned int count=0;
        for(int j=std::max<int>(0, i-radius); 
            j < std::min<int>(function.size()-1, i+radius); ++j)
        {
            value += function[j];
            ++count;
        }
        if(count != 0)
        {
            value /= count;
        }
        
        smoothed[i] = value;
    }
    
    return smoothed;
}

template<typename THistogram>
typename JointHistogramTransferFunctionCalculator<THistogram>::TransferFunctionPointer
JointHistogramTransferFunctionCalculator<THistogram>
::Smooth(TransferFunctionPointer function, int factor) const
{
    std::vector<MeasurementType> vector(function->GetSize());
    for(int i=0; i<function->GetSize(); ++i)
    {
        double node_value[4];
        function->GetNodeValue(i, node_value);
        vector[i] = node_value[1];
    }
    
    vector = this->Smooth(vector, factor);
    
    TransferFunctionPointer smoothed = TransferFunctionPointer::New();
    for(int i=0; i<function->GetSize(); ++i)
    {
        double node_value[4];
        function->GetNodeValue(i, node_value);
        
        smoothed->AddPoint(node_value[0], vector[i]);
    }
    
    return smoothed;
}

}

#endif // _8aac2d03_ab63_42a4_8848_6c4955d19365
