/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _MeanDiffusivityImageFilter_txx
#define _MeanDiffusivityImageFilter_txx

#include "itkFractionalAnisotropyImageFilter.h"

#include "itkVectorImage.h"
#include "itkImageRegionIteratorWithIndex.h"

#include <math.h>

namespace itk
{


template<typename TInputImage, typename TOutputImage>
MeanDiffusivityImageFilter<TInputImage, TOutputImage>
::MeanDiffusivityImageFilter()
{
    this->m_EigenSystem = false;
}

template<typename TInputImage, typename TOutputImage>
void
MeanDiffusivityImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
    Superclass::PrintSelf(os,indent);
}

template<typename TInputImage, typename TOutputImage>
void
MeanDiffusivityImageFilter<TInputImage, TOutputImage>
::SetEigenValue(EigenValueImageType *val)
{
    m_EigVal = val;
    this->m_EigenSystem = true;
}

template<typename TInputImage, typename TOutputImage>
void
MeanDiffusivityImageFilter<TInputImage, TOutputImage>
::SetEigenVector(EigenVectorImageType *vec)
{
    m_EigVec = vec;
}

template<typename TInputImage, typename TOutputImage>
typename MeanDiffusivityImageFilter<TInputImage, TOutputImage>::EigenValueImageType*
MeanDiffusivityImageFilter<TInputImage, TOutputImage>
::GetEigenValue()
{
    return m_EigVal;
}

template<typename TInputImage, typename TOutputImage>
typename MeanDiffusivityImageFilter<TInputImage, TOutputImage>::EigenVectorImageType*
MeanDiffusivityImageFilter<TInputImage, TOutputImage>
::GetEigenVector()
{
    return m_EigVec;
}

template<typename TInputImage, typename TOutputImage>
void 
MeanDiffusivityImageFilter<TInputImage, TOutputImage>
::BeforeThreadedGenerateData()
{
    if (!this->m_EigenSystem) {
        typedef itk::SymmetricSpectralAnalysisImageFilter<TInputImage, TInputImage> SpectralFilter;
        typename SpectralFilter::Pointer filter = SpectralFilter::New();
        filter->SetInput(0,this->GetInput(0));
        filter->Update();
        m_EigVal = filter->GetOutput(0);
        m_EigVec = filter->GetOutput(1);
        this->m_EigenSystem = true;
    }
}

template<typename TInputImage, typename TOutputImage>
void 
MeanDiffusivityImageFilter<TInputImage, TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int )
{
    typename OutputImageType::Pointer output = static_cast< OutputImageType * >(this->ProcessObject::GetOutput(0));

    ImageRegionConstIteratorWithIndex< OutputImageType > it(output, outputRegionForThread);
    it.GoToBegin();

    while( !it.IsAtEnd() ) {
        InputPixelType eigenvalues = m_EigVal->GetPixel(it.GetIndex());
        double ev1 = (double) eigenvalues[0];
        double ev2 = (double) eigenvalues[1];
        double ev3 = (double) eigenvalues[2];
        double md = (ev1+ev2+ev3)/3.0;
        this->GetOutput(0)->SetPixel(it.GetIndex(),(OutputPixelType)md);

        ++it;
    }

}

}

#endif

							
