/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _FractionalAnisotropyImageFilter_txx
#define _FractionalAnisotropyImageFilter_txx

#include "itkFractionalAnisotropyImageFilter.h"

#include "itkVectorImage.h"
#include "itkImageRegionIteratorWithIndex.h"

#include <math.h>

namespace itk
{


template<typename TInputImage, typename TOutputImage>
FractionalAnisotropyImageFilter<TInputImage, TOutputImage>
::FractionalAnisotropyImageFilter()
{
    this->m_EigenSystem = false;
}

template<typename TInputImage, typename TOutputImage>
void
FractionalAnisotropyImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
    Superclass::PrintSelf(os,indent);
}

template<typename TInputImage, typename TOutputImage>
void
FractionalAnisotropyImageFilter<TInputImage, TOutputImage>
::SetEigenValue(EigenValueImageType val)
{
    m_EigVal = val;
    this->m_EigenSystem = true;
}

template<typename TInputImage, typename TOutputImage>
void
FractionalAnisotropyImageFilter<TInputImage, TOutputImage>
::SetEigenVector(EigenVectorImageType vec)
{
    m_EigVec = vec;
}

template<typename TInputImage, typename TOutputImage>
typename FractionalAnisotropyImageFilter<TInputImage, TOutputImage>::EigenValueImageType
FractionalAnisotropyImageFilter<TInputImage, TOutputImage>
::GetEigenValue()
{
    return m_EigVal;
}

template<typename TInputImage, typename TOutputImage>
typename FractionalAnisotropyImageFilter<TInputImage, TOutputImage>::EigenVectorImageType
FractionalAnisotropyImageFilter<TInputImage, TOutputImage>
::GetEigenVector()
{
    return m_EigVec;
}

template<typename TInputImage, typename TOutputImage>
void 
FractionalAnisotropyImageFilter<TInputImage, TOutputImage>
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
FractionalAnisotropyImageFilter<TInputImage, TOutputImage>
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
        double fa = sqrt(0.5 * ((ev1-ev2)*(ev1-ev2) + (ev2-ev3)*(ev2-ev3) + (ev3-ev1)*(ev3-ev1)) / (ev1*ev1 + ev2*ev2 + ev3*ev3));
        this->GetOutput(0)->SetPixel(it.GetIndex(),(OutputPixelType)fa);

        ++it;
    }

}

}

#endif

							
