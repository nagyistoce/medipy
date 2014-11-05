/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "spectral_analysis.h"
#include "itkSymmetricSpectralAnalysisImageFilter.h"

void spectral_analysis(itk::VectorImage<float, 3>::Pointer dt6,
                       itk::VectorImage<float, 3>::Pointer eigval,
                       itk::VectorImage<float, 3>::Pointer eigvec)
{
    typedef itk::VectorImage<float, 3> VectorImage;

    typedef itk::SymmetricSpectralAnalysisImageFilter<VectorImage, VectorImage> Filter;
    Filter::Pointer filter = Filter::New();
    filter->GraftNthOutput(0,eigval);
    filter->GraftNthOutput(1,eigvec);
    filter->SetInput(0,dt6);
    filter->Update();
}
