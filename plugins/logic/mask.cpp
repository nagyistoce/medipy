/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "mask.h"

#include <itkImage.h>

#include "itkMaskWithValueImageFilter.h"

void
apply_mask(itk::Image<float, 3>::Pointer input,
           itk::Image<float, 3>::Pointer mask,
           float background, float outside,
           itk::Image<float, 3>::Pointer output)
{
    typedef itk::Image<float, 3> ImageType;

    typedef
        itk::MaskWithValueImageFilter<ImageType, ImageType, ImageType>
        MaskWithValueFilterType;

    MaskWithValueFilterType::Pointer maskFilter = MaskWithValueFilterType::New();

    maskFilter->SetInput1(input);
    maskFilter->SetInput2(mask);
    maskFilter->SetBackgroundValue(background);
    maskFilter->SetOutsideValue(outside);
    maskFilter->Update();

    output->SetRegions(maskFilter->GetOutput()->GetRequestedRegion());
    output->SetPixelContainer(maskFilter->GetOutput()->GetPixelContainer());
    output->CopyInformation(input);
}
