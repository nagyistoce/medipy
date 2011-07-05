/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef medipy_components_logic_mask_h
#define medipy_components_logic_mask_h

#include <itkImage.h>

void
apply_mask(itk::Image<float, 3>::Pointer input,
           itk::Image<float, 3>::Pointer mask,
           float background, float outside,
           itk::Image<float, 3>::Pointer output);

#endif // medipy_components_logic_mask_h
