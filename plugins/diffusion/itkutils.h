/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef UTILS_h
#define UTILS_h

#include <itkVectorImage.h>
#include <itkImage.h>
  
void dtiInvMatrix( itk::VectorImage<float, 3>::Pointer im);
void gradient(itk::Image<float, 3>::Pointer im, itk::VectorImage<float, 3>::Pointer grad);

#endif
