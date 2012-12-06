/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _bec069a6_9b54_4beb_bb22_47ac2fa157ac
#define _bec069a6_9b54_4beb_bb22_47ac2fa157ac

#include <itkImage.h>
#include <itkVectorImage.h>

void spectral_analysis(itk::VectorImage<float, 3>::Pointer dt6,
                       itk::VectorImage<float, 3>::Pointer eigval,
                       itk::VectorImage<float, 3>::Pointer eigvec);

#endif // _bec069a6_9b54_4beb_bb22_47ac2fa157ac
