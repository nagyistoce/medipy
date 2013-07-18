/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef FONCTIONS_WRAPPEE_TRAITEMENT_H
#define FONCTIONS_WRAPPEE_TRAITEMENT_H


#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_3d.h"

DllExport void IntensityNormalisationForRegisteredImages(grphic3d *imsrc, grphic3d *imref, grphic3d *imres, grphic3d *mask_imsrc, grphic3d * mask_imref, int method);


#ifdef __cplusplus
}
#endif /*__cplusplus*/

#endif
