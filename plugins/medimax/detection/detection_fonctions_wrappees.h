/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef DETECTION_FONCTIONS_WRAPPEE_H
#define DETECTION_FONCTIONS_WRAPPEE_H

#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/
#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_3d.h"
#include "noyau/imx_2d.h"


DllExport int ComputeRocCurve(grphic3d *imdetection,grphic3d *maskdetection,grphic3d *groundtruth, int nbpoint,char *nomfichier);


#ifdef __cplusplus
}
#endif /*__cplusplus*/

#endif
