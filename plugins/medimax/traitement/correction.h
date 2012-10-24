/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef __COREECTION_H
#define __COREECTION_H

#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_3d.h"
#include "noyau/imx_2d.h"
#include "noyau/gui/imx_picture3d.h"


DllExport int imx_corr_sagittal_3d_p(grphic3d *imori, grphic3d *imres);

void corr_sagittal_3d(void);
void imx_corr_sagittal_3d(int im_ori, int im_res);



#ifdef __cplusplus
}
#endif /*__cplusplus*/

#endif
