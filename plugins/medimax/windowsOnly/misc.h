/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef __MISC_H
#define __MISC_H

#include <config.h>

#include "noyau/imagix.h"
#include "noyau/imx_3d.h"

extern int save_mri_ipb_3d_p(char *file, grphic3d *image) ;

#ifdef __cplusplus
extern "C" {
#endif

/* http://www.eecg.utoronto.ca/~aamodt/sourceware/MSVC.html */
double rint( double x);

int random();
double drand48();

#ifdef __cplusplus
}
#endif

#endif
