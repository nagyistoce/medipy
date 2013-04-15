/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/**-----------------------------------------------------------------------
***	
***	file:		otsu_3d.h
***
***	project:	Imagix 1.01
***			
***
***	description:	Header file for project imagix (1.01)
***	
***	
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***     Last user action by: Mr. ARMSPACH on Juny  1995
***
***---------------------------------------------------------------------*/

#ifndef OTSU_3D_H
#define OTSU_3D_H

#include <config.h>

#ifdef __cplusplus 
extern "C" 
{
#endif /*__cplusplus*/ 

#include "noyau/imx_3d.h"

extern void otsu_thresholding_3d(void) ;
extern void otsu_brain_3d(void) ;
extern void otsu_ventricules_3d(void) ;
extern void imx_otsu_thresholding_3d(int im_deb, int im_res, int thr_nb);
extern DllExport void imx_otsu_thresholding_3d_p(grphic3d *imdeb, grphic3d *imres, int thr_nb);
extern void imx_otsu_brain_3d(int im_deb, int im_res);
extern void imx_otsu_brain_3d_p(grphic3d *imdeb, grphic3d *imres);
extern void imx_otsu_ventricules_3d(int im_deb, int im_res);
extern void imx_otsu_ventricules_3d_p(grphic3d *imdeb, grphic3d *imres);
extern void find_threshold_3d(grphic3d *imdeb, int *threshold, int thr_nb, BOOL bWithZeroValue);
extern void calc_prob_table_3d(grphic3d *imdeb, double *prob, BOOL bWithZeroValue);

#ifdef __cplusplus 
}
#endif /*__cplusplus*/ 

#endif
