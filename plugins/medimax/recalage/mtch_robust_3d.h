/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*******************************************************************************/	
/*!	   \file:		mtch_robust_3d.h
***
***		project:	Imagix 2.01
***			
***
***		\brief description:    Fichier pourle matching robust des images 3D
***	
***	
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***
********************************************************************************/

#ifndef __MTCH_ROBUST_3D_H
#define __MTCH_ROBUST_3D_H

extern void calc_robust_stat_3d(grphic3d *imref, grphic3d *imreca, int energy_function, int robust_estimator);
extern double calc_moy_image_3d(grphic3d *image);
extern double calc_ecart_type_image_3d(grphic3d *image);
extern int  define_noise_param_3d(grphic3d *im1, grphic3d *im2, int energy_function, double *moy, double *sigma);

// Variables globales
extern  double _moy_imref_3d;
extern double _moy_imreca_3d;

//extern double _ect_imref_3d;
//extern double _ect_imreca_3d;


#endif /*__MTCH_ROBUST_3D_H*/

