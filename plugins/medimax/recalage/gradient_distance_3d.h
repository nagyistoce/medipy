/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*******************************************************************************/
/*!	   \file:		gradient_distance_3d.h
***
***		project:	Imagix 2.01
***
***
***		\brief description:    
***
***
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***
********************************************************************************/

#ifndef __GRADIENT_DISTANCE_3D_H__
#define __GRADIENT_DISTANCE_3D_H__

#include "recalage/definitions_types_recalage.h"

/*********** GRADIENT DES TRANSFORMATIONS ERREUR QUADRATIQUE ******************/
extern int	gradient_rigid_quad_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP);
extern int	gradient_rigid_quad2_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP);
extern int	gradient_rigidz_quad_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP);
extern int	gradient_rigidz_quad2_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP);
extern int	gradient_affine_quad_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP);
extern int	gradient_affine_quad2_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP);
extern int	gradient_base_quad_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP);
extern int  gradient_base_quad2_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP);
/******************** GRADIENT INFORMATION MUTUELLE ***************************/
extern int	gradient_rigid_entropie_conjointe_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP);
extern int	gradient_rigid_IM_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP);
extern int	gradient_rigid_IMNS_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP);
extern int	gradient_rigid_IMNM1_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP);
extern int	gradient_rigid_IMNM2_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP);
extern int	gradient_rigidz_entropie_conjointe_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP);
extern int	gradient_rigidz_IM_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP);
extern int	gradient_rigidz_IMNS_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP);
extern int	gradient_rigidz_IMNM1_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP);
extern int	gradient_rigidz_IMNM2_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP);
extern int	gradient_affine_entropie_conjointe_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP);
extern int	gradient_affine_IM_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP);
extern int	gradient_affine_IMNS_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP);
extern int	gradient_affine_IMNM1_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP);
extern int	gradient_affine_IMNM2_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP);

int gradient_entropie_conjointe_3d(double nbVoxelsOverlap, VP_histos *donnees_VP, double **dEnergie);
int gradient_IM_3d(double nbVoxelsOverlap, VP_histos *donnees_VP, double **dEnergie);
int gradient_IMNS_3d(double nbVoxelsOverlap, VP_histos *donnees_VP, double **dEnergie);
int gradient_IMNM1_3d(double nbVoxelsOverlap, VP_histos *donnees_VP, double **dEnergie);
int gradient_IMNM2_3d(double nbVoxelsOverlap, VP_histos *donnees_VP, double **dEnergie);

int calcul_gradient_VP_3d(IM3DTYPE im_type, TRANSFO3D transfo, grphic3d *Imreca, grphic3d *Imref, int nb_param, double *param, double *grad, VP_histos *donnees_VP);
void calcul_transformation_derivee_rigid_3d(double *param, double *** DTRANS);
void calcul_transformation_derivee_rigid_zoom_3d(double *param, double *** DTRANS);
void calcul_transformation_derivee_affine_3d(double *param, double *** DTRANS);

extern gradient_func_t imx_choose_gradient_fct(TRANSFO3D transfo_type, dist_func_t dist);

#endif /*__GRADIENT_DISTANCE_3D_H__*/
