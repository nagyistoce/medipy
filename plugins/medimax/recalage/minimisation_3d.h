/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*******************************************************************************/
/*!	   \file:		minimisation_3d.h
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

#ifndef __MINIMISATION_3D_H__
#define __MINIMISATION_3D_H__

#include "recalage/definitions_types_recalage.h"
                              
/************************ MINIMISATIONS ***************************************/
extern double fonc_transfo_3d(grphic3d *imref, grphic3d *imreca, grphic3d *imres, transf_func_t transf, InterpolationFct inter, dist_func_t dist, contrainte_func_t contrainte, field3d *champ_ini, field3d *champ_fin, int nb, double *param, VP_histos *donnees_VP, field3d *points_calcules, double ***weights, ERREUR_RECALAGE *err);
extern double energie_transfo_3d(ptr_minimisation_param min_parameters);
extern double linemin_3d(ptr_minimisation_param min_parameters, double *direc);

extern double min_icm_3d(ptr_minimisation_param min_parameters);
extern double min_powel_3d(ptr_minimisation_param min_parameters);
extern double min_simplex_3d(ptr_minimisation_param min_parameters);
extern double min_desc_grad_3d(ptr_minimisation_param min_parameters);
extern double min_grad_conj_3d(ptr_minimisation_param min_parameters);
extern double min_desc_grad_mod_3d(ptr_minimisation_param min_parameters);
extern double min_quasi_newton_mod_3d(ptr_minimisation_param min_parameters);
extern double min_flux_local_3d(ptr_minimisation_param min_parameters);
extern double min_powel2_3d(ptr_minimisation_param min_parameters);
extern double min_simplex_robuste_3d(ptr_minimisation_param min_parameters);

/************************ Powel ToolBox ***************************************/
extern double PowellOptimization (ptr_minimisation_param min_parameters, double **dTabParam,double **dTabDirection);
extern void PowellVectorProduct (double ** dVectSrc, int nNbVect, double * dVectResult);
/************************ toolbox pour le simplexe *********************/
void get_simplexe_psum(double **simplexe, double *psum, int Nest, int nb_param);
double reflexion_simplex_3d(ptr_minimisation_param min_parameters,
                        double ** simplexe, double *energies,
                        double *psum, int ihi, int Nest, double fac);
/************************ toolbox pour le simplexe robuste *********************/
void compute_robust_weights_quad(ptr_minimisation_param min_parameters);
void compute_robust_weights_woods(ptr_minimisation_param min_parameters);

extern void tri_double_index(double *arr, int taille, int *index);

extern int imx_query_minimisation_type_3D(int dist_type);
extern min_func_t imx_choose_minimisation_fct(int min_type,int dist_type);
bool arret_selon_parametres(double *new_param, double *old_param, double *prec_param, int nb_param);
bool arret_selon_energie(double E1, double E2);

#endif /*__MINIMISATION_3D_H__*/
