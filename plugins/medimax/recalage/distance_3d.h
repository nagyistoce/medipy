/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*******************************************************************************/	
/*!	   \file:		distance_3d.h
***
***		project:	Imagix 2.01
***			
***
***		\brief description:    Calcul de la distance interimage
***	
***	
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***
********************************************************************************/


#ifndef __DISTANCE_3D_H__
#define __DISTANCE_3D_H__

#include "recalage/definitions_types_recalage.h"

/****************************** DISTANCES *************************************/
extern double   erreur_quad_3d(ptr_distance_param dist_param);
extern double   erreur_quad_robust_geman_3d(ptr_distance_param dist_param);
extern double   old_erreur_woods_3d(ptr_distance_param dist_param);
extern double   erreur_woods_3d(ptr_distance_param dist_param);
extern double   erreur_woods_robust_3d(ptr_distance_param dist_param);
extern double   erreur_IM_3d(ptr_distance_param dist_param);
extern double   erreur_entropie_conjointe_simple_3d(ptr_distance_param dist_param);
extern double   erreur_IM_simple_3d(ptr_distance_param dist_param);
extern double   erreur_IM_Normalisee_Studholme_simple_3d(ptr_distance_param dist_param);
extern double   erreur_IM_Normalisee_Maes1_simple_3d(ptr_distance_param dist_param);
extern double   erreur_IM_Normalisee_Maes2_simple_3d(ptr_distance_param dist_param);
extern double   erreur_CR_3d(ptr_distance_param dist_param);
extern double   erreur_CR_robust_3d(ptr_distance_param dist_param);
extern double   erreur_Woods_vraie_3d(ptr_distance_param dist_param);
extern double   erreur_woods_robust2_3d(ptr_distance_param dist_param);
extern double   erreur_entropie_conjointe_VP_3d(ptr_distance_param dist_param);
extern double   erreur_IM_VP_3d(ptr_distance_param dist_param);
extern double   erreur_IMNS_VP_3d(ptr_distance_param dist_param);
extern double   erreur_IMNM1_VP_3d(ptr_distance_param dist_param);
extern double   erreur_IMNM2_VP_3d(ptr_distance_param dist_param);
extern double   erreur_quad_robust2_3d(ptr_distance_param dist_param);
extern double   erreur_IM_VP_stochastique_3d(ptr_distance_param dist_param);
extern double   erreur_IM_BSpline_3d(ptr_distance_param dist_param);
extern double   erreur_IM_Region_3d(ptr_distance_param dist_param);
extern double   erreur_coeff_correlation_d(ptr_distance_param dist_param);
extern double	  erreur_ICP_3d(ptr_distance_param dist_param);
extern double	  erreur_ICP_sym_3d(ptr_distance_param dist_param);

/****fonctions auxilaires pour le calculs bases sur l'information mutuelle simple****/
double calcul_info_mutuelle_3d(IM3DTYPE im_type, ptr_distance_param dist_param);
extern double calcul_IM(double entropie_conjointe, double entropie1, double entropie2);
extern double calcul_IMNS(double entropie_conjointe, double entropie1, double entropie2);
extern double calcul_IMNM1(double entropie_conjointe, double entropie1, double entropie2);
extern double calcul_IMNM2(double entropie_conjointe, double entropie1, double entropie2);

/****fonctions auxilaires pour le calculs bases sur l'interpolation volume partiel****/
double calcul_info_mutuelle_VP_3d(IM3DTYPE im_type, ptr_distance_param dist_param);
extern int calcul_histogrammes_VP(ptr_distance_param dist_param);
extern int calcul_histogrammes_gradient_VP(grphic3d *Imreca, grphic3d *Imref, field3d *champ, VP_histos * donnees_VP, ERREUR_RECALAGE *err);
void update_histogrammes_VP(int pixelRef, VP_histos * donnees_VP, TYPEMRI3D *v, double *w, double xrel, double yrel, double zrel);
void update_histogrammes_gradient_VP(int pixelRef, VP_histos * donnees_VP, TYPEMRI3D *v, double dw[8][3], int *p, double xrel, double yrel, double zrel);
extern bool dist_utilisant_VP(dist_func_t distance);

/***fonctions auxilaires pour le calculs bases sur le calcul stochastique***/
double calcul_info_mutuelle_VP_stochastique_3d(IM3DTYPE im_type, ptr_distance_param dist_param);
extern field3d *choix_points_aleatoires(int width, int height, int depth, ERREUR_RECALAGE *err);
extern int calcul_histogrammes_VP_stochastique(ptr_distance_param dist_param);
extern void update_histogrammes_VP_stochastique(VP_histos * donnees_VP, TYPEMRI3D *vRef, TYPEMRI3D *vReca, double *wRef, double *wReca, double xrelRef, double yrelRef, double zrelRef, double xrelReca, double yrelReca, double zrelReca);

/***fonctions auxilaires pour le calculs bases sur le calcul IM BSpline***/
double calcul_info_mutuelle_BSpline_3d(IM3DTYPE im_type, ptr_distance_param dist_param);
extern int calcul_histogramme_BSpline(grphic3d *Imreca, grphic3d *Imref, field3d *champ, double **histo_conjoint, int nbclReca, int nbclRef, int maxReca, int maxRef);

//---calcul de l'information par la methode regionelle---//
double calcul_info_mutuelle_region_3d(IM3DTYPE im_type, ptr_distance_param dist_param);
int calcul_covar_region_3d(ptr_distance_param dist_param, const int N_R, double **covar_mat);
int calcul_entropie_covar_region_3d(const int N_R, double **covar_mat, double *entropie_conjointe, double *entropie1, double *entropie2);

//---calcul de l'image edgeness---//
void calcul_edgeness_3d();
int im_calcul_edgeness_3d(int im_deb, int im_res, int taille_voisinage);
int im_calcul_edgeness_3d_p(grphic3d *imdeb, grphic3d *imres, int taille_voisinage);

/****fonctions generales utilisees par les distances se basant sur l'information mutuelle****/
extern void calcul_histo_conjoint_simple(grphic3d *im1, grphic3d *im2, double **histo_conjoint, int nbcl1, int nbcl2, int max1, int max2, ERREUR_RECALAGE *err);
extern void calcul_histo_conjoint_lineaire(grphic3d *im1, grphic3d *im2, double **histo_conjoint, int nbcl1, int nbcl2, int max1, int max2, ERREUR_RECALAGE *err);
extern void calcul_histos_marginaux(double **histo_conjoint, double *histo1, double *histo2, int nbcl1, int nbcl2, double nbVoxelsHorsOverlap);
extern double calcul_entropie_conjointe(double **histo_conjoint, int nbcl1, int nbcl2, double nbVoxelsOverlap);
extern double calcul_entropie(double *histo, int nbcl, double nbVoxelsOverlap);
double max_erreur_info_mutuelle(IM3DTYPE im_type, int nbcl1, int nbcl2);
double calc_info_mutuelle_histos(IM3DTYPE im_type, double **histo_conjoint, double *histo1, double *histo2, int nbcl1, int nbcl2, double nbVoxelsOverlap, double nbVoxelsHorsOverlap, double *EntropieConjointe, double *Entropie1, double *Entropie2);

extern int calcul_nb_classes(int wdth, int hght, int dpth, int max1, int max2, int *nbcl1, int *nbcl2);
extern int imx_calculer_nombre_de_classes_3d_p(grphic3d *image, int *max_image);

extern VP_histos * cr_VP_histos(grphic3d *Imreca, grphic3d *Imref,bool minimisation_par_gradient);
extern void free_VP_histos(VP_histos **donnees_VP_ptr);

extern int imx_query_erreur_type_3D();
extern dist_func_t imx_choose_distance_fct(int dist_type);

#endif /*__DISTANCE_3D_H__*/
