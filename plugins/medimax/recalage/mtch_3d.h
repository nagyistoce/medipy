/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef MTCH_3D_H
#define MTCH_3D_H

/**-----------------------------------------------------------------------
***	
***	file:		mtch_3d
***
***	project:	Imagix 1.01
***			
***
***	description:    Fichier pour le matching des images 3D
***	
***	
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***
***---------------------------------------------------------------------*/


#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/


#include "recalage/definitions_types_recalage.h"

/******* VARIABLES GLOBALES*********/
extern base3d	BASE3D;	/* effectivement utilise dans chps_3d.c */

/****************PROTOTYPAGE DES FONCTIONS DU FICHIER MTCH_3D.C***************/

/************************* TRANSFORMATIONS  ***********************************/
extern int rigid_to_field_3d(int nb_param,  const double *param, field3d *champ,grphic3d *imref, grphic3d *imreca);
extern int	rigidz_to_field_3d(int nb_param, const double *param, field3d *champ,grphic3d *imref, grphic3d *imreca);
extern int	affine_to_field_3d(int nb_param, const double *param, field3d *champ,grphic3d *imref, grphic3d *imreca);
extern int	rigid_global_zoom_to_field_3d(int nb_param, const double *param, field3d *champ,grphic3d *imref, grphic3d *imreca);
extern int	affine_decouple_to_field_3d(int nb_param, const double *param, field3d *champ,grphic3d *imref, grphic3d *imreca);
extern int  champ_to_field_3d(transf3d *transfo,field3d *champ, grphic3d *imref,grphic3d *imreca);
//extern int	affinesscg_to_field_3d(int nb_param, const double *param, field3d *champ,grphic3d *imref, grphic3d *imreca);
extern int	base_to_field_3d(int nb_param,   const double *param, field3d *champ,grphic3d *imref, grphic3d *imreca);
/*************** PASSAGE TRANSFORMATIONS **************************************/
extern int	gravite_to_rigid_3d(double *param);
extern int	rigid_to_rigidz_3d(double *param);
extern int	rigid_to_rigid_global_zoom_3d(double *param);
extern int	rigid_global_zoom_to_rigidz_3d(double *param);
extern int	rigidz_to_affine_decouple_3d(double *param);
extern int	rigidz_to_affine_3d(double *param);
extern int	affine_decouple_to_affine_3d(double *param);
extern int	affine_to_affinesscg_3d(double *param);

/**************** GESTION DE LA BASE DE FONCTIONS *****************************/
/*extern int	Bsplined0(int pas, double **f, double **df);
extern int	Bsplined1(int pas, double **f, double **df);
extern int	Bsplined2(int pas, double **f, double **df);
extern int	Bsplineod1(int pas, double **f, double **df);*/
extern int	init_base_3d(int wdth, int hght, int dpth, scal_func_t scal_func);
extern int	init_base_resol_3d(int wdth, int hght, int dpth, scal_func_t scal_func, int resolution);
extern int	base_resol_up_3d(double *param, int nb_param);
extern int	reduc_base_3d(double *param, int nb_param, grphic3d *imref, grphic3d *imreca, grphic3d *imres);
extern int	end_base_3d(void);
/***************** ALIGNEMENT DES CENTRES DE GRAVITE **************************/
extern double 	*gravite_to_field_3d(grphic3d *imref, grphic3d *imreca, field3d *champ);
extern TYPEMRI3D imx_find_cg_with_seuil_3d_p(grphic3d *img,double *cgX,double *cgY,double *cgZ);
extern int       imx_find_cg_geometric_3d_p(grphic3d *img,double *cgX,double *cgY,double *cgZ);

//---toolbox pour les demande de matching---//
void query_simple_matching_parameters_3d( int (*matching_3d_func)(int, int, int, int, int, int, int, char *, char *,eResearchInterv, eMatchPrecision));
void query_simple_precision_matching_parameters_3d( int (*matching_3d_func)(int, int, int, int, int, int, int, char *, char *,eResearchInterv, eMatchPrecision));
void query_multires_matching_parameters_3d( int (*matching_3d_func)(int, int, int, int, int, int, int, char *, char *, eMatchPrecision, int, int));
void query_multires_precision_matching_parameters_3d( int (*matching_3d_func)(int, int, int, int, int, int, int, char *, char *, eMatchPrecision, int, int));
/************************ MATCHING  RIGIDE ************************************/
extern void 	matching_rigide_3d();
extern void 	matching_rigide_precision_3d();
extern int 	imx_matching_rigide_3d(int im_ref, int im_reca, int im_res, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres, char *nomfichtrf,eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision);
extern int 	imx_matching_rigide_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int dist_type, int inter_type, int min_type, int save_type,transf3d *transfres, transf3d *initrf,eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision);
/*-----------------MATCHING  RIGIDE MULTIRESOLUTION------------------------------*/
extern void 	matching_rigide_multi_3d();
extern void 	matching_rigide_multi_precision_3d();
extern int  imx_matching_rigide_multi_3d(int im_ref, int im_reca, int im_res, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres, char *nomfichtrf, eMatchPrecision matchPrecision, int start_resol, int end_resol);
extern int 	imx_matching_rigide_multi_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int dist_type, int inter_type, int min_type, int save_type,transf3d *transfres, transf3d *initrf, eMatchPrecision matchPrecision, int start_resol, int end_resol);
/************************ MATCHING  RIGIDE+ZOOM ************************************/
extern void 	matching_rigide_zoom_3d();
extern void 	matching_rigide_zoom_precision_3d();
extern int 	imx_matching_rigide_zoom_3d(int im_ref, int im_reca, int im_res, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres, char *nomfichtrf, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision);
extern int 	imx_matching_rigide_zoom_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int dist_type, int inter_type, int min_type, int save_type, transf3d *transfres, transf3d *inittrf,eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision);
/*-----------------MATCHING  RIGIDE+ZOOM MULTIRESOLUTION------------------------------*/
extern void 	matching_rigide_zoom_multi_3d();
extern void 	matching_rigide_zoom_multi_precision_3d();
extern int  imx_matching_rigide_zoom_multi_3d(int im_ref, int im_reca, int im_res, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres, char *nomfichtrf, eMatchPrecision matchPrecision, int start_resol,int end_resol);
extern int 	imx_matching_rigide_zoom_multi_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int dist_type, int inter_type, int min_type, int save_type,transf3d *transfres, transf3d *initrf, eMatchPrecision matchPrecision, int start_resol,int end_resol);
/************************ MATCHING  AFFINE ************************************/
extern void 	matching_affine_3d();
extern void 	matching_affine_precision_3d();
extern int 	imx_matching_affine_3d(int im_ref, int im_reca, int im_res, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres, char *nomfichtrf, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision);
extern int 	imx_matching_affine_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int dist_type, int inter_type, int min_type, int save_type, transf3d *transfres, transf3d *inittrf, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision);
/************************ MATCHING  AFFINE MULTIRESOLUTION************************************/
extern void 	matching_affine_multi_3d();
extern void 	matching_affine_multi_precision_3d();
extern int 	imx_matching_affine_multi_3d(int im_ref, int im_reca, int im_res, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres, char *nomfichtrf, eMatchPrecision matchPrecision, int start_resol,int end_resol);
extern int 	imx_matching_affine_multi_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int dist_type, int inter_type, int min_type, int save_type,transf3d *transfres, transf3d *initrf, eMatchPrecision matchPrecision, int start_resol, int end_resol);
/************************ MATCHING  RIGIDE MULTISTART MULTIRESOLUTION************************************/
extern void 	matching_rigide_multistart_multires_3d();
extern void 	matching_rigide_multistart_multires_precision_3d();
extern int 	imx_matching_rigide_multistart_multires_3d(int im_ref, int im_reca, int im_res, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres, char *nomfichtrf, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision);
extern int 	imx_matching_rigide_multistart_multires_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int dist_type, int inter_type, int min_type, int save_type,transf3d *transfres, transf3d *initrf,eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision);
/************************ MATCHING  RIGIDE ZOOM MULTISTART MULTIRESOLUTION************************************/
extern void 	matching_rigide_zoom_multistart_multires_3d();
extern void 	matching_rigide_zoom_multistart_multires_precision_3d();
extern int 	imx_matching_rigide_zoom_multistart_multires_3d(int im_ref, int im_reca, int im_res, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres, char *nomfichtrf, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision);
extern int 	imx_matching_rigide_zoom_multistart_multires_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int dist_type, int inter_type, int min_type, int save_type,transf3d *transfres, transf3d *initrf,eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision);
/************************ MATCHING  AFFINE MULTISTART MULTIRESOLUTION************************************/
extern void 	matching_affine_multistart_multires_3d();
extern void 	matching_affine_multistart_multires_precision_3d();
extern int 	imx_matching_affine_multistart_multires_3d(int im_ref, int im_reca, int im_res, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres, char *nomfichtrf, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision);
extern int 	imx_matching_affine_multistart_multires_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int dist_type, int inter_type, int min_type, int save_type,transf3d *transfres, transf3d *initrf,eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision);
/************************ MATCHING  BSPLINE ***********************************/
extern void 	matching_Bspline_3d();
extern int 	imx_matching_Bspline_3d(int im_ref, int im_reca, int im_res, int func_type, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres, int resolf, char *nomfichiertrf, int renormalisation);
extern int 	imx_matching_Bspline_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int func_type, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres, int resolf, char *nomfichiertrf, int renormalisation);
//extern int imx_matching_Bspline_pyramide_3d_p(grphic3d **pyr_imref, grphic3d **pyr_imreca, int nb_niv_pyr, grphic3d *imres, int func_type, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres, int resolf, char *nomfichiertrf);
/************************** MATCHING  FINE ************************************/
extern void matching_fine_3d(/*int renormalisation MODIF. BERST 15/10/2002*/);
extern int imx_matching_fine_3d(int im_ref, int im_reca, int im_res, int dist_type, int inter_type, char *nomfichres, int renormalisation);
extern int imx_matching_fine_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int dist_type, int inter_type, char *nomfichres, int renormalisation);
/************************** MATCHING  MANUEL ************************************/
extern void matching_manuel_3d(void);
extern int imx_matching_manuel_3d(int im_ref, int im_reca, int im_res, int inter_type, char *nomfichres);
extern int imx_matching_manuel_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int inter_type, char *nomfichres);
/************************ MATCHING  ICP ************************************/
extern void 	matching_rigide_ICP_3d();
extern int imx_matching_rigide_ICP_3d(int im_ref, int im_reca, int im_res, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres, char *nomfichtrf,eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision);
extern int imx_matching_rigide_ICP_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int dist_type, int inter_type, int min_type, int save_type, transf3d *transfres, transf3d *inittrf,eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision);

extern void 	matching_rigide_zoom_ICP_3d();
extern int imx_matching_rigide_zoom_ICP_3d(int im_ref, int im_reca, int im_res, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres, char *nomfichtrf,eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision);
extern int imx_matching_rigide_zoom_ICP_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int dist_type, int inter_type, int min_type, int save_type, transf3d *transfres, transf3d *inittrf,eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision);

extern void 	matching_affine_ICP_3d();
extern int imx_matching_affine_ICP_3d(int im_ref, int im_reca, int im_res, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres, char *nomfichtrf,eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision);
extern int imx_matching_affine_ICP_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int dist_type, int inter_type, int min_type, int save_type, transf3d *transfres, transf3d *inittrf,eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision);


//---creation d'une image isotrope---//
extern void creer_image_isotrope_3d();
extern int imx_creer_image_isotrope_3d_p(grphic3d *im, InterpolationFct interpol, grphic3d *imres);
 
//------INITIALISATION DU RECALAGE-----//

int imx_matching_preprocessing_3d_p(grphic3d *imref,grphic3d *imreca, grphic3d **imtref, grphic3d **imtreca, grphic3d **imtres, transf3d **preProcessingRefTrf, transf3d **preProcessingRecaTrf, dist_func_t  distance);
transf3d *imx_trouver_transf_initiale_matching_3d_p(grphic3d *imref, grphic3d *imreca, transf3d *inittrf);
//---recadrage du cerveau---//
transf3d * imx_matching_recadrage_cerveau_3d_p(grphic3d *imdeb,grphic3d **imres);
void imx_matching_recadrage_cerveau_3d(int im_deb, int im_res);
void matching_recadrage_cerveau_3d();


//fonctions d'affichage des parametres en fonction du type de transformation
void imx_aff_param(TRANSFO3D transfo, double *param);
void imx_aff_param_rigid(double *param);
void imx_aff_param_rigid_zoom(double *param);
void imx_aff_param_affine(double *param);
void imx_aff_param_rigid_global_zoom(double *param);
void imx_aff_param_affine_decouple(double *param);

//fonctions d'init des bornes des param en fonction du type de transformation
double get_facteur_precision(eMatchPrecision matchPrecision);
double get_facteur_intervalle_recherche(eResearchInterv FieldOfResearch);
void init_bornes_matching_centre(double *min_param, double *max_param, double *prec_param, double *param, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision);
void init_bornes_matching_rot(double *min_param, double *max_param, double *prec_param, double *param,  eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision, double dx, double dy, double dz);
void init_bornes_matching_trans(double *min_param, double *max_param, double *prec_param, double *param, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision, double dx, double dy, double dz);
void init_bornes_matching_zoom(double *min_param, double *max_param, double *prec_param, double *param, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision, double dx, double dy, double dz);
void init_bornes_matching_affine(double *min_param, double *max_param, double *prec_param, double *param, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision, double dx, double dy, double dz);
void init_bornes_matching_skew(double *min_param, double *max_param, double *prec_param, double *param, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision, double dx, double dy, double dz);
void init_bornes_matching_global_zoom(double *min_param, double *max_param, double *prec_param, double *param, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision, double dx, double dy, double dz);
void init_bornes_matching(TRANSFO3D transfo, double *min_param, double *max_param, double *prec_param, double *param, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision, grphic3d *imref, grphic3d *imreca);

extern int rigid_to_field_champ_3d(int nb_param,  const double *param, field3d* points_choisis, field3d *champ,grphic3d *imref, grphic3d *imreca);
extern int rigidz_to_field_champ_3d(int nb_param,  const double *param, field3d* points_choisis, field3d *champ,grphic3d *imref, grphic3d *imreca);
extern int affine_to_field_champ_3d(int nb_param,  const double *param, field3d* points_choisis, field3d *champ,grphic3d *imref, grphic3d *imreca);
extern int rigid_global_zoom_to_field_champ_3d(int nb_param,  const double *param, field3d* points_choisis, field3d *champ,grphic3d *imref, grphic3d *imreca);
extern int affine_decouple_to_field_champ_3d(int nb_param,  const double *param, field3d* points_choisis, field3d *champ,grphic3d *imref, grphic3d *imreca);

//sous fonctions utilisees dans les methode de matching
int imx_calc_image_resultat_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, field3d *champ, transf3d *transfres);
int imx_calc_transfo_sav_image_res_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, transf3d *transfres, int save_type);
int imx_aff_temps_calcul(int tdebut, int tfin);

//fonctions misc ayant rapport aux transformations
extern transf_func_t imx_query_transformation_fct(TRANSFO3D transfo);
bool transfo_type_est_parametree(TRANSFO3D transfo);
int imx_copie_dim_param_to_transf3d_p(grphic3d *imsrc, transf3d *transfo);

//recalage de base sans pretraitement et qui minimise selon le meme type de transformation que inittrf
double imx_matching_simple_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, dist_func_t dist, InterpolationFct inter, min_func_t minimisation, transf3d *transfres, transf3d *inittrf, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision, contrainte_func_t contrainte, field3d *champ_ini, field3d *champ_fin);
//recalage en multiresolution de base sans pretraitement et qui minimise selon le meme type de transformation que inittrf
double imx_matching_multires_simple_3d_p(grphic3d *imref, grphic3d *imreca, dist_func_t dist, InterpolationFct inter, min_func_t minimisation, transf3d *transfres, transf3d *inittrf, eMatchPrecision matchPrecision, contrainte_func_t contrainte, field3d *champ_ini, int start_resol,int end_resol);

//--------------- toolbox pour le multiresoltion multistart ----------------------/
int determination_minima_grille_angulaire(double *energies, int *indices, int nb_angles);
int evaluer_parametres_grille_fine(double **param_calc_fin, double **param_calc_coarse, int nb_angles_fin, int nb_angles_grossier, TRANSFO3D transfoType);
int initialiser_parametres_grossier(double **param_calc_coarse, int nb_angles_grossier, double *max_param, double *min_param, double *param, TRANSFO3D transfoType);
double imx_matching_multistart_multires_8mm_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, dist_func_t dist, InterpolationFct inter, min_func_t minimisation, transf3d *transfres, transf3d *inittrf, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision, contrainte_func_t contrainte, field3d *champ_ini, field3d *champ_fin, int *nb_minima, double **minima);
double imx_matching_multistart_multires_4mm_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, dist_func_t dist, InterpolationFct inter, min_func_t minimisation, transf3d *transfres, transf3d *inittrf, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision, contrainte_func_t contrainte, field3d *champ_ini, field3d *champ_fin, int nb_minima, double **minima);

//fonctions ayant trait a la demande de precision
extern eMatchPrecision imx_query_matching_precision();

extern ptr_distance_param cr_distParam_from_minParam(ptr_minimisation_param min_parameters, min_func_t minimisation);
extern int free_dist_param(ptr_distance_param dist_param);




typedef int (*matching_func_2d)(grphic *imref, grphic *imreca, grphic *imres, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres);


#ifdef __cplusplus
}
#endif /*__cplusplus*/

#endif

