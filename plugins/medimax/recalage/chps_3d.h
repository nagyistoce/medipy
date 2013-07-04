/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/


#ifndef chps_3d_h
#define chps_3d_h
/**-----------------------------------------------------------------------------
*** 
*** file:       chps_3d.h
***
*** project:    Imagix 1.01
***         
***
*** description:    Fichier de gestion et de traitement des champs 3D
***  
*** 
*** Copyright (c) 1993, ULP-IPB Strasbourg.
*** All rights are reserved.
***
***
***---------------------------------------------------------------------------*/


#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/

#include "recalage/definitions_types_recalage.h"

/****************PROTOTYPAGE DES FONCTIONS DU FICHIER CHPS_3D.C***************/

/***************** ALLOCATION ET LIBERATION MEMOIRE ***************************/
extern field3d *cr_field3d(int wdth, int hght, int dpth);
extern int  free_field3d(field3d *ch);
extern transf3d *cr_transf3d(int wdth, int hght, int dpth, char *t_ini);
extern transf3d *cr_transf3d_p(grphic3d *imref, char *t_ini);
extern transf3d *cr_transf_anisotrope_3d(transf3d *transfOriginal,grphic3d *imref, grphic3d *imreca);
extern transf3d *cr_copy_transf3d(transf3d *transfo_src);
extern int realloc_transf_param_3d(transf3d *transf,int nb_param);
extern int copie_transf_3d(transf3d *src, transf3d *dst);
extern int free_transf3d(transf3d *tr);
/************** OPERATION MATHEMATIQUES SUR LES CHAMPS ************************/
extern void     oper_field_3d(void);
extern int  add_field_3d(field3d *champ1, field3d *champ2, field3d *champres);
extern int  sub_field_3d(field3d *champ1, field3d *champ2, field3d *champres);
extern int  div_field_3d(field3d *champ1, field3d *champres, float coeff);
extern int  mul_field_3d(field3d *champ1, field3d *champres, float coeff);
extern int  square_field_3d(field3d *champ1, field3d *champres);
extern int  rootsquare_field_3d(field3d *champ1, field3d *champres);

//---combinaison de transformations---//
extern int decomposer_rotation(double ** matrice_rot, double *angles_rotation);
extern int composer_rotation (double *rotation1, double *rotation2, double *rotation_res);
extern int comb_transf_3d (transf3d *transfo1,transf3d *transfo2,transf3d *transfores, int interpolation);

//---inversion de transformation---//
extern void inverser_transformation_3d();
extern int imx_inverser_transformation_3d(char *nom_fich_transfo_deb, char *nom_fich_transfo_res);
extern int imx_inverser_transformation_3d_p(transf3d *transfo_deb, transf3d *transfo_res);
int imx_inverser_transformation_rigide_3d_p(transf3d *transfo_deb, transf3d *transfo_res);
int imx_inverser_transformation_affine_3d_p(transf3d *transfo_deb, transf3d *transfo_res);
int imx_inverser_transformation_rigidezoom_3d_p(transf3d *transfo_deb, transf3d *transfo_res);

/************************* CHAMPS -> IMAGES ***********************************/
int module_field_to_image_3d(field3d *ch, grphic3d *im);
int composante_field_to_image_3d(field3d *ch, grphic3d *im,int type);
int images_to_field_3d(grphic3d *ux,grphic3d *uy,grphic3d *uz,field3d *ch);

/********************  GRADIENT D'UNE IMAGE   *********************************/
int imx_gradient_3d_p(grphic3d *im, field3d *grad, int method, int t);
/*********************** CARACTERISTIQUES CHAMP *******************************/
double  erreur_field_3d(field3d *champ1, field3d *champ2); 
int carac_field_3d(field3d *champ);

/********************* METHODES D'INTERPOLATIONS ******************************/
extern int  inter_nearest_3d(grphic3d *imdeb, field3d *champ, grphic3d *imres);
extern int  inter_linear_3d(grphic3d *imdeb, field3d *champ, grphic3d *imres);
extern int  inter_sinc_3d(grphic3d *imdeb, field3d *champ, grphic3d *imres);
extern int  inter_qsinc2_3d(grphic3d *imdeb, field3d *champ, grphic3d *imres);
extern int  inter_qsinc3_3d(grphic3d *imdeb, field3d *champ, grphic3d *imres);
int inter_qsinc_3d(grphic3d *imdeb, field3d *champ, grphic3d *imres, int tf);
extern int  inter_Bspline_3d_2(grphic3d *imdeb, field3d *champ, grphic3d *imres);
extern int  inter_Bspline_3d_3(grphic3d *imdeb, field3d *champ, grphic3d *imres);
extern int  inter_Bspline_3d_4(grphic3d *imdeb, field3d *champ, grphic3d *imres);
extern int  inter_Bspline_3d_5(grphic3d *imdeb, field3d *champ, grphic3d *imres);
int inter_Bspline_3d(grphic3d *imdeb, field3d *champ, grphic3d *imres,int degre);
extern void inter_linear_chps_3d(field3d *champ1,double x,double y,double z,vector3d *p);
extern void Inter_Chps_Bspline_3d(field3d *ch1, field3d *ch2,field3d *chres,int degre,double dx1, double dy1, double dz1, double dx2, double dy2, double dz2);
extern int  inter_labeled_3d(grphic3d *imdeb, field3d *champ, grphic3d *imres);
extern int  inter_VP_3d(grphic3d *imdeb, field3d *champ, grphic3d *imres);

//---toolbox pour l'interpolation BSpline---//
double Bspline(int degre,double valeur);
void GetInterpolationCoefficients (double c[], int DataLength, double z[], int NbPoles, double Tolerance);
double InitialCausalCoefficient (double c[], int DataLength, double z, double tolerance);
double InitialAntiCausalCoefficient (double c[], int DataLength, double z);

extern int normaliser_3d(grphic3d *imdeb, grphic3d *imres, double ***tmp_res, int wdth, int hght, int dpth);
extern int normaliser_float_3d(grphic3d *imdeb, grphic3d *imres, float ***tmp_res, int wdth, int hght, int dpth);
/************************ Interpolation **************************************/
extern InterpolationFct imx_choose_interpolation_fct(int inter_type);
extern int imx_query_interpolation_type_3D(int dist_type);
/************* ENREGISTREMENT ET LECTURE SUR FICHIER **************************/
int test_fichier_transf_3d(char *nomfichier);
char    *quest_save_result_3d(int *save_type);
int     save_transf_3d(transf3d *transfo, char *nomfichier);
transf3d *load_transf_3d(char *nomfichier);
/****************************** DIVERS  ***************************************/
//---applications de transformation---//
extern void apply_transf_3d(void);
extern int imx_apply_transf_3d(int im_res, int im_deb, char *nomfichier, int inter_type);
extern int imx_apply_transf_3d_p(grphic3d *imres, grphic3d *imdeb, transf3d *transfo, InterpolationFct inter_fct);

extern void apply_transf_seuil_3d(void);
extern int imx_apply_transf_seuil_3d(int im_res, int im_deb, char *nomfichier, int inter_type);

//---visualisation d'une transformation---//
extern void     visu_transf_3d(void);  
int jacobien_transf3d_to_image_3d (transf3d *transfo, grphic3d *imres);
int imx_visu_field_3d(field3d *ch, grphic3d *imres, int type);
int jacobien_champ_interpolation_bspline1(transf3d *ch, grphic3d *imres);
int integrale_jacobien_bspline1(transf3d *transfo_src, grphic3d *imres);

/************************ DEFORMATION DE SYNTHESE *****************************/
extern int imx_defor_synth_3d_point(char *nomfichres, int wdth, int hght, int dpth,int,vector3d *Pts,vector3d *Dep);
extern void defor_synth_3d_point(void);
extern void  defor_synth_3d(void);
extern int imx_defor_synth_3d(char *nomfichres, int wdth, int hght, int dpth);
/************************* INTERPOLATION D'UN CHAMP ****************************
*********************A PARTIR DE VECTEURS SPECIFIE ****************************/
extern int inter_field_thinplate_3d(field3d *ch, vector3d *points, vector3d *vecteurs, int nb_points);


extern void  compare_field_3d(void) ;
extern int imx_compare_field_3d_p(field3d *ch1, field3d *ch2, grphic3d *immod, grphic3d *imphase) ;
extern int jacobien_field_to_image_3d(field3d *ch, grphic3d *im) ;
extern field3d *transf_to_field_3d(transf3d *transfo,grphic3d *imref, grphic3d *imtransfo) ;
extern transf3d *field_to_transf_3d(field3d *champ,grphic3d *imref, grphic3d *imtransfo) ;
extern void field_to_transf_3d_noalloc(field3d *champ,transf3d* transfo,grphic3d *imref, grphic3d *imreca);

extern int resize_transf_3d(transf3d *transfo, int wdth, int hght, int dpth);
extern int transf_to_field_3d_noalloc(transf3d *transfo, field3d *champ,grphic3d *imref, grphic3d *imtransf);
extern int transf_to_field_champ_3d_noalloc(transf3d *transfo, field3d* points_choisis, field3d *champ,grphic3d *imref, grphic3d *imreca);

extern void apply_inverse_primitives_transf_3d(void);
extern int imx_apply_inverse_primitives_transf_3d(int im_res, int im_deb, char *nomfichier);
extern int imx_apply_inverse_primitives_transf_3d_p(grphic3d *imres, grphic3d *imdeb, transf3d *transfo);
extern void convert_trf_to_chp(void);
extern void subsample_chp(void);

extern void warpDeformationField_menu(void);
extern void warpDeformationField_p(field3d* chRef, transf3d* transfo, field3d* chRes, int inter_type);
extern void eval_matrice_jacobienne_3d(field3d *ch, int i, int j, int k, double** J);
extern void reorientDeformationField_p(field3d* chRef, transf3d* transfo, field3d* chRes);

#ifdef __cplusplus
}
#endif /*__cplusplus*/


#endif /* chps_3d_h */
 
