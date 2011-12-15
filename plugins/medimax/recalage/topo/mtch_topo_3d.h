/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef mtch_topo_3d_h
#define mtch_topo_3d_h
#include "recalage/mtch_3d.h"
#include "recalage/topo/analyse_intervalle.h"
#include "recalage/topo/normalisation_topo.h"
#include "recalage/topo/divers_topo.h"


#define TOPO_DEBUG      	/* On commente la ligne pour passer en mode Debug */

#define TOPO_VERBOSE 		/* On commente la ligne pour passer en mode Verbose */


#define SAVE_INTERMEDIAIRE  /* On commente la ligne pour enregistrer les transformations des resolutions intermediaires (version qui ne marche plus )*/

#define SAVE_INTERMEDIAIRE_2  /* On commente la ligne pour enregistrer les transformations des resolutions intermediaires (version qui marche) */


#define NB_PAS_LINEMIN_TOPO 6

#define TOPO_COMPARE_CRITERE /* on commente la ligne pour enregistrer plein de critere dans un fichier */

/*********** Parametre intervenant dans les fonctions de cout **********/
#define TOPO_EPSILON 0.01  /* parametre de la fonction de cout L1L2 qui definit le seuil a partir du quel la fonction de cout passe d'un comportement L2 a L1 */
#define ORDRE_NORME_LP 1.2  /* 1 < P < 2 */
#define ORDRE_REGULARISATION_LP 1.2  /* 1 < P < 2 */
#define TOPO_SIZE_HISTOJOINT 64  
#define TOPO_SIZE_NORM_HISTO 256  

#define SIZE_MIN_REGULARISATION 1000
/****** conversion d'un triplet (i,j,k) en l'indice correspondant dans le tableau param  ***********/
#define TOP_D(nb_param)  ( floor (pow(1.0*(nb_param)/3.0,1.0/3.0) + 0.1 ) )
#define TOP_conv_ind(i,j,k,nb_param) ( 3*( (i)*(TOP_D(nb_param))*(TOP_D(nb_param)) + (j)*(TOP_D(nb_param)) + (k)) )
#define TOP_conv_ind_topD(i,j,k,topD) ( 3*( (i)*(topD)*(topD) + (j)*(topD) + (k)) )

#define CRITERE_ARRET_PARAM 0.05 /* exprime en pourcentage de la largeur d'une fonction Bspline */


#define TOPO_QUANTIFICATION_PARAM 0.0001 /* defini le pas de quantification (en voxel) pour la recherche des parametres du modele de deformation Bspline */

/***** Param pour la regularisation avec patch ******/
//#define TOPO_PATCH_SIZE 1 // le patch est de taille 2*TOPO_PATCH_SIZE+1
//#define TOPO_PATCH_NEIGHBOR 5 // la zone de recherche est de taille 2*TOPO_PATCH_NEIGHBOR+1
#define TOPO_PATCH_BETA 1 // coefficient de regularite

extern int TOPO_PATCH_SIZE;
extern int TOPO_PATCH_NEIGHBOR;


/***** DEFINITION DU TYPE hessien3d *******/
typedef struct s_hessien3d { 
  int		width,height,depth; 			      	
  mat3d	***raw;
} hessien3d, *ptr_hessien3d;


typedef struct s_ParamRecalageBspline {
field3d *chref;
grphic3d** serie_groupwise;
grphic3d* reca_groupwise,*ref_groupwise,*ref2_groupwise;
double ***dbl_ref_groupwise,***dbl_ref2_groupwise;
field3d  *grad_reca_groupwise, *grad_ref_groupwise, *grad_ref2_groupwise;
hessien3d *hessien_reca_groupwise, *hessien_ref_groupwise, *hessien_ref2_groupwise; 
int nb_tot_groupwise;
int nb_reca_groupwise;
} ParamRecalageBspline, *ptr_ParamRecalageBspline;



typedef struct s_HistoJoint {
int nbI,nbJ; // taille de l'histogramme joint
int Ntot; //nombre d'�l�ments total
double maxI,maxJ; // valeurs maximales dans les images I (image reca) et J (image r�f�rence)
double normI,normJ; // coefficient de normalisation pour obtenir la valeur correspondante dans l'histojoint : normI = nbI/maxI
double **hij; // proba jointe
double *hi,*hj; // proba marginales
double IM;
double **aux_hij,*aux_hi,*aux_hj; // variables auxiliaires
double aux_IM;

} HistoJoint, *ptr_HistoJoint;

extern HistoJoint HISTOJOINT;

extern int LOG_ENERGY;  // variable globale pour identifier un fichier de log dans la fonction compare_critere  
extern double  TOPO_ALPHA_ROBUST; /* parametre de la fonction de cout robuste de Geman Mclure*/
extern double TOPO_REGULARISATION_SIGMA_GAUSSIEN;
extern double TOPO_REGULARISATION_MEMBRANE_ELASTIQUE;
extern double TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE;
extern double TOPO_SEUIL_SIGNIFICATIVITE;
extern melange_gaussienne_2d MELANGE_GAUSSIENNE_2D;
extern double TOPO_PONDERATION_RECALAGE_SYMETRIQUE; /* parametre de pnderation dans le recalage symetrique. Il doit etre compris entre 0 et 1 */
extern int PRECISION_SIMULATION;
extern ParamRecalageBspline _ParamRecalageBspline;

extern double PATCH_REF;
extern double PATCH_RECA;

extern double TOPO_PATCH_COEFF;
extern int TOPO_NLMhwnx;
extern int TOPO_NLMhwny;
extern int TOPO_NLMhwnz;
extern int TOPO_NLMhwvsx;
extern int TOPO_NLMhwvsy;
extern int TOPO_NLMhwvsz;

/****************** Definition des pointeurs de fonction  ****************/
typedef double 	(*reg_func_locale_t)	(int nb_param, double *param, double *param_norm, int topi, int topj, int topk, grphic3d *mask,TSlpqr *Slpqr);
typedef int 	(*reg_grad_locale_t)	(int nb_param, double *param, double *param_norm, double *grad, int topi, int topj, int topk, grphic3d *mask,TSlpqr *Slpqr);
typedef int 	(*reg_hessien_locale_t)	(int nb_param, double *param, double *param_norm, mat3d *hessien, int topi, int topj, int topk, grphic3d *mask,TSlpqr *Slpqr);
typedef double 	(*dist_func_locale_t)	(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, double *param_norm, int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation);
typedef double 	(*min_func_locale_t)	(grphic3d *imref, grphic3d *imreca, grphic3d *imres, field3d *champ_fin, transf_func_t transf, InterpolationFct inter, dist_func_locale_t dist,reg_func_locale_t regularisation, int nb_param, double *param, int ***masque_param,double Jmin, double Jmax,char *nomfichiertrf);
typedef void 		(*norm_func_locale_t)	(grphic3d *im1, grphic3d *im2, grphic3d *imres, int nb_classe);


/************************ MATCHING  BSPLINE ************************************/
extern void 	matching_Bspline_topo_3d();
extern void matching_Bspline_topo_norm_3d();

extern int 	imx_matching_Bspline_topo_3d(int im_ref, int im_reca, int im_res, int func_type,
																				 int dist_type, int regularisation, int inter_type, int min_type, int save_type, char *nomfichres,int resolf, 
																				 double Jmin, double Jmax,int normalisation_type, int nb_classe ,int adaptatif, int biais);

extern int 	imx_matching_Bspline_topo_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int func_type, 
																					int dist_type, int regularisation, int inter_type, int min_type, int save_type, char *nomfichres, int resolf,
																				  double Jmin, double Jmax, int normalisation_type, int nb_classe,int adaptatif, int biais);


scal_func_t topo_choose_bspline(int func_type);
dist_func_locale_t topo_choose_distance(int dist_type);
min_func_locale_t topo_choose_optimisation(int min_type);
norm_func_locale_t topo_choose_normalisation(int normalisation_type, int nb_classe);
reg_func_locale_t topo_choose_regularisation(int reg_type);
void topo_masque_param (grphic3d *imtref,grphic3d *imtreca,grphic3d *imtres,int ***masque_param, int D);
void annule_param_avec_masque_param (double *param, int nb_param, int ***masque_param);
void update_TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE(grphic3d *imtref,grphic3d *imtreca, int nb_param, double *param, double *param_norm, dist_func_locale_t distance,reg_func_locale_t regularisation, int ***masque_param);
void update_TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE2(grphic3d *imref, grphic3d *imreca, grphic3d *imres, field3d *champ_fin, min_func_locale_t minimisation,
                          					transf_func_t the_transf, InterpolationFct inter, dist_func_locale_t dist, reg_func_locale_t regularisation,
																	 	int nb_param,  double *param, int ***masque_param,double Jmin, double Jmax,char *nomfichiertrf);

void init_melange_gaussienne_2d(int nb_classe);
void free_melange_gaussienne_2d(void);

void init_HistoJoint(int nbI,int nbJ);
void free_HistoJoint(void);
void update_HistoJoint_global(grphic3d *imtref,grphic3d *imtres);

void init_IM_locale_before_optimization_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, int topi, int topj, int topk);
void init_IM_locale_after_optimization_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, int topi, int topj, int topk);

void imx_realloc_mri_pow2_centered_p (grphic3d *im);
void imx_undo_mri_pow2_centered_p (grphic3d *im, int wdth, int hght, int dpth);
void update_maxnorm_HistoJoint(grphic3d *imref,grphic3d *imreca);


extern void 	matching_Bspline_topo_primitive_3d();
extern int imx_matching_Bspline_topo_primitive_3d(int im_ref, int im_reca, int im_res, int func_type, int dist_type, int reg_type, int inter_type, int min_type, int save_type, char *nomfichres, int resolf, double Jmin, double Jmax);
int imx_matching_Bspline_topo_primitive_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int func_type, int dist_type, int reg_type, int inter_type, int min_type, int save_type, char *nomfichres, int resolf, double Jmin, double Jmax);
void topo_masque_param_primitive (grphic3d *imtref,int ***masque_param, int nb_param);


/************************ Simulation d'atrophie ************************************/
extern void 	simul_atrophie_topo_3d();
int imx_simul_atrophie_topo_3d(int im_ref, int im_mask,int func_type, int dist_type, int reg_type, int inter_type, int min_type, int save_type, char *nomfichres, int resolf,int adaptatif,double Jmin, double Jmax);
int imx_simul_atrophie_topo_3d_p(grphic3d *imref,grphic3d *mask, int func_type, int dist_type, int reg_type, int inter_type, int min_type, int save_type, char *nomfichres, int resolf,int adaptatif,double Jmin, double Jmax);





/************************ Recalage symatrique ************************************/

extern void matching_Bspline_topo_norm_symetrique_3d();

extern int 	imx_matching_Bspline_topo_symetrique_3d(int im_ref, int im_reca, int im_res, int func_type,
																				 int dist_type, int regularisation, int inter_type, int min_type, int save_type, char *nomfichres,int resolf, 
																				 double Jmin, double Jmax,int normalisation_type, int nb_classe ,int adaptatif);

extern int 	imx_matching_Bspline_topo_symetrique_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int func_type, 
																					int dist_type, int regularisation, int inter_type, int min_type, int save_type, char *nomfichres, int resolf,
																				  double Jmin, double Jmax, int normalisation_type, int nb_classe,int adaptatif);
/************************ Projection d'un champ pour qu'il preserve la topologie  ************************************/
void Project_Field_Bspline_topo_3d();
int imx_Project_Field_Bspline_topo_3d(char *nomfichier,char* nomfichres,int min_type,int resolf,double Jmin,double Jmax,int reg_type,double prec, int im_invariant);
void imx_Project_Field_Bspline_topo_3d_p(field3d* ch,transf3d *transfores,transf3d *transfo_opp,int min_type,int resolf,double Jmin,double Jmax,int reg_type, double prec,grphic3d *mask);
int imx_dti_matching_Bspline_topo_3d_p(int func_type, int dist_type, int reg_type, int inter_type, int min_type, int save_type, char *nomfichres, int resolf, double Jmin, double Jmax, int normalisation_type, int nb_classe);
extern void dti_matching_Bspline_topo_3d();
void groupwise_matching_Bspline_topo_3d();
int imx_groupwise_matching_Bspline_topo_3d(char* file, int nb_img, int num_ref, int func_type, int dist_type, int reg_type, int inter_type, int min_type, int save_type, char* nomfichres, int resolf, double Jmin, double Jmax, int normalisation_type, int nb_classe,int adaptatif,char *nomfichiertrf);
int imx_groupwise_matching_Bspline_topo_3d_p(int func_type,int dist_type, int reg_type, int inter_type,int min_type,int save_type, char* nomfichres,int resolf, double Jmin, double Jmax, int normalisation_type, int nb_classe,int adaptatif,char *nomfichiertrf);

/***** Recalage moment geometrique ******/
void matching_Bspline_topo_momentgeometrique_3d();
int imx_matching_Bspline_topo_momentgeometrique_3d(int im_ref, int im_reca, int im_res, int func_type, int dist_type, int reg_type, int inter_type, int min_type, int save_type, char *nomfichres, int resolf, double Jmin, double Jmax);
int imx_matching_Bspline_topo_momentgeometrique_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int func_type, int dist_type, int reg_type, int inter_type, int min_type, int save_type, char *nomfichres, int resolf, double Jmin, double Jmax);


#endif

