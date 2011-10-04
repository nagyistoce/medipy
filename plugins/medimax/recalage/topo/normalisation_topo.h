/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef normalisation_topo_h
#define normalisation_topo_h

#include "noyau/imx_2d.h"

#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/








/********* Definition de la structure Param_gaussienne_2d *******/

typedef struct s_Param_Gaussienne_2d { 
double mx,my, 			/*moyennes*/
			 sx,sy,sxy,		/*covariance*/
			 l;						/* facteur multiplicatif */ 
} Param_Gaussienne_2d, *ptr_Param_Gaussienne_2d;

/********* Definition de la structure Melange_gaussienne_2d *******/

typedef struct s_melange_gaussienne_2d { 
Param_Gaussienne_2d *param;
int nb_classe; 
grphic *histo;
} melange_gaussienne_2d, *ptr_melange_gaussienne_2d;



/******* Fonctions calculant l'histogramme joint par differentes methodes ********/

extern void 	histo_joint_3d								(void);
extern void 	histo_joint_3d_p							(grphic3d *im1, grphic3d *im2, grphic *implot);
extern void 	histo_joint_linear_3d					(void);
extern void 	histo_joint_linear_3d_p				(grphic3d *im1, grphic3d *im2, grphic *implot);
extern void 	histo_joint_linear_norm_3d		(void);
extern void histo_joint_linear_norm_3d_p(grphic3d *im1, grphic3d *im2, grphic *implot, double **confiance);
extern void 	histo_joint_maxentropie_3d_p	(grphic3d *im1, grphic3d *im2, grphic *implot);



/******* Normalisation par m�lange de gaussiennes sans calcul de l'histogramme joint ********/

extern void  		normalisation_gaussienne2d				(void);
extern void  		normalisation_gaussienne2d_p			(grphic3d *im1,grphic3d *im2,grphic3d *imres,int nb_classe);
extern void 		ini_kmeans_histo_joint_3d_p				(grphic3d *im1,grphic3d *im2, Param_Gaussienne_2d *param, int nb_classe);
extern void 		mean_std_classe										(grphic *im,grphic *imcl,int nb_classe, Param_Gaussienne_2d *param);
extern void 		fit_EM_gaussienne2d								(grphic3d *im1,grphic3d *im2, Param_Gaussienne_2d *param, int nb_classe);
extern void 		estime_EM_gaussienne2d						(grphic3d *im1,grphic3d *im2,Param_Gaussienne_2d *param,Param_Gaussienne_2d *new_param,int nb_classe);
extern double 	vraisemblance_gaussienne_2d				(grphic3d *im1,grphic3d *im2, Param_Gaussienne_2d *param, int nb_classe);
extern double 	eval_proba_apriori								(double x, double y, Param_Gaussienne_2d *param,int classe,int nb_classe);
extern void 		apply_normalisation_gaussienne2d	(grphic3d *imreca,grphic3d *imref,grphic3d *imres,int nb_classe, Param_Gaussienne_2d *param);



/******* Normalisation par m�lange de gaussiennes avec calcul de l'histogramme joint ********/
extern void  		normalisation_gaussienne2d_histo_joint								(void);
extern void  		imx_normalisation_gaussienne2d_histo_joint_standart_p	(grphic3d *im1,grphic3d *im2,grphic3d *imres,int nb_classe);
extern void  		normalisation_gaussienne2d_histo_joint_standart_p			(grphic3d *im1,grphic3d *im2,grphic3d *imres,grphic *histo,int nb_classe);
extern void  		imx_normalisation_gaussienne2d_histo_joint_norm_p	(grphic3d *im1,grphic3d *im2,grphic3d *imres,int nb_classe);
extern void  		normalisation_gaussienne2d_histo_joint_norm_p			(grphic3d *im1,grphic3d *im2,grphic3d *imres,grphic *histo,int nb_classe);
extern void  		imx_normalisation_gaussienne2d_histo_joint_reclasse_p	(grphic3d *im1,grphic3d *im2,grphic3d *imres,int nb_classe);
extern void  		normalisation_gaussienne2d_histo_joint_reclasse_p			(grphic3d *im1,grphic3d *im2,grphic3d *imres,grphic *histo,int nb_classe);


extern void 		fit_gaussienne2d_EM													(grphic *histo, Param_Gaussienne_2d *param, int nb_classe);
extern void 		estime_gaussienne2d_EM											(grphic *histo,Param_Gaussienne_2d *param,Param_Gaussienne_2d *new_param,int nb_classe);
extern double 	vraisemblance_gaussienne2d									(grphic *histo, Param_Gaussienne_2d *param, int nb_classe);
extern double 	eval_gaussienne_2d													(double x, double y, Param_Gaussienne_2d *param);
extern void 		apply_normalisation_gaussienne_histo_joint	(grphic3d *imreca,grphic3d *imref,grphic3d *imres,grphic *histo, int nb_classe, Param_Gaussienne_2d *param, double maxreca, double maxref, double alpha);
extern void 		apply_normalisation_gaussienne_histo_joint_multimodal	(grphic3d *imreca,grphic3d *imref,grphic3d *imres,grphic *histo, int nb_classe, Param_Gaussienne_2d *param);



/******* Normalisation par m�lange de gaussiennes (cas multimodal)  ********/

extern void  	normalisation_gaussienne2d_histo_joint_multimodal				(void);
extern void  	imx_normalisation_gaussienne2d_histo_joint_multimodal_p	(grphic3d *im1,grphic3d *im2,grphic3d *imres,int nb_classe);
extern void  	normalisation_gaussienne2d_histo_joint_multimodal_p			(grphic3d *im1,grphic3d *im2,grphic3d *imres,grphic *histo,int nb_classe);
extern void  	normalisation_gaussienne2d_histo_joint_multimodal2_p			(grphic3d *im1,grphic3d *im2,grphic3d *imres,grphic *histo,int nb_classe);
extern void 	apply_normalisation_gaussienne_histo_joint2							(grphic3d *imreca,grphic3d *imref,grphic3d *imres,grphic *histo, int nb_classe, Param_Gaussienne_2d *param);
extern void 	fit_gaussienne2d_EM_multi																(grphic *histo, Param_Gaussienne_2d *param, int nb_classe);
void init_segment_gaussienne2d(grphic3d *imreca,grphic3d *imref,grphic3d *imreca_cl,grphic3d *imref_cl, grphic *histo, int nb_classe, Param_Gaussienne_2d *param);
void raffine_segment_croissance_region(grphic3d *im,grphic3d *im_cl, double *moyenne, double *std  , int nb_classe);
void raffine_segment_croissance_region2(grphic3d *im,grphic3d *im_cl, double *moyenne, double *std  , int nb_classe);
void raffine_segment_croissance_region3(grphic3d *im,grphic3d *im_cl, double *moyenne, double *std  , int nb_classe);
void segment_mean_std_classe(grphic3d *im1,grphic3d *im1_cl,double *moyenne1,double *std1,int *nb1,int nb_classe);
void segment_mean_std_classe_robust(grphic3d *im,grphic3d *im_cl,double *moyenne,double *std,int nb_classe);
void apply_norm_mean_std_classe(grphic3d *im1,grphic3d *im1_cl,grphic3d *imres,double *moyenne1,double *moyenne2,double *std1,double *std2,double *covariance, Param_Gaussienne_2d *param,int
nb_classe);
void segment_covariance_classe(grphic3d *im1,grphic3d *im1_cl,grphic3d *im2,grphic3d *im2_cl,double *covariance,double *moyenne1, double
*moyenne2,double *std1, double *std2,int nb_classe);
void apply_norm_mean_std_classe_vp(grphic3d *im1,grphic3d *im1_cl,grphic3d *imres,double *moyenne1,double *moyenne2,double *std1,double *std2,double *covariance, Param_Gaussienne_2d *param,int
nb_classe);
void apply_norm_mean_std_classe_reg(grphic3d *im1,grphic3d *im1_cl,grphic3d *imres,double *moyenne1,double *moyenne2,double *std1,double *std2,double *covariance, Param_Gaussienne_2d *param,int
nb_classe);
void apply_norm_quantile_classe(grphic3d *im1,grphic3d *im1_cl,grphic3d *im2,grphic3d *im2_cl,grphic3d *imres,double *covariance,int nb_classe,int nb_quantile);
void raffine_segment_ML(grphic3d *im,grphic3d *im_cl, double *moyenne, double *std , int nb_classe );
double eval_mahalanobis_2d(double x, double y, Param_Gaussienne_2d *param);
int merge_classe(grphic3d *im1,grphic3d *im1_cl,grphic3d *im2,grphic3d *im2_cl,double *moyenne1,double *moyenne2,double *std1,double *std2,int *nb1,int *nb2,int nb_classe);
void raffine_segment_MAP(grphic3d *im,grphic3d *im_cl, double *moyenne, double *std , int *nb, int nb_classe );
void remplissage_carte_de_probabilite(grphic3d *im1,grphic3d *im1_cl,grphic3d *im2,grphic3d *im2_cl,grphic3d *imres,grphic *histo,double *moyenne1,double *moyenne2,double *std1,double *std2,int *nb1,int *nb2,int nb_classe);
void raffine_segment_potts(grphic3d *im,grphic3d *im_cl,double *moyenne,double *std,int nb_classe);
void raffine_segment_potts_phi_variable(grphic3d *im,grphic3d *im_cl,double *moyenne,double *std,int* nb,int nb_classe);
void calcul_fonction_partition_potts (void);
void swendsen_wang(grphic3d * im,double phi,int nb_classe);
int calcul_Uz (grphic3d *im,grphic3d *mask);
void calcul_fonction_partition_potts_p(double *Z, double incr, double phi_max ,int nb_classe);
void apply_norm_mean_std_classe_proba_reg(grphic3d *im1,grphic3d *im1_cl,grphic3d *imres,double *moyenne1,double *moyenne2,double *std1,double *std2,double *covariance, Param_Gaussienne_2d *param,int
nb_classe);

void  imx_normalisation_gaussienne2d_histo_joint_norm_symetrique_p(grphic3d *im1,grphic3d *im2, grphic3d *im1d,grphic3d *im2d, grphic3d *im1res,grphic3d *im2res,int nb_classe);
void  imx_normalisation_gaussienne2d_histo_joint_norm_symetrique_nonsymetrise_p(grphic3d *im1,grphic3d *im2, grphic3d *im1d,grphic3d *im2d, grphic3d *im1res,grphic3d *im2res,int nb_classe);
void  imx_normalisation_gaussienne2d_histo_joint_norm_symetrique_iteratif_p(grphic3d *im1,grphic3d *im2, grphic3d *im1d,grphic3d *im2d, grphic3d *im1res,grphic3d *im2res,int nb_classe);
void Update_image_from_histojointnorm_3d_p(grphic3d *im1, grphic3d *im2, grphic *implot);
void  imx_normalisation_gaussienne2d_histo_joint_norm_symetrique(void);

/*******  Quelques variantes dans l'estimation des gaussiennes  ********/

extern void 	estime_gaussienne2d_isosigma_EM			(grphic *histo,Param_Gaussienne_2d *param,Param_Gaussienne_2d *new_param,int nb_classe);
extern void 	estime_gaussienne2d_isosxsy_EM			(grphic *histo,Param_Gaussienne_2d *param,Param_Gaussienne_2d *new_param,int nb_classe);
extern void  	trimmed_histo_gaussienne_2d					(grphic *histo,grphic *trimmedhisto, Param_Gaussienne_2d *param,int nb_classe, double percentage);
extern void  	trimmed_histo_global_gaussienne_2d	(grphic *histo,grphic *trimmedhisto, Param_Gaussienne_2d *param,int nb_classe, double percentage);
extern void  	trimmed_histo_local_gaussienne_2d		(grphic *histo,grphic *trimmedhisto, Param_Gaussienne_2d *param,int nb_classe, double percentage);
extern void 	robust_estime_gaussienne2d_EM				(grphic *histo,Param_Gaussienne_2d *param,Param_Gaussienne_2d *new_param,int nb_classe);


/*******  Suppression de points non significatifs dans l'histogramme joint  ********/

extern void  trimmed_histo												(grphic *histo,grphic *trimmedhisto,double percentage);
extern void  trimmed_info_histo										(grphic *histo,grphic *trimmedhisto,double nb_sigma);
extern void  trimmed_info_xy_histo								(grphic *histo,grphic *trimmedhisto,double nb_sigma);
extern void  trimmed_robust_histo									(grphic *histo,grphic
*trimmedhisto,double nb_sigma);


/*******  Reclassement de points non significatifs de l'histogramme joint  ********/

extern void 	histo_joint_maxinfo_3d_p						(grphic3d *im1, grphic3d *im2, grphic *implot);
extern void 	histo_joint_maxinfo_multimodal_3d_p	(grphic3d *im1, grphic3d *im2, grphic *implot);


/*******  Divers fonctions relatives a l'EM  ********/

extern void  		plot_iso_gaussienne_2d	(grphic *histo,Param_Gaussienne_2d *param,double n_iso, int value);
extern double 	distance_mahalanobis_2d	(double x, double y, Param_Gaussienne_2d *param);
extern void 		verif_EM								(void);



/*******  Fonctions interfacant les differentes methodes de normalisation avec le recalage non-rigide  ********/

extern void 	TOP_norm_mean 												(grphic3d *im1, grphic3d *im2, grphic3d *imres, int nb_classe);
extern void 	TOP_norm_meanecty 										(grphic3d *im1, grphic3d *im2, grphic3d *imres, int nb_classe);
extern void 	TOP_norm_rl 													(grphic3d *im1, grphic3d *im2, grphic3d *imres, int nb_classe);
extern void 	TOP_norm_mean_ratio 									(grphic3d *im1, grphic3d *im2, grphic3d *imres, int nb_classe);
extern void 	TOP_norm_minmax 											(grphic3d *im1, grphic3d *im2, grphic3d *imres, int nb_classe);
extern void 	TOP_norm_spline 											(grphic3d *im1, grphic3d *im2, grphic3d *imres, int nb_classe);
extern void 	TOP_norm_GM 													(grphic3d *im1, grphic3d *im2, grphic3d *imres, int nb_classe);
extern void 	TOP_HistogramMatchingImageFilter_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres, int nb_classe);



/*******  Differentes methodes de normalisation   ********/

extern void 	normalisation_histo_joint					(void);
extern void 	normalise_histo_joint							(grphic3d *imreca,grphic3d *imref,grphic3d *imres,int nb_classe);

extern void 	normalisation_kmeans							(void);
extern void 	normalise_kmeans_p								(grphic3d *imreca,grphic3d *imref,grphic3d *imres,int nb_classe);

extern void 	kmeans_histo_joint_3d							(void);
extern void 	kmeans_histo_joint_3d_p						(grphic *im1,grphic *imres,int nb_classe);
extern void 	kmeans_L1_histo_joint_3d_p				(grphic *im1,grphic *imres,int nb_classe);

extern void 	normalisation_kmeans_histo_joint	(void);
extern void 	normalise_kmeans_histo_joint_p		(grphic3d *im1,grphic3d *im2,grphic3d *imres,int nb_classe);

extern void 	normalisation_markov							(void);
extern void 	normalise_markov_p 								(grphic3d *im1, grphic3d *im2, grphic3d *imres, int nb_classe);

extern void 	normalisation_local_3d_p					(grphic3d *im1, grphic3d *im2, grphic3d *imres);

extern void  normalisation_equalisation_histo		(void);
extern void  normalisation_equalisation_histo_p	(grphic3d *im1,grphic3d *im2,grphic3d *imres,int nb_classe);


/*******  normalisation avec segmentation sous-jacente  ********/
extern void  normalisation_segmentation				(void);
extern void  normalisation_segmentation_p			(grphic3d *im1,grphic3d *im2,grphic3d *imres,grphic *histo,int nb_classe);
extern void  imx_normalisation_segmentation_p	(grphic3d *im1,grphic3d *im2,grphic3d *imres,int nb_classe);
extern void  segmentation_histo								(grphic *histo,Param_Gaussienne_2d *param, int nb_classe);

extern void  imx_norm_seuil_meanecty_robust_3d(void);
extern void imx_norm_seuil_meanecty_robust_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres, int nb_classe);
void imx_norm_seuil_meanecty_symetrique_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *im1d, grphic3d *im2d, grphic3d *imres1, grphic3d *imres2);

double calcul_IM_histo(grphic *histo);

void equalisation_histogramme(grphic3d *imsrc,grphic3d *imres,double *ft,int size_ft);
void inv_equalisation_histogramme(grphic3d *imsrc,grphic3d *imres,double *ft,int size_ft);

int double_compare_function2(const void *a, const void *b);
double eval_binomiale(int N,double p,int k);

extern void segment_gaussienne2d_recherche_voisinage(grphic3d *imreca,grphic3d *imref,grphic3d *imreca_cl,grphic3d *imref_cl, grphic *histo, int nb_classe, Param_Gaussienne_2d *param,double ** confiance);
void segment_signe_covariance_classe_robuste(grphic3d *im1,grphic3d *im1_cl,grphic3d *im2,grphic3d *im2_cl,double *covariance,double *moyenne1, double
*moyenne2,double *std1, double *std2,int nb_classe,grphic *histo);
void apply_norm_mean_std_classe_optcov(grphic3d *im1,grphic3d *im1_cl,grphic3d *im2,grphic3d *imres,double *moyenne1,double *moyenne2,double *std1,double *std2,double *covariance, Param_Gaussienne_2d *param, int
nb_classe);

void apply_norm_regression_lineaire(grphic3d *im1,grphic3d *im1_cl,grphic3d *im2,grphic3d *imres,grphic *histo,int nb_classe);

extern void norm_histojoint_spline_Bosc(void);
extern void  normalisation_gaussienne2d_histo_joint_standart(void);
extern void  normalisation_gaussienne_histo(void);
extern void  normalisation_quantile_histo(void);
extern void  normalisation_gaussienne2D_histojoint_norm(void);
extern void  normalisation_moyenne_ecart_type_robust(void);
extern void Troncature_intensite_nsigma_grphic3d_p(grphic3d *im, grphic3d *imres, int nb_sigma);
extern void equalisation_histo_3d(void);
void equalisation_histo_3d_p(grphic3d *im1,grphic3d *imres);

#ifdef __cplusplus
}
#endif /*__cplusplus*/


#endif
