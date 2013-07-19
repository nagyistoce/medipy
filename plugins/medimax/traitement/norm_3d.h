/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*!**********************************************************************
***	
***	\file		norm_3d.h
***
***	project:	Imagix 2.01
***			
***
***	\brief description:	normalisation d'images
***	
***	
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***     Last user action by: Mr. ARMSPACH on Dec 15th 1996
***
*************************************************************************/
#ifndef __NORM_3D_H
#define __NORM_3D_H

#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/

//prototype d'une fonction de normalisation
typedef void (*norm_func_t)(grphic3d *im1, grphic3d *im2, grphic3d *imres);

//---methodes de normalisation d'images masquees---//
//normalisation d'images non recalees
extern void norm_roi_3d(void) ;
extern void imx_norm_roi_3d (int im_1, int im_2, int im_res, int type_norm);
//normalisation d'images recalees
extern void norm_roi_reca_3d(void) ;
extern void imx_norm_roi_reca_3d (int im_1, int im_2, int im_res, int type_norm);
//methodes de normalisation
extern void imx_norm_roi_mean_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern void imx_norm_roi_meanecty_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern void imx_norm_roi_rl_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern void imx_norm_roi_rl_1_param_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern void imx_norm_roi_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern void imx_norm_roi_minmax_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
//toolbox pour les methodes de normalisation
extern int imx_calc_taille_mask_3d_p(grphic3d *im);
extern int imx_calc_moy_im_mask_3d_p(grphic3d *im, double *moyenne);
extern int imx_calc_sigma_im_mask_3d_p(grphic3d *im, int taille, double moyenne, double *sigma);
extern int imx_calc_covar_im_mask_3d_p(grphic3d *im1, grphic3d *im2, double moyenne1, double moyenne2, double *covariance);
extern int imx_norm_lin_im_3d_p(grphic3d *im1, grphic3d *imres, double alpha, double beta);

//autres normalisations
extern void norm_im_3d(void) ;
extern void imx_norm_im_3d (int im_1, int im_2, int im_res, int type_norm);
extern void imx_norm_im_mean_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern void imx_norm_im_meanecty_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern void imx_norm_im_rl_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern void imx_norm_im_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern void imx_norm_im_minmax_3d_p(grphic3d *im1,grphic3d *im2,grphic3d *imres);
extern void norm_seuil_3d(void) ;
extern void imx_norm_seuil_3d (int im_1, int im_2, int im_res, int type_norm);
extern void imx_norm_seuil_mean_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern void imx_norm_seuil_meanecty_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern void imx_norm_seuil_rl_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern void imx_norm_seuil_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern void imx_norm_seuil_minmax_3d_p(grphic3d *im1,grphic3d *im2,grphic3d *imres);
extern int  imx_norm_rcoeff_3d_p (grphic3d *im1, grphic3d *im2);
extern void norm_local_3d(void) ;
extern void imx_norm_local_3d (int im_1, int im_2, int im_res, int N);
extern int  imx_norm_local_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres, int N);
extern void norm_coronal_3d(void) ;
extern void imx_norm_coronal_3d (int im_1, int im_2, int im_res);
extern void imx_norm_coronal_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern void norm_sagittal_3d(void) ;
extern void imx_norm_sagittal_3d (int im_1, int im_2, int im_res);
extern void imx_norm_sagittal_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern void norm_transverse_3d(void) ;
extern void imx_norm_transverse_3d (int im_1, int im_2, int im_res);
extern void imx_norm_transverse_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern void norm_tarik_3d(void) ;
extern void imx_norm_tarik_3d (int im_1, int im_2, int im_res, int N);
extern int  imx_norm_tarik_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres, int N);
extern int imx_norm_rcoeff_3d_p(grphic3d *im1, grphic3d *im2);
void normalise_mixture_gaussienne(void);
void CPP_normalisation_histo(grphic3d *im1,grphic3d *im2,grphic3d *imres);
void CPP_normalisation_GM(grphic3d *im1,grphic3d *im2,grphic3d *imres,int nb_gaussienne, double facteur);
void norm_seuil_coeff_3d(void);
void imx_norm_seuil_coeff_3d(int im_1,int im_res,float mean,float std);
void imx_norm_seuil_coeff_3d_p(grphic3d * im,grphic3d * imres,float mean,float std);

//---divers
extern int imx_query_normalisation_type_noReca_3d();
extern int imx_query_normalisation_type_reca_3d();
extern norm_func_t imx_choose_normalisation_roi_noReca_fct(int norm_type);
extern norm_func_t imx_choose_normalisation_roi_reca_fct(int norm_type);

extern void imx_corr_sagittal_3d(int im_ori, int im_res);
extern void imx_corr_sagittal_3d_p(grphic3d *imori, grphic3d *imres);

#ifdef __cplusplus
}
#endif /*__cplusplus*/

#endif /* __NORM_3D_H*/
