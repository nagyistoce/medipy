/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef hessien_distance_topo_h
#define hessien_distance_topo_h

#include "recalage/topo/divers_topo.h"
#include "recalage/topo/mtch_topo_3d.h"






/***** ALLOCATION ET LIBERATION MEMOIRE ********/
hessien3d *cr_hessien3d(int wdth, int hght, int dpth);
int	free_hessien3d(hessien3d *ch);
int	mul_hessien_3d(hessien3d *hessien, hessien3d *hessienres, float coeff);

int imx_hessien_3d_p(grphic3d *im, hessien3d *hessien, int method, int t);


int hessien_base_quad_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2,hessien3d *hessienI2,hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg);
int hessien_base_quad_sym_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, hessien3d *hessienI2,hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg);

int hessien_base_L1L2_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, hessien3d *hessienI2,hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg);
int hessien_base_L1L2_sym_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, hessien3d *hessienI2,hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg);

int hessien_base_L1norm_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, hessien3d *hessienI2,hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg);

int hessien_base_geman_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, hessien3d *hessienI2,hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg);

int hessien_regularisation_energie_membrane_local(int nb_param, double *param, double *param_norm, mat3d *hessien, int topi, int topj, int topk,grphic3d *mask, TSlpqr *Slpqr);
int hessien_regularisation_energie_membrane_Lp_local(int nb_param, double *param,double *param_norm, mat3d *hessien, int topi, int topj, int topk,grphic3d *mask, TSlpqr *Slpqr);
int hessien_regularisation_log_jacobien_local(int nb_param, double *param,double *param_norm, mat3d *hessien, int topi, int topj, int topk,grphic3d *mask, TSlpqr *Slpqr);

int hessien_base_quad_locale_pondere_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, hessien3d *hessienI2,hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg);

int hessien_base_IM_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, hessien3d *hessienI2,hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg);

int hessien_base_ICP_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, hessien3d *hessienI2,hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg);
int hessien_base_ICP_sym_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, hessien3d *hessienI2,hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg);

int hessien_atrophie_jacobien_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, hessien3d *hessienI2,hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg);
int hessien_atrophie_log_jacobien_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, hessien3d *hessienI2,hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg);

int hessien_quad_locale_symetrique_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2,hessien3d *hessienI2,hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg);
int hessien_quad_locale_symetrique_coupe_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2,hessien3d *hessienI2,hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg);

int hessien_quad_sous_champ_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, hessien3d *hessienI2,
hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg);

int hessien_groupwise_variance_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2,hessien3d *hessienI2,hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg);
int hessien_groupwise_variance_locale_nonsym_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, hessien3d *hessienI2,
hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg);


int hessien_regularisation_dist_identite_local(int nb_param, double *param, double *param_norm, mat3d *hessien, int topi, int topj, int topk,grphic3d *mask, TSlpqr *Slpqr);
int hessien_regularisation_patch_imagebased_local(int nb_param, double *param, double *param_norm, mat3d *hessien, int topi, int topj, int topk,grphic3d *mask, TSlpqr *Slpqr);

#endif
