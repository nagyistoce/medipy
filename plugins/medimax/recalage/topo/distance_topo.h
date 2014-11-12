/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef distance_topo_h
#define distance_topo_h



extern double Energie_quad_locale_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, double *param_norm, int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation);
extern double Energie_quad_sym_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation);

extern double Energie_Lp_locale_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, double *param_norm, int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation);
extern double Energie_Lp_sym_locale_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, double *param_norm, int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation);

extern double Energie_L1L2_locale_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, double *param_norm, int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation);
extern double Energie_L1L2_sym_locale_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, double *param_norm, int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation);

extern double Energie_quad_topo_locale_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, double *param_norm, int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation);
extern double Energie_quad_sym_topo_locale_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, double *param_norm, int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation);

extern double Energie_L1norm_locale_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, double *param_norm, int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation);

extern double Energie_geman_locale_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, double *param_norm, int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation);
extern double Energie_geman_sym_locale_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, double *param_norm, int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation);

extern double Energie_globale_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, double *param_norm, dist_func_locale_t dist,reg_func_locale_t regularisation,int ***masque_param);
extern double Energie_globale_sym_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, double *param_norm, dist_func_locale_t dist,int ***masque_param,reg_func_locale_t regularisation);


extern double Energie_groupwise_variance_locale_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, double *param_norm, int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation);
extern double Energie_groupwise_variance_locale_nonsym_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation);

void update_TOPO_ALPHA_ROBUST_global(grphic3d *imref, grphic3d *imreca, int nb_param, double *param);
void update_TOPO_ALPHA_ROBUST_local(grphic3d *imref, grphic3d *imreca, int nb_param, double *param,int topi, int topj, int topk);

double regularisation_energie_membrane_local(int nb_param, double *param,double *param_norm,int topi, int topj, int topk,grphic3d *mask, TSlpqr *Slpqr);
double regularisation_energie_membrane_Lp_local(int nb_param, double *param,double *param_norm,int topi, int topj, int topk, grphic3d *mask, TSlpqr *Slpqr);
double regularisation_globale_3d(int nb_param, double *param,double *param_norm, int ***masque_param, grphic3d *mask, reg_func_locale_t regularisation);
double regularisation_log_jacobien_local(int nb_param, double *param,double *param_norm,int topi, int topj, int topk,grphic3d *mask, TSlpqr *Slpqr);
double regularisation_energie_membrane_jacobien_local(int nb_param, double *param,double *param_norm, int topi, int topj, int topk,grphic3d *mask, TSlpqr *Slpqr);
double regularisation_log_jacobien_centre_local(int nb_param, double *param,double *param_norm, int topi, int topj, int topk,grphic3d *mask, TSlpqr *Slpqr);
double regularisation_dist_identite_local(int nb_param, double *param,double *param_norm, int topi, int topj, int topk,grphic3d *mask, TSlpqr *Slpqr);

extern double Energie_quad_locale_pondere_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, double *param_norm, int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation);
extern double Energie_ML_locale_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, double *param_norm, int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation);


double Energie_IM_locale_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, double *param_norm, int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation);
double Calcul_IM_global_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param);

double Energie_ICP_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation);
double Energie_ICP_sym_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation);

double Energie_atrophie_jacobien_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation);
double Energie_atrophie_log_jacobien_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation);


extern double Energie_quad_locale_symetrique_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, double *param_norm, int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation);
extern double Energie_quad_locale_symetrique_coupe_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, double *param_norm, int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation);


extern double Energie_quad_sous_champ_locale_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, double *param_norm, int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation);

extern double Energie_DTI_quad_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation);

double Energie_patch_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation);

double image_trilinear_value(grphic3d *im, double x, double y, double z);

extern double regularisation_patch_imagebased_local(int nb_param, double *param,double *param_norm, int topi, int topj, int topk,grphic3d *mask, TSlpqr *Slpqr);

#endif
