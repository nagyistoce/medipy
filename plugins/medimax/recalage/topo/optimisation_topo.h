/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef optimisation_topo_h
#define optimisation_topo_h

extern double TOP_linemin_locale_3d(grphic3d *imref, grphic3d *imreca, grphic3d *imres, field3d *champ_fin,transf_func_t transf, InterpolationFct inter, dist_func_locale_t dist, reg_func_locale_t regularisation,int nb_param, double *param, double *param_norm, double *direc, int topi, int topj, int topk, double Jm, double JM,TSlpqr *Slpqr);

extern double TOP_min_desc_grad_locale_3d(grphic3d *imref, grphic3d *imreca, grphic3d *imres, field3d *champ_fin, transf_func_t transf, InterpolationFct inter, dist_func_locale_t dist, reg_func_locale_t regularisation,int nb_param, double *param, int ***masque_param,double Jmin, double Jmax,char *nomfichiertrf);

extern double TOP_min_icm_locale_3d(grphic3d *imref, grphic3d *imreca, grphic3d *imres, field3d *champ_fin, transf_func_t transf, InterpolationFct inter, dist_func_locale_t dist, reg_func_locale_t regularisation, int nb_param, double *param, int ***masque_param,double Jmin, double Jmax,char *nomfichiertrf);

extern double TOP_min_marquardt_locale_3d(grphic3d *imref, grphic3d *imreca, grphic3d *imres, field3d *champ_fin, transf_func_t transf, InterpolationFct inter, dist_func_locale_t dist, reg_func_locale_t regularisation, int nb_param, double *param, int ***masque_param,double Jmin, double Jmax,char *nomfichiertrf);
extern double TOP_min_marquardt_locale_lent_3d(grphic3d *imref, grphic3d *imreca, grphic3d *imres, field3d *champ_fin, transf_func_t transf, InterpolationFct inter, dist_func_locale_t dist, reg_func_locale_t regularisation, int nb_param, double *param, int ***masque_param,double Jmin, double Jmax,char *nomfichiertrf);

extern double TOP_min_marquardt_penal_locale_3d(grphic3d *imref, grphic3d *imreca, grphic3d *imres, field3d *champ_fin, transf_func_t transf, InterpolationFct inter, dist_func_locale_t dist,reg_func_locale_t regularisation,  int nb_param, double *param, int ***masque_param,double Jmin, double Jmax,char *nomfichiertrf);

extern double TOP_min_marquardt_adapt_penal_locale_3d(grphic3d *imref, grphic3d *imreca, grphic3d *imres, field3d *champ_fin, transf_func_t transf, InterpolationFct inter, dist_func_locale_t dist, reg_func_locale_t regularisation, int nb_param, double *param, int ***masque_param,double Jmin, double Jmax,char *nomfichiertrf);

extern double TOP_min_marquardt_penal_locale_gradient_3d(grphic3d *imref, grphic3d *imreca, grphic3d *imres, field3d *champ_fin, transf_func_t transf, InterpolationFct inter, dist_func_locale_t dist, reg_func_locale_t regularisation, int nb_param, double *param, int ***masque_param,double Jmin, double Jmax,char *nomfichiertrf);

#endif
