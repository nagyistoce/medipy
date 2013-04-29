/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef validation_topo_3d_h
#define validation_topo_3d_h

void superposition_damier_3d(void);
void superposition_damier_3d_p(int im_reca,int im_ref,int im_res,int step);


void synthese_champ_topo_3d(void);
void rand_champ_topo(char *nomfichres,int width, int height, int depth, float dx, float dy, float dz, int size_min);
void  id_champ(char *nomfichres,int width, int height, int depth, float dx, float dy, float dz);
void topo_synth_3d(void);
void rand_topo_field(char *nomfichres,int taille, int resol, double Jm, double JM);                                              
void chps_test_3d(char *nomfichres);


void distance_landmark(void);
void distance_landmark_p(grphic3d *im1,grphic3d *im2,grphic3d *imres,field3d *champ,char* nomfichres);

void info_landmark(void);
void info_landmark_p(grphic3d *im1, field3d* champ);

void LVV_operator(void);
void LVV_operator_p(grphic3d *im1,grphic3d *imres,float sigma);

void verif_topo(void);

void compare_critere(grphic3d *imref, grphic3d *imreca, grphic3d *imres,field3d *champ_fin,double *param, double *param_norm, int nb_param, int ***masque_param, int nb_iter);
extern void distances_inter_images_3d(void);
void imx_distances_inter_images_3d(int im_1, int im_2);

void correlation_inter_images_3d(void);
void imx_correlation_inter_images_3d(int im_1, int im_2);
#endif
