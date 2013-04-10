/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/**-----------------------------------------------------------------------
***	
***	file:		tala_3d.h
***
***	project:	Imagix 1.01
***			
***
***	description:    Fichier relatif a l'atlas de Talairach.
***	
***	
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***
***---------------------------------------------------------------------*/

#ifndef _TALA_3D_H
#define _TALA_3D_H

/* chemins d'acces des fichiers de donn�s pour labellisation */
#define AREA_LABEL_DIR "/opt/medimaxlibs/v3.1/talairach/TTareas.txt"
#define GYRAL_LABEL_DIR "/opt/medimaxlibs/v3.1/talairach/TTgyral.txt"

//#include "noyau/imx_3d.h"

typedef struct s_dmatrice {
	int li;
	int co;
	double coeff[5][5];
	} dmatrice;


/* AUTOMATICALLY EXTRACTED PROTOTYPES (MB 6/12/2000) */
#ifdef __cplusplus 
extern "C" 
{
#endif /*__cplusplus*/ 
extern dmatrice somme_matrices(dmatrice a, dmatrice b) ;
extern dmatrice difference_matrices(dmatrice a, dmatrice b) ;
extern dmatrice produit_matrices(dmatrice a, dmatrice b) ;
extern dmatrice transposee(dmatrice a) ;
extern double norme(dpoint3d u) ;
extern double prod_scal(dpoint3d u, dpoint3d v) ;
extern dpoint3d prod_vect(dpoint3d u, dpoint3d v) ;
extern dpoint3d vecteur(dpoint3d A, dpoint3d B) ;
extern dpoint3d projection_ortho(dpoint3d MM, dplan plan) ;
extern dpoint3d centre_de_masse_p(grphic3d *im1) ;
extern double * cliquez_sur(char *nom_point, int im_1) ;
extern int determine_le_repere_p(grphic3d *im1) ;
extern void saisie_manuelle_p(grphic3d *im1, int im_1) ;
extern int simplex(grphic3d *im1, float **vertice, float *y, int dim, float (*distance) (/* ??? */), float ftol) ;
extern float contraction(grphic3d *im1, float coeff, float **vertice, float *y, float *psum, int dim, int ihi, float (*distance) (/* ??? */)) ;
extern void get_psum(float **vertice, int dim, float *psum) ;
extern void minimise_simplex(grphic3d *im1, float **vertice, float *y, int dim, float (*distance) (/* ??? */)) ;
extern void recherche_plan_simplex(grphic3d *im1, dplan *plan, dpoint3d *point, float (*distance) (/* ??? */)) ;
extern dplan defini_plan_param_simplex(float *param) ;
extern dplan defini_plan_param_ICM(dpoint3d point, dpoint3d u, dpoint3d v) ;
extern void recherche_plan_ICM(grphic3d *im1, dplan *plan, dpoint3d *P, float (*distance) (/* ??? */)) ;
extern dpoint3d rotation_vect(dpoint3d axe, double angle, dpoint3d vect1) ;
extern float dist_m_c(grphic3d *im1, dplan plan) ;
extern float dist_McClure(grphic3d *im1, dplan plan) ;
extern float dist_cross(grphic3d *im1, dplan plan) ;
extern void trace_plan(grphic3d *im1, dplan plan) ;
extern void determine_rotation(grphic3d *im1, dplan plan, dmatrice *A, dmatrice *b, dmatrice *invA) ;
extern void transfo_image3d_proche_voisin(grphic3d *im1, grphic3d *im2, dmatrice A, dmatrice b, dmatrice invA) ;
extern void transfo_image3d_trilineaire(grphic3d *im1, grphic3d *im2, dmatrice A, dmatrice b, dmatrice invA) ;
extern dpoint3d point_antecedent_transfo(dpoint3d point, dmatrice A, dmatrice b, dmatrice invA) ;
extern void repere_antecedent_transfo_p(grphic3d *im1, grphic3d *im2, dmatrice A, dmatrice b, dmatrice invA) ;
extern int position_relative(double a, double b, double x, int n) ;
extern int coord_ds_grille(grphic3d *im1, dpoint3d point) ;
extern int image_seuillee_p(grphic3d *im1) ;
extern void recherche_bords_cerveau_p(grphic3d *im1) ;
extern void tala_manu(void) ;
extern int recale_plan_sagittal_p(grphic3d *im1, int im_1) ;
extern void tala_semi2(void) ;
extern void tala_semi2bis(grphic3d *im1, int im_1) ;
extern void tala_semi8(void) ;
extern void tala_auto(void) ;
extern void tala_coord(void) ;
extern void tala_init(void) ;
extern void tala_recale_images_p(int choix, grphic3d *im1, grphic3d *im2, grphic3d *im3) ;
extern void determine_transfo_cube(dpoint3d O1, dpoint3d u1, dpoint3d v1, dpoint3d w1, dpoint3d O2, dpoint3d u2, dpoint3d v2, dpoint3d w2, dmatrice *A, dmatrice *b, dmatrice *invA) ;
extern void transfo_cube_p(int choix, grphic3d *im1, dpoint3d O1, dpoint3d u1, dpoint3d v1, dpoint3d w1, grphic3d *im2, dpoint3d O2, dpoint3d u2, dpoint3d v2, dpoint3d w2) ;
extern void transfo_cube_proche_voisin_p(grphic3d *im1, dpoint3d O1, dpoint3d u1, dpoint3d v1, dpoint3d w1, grphic3d *im2, dpoint3d O2, dpoint3d u2, dpoint3d v2, dpoint3d w2) ;
extern void transfo_cube_trilineaire_p(grphic3d *im1, dpoint3d O1, dpoint3d u1, dpoint3d v1, dpoint3d w1, grphic3d *im2, dpoint3d O2, dpoint3d u2, dpoint3d v2, dpoint3d w2) ;
extern void creer_cube(dpoint3d O, dpoint3d *u, dpoint3d *v, dpoint3d *w, dpoint3d Tu, dpoint3d Tv, dpoint3d Tw, double x, double y, double z) ;
extern void tala_recale_images(void) ;
extern void tala_grid_on(void) ;
extern void tala_grid_off(void) ;
extern void tala_rm_grid(void) ;
extern void tala_copy_grid(void) ;
extern void tala_copy_grid_p(grphic3d *im1, grphic3d *im2) ;
extern void Tala_Save_3d(void) ;
extern void position_ds_talairach(grphic3d *im, dpoint3d point, char *reponse);
extern void position_ds_talairach_mm(grphic3d *im, dpoint3d point, char *reponse, char *lblArea, char *lblGyrus);
extern int loadLabels(grphic3d *im);
extern void extraction_zones_talairach(void);

typedef struct _LabelsBAGyri
{
	int BA;
	int Gyri;
}LabelsBAGyri;


#ifdef __cplusplus 
}
#endif /*__cplusplus*/ 

#endif /* Fin _TALA_3D*/
