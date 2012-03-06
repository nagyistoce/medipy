/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef divers_topo_h
#define divers_topo_h

#include "recalage/topo/analyse_intervalle.h"
#include "recalage/definitions_types_recalage.h"


#define PRECISION_INVERSION 0.1



/***************************** DEFINITION DU TYPE boite_inv_bspline **********************/
typedef struct s_BoiteInvBspline {
int resol; 									// Resolution du modele de deformation Bspline
int topi,topj,topk; 				// Indices de la boite considere
double x0,y0,z0,xm,ym,zm;		// Coordonnées des bornes de la boite
double ax000,ax001,ax010,ax100,ax110,ax101,ax011,ax111;  	// Coefficients suivant x des fonctions Bsplines intervenant sur la boite
double ay000,ay001,ay010,ay100,ay110,ay101,ay011,ay111;		// Coefficients suivant y des fonctions Bsplines intervenant sur la boite
double az000,az001,az010,az100,az110,az101,az011,az111;		// Coefficients suivant z des fonctions Bsplines intervenant sur la boite
double norm;							// Coefficient de normalisation des splines
} BoiteInvBspline, *ptr_BoiteInvBspline;


/***************************** DEFINITION DU TYPE mat3d **********************/
typedef struct s_mat3d {
 float H[3][3];
} mat3d, *ptr_mat3d;

int inversion_mat3d	(mat3d *A, mat3d* invA);
double determinant_mat3d(mat3d *A);

int	jacobien_bspline1_3d							(transf3d *transfo, grphic3d *imres);
int	borne_jacobien_bspline1_local_3d	(int nb_param, double *p, double *param_norm,int topi, int topj, int topk,double *min,double *max);

void test_correlation(grphic3d *imref, grphic3d *imres,int nb_param,int ***masque_param);
void test_correlation2(grphic3d *imref, grphic3d *imres,int nb_param,int ***masque_param);

double 	dim_min_box	(TOPTbox b);
TOPTbox find_box		(double x, double y, double z, TOPTliste myListe, int *status);

void 	save_histo		(void);
int		save_histo_3d	(int im_1,float xmax,float xmin,int Nbpas,char *nomfichres);

double my_funct	(double x, double y, double z, TOPTliste liste1, TOPTliste liste2, int i, int j,int k  );
double my_funct2(double x, double y, double z, TOPTliste liste1, int i, int j,int k );


void	inverse_transf_3d	(void);
void	inverse_transf_3d_p(transf3d *transfo1,transf3d *transfores,double prec);
int 	inv_affine_3d			(transf3d*,transf3d*);
int 	inv_field_3d			(transf3d *transfo,field3d *chres, double prec);
int 	inv_bspline_3d			(transf3d* transfo,field3d *chres,double prec);
transf3d* ConvertTransfoField3dToBase3d(transf3d *transfo_field3d);
int ParamBspline_to_field_3d(int nb_param, double *param, field3d *champ, grphic3d *imref, grphic3d *imreca);

void deform_anim_3d(void);

void filtre_champ_diffusion_anisotrope_3d();
int imx_filtre_champ_diffusion_anisotrope_3d(int im_src, char* nomfich1, char* nomfichres, int nb_iter,double Kd);
int imx_filtre_champ_diffusion_anisotrope_3d_p(grphic3d *imsrc,field3d* ch,field3d* chres,int nb_iter,double Kd);

void filtre_champ_gaussien_3d(void);
int filtre_champ_gaussien_3d_p(char* nomfich1, char* nomfichres, int wnd,double sigma);
void imx_filtre_champ_gaussien_3d_p(field3d *imdeb, field3d *imres, int wnd, double sigma);
void imx_filtre_param_gaussien_3d_p(double *param_deb, double *param_res, int nb_param, double sigma);

double raffinement_bspline_inv_3d(double *param,int nb_param,field3d *chdepart, field3d *chres,field3d *chmin,field3d *chmax,int i,int j,int k,double prec);
void eval_deplacement_bspline1_3d(double *param, int nb_param, int wdth, int hght, int dpth, double x, double y, double z, dvector3d* u);
void eval_transfo_bspline1_3d(double *param, int nb_param, int wdth, int hght, int dpth, double x, double y, double z, dvector3d* u);
void eval_transfo_bspline1_BoiteInvBspline_3d(BoiteInvBspline *ParamInvBspline, double x, double y, double z, dvector3d* u);
void eval_deplacement_bspline1_BoiteInvBspline_3d(BoiteInvBspline *ParamInvBspline, double x, double y, double z, dvector3d* u);
void eval_jacobien_transfo_bspline1_3d(double *param, int nb_param, int wdth, int hght, int dpth, double x, double y, double z, double **J);
void eval_jacobien_transfo_bspline1_BoiteInvBspline_3d(BoiteInvBspline *ParamInvBspline, double x, double y, double z, double** J);
double inv_bspline_reduction_boite_3d(double *param, int nb_param, int wdth, int hght, int dpth, double precis, double i, double j, double k, TOPTbox* myB);
void grignote_eval_proposition_x(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double y, double z,double i,double j,double k, double *x1, double *x2, double *x3);
void grignote_eval_proposition_y(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double x, double z,double i,double j,double k, double *y1, double *y2, double *y3);
void grignote_eval_proposition_z(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double x, double y,double i,double j,double k, double *z1, double *z2, double *z3);
double grignote_boite_suivant_x(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i,double j,double k);
double grignote_boite_suivant_y(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i,double j,double k);
double grignote_boite_suivant_z(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i,double j,double k);
double raffinement_bspline_inv_3d2(double *param,int nb_param,field3d *chdepart,field3d *chres,field3d *chmin,field3d *chmax,int i,int j,int k,double prec);
void update_bloc_transfo_bspline1_3d(double *param, int nb_param, int wdth, int hght, int dpth, TOPTbox* myB);
void inv_bspline_casser_boite_3d(double *param, int nb_param, int width, int height, int depth, TOPTbox *myB, TOPTliste *myQ, double xmid, double ymid, double zmid, dvector3d hmid);

int grignote_boite_update_xm(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double x , double i,double j,double k);
int grignote_boite_update_xM(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double x , double i,double j,double k);
int grignote_boite_update_ym(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double x , double i,double j,double k);
int grignote_boite_update_yM(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double x , double i,double j,double k);
int grignote_boite_update_zm(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double x , double i,double j,double k);
int grignote_boite_update_zM(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double x , double i,double j,double k);

int comp_double(const void *num1, const void *num2);

void	itkInverseDeformationFieldImageFilter_3d(void);
int inv_bspline_find_solution_with_fixed_point_algo(double *param, int nb_param, int wdth, int hght, int dpth, double precis, double i, double j, double k, TOPTbox* myB, dvector3d* res, double *distance);

int test_if_x_suivant_hx_is_updatable(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i, double x);
int test_if_x_suivant_hy_is_updatable(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i, double x);
int test_if_x_suivant_hz_is_updatable(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i, double x);

int test_if_y_suivant_hx_is_updatable(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i, double x);
int test_if_y_suivant_hy_is_updatable(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i, double x);
int test_if_y_suivant_hz_is_updatable(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i, double x);

int test_if_z_suivant_hx_is_updatable(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i, double x);
int test_if_z_suivant_hy_is_updatable(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i, double x);
int test_if_z_suivant_hz_is_updatable(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i, double x);

double solve_x_suivant_hx(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double y, double z, double i);
double solve_x_suivant_hy(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double y, double z, double j);
double solve_x_suivant_hz(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double y, double z, double k);

double solve_y_suivant_hx(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double x, double z, double i);
double solve_y_suivant_hy(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double x, double z, double j);
double solve_y_suivant_hz(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double x, double z, double k);

double solve_z_suivant_hx(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double x, double y, double i);
double solve_z_suivant_hy(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double x, double y, double j);
double solve_z_suivant_hz(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double x, double y, double k);


int update_x_suivant_hx(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i, int is_xm, int is_xM);
int update_x_suivant_hy(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i, int is_xm, int is_xM);
int update_x_suivant_hz(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i, int is_xm, int is_xM);

int update_y_suivant_hx(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i, int is_xm, int is_xM);
int update_y_suivant_hy(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i, int is_xm, int is_xM);
int update_y_suivant_hz(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i, int is_xm, int is_xM);

int update_z_suivant_hx(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i, int is_xm, int is_xM);
int update_z_suivant_hy(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i, int is_xm, int is_xM);
int update_z_suivant_hz(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i, int is_xm, int is_xM);


double eval_deplacement_bspline1_BoiteInvBspline_3d_x(BoiteInvBspline *ParamInvBspline, double x, double y, double z);
double eval_transfo_bspline1_BoiteInvBspline_3d_x(BoiteInvBspline *ParamInvBspline, double x, double y, double z);
double eval_deplacement_bspline1_BoiteInvBspline_3d_y(BoiteInvBspline *ParamInvBspline, double x, double y, double z);
double eval_transfo_bspline1_BoiteInvBspline_3d_y(BoiteInvBspline *ParamInvBspline, double x, double y, double z);
double eval_deplacement_bspline1_BoiteInvBspline_3d_z(BoiteInvBspline *ParamInvBspline, double x, double y, double z);
double eval_transfo_bspline1_BoiteInvBspline_3d_z(BoiteInvBspline *ParamInvBspline, double x, double y, double z);

int inv_bspline_find_solution_with_fixed_point_algo_global(double *param, int nb_param, int wdth, int hght, int dpth, double precis, double i, double j, double k, TOPTbox* myB, dvector3d* res, double *distance);

int acknowledge_invbspline(double *param, int nb_param, int width, int height, int depth, double i, double j, double k, TOPTbox myB);			
int ProjectBase3dintoTopologyPreservingTransformation(transf3d *transfo_base3d );
int	base_resol_up_light_3d(double *param, int nb_param);


void fusion_transf_3d();
void imx_fusion_transf_3d(int im_1, char* nomfichier1, int im_2, char* nomfichier2, char* nomfichres);
void imx_fusion_transf_3d_p(grphic3d *im1, field3d *ch1, grphic3d *im2, field3d *ch2, field3d *chres);


void convert_IPB_to_trf();
void imx_convert_IPB_to_trf(int im_x, int im_y, int im_z, char* nomfichres);
void imx_convert_IPB_to_trf_p(grphic3d *imx, grphic3d *imy, grphic3d *imz, field3d* chres);

void ComputeMoment3D(void);
void imx_ComputeMoment3D(int im_deb, int im_res, int p, int q, int r, float rayon);
void imx_ComputeMoment3D_p(grphic3d * imdeb, grphic3d * imres, int p, int q, int r, float rayon);
void ComputeRotationInvariantMoment3D(void);
void imx_ComputeRotationInvariantMoment3D(int im_deb, int im_res, int num, float rayon);
void imx_ComputeRotationInvariantMoment3D_p(grphic3d * imdeb, grphic3d * imres, int num, float rayon);

void stat_fich_histo(void);
void Compute_energie_groupwise_variance(void);

void NLMWeightOnePoint3D(grphic3d* Image, double *** weights, float NLMsmooth, int NLMhwnx, int NLMhwny, int NLMhwnz, int NLMhwvsx, int NLMhwvsy, int NLMhwvsz,  int x, int y, int z);

float NLMSmoothComputation3D(grphic3d* Image, int NLMhwnx, int NLMhwny, int NLMhwnz, float NLMbeta, int padding);

int	integrale_jacobien_bspline1_compute(transf3d *transfo, grphic3d *imres);

#endif
