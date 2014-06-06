/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef analyse_intervalle_h
#define analyse_intervalle_h
#include "recalage/definitions_types_recalage.h"


#define MIN2(x,y) (((x)<(y)) ? (x) : (y))
#define MAX2(x,y) (((x)>(y)) ? (x) : (y))

extern int RENVERSEMENT_SEUIL;

/***************************** DEFINITION DU TYPE Tdroite **********************/
typedef struct s_Tdroite { 
double pente, ordorg; 
} Tdroite, *ptr_Tdroite;

/***************************** DEFINITION DU TYPE Tparall **********************/
typedef struct s_Tparall {
double xm, xM, ym, yM, zm, zM; 
} Tparall, *ptr_Tparall;

/***************************** DEFINITION DU TYPE TSlpqr ***********************/
typedef struct s_TSlpqr { 
int ind;												// indice des coefficients de la deformation que l'on modifie
double bdeltam, bdeltaM;    		// borne min et max
double Jm, JM;									// bornes sur le jacobien
Tparall b;											// coordonnees de la boite
Tdroite alpha[20];							// expression affine des differents coefficient du jacobien
double alx[8], aly[8], alz[8];	// valeurs des coefficients de la deformation sur la boite
} TSlpqr, *ptr_TSlpqr;

/***************************** DEFINITION DU TYPE TOPTbox ***********************/
typedef struct s_TOPTbox
{ double xm,xM,ym,yM,zm,zM;
  double fm,fM;
	double ax,ay,az;
  struct s_TOPTbox *next;
	dvector3d hx[8];
} TOPTbox;

/***************************** DEFINITION DU TYPE TOPTbox ***********************/
typedef struct s_TOPTliste
{ TOPTbox *tete;
  unsigned long nelem;
} TOPTliste;


/************************** FONCTION RALATIVE A LA CONSERVATION TOPOLOGIE (ANALYSE PAR INTERVALLES  *************************************/
extern double TOP_linemin_maxpas_3d			(int nb_param,double *param,double *direc,int i,int j,int k, double Jm, double JM, TSlpqr *Slpqr);
extern void 	TOP_poly_pente_coeff			(double *dir_des, int indice, double *alx, double *aly, double *alz, int l, double *pente);
extern void 	TOP_initialization_delta	(TSlpqr *box);
extern double TOP_eval_poly							(double xt, double yt, double zt, Tdroite *droite, double cpente, double cordorg);
extern void 	TOP_init_slpqr						(TSlpqr *box,int p,int q, int r, int D, int nb_param, double *param);
void 					TOP_poly_init							(double *alx, double *aly, double *alz, int l, double *alpha);
void 					poly_rescale							(double *a, double a1, double b1, double a2, double b2, double a3, double b3, double *alpha);


void     master_minim_main   (double seuil, double *a, Tparall b, double *binf, double *bsup, double *bsupx, double *bsupy, double *bsupz);
void     minim_main          (double seuil, double precis, double *a, double xinf, double xsup, double yinf, double ysup, double zinf, double zsup, double *binf, double *bsup, double *bsupx, double *bsupy, double *bsupz,double J);
void     pave_grignote       (double precis,double *a,TOPTbox *b, int *status);
void     poly_grad           (double a[], double mygrad[][4]);
double   dim_max_parall      (Tparall p);
double   dim_max_box         (TOPTbox b);
void     poly_marg           (int indice, double *a, double xt, double yt, double zt, double *marg);
void     poly_marg_poly_fast (int param, double my_grad[][4], double big_mat[][10]);
void     if_natural          (double *a, double xinf, double xsup, double yinf, double ysup, double zinf, double zsup, double *auxdbl2);
void     poly_transl         (double *a, double b1, double b2, double b3, double *alpha);
void     eval_fun            (double *a, double xt, double yt, double zt, double *aux, double *mymono);
void 		 fast_eval_fun			 (double *a, double xt, double yt, double zt, double *aux, double *mymono);
void eval_integrale_fun(double *a, double x0, double y0, double z0, double x1, double y1, double z1, double *aux);
TOPTbox  depile_end          (TOPTliste *myListe);
void     empile              (TOPTliste *myL, TOPTbox myB);
void     empile_front        (TOPTliste *myL, TOPTbox myB);
void     minim_iter          (double precis, double *a, TOPTbox myB, double *x_hat, double *y_hat, double *z_hat);
double   minim_min_ligne     (double *poly_marg,double vinf, double vsup);
void     pave_casser         (double precis, TOPTbox myB, double *xcass, double *ycass, double *zcass, TOPTliste *myL);
void     vecteur_insere      (int *ncomp, double *lin, double val, double precis);
int      acknowledge         (TSlpqr *Slpqr, double delta_s, int minimiser);
void     make_propo          (TSlpqr *Slpqr, double delta_s, int minimiser, double *aux_i, double *aux_s);
void     TOP_refinement      (TSlpqr *Slpqr, double *delta_i, double *delta_s);


extern void 	TOP_verification_delta	(TSlpqr *Slpqr, double coeff, double *result);
extern int 		TOP_verification_all		(int nb_param,double *param,int ***masque_param,double Jm,double JM);
int 					TOP_verification_locale	(int nb_param,double *param,double *param_norm,double Jm,double JM,int i,int j,int k,TSlpqr *Slpqr);
void 					verif_poly							(int i,int j,int k,int ind,int echelle,int nb_param,double *param,TSlpqr *Slpqr);

#endif
