/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef IMX_MATH_H
#define IMX_MATH_H

#ifdef __cplusplus
extern "C" {
#endif


#ifndef M_PI
# define M_PI	3.14159265358979323846264338327950288
#endif

#ifndef SIN60
# define SIN60	0.86602540378443865	/* sin(60 deg) */
# define COS72	0.30901699437494742	/* cos(72 deg) */
# define SIN72	0.95105651629515357	/* sin(72 deg) */
#endif

# define REAL		double

float _t_swap;
#define SWAPF(f1, f2) \
{ _t_swap = f1; f1 = f2; f2 = _t_swap; }
float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

/* parametres de la fonction a minimiser */
typedef struct s_param_simplex {

	int nb_point;       /* nb points des tableaux */
	float *tab_val_in;  /* entrees pour la fct de transfert */
	float *tab_val_out; /* sortie de la fct transfert */
	float *x;           /* abscisse de la fct */

} Param_simplex;

/* Represente 1 point a traiter */
typedef struct vtx{

	float *valeur_estimer;  /* tableau des valeurs a tester */
	int nb_valeur_estimer;
	float r;                /* residu */

} Vertex;



extern void init_rand_seed(void);

void nrerror(char *msg);
float *vectorF(long int nl, long int nh);
void free_vectorF(float *v, long int nl, long int nh);


void four1(float *data, long unsigned int nn, int isign);
void twofft(float *data1, float *data2, float *fft1, float *fft2, long unsigned int n);
void convlv(float *data, long unsigned int n, float *respns, long unsigned int m, int isign, float *ans);
void realft(float *data, long unsigned int n, int isign);

void simplex_gen(void (*function) (/* ??? */), Param_simplex f, Vertex *vtx, int nbv);
void swap_vtx(Vertex *v1, Vertex *v2);
void init_simplex(Param_simplex *fun, int nb_pts, Vertex *vtx, int nb_vtx, int nb_param);
void free_simplex(Param_simplex *fun, Vertex *vtx, int nb_vtx);

extern void nrerror(char *msg) ;
extern float *vectorF(long int nl, long int nh) ;
extern void free_vectorF(float *v, long int nl, long int nh) ;
extern void four1(float *data, long unsigned int nn, int isign) ;
extern void twofft(float *data1, float *data2, float *fft1, float *fft2, long unsigned int n) ;
extern void convlv(float *data, long unsigned int n, float *respns, long unsigned int m, int isign, float *ans) ;
extern void realft(float *data, long unsigned int n, int isign) ;
extern void fft_free (void) ;
extern int fftn (int ndim, int *dims, double *Re, double *Im, int iSign) ;
/*  extern static int FFTRADIX (double *Re, double *Im, unsigned int nTotal, unsigned int nPass, unsigned int nSpan, int iSign, int max_factors, int max_perm) ; */
extern void test_conv(void) ;
extern void init_simplex (Param_simplex *fun, int nb_pts, Vertex *vtx, int nb_vtx, int nb_param) ;
extern void free_simplex(Param_simplex *fun, Vertex *vtx, int nb_vtx) ;
extern void swap_vtx(Vertex *v1, Vertex *v2) ;
extern void simplex_gen(void (*function) (/* ??? */), Param_simplex f, Vertex *vtx, int nbv) ;
extern double calc_omega(double *prob, int k1, int k2);
extern double calc_mu(double *prob, int k1, int k2);
extern float   ynist(float aminh, float amaxh, long int Nbpas) ;
extern int  ystog(float *ibuf, long int Nbpas, float aminh, float reso, float valeur, float duree) ;
extern void get_gaussian(float s, float *y, int *len);
extern void get_derigaussian(float s, float *y, int *len);
#ifdef __cplusplus 
}
#endif /*__cplusplus*/ 

#endif
