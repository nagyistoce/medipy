/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/**-----------------------------------------------------------------------
***	
***	file:       imx_matrix.c
***
***	project:    Imagix 1.01 
***			
***
***	description:	Tools.
***	
***	
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***     Last user action by: Mr. ARMSPACH on Aug 17th 1994
***     Last user action by: Mr. NIKOU    on Dec 13th 1995
***
***---------------------------------------------------------------------*/

#ifndef __IMX_MATRIX_H
#define __IMX_MATRIX_H

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_3d.h"

#ifdef __cplusplus
extern "C" {
#endif


//------------- VECTEURS -----------------//
//---double
extern double *alloc_dvector(int n);
extern double *realloc_dvector(double *vector, int n);
extern void free_dvector(double *v, int n);
//---entiers
extern int *alloc_ivector(int n);
extern int *realloc_ivector(int *vector, int n);
extern void free_ivector(int *v, int n);
//---flottants
extern float *alloc_vector(int n);
extern void free_vector(float *v, int n);

//------------- MATRICES 2D  -------------//
//---double
extern double **alloc_dmatrix(int nb_lig, int nb_col);
extern void free_dmatrix(double **m, int co, int li);
//--entiers
extern int **alloc_imatrix(int nb_lig, int nb_col);
extern void free_imatrix(int **m, int co, int li);
//---flottants
extern float **alloc_matrix(int nb_lig, int nb_col);
extern void free_matrix(float **m, int co, int li);

extern double **multiply_matrix(double **a, int l1, int c1, double **b, int l2, int c2);
extern double **add_matrix(double **a, int l, int c, double **b);
extern double **sub_matrix(double **a, int l, int c, double **b);
extern double *matrix_multiply_vector(double **a, int l, int c, double *b);
extern double  interpolant(double x, int type);
extern double **matrix_transp(double **a, int li, int co);
extern double **matrix_inversion(double **a, int n);
extern double **matrix_inversion_sym(double **a, int n);
extern int  update_matrix_BFGS(double **Hk1, double **Hk2, double *paramk1, double *paramk2, double *gradk1, double *gradk2, int nb_param);
extern double **eigen_values_jacobi(double **a, int n, int sort);
extern void back_subst (double **a, int n, int *index, double *b);

//---decomposition de Choleski d'une matrice symetrique definie positive---//
extern int choldc(double **a, double *diag, int n);
//---decomposition LU d'une matrice---//
extern void lu_decomposition(double **a, int n, int *index, double *d);

extern double Determinant(double **a,int n);

//------------- MATRICES 3D  -------------//
//---double
extern double 	***alloc_dmatrix_3d(int nb_lig, int nb_col, int nb_dep);
extern int	free_dmatrix_3d(double ***m);
//--entiers
extern int ***alloc_imatrix_3d(int nb_lig, int nb_col, int nb_dep);
extern int free_imatrix_3d(int ***m);
//---flottants
extern float ***alloc_matrix_3d(int nb_lig, int nb_col, int nb_dep);
extern int free_matrix_3d(float ***m);

extern void arrange_vector(double *arr, int n, int *index);
extern double *vecteur_moyen(double **vectors, int t, int n);
extern double **covariance_of_small_training_set(double **vectors, double *vm, int t, int n);
extern double *add_vectors(double *a, double *b, int n);
double **cons_multiply_matrix(double c, double **a, int li, int co);


int imx_convert_dmatrix3d_grphic3d(double ***mat,grphic3d* im);

#ifdef WIN32
extern int truncate( const char *path, int length );
#endif

#ifdef __cplusplus
}
#endif

#endif  // #ifndef __IMX_MATRIX_H
