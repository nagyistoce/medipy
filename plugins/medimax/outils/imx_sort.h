/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*!*************************************************************************
**
**	\file		imx_sort.c
**
**	project:	Imagix 1.01
**
**
**	\brief description:	Outils de tri rapide
**
**
**	Copyright (c) 1993, ULP-IPB Strasbourg.
**	All rights are reserved.
**
**     Last user action by: Mr. ARMSPACH on Aug 17th 1994
**     Last user action by: Mr. NIKOU    on Dec 13th 1995
**
***************************************************************************/
#ifndef _IMX_SORT_H
#define _IMX_SORT_H


extern double quick_select(int k, double *arr, int n);   /*** arr[0....n-1] ****/
extern int quick_select_integer(int k, int *arr, int n);   /*** arr[0....n-1] ****/
extern void indexed_quicksort(double *arr, int n, int *index);
extern void d_indexed_quicksort(double *arr, int n, int *index);
extern void simultaneous_quicksort(double *arr, int n, double **matrix, int m);
extern int tri_rapide(long int *tab, int bas, int haut) ;
extern int tri_rapide_double(double *tab, int bas, int haut) ;

#endif
