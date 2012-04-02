/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/**----------------------------------------------------------------------------
***
***		File:           imx_log.h
***
***		Project:		Imagix 2.01
***		
***     Description:
***     Fichier ANSI pour la gestion du fichier LOG
***
***     Copyright (c) 1993, ULP-IPB Strasbourg.
***     All rights are reserved.
***
***     - When changing please issue your name... -
***     Last user action by: on July th 1999
***
***------------------------------------------------------------------------- */

#ifndef __IMX_LOG_H
#define __IMX_LOG_H

#include <stdio.h> 
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/io/imx__atm.h"

extern FILE	*FICHIER_LOG;

typedef struct s_flog {
int indice;
int nb_element;
FILE *ptr_fichier;
} flog, *ptr_flog;

#endif
 
/* AUTOMATICALLY EXTRACTED PROTOTYPES (MB 6/12/2000) */
#ifdef __cplusplus 
extern "C" 
{
#endif /*__cplusplus*/ 
extern int init_log(char *strFile) ;
extern int aff_log(char *s, ...) ;
extern int end_log(void) ;

extern int init_flog(char *strFile) ;
extern int aff_flog(int indice,char *s, ...) ;
extern int end_flog(int indice) ;

extern int ipmpar(int *i) ;
extern double algdiv(double *a,double *b) ;
extern double alngam(double *x) ;
extern double alnrel(double *a) ;
extern double apser(double *a,double *b,double *x,double *eps) ;
extern double basym(double *a,double *b,double *lambda,double *eps) ;
extern double bcorr(double *a0,double *b0) ;
extern double betaln(double *a0,double *b0) ;
extern double bpser(double *a,double *b,double *x,double *eps) ;
extern double brcmp1(int *mu,double *a,double *b,double *x,double *y) ;
extern double brcomp(double *a,double *b,double *x,double *y) ;
extern double bup(double *a,double *b,double *x,double *y,int *n,double *eps) ;
extern void cumchi(double *x,double *df,double *cum,double *ccum) ;
extern void cumf(double *f,double *dfn,double *dfd,double *cum,double *ccum) ;
extern void cumgam(double *x,double *a,double *cum,double *ccum) ;
extern void cumnor(double *arg,double *result,double *ccum) ;
extern void cumpoi(double *s,double *xlam,double *cum,double *ccum) ;
extern void cumt(double *t,double *df,double *cum,double *ccum) ;
extern double devlpl(double a[],int *n,double *x) ;
extern double dinvnr(double *p,double *q) ;
extern double dt1(double *p,double *q,double *df) ;
extern void dstzr(double *zxlo,double *zxhi,double *zabstl,double *zreltl) ;
extern double erf1(double *x) ;
extern double erfc1(int *ind,double *x) ;
extern double esum(int *mu,double *x) ;
extern double exparg(int *l) ;
extern double fpser(double *a,double *b,double *x,double *eps) ;
extern double gam1(double *a) ;
extern double gamln(double *a) ;
extern double gamln1(double *a) ;
extern double Xgamm(double *a) ;
extern void gratio(double *a,double *x,double *ans,double *qans,int *ind) ;
extern double gsumln(double *a,double *b) ;
extern double psi(double *xx) ;
extern double rcomp(double *a,double *x) ;
extern double rexp(double *x) ;
extern double rlog(double *x) ;
extern double rlog1(double *x) ;
extern double spmpar(int *i) ;
extern double stvaln(double *p) ;
extern double fifdint(double a) ;
extern double fifdmax1(double a,double b) ;
extern double fifdmin1(double a,double b) ;
extern double fifdsign(double mag,double sign) ;
extern long fifidint(double a) ;
extern long fifmod(long a,long b) ;
extern void ftnstop(char* msg) ;
#ifdef __cplusplus 
}
#endif /*__cplusplus*/ 
