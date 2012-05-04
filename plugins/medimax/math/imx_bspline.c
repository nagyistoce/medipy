/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/********************************************************************************
**/	
/*!	\file:		imx_bspline.c
***
***	project:	Imagix 2.01
***			
***
***	\brief description:    Fichier fonction spline
***	
***	
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***
********************************************************************************/

#include <config.h>
#include	<stdlib.h> 	/* Needed when you use GET_FLOAT	*/

#include "noyau/imagix.h"
#include "math/imx_matrix.h"
#include "math/imx_bspline.h"

/********* VARIABLES GLOBALES *********/
/*base de decomposition des deformation*/
base	BASE;




/*******************************************************************************
********************************************************************************
**************** GESTION DE LA BASE DE FONCTIONS *******************************
********************************************************************************
*******************************************************************************/

/**************************************************************************
**      Bsplined0(pas,f,df)
*/	
/*!      Fonction d'echelle Bspline de degre 0 \n
**	f (E/S) et (E/S) df sont alloues si il ne le sont pas et realloues sinon
**	retourne le nombre d'intervalles occupes pas la fonction
**************************************************************************/
int	Bsplined0(int pas, double **f, double **df)
{
 int i;
 double *ff,*dff;
 	
 ff=*f;dff=*df;	
 if (f==NULL) ff=CALLOC(pas,double); else ff=REALLOC(ff,pas,double);	
 if (df==NULL) dff=CALLOC(pas,double); else dff=REALLOC(dff,pas,double);	
 
 if (ff==NULL || dff==NULL) {printf("ERREUR d'allocation dans Bsplined0\n");exit(1);}

 for (i=0;i<pas;i++) {ff[i]=1.0;dff[i]=0.0;}
 
 *f=ff;*df=dff;
 return(1);
}

/**************************************************************************
**      Bsplined1(pas,f,df)
*/   
/*! 	 Fonction d'echelle Bspline de degre 1 \n
**	f (E/S) et df (E/S) sont alloues si il ne le sont pas et realloues sinon
**	retourne le nombre d'intervalles occupes pas la fonction
**************************************************************************/
int	Bsplined1(int pas, double **f, double **df)
{
 int i;
 double *ff,*dff;
 	
 ff=*f;dff=*df;	
 if (f==NULL) ff=CALLOC(2*pas,double); else ff=REALLOC(ff,2*pas+1,double);	
 if (df==NULL) dff=CALLOC(2*pas,double); else dff=REALLOC(dff,2*pas+1,double);	
 
 if (ff==NULL || dff==NULL) {printf("ERREUR d'allocation dans Bsplined0\n");exit(1);}

 for (i=0;i<pas;i++) {ff[i]=(double)(i)/(double)pas;dff[i]=1.0/(double)pas;}
 for (i=pas;i<2*pas;i++) {ff[i]=(double)(2*pas-i)/(double)pas;dff[i]=-1.0/(double)pas;}
 
 *f=ff;*df=dff;
 return(2);
}

/**************************************************************************
**      Bsplined2(pas,f,df)
*/   
/*! 	 Fonction d'echelle Bspline de degre 1 \n
**	f (E/S) et df (E/S) sont alloues si il ne le sont pas et realloues sinon
**	retourne le nombre d'intervalles occupes pas la fonction
**************************************************************************/
int	Bsplined2(int pas, double **f, double **df)
{
 int i;
 double *ff,*dff;
 	
 ff=*f;dff=*df;	
 if (f==NULL) ff=CALLOC(3*pas,double); else ff=REALLOC(ff,3*pas,double);	
 if (df==NULL) dff=CALLOC(3*pas,double); else dff=REALLOC(dff,3*pas,double);	
 
 if (ff==NULL || dff==NULL) {printf("ERREUR d'allocation dans Bsplined0\n");exit(1);}

 for (i=0;i<pas;i++) 
  {ff[i]=(double)(i*i)/(double)(2*pas*pas);dff[i]=(double)i/(double)(pas*pas);}
 for (i=pas;i<2*pas;i++) 
  {ff[i]=(double)(3*i)/(double)pas-(double)(i*i)/(double)(pas*pas)-3.0/2.0;
   dff[i]=3.0/(double)pas-(double)(2*i)/(double)(pas*pas);}
 for (i=2*pas;i<3*pas;i++) 
  {ff[i]=(double)(i*i)/(double)(2*pas*pas)-(double)(3*i)/(double)pas+9.0/2.0;
   dff[i]=(double)i/(double)(pas*pas)-3.0/(double)pas;}
  
 *f=ff;*df=dff;
 return(3);
}

/**************************************************************************
**      Bsplineod1(pas,f,df)
*/   
/*! 	 Fonction d'echelle Bspline de degre 1 \n
**	f (E/S) et df (E/S) sont alloues si il ne le sont pas et realloues sinon
**	retourne le nombre d'intervalles occupes pas la fonction
**************************************************************************/
int	Bsplineod1(int pas, double **f, double **df)
{
 int i;
 double *ff,*dff;
 int dpas;
 
 dpas=(int)(pas/2);
 	
 ff=*f;dff=*df;	
 if (f==NULL) ff=CALLOC(2*pas,double); else ff=REALLOC(ff,2*pas,double);	
 if (df==NULL) dff=CALLOC(2*pas,double); else dff=REALLOC(dff,2*pas,double);	
 
 if (ff==NULL || dff==NULL) {printf("ERREUR d'allocation dans Bsplineod1\n");exit(1);}

 for (i=0;i<dpas;i++) 
  {ff[i]=-2.0*(double)i/(double)pas;dff[i]=-2.0/(double)pas;}
 for (i=dpas;i<pas;i++) 
  {ff[i]=(double)(4*i)/(double)pas-3.0;
   dff[i]=4.0/(double)pas;}
 for (i=pas;i<pas+dpas;i++) 
  {ff[i]=5.0-(double)(4*i)/(double)pas;
   dff[i]=-4.0/(double)pas;}
 for (i=pas+dpas;i<2*pas;i++) 
  {ff[i]=(double)(2*i)/(double)pas-4.0;
   dff[i]=2.0/(double)pas;}
  
 *f=ff;*df=dff;
 return(2);
}

/**************************************************************************
**      init_base(wdth,hght,scal_func)
**	
**      Initialisation du tableau base en fonction de la taille de
**	l'image, de la resolution et de la fonction d'echelle utilisee
**	la fonction retourne un entier correspondant au nombre de parametres
**	La base est allouee pendant cette initialisation et la resolution de
**	depart est calculee automatiquement en fonction du support de la fonction
**	d'echelle
**************************************************************************/
int	init_base(int wdth, int hght, int (*scal_func) (/* ??? */))
{
 int	*x0,*x1,*y0,*y1,*idx,*idy,tfilt;
 int	nbint,pasx,pasy,nbfx,nbfy,supx,supy;
 int	i,j,l;
 int	resol,nb_func;
 double *fx,*fy,*dfx,*dfy,*px,*py,*filt;
 
 fx=fy=dfx=dfy=NULL;
 /*calcul de la resolution initiale et d'autres parametres*/
 resol=-1;
 do
 {
   resol++;
   /*nombre d'intervalle de decoupage*/
    nbint=(int)(pow(2.0,(double)resol));
   /*pas suivant x et y*/ 
    pasx=wdth/nbint;pasy=hght/nbint;
   /*calcul du nombre de fonctions suivant x et y*/
    supx=(*scal_func)(pasx,&fx,&dfx);nbfx=nbint+1-supx;
    supy=(*scal_func)(pasy,&fy,&dfy);nbfy=nbint+1-supy;
    nb_func=nbfx*nbfy;
 } while (nbfx<=0 || nbfy<=0);
 
 /*calcul du filtre en fonction de la fonction d'echelle*/
  tfilt=2;filt=CALLOC(tfilt,double);filt[0]=filt[1]=1.0;
  if (scal_func==Bsplined1)
    {tfilt=3;filt=CALLOC(tfilt,double);filt[0]=filt[2]=0.5;filt[1]=1.0;}
  if (scal_func==Bsplined2)
    {tfilt=4;filt=CALLOC(tfilt,double);filt[0]=filt[3]=1.0/4.0;filt[1]=filt[2]=3.0/4.0;}
  if (scal_func==Bsplineod1)
    {tfilt=3;filt=CALLOC(tfilt,double);filt[0]=filt[2]=-1.0;filt[1]=1.0;}
 BASE.tfilt=tfilt;BASE.filt=filt;

   
 /*quelques caracteristique de la base*/
 BASE.width=wdth;BASE.height=hght;
 BASE.nb_func=nb_func;BASE.nbfx=nbfx;BASE.nbfy=nbfy;BASE.resol=resol;
 BASE.fx=fx;BASE.fy=fy;BASE.dfx=dfx;BASE.dfy=dfy;
 BASE.scal_func=scal_func;
  
 /*allocation des element*/
 x0=CALLOC(tfilt,int); y0=CALLOC(nb_func,int); 
 x1=CALLOC(nb_func,int);y1=CALLOC(nb_func,int); 
 idx=CALLOC(nb_func,int);idy=CALLOC(nb_func,int); 
 px=CALLOC(nb_func,double);py=CALLOC(nb_func,double); 
 if (x0==NULL || x1==NULL || y0==NULL || y1==NULL || idx==NULL || idy==NULL || px==NULL || py==NULL)
  {printf("ERREUR d'allocation dans init_base\n");exit(1);}
 
  /*remplissage de la base*/
 l=0;
 for (i=0;i<nbfx;i++)
  for (j=0;j<nbfy;j++) 
   {
    x0[l]=i*pasx;x1[l]=i*pasx+supx*pasx;idx[l]=i;px[l]=0.0;
    y0[l]=j*pasy;y1[l]=j*pasy+supy*pasy;idy[l]=j;py[l]=0.0;
    l++;
   }
 BASE.x0=x0;BASE.x1=x1;BASE.y0=y0;BASE.y1=y1;BASE.idx=idx;BASE.idy=idy;
 BASE.px=px;BASE.py=py;
      
 return(2*nb_func);
}

/**************************************************************************
**      init_base_resol(wdth,hght,dpth,scal_func,resolution)
**	
**      Initialisation du tableau base en fonction de la taille de
**	l'image, de la resolution et de la fonction d'echelle utilisee
**	la fonction retourne un entier correspondant au nombre de parametres
**	La base est allouee pendant cette initialisation et la resolution est
**	passee l'utilisateur
**************************************************************************/
int	init_base_resol(int wdth, int hght, int (*scal_func) (/* ??? */), int resolution)
{
 int	*x0,*x1,*y0,*y1,*idx,*idy,tfilt;
 int	nbint,pasx,pasy,nbfx,nbfy,supx,supy;
 int	i,j,l;
 int	resol,nb_func;
 double *fx,*fy,*dfx,*dfy,*px,*py,*filt;
 
 fx=fy=dfx=dfy=NULL;
 
 resol=resolution;
  /*nombre d'intervalle de decoupage*/
   nbint=(int)(pow(2.0,(double)resol));
  /*pas suivant x et y*/ 
    pasx=wdth/nbint;pasy=hght/nbint;
  /*calcul du nombre de fonctions suivant x et y*/
    supx=(*scal_func)(pasx,&fx,&dfx);nbfx=nbint+1-supx;
    supy=(*scal_func)(pasy,&fy,&dfy);nbfy=nbint+1-supy;
    nb_func=nbfx*nbfy;
 
 /*calcul du filtre en fonction de la fonction d'echelle*/
  tfilt=2;filt=CALLOC(tfilt,double);filt[0]=filt[1]=1.0;
  if (scal_func==(int (*)()) Bsplined1)
    {tfilt=3;filt=CALLOC(tfilt,double);filt[0]=filt[2]=0.5;filt[1]=1.0;}
  if (scal_func==(int (*)()) Bsplined2)
    {tfilt=4;filt=CALLOC(tfilt,double);filt[0]=filt[3]=1.0/4.0;filt[1]=filt[2]=3.0/4.0;}
  if (scal_func==(int (*)()) Bsplineod1)
    {tfilt=3;filt=CALLOC(tfilt,double);filt[0]=filt[2]=-1.0;filt[1]=1.0;}
 BASE.tfilt=tfilt;BASE.filt=filt; 

   
 /*quelques caracteristique de la base*/
 BASE.width=wdth;BASE.height=hght;
 BASE.nb_func=nb_func;BASE.nbfx=nbfx;BASE.nbfy=nbfy;BASE.resol=resol;
 BASE.fx=fx;BASE.fy=fy;BASE.dfx=dfx;BASE.dfy=dfy;
 BASE.scal_func=scal_func;
  
 /*allocation des element*/
 x0=CALLOC(nb_func,int); y0=CALLOC(nb_func,int);  
 x1=CALLOC(nb_func,int);y1=CALLOC(nb_func,int); 
 idx=CALLOC(nb_func,int);idy=CALLOC(nb_func,int);  
 px=CALLOC(nb_func,double);py=CALLOC(nb_func,double); 
 if (x0==NULL || x1==NULL || y0==NULL || y1==NULL || idx==NULL || idy==NULL || px==NULL ||
 py==NULL)
  {printf("ERREUR d'allocation dans init_base\n");exit(1);}
 
  /*remplissage de la base*/
 l=0;
 for (i=0;i<nbfx;i++)
  for (j=0;j<nbfy;j++) 
   {
     x0[l]=i*pasx;x1[l]=i*pasx+supx*pasx;idx[l]=i;px[l]=0.0;
     y0[l]=j*pasy;y1[l]=j*pasy+supy*pasy;idy[l]=j;py[l]=0.0;
     l++;
    }
 BASE.x0=x0;BASE.x1=x1;BASE.y0=y0;BASE.y1=y1;BASE.idx=idx;BASE.idy=idy;
 BASE.px=px;BASE.py=py;
      
 return(3*nb_func);
}

/**************************************************************************
**      base_resol_up(param,nb_param)
**	
**      Passage de la base de fonction a la resolution superieure
**	Modifie les parametres en consequence :
**		nouvelle taille en resultat de la fonction
**		nouvelles valeurs dans le tableau param
**************************************************************************/
int	base_resol_up(double *param, int nb_param)
{
 int	i,j,l;
 int	pasx,pasy,wdth,hght;
 int	resol,nbint,tfilt,nb_func,nbfx,nbfy,supx,supy;
 double *fx,*fy,*dfx,*dfy,*px,*py,*filt;
 int	*x0,*x1,*y0,*y1,*idx,*idy;
 double **Mx,**My,**F2D;
 int	i0,j0;
 int	(*scal_func)();
 
  wdth=BASE.width;hght=BASE.height;
  
 /*parametre de la base dans les variables locales*/
  resol=BASE.resol;resol++; /*passage a la resolution superieure*/
  fx=BASE.fx;fy=BASE.fy;dfx=BASE.dfx;dfy=BASE.dfy;
  filt=BASE.filt;tfilt=BASE.tfilt;
  x0=BASE.x0;x1=BASE.x1;y0=BASE.y0;y1=BASE.y1;
  px=BASE.px;py=BASE.py;idx=BASE.idx;idy=BASE.idy;
  nb_func=BASE.nb_func;nbfx=BASE.nbfx;nbfy=BASE.nbfy;
  scal_func=BASE.scal_func;
 /*les parametres du tableau param sont mis dans la base*/
  if (nb_param!=2*nb_func) {printf("ERREUR dans base_resol_up \n");exit(1);} 
  for (l=0;l<nb_func;l++) {px[l]=param[2*l];py[l]=param[2*l+1];}
  
 /*allocation et remplissage des matrices de parametres Mx et My*/
 /*c'est sur ces matrices que l'on appliquera le filtre de passage d'une resolution a une
 superieure*/ 
   /*calcul de la taille des matrices apres zero filling, allocation memoire et mise a zero*/
    nbfx=tfilt+2*nbfx-2;nbfy=tfilt+2*nbfy-2;
    Mx=alloc_dmatrix(nbfx,nbfy);My=alloc_dmatrix(nbfx,nbfy);
    for (i=0;i<nbfx;i++) for (j=0;j<nbfy;j++) Mx[i][j]=My[i][j]=0.0;
   /*allocation et calcul du filtre de convolution 2D*/
    F2D=alloc_dmatrix(tfilt,tfilt);
    for (i=0;i<tfilt;i++) 
     for (j=0;j<tfilt;j++)
      F2D[i][j]=filt[i]*filt[j];
      
   /*calcul de Mx et My par application du filtre de convolution 2D*/
    for (l=0;l<nb_func;l++)
     {
      i0=2*idx[l];j0=2*idy[l];
      for (i=i0;i<i0+tfilt;i++)
       for (j=j0;j<j0+tfilt;j++)
        {
	 Mx[i][j]+=px[l]*F2D[i-i0][j-j0];
	 My[i][j]+=py[l]*F2D[i-i0][j-j0];
	}
     }   
 
 /*mise a jour des differents parametres de la base*/
   /*nombre d'intervalle de decoupage*/
    nbint=(int)(pow(2.0,(double)resol));
   /*pas suivant x et y*/ 
    pasx=wdth/nbint;pasy=hght/nbint;
   /*mise a jour de l'echantillonage de la fonction d'echelle*/
    supx=(*scal_func)(pasx,&fx,&dfx);
    supy=(*scal_func)(pasy,&fy,&dfy);
   /*verification que ca correspond avec ce qui est deja calcule*/
    if (nbfx!=nbint+1-supx || nbfy!=nbint+1-supy) printf("ERREUR 1 dans base_resol_up\n");
    nb_func=nbfx*nbfy;
  /*reallocation des element*/
    x0=REALLOC(x0,nb_func,int);y0=REALLOC(y0,nb_func,int); 
    x1=REALLOC(x1,nb_func,int);y1=REALLOC(y1,nb_func,int); 
    idx=REALLOC(idx,nb_func,int);idy=REALLOC(idy,nb_func,int); 
    px=REALLOC(px,nb_func,double);py=REALLOC(py,nb_func,double); 
    if (x0==NULL || x1==NULL || y0==NULL || y1==NULL || idx==NULL || idy==NULL || px==NULL || py==NULL)
     {printf("ERREUR d'allocation dans init_base\n");exit(1);}
  /*remplissage de la base*/
   l=0;
    for (i=0;i<nbfx;i++)
     for (j=0;j<nbfy;j++) 
      {
       x0[l]=i*pasx;x1[l]=i*pasx+supx*pasx;idx[l]=i;px[l]=Mx[i][j];
       y0[l]=j*pasy;y1[l]=j*pasy+supy*pasy;idy[l]=j;py[l]=My[i][j];
       l++;
      }
   BASE.x0=x0;BASE.x1=x1;BASE.y0=y0;BASE.y1=y1;BASE.idx=idx;BASE.idy=idy;
   BASE.px=px;BASE.py=py;
   BASE.resol=resol;BASE.nb_func=nb_func;BASE.nbfx=nbfx;BASE.nbfy=nbfy;
   BASE.fx=fx;BASE.fy=fy;BASE.dfx=dfx;BASE.dfy=dfy;
  
  /*nouvelles valeurs des parametres dans le tableau param*/
   for (l=0;l<nb_func;l++)
    {param[2*l]=px[l];param[2*l+1]=py[l];}
  
 free_dmatrix(Mx,nbfx,nbfy);free_dmatrix(My,nbfx,nbfy);free_dmatrix(F2D,tfilt,tfilt);	
 return(2*nb_func);
}

/**************************************************************************
**      end_base
**	
**      Terminaison de la base
**	Liberation memoire des different elements
**************************************************************************/
int	end_base(void)
{
 free(BASE.x0);free(BASE.x1);
 free(BASE.y0);free(BASE.y1);
 free(BASE.idx);free(BASE.idy);
 free(BASE.px);free(BASE.py);
 free(BASE.fx);free(BASE.fy);
 free(BASE.dfx);free(BASE.dfy);
 free(BASE.filt);
 return(1);
}
