/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*******************************************************************************/   
/*!    \file:       mtch_3d.c
***
***     project:    Imagix 1.01
***         
***
***     \brief description:    Fichier pourle matching des images 3D
*** 
*** 
*** Copyright (c) 1993, ULP-IPB Strasbourg.
*** All rights are reserved.
***
***
********************************************************************************/
#include <config.h>
#include <stdio.h> 
#include <stdlib.h>     /* Needed when you use GET_FLOAT    */
#include <math.h>   /* Needed when you use math function    */
#include <limits.h>
#include <float.h>
#include <time.h>

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_3d.h"
#include "noyau/imx_lang.h"
#include "noyau/io/imx_head.h"
#include "recalage/chps_3d.h"
#include "recalage/distance_3d.h"
#include "recalage/mtch_3d.h"
#include "math/imx_matrix.h"
#include "outils/imx_sort.h"
#include "noyau/gui/imx_picture3d.h"
#include "noyau/mani_3d.h"
#include "traitement/norm_3d.h"
#include "recalage/mtch_robust_3d.h"
#include "outils/imx_misc.h"
#include "recalage/minimisation_3d.h"
#include "segmentation/otsu_3d.h"

#include "math/imx_bspline.h"   /* pour les fonctions d'interpolation */

#include "noyau/io/imx_log.h"
#include "math/oper_3d.h"    
#include "traitement/trai_3d.h"
#include "recalage/erreur_recalage.h"
#include "recalage/multiresolution_3d.h"

/*static int  NBIM=0;*/
base3d  BASE3D;

/* variable pour restreindre l'optimisation suivant un certain nombre de parametre*/
double *_mtch_3d_hybride_param = NULL;
    
/* seuil a partir duquel on effectue le recalage en local. Apparemment, on a de
meilleurs resultats pour un seuil de 3, mais le recalage est nettement plus lent. */
#define RESOL_LOCAL /*3*/4

//definition du diametre du cerveau en mm
#ifndef BRAIN_RADIUS
#define BRAIN_RADIUS 100.0
#endif

/* DOXYGEN SECTION */
/*! \defgroup Recalage
    \brief module de recalage d'images
*/

/*! \defgroup InterpolationFonction
    \ingroup Recalage
    \brief fonction d'interpolation
*/

/*! \defgroup DistanceFonction
    \ingroup Recalage
    \brief fonction de calcul de distance
*/

/*! \defgroup MinimisationFonction
    \ingroup Recalage
    \brief fonction de minimisation
*/


/*! \defgroup GradientFonction
    \ingroup Recalage
    \brief Calcul de gradient
*/

/*! \defgroup TransfoFonction
    \ingroup Recalage
    \brief Fonction de transformation (params - champ)
*/
/* END DOXYGEN SECTION */

/*********************************************************************/
/*     cr_transf_anisotrope_3d(m)                                    */
/*                                                                   */
/*!     Creation d'une transformation anisotrope a partir
**      d'une transfo isotrope (utilise dx,dy,dz des images pour obtenir la transfo
**      de imreca dans l'espace de imref).
**      \param transfOriginal : la transformation en mm (rigid,rigid_zoom,affine)
**      \param imref   : image reference (dx dy dz d'arrivee)
**      \param imreca  : image a transformer (dx dy dz de depart)
**      \retval transf3d * : une transformation qui amene de l'espace de depart vers l'espace d'arrivee
*********************************************************************/
transf3d *cr_transf_anisotrope_3d(transf3d *transfOriginal,grphic3d *imref, grphic3d *imreca)
{
  transf3d *transfS1=NULL,*transfS2=NULL,*transfTmp=NULL;
  transf3d *transfAnisotrope=NULL;
  double dx,dy,dz;
  
  if (!imref || !imreca)
    return NULL;                                  

  transfAnisotrope = cr_transf3d(imref->width,imref->height,imref->depth,NULL);
  transfAnisotrope->param = CALLOC(15,double);
  transfS1 = cr_transf3d(0,0,0,NULL);  transfS1->param = CALLOC(15,double);
  transfS1->nb_param = 15; transfS1->typetrans = AFFINE3D;
  transfS2 = cr_transf3d(0,0,0,NULL);  transfS2->param = CALLOC(15,double);
  transfS2->nb_param = 15; transfS2->typetrans = AFFINE3D;
  transfTmp = cr_transf3d(0,0,0,NULL); transfTmp->param = CALLOC(15,double);

  
  dx = imreca->dx ; dy = imreca->dy ; dz = imreca->dz;
  transfS1->param[0] =1.0/dx; transfS1->param[1] =0.0   ;transfS1->param[2] =0.0;    //matrix
  transfS1->param[3] =0.0;    transfS1->param[4] =1.0/dy;transfS1->param[5] =0.0;    //matrix
  transfS1->param[6] =0.0;    transfS1->param[7] =0.0   ;transfS1->param[8] =1.0/dz; //matrix
  transfS1->param[9] =0.0;    transfS1->param[10]=0.0   ;transfS1->param[11]=0.0;    //translation
  transfS1->param[12]=0.0;    transfS1->param[13]=0.0   ;transfS1->param[14]=0.0;    //centre

  comb_transf_3d(transfS1,transfOriginal,transfTmp,0);

  dx = imref->dx ; dy = imref->dy ; dz = imref->dz;
  transfS2->param[0] =dx;     transfS2->param[1] =0.0   ;transfS2->param[2] =0.0;    //matrix
  transfS2->param[3] =0.0;    transfS2->param[4] =dy    ;transfS2->param[5] =0.0;    //matrix
  transfS2->param[6] =0.0;    transfS2->param[7] =0.0   ;transfS2->param[8] =dz;     //matrix
  transfS2->param[9] =0.0;    transfS2->param[10]=0.0   ;transfS2->param[11]=0.0;    //translation
  transfS2->param[12]=0.0;    transfS2->param[13]=0.0   ;transfS2->param[14]=0.0;    //centre
  
  comb_transf_3d(transfTmp,transfS2,transfAnisotrope,0);


  
  if (transfS1)  free_transf3d(transfS1);  
  if (transfS2)  free_transf3d(transfS2);
  if (transfTmp) free_transf3d(transfTmp);

  return transfAnisotrope;  
}

 
/*******************************************************************************
********************************************************************************
************************* TRANSFORMATIONS  *************************************
********************************************************************************
*******************************************************************************/
/*! \ingroup TransfoFonction @{ */
/*******************************************************************************
**     rigid_to_field_3d(nb_param,param,champ)                          
*/                                                                   
/*!     remplit un champ dense 3D pour une transformation rigide                           
**  applique autour du centre (xc,yc) et avec une translation (tx,ty)
**      Les parametres dans le tableau param, sont ordonnes ainsi
**  tetax tetay tetaz tx ty tz xc yc zc    
**  \param nb_param : inutilise                            
**  \param param  : tableau de param
**  \param champ : le champ (E/S)
**  \param imref : utile pour creer un champ d'une transfo anisotrope (espace imreca -> espace imref)
**  \param imref : utile pour creer un champ d'une transfo anisotrope (espace imreca -> espace imref)
**  \note : si imref et imreca == NULL : la champ est genere tel quel a partir de param.
**  \retval : 1                            
*******************************************************************************/
int rigid_to_field_3d(int nb_param, const double *param, field3d *champ, grphic3d *imref, grphic3d *imreca)
{   
  int i,j,k,wdth,hght,dpth;
  double x,y,z;
  double tetax,tetay,tetaz,tx,ty,tz,xc,yc,zc;
  double cxx,sxx,cyy,syy,czz,szz;
  double a11,a12,a13,a21,a22,a23,a31,a32,a33;
  double tpx11,tpx21,tpx31,tpy12,tpy22,tpy32;
  vector3d ***data;
  transf3d* transfOriginal   = NULL;
  transf3d* transfAnisotrope = NULL;

  transfOriginal = cr_transf3d(0,0,0,NULL);
  transfOriginal->param = CALLOC(15,double);
  transfOriginal->typetrans = RIGID3D;
  transfOriginal->nb_param = nb_param;
  for (i=0;i<nb_param;i++)
    transfOriginal->param[i] = param[i];


  if (imref && imreca)
    transfAnisotrope = cr_transf_anisotrope_3d(transfOriginal,imref,imreca);  
  
  wdth=champ->width;hght=champ->height;dpth=champ->depth;
  data=champ->raw;  

  if (transfAnisotrope && transfAnisotrope->typetrans==AFFINE3D)
  {
    a11=transfAnisotrope->param[0];a12=transfAnisotrope->param[1];a13=transfAnisotrope->param[2];
    a21=transfAnisotrope->param[3];a22=transfAnisotrope->param[4];a23=transfAnisotrope->param[5];
    a31=transfAnisotrope->param[6];a32=transfAnisotrope->param[7];a33=transfAnisotrope->param[8];
    tx=transfAnisotrope->param[9];ty=transfAnisotrope->param[10];tz=transfAnisotrope->param[11];
    xc=transfAnisotrope->param[12];yc=transfAnisotrope->param[13];zc=transfAnisotrope->param[14];
  }
  else
  {
    tetax=param[0];tetay=param[1];tetaz=param[2];
    tx=param[3];ty=param[4];tz=param[5];
    xc=param[6];yc=param[7];zc=param[8];
    cxx=cos(tetax);sxx=sin(tetax);cyy=cos(tetay);syy=sin(tetay);czz=cos(tetaz);szz=sin(tetaz);
    a11=cyy*czz;a12=-cyy*szz;a13=-syy;
    a21=cxx*szz-sxx*syy*czz;
    a22=cxx*czz+sxx*syy*szz;
    a23=-sxx*cyy;
    a31=sxx*szz+cxx*syy*czz;
    a32=sxx*czz-cxx*syy*szz;
    a33=cxx*cyy;  
  }
  
  a11=a11-1.0;a22=a22-1.0;a33=a33-1.0;
  for (i=0;i<wdth;i++)
   {
    x=i-xc;tpx11=a11*x;tpx21=a21*x;tpx31=a31*x;
    for (j=0;j<hght;j++)
    {
     y=j-yc;tpy12=a12*y;tpy22=a22*y;tpy32=a32*y;
     for (k=0;k<dpth;k++)
      {
       z=k-zc;
       data[i][j][k].x=(float)(tpx11+tpy12+a13*z+tx);
       data[i][j][k].y=(float)(tpx21+tpy22+a23*z+ty);
       data[i][j][k].z=(float)(tpx31+tpy32+a33*z+tz);
      }
    }
   }
  if (transfOriginal) free_transf3d(transfOriginal);
  if (transfAnisotrope) free_transf3d(transfAnisotrope);
  
  return(1);
}

/*! @} */

/*! \ingroup TransfoFonction @{ */
/*******************************************************************************
**     rigidz_to_field_3d(nb_param,param,champ)                          
*/                                                                   
/*!     rempli un champ dense 3D pour une transformation rigide avec zoom                        
**  applique autour du centre (xc,yc,zc) et avec une translation (tx,ty,tz)
**      Les parametres dans le tableau param, sont ordonnes ainsi
**  tetax tetay tetaz sx sy sz tx ty tz xc yc zc                                   
**  \param nb_param : inutilise                            
**  \param param  : tableau de param
**  \param champ : le champ (E/S)
**  \param imref : utile pour creer un champ d'une transfo anisotrope (espace imreca -> espace imref)
**  \param imref : utile pour creer un champ d'une transfo anisotrope (espace imreca -> espace imref)
**  \note : si imref et imreca == NULL : la champ est genere tel quel a partir de param.
**  \retval : 1                            
*******************************************************************************/
int rigidz_to_field_3d(int nb_param, const double *param, field3d *champ, grphic3d *imref, grphic3d *imreca)
{   
  int i,j,k,wdth,hght,dpth;
  double x,y,z;
  double tetax,tetay,tetaz,tx,ty,tz,xc,yc,zc,sx,sy,sz;
  double cxx,sxx,cyy,syy,czz,szz;
  double a11,a12,a13,a21,a22,a23,a31,a32,a33;
  double tpx11,tpx21,tpx31,tpy12,tpy22,tpy32;
  vector3d ***data;
  transf3d* transfOriginal   = NULL;
  transf3d* transfAnisotrope = NULL;

  transfOriginal = cr_transf3d(0,0,0,NULL);
  transfOriginal->param = CALLOC(15,double);
  transfOriginal->typetrans = RIGIDZOOM3D;
  transfOriginal->nb_param = nb_param;
  for (i=0;i<nb_param;i++)
  {
    transfOriginal->param[i] = param[i];
  }  

  if (imref && imreca)
    transfAnisotrope = cr_transf_anisotrope_3d(transfOriginal,imref,imreca);

  wdth=champ->width;hght=champ->height;dpth=champ->depth;
  data=champ->raw;  

  if (transfAnisotrope && transfAnisotrope->typetrans==AFFINE3D)
  {
    a11=transfAnisotrope->param[0];a12=transfAnisotrope->param[1];a13=transfAnisotrope->param[2];
    a21=transfAnisotrope->param[3];a22=transfAnisotrope->param[4];a23=transfAnisotrope->param[5];
    a31=transfAnisotrope->param[6];a32=transfAnisotrope->param[7];a33=transfAnisotrope->param[8];
    tx=transfAnisotrope->param[9];ty=transfAnisotrope->param[10];tz=transfAnisotrope->param[11];
    xc=transfAnisotrope->param[12];yc=transfAnisotrope->param[13];zc=transfAnisotrope->param[14];
  }
  else
  {  
    tetax=param[0];tetay=param[1];tetaz=param[2];
    sx=param[3];sy=param[4];sz=param[5];
    tx=param[6];ty=param[7];tz=param[8];
    xc=param[9];yc=param[10];zc=param[11];
    
    cxx=cos(tetax);sxx=sin(tetax);cyy=cos(tetay);syy=sin(tetay);czz=cos(tetaz);szz=sin(tetaz);
    a11=sx*cyy*czz;a12=-sy*cyy*szz;a13=-sz*syy;
    a21=sx*(cxx*szz-sxx*syy*czz);
    a22=sy*(cxx*czz+sxx*syy*szz);
    a23=-sz*sxx*cyy;
    a31=sx*(sxx*szz+cxx*syy*czz);
    a32=sy*(sxx*czz-cxx*syy*szz);
    a33=sz*cxx*cyy;
  }
    
  a11=a11-1.0;a22=a22-1.0;a33=a33-1.0;  
  for (i=0;i<wdth;i++)
   {
    x=(double)i-xc;tpx11=a11*x;tpx21=a21*x;tpx31=a31*x;
    for (j=0;j<hght;j++)
    {
     y=(double)j-yc;tpy12=a12*y;tpy22=a22*y;tpy32=a32*y;
     for (k=0;k<dpth;k++)
      {
       z=(double)k-zc;
       data[i][j][k].x=(float)(tpx11+tpy12+a13*z+tx);
       data[i][j][k].y=(float)(tpx21+tpy22+a23*z+ty);
       data[i][j][k].z=(float)(tpx31+tpy32+a33*z+tz);
      }
    }
   }
  if (transfAnisotrope)
      free_transf3d(transfAnisotrope);
  if (transfOriginal)
      free_transf3d(transfOriginal);

  return(1);
 }
/*! @} */

/*! \ingroup TransfoFonction @{ */
/*******************************************************************************
**     affine_to_field_3d(nb_param,param,champ)                          
*/                                                                   
/*!     rempli un champ dense 3D pour une transformation affine                           
**  applique autour du centre (xc,yc,zc) et avec une translation (tx,ty,tz)
**      Les parametres dans le tableau param, sont ordonnes ainsi
**  a11 a12 a13 a21 a22 a23 a31 a32 a33 tx ty tz xc yc zc                              
**  \param nb_param : inutilise                            
**  \param param  : tableau de param
**  \param champ : le champ (E/S)
**  \retval : 1                            
*******************************************************************************/
int affine_to_field_3d(int nb_param, const double *param, field3d *champ, grphic3d *imref, grphic3d *imreca)
{   
  int i,j,k,wdth,hght,dpth;
  double x,y,z;
  double /*tetax,tetay,tetaz,*/tx,ty,tz,xc,yc,zc/*,sx,sy,sz*/;
  double a11,a12,a13,a21,a22,a23,a31,a32,a33;
  double tpx11,tpx21,tpx31,tpy12,tpy22,tpy32;
  vector3d ***data;
  transf3d* transfOriginal   = NULL;
  transf3d* transfAnisotrope = NULL;

  transfOriginal = cr_transf3d(0,0,0,NULL);
  transfOriginal->param = CALLOC(15,double);
  transfOriginal->typetrans =AFFINE3D;
  transfOriginal->nb_param = nb_param;
  for (i=0;i<nb_param;i++)
    transfOriginal->param[i] = param[i];


  if (imref && imreca)
    transfAnisotrope = cr_transf_anisotrope_3d(transfOriginal,imref,imreca);

  wdth=champ->width;hght=champ->height;dpth=champ->depth;
  data=champ->raw;  
    
  if (transfAnisotrope && transfAnisotrope->typetrans==AFFINE3D)
  {
    a11=transfAnisotrope->param[0];a12=transfAnisotrope->param[1];a13=transfAnisotrope->param[2];
    a21=transfAnisotrope->param[3];a22=transfAnisotrope->param[4];a23=transfAnisotrope->param[5];
    a31=transfAnisotrope->param[6];a32=transfAnisotrope->param[7];a33=transfAnisotrope->param[8];
    tx=transfAnisotrope->param[9];ty=transfAnisotrope->param[10];tz=transfAnisotrope->param[11];
    xc=transfAnisotrope->param[12];yc=transfAnisotrope->param[13];zc=transfAnisotrope->param[14];
  }
  else
  {
    a11=param[0];a12=param[1];a13=param[2]; 
    a21=param[3];a22=param[4];a23=param[5];
    a31=param[6];a32=param[7];a33=param[8];
    tx=param[9];ty=param[10];tz=param[11];
    xc=param[12];yc=param[13];zc=param[14]; 
  }
      
  a11=a11-1.0;a22=a22-1.0;a33=a33-1.0;  
  for (i=0;i<wdth;i++)
   {
    x=(double)i-xc;tpx11=a11*x;tpx21=a21*x;tpx31=a31*x;
    for (j=0;j<hght;j++)
    {
     y=(double)j-yc;tpy12=a12*y;tpy22=a22*y;tpy32=a32*y;
     for (k=0;k<dpth;k++)
      {
       z=(double)k-zc;
       data[i][j][k].x=(float)(tpx11+tpy12+a13*z+tx);
       data[i][j][k].y=(float)(tpx21+tpy22+a23*z+ty);
       data[i][j][k].z=(float)(tpx31+tpy32+a33*z+tz);
      }
    }
   }
  free_transf3d(transfAnisotrope);
  free_transf3d(transfOriginal);
  return(1);
 }
/*! @} */

/*! \ingroup TransfoFonction @{ */
/*******************************************************************************
**     rigid_global_zoom_to_field_3d(nb_param,param,champ)
*/
/*!     rempli un champ dense 3D pour une transformation rigide avec zoom global isotrope
**  applique autour du centre (xc,yc,zc) et avec une translation (tx,ty,tz)
**      Les parametres dans le tableau param, sont ordonnes ainsi
**  tetax tetay tetaz sx sy sz tx ty tz xc yc zc
**  \param nb_param : inutilise
**  \param param  : tableau de param
**  \param champ : le champ (E/S)
**  \param imref : utile pour creer un champ d'une transfo anisotrope (espace imreca -> espace imref)
**  \param imref : utile pour creer un champ d'une transfo anisotrope (espace imreca -> espace imref)
**  \note : si imref et imreca == NULL : la champ est genere tel quel a partir de param.
**  \retval : 1
*******************************************************************************/
int rigid_global_zoom_to_field_3d(int nb_param, const double *param, field3d *champ,grphic3d *imref, grphic3d *imreca)
{
 double *param_zoom;
 int nb_param_zoom=12;
 int i, res;

 //on met tout comme il faut en zoom anisotrope
 param_zoom=CALLOC(nb_param_zoom,double);
 //3 premiers param=angles d'euler
 for (i=0;i<3;i++) param_zoom[i]=param[i];
 //troisieme=scale global
 for (i=3;i<6;i++) param_zoom[i]=param[3];
 //ensuite translations et centre
 for (i=6;i<12;i++) param_zoom[i]=param[i-2];
 //on utilise ce qui est deja fait
 res=rigidz_to_field_3d(nb_param_zoom, param_zoom, champ, imref, imreca);

 FREE(param_zoom);

 return res;
}
/*! @} */

/*! \ingroup TransfoFonction @{ */
/*******************************************************************************
**     affine_decouple_to_field_3d(nb_param,param,champ)
*/
/*!     rempli un champ dense 3D pour une transformation affine ou les parametres sont decouples
**  applique autour du centre (xc,yc,zc) et avec une translation (tx,ty,tz)
**      Les parametres dans le tableau param, sont ordonnes ainsi
**  tetax tetay tetaz sx sy sz tx ty tz xc yc zc
**  \param nb_param : inutilise
**  \param param  : tableau de param

**  \param champ : le champ (E/S)
**  \param imref : utile pour creer un champ d'une transfo anisotrope (espace imreca -> espace imref)
**  \param imref : utile pour creer un champ d'une transfo anisotrope (espace imreca -> espace imref)
**  \note : si imref et imreca == NULL : la champ est genere tel quel a partir de param.
**  \retval : 1
*******************************************************************************/
int affine_decouple_to_field_3d(int nb_param, const double *param, field3d *champ,grphic3d *imref, grphic3d *imreca)
{
 double matrice_affine[15];
 int i, res=0;

 //on utilise ce qui existe deja
 for (i=0;i<nb_param;i++) matrice_affine[i]=param[i];
 affine_decouple_to_affine_3d(matrice_affine);
 res=affine_to_field_3d(nb_param, matrice_affine, champ, imref, imreca);

 return res;
}
/*! @} */

/*! \ingroup TransfoFonction @{ */
/*******************************************************************************
**     affine_to_field_3d(nb_param,param,champ)                          
*/                                                                   
/*!     rempli un champ dense 3D pour une transformation affine                           
**  applique autour du centre (xc,yc,zc) et avec une translation (tx,ty,tz)
**      Les parametres dans le tableau param, sont ordonnes ainsi
**  a11 a12 a13 a21 a22 a23 a31 a32 a33 tx ty tz xc yc zc                              
**  \param nb_param : inutilise                            
**  \param param  : tableau de param
**  \param champ : le champ (E/S)
**  \retval : 1                            
*******************************************************************************/
/*
int affinesscg_to_field_3d(int nb_param, const double *param, field3d *champ, grphic3d *imref, grphic3d *imreca)
{   
  int i,j,k,wdth,hght,dpth;
  double tx,ty,tz;
  double a11,a12,a13,a21,a22,a23,a31,a32,a33;
  double tpx11,tpx21,tpx31,tpy12,tpy22,tpy32;
  vector3d ***data;

  wdth=champ->width;hght=champ->height;dpth=champ->depth;
  data=champ->raw;  
    
  a11=param[0];a12=param[1];a13=param[2]; 
  a21=param[3];a22=param[4];a23=param[5];
  a31=param[6];a32=param[7];a33=param[8];
  tx=param[9];ty=param[10];tz=param[11];
    
  a11=a11-1.0;a22=a22-1.0;a33=a33-1.0;  
  for (i=0;i<wdth;i++)
   {
    tpx11=a11*i;tpx21=a21*i;tpx31=a31*i;
    for (j=0;j<hght;j++)
    {
     tpy12=a12*j;tpy22=a22*j;tpy32=a32*j;
     for (k=0;k<dpth;k++)
      {
       data[i][j][k].x=(float)(tpx11+tpy12+a13*k+tx);
       data[i][j][k].y=(float)(tpx21+tpy22+a23*k+ty);
       data[i][j][k].z=(float)(tpx31+tpy32+a33*k+tz);
      }
    }
   }
  
  return(1);
 }
*/
/*! @} */

/*! \ingroup TransfoFonction @{ */
/*******************************************************************************
**     base_to_field_3d(nb_param,param,champ)                          
*/
/*!                                                                   
**     transformation liee a la base de fonctions ondelettes
**     \param nb_param : nbre d'element dans param
**     \param param : tableau de parametres
**     \param champ : (E/S)
**     \param imref,imreca : pas utilise pour le moment (a voir pour les transfos anisotropiques)
**     \retval : 1
*******************************************************************************/
int base_to_field_3d(int nb_param, const double *param, field3d *champ, grphic3d *imref, grphic3d *imreca)
{   
  int i,j,k,l;
  int wdth,hght,dpth;
  int *x0,*x1,*y0,*y1,*z0,*z1/*,*idx,*idy,*idz*/;
  int x00,x11,y00,y11,z00,z11;
  double *fx,*fy,*fz,px,py,pz,f;
  double dxreca,dyreca,dzreca;
  double dxref,dyref,dzref;
    
  vector3d ***data;
 
  wdth=champ->width;hght=champ->height;dpth=champ->depth;
  data=champ->raw;
  /*mise a zero du champ*/
  for (i=0;i<wdth;i++)
   for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
     {data[i][j][k].x=data[i][j][k].y=data[i][j][k].z=0.0;}
    
  x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
  fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;
 
  if (imreca)
  {
    dxreca = imreca->dx;
    dyreca = imreca->dy;
    dzreca = imreca->dz;
  }
  else
  {
    dxreca = 1.0;
    dyreca = 1.0;
    dzreca = 1.0;
  }
    if (imref)
    {
    if (wdth != (int)imref->width || hght != (int)imref->height || dpth != (int)imref->depth)
      return 0;
        dxref = imref->dx;    
        dyref = imref->dy;    
        dzref = imref->dz;    
    }
    else
    {
           dxref = 1.0;
           dyref = 1.0;
           dzref = 1.0;
    }

if ((dxref!=dxreca)||(dyref!=dyreca)||(dzref!=dzreca))
        {   
        PUT_ERROR("ERREUR dans base_to_field_3d : le champ est l'image n'ont pas les meme dx,dy,dz\n");
       printf("ERREUR dans base_to_field_3d : le champ est l'image n'ont pas les meme dx,dy,dz\n");
         return(-1);
    }
 
 for (l=0;l<nb_param/3;l++)
        {
        x00=x0[l];x11=x1[l];y00=y0[l];y11=y1[l];z00=z0[l];z11=z1[l];
        px=param[3*l]/dxreca;py=param[3*l+1]/dyreca;pz=param[3*l+2]/dzreca;
        for (i=x00;i<x11;i++)
            for (j=y00;j<y11;j++)
                for (k=z00;k<z11;k++) 
            {
            f=fx[i-x00]*fy[j-y00]*fz[k-z00];

            data[i][j][k].x=(float)(data[i][j][k].x+px*f);
            data[i][j][k].y=(float)(data[i][j][k].y+py*f);
            data[i][j][k].z=(float)(data[i][j][k].z+pz*f);       
                }
        }
    
    
  
 
  return(1);
}


/*! @} */
/*! \ingroup TransfoFonction @{ */
/*******************************************************************************
**     champ_to_field_3d(nb_param,param,champ)
*/
/*!
**     transforme un champ millimetrique en champ voxelique
**
**     \retval : 1
*******************************************************************************/
int champ_to_field_3d(transf3d *transfo, field3d *champ, grphic3d *imref, grphic3d *imreca)
{
  int i,j,k,l,wdth,hght,dpth;
  vector3d ***data,*ch;
  float dxreca,dyreca,dzreca;
  float dxref,dyref,dzref;

  wdth = transfo->width; hght = transfo->height; dpth = transfo->depth;

  if (imreca)
  {
    dxreca = imreca->dx; dyreca = imreca->dy; dzreca = imreca->dz;
  }
  else
        {
    dxreca = 1.0; dyreca = 1.0; dzreca = 1.0;
    }

    if (imref)
    {
    if (wdth != (int)imref->width || hght != (int)imref->height || dpth != (int)imref->depth)
      return 0;
        dxref = imref->dx; dyref = imref->dy; dzref = imref->dz;
    }
    else
    {
        if (transfo->dx && transfo->dy && transfo->dz && imreca != NULL)
        {
           dxref = transfo->dx; dyref = transfo->dy; dzref = transfo->dz;
        }
        else
        {
           dxref = 1.0; dyref = 1.0; dzref = 1.0;
        }
    }


  ch = transfo->ch;
  data = champ->raw;
  l=0;
  for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
  for (k=0;k<dpth;k++)
  {
    data[i][j][k].x = (float)((i*dxref+ch[l].x)/dxreca-i);
    data[i][j][k].y = (float)((j*dyref+ch[l].y)/dyreca-j);
    data[i][j][k].z = (float)((k*dzref+ch[l].z)/dzreca-k);
    l++;
  }
  return(1);
}
/*! @} */

/*******************************************************************************
********************************************************************************
*************** PASSAGE TRANSFORMATIONS  ***************************************
********************************************************************************
*******************************************************************************/

/*******************************************************************************
**     gravite_to_rigid_3d(param)                  
*/
/*!                                                                   
**      met a jour les valeurs du tableau param pour le passage de la mise en
**  commun des centre de gravite vers une transfo rigide\n
**  retourne le nouveau nombre de parametres
**  \param param : le tableau de param 
**  \retval : 9 (nbre de param)
*******************************************************************************/
int gravite_to_rigid_3d(double *param)
{
 double tx,ty,tz,xg,yg,zg;
 
 tx=param[0];ty=param[1];tz=param[2];
 xg=param[3];yg=param[4];zg=param[5];
 
 param[0]=param[1]=param[2]=0.0;
 param[3]=tx;param[4]=ty;param[5]=tz;
 param[6]=xg;param[7]=yg;param[8]=zg;
 
 return(9); 
}

/*******************************************************************************
**     rigid_to_rigidz_3d(param)                  
*/
/*!                                                                   
**      met a jour les valeurs du tableau param pour le passage d'une transfo
**  rigide vers une transfo rigide+zoom\n
**  retourne le nouveau nombre de parametres
**  \param param : le tableau de param 
**  \retval 12 (nbre de param)
*******************************************************************************/
int rigid_to_rigidz_3d(double *param)
{
 double tetax,tetay,tetaz,tx,ty,tz,xg,yg,zg;
 
 tetax=param[0];tetay=param[1];tetaz=param[2];
 tx=param[3];ty=param[4];tz=param[5];
 xg=param[6];yg=param[7];zg=param[8];
 
 param[0]=tetax;param[1]=tetay;param[2]=tetaz;
 param[3]=param[4]=param[5]=1.0;
 param[6]=tx;param[7]=ty;param[8]=tz;
 param[9]=xg;param[10]=yg;param[11]=zg;
 
 return(12); 
}

/*******************************************************************************
**     rigid_to_rigid_global_zoom_3d(param)
*/
/*!
**      met a jour les valeurs du tableau param pour le passage d'une transfo
**  rigide vers une transfo rigide+zoom global isotrope\n
**  retourne le nouveau nombre de parametres
**  \param param : le tableau de param
**  \retval 12 (nbre de param)
*******************************************************************************/
int rigid_to_rigid_global_zoom_3d(double *param)
{
 double tetax,tetay,tetaz,tx,ty,tz,xg,yg,zg;

 tetax=param[0];tetay=param[1];tetaz=param[2];
 tx=param[3];ty=param[4];tz=param[5];
 xg=param[6];yg=param[7];zg=param[8];

 param[0]=tetax;param[1]=tetay;param[2]=tetaz;
 param[3]=1.0;
 param[4]=tx;param[5]=ty;param[6]=tz;
 param[7]=xg;param[8]=yg;param[9]=zg;

 return(10);
}

/*******************************************************************************
**     rigid_global_zoom_to_rigidz_3d(param)
*/
/*!
**      met a jour les valeurs du tableau param pour le passage d'une transfo
**  rigide+zoom global isotrope vers une transfo rigide + zoom anisotrope \n
**  retourne le nouveau nombre de parametres
**  \param param : le tableau de param
**  \retval 12 (nbre de param)
*******************************************************************************/
int rigid_global_zoom_to_rigidz_3d(double *param)
{
 double tetax,tetay,tetaz,tx,ty,tz,xg,yg,zg,scale;

 tetax=param[0];tetay=param[1];tetaz=param[2];
 scale=param[3];
 tx=param[4];ty=param[5];tz=param[6];
 xg=param[7];yg=param[8];zg=param[9];

 param[0]=tetax;param[1]=tetay;param[2]=tetaz;
 param[3]=param[4]=param[5]=scale;
 param[6]=tx;param[7]=ty;param[8]=tz;
 param[9]=xg;param[10]=yg;param[11]=zg;

 return(12);
}

/*******************************************************************************
**     rigidz_to_affine_decouple_3d(param)
*/
/*!
**      met a jour les valeurs du tableau param pour le passage d'une transfo
**  rigide+zoom anisotrope vers une transfo affinne avec parametres decouples\n
**  retourne le nouveau nombre de parametres
**  \param param : le tableau de param
**  \retval 12 (nbre de param)
*******************************************************************************/
int rigidz_to_affine_decouple_3d(double *param)
{
 double tetax,tetay,tetaz,sx,sy,sz,tx,ty,tz,xg,yg,zg;

 tetax=param[0];tetay=param[1];tetaz=param[2];
 sx=param[3];sy=param[4];sz=param[5];
 tx=param[6];ty=param[7];tz=param[8];
 xg=param[9];yg=param[10];zg=param[11];

 //angles d'euler
 param[0]=tetax;param[1]=tetay;param[2]=tetaz;
 //translations
 param[3]=tx; param[4]=ty;param[5]=tz;
 //scales
 param[6]=sx;param[7]=sy;param[8]=sz;
 //skews
 param[9]=param[10]=param[11]=0.0;
 //centres
 param[12]=xg;param[13]=yg;param[14]=zg;

 return(15);
}

/*******************************************************************************
**     rigidz_to_affine_3d(param)                  
*/
/*!                                                                   
**     met a jour les valeurs du tableau param pour le passage d'une transfo
**  rigide+zoom vers une transfo affine\n
**  retourne le nouveau nombre de parametres
**  \param param : le tableau de param 
**  \retval 15 (nbre de param)
*******************************************************************************/
int rigidz_to_affine_3d(double *param)
{
 double tetax,tetay,tetaz,sx,sy,sz,tx,ty,tz,xg,yg,zg;
 double cxx,sxx,cyy,syy,czz,szz;
 double a11,a12,a13,a21,a22,a23,a31,a32,a33;
 
 tetax=param[0];tetay=param[1];tetaz=param[2];
 sx=param[3];sy=param[4];sz=param[5];
 tx=param[6];ty=param[7];tz=param[8];
 xg=param[9];yg=param[10];zg=param[11];
 
 cxx=cos(tetax);sxx=sin(tetax);cyy=cos(tetay);syy=sin(tetay);czz=cos(tetaz);szz=sin(tetaz);
 a11=sx*cyy*czz;a12=-sy*cyy*szz;a13=-sz*syy;
 a21=sx*(cxx*szz-sxx*syy*czz);
 a22=sy*(cxx*czz+sxx*syy*szz);
 a23=-sz*sxx*cyy;
 a31=sx*(sxx*szz+cxx*syy*czz);
 a32=sy*(sxx*czz-cxx*syy*szz);
 a33=sz*cxx*cyy;
 param[0]=a11;param[1]=a12;param[2]=a13;
 param[3]=a21;param[4]=a22;param[5]=a23;
 param[6]=a31;param[7]=a32;param[8]=a33;
 param[9]=tx;param[10]=ty;param[11]=tz;
 param[12]=xg;param[13]=yg;param[14]=zg;
 
 return(15); 
}

/*******************************************************************************
**     affine_decouple_to_affine_3d(param)
*/
/*!
**     met a jour les valeurs du tableau param pour le passage d'une transfo
**  affine avec parametres decouples vers une transfo affine\n
**  retourne le nouveau nombre de parametres
**  \param param : le tableau de param
**  \retval 15 (nbre de param)
*******************************************************************************/
int affine_decouple_to_affine_3d(double *param)
{
  double matrice_affine[15]={0.0};
  //double **matrice_rigid_zoom, **matrice_skew;
  double matrice_rigid_zoom[3][3]={{0.0}}, matrice_skew[3][3]={{0.0}};
  int i,j,k;
  double tetax,tetay,tetaz,tx,ty,tz,xc,yc,zc,sx,sy,sz,syx,szx,szy;
  double cxx,sxx,cyy,syy,czz,szz;
  double val;

  //angle d'euler
  tetax=param[0];tetay=param[1];tetaz=param[2];
  //translations
  tx=param[3];ty=param[4];tz=param[5];
  //scales
  sx=param[6];sy=param[7];sz=param[8];
  //skews
  syx=param[9]; szx=param[10]; szy=param[11];
  //centre
  xc=param[12];yc=param[13];zc=param[14];

  //construction de la mtrice rigid+zoom
  cxx=cos(tetax);sxx=sin(tetax);cyy=cos(tetay);syy=sin(tetay);czz=cos(tetaz);szz=sin(tetaz);
  matrice_rigid_zoom[0][0]=sx*cyy*czz; matrice_rigid_zoom[0][1]=-sy*cyy*szz; matrice_rigid_zoom[0][2]=-sz*syy;
  matrice_rigid_zoom[1][0]=sx*(cxx*szz-sxx*syy*czz); matrice_rigid_zoom[1][1]=sy*(cxx*czz+sxx*syy*szz); matrice_rigid_zoom[1][2]=-sz*sxx*cyy;
  matrice_rigid_zoom[2][0]=sx*(sxx*szz+cxx*syy*czz); matrice_rigid_zoom[2][1]=sy*(sxx*czz-cxx*syy*szz); matrice_rigid_zoom[2][2]=sz*cxx*cyy;

  //matrice de skew
  matrice_skew[0][0]=1; matrice_skew[0][1]=syx; matrice_skew[0][2]=szx;
  matrice_skew[1][0]=-syx; matrice_skew[1][1]=1; matrice_skew[1][2]=szy;
  matrice_skew[2][0]=-szx; matrice_skew[2][1]=-szy; matrice_skew[2][2]=1;

  //multiplication des 2 matrices (rigide+zoom et skew) pour obtenir la globale
  for (i=0;i<3;i++)
  {
   for (j=0;j<3;j++)
   {
    val=0.0;
    for (k=0;k<3;k++)
    {
     val += matrice_rigid_zoom[i][k]*matrice_skew[k][j];
    }
    matrice_affine[3*i+j]=val;
   }
  }

  matrice_affine[9]=tx; matrice_affine[10]=ty; matrice_affine[11]=tz;
  matrice_affine[12]=xc; matrice_affine[13]=yc; matrice_affine[14]=zc;

  //on met a jour les parametres
  for (i=0;i<15;i++) param[i]=matrice_affine[i];

  return 15;
}

/*******************************************************************************
**     affine_to_affinesscg_3d(param)                  
*/
/*!                                                                   
**     met a jour les valeurs du tableau param pour le passage d'une transfo
**  affine vers une transfo affine sans le centre de gravite\n
**  retourne le nouveau nombre de parametres
**  \param param : le tableau de param 
**  \retval 12 (nbre de param)
*******************************************************************************/
int affine_to_affinesscg_3d(double *param)
{
 double tx,ty,tz,xg,yg,zg;
 double tx_aff,ty_aff,tz_aff;
 double a11,a12,a13,a21,a22,a23,a31,a32,a33;
 
 a11=param[0];a12=param[1];a13=param[2];
 a21=param[3];a22=param[4];a23=param[5];
 a31=param[6];a32=param[7];a33=param[8];
 tx=param[9];ty=param[10];tz=param[11];
 xg=param[12];yg=param[13];zg=param[14];
 tx_aff = (1.0-a11)*xg + (-a12)*yg     + (-a13)*zg         +tx;
 ty_aff = (-a21)*xg     + (1.0-a22)*yg + (-a23)*zg         +ty;
 tz_aff = (-a31)*xg     + (-a32)*yg       + (1.0-a33)*zg   +tz;
 
 param[9]  = tx_aff;
 param[10] = ty_aff;
 param[11] = tz_aff;
 param[12]=0.0;param[13]=0.0;param[14]=0.0;
 
 return(12); 
}

/*******************************************************************************
********************************************************************************
**************** GESTION DE LA BASE DE FONCTIONS *******************************
********************************************************************************
*******************************************************************************/


/**************************************************************************
**      init_base_3d(wdth,hght,dpth,scal_func)
*/
/*! 
**      Initialisation du tableau base en fonction de la taille de
**  l'image, de la resolution et de la fonction d'echelle utilisee
**  la fonction retourne un entier correspondant au nombre de parametres
**  La base est allouee pendant cette initialisation et la resolution de
**  depart est calculee automatiquement en fonction du support de la fonction

**  d'echelle
**  \param wdth, hght,dpth : la taille
**  \param scal_func : fonction d'echelle
**  \retval nbre de parametres
**  \remark la variable global BASE3D est modifiee
**************************************************************************/
int init_base_3d(int wdth, int hght, int dpth, scal_func_t scal_func)
{
 int    *x0,*x1,*y0,*y1,*z0,*z1,*idx,*idy,*idz,tfilt;
 int    nbint,pasx,pasy,pasz,nbfx,nbfy,nbfz,supx,supy,supz;
 int    i,j,k,l;
 int    resol,nb_func;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz,*px,*py,*pz,*filt;
 
 fx=fy=fz=dfx=dfy=dfz=NULL;
 /*calcul de la resolution initiale et d'autres parametres*/
 resol=-1;
 do
 {
   resol++;
   /*nombre d'intervalle de decoupage*/
    nbint=(int)(pow(2.0,(double)resol));
   /*pas suivant x et y*/ 
    pasx=wdth/nbint;pasy=hght/nbint;pasz=dpth/nbint;
   /*calcul du nombre de fonctions suivant x et y*/
    supx=(*scal_func)(pasx,&fx,&dfx);nbfx=nbint+1-supx;
    supy=(*scal_func)(pasy,&fy,&dfy);nbfy=nbint+1-supy;
    supz=(*scal_func)(pasz,&fz,&dfz);nbfz=nbint+1-supz;
    nb_func=nbfx*nbfy*nbfz;
 } while (nbfx<=0 || nbfy<=0 || nbfz<=0);
 
 /*calcul du filtre en fonction de la fonction d'echelle*/
  tfilt=2;filt=NULL;
  if (scal_func==Bsplined0) 
    {tfilt=2;filt=CALLOC(tfilt,double);filt[0]=filt[1]=1.0;}
  if (scal_func==Bsplined1)
    {tfilt=3;filt=CALLOC(tfilt,double);filt[0]=filt[2]=0.5;filt[1]=1.0;}
  if (scal_func==Bsplined2)
    {tfilt=4;filt=CALLOC(tfilt,double);filt[0]=filt[3]=1.0/4.0;filt[1]=filt[2]=3.0/4.0;}
  if (scal_func==Bsplineod1)
    {tfilt=3;filt=CALLOC(tfilt,double);filt[0]=filt[2]=-1.0;filt[1]=1.0;}
 BASE3D.tfilt=tfilt;BASE3D.filt=filt; 

   
 /*quelques caracteristique de la base*/
 BASE3D.width=wdth;BASE3D.height=hght;BASE3D.depth=dpth;
 BASE3D.nb_func=nb_func;BASE3D.nbfx=nbfx;BASE3D.nbfy=nbfy;BASE3D.nbfz=nbfz;BASE3D.resol=resol;
 BASE3D.fx=fx;BASE3D.fy=fy;BASE3D.fz=fz;BASE3D.dfx=dfx;BASE3D.dfy=dfy;BASE3D.dfz=dfz;
 BASE3D.scal_func=scal_func;
  
 /*allocation des element*/
 x0=CALLOC(nb_func,int); y0=CALLOC(nb_func,int);z0=CALLOC(nb_func,int);  
 x1=CALLOC(nb_func,int);y1=CALLOC(nb_func,int);z1=CALLOC(nb_func,int); 
 idx=CALLOC(nb_func,int);idy=CALLOC(nb_func,int);idz=CALLOC(nb_func,int);  
 px=CALLOC(nb_func,double);py=CALLOC(nb_func,double);pz=CALLOC(nb_func,double); 
 if (x0==NULL || x1==NULL || y0==NULL || y1==NULL || z0==NULL || z1==NULL|| idx==NULL || idy==NULL || idz==NULL|| px==NULL ||
 py==NULL || pz==NULL)
  {printf("ERREUR d'allocation dans init_base_3d\n");exit(1);}
 
  /*remplissage de la base*/
 l=0;
 for (i=0;i<nbfx;i++)
  for (j=0;j<nbfy;j++) 
   for (k=0;k<nbfz;k++)
    {
     x0[l]=i*pasx;x1[l]=i*pasx+supx*pasx;idx[l]=i;px[l]=0.0;
     y0[l]=j*pasy;y1[l]=j*pasy+supy*pasy;idy[l]=j;py[l]=0.0;
     z0[l]=k*pasy;z1[l]=k*pasz+supz*pasz;idz[l]=k;pz[l]=0.0;
     l++;
    }
 BASE3D.x0=x0;BASE3D.x1=x1;BASE3D.y0=y0;BASE3D.y1=y1;BASE3D.z0=z0;BASE3D.z1=z1;BASE3D.idx=idx;BASE3D.idy=idy;BASE3D.idz=idz;
 BASE3D.px=px;BASE3D.py=py;BASE3D.pz=pz;
      
 return(3*nb_func);
}


/**************************************************************************
**      init_base_resol_3d(wdth,hght,dpth,scal_func,resolution)
**  
**      Initialisation du tableau base en fonction de la taille de
**  l'image, de la resolution et de la fonction d'echelle utilisee
**  la fonction retourne un entier correspondant au nombre de parametres
**  La base est allouee pendant cette initialisation et la resolution est
**  passee l'utilisateur
**************************************************************************/
int init_base_resol_3d(int wdth, int hght, int dpth, scal_func_t scal_func, int resolution)
{
 int    *x0,*x1,*y0,*y1,*z0,*z1,*idx,*idy,*idz,tfilt;
 int    nbint,pasx,pasy,pasz,nbfx,nbfy,nbfz,supx,supy,supz;
 int    i,j,k,l;
 int    resol,nb_func;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz,*px,*py,*pz,*filt;
 
 fx=fy=fz=dfx=dfy=dfz=NULL;
 
 resol=resolution;
  /*nombre d'intervalle de decoupage*/
   nbint=(int)(pow(2.0,(double)resol));
  /*pas suivant x et y*/ 
    pasx=wdth/nbint;pasy=hght/nbint;pasz=dpth/nbint;
  /*calcul du nombre de fonctions suivant x et y*/
    supx=(*scal_func)(pasx,&fx,&dfx);nbfx=nbint+1-supx;
    supy=(*scal_func)(pasy,&fy,&dfy);nbfy=nbint+1-supy;
    supz=(*scal_func)(pasz,&fz,&dfz);nbfz=nbint+1-supz;
    nb_func=nbfx*nbfy*nbfz;
 
 /*calcul du filtre en fonction de la fonction d'echelle*/
  tfilt=2;filt=NULL;
  if (scal_func==Bsplined0)
    {tfilt=2;filt=CALLOC(tfilt,double);filt[0]=filt[1]=1.0;}
  if (scal_func==Bsplined1)
    {tfilt=3;filt=CALLOC(tfilt,double);filt[0]=filt[2]=0.5;filt[1]=1.0;}
  if (scal_func==Bsplined2)
    {tfilt=4;filt=CALLOC(tfilt,double);filt[0]=filt[3]=1.0/4.0;filt[1]=filt[2]=3.0/4.0;}
  if (scal_func==Bsplineod1)
    {tfilt=3;filt=CALLOC(tfilt,double);filt[0]=filt[2]=-1.0;filt[1]=1.0;}
 BASE3D.tfilt=tfilt;BASE3D.filt=filt; 

   
 /*quelques caracteristique de la base*/
 BASE3D.width=wdth;BASE3D.height=hght;BASE3D.depth=dpth;
 BASE3D.nb_func=nb_func;BASE3D.nbfx=nbfx;BASE3D.nbfy=nbfy;BASE3D.nbfz=nbfz;BASE3D.resol=resol;
 BASE3D.fx=fx;BASE3D.fy=fy;BASE3D.fz=fz;BASE3D.dfx=dfx;BASE3D.dfy=dfy;BASE3D.dfz=dfz;
 BASE3D.scal_func=scal_func;
  
 /*allocation des element*/
 x0=CALLOC(nb_func,int); y0=CALLOC(nb_func,int);z0=CALLOC(nb_func,int);  
 x1=CALLOC(nb_func,int);y1=CALLOC(nb_func,int);z1=CALLOC(nb_func,int); 
 idx=CALLOC(nb_func,int);idy=CALLOC(nb_func,int);idz=CALLOC(nb_func,int);  
 px=CALLOC(nb_func,double);py=CALLOC(nb_func,double);pz=CALLOC(nb_func,double); 
 if (x0==NULL || x1==NULL || y0==NULL || y1==NULL || z0==NULL || z1==NULL|| idx==NULL || idy==NULL || idz==NULL|| px==NULL ||
 py==NULL || pz==NULL)
  {printf("ERREUR d'allocation dans init_base_3d\n");exit(1);}
 
  /*remplissage de la base*/
 l=0;
 for (i=0;i<nbfx;i++)
  for (j=0;j<nbfy;j++) 
   for (k=0;k<nbfz;k++)
    {
     x0[l]=i*pasx;x1[l]=i*pasx+supx*pasx;idx[l]=i;px[l]=0.0;
     y0[l]=j*pasy;y1[l]=j*pasy+supy*pasy;idy[l]=j;py[l]=0.0;
     z0[l]=k*pasy;z1[l]=k*pasz+supz*pasz;idz[l]=k;pz[l]=0.0;
     l++;
    }
 BASE3D.x0=x0;BASE3D.x1=x1;BASE3D.y0=y0;BASE3D.y1=y1;BASE3D.z0=z0;BASE3D.z1=z1;BASE3D.idx=idx;BASE3D.idy=idy;BASE3D.idz=idz;
 BASE3D.px=px;BASE3D.py=py;BASE3D.pz=pz;
      
 return(3*nb_func);
}




/*
  resize the base (the resolution and other parameters are not altered)

  This enables us to change the image resolution during matching
*/

int resize_base_3d(int wdth, int hght, int dpth)
{
  int pasx, pasy, pasz;
  int resol, nbint;
  int nbfx, nbfy, nbfz, supx, supy, supz;
  double *fx, *fy, *fz, *dfx, *dfy, *dfz;
  int (*scal_func)();
  int *x0, *x1, *y0, *y1, *z0, *z1;
  int i, j, k, l;


  /*parametre de la base dans les variables locales*/
  resol=BASE3D.resol;
  fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;
  scal_func=BASE3D.scal_func;
  /*nb_func=BASE3D.nb_func;*/nbfx=BASE3D.nbfx;nbfy=BASE3D.nbfy;nbfz=BASE3D.nbfz;
  x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;

  /*mise a jour des differents parametres de la base*/
  /*nombre d'intervalle de decoupage*/
  nbint=/*(int)(pow(2.0,(double)resol))*/ 1 << resol;
  /*pas suivant x y et z*/ 
  pasx=wdth/nbint;pasy=hght/nbint;pasz=dpth/nbint;
  /*mise a jour de l'echantillonage de la fonction d'echelle*/
  /* on libere d'abord tout pour empecher la fragmentation du tas */
  free(fx);
  free(fy);
  free(fz);
  fx = NULL;
  fy = NULL;
  fz = NULL;
  free(dfx);
  free(dfy);
  free(dfz);
  dfx = NULL;
  dfy = NULL;
  dfz = NULL;
  supx=(*scal_func)(pasx,&fx,&dfx);
  supy=(*scal_func)(pasy,&fy,&dfy);
  supz=(*scal_func)(pasz,&fz,&dfz);

  l=0;
  for (i=0; i<nbfx; i++)
   for (j=0; j<nbfy; j++) 
    for (k=0; k<nbfz; k++)
    {
      x0[l] = i*pasx;
      x1[l] = i*pasx + supx*pasx;
      y0[l] = j*pasy;
      y1[l] = j*pasy + supy*pasy;
      z0[l] = k*pasy;
      z1[l] = k*pasz + supz*pasz;
      l++;
    }


  BASE3D.fx=fx;BASE3D.fy=fy;BASE3D.fz=fz;BASE3D.dfx=dfx;BASE3D.dfy=dfy;BASE3D.dfz=dfz;
  BASE3D.width=wdth;BASE3D.height=hght;BASE3D.depth=dpth;


  return 0;
}



/**************************************************************************
**      base_resol_up_3d(param,nb_param)
*/  
/*!      Passage de la base de fonction a la resolution superieure
**  Modifie les parametres en consequence :
**      nouvelle taille en resultat de la fonction
**      nouvelles valeurs dans le tableau param
**  \retval le nombre de parametres
**  \remark modifie la variable globale BASE
**************************************************************************/

int base_resol_up_3d(double *param, int nb_param)
{
 int    i,j,k,l;
 int    pasx,pasy,pasz,wdth,hght,dpth;
 int    resol,nbint,tfilt,nb_func,nbfx,nbfy,nbfz,supx,supy,supz;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz,*px,*py,*pz,*filt;
 int    *x0,*x1,*y0,*y1,*z0,*z1,*idx,*idy,*idz;
 double ***Mx,***My,***Mz,***F3D;
 int    i0,j0,k0;
 int    (*scal_func)();
 
  wdth=BASE3D.width;hght=BASE3D.height;dpth=BASE3D.depth;
  
 /*parametre de la base dans les variables locales*/
  resol=BASE3D.resol;resol++; /*passage a la resolution superieure*/
  fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;
  filt=BASE3D.filt;tfilt=BASE3D.tfilt;
  x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;
  px=BASE3D.px;py=BASE3D.py;pz=BASE3D.pz;idx=BASE3D.idx;idy=BASE3D.idy;idz=BASE3D.idz;
  nb_func=BASE3D.nb_func;nbfx=BASE3D.nbfx;nbfy=BASE3D.nbfy;nbfz=BASE3D.nbfz;
  scal_func=BASE3D.scal_func;
 /*les parametres du tableau param sont mis dans la base*/
  if (nb_param!=3*nb_func) {printf("ERREUR dans base_resol_up \n");return(1);} 
  for (l=0;l<nb_func;l++) {px[l]=param[3*l];py[l]=param[3*l+1];pz[l]=param[3*l+2];}
  
 /*allocation et remplissage des matrices de parametres Mx, My et Mz*/
 /*c'est sur ces matrices que l'on appliquera le filtre de passage d'une resolution a une
 superieure*/ 
   /*calcul de la taille des matrices apres zero filling, allocation memoire et mise a zero*/
    nbfx=tfilt+2*nbfx-2;nbfy=tfilt+2*nbfy-2;nbfz=tfilt+2*nbfz-2;
    Mx=alloc_dmatrix_3d(nbfx,nbfy,nbfz);My=alloc_dmatrix_3d(nbfx,nbfy,nbfz);Mz=alloc_dmatrix_3d(nbfx,nbfy,nbfz);
    for (i=0;i<nbfx;i++) for (j=0;j<nbfy;j++) for (k=0;k<nbfz;k++) Mx[i][j][k]=My[i][j][k]=Mz[i][j][k]=0.0;
   /*allocation et calcul du filtre de convolution 2D*/
    F3D=alloc_dmatrix_3d(tfilt,tfilt,tfilt);
    for (i=0;i<tfilt;i++) 
     for (j=0;j<tfilt;j++)
      for (k=0;k<tfilt;k++)
       F3D[i][j][k]=filt[i]*filt[j]*filt[k];
      
   /*calcul de Mx My et Mz par application du filtre de convolution 2D*/
    for (l=0;l<nb_func;l++)
     {
      i0=2*idx[l];j0=2*idy[l];k0=2*idz[l];
      for (i=i0;i<i0+tfilt;i++)
       for (j=j0;j<j0+tfilt;j++)
        for (k=k0;k<k0+tfilt;k++)
         {
      Mx[i][j][k]+=px[l]*F3D[i-i0][j-j0][k-k0];
      My[i][j][k]+=py[l]*F3D[i-i0][j-j0][k-k0];
      Mz[i][j][k]+=pz[l]*F3D[i-i0][j-j0][k-k0];
     }
     }   
 
 /*mise a jour des differents parametres de la base*/
   /*nombre d'intervalle de decoupage*/
    nbint=(int)(pow(2.0,(double)resol));
   /*pas suivant x y et z*/ 
    pasx=wdth/nbint;pasy=hght/nbint;pasz=dpth/nbint;
   /*mise a jour de l'echantillonage de la fonction d'echelle*/
    supx=(*scal_func)(pasx,&fx,&dfx);
    supy=(*scal_func)(pasy,&fy,&dfy);
    supz=(*scal_func)(pasz,&fz,&dfz);
   /*verification que ca correspond avec ce qui est deja calcule*/
    if (nbfx!=nbint+1-supx || nbfy!=nbint+1-supy || nbfz!=nbint+1-supz) printf("ERREUR 1 dans base_resol_up_3d\n");
    nb_func=nbfx*nbfy*nbfz;
  /*reallocation des element*/
    x0=REALLOC(x0,nb_func,int);
    y0=REALLOC(y0,nb_func,int);
    z0=REALLOC(z0,nb_func,int);  
    x1=REALLOC(x1,nb_func,int);y1=REALLOC(y1,nb_func,int);z1=REALLOC(z1,nb_func,int); 
    idx=REALLOC(idx,nb_func,int);idy=REALLOC(idy,nb_func,int);idz=REALLOC(idz,nb_func,int); 
    px=REALLOC(px,nb_func,double);py=REALLOC(py,nb_func,double);pz=REALLOC(pz,nb_func,double); 
    if (x0==NULL || x1==NULL || y0==NULL || y1==NULL || z0==NULL || z1==NULL || idx==NULL || idy==NULL || idz==NULL 
    || px==NULL || py==NULL || pz==NULL)
     {printf("ERREUR d'allocation dans base_resol_up_3d\n");exit(1);}
  /*remplissage de la base*/
   l=0;
    for (i=0;i<nbfx;i++)
     for (j=0;j<nbfy;j++) 
      for (k=0;k<nbfz;k++)
       {
        x0[l]=i*pasx;x1[l]=i*pasx+supx*pasx;idx[l]=i;px[l]=Mx[i][j][k];
        y0[l]=j*pasy;y1[l]=j*pasy+supy*pasy;idy[l]=j;py[l]=My[i][j][k];
    z0[l]=k*pasz;z1[l]=k*pasz+supz*pasz;idz[l]=k;pz[l]=Mz[i][j][k];
        l++;
       }
   BASE3D.x0=x0;BASE3D.x1=x1;BASE3D.y0=y0;BASE3D.y1=y1;BASE3D.z0=z0;BASE3D.z1=z1;BASE3D.idx=idx;BASE3D.idy=idy;BASE3D.idz=idz;
   BASE3D.px=px;BASE3D.py=py;BASE3D.pz=pz;
   BASE3D.resol=resol;BASE3D.nb_func=nb_func;BASE3D.nbfx=nbfx;BASE3D.nbfy=nbfy;BASE3D.nbfz=nbfz;
   BASE3D.fx=fx;BASE3D.fy=fy;BASE3D.fz=fz;BASE3D.dfx=dfx;BASE3D.dfy=dfy;BASE3D.dfz=dfz;
  
  /*nouvelles valeurs des parametres dans le tableau param*/
   for (l=0;l<nb_func;l++)
    {param[3*l]=px[l];param[3*l+1]=py[l];param[3*l+2]=pz[l];}
  
 free_dmatrix_3d(Mx);free_dmatrix_3d(My);free_dmatrix_3d(Mz);free_dmatrix_3d(F3D);  
 return(3*nb_func);
}

/**************************************************************************
**      reduc_base_3d(param,nb_param,imref,imreca,imres)
*/  
/*!      Reduit la base de fonction en tenant compte des info contenu dans
**  imref, imreca et imres
**  \retval retourne le nouveau nombre de parametre et met a jour la valeur
**  du tableau param, 1 si echec (pas de modif)
**************************************************************************/
int reduc_base_3d(double *param, int nb_param, grphic3d *imref, grphic3d *imreca, grphic3d *imres)
{
 int    wdth,hght,dpth,i,j,k,l,m,nb_func,i0,i1,j0,j1,k0,k1,flag;
 double *px,*py,*pz;  
 int    *idx,*idy,*idz,*x0,*y0,*z0,*x1,*y1,*z1;   

 wdth=BASE3D.width;hght=BASE3D.height;dpth=BASE3D.depth;

 if ((int)imref->width!=wdth || (int)imref->height!=hght || (int)imref->depth!=dpth ||
     (int)imreca->width!=wdth || (int)imreca->height!=hght || (int)imreca->depth!=dpth ||
     (int)imres->width!=wdth || (int)imres->height!=hght || (int)imres->depth!=dpth ) {printf("ERREUR dans reduc_base_3d\n");return(1);} 
     
 nb_func=BASE3D.nb_func;
 px=BASE3D.px;py=BASE3D.py;pz=BASE3D.pz;
 idx=BASE3D.idx;idy=BASE3D.idy;idz=BASE3D.idz;x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;
 
 for (l=0;l<nb_func;l++) {px[l]=param[3*l];py[l]=param[3*l+1];pz[l]=param[3*l+2];}
 
 m=0;
 for (l=0;l<nb_func;l++)
 {
  i0=x0[l];i1=x1[l];j0=y0[l];j1=y1[l];k0=z0[l];k1=z1[l];
  i=i0;flag=0;
  do 
   { 
    j=j0;
    do 
     {
      k=k0;
      do
       {
        if (imref->mri[i][j][k]!=0 || imreca->mri[i][j][k]!=0 || imres->mri[i][j][k]!=0) flag=1;
        k++;
       } while (k<k1 && flag==0);
      j++;
     } while (j<j1 && flag==0);
    i++;
   } while (i<i1 && flag==0);

  if (flag!=0) /*cette fonction de la base a lieu d'exister*/
  {
   x0[m]=x0[l];x1[m]=x1[l];y0[m]=y0[l];y1[m]=y1[l];z0[m]=z0[l];z1[m]=z1[l];
   idx[m]=idx[l];idy[m]=idy[l];idz[m]=idz[l];px[m]=px[l];py[m]=py[l];pz[m]=pz[l];
   m++;
  }
 }
 nb_func=m;
 x0=REALLOC(x0,nb_func,int);x1=REALLOC(x1,nb_func,int);
 y0=REALLOC(y0,nb_func,int);y1=REALLOC(y1,nb_func,int);
 z0=REALLOC(z0,nb_func,int);z1=REALLOC(z1,nb_func,int);
 idx=REALLOC(idx,nb_func,int);idy=REALLOC(idy,nb_func,int);idz=REALLOC(idz,nb_func,int);
 px=REALLOC(px,nb_func,double);py=REALLOC(py,nb_func,double);pz=REALLOC(pz,nb_func,double);
 
 BASE3D.nb_func=nb_func;
 BASE3D.x0=x0;BASE3D.x1=x1;BASE3D.y0=y0;BASE3D.y1=y1;BASE3D.z0=z0;BASE3D.z1=z1;
 BASE3D.idx=idx;BASE3D.idy=idy;BASE3D.idz=idz;BASE3D.px=px;BASE3D.py=py;BASE3D.pz=pz;
 
 for (l=0;l<nb_func;l++)
  {param[3*l]=px[l];param[3*l+1]=py[l];param[3*l+2]=pz[l];}  
   
 return(3*nb_func);
}
/**************************************************************************
**      end_base_3d
*/
/*! 
**      Terminaison de la base
**  Liberation memoire des different elements
**  \remark la variable globale BASE3D est modifie
**************************************************************************/
int end_base_3d(void)
{
 if (BASE3D.x0)free(BASE3D.x0);
 if (BASE3D.x1) free(BASE3D.x1);
 if (BASE3D.y0) free(BASE3D.y0);
 if (BASE3D.y1) free(BASE3D.y1);
 if (BASE3D.z0) free(BASE3D.z0);
 if (BASE3D.z1) free(BASE3D.z1);
 if (BASE3D.idx) free(BASE3D.idx);
 if (BASE3D.idy) free(BASE3D.idy);
 if (BASE3D.idz) free(BASE3D.idz);
 if (BASE3D.px) free(BASE3D.px);
 if (BASE3D.py) free(BASE3D.py);
 if (BASE3D.pz) free(BASE3D.pz);
 if (BASE3D.fx) free(BASE3D.fx);
 if (BASE3D.fy) free(BASE3D.fy);
 if (BASE3D.fz) free(BASE3D.fz);
 if (BASE3D.dfx) free(BASE3D.dfx);
 if (BASE3D.dfy) free(BASE3D.dfy);
 if (BASE3D.dfz) free(BASE3D.dfz);
 if (BASE3D.filt) free(BASE3D.filt);
 return(1);
}

/*******************************************************************************
********************************************************************************
*************** ALIGNEMENT DES CENTRES DE GRAVITE ******************************
********************************************************************************
*******************************************************************************/

/*******************************************************************************
**     imx_find_cg_with_seuil_3d_p(img,cgX,cgY,cgZ)
*/
/*!
**      trouve le centre de gravite d'une image a l'aide d'un seuillage (otsu)
**  \param img : l'image en question
**  \param cgX : (out) coordX du centre de gravite (en voxel)
**  \param cgY : (out) coordY du centre de gravite (en voxel)
**  \param cgZ : (out) coordZ du centre de gravite (en voxel)
**  \return la valeur du seuil qui a ete retenue pour le calcul
**
*******************************************************************************/
TYPEMRI3D imx_find_cg_with_seuil_3d_p(grphic3d *img,double *cgX,double *cgY,double *cgZ)
{
  int threshold;
  unsigned int width,height,depth;
  unsigned int i,j,k;
  double val;
  double somme = 0.0;
  double seuil;
  double rcoeff;
  TYPEMRI3D ***mri;

  if (!img)
  {
      PUT_WARN("ERROR in imx_find_cg_with_seuil_3d_p\n"); return 0;
  }

  width = img->width;
  height= img->height;
  depth = img->depth;

  find_threshold_3d(img,&threshold,1,FALSE);
  rcoeff = img->rcoeff;
  seuil = threshold*rcoeff;
  mri = img->mri;
  *cgX = 0.;*cgY = 0.;*cgZ = 0.;

  for(i=0;i<width;i++)
  for(j=0;j<height;j++)
  for(k=0;k<depth;k++)
  {
      val = mri[i][j][k]*rcoeff;
     if (val>seuil) {(*cgX)+=(double)i*val;(*cgY)+=(double)j*val;(*cgZ)+=(double)k*val;somme+=val;}
  }

  if (somme)
  {
    *cgX=(*cgX)/somme;
    *cgY=(*cgY)/somme;
    *cgZ=(*cgZ)/somme;
  }

  return (TYPEMRI3D)threshold;
}

/*******************************************************************************
**     imx_find_cg_geometric_3d_p(img,cgX,cgY,cgZ)
*/
/*!
**      trouve le centre geometrique d'une image (il vaut mieux qu'elle soit seuille!)
**  \param img : l'image en question
**  \param cgX : (out) coordX du centre de gravite (en voxel)
**  \param cgY : (out) coordY du centre de gravite (en voxel)
**  \param cgZ : (out) coordZ du centre de gravite (en voxel)
**  \return 0
**
*******************************************************************************/
int imx_find_cg_geometric_3d_p(grphic3d *img,double *cgX,double *cgY,double *cgZ)
{
  unsigned int width,height,depth;
  unsigned int i,j,k;
  int somme = 0;
  TYPEMRI3D ***mri;

  if (!img)
  {
      PUT_WARN("ERROR in imx_find_cg_geometric_3d_p\n"); 
      return 1;
  }

  width = img->width;
  height= img->height;
  depth = img->depth;

  mri = img->mri;
  *cgX = 0.;*cgY = 0.;*cgZ = 0.;

  for(i=0;i<width;i++)
  for(j=0;j<height;j++)
  for(k=0;k<depth;k++)
  {
      if  (mri[i][j][k])
      {
        (*cgX)+=(double)i;
        (*cgY)+=(double)j;
        (*cgZ)+=(double)k;
        somme++;
      }
  }

  if (somme)
  {
    *cgX=(*cgX)/somme;
    *cgY=(*cgY)/somme;
    *cgZ=(*cgZ)/somme;
  }

  return 0;
}



/*******************************************************************************
**     gravite_to_field_3d(imref,imreca,champ)                    
*/
/*!                                                                   
**     Champ permettant de superposer les centre de gravite de imref et imreca 
**     ATTENTION : les images sont supposees deja seuillees
**    \param imref, imreca : image reference et recale
**    \retval En resultat, la fonction retourne un tableau (alloue ici) contenant
**  tx, ty, tz,  xg, yg, zg 
*******************************************************************************/
double *gravite_to_field_3d(grphic3d *imref, grphic3d *imreca, field3d *champ)
{
 double xgref,ygref,zgref,xgreca,ygreca,zgreca,tx,ty,tz;
 double *res;
 double tmp[16];
 int i;
  
 //calcul des centres de gravites
 xgref=ygref=zgref=0.0;xgreca=ygreca=zgreca=0.0;

 imx_find_cg_with_seuil_3d_p(imref,&xgref,&ygref,&zgref);
 imx_find_cg_with_seuil_3d_p(imreca,&xgreca,&ygreca,&zgreca);

 tx=xgreca*imreca->dx-xgref*imref->dx;
 ty=ygreca*imreca->dy-ygref*imref->dy;
 tz=zgreca*imreca->dz-zgref*imref->dz; 

 //resultats
 res=CALLOC(7,double);
 res[0]=tx;res[1]=ty;res[2]=tz;
 res[3]=xgref*imref->dx;res[4]=ygref*imref->dy;res[5]=zgref*imref->dz;

 for (i=0;i<6;i++)
  tmp[i] = res[i];
 gravite_to_rigid_3d(tmp);
 if (champ) rigid_to_field_3d(9,tmp,champ,imref,imreca);

 
 aff_log("translation: (%.2f %.2f %.2f) \n",tx,ty,tz);
 aff_log("centre: (%.2f %.2f %.2f) \n",res[3],res[4],res[5]);
 return(res);
}

/*******************************************************************************
********************************************************************************
************************** MATCHING  RIGIDE ************************************
********************************************************************************
*******************************************************************************/

void query_simple_matching_parameters_3d( int (*matching_3d_func)(int, int, int, int, int, int, int, char *, char *,eResearchInterv, eMatchPrecision))
{
 int im_ref,im_reca,im_res;
 char *nomfichres=NULL;
 char *nomfichtrf=NULL;
 int dist_type,inter_type,min_type,save_type;
 eResearchInterv FieldOfResearch;

 //demande des images
 im_reca=GET_PLACE3D(TEXT0030);
 im_ref=GET_PLACE3D(TEXT0031);
 im_res=GET_PLACE3D(TEXT0006);

 //distance
 dist_type=imx_query_erreur_type_3D();
 //interpolation
 inter_type=imx_query_interpolation_type_3D(dist_type);
 //minimisation
 min_type=imx_query_minimisation_type_3D(dist_type);
 //fichier d'enregistrement du resultat
 nomfichres=quest_save_result_3d(&save_type);
 //intervalle de recherche
 FieldOfResearch = MEDIUM;

 //matching
 matching_3d_func(im_ref,im_reca,im_res,dist_type,inter_type,min_type,save_type,nomfichres,nomfichtrf,FieldOfResearch,NORMAL);

 if (nomfichres) free(nomfichres);
 if (nomfichtrf) free(nomfichtrf);

 show_picture_3d(im_res);
}

void query_simple_precision_matching_parameters_3d( int (*matching_3d_func)(int, int, int, int, int, int, int, char *, char *,eResearchInterv, eMatchPrecision))
{
 int im_ref,im_reca,im_res;
 char *nomfichres=NULL;
 char *nomfichtrf=NULL;
 int dist_type,inter_type,min_type,save_type;
 eResearchInterv FieldOfResearch;
 eMatchPrecision matchPrecision;

 //demande des images
 im_reca=GET_PLACE3D(TEXT0030);
 im_ref=GET_PLACE3D(TEXT0031);
 im_res=GET_PLACE3D(TEXT0006);

 //distance
 dist_type=imx_query_erreur_type_3D();
 //interpolation
 inter_type=imx_query_interpolation_type_3D(dist_type);
 //minimisation
 min_type=imx_query_minimisation_type_3D(dist_type);
 //precision voulue
 matchPrecision=imx_query_matching_precision();
 //fichier d'enregistrement du resultat
 nomfichres=quest_save_result_3d(&save_type);
 //intervalle de recherche
 FieldOfResearch = MEDIUM;

 //matching
 matching_3d_func(im_ref,im_reca,im_res,dist_type,inter_type,min_type,save_type,nomfichres,nomfichtrf,FieldOfResearch,matchPrecision);

 if (nomfichres) free(nomfichres);
 if (nomfichtrf) free(nomfichtrf);

 show_picture_3d(im_res);
}

void query_multires_matching_parameters_3d( int (*matching_3d_func)(int, int, int, int, int, int, int, char *, char *, eMatchPrecision, int, int))
{
 int im_ref,im_reca,im_res;
 char *nomfichres=NULL;
 char *nomfichtrf=NULL;
 int dist_type,inter_type,min_type,save_type,start_resol,end_resol;
 int err=0;

 //demande des images
 im_reca=GET_PLACE3D(TEXT0030);
 im_ref=GET_PLACE3D(TEXT0031);
 im_res=GET_PLACE3D(TEXT0006);

 //distance
 dist_type=imx_query_erreur_type_3D();
 //interpolation
 inter_type=imx_query_interpolation_type_3D(dist_type);
 //minimisation
 min_type=imx_query_minimisation_type_3D(dist_type);
 //nbre de resolution
 start_resol = GET_INT("resolution large ?", 2, &err);
 if (err)
 { PUT_WARN("CANCELED\n"); return;}
 end_resol = GET_INT("resolution fine ? ", 0, &err);
 if (err)
 { PUT_WARN("CANCELED\n"); return;}

 //question sur l'enregistrement du champ resultat
 nomfichres=quest_save_result_3d(&save_type);

 matching_3d_func(im_ref,im_reca,im_res,dist_type,inter_type,min_type,save_type,nomfichres,nomfichtrf,NORMAL,start_resol,end_resol);

 if (nomfichres) free(nomfichres);
 if (nomfichtrf) free(nomfichtrf);

 show_picture_3d(im_res);
}

void query_multires_precision_matching_parameters_3d( int (*matching_3d_func)(int, int, int, int, int, int, int, char *, char *, eMatchPrecision, int, int))
{
 int im_ref,im_reca,im_res;
 char *nomfichres=NULL;
 char *nomfichtrf=NULL;
 int dist_type,inter_type,min_type,save_type,start_resol,end_resol;
 int err=0;
 eMatchPrecision matchPrecision;

 //demande des images
 im_reca=GET_PLACE3D(TEXT0030);
 im_ref=GET_PLACE3D(TEXT0031);
 im_res=GET_PLACE3D(TEXT0006);

 //distance
 dist_type=imx_query_erreur_type_3D();
 //interpolation
 inter_type=imx_query_interpolation_type_3D(dist_type);
 //minimisation
 min_type=imx_query_minimisation_type_3D(dist_type);
 //nbre de resolution
 start_resol = GET_INT("resolution large ?", 2, &err);
 if (err)
 { PUT_WARN("CANCELED\n"); return;}
 end_resol = GET_INT("resolution fine ? ", 0, &err);
 if (err)
 { PUT_WARN("CANCELED\n"); return;}
 //precision voulue
 matchPrecision=imx_query_matching_precision();

 //question sur l'enregistrement du champ resultat
 nomfichres=quest_save_result_3d(&save_type);

 matching_3d_func(im_ref,im_reca,im_res,dist_type,inter_type,min_type,save_type,nomfichres,nomfichtrf,matchPrecision,start_resol,end_resol);

 if (nomfichres) free(nomfichres);
 if (nomfichtrf) free(nomfichtrf);

 show_picture_3d(im_res);
}

/*******************************************************************************
**        matching_rigide_3d                                                 
**                                                                        
**                             
*******************************************************************************/
/*! 
     recalage rigide de deux images 3D (sans zoom)        
  
*/
void matching_rigide_3d()
{
  query_simple_matching_parameters_3d(imx_matching_rigide_3d);
}

/*!
     recalage rigide de deux images 3D (sans zoom)

*/
void matching_rigide_precision_3d()
{
 query_simple_precision_matching_parameters_3d(imx_matching_rigide_3d);
}

/*******************************************************************************
**        imx_matching_rigide_3d                                             
*/                                                                       
/*!       recalage rigide sans zoom de 2 images 3D                                   
** 
**  \param im_ref : numero de l'image de reference, 
**  \param im_reca  : numero de l'image a recaler,
**  \param im_res : numero de l'image resultat
**  les autres parametres sont expliques dans ::imx_matching_rigide_3d_p    
**  \retval 1
**  
*******************************************************************************/
int imx_matching_rigide_3d(int im_ref, int im_reca, int im_res, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres, char *nomfichtrf,eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision)
{
  grphic3d *imref,*imreca,*imres;
  transf3d *inittrf=NULL;
  transf3d *transfres=NULL;
 
  imref=ptr_img_3d(im_ref);
  imreca=ptr_img_3d(im_reca);
  imres=ptr_img_3d(im_res);
   
  if (nomfichtrf !=NULL)
     inittrf=load_transf_3d(nomfichtrf);
  transfres = cr_transf3d_p(imref,NULL); 

  init_log(NULL);
  aff_log("\n RECALAGE RIGIDE SANS ZOOM ");
  aff_log("\n");
  aff_log("recalage de %s sur %s \n",(imreca->patient).name,(imref->patient).name);

  imx_inimaxminpixel_3d_p(imref);
  imx_inimaxminpixel_3d_p(imreca);
  if ((imref->min_pixel<0)||(imreca->min_pixel<0))
  { PUT_WARN("L'image contient des valeurs negatives!!\nOn ajoute -val_min pour la recaler\n"); }

  imx_matching_rigide_3d_p(imref,imreca,imres,dist_type,inter_type,min_type,save_type,transfres,inittrf,FieldOfResearch,matchPrecision);

  /*enregistrement du champ resultat*/
  if (nomfichres!=NULL)
  { 
    save_transf_3d(transfres,nomfichres);
  }

  end_log();

  if (inittrf)
    free_transf3d(inittrf);

  free_transf3d(transfres);
 
  return(1);
}

/*! \ingroup Recalage  @{ */
/*******************************************************************************
**        imx_matching_rigide_3d_p                                            
*/
/*!                                                                        
**       \brief recalage rigide entre 2 images 3D                                    
**
**  \param imref : images de reference, 
**  \param imreca  : image a recaler,
**  \param imres : images resultat (E/S)
**  \param dist_type : numero de la fonction de distance 
**                          \li 0 : Quadratic
**                          \li 1 : Quadratic Robust
**                          \li 2 : woods
**                          \li 3 : woods robust
**                          \li 4 : IM Old
**                          \li 5 : IM
**                          \li 6 : IM NORM STUDHOLME
**                          \li 7 : IM NORM MAES1
**                          \li 8 : IM NORM MAES2
**                          \li 9 : Entropie conjointe
**                          \li 10 : Correlation ratio
**  \param inter_type : numero de la fonction d'interpolation des images
**                          \li 0 : nearest
**                          \li 1 : linear
**                          \li 2 : sin card
**                          \li 3 : quick sin card2
**                          \li 4 : quick sin card3 
**                          \li 5 : bspline 2
**                          \li 6 : bspline 3
**                          \li 7 : bspline 4
**                          \li 8 : bspline 5
**  \param min_type : numero de la fonction de minimisation utilisee aux plus faibles
    echelles de la transformation 
**                          \li 0 : ICM
**                          \li 1 : simplex
**                          \li 2 : Descente gradient
**                          \li 3 : Gradient conjugue
**                          \li 4 : Descente gradient modifiee
**                          \li 5 : Quasi Newton modifiee
**  \param save_type : Type de sauvegarde
**                          \li 0 : champ
**                          \li 1 : parametres
**  \param nomfichres :  nom du fichier de sauvegarde des transformations (*.trf)
**  \param nomfichtrf : initialisation du champ de depart a la valeur d'un fichier (*.trf)\n
**                      si NULL -> alignement centre de gravite
**  \param FieldOfResearch :  voir init_bornes_matching
**  \param matchPrecision  :  voir init_bornes_matching
**  \retval 1
**
*******************************************************************************/
int imx_matching_rigide_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int dist_type, int inter_type, int min_type, int save_type, transf3d *transfres, transf3d *inittrf,eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision)
{
 grphic3d *imtref,*imtreca,*imtres;
 field3d *champ1;
 InterpolationFct interpol;
 dist_func_t  distance;
 min_func_t minimisation;
 double  mini;
 int    tdebut,tfin;
 transf3d *preProcessingRefTrf=NULL, *preProcessingRecaTrf=NULL, *transfoini=NULL;
 
#ifndef WIN32
 setvbuf(stdout, (char *)NULL, _IONBF, 0); //setbuf(stdout,NULL);
#endif
 tdebut=time(NULL);
 
 //---------mise en place des parametres du recalage-------/
  
  distance=imx_choose_distance_fct(dist_type);
  
  interpol=imx_choose_interpolation_fct(inter_type);
     
  minimisation=imx_choose_minimisation_fct(min_type,dist_type);  

 // ------------------PREPROCESS----------------------
 {
  imx_matching_preprocessing_3d_p(imref, imreca, &imtref, &imtreca, &imtres, &preProcessingRefTrf, &preProcessingRecaTrf, distance);
 }

 // ------------------INITIALISATION PAR TRANSFO DE DEPART---------------------- 
 {
  //on veut une transformation de depart compatible avec le recalage courant
  if ((inittrf)&&(inittrf->typetrans==RIGID3D))
  {
   transfoini=imx_trouver_transf_initiale_matching_3d_p(imtref, imtreca, inittrf);
  }
  else transfoini=imx_trouver_transf_initiale_matching_3d_p(imtref, imtreca, NULL);
 }

 //------------------RECALAGE RIGIDE----------------------------------------/
 aff_log("\n TRANSFO RIGIDE \n");
 {
  champ1=cr_field3d(imtref->width,imtref->height,imtref->depth);

  mini=imx_matching_simple_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, FieldOfResearch, matchPrecision, NULL, NULL, champ1);

  free_field3d(champ1);
 }

 //on remet le tout dans le repere des images originales
 imx_copie_dim_param_to_transf3d_p(imref, transfres);
 if (preProcessingRecaTrf) comb_transf_3d(preProcessingRecaTrf,transfres,transfres,0);
 if (preProcessingRefTrf) comb_transf_3d(transfres,preProcessingRefTrf,transfres,0);

 if (preProcessingRecaTrf) free_transf3d(preProcessingRecaTrf);
 if (preProcessingRefTrf) free_transf3d(preProcessingRefTrf);
 if (imtres) free_grphic3d(imtres);
 if (imtref!=imref) free_grphic3d(imtref);
 if (imtreca!=imreca) free_grphic3d(imtreca);
 
 // ------------------CALCUL DE L'IMAGE A AFFICHER ET DE LA TRANSFO A SAUVEGARDER----------------------
 imx_aff_param(transfres->typetrans, transfres->param);
 imx_calc_transfo_sav_image_res_3d_p(imref, imreca, imres, transfres, save_type);

 //liberation memoire des variables dynamiques 
 if (transfoini) free_transf3d(transfoini);
 
 tfin=time(NULL);
 //affichage du temps de calcul total 
 {
  imx_aff_temps_calcul(tdebut, tfin);
 }

 return(1);  
}
/*! @} */

/*! 
     recalage rigide de deux images 3D (sans zoom)  multiresolution      
*/
void matching_rigide_multi_3d()
{
 query_multires_matching_parameters_3d(imx_matching_rigide_multi_3d);
}

/*!
     recalage rigide de deux images 3D (sans zoom)  multiresolution
*/
void matching_rigide_multi_precision_3d()
{
 query_multires_precision_matching_parameters_3d(imx_matching_rigide_multi_3d);
}

/*******************************************************************************
**        imx_matching_rigide_multi_3d                                             
*/                                                                       
/*!       recalage rigide multiresolution sans zoom de 2 images 3D                                  
** 
**  \param im_ref : numero de l'image de reference, 
**  \param im_reca  : numero de l'image a recaler,
**  \param im_res : numero de l'image resultat
**  les autres parametres sont expliques dans ::imx_matching_rigide_3d_p    
**  \retval 1
**  
*******************************************************************************/
int imx_matching_rigide_multi_3d(int im_ref, int im_reca, int im_res, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres, char *nomfichtrf, eMatchPrecision matchPrecision, int start_resol, int end_resol)
{
   grphic3d *imref=NULL, *imreca=NULL, *imres=NULL;
   transf3d *inittrf=NULL;
   transf3d *transfres=NULL;
   
   imref=ptr_img_3d(im_ref);
   imreca=ptr_img_3d(im_reca);
   imres=ptr_img_3d(im_res);

   if (nomfichtrf !=NULL) // initialisation de la transfo par un fichier.
     inittrf=load_transf_3d(nomfichtrf);
   transfres=cr_transf3d_p(imref,NULL);

  imx_inimaxminpixel_3d_p(imref);
  imx_inimaxminpixel_3d_p(imreca);
  if ((imref->min_pixel<0)||(imreca->min_pixel<0))
  { PUT_WARN("L'image contient des valeurs negatives!!\nOn ajoute -val_min pour la recaler\n"); }

   init_log(NULL);
   aff_log("\n RECALAGE RIGIDE SANS ZOOM EN MULTIRESOLUTION");
   aff_log("\n");
   aff_log("recalage de %s sur %s \n",(imreca->patient).name,(imref->patient).name);

   imx_matching_rigide_multi_3d_p(imref,imreca,imres,dist_type,inter_type,min_type,save_type,transfres,inittrf,matchPrecision,start_resol,end_resol);

   /*enregistrement du champ resultat*/
   if (nomfichres!=NULL)
   { 
     save_transf_3d(transfres,nomfichres);
   } 

   if (inittrf)
     free_transf3d(inittrf);
   free_transf3d(transfres);

   end_log();
   return(1);
}



/*******************************************************************************
**        imx_matching_rigide_mulit_3d_p                                             
*/                                                                       
/*!       recalage rigide multiresolution sans zoom de 2 images 3D                                  
** 
**  \param im_ref : numero de l'image de reference, 
**  \param im_reca  : numero de l'image a recaler,
**  \param im_res : numero de l'image resultat
**  les autres parametres sont expliques dans ::imx_matching_rigide_3d_p    
**  \retval 1
**  
*******************************************************************************/
int imx_matching_rigide_multi_3d_p( grphic3d *imref, 
                                    grphic3d *imreca, 
                                    grphic3d *imres, 
                                    int dist_type, 
                                    int inter_type, 
                                    int min_type, 
                                    int save_type, 
                                    transf3d *transfres, 
                                    transf3d *inittrf,
                                    eMatchPrecision matchPrecision,
                                    int start_resol,
                                    int end_resol)
{
  grphic3d *imref_preprocessed=NULL, *imreca_preprocessed=NULL;
  transf3d *preProcessingRefTrf=NULL, *preProcessingRecaTrf=NULL, *transfoini=NULL;
  InterpolationFct interpol;
  dist_func_t  distance;
  min_func_t minimisation;
  int    tdebut,tfin;
  tdebut=time(NULL);
  
  if (end_resol>start_resol)
  {
     PUT_WARN("ERROR in imx_matching_rigide_multi_3d_p\n"); return 0;
  }

 //---------mise en place des parametres du recalage-------/

 distance=imx_choose_distance_fct(dist_type);

 interpol=imx_choose_interpolation_fct(inter_type);

 minimisation=imx_choose_minimisation_fct(min_type,dist_type);

 // ------------------PREPROCESS----------------------
 {
  imx_matching_preprocessing_3d_p(imref, imreca, &imref_preprocessed, &imreca_preprocessed, NULL, &preProcessingRefTrf, &preProcessingRecaTrf, distance);
 }

 // ------------------INITIALISATION PAR TRANSFO DE DEPART----------------------
 {
  //on veut une transformation de depart compatible avec le recalage courant
  if ((inittrf)&&(inittrf->typetrans==RIGID3D))
  {
   transfoini=imx_trouver_transf_initiale_matching_3d_p(imref_preprocessed, imreca_preprocessed, inittrf);
  }
  else transfoini=imx_trouver_transf_initiale_matching_3d_p(imref_preprocessed, imreca_preprocessed, NULL);
 }

 //------------------ RECALAGE -----------------------/
 //recalage en multiresolution sans changement de degre de liberte
 imx_matching_multires_simple_3d_p(imref_preprocessed, imreca_preprocessed, distance, interpol, minimisation, transfres, transfoini, matchPrecision, NULL, NULL, start_resol, end_resol);

 //on remet le tout dans le repere des images originales
 imx_copie_dim_param_to_transf3d_p(imref, transfres);
 if (preProcessingRecaTrf) comb_transf_3d(preProcessingRecaTrf,transfres,transfres,0);
 if (preProcessingRefTrf) comb_transf_3d(transfres,preProcessingRefTrf,transfres,0);

 if (preProcessingRecaTrf) free_transf3d(preProcessingRecaTrf);
 if (preProcessingRefTrf) free_transf3d(preProcessingRefTrf);
 if (imref_preprocessed!=imref) free_grphic3d(imref_preprocessed);
 if (imreca_preprocessed!=imreca) free_grphic3d(imreca_preprocessed);

 // ------------------CALCUL DE L'IMAGE A AFFICHER ET DE LA TRANSFO A SAUVEGARDER----------------------
 imx_aff_param(transfres->typetrans, transfres->param);
 imx_calc_transfo_sav_image_res_3d_p(imref, imreca, imres, transfres, save_type);

 //liberations memoire
 if (transfoini) free_transf3d(transfoini);

 tfin=time(NULL);
 //affichage du temps de calcul total
 {
  imx_aff_temps_calcul(tdebut, tfin);
 }

 return(1);
}

/*******************************************************************************
**        matching_rigide_zoom_3d                                                 
**                                                                        
**       recalage rigide + zoom de deux images 3D (avec zoom)                                
*******************************************************************************/
void matching_rigide_zoom_3d()
{
  query_simple_matching_parameters_3d(imx_matching_rigide_zoom_3d);
}

/*******************************************************************************
**        matching_rigide_zoom_precision_3d
**
**       recalage rigide de deux images 3D (avec zoom)
*******************************************************************************/
void matching_rigide_zoom_precision_3d()
{
 query_simple_precision_matching_parameters_3d(imx_matching_rigide_zoom_3d);
}

/*******************************************************************************
**        imx_matching_rigide_zoom_3d                                             
**                                                                       
**       recalage rigide avec zoom de 2 images 3D                                   
*******************************************************************************/
int imx_matching_rigide_zoom_3d(int im_ref, int im_reca, int im_res, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres, char *nomfichtrf, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision)
{
 grphic3d *imref,*imreca,*imres;
 transf3d *transfres=NULL;
 transf3d *inittrf=NULL;

 imref=ptr_img_3d(im_ref);
 imreca=ptr_img_3d(im_reca);
 imres=ptr_img_3d(im_res);

 if (nomfichtrf !=NULL)
     inittrf=load_transf_3d(nomfichtrf);
 transfres = cr_transf3d_p(imref, NULL);

 imx_inimaxminpixel_3d_p(imref);
 imx_inimaxminpixel_3d_p(imreca);
 if ((imref->min_pixel<0)||(imreca->min_pixel<0))
 { PUT_WARN("L'image contient des valeurs negatives!!\nOn ajoute -val_min pour la recaler\n"); }

 init_log(NULL);
 aff_log("\n RECALAGE RIGIDE AVEC ZOOM ");
 aff_log("\n");
 aff_log("recalage de %s sur %s \n",(imreca->patient).name,(imref->patient).name);

 imx_matching_rigide_zoom_3d_p(imref,imreca,imres,dist_type,inter_type,min_type,save_type,transfres,inittrf,FieldOfResearch,matchPrecision);

 /*enregistrement du champ resultat*/
 if (nomfichres!=NULL)
    save_transf_3d(transfres,nomfichres);

 end_log();

 free_transf3d(transfres);
 
  return(1);
}


/*! \ingroup     Recalage  @{ */
/*******************************************************************************
**        imx_matching_rigide_zoom_3d_p
**                                         
**       \param :pour les parametres voir imx_matching_rigide_3d_p
**                                                              
**       recalage rigide avec zoom de 2 images 3D                                    
*******************************************************************************/
int imx_matching_rigide_zoom_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int dist_type, int inter_type, int min_type, int save_type, transf3d *transfres, transf3d *inittrf,eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision)
{
 int wdth,hght,dpth;
 grphic3d *imtref,*imtreca,*imtres;
 field3d *champ1;
 InterpolationFct interpol;
 dist_func_t distance;
 min_func_t minimisation;
 double  mini;
 int    tdebut,tfin;
 transf3d *preProcessingRefTrf=NULL, *preProcessingRecaTrf=NULL, *transfoini=NULL;
 
 wdth=imref->width;hght=imref->height;dpth=imref->depth;
#ifndef WIN32
 setvbuf(stdout, (char *)NULL, _IONBF, 0); //setbuf(stdout,NULL);
#endif
 tdebut=time(NULL);

 //---------mise en place des parametres du recalage-------/

  distance=imx_choose_distance_fct(dist_type);

  interpol=imx_choose_interpolation_fct(inter_type);

  minimisation=imx_choose_minimisation_fct(min_type,dist_type);

 // ------------------PREPROCESS----------------------
 {
  imx_matching_preprocessing_3d_p(imref, imreca, &imtref, &imtreca, &imtres, &preProcessingRefTrf, &preProcessingRecaTrf, distance);
 }

 // ------------------INITIALISATION PAR TRANSFO DE DEPART----------------------
 {
  //on veut une transformation de depart compatible avec le recalage courant
  if ((inittrf)&&(inittrf->typetrans==RIGID3D))
  {
   transfoini=imx_trouver_transf_initiale_matching_3d_p(imtref, imtreca, inittrf);
  }
  else transfoini=imx_trouver_transf_initiale_matching_3d_p(imtref, imtreca, NULL);
 }

 //------------------RECALAGE ----------------------------------------/
 {
 wdth=imtref->width;hght=imtref->height;dpth=imtref->depth;
 champ1=cr_field3d(wdth,hght,dpth);

 //recalage rigide
 aff_log("\n TRANSFO RIGIDE \n");
 {
  mini=imx_matching_simple_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, FieldOfResearch, matchPrecision, NULL, NULL, champ1);
 }

 //construction de la transfo initiale a injecter en debut du recalage rigide+zoom
 imx_aff_param(RIGID3D, transfres->param);
 copie_transf_3d(transfres, transfoini);
 transfoini->param=REALLOC(transfoini->param,transfoini->nb_param+3,double);
 transfoini->nb_param=rigid_to_rigidz_3d(transfoini->param);
 transfoini->typetrans=RIGIDZOOM3D;

 //recalage rigide+zoom
 aff_log("\n TRANSFO RIGIDE+ZOOM \n");
 {
  mini=imx_matching_simple_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, FieldOfResearch, matchPrecision, NULL, NULL, champ1);
 }

 free_field3d(champ1);
 }

 //on remet le tout dans le repere des images originales
 imx_copie_dim_param_to_transf3d_p(imref, transfres);
 if (preProcessingRecaTrf) comb_transf_3d(preProcessingRecaTrf,transfres,transfres,0);
 if (preProcessingRefTrf) comb_transf_3d(transfres,preProcessingRefTrf,transfres,0);

 if (preProcessingRecaTrf) free_transf3d(preProcessingRecaTrf);
 if (preProcessingRefTrf) free_transf3d(preProcessingRefTrf);
 if (imtres) free_grphic3d(imtres);
 if (imtref!=imref) free_grphic3d(imtref);
 if (imtreca!=imreca) free_grphic3d(imtreca);
 
 // ------------------CALCUL DE L'IMAGE A AFFICHER ET DE LA TRANSFO A SAUVEGARDER----------------------
 imx_aff_param(transfres->typetrans, transfres->param);
 imx_calc_transfo_sav_image_res_3d_p(imref, imreca, imres, transfres, save_type);

 //liberation memoire des variables dynamiques
 if (transfoini) free_transf3d(transfoini);

 tfin=time(NULL);
 //affichage du temps de calcul total
 {
  imx_aff_temps_calcul(tdebut, tfin);
 }
 end_log();

 return(1);  
}
/*! @} */

/*! 
     recalage rigide de deux images 3D (sans zoom)  multiresolution      
*/
void matching_rigide_zoom_multi_3d()
{
 query_multires_matching_parameters_3d(imx_matching_rigide_zoom_multi_3d);
}

/*!
     recalage rigide de deux images 3D (sans zoom)  multiresolution
*/
void matching_rigide_zoom_multi_precision_3d()
{
 query_multires_precision_matching_parameters_3d(imx_matching_rigide_zoom_multi_3d);
}

/*******************************************************************************
**        imx_matching_rigide_zoom_multi_3d                                             
*/                                                                       
/*!       recalage rigide multiresolution avec zoom de 2 images 3D                                  
** 
**  \param im_ref : numero de l'image de reference, 
**  \param im_reca  : numero de l'image a recaler,
**  \param im_res : numero de l'image resultat
**  les autres parametres sont expliques dans ::imx_matching_rigide_3d_p    
**  \retval 1
**  
*******************************************************************************/
int imx_matching_rigide_zoom_multi_3d(int im_ref, int im_reca, int im_res, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres,char *nomfichtrf, eMatchPrecision matchPrecision, int start_resol,int end_resol)
{
   grphic3d *imref=NULL, *imreca=NULL, *imres=NULL;
   transf3d *inittrf=NULL;
   transf3d *transfres=NULL;
   
   imref=ptr_img_3d(im_ref);
   imreca=ptr_img_3d(im_reca);
   imres=ptr_img_3d(im_res);

  if (nomfichtrf !=NULL) // initialisation de la transfo par un fichier.
     inittrf=load_transf_3d(nomfichtrf);
  transfres=cr_transf3d_p(imref,NULL);

  imx_inimaxminpixel_3d_p(imref);
  imx_inimaxminpixel_3d_p(imreca);
  if ((imref->min_pixel<0)||(imreca->min_pixel<0))
  { PUT_WARN("L'image contient des valeurs negatives!!\nOn ajoute -val_min pour la recaler\n"); }

  init_log(NULL);
  aff_log("\n RECALAGE RIGIDE AVEC ZOOM EN MULTIRESOLUTION");
  aff_log("\n");
  aff_log("recalage de %s sur %s \n",(imreca->patient).name,(imref->patient).name);
 
  imx_matching_rigide_zoom_multi_3d_p(imref,imreca,imres,dist_type,inter_type,min_type,save_type,transfres,inittrf,matchPrecision, start_resol, end_resol);

   /*enregistrement du champ resultat*/
   if (nomfichres!=NULL)
   {
     save_transf_3d(transfres,nomfichres);
   } 

   if (inittrf)
     free_transf3d(inittrf);

   free_transf3d(transfres);

   end_log();

   return(1);
}

/*******************************************************************************
**        imx_matching_rigide_zoom_multi_3d                                             
*/                                                                       
/*!       recalage rigide multiresolution avec zoom de 2 images 3D                                  
** 
**  \param im_ref : numero de l'image de reference, 
**  \param im_reca  : numero de l'image a recaler,
**  \param im_res : numero de l'image resultat
**  les autres parametres sont expliques dans ::imx_matching_rigide_3d_p    
**  \retval 1
**  
*******************************************************************************/
int imx_matching_rigide_zoom_multi_3d_p( grphic3d *imref, 
                                    grphic3d *imreca, 
                                    grphic3d *imres, 
                                    int dist_type, 
                                    int inter_type, 
                                    int min_type, 
                                    int save_type, 
                                    transf3d *transfres, 
                                    transf3d *inittrf,
                                    eMatchPrecision matchPrecision, 
                                    int start_resol,
                                    int end_resol)
{
 int j;
 grphic3d *imref_preprocessed=NULL, *imreca_preprocessed=NULL;
 transf3d *preProcessingRefTrf=NULL, *preProcessingRecaTrf=NULL, *transfoini=NULL;
 InterpolationFct interpol;
 dist_func_t  distance;
 min_func_t minimisation;
 int    tdebut,tfin;
 tdebut=time(NULL);

 if (end_resol>start_resol)
 {
   PUT_WARN("ERROR in imx_matching_rigide_zoom_multi_3d_p\n"); return 0;
 }

 //---------mise en place des parametres du recalage-------/

 distance=imx_choose_distance_fct(dist_type);

 interpol=imx_choose_interpolation_fct(inter_type);

 minimisation=imx_choose_minimisation_fct(min_type,dist_type);

 // ------------------PREPROCESS----------------------
 {
  imx_matching_preprocessing_3d_p(imref, imreca, &imref_preprocessed, &imreca_preprocessed, NULL, &preProcessingRefTrf, &preProcessingRecaTrf, distance);
 }

 // ------------------INITIALISATION PAR TRANSFO DE DEPART----------------------
 {
  //on veut une transformation de depart compatible avec le recalage courant
  if ((inittrf)&&(inittrf->typetrans==RIGID3D))
  {
   transfoini=imx_trouver_transf_initiale_matching_3d_p(imref_preprocessed, imreca_preprocessed, inittrf);
  }
  else transfoini=imx_trouver_transf_initiale_matching_3d_p(imref_preprocessed, imreca_preprocessed, NULL);
 }

 //------------------ RECALAGE -----------------------/
 for (j=0;j<2;j++)
 {//boucle sur les methodes de recalages
  switch (j)
  {
      case 0 : aff_log("\nRECALAGE RIGIDE\n");break;
      case 1 : aff_log("\nRECALAGE RIGIDE + ZOOM\n");break;
      default : fprintf(stderr,"not possible\n"); return 0;
  }

  //recalage en multiresolution sans changement de degre de liberte
  imx_matching_multires_simple_3d_p(imref_preprocessed, imreca_preprocessed, distance, interpol, minimisation, transfres, transfoini, matchPrecision, NULL, NULL, start_resol, end_resol);

  //construction de la transfo initiale a injecter en debut du recalage suivant
  if (j==0)
  {
   transfoini->param=REALLOC(transfoini->param,transfoini->nb_param+3,double);
   transfoini->nb_param=rigid_to_rigidz_3d(transfoini->param);
   transfoini->typetrans=RIGIDZOOM3D;
  }
 }

 //on remet le tout dans le repere des images originales
 imx_copie_dim_param_to_transf3d_p(imref, transfres);
 if (preProcessingRecaTrf) comb_transf_3d(preProcessingRecaTrf,transfres,transfres,0);
 if (preProcessingRefTrf) comb_transf_3d(transfres,preProcessingRefTrf,transfres,0);

 if (preProcessingRecaTrf) free_transf3d(preProcessingRecaTrf);
 if (preProcessingRefTrf) free_transf3d(preProcessingRefTrf);
 if (imref_preprocessed!=imref) free_grphic3d(imref_preprocessed);
 if (imreca_preprocessed!=imreca) free_grphic3d(imreca_preprocessed);
 
 // ------------------CALCUL DE L'IMAGE A AFFICHER ET DE LA TRANSFO A SAUVEGARDER----------------------
 imx_aff_param(transfres->typetrans, transfres->param);
 imx_calc_transfo_sav_image_res_3d_p(imref, imreca, imres, transfres, save_type);

 //liberations memoire
 if (transfoini) free_transf3d(transfoini);

 tfin=time(NULL);
 //affichage du temps de calcul total
 {
  imx_aff_temps_calcul(tdebut, tfin);
 }

  return(1);
}

/*******************************************************************************
********************************************************************************
************************** MATCHING  AFFINE ************************************
********************************************************************************
*******************************************************************************/

/*******************************************************************************
**        matching_affine_3d                                                 
*/
/*!                                                                        
**      \brief recalage affine de deux images 3D                                
*******************************************************************************/
void matching_affine_3d()
{
  query_simple_matching_parameters_3d(imx_matching_affine_3d);
}

/*******************************************************************************
**        matching_affine_precision_3d
*/
/*!
**      \brief recalage affine de deux images 3D
*******************************************************************************/
void matching_affine_precision_3d()
{
 query_simple_precision_matching_parameters_3d(imx_matching_affine_3d);
}


/*******************************************************************************
**        imx_matching_affine_3d                                             
*/
/*!     \brief  recalage affine de 2 images 3D                                    
**
**      \param  imref : numero de l'image de ref
**      \param  imreca : numero de l'image a recaler
**      \param  imres : numero de l'image ou l'on stocke le resultat du recalage de imreca
    (peut etre NULL)
**      \param  dist_type : numero de la fonction de distance (0 pour distance quadratique)
                            \li 0 : QUAD
                            \li 1 : QUAD ROBUST
                            \li 2 : WOODS
                            \li 3 : WOODS ROBUST
                            \li 4 : IM
                            \li 
                            
**      \param  inter_type : numero de la fonction d'interpolation des images (1 pour lineaire)
                            \li 0 : NEAREST
                            \li 1 : LINEAR
                            \li 2 : SINC
                            \li 3 : QSINC2
                            \li 4 : QSINC3
                            \li 5 : BSPLINE2
                            \li 6 : BSPLINE3 
                            \li 7 : BSPLINE4
                            \li 8 : BSPLINE5


**      \param  min_type : numero de la fonction de minimisation utilisee aux plus faibles
    echelles de la transformation (5 pour quasi-newton)
                            \li 0 : ICM
                            \li 1 : SIMPLEX 
                            \li 2 : DESCGRAD 
                            \li 3 : GRADCONJ
                            \li 4 : DESCGRADMOD
                            \li 5 : QUASINEWTONMOD

**      \param save_type : 
**                          \li 0 : champ
**                          \li 1 : parametres
**      \param  nomfichres :  nom du fichier de sauvegarde de la transformation
**      \retval 1
*******************************************************************************/
int imx_matching_affine_3d(int im_ref, int im_reca, int im_res, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres, char *nomfichtrf, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision)
{
 grphic3d *imref,*imreca,*imres;
 transf3d *inittrf=NULL;
 transf3d *transfres=NULL;
 imref=ptr_img_3d(im_ref);
 imreca=ptr_img_3d(im_reca);
 imres=ptr_img_3d(im_res);

 if (nomfichtrf !=NULL)
     inittrf=load_transf_3d(nomfichtrf);  
 transfres=cr_transf3d_p(imref,NULL);

 imx_inimaxminpixel_3d_p(imref);
 imx_inimaxminpixel_3d_p(imreca);
 if ((imref->min_pixel<0)||(imreca->min_pixel<0))
 { PUT_WARN("L'image contient des valeurs negatives!!\nOn ajoute -val_min pour la recaler\n"); }

 init_log(NULL);
 aff_log("\n RECALAGE AFFINE ");
 aff_log("\n");
 aff_log("recalage de %s sur %s \n",(imreca->patient).name,(imref->patient).name);

 imx_matching_affine_3d_p(imref,imreca,imres,dist_type,inter_type,min_type,save_type,transfres,inittrf,FieldOfResearch,matchPrecision);

 /*enregistrement du champ resultat*/
 if (nomfichres!=NULL)
 { 
   save_transf_3d(transfres,nomfichres);
 } 

 end_log();
/*   if (inittrf)
     free_transf3d(inittrf);
*/
  free_transf3d(transfres);
  return(1);
}

/*! \ingroup     Recalage @{ */
/*******************************************************************************
**        imx_matching_affine_3d_p                                            
**                                                                        
*/
/*!     \brief  recalage affine de 2 images 3D                                    
**
**      \param  imref : image de reference
**      \param  imreca : image a recaler
**      \param  imres : image resultat (E/S)
**      \param  dist_type : numero de la fonction de distance 
                            \li 0 : QUAD
                            \li 1 : QUAD ROBUST
                            \li 2 : WOODS
                            \li 3 : WOODS ROBUST
                            \li 4 : IM
                            \li 
                            
**      \param inter_type : numero de la fonction d'interpolation des images 
                            \li 0 : NEAREST
                            \li 1 : LINEAR
                            \li 2 : SINC
                            \li 3 : QSINC2
                            \li 4 : QSINC3
                            \li 5 : BSPLINE2
                            \li 6 : BSPLINE3 
                            \li 7 : BSPLINE4
                            \li 8 : BSPLINE5

**      \param  min_type : numero de la fonction de minimisation utilisee aux plus faibles
    echelles de la transformation 
                            \li 0 : ICM
                            \li 1 : SIMPLEX 
                            \li 2 : DESCGRAD 
                            \li 3 : GRADCONJ
                            \li 4 : DESCGRADMOD
                            \li 5 : QUASINEWTONMOD

**      \param save_type : 
**                          \li 0 : champ
**                          \li 1 : parametres
**      \param  nomfichres :   nom du fichier de sauvegarde de la transfo (*.trf)
**  \param FieldOfResearch :  voir init_bornes_matching
**  \param matchPrecision  :  voir init_bornes_matching
**      \retval 1
*******************************************************************************/
int imx_matching_affine_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int dist_type, int inter_type, int min_type, int save_type, transf3d *transfres,transf3d *inittrf, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision)
{
 int    wdth,hght,dpth;
 grphic3d *imtref,*imtreca,*imtres;
 field3d *champ1;
 InterpolationFct interpol;
 dist_func_t distance;
 min_func_t minimisation;
 double  mini;
 int    tdebut,tfin;
 transf3d *preProcessingRefTrf=NULL, *preProcessingRecaTrf=NULL, *transfoini=NULL;

 wdth=imref->width;hght=imref->height;dpth=imref->depth;
#ifndef WIN32
 setvbuf(stdout, (char *)NULL, _IONBF, 0); //setbuf(stdout,NULL);
#endif
 tdebut=time(NULL);

 //---------mise en place des parametres du recalage-------/

  distance=imx_choose_distance_fct(dist_type);

  interpol=imx_choose_interpolation_fct(inter_type);

  minimisation=imx_choose_minimisation_fct(min_type,dist_type);

 // ------------------PREPROCESS----------------------
 {
  imx_matching_preprocessing_3d_p(imref, imreca, &imtref, &imtreca, &imtres, &preProcessingRefTrf, &preProcessingRecaTrf, distance);
 }

 // ------------------INITIALISATION PAR TRANSFO DE DEPART----------------------
 {
  //on veut une transformation de depart compatible avec le recalage courant
  if ((inittrf)&&(inittrf->typetrans==RIGID3D))
  {
   transfoini=imx_trouver_transf_initiale_matching_3d_p(imtref, imtreca, inittrf);
  }
  else transfoini=imx_trouver_transf_initiale_matching_3d_p(imtref, imtreca, NULL);
 }

 //------------------RECALAGE ----------------------------------------/
 {
 wdth=imtref->width;hght=imtref->height;dpth=imtref->depth;
 champ1=cr_field3d(wdth,hght,dpth);

 //recalage rigide
 aff_log("\n TRANSFO RIGIDE \n");
 {
  mini=imx_matching_simple_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, FieldOfResearch, matchPrecision, NULL, NULL, champ1);
 }

 //construction de la transfo initiale a injecter en debut du recalage rigide+zoom
 imx_aff_param(RIGID3D, transfres->param);
 copie_transf_3d(transfres, transfoini);
 transfoini->param=REALLOC(transfoini->param,transfoini->nb_param+3,double);
 transfoini->nb_param=rigid_to_rigidz_3d(transfoini->param);
 transfoini->typetrans=RIGIDZOOM3D;

 //recalage rigide+zoom
 aff_log("\n TRANSFO RIGIDE+ZOOM \n");
 {
  mini=imx_matching_simple_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, FieldOfResearch, matchPrecision, NULL, NULL, champ1);
 }

 //construction de la transfo initiale a injecter en debut du recalage affine
 imx_aff_param(RIGIDZOOM3D, transfres->param);
 copie_transf_3d(transfres, transfoini);
 transfoini->param=REALLOC(transfoini->param,transfoini->nb_param+3,double);
 //transfoini->nb_param=rigidz_to_affine_3d(transfoini->param);
 //transfoini->typetrans=AFFINE3D;
 transfoini->nb_param=rigidz_to_affine_decouple_3d(transfoini->param);
 transfoini->typetrans=AFFINEDECOUPLE3D;

 //recalage affine
 aff_log("\n TRANSFO AFFINE \n");
 {
  mini=imx_matching_simple_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, FieldOfResearch, matchPrecision, NULL, NULL, champ1);
 }

 free_field3d(champ1);
 }

 //on remet le tout dans le repere des images originales
 transfres->nb_param=affine_decouple_to_affine_3d(transfres->param);
 transfres->typetrans=AFFINE3D;
 imx_copie_dim_param_to_transf3d_p(imref, transfres);
 if (preProcessingRecaTrf) comb_transf_3d(preProcessingRecaTrf,transfres,transfres,0);
 if (preProcessingRefTrf) comb_transf_3d(transfres,preProcessingRefTrf,transfres,0);

 if (preProcessingRecaTrf) free_transf3d(preProcessingRecaTrf);
 if (preProcessingRefTrf) free_transf3d(preProcessingRefTrf);
 if (imtres) free_grphic3d(imtres);
 if (imtref!=imref) free_grphic3d(imtref);
 if (imtreca!=imreca) free_grphic3d(imtreca);
 
 // ------------------CALCUL DE L'IMAGE A AFFICHER ET DE LA TRANSFO A SAUVEGARDER----------------------
 imx_aff_param(transfres->typetrans, transfres->param);
 imx_calc_transfo_sav_image_res_3d_p(imref, imreca, imres, transfres, save_type);

 //liberation memoire des variables dynamiques
 if (transfoini) free_transf3d(transfoini);

 tfin=time(NULL);
 //affichage du temps de calcul total
 {
  imx_aff_temps_calcul(tdebut, tfin);
 }
 end_log();

 return(1);  
}
/*! @} */
 /*!
     recalage affine de deux images 3D  multiresolution
*/
void matching_affine_multi_3d()
{
 query_multires_matching_parameters_3d(imx_matching_affine_multi_3d);
}

 /*!
     recalage affine de deux images 3D  multiresolution
*/
void matching_affine_multi_precision_3d()
{
 query_multires_precision_matching_parameters_3d(imx_matching_affine_multi_3d);
}

/*******************************************************************************
**        imx_matching_affine_multi_3d
*/
/*!       recalage affine multiresolution  de 2 images 3D
**

**  \param im_ref : numero de l'image de reference,
**  \param im_reca  : numero de l'image a recaler,
**  \param im_res : numero de l'image resultat
**  les autres parametres sont expliques dans ::imx_matching_affine_3d_p
**  \retval 1
**
*******************************************************************************/
int imx_matching_affine_multi_3d(int im_ref, int im_reca, int im_res, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres,char *nomfichtrf, eMatchPrecision matchPrecision, int start_resol, int end_resol)
{
   grphic3d *imref=NULL, *imreca=NULL, *imres=NULL;
   transf3d *inittrf=NULL;
   transf3d *transfres=NULL;

   imref=ptr_img_3d(im_ref);
   imreca=ptr_img_3d(im_reca);

   imres=ptr_img_3d(im_res);

   if (nomfichtrf !=NULL) // initialisation de la transfo par un fichier.
     inittrf=load_transf_3d(nomfichtrf);
   transfres=cr_transf3d_p(imref,NULL);

   imx_inimaxminpixel_3d_p(imref);
   imx_inimaxminpixel_3d_p(imreca);
   if ((imref->min_pixel<0)||(imreca->min_pixel<0))
   { PUT_WARN("L'image contient des valeurs negatives!!\nOn ajoute -val_min pour la recaler\n"); }

   init_log(NULL);
   aff_log("\n RECALAGE AFFINE EN MULTIRESOLUTION");
   aff_log("\n");
   aff_log("recalage de %s sur %s \n",(imreca->patient).name,(imref->patient).name);

   imx_matching_affine_multi_3d_p(imref,imreca,imres,dist_type,inter_type,min_type,save_type,transfres,inittrf,matchPrecision, start_resol, end_resol);

   /*enregistrement du champ resultat*/
   if (nomfichres!=NULL)
   {
     save_transf_3d(transfres,nomfichres);
   }

   if (inittrf)
     free_transf3d(inittrf);
   free_transf3d(transfres);

   end_log();
   return(1);
}

/*******************************************************************************
**        imx_matching_affine_multi_3d
*/
/*!       recalage affine multiresolution  de 2 images 3D
**
**  \param im_ref : numero de l'image de reference,
**  \param im_reca  : numero de l'image a recaler,
**  \param im_res : numero de l'image resultat
**  les autres parametres sont expliques dans ::imx_matching_rigide_3d_p
**  \retval 1
**
*******************************************************************************/
int imx_matching_affine_multi_3d_p( grphic3d *imref,
                                    grphic3d *imreca,
                                    grphic3d *imres,
                                    int dist_type,
                                    int inter_type,
                                    int min_type,
                                    int save_type,
                                    transf3d *transfres,
                                    transf3d *inittrf,
                                    eMatchPrecision matchPrecision,
                                    int start_resol, 
                                    int end_resol
                                    )
{
 int j;
 grphic3d *imref_preprocessed=NULL, *imreca_preprocessed=NULL;
 transf3d *preProcessingRefTrf=NULL, *preProcessingRecaTrf=NULL, *transfoini=NULL;
 InterpolationFct interpol;
 dist_func_t  distance;
 min_func_t minimisation;
 int    tdebut,tfin;
 tdebut=time(NULL);

 if (end_resol>start_resol)
 {
   PUT_WARN("ERROR in imx_matching_affine_multi_3d_p\n"); return 0;
 }

 //---------mise en place des parametres du recalage-------/

 distance=imx_choose_distance_fct(dist_type);

 interpol=imx_choose_interpolation_fct(inter_type);

 minimisation=imx_choose_minimisation_fct(min_type,dist_type);

 // ------------------PREPROCESS----------------------
 {
  imx_matching_preprocessing_3d_p(imref, imreca, &imref_preprocessed, &imreca_preprocessed, NULL, &preProcessingRefTrf, &preProcessingRecaTrf, distance);
 }

 // ------------------INITIALISATION PAR TRANSFO DE DEPART----------------------
 {
  //on veut une transformation de depart compatible avec le recalage courant
  if ((inittrf)&&(inittrf->typetrans==RIGID3D))
  {
   transfoini=imx_trouver_transf_initiale_matching_3d_p(imref_preprocessed, imreca_preprocessed, inittrf);
  }
  else transfoini=imx_trouver_transf_initiale_matching_3d_p(imref_preprocessed, imreca_preprocessed, NULL);
 }

 //------------------ RECALAGE -----------------------/
 for (j=0;j<3;j++)
 {//boucle sur les methodes de recalages
  switch (j)
  {
      case 0 : aff_log("\nRECALAGE RIGIDE\n");break;
      case 1 : aff_log("\nRECALAGE RIGIDE + ZOOM\n");break;
      case 2 : aff_log("\nRECALAGE AFFINE\n");break;
      default : fprintf(stderr,"not possible\n"); return 0;
  }

  //recalage en multiresolution sans changement de degre de liberte
  imx_matching_multires_simple_3d_p(imref_preprocessed, imreca_preprocessed, distance, interpol, minimisation, transfres, transfoini, matchPrecision, NULL, NULL, start_resol, end_resol);

  //construction de la transfo initiale a injecter en debut du recalage suivant
  switch (j)
  {
      case 0 : transfoini->param=REALLOC(transfoini->param,transfoini->nb_param+3,double);
               transfoini->nb_param=rigid_to_rigidz_3d(transfoini->param);
               transfoini->typetrans=RIGIDZOOM3D;
               break;
      case 1 : transfoini->param=REALLOC(transfoini->param,transfoini->nb_param+3,double);
               //transfoini->nb_param=rigidz_to_affine_3d(transfoini->param);
               //transfoini->typetrans=AFFINE3D;
               transfoini->nb_param=rigidz_to_affine_decouple_3d(transfoini->param);
               transfoini->typetrans=AFFINEDECOUPLE3D;
               break;
      case 2 : break;
      default : fprintf(stderr,"not possible\n"); return 0;
  }
  
 }

 //on remet le tout dans le repere des images originales
 transfres->nb_param=affine_decouple_to_affine_3d(transfres->param);
 transfres->typetrans=AFFINE3D;
 imx_copie_dim_param_to_transf3d_p(imref, transfres);
 if (preProcessingRecaTrf) comb_transf_3d(preProcessingRecaTrf,transfres,transfres,0);
 if (preProcessingRefTrf) comb_transf_3d(transfres,preProcessingRefTrf,transfres,0);

 if (preProcessingRecaTrf) free_transf3d(preProcessingRecaTrf);
 if (preProcessingRefTrf) free_transf3d(preProcessingRefTrf);
 if (imref_preprocessed!=imref) free_grphic3d(imref_preprocessed);
 if (imreca_preprocessed!=imreca) free_grphic3d(imreca_preprocessed);

 // ------------------CALCUL DE L'IMAGE A AFFICHER ET DE LA TRANSFO A SAUVEGARDER----------------------
 imx_aff_param(transfres->typetrans, transfres->param);
 imx_calc_transfo_sav_image_res_3d_p(imref, imreca, imres, transfres, save_type);

 //liberations memoire
 if (transfoini) free_transf3d(transfoini);
 
 tfin=time(NULL);
 //affichage du temps de calcul total
 {
  imx_aff_temps_calcul(tdebut, tfin);
 }

  return(1);
}


/*******************************************************************************
********************************************************************************
****************** MATCHING  RIGIDE MULTISTART MULTIRES ************************
********************************************************************************

*******************************************************************************/

/*******************************************************************************
**        matching_rigide_multistart_multires_3d
*/
/*!
**      \brief recalage rigide multistart et multiresolution de deux images 3D
*******************************************************************************/
extern void     matching_rigide_multistart_multires_3d()
{
  query_simple_matching_parameters_3d(imx_matching_rigide_multistart_multires_3d);
}

/*******************************************************************************
**        matching_rigide_multistart_multires_precision_3d
*/
/*!
**      \brief recalage rigide multistart et multiresolution de deux images 3D
*******************************************************************************/
extern void     matching_rigide_multistart_multires_precision_3d()
{
 query_simple_precision_matching_parameters_3d(imx_matching_rigide_multistart_multires_3d);
}

/*******************************************************************************
**        imx_matching_rigide_multistart_multires_3d
*/
/*!       recalage rigide multistart et multiresolution  de 2 images 3D
**
**  \param im_ref : numero de l'image de reference,
**  \param im_reca  : numero de l'image a recaler,
**  \param im_res : numero de l'image resultat
**  les autres parametres sont expliques dans ::imx_matching_affine_3d_p
**  \retval 1
**
*******************************************************************************/
extern int  imx_matching_rigide_multistart_multires_3d(int im_ref, int im_reca, int im_res, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres, char *nomfichtrf,eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision)
{
  grphic3d *imref,*imreca,*imres;
  transf3d *inittrf=NULL;
  transf3d *transfres=NULL;

  imref=ptr_img_3d(im_ref);
  imreca=ptr_img_3d(im_reca);
  imres=ptr_img_3d(im_res);

  if (nomfichtrf !=NULL)
     inittrf=load_transf_3d(nomfichtrf);
  transfres = cr_transf3d_p(imref,NULL);

  imx_inimaxminpixel_3d_p(imref);
  imx_inimaxminpixel_3d_p(imreca);
  if ((imref->min_pixel<0)||(imreca->min_pixel<0))
  { PUT_WARN("L'image contient des valeurs negatives!!\nOn ajoute -val_min pour la recaler\n"); }

  init_log(NULL);
  aff_log("\n RECALAGE RIGIDE EN MULTIRESOLUTION AVEC MULTISTART (methode Jenkison&Smith) ");
  aff_log("\n");
  aff_log("recalage de %s sur %s \n",(imreca->patient).name,(imref->patient).name);


  imx_matching_rigide_multistart_multires_3d_p(imref,imreca,imres,dist_type,inter_type,min_type,save_type,transfres,inittrf,FieldOfResearch,matchPrecision);

  /*enregistrement du champ resultat*/
  if (nomfichres!=NULL)
  {
    save_transf_3d(transfres,nomfichres);
  }

  end_log();

  if (inittrf)
    free_transf3d(inittrf);

  free_transf3d(transfres);

  return(0);
}

/*******************************************************************************
**        imx_matching_rigide_multistart_multires_3d_p
*/
/*!       recalage rigide multistart et multiresolution  de 2 images 3D
**
**  \param im_ref : numero de l'image de reference,
**  \param im_reca  : numero de l'image a recaler,
**  \param im_res : numero de l'image resultat
**  les autres parametres sont expliques dans ::imx_matching_rigide_3d_p
**  \retval 1
**
*******************************************************************************/
extern int imx_matching_rigide_multistart_multires_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int dist_type, int inter_type, int min_type, int save_type, transf3d *transfres, transf3d *inittrf,eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision)
{
 int k;
 grphic3d *imtref,*imtreca,*imtres, *imref_preprocessed, *imreca_preprocessed, *imref_isotrope;
 field3d *champ_fin;
 InterpolationFct interpol;
 dist_func_t  distance;
 min_func_t minimisation;
 int    tdebut,tfin;
 int nb_minima=3, nb_minima_max=3;
 double **minima;
 transf3d *preProcessingRefTrf=NULL, *preProcessingRecaTrf=NULL, *transfoini=NULL;

 minima=CALLOC(nb_minima_max, double*);
 for (k=0;k<nb_minima_max;k++) minima[k]=CALLOC(15,double);

#ifndef WIN32
 setvbuf(stdout, (char *)NULL, _IONBF, 0); //setbuf(stdout,NULL);
#endif
 tdebut=time(NULL);

 //---------mise en place des parametres du recalage-------/

 distance=imx_choose_distance_fct(dist_type);

 interpol=imx_choose_interpolation_fct(inter_type);

 minimisation=imx_choose_minimisation_fct(min_type,dist_type);

 // ------------------PREPROCESS----------------------
 {
  imx_matching_preprocessing_3d_p(imref, imreca, &imref_preprocessed, &imreca_preprocessed, NULL, &preProcessingRefTrf, &preProcessingRecaTrf, distance);
 }

 // ------------------INITIALISATION PAR TRANSFO DE DEPART----------------------
 {
  //on veut une transformation de depart compatible avec le recalage courant
  if ((inittrf)&&(inittrf->typetrans==RIGID3D))
  {
   transfoini=imx_trouver_transf_initiale_matching_3d_p(imref_preprocessed, imreca_preprocessed, inittrf);
  }
  else transfoini=imx_trouver_transf_initiale_matching_3d_p(imref_preprocessed, imreca_preprocessed, NULL);
 }

 //calcule une image isotrope en resolution 2 mm a partir de l'image de depart
 imref_isotrope=cr_grphic3d(imref_preprocessed);
 imx_creer_image_isotrope_3d_p(imref_preprocessed, inter_qsinc3_3d, imref_isotrope);

 imtref=cr_grphic3d(imref_isotrope);
 imtreca=cr_grphic3d(imreca_preprocessed);

 //------------------RECALAGE EN RESOLUTION 8 MM----------------------------------------/

 pyramide_gaussienne_recalage_3d_p(imref_isotrope, imreca_preprocessed, imtref, imtreca, 4.0);

 aff_log("\n minimisation sur image 32x32x32 \n");
 {
  champ_fin=cr_field3d(imtref->width,imtref->height,imtref->depth);
  imtres=cr_grphic3d(imtref);
  imx_copie_dim_param_to_transf3d_p(imtref, transfoini);

  imx_matching_multistart_multires_8mm_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, LARGE, matchPrecision, NULL, NULL, champ_fin, &nb_minima, minima);

  free_grphic3d(imtres);
  free_field3d(champ_fin);
 }

 imx_aff_param(RIGID3D, transfres->param);
 copie_transf_3d(transfres, transfoini);

 //------------------RECALAGE EN RESOLUTION 4 MM----------------------------------------/

 pyramide_gaussienne_recalage_3d_p(imref_isotrope, imreca_preprocessed, imtref, imtreca, 2.0);

 aff_log("\n minimisation sur image 64x64x64 \n");
 {
  champ_fin=cr_field3d(imtref->width,imtref->height,imtref->depth);
  imtres=cr_grphic3d(imtref);
  imx_copie_dim_param_to_transf3d_p(imtref, transfoini);

  imx_matching_multistart_multires_4mm_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, LARGE, matchPrecision, NULL, NULL, champ_fin, nb_minima, minima);

  free_grphic3d(imtres);
  free_field3d(champ_fin);
 }

 imx_aff_param(RIGID3D, transfres->param);
 copie_transf_3d(transfres, transfoini);

 //------------------RECALAGE EN RESOLUTION 2 MM----------------------------------------/

 pyramide_gaussienne_recalage_3d_p(imref_isotrope, imreca_preprocessed, imtref, imtreca, 1.0);

 aff_log("\n minimisation sur image 128x128x128 \n");
 {
  champ_fin=cr_field3d(imtref->width,imtref->height,imtref->depth);
  imtres=cr_grphic3d(imtref);
  imx_copie_dim_param_to_transf3d_p(imtref, transfoini);

  imx_matching_simple_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, SMALL, matchPrecision, NULL, NULL, champ_fin);

  free_grphic3d(imtres);
  free_field3d(champ_fin);
 }

 imx_aff_param(RIGID3D, transfres->param);
 copie_transf_3d(transfres, transfoini);

 free_grphic3d(imtref); free_grphic3d(imtreca);

 //------------------RECALAGE EN RESOLUTION INITIALE----------------------------------------/

 aff_log("\n minimisation sur image originale \n");
 {
  imtref=imref_preprocessed; imtreca=imreca_preprocessed; imtres=cr_grphic3d(imtref);
  champ_fin=cr_field3d(imtref->width,imtref->height,imtref->depth);
  imx_copie_dim_param_to_transf3d_p(imtref, transfoini);

  imx_matching_simple_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, SMALL, matchPrecision, NULL, NULL, champ_fin);

  free_field3d(champ_fin);
  free_grphic3d(imtres);
 }

 //on remet le tout dans le repere des images originales
 imx_copie_dim_param_to_transf3d_p(imref, transfres);
 if (preProcessingRecaTrf) comb_transf_3d(preProcessingRecaTrf,transfres,transfres,0);
 if (preProcessingRefTrf) comb_transf_3d(transfres,preProcessingRefTrf,transfres,0);

 if (preProcessingRecaTrf) free_transf3d(preProcessingRecaTrf);
 if (preProcessingRefTrf) free_transf3d(preProcessingRefTrf);
 if (imref_preprocessed!=imref) free_grphic3d(imref_preprocessed);
 if (imreca_preprocessed!=imreca) free_grphic3d(imreca_preprocessed);
 free_grphic3d(imref_isotrope);
  
 // ------------------CALCUL DE L'IMAGE A AFFICHER ET DE LA TRANSFO A SAUVEGARDER----------------------
 imx_aff_param(transfres->typetrans, transfres->param);
 imx_calc_transfo_sav_image_res_3d_p(imref, imreca, imres, transfres, save_type);

 //liberations memoire
 if (transfoini) free_transf3d(transfoini);
 for (k=0;k<nb_minima_max;k++) FREE(minima[k]);
 FREE(minima);

 tfin=time(NULL);
 //affichage du temps de calcul total
 {
  imx_aff_temps_calcul(tdebut, tfin);
 }

 return 0;
}

/*******************************************************************************
********************************************************************************
****************** MATCHING  RIGIDE ZOOM MULTISTART MULTIRES ************************
********************************************************************************

*******************************************************************************/

/*******************************************************************************
**        matching_rigide_zoom_multistart_multires_3d
*/
/*!
**      \brief recalage rigide + zoom multistart et multiresolution de deux images 3D
*******************************************************************************/
extern void     matching_rigide_zoom_multistart_multires_3d()
{
  query_simple_matching_parameters_3d(imx_matching_rigide_zoom_multistart_multires_3d);
}

/*******************************************************************************
**        matching_rigide_zoom_multistart_multires_precision_3d
*/
/*!
**      \brief recalage rigide + zoom multistart et multiresolution de deux images 3D
*******************************************************************************/
extern void     matching_rigide_zoom_multistart_multires_precision_3d()
{
 query_simple_precision_matching_parameters_3d(imx_matching_rigide_zoom_multistart_multires_3d);
}

/*******************************************************************************
**        imx_matching_rigide_zoom_multistart_multires_3d
*/
/*!       recalage rigide +zoom multistart et multiresolution  de 2 images 3D
**
**  \param im_ref : numero de l'image de reference,
**  \param im_reca  : numero de l'image a recaler,
**  \param im_res : numero de l'image resultat
**  les autres parametres sont expliques dans ::imx_matching_affine_3d_p
**  \retval 1
**
*******************************************************************************/
extern int  imx_matching_rigide_zoom_multistart_multires_3d(int im_ref, int im_reca, int im_res, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres, char *nomfichtrf,eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision)
{
  grphic3d *imref,*imreca,*imres;
  transf3d *inittrf=NULL;
  transf3d *transfres=NULL;

  imref=ptr_img_3d(im_ref);
  imreca=ptr_img_3d(im_reca);
  imres=ptr_img_3d(im_res);

  if (nomfichtrf !=NULL)
  inittrf=load_transf_3d(nomfichtrf);
  transfres = cr_transf3d_p(imref,NULL);

  imx_inimaxminpixel_3d_p(imref);
  imx_inimaxminpixel_3d_p(imreca);
  if ((imref->min_pixel<0)||(imreca->min_pixel<0))
  { PUT_WARN("L'image contient des valeurs negatives!!\nOn ajoute -val_min pour la recaler\n"); }

  init_log(NULL);
  aff_log("\n RECALAGE RIGIDE+ZOOM EN MULTIRESOLUTION AVEC MULTISTART (methode Jenkison&Smith) ");
  aff_log("\n");
  aff_log("recalage de %s sur %s \n",(imreca->patient).name,(imref->patient).name);

  imx_matching_rigide_zoom_multistart_multires_3d_p(imref,imreca,imres,dist_type,inter_type,min_type,save_type,transfres,inittrf,FieldOfResearch,matchPrecision);

  /*enregistrement du champ resultat*/
  if (nomfichres!=NULL)
  {
    save_transf_3d(transfres,nomfichres);
  }

  end_log();

  if (inittrf)
    free_transf3d(inittrf);

  free_transf3d(transfres);

  return(0);
}

/*******************************************************************************
**        imx_matching_rigide_zoom_multistart_multires_3d_p
*/
/*!       recalage rigide + zoom multistart et multiresolution  de 2 images 3D
**
**  \param im_ref : numero de l'image de reference,
**  \param im_reca  : numero de l'image a recaler,
**  \param im_res : numero de l'image resultat
**  les autres parametres sont expliques dans ::imx_matching_rigide_3d_p
**  \retval 1
**

*******************************************************************************/
extern int imx_matching_rigide_zoom_multistart_multires_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int dist_type, int inter_type, int min_type, int save_type, transf3d *transfres, transf3d *inittrf,eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision)
{
 int k;
 grphic3d *imtref,*imtreca,*imtres, *imref_preprocessed, *imreca_preprocessed, *imref_isotrope;
 field3d *champ_fin;
 InterpolationFct interpol;
 dist_func_t  distance;
 min_func_t minimisation;
 int tdebut,tfin;
 int nb_minima=3, nb_minima_max=3;
 double **minima;
 transf3d *preProcessingRefTrf=NULL, *preProcessingRecaTrf=NULL, *transfoini=NULL;

 minima=CALLOC(nb_minima_max, double*);
 for (k=0;k<nb_minima_max;k++) minima[k]=CALLOC(15,double);

#ifndef WIN32
 setvbuf(stdout, (char *)NULL, _IONBF, 0); //setbuf(stdout,NULL);
#endif
 tdebut=time(NULL);

 //---------mise en place des parametres du recalage-------/

  distance=imx_choose_distance_fct(dist_type);

  interpol=imx_choose_interpolation_fct(inter_type);

  minimisation=imx_choose_minimisation_fct(min_type,dist_type);

 // ------------------PREPROCESS----------------------
 {
  imx_matching_preprocessing_3d_p(imref, imreca, &imref_preprocessed, &imreca_preprocessed, NULL, &preProcessingRefTrf, &preProcessingRecaTrf, distance);
 }

 // ------------------INITIALISATION PAR TRANSFO DE DEPART----------------------
 {
  //on veut une transformation de depart compatible avec le recalage courant
  if ((inittrf)&&(inittrf->typetrans==RIGID3D))
  {
   transfoini=imx_trouver_transf_initiale_matching_3d_p(imref_preprocessed, imreca_preprocessed, inittrf);
  }
  else transfoini=imx_trouver_transf_initiale_matching_3d_p(imref_preprocessed, imreca_preprocessed, NULL);
 }
 //on se met en rigide+zoom global directement
 transfoini->param=REALLOC(transfoini->param,transfoini->nb_param+1,double);
 transfoini->nb_param=rigid_to_rigid_global_zoom_3d(transfoini->param);
 transfoini->typetrans=RIGIDGLOBALZOOM3D;

 //calcule une image isotrope en resolution 2 mm a partir de l'image de depart
 imref_isotrope=cr_grphic3d(imref_preprocessed);
 imx_creer_image_isotrope_3d_p(imref_preprocessed, inter_qsinc3_3d, imref_isotrope);

 imtref=cr_grphic3d(imref_isotrope);
 imtreca=cr_grphic3d(imreca_preprocessed);

 //------------------RECALAGE EN RESOLUTION 8 MM----------------------------------------/

 pyramide_gaussienne_recalage_3d_p(imref_isotrope, imreca_preprocessed, imtref, imtreca, 4.0);

 aff_log("\n minimisation sur image 32x32x32 \n");
 {
  champ_fin=cr_field3d(imtref->width,imtref->height,imtref->depth);
  imtres=cr_grphic3d(imtref);
  imx_copie_dim_param_to_transf3d_p(imtref, transfoini);

  imx_matching_multistart_multires_8mm_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, LARGE, matchPrecision, NULL, NULL, champ_fin, &nb_minima, minima);

  free_grphic3d(imtres);
  free_field3d(champ_fin);
 }

 imx_aff_param(transfres->typetrans, transfres->param);
 copie_transf_3d(transfres, transfoini);

 //------------------RECALAGE EN RESOLUTION 4 MM----------------------------------------/

 pyramide_gaussienne_recalage_3d_p(imref_isotrope, imreca_preprocessed, imtref, imtreca, 2.0);

 aff_log("\n minimisation sur image 64x64x64 \n");
 {
  champ_fin=cr_field3d(imtref->width,imtref->height,imtref->depth);
  imtres=cr_grphic3d(imtref);
  imx_copie_dim_param_to_transf3d_p(imtref, transfoini);

  imx_matching_multistart_multires_4mm_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, LARGE, matchPrecision, NULL, NULL, champ_fin, nb_minima, minima);

  free_grphic3d(imtres);
  free_field3d(champ_fin);
 }

 imx_aff_param(transfres->typetrans, transfres->param);
 copie_transf_3d(transfres, transfoini);

 //------------------RECALAGE EN RESOLUTION 2 MM----------------------------------------/

 pyramide_gaussienne_recalage_3d_p(imref_isotrope, imreca_preprocessed, imtref, imtreca, 1.0);

 aff_log("\n minimisation sur image 128x128x128 \n");
 {
  champ_fin=cr_field3d(imtref->width,imtref->height,imtref->depth);
  imtres=cr_grphic3d(imtref);
  imx_copie_dim_param_to_transf3d_p(imtref, transfoini);

  //recalage rigide+zoom global
  imx_matching_simple_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, SMALL, matchPrecision, NULL, NULL, champ_fin);

  //construction de la transfo initiale a injecter en debut du recalage rigide+zoom anisotrope
  imx_aff_param(transfres->typetrans, transfres->param);
  copie_transf_3d(transfres, transfoini);
  transfoini->param=REALLOC(transfoini->param,transfoini->nb_param+2,double);
  transfoini->nb_param=rigid_global_zoom_to_rigidz_3d(transfoini->param);
  transfoini->typetrans=RIGIDZOOM3D;

  //recalage rigide+zoom anisotrope
  imx_matching_simple_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, SMALL, matchPrecision, NULL, NULL, champ_fin);

  free_grphic3d(imtres);
  free_field3d(champ_fin);
 }

 imx_aff_param(transfres->typetrans, transfres->param);
 copie_transf_3d(transfres, transfoini);

 free_grphic3d(imtref); free_grphic3d(imtreca);

 //------------------RECALAGE EN RESOLUTION INITIALE----------------------------------------/

 aff_log("\n minimisation sur image originale \n");
 {
  imtref=imref_preprocessed; imtreca=imreca_preprocessed; imtres=cr_grphic3d(imtref);
  champ_fin=cr_field3d(imtref->width,imtref->height,imtref->depth);
  imx_copie_dim_param_to_transf3d_p(imtref, transfoini);

  imx_matching_simple_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, SMALL, matchPrecision, NULL, NULL, champ_fin);

  free_field3d(champ_fin);
  free_grphic3d(imtres);
 }

 //on remet le tout dans le repere des images originales
 imx_copie_dim_param_to_transf3d_p(imref, transfres);
 if (preProcessingRecaTrf) comb_transf_3d(preProcessingRecaTrf,transfres,transfres,0);
 if (preProcessingRefTrf) comb_transf_3d(transfres,preProcessingRefTrf,transfres,0);

 //liberations memoire
 if (preProcessingRecaTrf) free_transf3d(preProcessingRecaTrf);
 if (preProcessingRefTrf) free_transf3d(preProcessingRefTrf);
 if (imref_preprocessed!=imref) free_grphic3d(imref_preprocessed);
 if (imreca_preprocessed!=imreca) free_grphic3d(imreca_preprocessed);
  if (transfoini) free_transf3d(transfoini);
 for (k=0;k<nb_minima_max;k++) FREE(minima[k]);
 FREE(minima);
 free_grphic3d(imref_isotrope);
 
 // ------------------CALCUL DE L'IMAGE A AFFICHER ET DE LA TRANSFO A SAUVEGARDER----------------------
 imx_aff_param(transfres->typetrans, transfres->param);
 imx_calc_transfo_sav_image_res_3d_p(imref, imreca, imres, transfres, save_type); 

 tfin=time(NULL);
 //affichage du temps de calcul total
 {
  imx_aff_temps_calcul(tdebut, tfin);
 }

 return 0;
}

/*******************************************************************************
********************************************************************************
****************** MATCHING  AFFINE MULTISTART MULTIRES ************************
********************************************************************************

*******************************************************************************/

/*******************************************************************************
**        matching_affine_multistart_multires_3d
*/
/*!
**      \brief recalage affine multistart et multiresolution de deux images 3D
*******************************************************************************/
extern void     matching_affine_multistart_multires_3d()
{
  query_simple_matching_parameters_3d(imx_matching_affine_multistart_multires_3d);
}

/*******************************************************************************
**        matching_affine_multistart_multires_3d
*/
/*!
**      \brief recalage affine multistart et multiresolution de deux images 3D
*******************************************************************************/
void matching_affine_multistart_multires_precision_3d()
{
 query_simple_precision_matching_parameters_3d(imx_matching_affine_multistart_multires_3d);
}

/*******************************************************************************
**        imx_matching_affine_multistart_multires_3d
*/
/*!       recalage affine multistart et multiresolution  de 2 images 3D
**
**  \param im_ref : numero de l'image de reference,
**  \param im_reca  : numero de l'image a recaler,
**  \param im_res : numero de l'image resultat
**  les autres parametres sont expliques dans ::imx_matching_affine_3d_p
**  \retval 1
**
*******************************************************************************/
extern int  imx_matching_affine_multistart_multires_3d(int im_ref, int im_reca, int im_res, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres, char *nomfichtrf,eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision)
{
  grphic3d *imref,*imreca,*imres;
  transf3d *inittrf=NULL;
  transf3d *transfres=NULL;

  imref=ptr_img_3d(im_ref);
  imreca=ptr_img_3d(im_reca);
  imres=ptr_img_3d(im_res);

  if (nomfichtrf !=NULL)
  inittrf=load_transf_3d(nomfichtrf);
  transfres = cr_transf3d_p(imref,NULL);


  imx_inimaxminpixel_3d_p(imref);
  imx_inimaxminpixel_3d_p(imreca);
  if ((imref->min_pixel<0)||(imreca->min_pixel<0))
  { PUT_WARN("L'image contient des valeurs negatives!!\nOn ajoute -val_min pour la recaler\n"); }

  init_log(NULL);
  aff_log("\n RECALAGE AFFINE EN MULTIRESOLUTION AVEC MULTISTART (methode Jenkison&Smith) ");
  aff_log("\n");
  aff_log("recalage de %s sur %s \n",(imreca->patient).name,(imref->patient).name);

  imx_matching_affine_multistart_multires_3d_p(imref,imreca,imres,dist_type,inter_type,min_type,save_type,transfres,inittrf,FieldOfResearch,matchPrecision);

  /*enregistrement du champ resultat*/
  if (nomfichres!=NULL)
  {
    save_transf_3d(transfres,nomfichres);
  }

  end_log();

  if (inittrf)
    free_transf3d(inittrf);

  free_transf3d(transfres);

  return(0);
}

/*******************************************************************************
**        imx_matching_affine_multistart_multires_3d_p
*/
/*!       recalage affine multistart et multiresolution  de 2 images 3D
**
**  \param im_ref : numero de l'image de reference,
**  \param im_reca  : numero de l'image a recaler,
**  \param im_res : numero de l'image resultat
**  les autres parametres sont expliques dans ::imx_matching_rigide_3d_p
**  \retval 1
**
*******************************************************************************/
extern int imx_matching_affine_multistart_multires_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int dist_type, int inter_type, int min_type, int save_type, transf3d *transfres, transf3d *inittrf,eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision)
{
 int k;
 grphic3d *imtref,*imtreca,*imtres, *imref_preprocessed, *imreca_preprocessed, *imref_isotrope;
 field3d *champ_fin;
 InterpolationFct interpol;
 dist_func_t  distance;
 min_func_t minimisation;
 int    tdebut,tfin;
 int nb_minima=3, nb_minima_max=3;
 double **minima;
 transf3d *preProcessingRefTrf=NULL, *preProcessingRecaTrf=NULL, *transfoini=NULL;

 minima=CALLOC(nb_minima_max, double*);
 for (k=0;k<nb_minima_max;k++) minima[k]=CALLOC(15,double);

#ifndef WIN32
 setvbuf(stdout, (char *)NULL, _IONBF, 0); //setbuf(stdout,NULL);
#endif
 tdebut=time(NULL);

 //---------mise en place des parametres du recalage-------/

  distance=imx_choose_distance_fct(dist_type);

  interpol=imx_choose_interpolation_fct(inter_type);

  minimisation=imx_choose_minimisation_fct(min_type,dist_type);

 // ------------------PREPROCESS----------------------
 {
  imx_matching_preprocessing_3d_p(imref, imreca, &imref_preprocessed, &imreca_preprocessed, NULL, &preProcessingRefTrf, &preProcessingRecaTrf, distance);
 }

 // ------------------INITIALISATION PAR TRANSFO DE DEPART----------------------
 {
  //on veut une transformation de depart compatible avec le recalage courant
  if ((inittrf)&&(inittrf->typetrans==RIGID3D))
  {
   transfoini=imx_trouver_transf_initiale_matching_3d_p(imref_preprocessed, imreca_preprocessed, inittrf);
  }
  else transfoini=imx_trouver_transf_initiale_matching_3d_p(imref_preprocessed, imreca_preprocessed, NULL);
 }
 //on se met en rigide+zoom global directement
 transfoini->param=REALLOC(transfoini->param,transfoini->nb_param+1,double);
 transfoini->nb_param=rigid_to_rigid_global_zoom_3d(transfoini->param);
 transfoini->typetrans=RIGIDGLOBALZOOM3D;

 //calcule une image isotrope en resolution 2 mm a partir de l'image de depart
 imref_isotrope=cr_grphic3d(imref_preprocessed);
 imx_creer_image_isotrope_3d_p(imref_preprocessed, inter_qsinc3_3d, imref_isotrope);

 imtref=cr_grphic3d(imref_isotrope);
 imtreca=cr_grphic3d(imreca_preprocessed);

 //------------------RECALAGE EN RESOLUTION 8 MM----------------------------------------/

 pyramide_gaussienne_recalage_3d_p(imref_isotrope, imreca_preprocessed, imtref, imtreca, 4.0);

 aff_log("\n minimisation sur image 32x32x32 \n");
 {
  champ_fin=cr_field3d(imtref->width,imtref->height,imtref->depth);
  imtres=cr_grphic3d(imtref);
  imx_copie_dim_param_to_transf3d_p(imtref, transfoini);

  imx_matching_multistart_multires_8mm_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, LARGE, matchPrecision, NULL, NULL, champ_fin, &nb_minima, minima);

  free_grphic3d(imtres);
  free_field3d(champ_fin);
 }

 imx_aff_param(transfres->typetrans, transfres->param);
 copie_transf_3d(transfres, transfoini);

 //------------------RECALAGE EN RESOLUTION 4 MM----------------------------------------/

 pyramide_gaussienne_recalage_3d_p(imref_isotrope, imreca_preprocessed, imtref, imtreca, 2.0);

 aff_log("\n minimisation sur image 64x64x64 \n");
 {
  champ_fin=cr_field3d(imtref->width,imtref->height,imtref->depth);
  imtres=cr_grphic3d(imtref);
  imx_copie_dim_param_to_transf3d_p(imtref, transfoini);

  imx_matching_multistart_multires_4mm_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, LARGE, matchPrecision, NULL, NULL, champ_fin, nb_minima, minima);

  free_grphic3d(imtres);
  free_field3d(champ_fin);
 }

 imx_aff_param(transfres->typetrans, transfres->param);
 copie_transf_3d(transfres, transfoini);

 //------------------RECALAGE EN RESOLUTION 2 MM----------------------------------------/

 pyramide_gaussienne_recalage_3d_p(imref_isotrope, imreca_preprocessed, imtref, imtreca, 1.0);

 aff_log("\n minimisation sur image 128x128x128 \n");
 {
  champ_fin=cr_field3d(imtref->width,imtref->height,imtref->depth);
  imtres=cr_grphic3d(imtref);
  imx_copie_dim_param_to_transf3d_p(imtref, transfoini);

  //recalage rigide+zoom global
imx_matching_simple_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, SMALL, matchPrecision, NULL, NULL, champ_fin);

  //construction de la transfo initiale a injecter en debut du recalage rigide+zoom anisotrope
  imx_aff_param(transfres->typetrans, transfres->param);
  copie_transf_3d(transfres, transfoini);
  transfoini->param=REALLOC(transfoini->param,transfoini->nb_param+2,double);
  transfoini->nb_param=rigid_global_zoom_to_rigidz_3d(transfoini->param);
  transfoini->typetrans=RIGIDZOOM3D;

  //recalage rigide+zoom anisotrope
  imx_matching_simple_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, SMALL, matchPrecision, NULL, NULL, champ_fin);

  //construction de la transfo initiale a injecter en debut du recalage affine decouple
  imx_aff_param(transfres->typetrans, transfres->param);
  copie_transf_3d(transfres, transfoini);
  transfoini->param=REALLOC(transfoini->param,transfoini->nb_param+3,double);
  transfoini->nb_param=rigidz_to_affine_decouple_3d(transfoini->param);
  transfoini->typetrans=AFFINEDECOUPLE3D;

  //recalage affine decouple
  imx_matching_simple_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, SMALL, matchPrecision, NULL, NULL, champ_fin);

  free_grphic3d(imtres);
  free_field3d(champ_fin);
 }

 imx_aff_param(transfres->typetrans, transfres->param);
 copie_transf_3d(transfres, transfoini);

 free_grphic3d(imtref); free_grphic3d(imtreca);

 //------------------RECALAGE EN RESOLUTION INITIALE----------------------------------------/

 aff_log("\n minimisation sur image originale \n");
 {
  imtref=imref_preprocessed; imtreca=imreca_preprocessed; imtres=cr_grphic3d(imtref);
  champ_fin=cr_field3d(imtref->width,imtref->height,imtref->depth);
  imx_copie_dim_param_to_transf3d_p(imtref, transfoini);

  imx_matching_simple_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, SMALL, matchPrecision, NULL, NULL, champ_fin);

  free_field3d(champ_fin);
  free_grphic3d(imtres);
 }

 //on remet le tout dans le repere des images originales
 transfres->nb_param=affine_decouple_to_affine_3d(transfres->param);
 transfres->typetrans=AFFINE3D;
 imx_copie_dim_param_to_transf3d_p(imref, transfres);
 if (preProcessingRecaTrf) comb_transf_3d(preProcessingRecaTrf,transfres,transfres,0);
 if (preProcessingRefTrf) comb_transf_3d(transfres,preProcessingRefTrf,transfres,0);

 //liberations memoire
 if (preProcessingRecaTrf) free_transf3d(preProcessingRecaTrf);
 if (preProcessingRefTrf) free_transf3d(preProcessingRefTrf);
 if (imref_preprocessed!=imref) free_grphic3d(imref_preprocessed);
 if (imreca_preprocessed!=imreca) free_grphic3d(imreca_preprocessed);
 for (k=0;k<nb_minima_max;k++) FREE(minima[k]);
 FREE(minima);
 free_grphic3d(imref_isotrope);
 if (transfoini) free_transf3d(transfoini);

 // ------------------CALCUL DE L'IMAGE A AFFICHER ET DE LA TRANSFO A SAUVEGARDER----------------------
 imx_aff_param(transfres->typetrans, transfres->param);
 imx_calc_transfo_sav_image_res_3d_p(imref, imreca, imres, transfres, save_type);
 
 

 tfin=time(NULL);
 //affichage du temps de calcul total
 {
  imx_aff_temps_calcul(tdebut, tfin);
 }

 return 0;
}

/*******************************************************************************
********************************************************************************
************************** MATCHING  BSPLINE ************************************
********************************************************************************
*******************************************************************************/

/*******************************************************************************
**        matching_Bspline_3d                                                
*/                                                                        
/*!       recalage Bspline de deux images 2D                                
*******************************************************************************/
void matching_Bspline_3d()
{
 int im_ref,im_reca,im_res;
 char *quest[18],*nomfichres;
 int dist_type,inter_type,min_type,save_type,func_type,i;
 int renormalisation;

 int resolf, e=0;
 char *nomfichiertrf;
 
 /*question sur les images a recaler*/
 im_reca=GET_PLACE3D(TEXT0030);
 im_ref=GET_PLACE3D(TEXT0031);
 im_res=GET_PLACE3D(TEXT0006);

 for(i=0;i<18;i++)
    quest[i]=CALLOC(80,char);
 
 /*question sur les parametres de la methode*/
 /*type de fonction*/
  strcpy(quest[0],"B-Spline d0");
  strcpy(quest[1],"B-Spline d1");
  strcpy(quest[2],"B-Spline d2");
  strcpy(quest[3],"B-Spline o d1");
  strcpy(quest[4],"\0");
  func_type=GETV_QCM("Fonction",(char **)quest);
 /*distance*/
  dist_type=imx_query_erreur_type_3D();
 /*interpolation*/
 inter_type=imx_query_interpolation_type_3D(dist_type);
/*minimisation*/
  min_type=imx_query_minimisation_type_3D(dist_type);
   
 for(i=0;i<7;i++)
    free(quest[i]);

 /*question sur l'enregistrement du champ resultat*/
 nomfichres=quest_save_result_3d(&save_type);

 {
  char  s[256];
  int   ans,err=0;
  

  sprintf(s,"Use a transformation as initialisation?\n");
  ans=GET_EXIST(s,&err);
  if (ans==1)
  {
   char nomfichier[250];
   
   strcpy(nomfichier,GET_FILE("*.trf",&err));
   put_file_extension(nomfichier,".trf",nomfichier);
   nomfichiertrf=CALLOC(strlen(nomfichier),char);
   strcpy(nomfichiertrf,nomfichier);
  }
  else
    nomfichiertrf = NULL;
 }

 /*on demande la resolution finale*/
 resolf= GET_INT("resolution", 0, &e);

 /* renormalisation : TRUE si l'on veut effectuer une renormalisation des
    images (recommande) */

 renormalisation=FALSE;

 imx_matching_Bspline_3d(im_ref,im_reca,im_res,func_type,dist_type,inter_type,min_type,save_type,nomfichres,resolf,nomfichiertrf,renormalisation);

 if (nomfichiertrf)
  free(nomfichiertrf);
 
 if (nomfichres)
  free(nomfichres);
  
 show_picture_3d(im_res);
 
}

/*******************************************************************************
**        imx_matching_Bspline_3d                                             
*/
/*!                                                                       
**       recalage Bspline de 2 images 3D                                   
*******************************************************************************/
/*!
  \param im_ref : numero de l'image sur laquelle on recale imreca
  \param im_reca : numero de l'image a recaler
  \param im_res : numero de l'image ou l'on stocke le resultat du recalage de imreca
    (peut etre NULL)
  \param func_type : numero de la fonction d'interpolation du champ (1 pour Bspline
    degre 1)
  \param dist_type : numero de la fonction de distance (0 pour distance quadratique)

  \param inter_type : numero de la fonction d'interpolation des images (1 pour lineaire)
  \param min_type : numero de la fonction de minimisation utilisee aux plus faibles
    echelles de la transformation (5 pour quasi-newton)
  \param save_type : 0 si on ne veut rien sauver, 1 si on veut sauver les parametres
    de la transformation (fichier trf), 2 si on veut sauver le champ
  \param nomfichres : fichier ou le resultat est sauve, si (save_type != 0)
  \param resolf : resolution finale du recalage ; typiquement 4 pour des images
    128*128*128
  \param  nomfichiertrf : nom du fichier trf a partir duquel on commence le recalage
    (peut etre NULL, mais ce n'est pas recommande ; il convient d'effectuer
    un recalage affine et de donner la transformation resultante comme etat
    initial)
  \param renormalisation : TRUE si l'on veut effectuer une renormalisation des
    images (recommande)
  \retval 1
*/
int imx_matching_Bspline_3d(int im_ref, int im_reca, int im_res, int func_type, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres, int resolf, char *nomfichiertrf, int renormalisation)
{
 grphic3d *imref,*imreca,*imres;
 
 imref=ptr_img_3d(im_ref);
 imreca=ptr_img_3d(im_reca);
 imres=ptr_img_3d(im_res);
 imx_matching_Bspline_3d_p(imref,imreca,imres,func_type,dist_type,inter_type,min_type,save_type,nomfichres,resolf,nomfichiertrf,renormalisation);
 
  return(1);
}

/*! \ingroup     Recalage  @{ */
/*******************************************************************************
**        imx_matching_Bspline_3d_p                                            
**                                                                        
**       recalage Bspline de 2 images 3D                                    
*******************************************************************************/
/*!  \brief recalage Bspline de 2 images 3D     
  \param imref : pointeur sur l'image sur laquelle on recale imreca
  \param imreca : pointeur sur l'image a recaler
  \param imres : pointeur sur l'image ou l'on stocke le resultat du recalage de imreca
    (peut etre NULL)
  \param func_type : degre de la fonction bspline (1 pour Bspline
    degre 1)
  \param dist_type : numero de la fonction de distance (0 pour distance quadratique)
  \param inter_type : numero de la fonction d'interpolation des images (1 pour lineaire)
  \param min_type : numero de la fonction de minimisation utilisee aux plus faibles
    echelles de la transformation (5 pour quasi-newton)
  \param save_type : 0 si on ne veut rien sauver, 1 si on veut sauver les parametres
    de la transformation (fichier trf), 2 si on veut sauver le champ
  \param nomfichres : fichier ou le resultat est sauve, si (save_type != 0)
  \param resolf : resolution finale du recalage ; typiquement 4 pour des images
    128*128*128
  \param nomfichiertrf : nom du fichier trf a partir duquel on commence le recalage
    (peut etre NULL, mais ce n'est pas recommande ; il convient d'effectuer
    un recalage affine et de donner la transformation resultante comme etat
    initial)
  \param renormalisation : TRUE si l'on veut effectuer une renormalisation des
    images (recommande)
  \retval 1
*/
int imx_matching_Bspline_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres,
                                int func_type, int dist_type, int inter_type,
                                int min_type, int save_type, char *nomfichres,
                                int resolf, char *nomfichiertrf, int renormalisation)
{
 TDimension wdth,hght,dpth;
 int nb_param;
 double *param,*min_param,*max_param,*prec_param;
 grphic3d *imtref,*imtreca,*imtres;
 field3d *champ,*champ_ini;
 InterpolationFct interpol;
 dist_func_t distance;
 min_func_t minimisation;
 scal_func_t scal_func;
 double  /*tx,ty,xg,yg,teta,sx,sy,*/mini;
 int    tdebut,tfin;
 int nbmax_param=760000; /* valable jusqu'a resolution 6 */
 ERREUR_RECALAGE err=NO_ERR;

 ptr_minimisation_param min_parameters=CALLOC(1, minimisation_param);
 
 wdth=imref->width;hght=imref->height;dpth=imref->depth;
#ifndef WIN32
 setvbuf(stdout, (char *)NULL, _IONBF, 0); // setbuf(stdout,NULL);
#endif
 tdebut=time(NULL);
 
   /* Allocation tableau param[],min_param[],max_param[],prec_param[];*/
 if(resolf==1)nbmax_param=3+10;
 if(resolf==2)nbmax_param=81+10;
 if(resolf==3)nbmax_param=1029+10;
 if(resolf==4)nbmax_param=10125+10;
 if(resolf==5)nbmax_param=89373+10;
 if(resolf==6)nbmax_param=750141+10;
// if((param = calloc(nbmax_param, sizeof(double))) == NULL) {
//   PUT_ERROR("[imx_matching_Bspline_3d_p] memory allocation error !\n"); return(0); }
 if((min_param = calloc(nbmax_param, sizeof(double))) == NULL) {
   PUT_ERROR("[imx_matching_Bspline_3d_p] memory allocation error !\n"); return(0); }
 if((max_param = calloc(nbmax_param, sizeof(double))) == NULL) {
   PUT_ERROR("[imx_matching_Bspline_3d_p] memory allocation error !\n"); return(0); }
 if((prec_param = calloc(nbmax_param, sizeof(double))) == NULL) {
   PUT_ERROR("[imx_matching_Bspline_3d_p] memory allocation error !\n"); return(0); }
   
 init_log(NULL);
 aff_log("\nRECALAGE ");
 /*************mise en place des parametres du recalage******/
  switch (func_type)
  {
    case 0: scal_func=Bsplined0;aff_log("BSPLINE DEGRE 0 ");break;
    case 1: scal_func=Bsplined1;aff_log("BSPLINE DEGRE 1 ");break;
    case 2: scal_func=Bsplined2;aff_log("BSPLINE DEGRE 2 ");break;
    case 3: scal_func=Bsplineod1;aff_log("BSPLINE O DEGRE 1 ");break;
    default: scal_func=Bsplined0;aff_log("BSPLINE DEGRE 0 ");break;
  }
  distance=imx_choose_distance_fct(dist_type);
  interpol=imx_choose_interpolation_fct(inter_type);  
  minimisation=imx_choose_minimisation_fct(min_type,dist_type);
 aff_log("\n");
 aff_log("recalage de %s sur %s \n",(imreca->patient).name,(imref->patient).name);
 
 //creation des images temporaires
 imtref=cr_grphic3d(imref);imx_copie_3d_p(imref,imtref);
 imtreca=cr_grphic3d(imreca);imx_copie_3d_p(imreca,imtreca);
 //initialisations specifiques a la distance utilisee 
 if ((distance==erreur_quad_3d)||(distance==erreur_quad_robust_geman_3d))
 {
  //normalisation de imreca par rapport a imref: resultat dans imtreca
  imx_norm_seuil_meanecty_3d_p(imreca,imtref,imtreca);
 }
 imtres=cr_grphic3d(imtref);

 /*initialisation du champ de depart a la valeur d'un fichier*/
 if (nomfichiertrf != NULL)
 {
//#if 0
//   /* arrgh, on fragmente le tas ! */
//   transf3d *transfo;
//
//   transfo = load_transf_3d(nomfichiertrf);
//   if ((transfo->width != wdth) || (transfo->height != hght)
//         || (transfo->depth != dpth))
//   {
//     /*err =*/ resize_transf_3d(transfo, wdth, hght, dpth);
//     /*if (err)
//       goto error1;*/
//   }
//   champ_ini = transf_to_field_3d(transfo);
//   free_transf3d(transfo);
//#else
   /* correction pour ce probleme */
   transf3d *transfo;

   champ_ini=cr_field3d(wdth,hght,dpth); 

   transfo = load_transf_3d(nomfichiertrf);
   if ((transfo->width != wdth) || (transfo->height != hght)
         || (transfo->depth != dpth))
   {
     /*err =*/ resize_transf_3d(transfo, wdth, hght, dpth);
     /*if (err)
       goto error1;*/
   }
   /*err =*/ transf_to_field_3d_noalloc(transfo, champ_ini,imref,imreca);
   /*if (err)
     goto error1;*/
   free_transf3d(transfo);
//#endif
 }
 else
   champ_ini = NULL;

  /*allocation memoire de champ*/
   champ=cr_field3d(wdth,hght,dpth); 
   
 /********************RECALAGE ONDELETTE******************************************/
 {
  int l,resol,r;
  double lar;
  lar=5.0;

  nb_param=init_base_3d(wdth,hght,dpth,scal_func);
  resol=BASE3D.resol;

//-----------------------------
 min_parameters->imref=imtref;
 min_parameters->imreca=imtreca;
 min_parameters->imres=imtres;
 min_parameters->champ_ini=champ_ini;
 min_parameters->champ_fin=champ;
 min_parameters->inter=interpol;
 min_parameters->dist=distance;
 min_parameters->contrainte=NULL;
 min_parameters->min_param=min_param;
 min_parameters->max_param=max_param;
 min_parameters->prec_param=prec_param;
 min_parameters->transfres=cr_transf3d_p(imtref, NULL);
 min_parameters->transfres->typetrans=BSPLINE3D;
 min_parameters->transfres->param=CALLOC(nbmax_param, double);
 min_parameters->err=&err;
 min_parameters->dist_param=cr_distParam_from_minParam(min_parameters, minimisation);

 param=min_parameters->transfres->param;
  
  /* Initialisation tableau */
  for (l=0;l<nb_param;l++) 
   {param[l]=0.0;prec_param[l]=0.001;}
  /* Debut traitement */
  for (r=resol;r<=resolf;r++)
  {
   aff_log("RESOLUTION %d ",r);
   for (l=0;l<nb_param;l++) {min_param[l]=param[l]-lar;max_param[l]=param[l]+lar;prec_param[l]=0.001;}

   min_parameters->transfres->nb_param=nb_param;

   if (r<RESOL_LOCAL) mini=minimisation(min_parameters);
   else mini=min_flux_local_3d(min_parameters);
   
   if (r!=resolf)
    {/*on passe a la resolution superieure*/
     nb_param=base_resol_up_3d(param,nb_param);printf("nb_param: %d ",nb_param);
     nb_param=reduc_base_3d(param,nb_param,imtref,imtreca,imtres);
     printf("   nb_param: %d \n",nb_param);
    }

    /*affichage du temps de calcul partiel*/
    tfin=time(NULL);
    {
     imx_aff_temps_calcul(tdebut, tfin);
    }

  }
   
  if (imres != NULL)
  { 
    interpol=inter_qsinc3_3d;
    interpol(imreca,champ,imtres);
    imx_copie_3d_p(imtres,imres);
  }

  /*enregistrement du champ resultat dans un fichier*/
   if (nomfichres!=NULL)
   { 

    transf3d *transfo; 
    
    if (save_type==0) transfo=field_to_transf_3d(champ,imtref,imtreca); 
    else
    {  
     double   ***Mx,***My,***Mz;
     int    i,j,k;  
       
     /*on met dans param tous les parametres, y compris ceux qui sont nul=> il faut annuler
     la troncature du a reduc_base_3d*/
     Mx=alloc_dmatrix_3d(BASE3D.nbfx,BASE3D.nbfy,BASE3D.nbfz);
     My=alloc_dmatrix_3d(BASE3D.nbfx,BASE3D.nbfy,BASE3D.nbfz);
     Mz=alloc_dmatrix_3d(BASE3D.nbfx,BASE3D.nbfy,BASE3D.nbfz);
     for (i=0;i<BASE3D.nbfx;i++) for (j=0;j<BASE3D.nbfy;j++) for (k=0;k<BASE3D.nbfz;k++)
        Mx[i][j][k]=My[i][j][k]=Mz[i][j][k]=0.0; 
     for (l=0;l<BASE3D.nb_func;l++) 
     {i=BASE3D.idx[l];j=BASE3D.idy[l];k=BASE3D.idz[l];
     Mx[i][j][k]=param[3*l];My[i][j][k]=param[3*l+1];Mz[i][j][k]=param[3*l+2];}
     nb_param=3*BASE3D.nbfx*BASE3D.nbfy*BASE3D.nbfz;
     l=0;
     for (i=0;i<BASE3D.nbfx;i++) for (j=0;j<BASE3D.nbfy;j++) for (k=0;k<BASE3D.nbfz;k++)
      {param[3*l]=Mx[i][j][k];param[3*l+1]=My[i][j][k];param[3*l+2]=Mz[i][j][k];l++;}
     free_dmatrix_3d(Mx);free_dmatrix_3d(My);free_dmatrix_3d(Mz);
     
     
     transfo=cr_transf3d_p(imref,nomfichiertrf);
     transfo->nb_param=nb_param;transfo->param=CALLOC(nb_param,double);
     transfo->typetrans=BSPLINE3D;
     transfo->resol=resolf;transfo->degre=func_type;
         transfo->dx=imres->dx;
         transfo->dy=imres->dy;
         transfo->dz=imres->dz;
         
     for (l=0;l<nb_param;l++) 
            {
            transfo->param[l]=param[l];
            } 
    }
    save_transf_3d(transfo,nomfichres);
    free_transf3d(transfo);
   } 
 }

  /*liberation de la base*/
   end_base_3d();

  /*caracteristique du champ*/
  if (champ_ini==NULL)
   {
    aff_log("CHAMP: ");carac_field_3d(champ);
   }
  else
   {
    aff_log("CHAMP TOTAL: ");carac_field_3d(champ);
    sub_field_3d(champ,champ_ini,champ);
    aff_log("CHAMP DEFOR: ");carac_field_3d(champ);
   }
 /*liberation memoire des variables dynamiques*/ 
  /*liberation memoire du champ*/
   free_field3d(champ);
   if (champ_ini!=NULL) free_field3d(champ_ini);
 
 //liberation memoire des images temporaires
 free_grphic3d(imtref); free_grphic3d(imtreca); free_grphic3d(imtres);
 //liberations memoire specifiques a la distance utilisee

 /*affichage du temps de calcul total*/
 tfin=time(NULL);
 {
  imx_aff_temps_calcul(tdebut, tfin);
 }
 end_log();
 
// free(param);
 free(min_param);free(max_param);free(prec_param);

 FREE(min_parameters);
 
 return(1);  
}
/*! @}  */

#if 0
/*
  imx_matching_Bspline_pyramide_3d_p() : semblable a imx_matching_Bspline_3d_p,
  mais on effectue le recalage en variant la resolution des images

pyr_imref : pyramide gaussienne
pyr_imreca : pyramide gaussienne
nb_niv_pyr : nombre de niveaux dans pyr_imref et pyr_imreca (on prendra le plus
petit des deux s'ils sont differents)
*/
static void change_niveau(int nouveau_niveau, int ancien_niveau,
                            grphic3d **pyr_imref, grphic3d **pyr_imreca,
                            grphic3d *imtref, grphic3d *imtreca, grphic3d *imtres,
                            int *wdth, int *hght, int *dpth,
                            transf3d *transfo_ini,
                            field3d **champ, field3d **champ_ini,
                            double *param, int *nb_param, int renormalisation)
{
  int i;
  double coeff;


  /*liberation memoire du champ*/
  free_field3d(*champ);
  if (*champ_ini != NULL)
    free_field3d(*champ_ini);

  /* mise a l'echelle des parametres de la tranformation */
  coeff = pow(2., ancien_niveau-nouveau_niveau);
  for (i=0; i<*nb_param; i++)
  {
    param[i] *= coeff;
  }

  /* et on refait tout le boulot */
  *wdth = pyr_imref[nouveau_niveau]->width;
  *hght = pyr_imref[nouveau_niveau]->height;
  *dpth = pyr_imref[nouveau_niveau]->depth;
  /*normalisation de imreca par rapport a imref: resultat dans imtreca*/ 
  if (renormaliation)
  {
    imx_copie_3d_p(pyr_imref[nouveau_niveau], imtref);
    imx_copie_3d_p(pyr_imreca[nouveau_niveau], imtreca);

    imx_norm_seuil_meanecty_3d_p(pyr_imreca[nouveau_niveau], imtref, imtreca);
  }
  else
  {
    imtref = pyr_imref[nouveau_niveau];
    imtreca = pyr_imreca[nouveau_niveau];
  }
  imx_copie_param_3d_p(pyr_imreca[nouveau_niveau], imtres);

  /* reallocation de la base a notre nouvelle resolution de travail */
  (void) resize_base_3d(*wdth, *hght, *dpth);
  *nb_param = reduc_base_3d(param, *nb_param, imtref, imtreca, imtres);


  /* reinitialisation du champ de depart */
  if (transfo_ini)
  {
    if ((transfo_ini->width != *wdth) || (transfo_ini->height != *hght)
          || (transfo_ini->depth != *dpth))
    {
      /*err =*/ resize_transf_3d(transfo_ini, *wdth, *hght, *dpth);
      /*if (err)
        goto error1;*/
    }
    *champ_ini = transf_to_field_3d(transfo_ini);
  }
  else
    *champ_ini = NULL;

  /* reallocation de champ */
  *champ = cr_field3d(*wdth, *hght, *dpth);
}


int imx_matching_Bspline_pyramide_3d_p(grphic3d **pyr_imref, grphic3d **pyr_imreca,
                                        int nb_niv_pyr,
                                        grphic3d *imres, int func_type,
                                        int dist_type, int inter_type,
                                        int min_type, int save_type,
                                        char *nomfichres, int resolf,
                                        char *nomfichiertrf, int renormalisation)
{
 int    /*k,*/wdth,hght,dpth,nb_param;
 double param[50000],min_param[50000],max_param[50000],prec_param[50000];
 grphic3d *imtref,*imtreca,*imtres;
 field3d *champ,*champ_ini;
 InterpolationFct interpol;
 dist_func_t distance;
 min_func_t minimisation;
 scal_func_t scal_func;
 double  /*tx,ty,xg,yg,teta,sx,sy,*/mini;
 int    tdebut,tfin;
 transf3d *transfo_ini;

 
#ifndef WIN32
 setvbuf(stdout, (char *)NULL, _IONBF, 0); //setbuf(stdout,NULL);
#endif
 tdebut=time(NULL);
 
 init_log(NULL);
 aff_log("\nRECALAGE ");

 /*************mise en place des parametres du recalage******/
  switch (func_type)
  {
    case 0: scal_func=Bsplined0;aff_log("BSPLINE DEGRE 0 ");break;
    case 1: scal_func=Bsplined1;aff_log("BSPLINE DEGRE 1 ");break;
    case 2: scal_func=Bsplined2;aff_log("BSPLINE DEGRE 2 ");break;
    case 3: scal_func=Bsplineod1;aff_log("BSPLINE O DEGRE 1 ");break;
    default: scal_func=Bsplined0;aff_log("BSPLINE DEGRE 0 ");break;
  }
  distance=imx_choose_distance_fct(dist_type);
  interpol=imx_choose_interpolation_fct(inter_type);  
  minimisation=imx_choose_minimisation_fct(min_type,dist_type);

 aff_log("\n");
 aff_log("recalage de %s sur %s \n",(pyr_imreca[0]->patient).name,(pyr_imref[0]->patient).name);
 
  /*creation, allocation et initialisation des variables dynamiques*/
  /*creation des images temporaires*/
  if (renormalisation)
  {
    imtref=cr_grphic3d(pyr_imref[0]);
    imtreca=cr_grphic3d(pyr_imreca[0]);
  }
  imtres=cr_grphic3d(pyr_imreca[0]);

  /*initialisation du champ de depart a la valeur d'un fichier*/
  if (nomfichiertrf != NULL)
  {

    transfo_ini = load_transf_3d(nomfichiertrf);
  }
  else
    transfo_ini = NULL;
   

  /********************RECALAGE ONDELETTE******************************************/
  {
    int l,resol,r;
    double lar;

    int niv_pyr_voulu, niv_pyr_reel;

    lar=5.0;

    switch (func_type)
    {
    case 0: /* BSPLINE DEGRE 0 */
      resol = 0;
      break;
    case 1: /* BSPLINE DEGRE 1 */
      resol = 1;
      break;
    case 2: /* BSPLINE DEGRE 2 */
      resol = 2;
      break;
    case 3: /* BSPLINE O DEGRE 1 */
      resol = 1;
      break;
    default: /* BSPLINE DEGRE 0 */
      resol = 0;
      break;
    }


#if 0
    niv_pyr_voulu = resolf - resol;
#else
    niv_pyr_voulu = resolf - resol - 1;
    if (niv_pyr_voulu < 0)
      niv_pyr_voulu = 0;
#endif
    niv_pyr_reel = (niv_pyr_voulu <= nb_niv_pyr) ? niv_pyr_voulu : nb_niv_pyr;

    wdth=pyr_imref[niv_pyr_reel]->width;
    hght=pyr_imref[niv_pyr_reel]->height;
    dpth=pyr_imref[niv_pyr_reel]->depth;

    /*normalisation de imreca par rapport a imref: resultat dans imtreca*/ 
    if (renormalisation)
    {
      imx_copie_3d_p(pyr_imref[niv_pyr_reel],imtref);
      imx_copie_3d_p(pyr_imreca[niv_pyr_reel],imtreca);
      imx_norm_seuil_meanecty_3d_p(pyr_imreca[niv_pyr_reel],imtref,imtreca);
    }
    else
    {
      imtref = pyr_imref[niv_pyr_reel];
      imtreca = pyr_imreca[niv_pyr_reel];
    }
    imx_copie_param_3d_p(pyr_imreca[niv_pyr_reel],imtres);

    /* initialisation du champ de depart */
    if (transfo_ini)
    {
      if ((transfo_ini->width != wdth) || (transfo_ini->height != hght)
           || (transfo_ini->depth != dpth))
      {
        /*err =*/ resize_transf_3d(transfo_ini, wdth, hght, dpth);
        /*if (err)
          goto error1;*/
      }
      champ_ini = transf_to_field_3d(transfo_ini);
    }
    else
      champ_ini = NULL;

    /*allocation memoire de champ*/
    champ=cr_field3d(wdth,hght,dpth); 

    /* allocation de la base */
    nb_param=init_base_3d(wdth,hght,dpth,scal_func);
    /* parametres initiaux a zero */
    for (l=0;l<nb_param;l++) 
    {param[l]=0.0;prec_param[l]=0.001;}



    for (r=resol;r<=resolf;r++)
    {
      aff_log("RESOLUTION %d ",r);
      for (l=0;l<nb_param;l++) {min_param[l]=param[l]-lar;max_param[l]=param[l]+lar;prec_param[l]=0.001;}
      if (r<RESOL_LOCAL)
        mini=minimisation(imtref,imtreca,imtres,champ_ini,champ,base_to_field_3d,interpol,distance,NULL,nb_param,min_param,max_param,prec_param,param);
      else
        mini=min_flux_local_3d(imtref,imtreca,imtres,champ_ini,champ,base_to_field_3d,interpol,distance,NULL,nb_param,min_param,max_param,prec_param,param);

      if (r!=resolf)
      {

        /* on passe a la resolution superieure pour les images */
        if (niv_pyr_voulu > 0)
        {
          niv_pyr_voulu--;
          if (niv_pyr_voulu < nb_niv_pyr)
          { /* actual resolution has changed */
#if 0
            niv_pyr_reel = niv_pyr_voulu;

            /*liberation memoire du champ*/
            free_field3d(champ);
            if (champ_ini!=NULL) free_field3d(champ_ini);

            /* mise a l'echelle des parametres de la tranformation */
            for (l=0;l<nb_param;l++)
            {
              param[l] *= 2.;
            }

            /* et on refait tout le boulot */
            wdth=pyr_imref[niv_pyr_reel]->width;
            hght=pyr_imref[niv_pyr_reel]->height;
            dpth=pyr_imref[niv_pyr_reel]->depth;
            /*normalisation de imreca par rapport a imref: resultat dans imtreca*/ 
            if (renormalisation)
            {
              imx_copie_3d_p(pyr_imref[niv_pyr_reel],imtref);
              imx_copie_3d_p(pyr_imreca[niv_pyr_reel],imtreca);
              imx_norm_seuil_meanecty_3d_p(pyr_imreca[niv_pyr_reel],imtref,imtreca);
            }
            else
            {
              imtref = pyr_imref[niv_pyr_reel];
              imtreca = pyr_imreca[niv_pyr_reel];
            }
            imx_copie_param_3d_p(pyr_imreca[niv_pyr_reel],imtres);

            /* reallocation de la base a notre nouvelle resolution de travail */
            (void)resize_base_3d(wdth, hght, dpth);
            nb_param=reduc_base_3d(param,nb_param,imtref,imtreca,imtres);

            /* reinitialisation du champ de depart */
            if (transfo_ini)
            {
              if ((transfo_ini->width != wdth) || (transfo_ini->height != hght)
                    || (transfo_ini->depth != dpth))
              {
                /*err =*/ resize_transf_3d(transfo_ini, wdth, hght, dpth);
                /*if (err)
                  goto error1;*/
              }
              champ_ini = transf_to_field_3d(transfo_ini);
            }
            else
              champ_ini = NULL;

            /* reallocation de champ */
            champ=cr_field3d(wdth,hght,dpth);

#elif 0

            change_niveau(niv_pyr_voulu, niv_pyr_reel, pyr_imref, pyr_imreca,
                            imtref, imtreca, imtres, &wdth, &hght, &dpth,
                            transfo_ini, &champ, &champ_ini, param, &nb_param);

            niv_pyr_reel = niv_pyr_voulu;

#if 0
            /* nouvelle minimisation */
            for (l=0;l<nb_param;l++) {min_param[l]=param[l]-lar;max_param[l]=param[l]+lar;prec_param[l]=0.001;}
            if (r<RESOL_LOCAL)
              mini=minimisation(imtref,imtreca,imtres,champ_ini,champ,base_to_field_3d,interpol,distance,NULL,nb_param,min_param,max_param,prec_param,param);
            else
              mini=min_flux_local_3d(imtref,imtreca,imtres,champ_ini,champ,base_to_field_3d,interpol,distance,NULL,nb_param,min_param,max_param,prec_param,param);
#endif

#else

            change_niveau(0, niv_pyr_reel, pyr_imref, pyr_imreca,
                            imtref, imtreca, imtres, &wdth, &hght, &dpth,
                            transfo_ini, &champ, &champ_ini, param, &nb_param);

            /* nouvelle minimisation */
            for (l=0;l<nb_param;l++) {min_param[l]=param[l]-lar;max_param[l]=param[l]+lar;prec_param[l]=0.001;}
            if (r<RESOL_LOCAL)
              mini=minimisation(imtref,imtreca,imtres,champ_ini,champ,base_to_field_3d,interpol,distance,NULL,nb_param,min_param,max_param,prec_param,param);
            else
              mini=min_flux_local_3d(imtref,imtreca,imtres,champ_ini,champ,base_to_field_3d,interpol,distance,NULL,nb_param,min_param,max_param,prec_param,param);

            change_niveau(niv_pyr_voulu, 0, pyr_imref, pyr_imreca,
                            imtref, imtreca, imtres, &wdth, &hght, &dpth,
                            transfo_ini, &champ, &champ_ini, param, &nb_param);


            niv_pyr_reel = niv_pyr_voulu;

#endif
          }
        }



      /*on passe a la resolution superieure pour la transformation*/
      nb_param=base_resol_up_3d(param,nb_param);printf("nb_param: %d ",nb_param);
      nb_param=reduc_base_3d(param,nb_param,imtref,imtreca,imtres);
      printf("   nb_param: %d \n",nb_param);
    }


    /*affichage du temps de calcul partiel*/
    tfin=time(NULL);
    {
      imx_aff_temps_calcul(tdebut, tfin);
    }
  }

  /*liberation de la base*/
  end_base_3d();
 
  free_transf3d(transfo_ini);


  if (imres != NULL)
  { /* des fois qu'on appelerait la fonction juste pour le champ, pas le resultat */
    /* application de la transformation aux images de depart
    Dans tous les cas interpolation qsinc3 */
   interpol=inter_qsinc3_3d;
   interpol(pyr_imreca[0],champ,imtres);
   imx_copie_3d_p(imtres,imres);
  }

  /*enregistrement du champ resultat dans un fichier*/
   if (nomfichres!=NULL)
   { 
    transf3d *transfo; 
    
    if (save_type==0) transfo=field_to_transf_3d(champ); 
    else
    {  
     double   ***Mx,***My,***Mz;
     int    i,j,k;  
       
     /*on met dans param tous les parametres, y compris ceux qui sont nul=> il faut annuler
     la troncature du a reduc_base_3d*/
     Mx=alloc_dmatrix_3d(BASE3D.nbfx,BASE3D.nbfy,BASE3D.nbfz);
     My=alloc_dmatrix_3d(BASE3D.nbfx,BASE3D.nbfy,BASE3D.nbfz);
     Mz=alloc_dmatrix_3d(BASE3D.nbfx,BASE3D.nbfy,BASE3D.nbfz);
     for (i=0;i<BASE3D.nbfx;i++) for (j=0;j<BASE3D.nbfy;j++) for (k=0;k<BASE3D.nbfz;k++)
        Mx[i][j][k]=My[i][j][k]=Mz[i][j][k]=0.0; 
     for (l=0;l<BASE3D.nb_func;l++) 
     {i=BASE3D.idx[l];j=BASE3D.idy[l];k=BASE3D.idz[l];
     Mx[i][j][k]=param[3*l];My[i][j][k]=param[3*l+1];Mz[i][j][k]=param[3*l+2];}
     nb_param=3*BASE3D.nbfx*BASE3D.nbfy*BASE3D.nbfz;
     l=0;
     for (i=0;i<BASE3D.nbfx;i++) for (j=0;j<BASE3D.nbfy;j++) for (k=0;k<BASE3D.nbfz;k++)
      {param[3*l]=Mx[i][j][k];param[3*l+1]=My[i][j][k];param[3*l+2]=Mz[i][j][k];l++;}
     free_dmatrix_3d(Mx);free_dmatrix_3d(My);free_dmatrix_3d(Mz);
     
     
     transfo=cr_transf3d_p(imref,nomfichiertrf);
     transfo->nb_param=nb_param;transfo->param=CALLOC(nb_param,double);
     transfo->typetrans=BSPLINE3D;
     transfo->resol=resolf;transfo->degre=func_type;
     for (l=0;l<nb_param;l++) transfo->param[l]=param[l]; 
    }
    save_transf_3d(transfo,nomfichres);
    free_transf3d(transfo);
   } 
 }


  /*caracteristique du champ*/
  if (champ_ini==NULL)
   {
    aff_log("CHAMP: ");carac_field_3d(champ);
   }
  else
   {
    aff_log("CHAMP TOTAL: ");carac_field_3d(champ);
    sub_field_3d(champ,champ_ini,champ);
    aff_log("CHAMP DEFOR: ");carac_field_3d(champ);
   }
 /*liberation memoire des variables dynamiques*/ 
  /*liberation memoire du champ*/
   free_field3d(champ);
   if (champ_ini!=NULL) free_field3d(champ_ini);
  /*liberation memoire des images temporaires*/
  if (renormalisation)
  {
   free_grphic3d(imtref);free_grphic3d(imtreca); 
  }
   free_grphic3d(imtres);
 
 
 /*affichage du temps de calcul total*/
 tfin=time(NULL);
 {
  imx_aff_temps_calcul(tdebut, tfin);
 }
 end_log();
 return(1);  
}
#endif



/*******************************************************************************
********************************************************************************
************************** MATCHING  FINE ************************************
********************************************************************************
*******************************************************************************/

/*******************************************************************************
**        matching_fine_3d                                                 
**                                                                        
**       recalage test de deux images 2D                                
*******************************************************************************/
/*
  renormalisation : TRUE si l'on veut effectuer une renormalisation des
    images (recommande)
*/
void matching_fine_3d()
{
 int im_ref,im_reca,im_res;
 char *quest[7],*nomfichres;
 int dist_type,inter_type,save_type,i;
 int renormalisation=TRUE;
 
 /*question sur les images a recaler*/
 im_reca=GET_PLACE3D(TEXT0030);
 im_ref=GET_PLACE3D(TEXT0031);
 im_res=GET_PLACE3D(TEXT0006);

 for(i=0;i<7;i++)
    quest[i]=CALLOC(80,char);
 
 /*question sur les parametres de la methode*/
 /*distance*/
  strcpy(quest[0],"Quadratic");
  strcpy(quest[1],"woods");
  strcpy(quest[2],"\0");
  dist_type=GETV_QCM("Distance",(char **)quest);
 /*interpolation*/
  inter_type=imx_query_interpolation_type_3D(dist_type);
   
 for(i=0;i<7;i++)
    free(quest[i]);

 /*question sur l'enregistrement du champ resultat*/
 nomfichres=quest_save_result_3d(&save_type);

 imx_matching_fine_3d(im_ref,im_reca,im_res,dist_type,inter_type,nomfichres,renormalisation);
 
 

 if (nomfichres)
  free(nomfichres);

 show_picture_3d(im_res);
 
}

/*******************************************************************************
**        imx_matching_fine_3d                                             
**                                                                       
**       recalage fine de 2 images 3D                                   
*******************************************************************************/
int imx_matching_fine_3d(int im_ref, int im_reca, int im_res, int dist_type, int inter_type, char *nomfichres, int renormalisation)
{
 grphic3d *imref,*imreca,*imres;
 
 imref=ptr_img_3d(im_ref);
 imreca=ptr_img_3d(im_reca);
 imres=ptr_img_3d(im_res);
 imx_matching_fine_3d_p(imref,imreca,imres,dist_type,inter_type,nomfichres,renormalisation);
 
  return(1);
}


/*******************************************************************************
**        imx_matching_fine_3d_p                                            
**                                                                        
**       recalage fine de 2 images 3D                                    
*******************************************************************************/
int imx_matching_fine_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int dist_type, int inter_type, char *nomfichres, int renormalisation)
{
 int    i,j,k,wdth,hght,dpth;
 grphic3d *imtref,*imtreca,*imtres;
 field3d *champ,*champt,*cht/*,*champr*/;
 vector3d ***data,***datat;
 InterpolationFct interpol;
 dist_func_t distance;
 char   *nomfichiertrf;
 int    tdebut,tfin;
 VP_histos * donnees_VP = NULL;
 ptr_distance_param dist_param=CALLOC(1, distance_param);

 wdth=imref->width;hght=imref->height;dpth=imref->depth;
#ifndef WIN32
 setvbuf(stdout, (char *)NULL, _IONBF, 0); //setbuf(stdout,NULL);
#endif
 tdebut=time(NULL);
 
 init_log(NULL);
 aff_log("\nRECALAGE FINE ");
 /*************mise en place des parametres du recalage******/
  distance=imx_choose_distance_fct(dist_type);
  interpol=imx_choose_interpolation_fct(inter_type); 

  aff_log("\n");
 aff_log("recalage de %s sur %s \n",(imreca->patient).name,(imref->patient).name);
 
 /*creation, allocation et initialisation des variables dynamiques*/
  /*creation des images temporaires*/
  if (renormalisation)
  {
   imtref=cr_grphic3d(imref);imx_copie_3d_p(imref,imtref);
   imtreca=cr_grphic3d(imreca);imx_copie_3d_p(imreca,imtreca);
  }
  else
  {
    imtref = imref;
    imtreca = imreca;
  }
   imtres=cr_grphic3d(imreca);imx_copie_param_3d_p(imreca,imtres);
  if (renormalisation)
  {
  /*normalisation de imreca par rapport a imref: resultat dans imtreca*/ 
   imx_norm_seuil_meanecty_3d_p(imreca,imtref,imtreca);
  }
   

 /*initialisation du champ de depart a la valeur d'un fichier*/
 {
  char  s[256];
  int   ans,err=0;
  
  nomfichiertrf=NULL;
  sprintf(s,"Use a transformation as initialisation?\n");
  ans=GET_EXIST(s,&err);
  if (ans==1)
  {
   char nomfichier[250];
   transf3d *transfo;
   
   strcpy(nomfichier,GET_FILE("*.trf",&err));
   transfo=load_transf_3d(nomfichier);
   champ=transf_to_field_3d(transfo,NULL,NULL);
   free_transf3d(transfo);
   nomfichiertrf=CALLOC(strlen(nomfichier),char);
   strcpy(nomfichiertrf,nomfichier);
  }
  else 
  {/*le champ est alloue et mis a zero*/
   champ=cr_field3d(wdth,hght,dpth);
   data=champ->raw;
   for (i=0;i<wdth;i++) for (j=0;j<hght;j++) for (k=0;k<dpth;k++) 
    data[i][j][k].x=data[i][j][k].y=data[i][j][k].z=0.0;
  }
 }

 data=champ->raw;
 champt=cr_field3d(wdth,hght,dpth);datat=champt->raw;
 
  /*recalage par fine tuning*/
  {
  double E1,E2;

  field3d  *grad;
  vector3d ***datag,***datac;
  double   gx,gy,gz,ng,dgx,dgy,dgz,ndg,diff;
  double   filt[10],sigma,norm,tfil;
  int      i1,i2,j1,j2,k1,k2,ic,jc,kc;
  int     nbiter,tfilt;
  
  nbiter=0;
  grad=cr_field3d(wdth,hght,dpth);
  datac=champ->raw;datag=grad->raw;  
  
  /*calcul de l'energie de depart*/
  interpol(imtreca,champ,imtres);
  dist_param->imreca=imtref;
  dist_param->imref=imtres;
  dist_param->champ=champ;
  dist_param->donnees_VP=donnees_VP;
//  E2=distance(imtref,imtres,champ,donnees_VP, NULL, NULL, NULL);
  E2=distance(dist_param);
  aff_log("Energie de depart: %f\n",E2);
  
  do
  {
   E1=E2;
   
   printf("calcul du gradient \r");
   /*calcul du gradient de imtres*/
    imx_gradient_3d_p(imtres,grad,2,4);
    
   printf("calcul du champ\r"); 
   /*calcul du champ residuel*/
   for (i=0;i<wdth;i++)
    for (j=0;j<hght;j++)
     for (k=0;k<dpth;k++)
     {
      diff=imtref->mri[i][j][k]-imtres->mri[i][j][k];
      gx=datag[i][j][k].x;gy=datag[i][j][k].y;gz=datag[i][j][k].z;
      ng=gx*gx+gy*gy+gz*gz;
      if (ng>0.0001)
      {
       dgx=diff*gx/ng;dgy=diff*gy/ng;dgz=diff*gz/ng;
       ndg=dgx*dgx+dgy*dgy+dgz*dgz;ndg=sqrt(ndg);
       if (ndg>1.0) {dgx=dgx/ndg;dgy=dgy/ndg;dgz=dgz/ndg;}
       data[i][j][k].x+=(float)dgx;
       data[i][j][k].y+=(float)dgy;
       data[i][j][k].z+=(float)dgz;
      }  
     } 
     
   printf("filtrage du champ \r");  
   /*filtrage gaussien du champ*/  
   tfilt=1;sigma=1.0;norm=0.0;
   for(i=0;i<=2*tfilt;i++)
    {
    filt[i]=(1/(2.0*PI*sigma))*exp(-1.0*((double)(i-tfilt))*((double)(i-tfilt))/(2.0*sigma*sigma));
    norm+=filt[i];
    }
   for(i=0;i<=2*tfilt;i++) filt[i]/=norm;
   
   for (i=tfilt;i<wdth-tfilt;i++)
    for (j=tfilt;j<hght-tfilt;j++)
     for (k=tfilt;k<dpth-tfilt;k++)
     {
      gx=0.0;gy=0.0;gz=0.0;
      i1=i-tfilt;j1=j-tfilt;k1=k-tfilt;
      i2=i+tfilt;j2=j+tfilt;k2=k+tfilt;
      for (ic=i1;ic<=i2;ic++)
       for (jc=j1;jc<=j2;jc++)
        for (kc=k1;kc<=k2;kc++)
    {
     tfil=filt[ic-i1]*filt[jc-j1]*filt[kc-k1];
     gx+=tfil*data[ic][jc][kc].x;
     gy+=tfil*data[ic][jc][kc].y;
     gz+=tfil*data[ic][jc][kc].z;
    }
      datat[i][j][k].x=(float)gx;
      datat[i][j][k].y=(float)gy;
      datat[i][j][k].z=(float)gz;   
     }
     
   printf("fin filtrage du champ \n");   
     cht=champ;champ=champt;champt=cht;
     datat=champt->raw;data=champ->raw;
  
     
   /*calcul de E2*/;
  interpol(imtreca,champ,imtres);
  dist_param->imreca=imtref;
  dist_param->imref=imtres;
  dist_param->champ=champ;
  dist_param->donnees_VP=donnees_VP;
//  E2=distance(imtref,imtres,champ,donnees_VP, NULL, NULL, NULL);
  E2=distance(dist_param);
  /*{
   char nomimagecoro[255],nomimagesagi[255],nomimageaxia[255],*txt;
   grphic3d  *imt2,*imt4;
   if (NBIM<10) {sprintf(nomimagecoro,"imcoro0%d.gif",NBIM);sprintf(nomimagesagi,"imsagi0%d.gif",NBIM);
        sprintf(nomimageaxia,"imaxia0%d.gif",NBIM);}
     else       {sprintf(nomimagecoro,"imcoro%d.gif",NBIM);sprintf(nomimagesagi,"imsagi%d.gif",NBIM);
        sprintf(nomimageaxia,"imaxia%d.gif",NBIM);}
   NBIM++;
   imx_copie_3d_p(imtres,ptr_img_3d(4));
   imt2=ptr_img_3d(2);
   imt4=ptr_img_3d(4);
   imt4->cutoff_max=imt2->cutoff_max;
   show_picture_3d(4);
   txt=CALLOC(MAX_TEXT_LEN,char);
   sprintf(txt,"RAFFINEMENT FINAL");

   affiche_text_obj(txt,imt4->img_coro,70,13,5);
   affiche_text_obj(txt,imt4->img_sagi,70,13,5);show_picture(imt4->img_tran);
   capture_picture(imt4->img_coro,nomimagecoro,GIF_FORMAT,1);
   capture_picture(imt4->img_sagi,nomimagesagi,GIF_FORMAT,1);
   
   imx_copie(imt4->img_tran,imt4->img_sagi);show_picture(imt4->img_sagi);
   affiche_text_obj(txt,imt4->img_sagi,70,13,5);show_picture(imt4->img_tran);
   capture_picture(imt4->img_sagi,nomimageaxia,GIF_FORMAT,1);
   free(txt);
  }*/    
  
  nbiter++;   
  printf("nbiter: %d  energie: %f  reduction: %f \r",nbiter,E2,(E1-E2)/E1*100);
  } while ((E1-E2)/E1*100>0.5);
  aff_log("nbiter: %d  energie: %f  reduction: %f \n",nbiter,E2,(E1-E2)/E1*100); 
   free_field3d(grad);
  } 
  
 /*application de la transformation aux images de depart
   Dans tous les cas interpolation qsinc3            */
   interpol=inter_qsinc3_3d;
   interpol(imreca,champ,imtres);
   imx_copie_3d_p(imtres,imres);
 
 aff_log("CHAMP: ");carac_field_3d(champ);
  
  /*enregistrement du champ resultat dans un fichier*/
   if (nomfichres!=NULL) 
   {
    transf3d *transfo;
    
    transfo=field_to_transf_3d(champ,imtref,imtreca);
    save_transf_3d(transfo,nomfichres);
    free_transf3d(transfo);
   } 

  /*liberation memoire des variables dynamiques*/ 
  /*liberation memoire du champ*/
   free_field3d(champ);free_field3d(champt);
  /*liberation memoire des images temporaires*/
  if (renormalisation)
  {
   free_grphic3d(imtref);free_grphic3d(imtreca); 
  }
   free_grphic3d(imtres);

 FREE(dist_param);

 tfin=time(NULL);
 /*affichage du temps de calcul total*/
 {
  imx_aff_temps_calcul(tdebut, tfin);
 }
 end_log();
 return(1);  

}

/*******************************************************************************
**        imx_creer_image_isotrope_3d_p(im, interpol, imres)
*/
/*!       calcule l'image isotrope en resolution 1.0 mm et de dimension max 256x256x256
**        qui contient l'image de depart
**
**  \param im : image de depart
**  \param interpol  : interpolation utilisee
**  \param imres : image resultat
**  \retval 0
**
*******************************************************************************/
int imx_creer_image_isotrope_3d_p(grphic3d *im, InterpolationFct interpol, grphic3d *imres)
{
 int wdth, hght, dpth;
 int new_wdth, new_hght, new_dpth;
 double *param=CALLOC(9,double);
 field3d *champ;
 float new_resol=0.0;
 int dim_max=128;
 grphic3d *im_temp=NULL;
 double fac_x=0.0, fac_y=0.0, fac_z=0.0;

 //on cherche une image en resolution 2 mm
 new_resol=2.0;
 wdth=im->width; hght=im->height; dpth=im->depth;
 new_wdth=MINI((int)ceil(wdth*im->dx/new_resol),dim_max);
 new_hght=MINI((int)ceil(hght*im->dy/new_resol),dim_max);
 new_dpth=MINI((int)ceil(dpth*im->dz/new_resol),dim_max);

 //creation de la nouvelle image
 imx_copie_param_3d_p(im, imres);
 free_mri_3d(imres->mri);
 imres->mri=cr_mri_3d(new_wdth,new_hght,new_dpth);
 imres->width=new_wdth; imres->height=new_hght; imres->depth=new_dpth;
 imres->dx=new_resol;
 imres->dy=new_resol;
 imres->dz=new_resol;

 //on va filtrer l'image dans les directions ou elle a une resolution superieure a 2 mm
 fac_x=2.0/im->dx; if (fac_x<=1.0) fac_x=0.0;
 fac_y=2.0/im->dy; if (fac_y<=1.0) fac_y=0.0;
 fac_z=2.0/im->dz; if (fac_z<=1.0) fac_z=0.0;
 if ((fac_x!=0.0)||(fac_y!=0.0)||(fac_z!=0.0))
 {
  im_temp=cr_grphic3d(im);
  imx_filtre_gaussien_anisotrope_3d_p(im, im_temp, 9, 0.425*fac_x, 0.425*fac_y, 0.425*fac_z);
 }
 else im_temp=im;

 //calcul de son contenu par interpolation
 param[6]=new_wdth/2; param[7]=new_hght/2; param[8]=new_dpth/2;
 champ=cr_field3d(new_wdth,new_hght,new_dpth);
 rigid_to_field_3d(9, param, champ, imres, im_temp);
 interpol(im_temp, champ, imres);

 //liberation memoire
 free_field3d(champ);
 FREE(param);
 if (im_temp!=im) free_grphic3d(im_temp);

 return 0;
}

void imx_creer_image_isotrope_3d(int im_deb, int inter_type, int im_res)
{
 grphic3d *imdeb=NULL, *imres=NULL;
 InterpolationFct interpol;

 imdeb=ptr_img_3d(im_deb);
 imres=ptr_img_3d(im_res);
 interpol=imx_choose_interpolation_fct(inter_type);

 imx_creer_image_isotrope_3d_p(imdeb,interpol,imres);
}

void creer_image_isotrope_3d()
{
 int im_deb, im_res, inter_type;

 /*question sur l'image*/
 im_deb=GET_PLACE3D(TEXT0030);
 /*question sur l'image*/
 im_res=GET_PLACE3D(TEXT0030);
 /*interpolation*/
 inter_type=imx_query_interpolation_type_3D(0);

 imx_creer_image_isotrope_3d(im_deb,inter_type,im_res);

 show_picture_3d(im_res);
}

/*------------------------------------------------*/
/*!
**  \brief affiche les parametres en fonction
**   du type de transformation
**
**  \param  transfo : type de transformation
**  \param  param : parametres
*/
/*------------------------------------------------*/
void imx_aff_param(TRANSFO3D transfo, double *param)
{
 switch (transfo)
 {
  case RIGID3D : imx_aff_param_rigid(param);break;
  case RIGIDZOOM3D : imx_aff_param_rigid_zoom(param);break;
  case AFFINE3D : imx_aff_param_affine(param);break;
  case RIGIDGLOBALZOOM3D : imx_aff_param_rigid_global_zoom(param);break;
  case AFFINEDECOUPLE3D : imx_aff_param_affine_decouple(param);break;
  default : imx_aff_param_rigid(param);break;
 }
}

/*------------------------------------------------*/
/*!
**  \brief affiche les parametres rigides
**
**  \param  param : parametres
*/
/*------------------------------------------------*/
void imx_aff_param_rigid(double *param)
{
 aff_log("rotation: %.2f  trans: %.2f   centre: %.2f \n",param[0]*180.0/PI,param[3],param[6]);
 aff_log("          %.2f         %.2f           %.2f \n",param[1]*180.0/PI,param[4],param[7]);
 aff_log("          %.2f         %.2f           %.2f \n",param[2]*180.0/PI,param[5],param[8]);
}

/*------------------------------------------------*/
/*!
**  \brief affiche les parametres rigides et de zoom
**
**  \param  param : parametres
*/
/*------------------------------------------------*/

void imx_aff_param_rigid_zoom(double *param)
{
 aff_log("rotation: %.2f zoom: %.2f  trans: %.2f   centre: %.2f \n",param[0]*180.0/PI,param[3],param[6],param[9]);
 aff_log("          %.2f zoom: %.2f  trans: %.2f   centre: %.2f \n",param[1]*180.0/PI,param[4],param[7],param[10]);
 aff_log("          %.2f zoom: %.2f  trans: %.2f   centre: %.2f \n",param[2]*180.0/PI,param[5],param[8],param[11]);
}

/*------------------------------------------------*/
/*!
**  \brief affiche les parametres affine
**
**  \param  param : parametres
*/
/*------------------------------------------------*/
void imx_aff_param_affine(double *param)
{
 aff_log("matrice: %.2f  %.2f  %.2f trans:  %.2f  centre:  %.2f \n",param[0],param[1],param[2],param[9],param[12]);
 aff_log("         %.2f  %.2f  %.2f         %.2f           %.2f \n",param[3],param[4],param[5],param[10],param[13]);
 aff_log("         %.2f  %.2f  %.2f         %.2f           %.2f \n",param[6],param[7],param[8],param[11],param[14]);
}

/*------------------------------------------------*/
/*!
**  \brief affiche les parametres rigide + zoom global isotrope
**
**  \param  param : parametres
*/
/*------------------------------------------------*/
void imx_aff_param_rigid_global_zoom(double *param)
{
 aff_log("rotation: %.2f zoom: %.2f  trans: %.2f   centre: %.2f \n",param[0]*180.0/PI,param[3],param[4],param[7]);
 aff_log("          %.2f             trans: %.2f   centre: %.2f \n",param[1]*180.0/PI,param[5],param[8]);
 aff_log("          %.2f             trans: %.2f   centre: %.2f \n",param[2]*180.0/PI,param[6],param[9]);
}


/*------------------------------------------------*/
/*!
**  \brief affiche les parametres affines decouples
**
**  \param  param : parametres
*/
/*------------------------------------------------*/
void imx_aff_param_affine_decouple(double *param)
{
 aff_log("rotation: %.2f trans: %.2f zoom: %.2f skewyx: %.2f centre: %.2f \n",param[0]*180.0/PI,param[3],param[6],param[9],param[12]);
 aff_log("          %.2f trans: %.2f zoom: %.2f skewzx: %.2f centre: %.2f \n",param[1]*180.0/PI,param[4],param[7],param[10],param[13]);
 aff_log("          %.2f trans: %.2f zoom: %.2f skewzy: %.2f centre: %.2f \n",param[2]*180.0/PI,param[5],param[8],param[11],param[14]);
}

/*------------------------------------------------*/
/*!
**  \brief definit le facteur a applique a la precision standard
**
**  \param  matchPrecision : enum definissant la precision voulue
**  \retval : la facteur determine
*/
/*------------------------------------------------*/
double get_facteur_precision(eMatchPrecision matchPrecision)
{
 double larPrec=1.0;

 switch(matchPrecision)
 {
   case NORMAL        : larPrec=1.0;break;
   case PRECIS        : larPrec = 1e-1;break;
   case TRES_PRECIS   : larPrec = 1e-2;break;
   case PRECISION_MAX : larPrec = 1e-4;break;
   default : fprintf(stderr,"SHOULD NOT SEE ME!!get_facteur_precision\n");
 }

 return larPrec;
}

/*------------------------------------------------*/
/*!
**  \brief definit le facteur a applique a la precision standard
**
**  \param  FieldOfResearch : enum definissant l'intervalle de recherche voulu
**  \retval : la facteur determine
*/
/*------------------------------------------------*/
double get_facteur_intervalle_recherche(eResearchInterv FieldOfResearch)
{
 double lar=1.0;

 switch(FieldOfResearch)
 {
   case SMALL  : break;
   case MEDIUM : lar *=2;break;
   case LARGE  : lar *=4;break;

   default : fprintf(stderr,"SHOULD NOT SEE ME!!get_facteur_intervalle_recherche\n");
 }

 return lar;
}

/*------------------------------------------------*/
/*!
**  \brief fixe les bornes de centre
**
**  \param  min_param : bornes min de centre
**  \param  max_param : bornes max de centre
**  \param  prec_param : precision de centre
**  \param  param : parametres de centre de depart
*/
/*------------------------------------------------*/
void init_bornes_matching_centre(double *min_param, double *max_param, double *prec_param, double *param, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision)
{
 double lar;
 int k;

 lar=0.0; //centre
 for (k=0;k<3;k++)
 {min_param[k]=param[k]-lar;max_param[k]=param[k]+lar;prec_param[k]=0.0;}

}

/*------------------------------------------------*/
/*!
**  \brief fixe les bornes de rotation
**
**  \param  min_param : bornes min de rotation
**  \param  max_param : bornes max de rotation
**  \param  prec_param : precision de rotation
**  \param  param : parametres de rotation de depart
*/
/*------------------------------------------------*/
void init_bornes_matching_rot(double *min_param, double *max_param, double *prec_param, double *param, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision, double dx, double dy, double dz)
{
 double lar, larPrec=1.0;
 int k;
 lar=0.05; //angles
 
 lar=lar*get_facteur_intervalle_recherche(FieldOfResearch);
 larPrec=get_facteur_precision(matchPrecision);

 for (k=0;k<3;k++)
 {min_param[k]=param[k]-lar;max_param[k]=param[k]+lar;} 

 //definition de la precision des parametres en fonction de la precision des images
 prec_param[0]=dx/(2.0*BRAIN_RADIUS)*larPrec;
 prec_param[1]=dy/(2.0*BRAIN_RADIUS)*larPrec;
 prec_param[2]=dz/(2.0*BRAIN_RADIUS)*larPrec;
 
}

/*------------------------------------------------*/
/*!
**  \brief fixe les bornes de translation
**
**  \param  min_param : bornes min de translation
**  \param  max_param : bornes max de translation
**  \param  prec_param : precision de translation
**  \param  param : parametres de translation de depart
*/
/*------------------------------------------------*/
void init_bornes_matching_trans(double *min_param, double *max_param, double *prec_param, double *param, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision, double dx, double dy, double dz)
{
 double lar, larPrec=1.0;
 int k;
 
 lar=5.0; //translations

 lar=lar*get_facteur_intervalle_recherche(FieldOfResearch);
 larPrec=get_facteur_precision(matchPrecision);

 for (k=0;k<3;k++)
 {min_param[k]=param[k]-lar;max_param[k]=param[k]+lar;} 

 //definition de la precision des parametres en fonction de la precision des images
 prec_param[0]=dx/2.0*larPrec;
 prec_param[1]=dy/2.0*larPrec;
 prec_param[2]=dz/2.0*larPrec;

}

/*------------------------------------------------*/
/*!
**  \brief fixe les bornes de zoom
**
**  \param  min_param : bornes min de zoom
**  \param  max_param : bornes max de zoom
**  \param  prec_param : precision de zoom
**  \param  param : parametres de zoom de depart 
*/
/*------------------------------------------------*/
void init_bornes_matching_zoom(double *min_param, double *max_param, double *prec_param, double *param, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision, double dx, double dy, double dz)
{
 double lar, larPrec=1.0;
 int k;

 lar=0.05; //facteurs de zoom

 lar=(1.0-lar*get_facteur_intervalle_recherche(FieldOfResearch));
 larPrec=get_facteur_precision(matchPrecision);

 for (k=0;k<3;k++)
 {min_param[k]=param[k]*lar;max_param[k]=param[k]/lar;} 


 //definition de la precision des parametres en fonction de la precision des images
 prec_param[0]=dx/(2.0*BRAIN_RADIUS)*larPrec;
 prec_param[1]=dy/(2.0*BRAIN_RADIUS)*larPrec;
 prec_param[2]=dz/(2.0*BRAIN_RADIUS)*larPrec;

}

/*------------------------------------------------*/
/*!
**  \brief determine les bornes du recalage
**   en fonction du recalage choisi
**
**  \param  transfo : le type transformation choisi
**  \param  min_param : bornes min des parametres
**  \param  max_param : bornes max des parametres
**  \param  prec_param : precision des parametres
**  \param  param : parametres de depart du recalage
*/
/*------------------------------------------------*/
void init_bornes_matching_affine(double *min_param, double *max_param, double *prec_param, double *param, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision, double dx, double dy, double dz)
{
 double lar, larPrec=1.0;
 int k;

 lar=0.3; //parametres de la matrice affine

 lar=lar*get_facteur_intervalle_recherche(FieldOfResearch);
 larPrec=get_facteur_precision(matchPrecision);

 for (k=0;k<9;k++)
 {min_param[k]=param[k]-lar;max_param[k]=param[k]+lar;} 

 //definition de la precision des parametres en fonction de la precision des images
 prec_param[0]=dx/(2.0*BRAIN_RADIUS)*larPrec;
 prec_param[1]=dy/(2.0*BRAIN_RADIUS)*larPrec;
 prec_param[2]=dz/(2.0*BRAIN_RADIUS)*larPrec;

}

/*------------------------------------------------*/
/*!
**  \brief determine les bornes du recalage
**   en fonction du recalage choisi
**
**  \param  transfo : le type transformation choisi
**  \param  min_param : bornes min des parametres
**  \param  max_param : bornes max des parametres
**  \param  prec_param : precision des parametres
**  \param  param : parametres de depart du recalage
*/
/*------------------------------------------------*/
void init_bornes_matching_global_zoom(double *min_param, double *max_param, double *prec_param, double *param, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision, double dx, double dy, double dz)
{
 double lar, larPrec=1.0;

 lar=0.05; //facteurs de zoom

 lar=(1.0-lar*get_facteur_intervalle_recherche(FieldOfResearch));
 larPrec=get_facteur_precision(matchPrecision);

 min_param[0]=param[0]*lar;max_param[0]=param[0]/lar; 

 //definition de la precision des parametres en fonction de la precision des images
 prec_param[0]=dx/(2.0*BRAIN_RADIUS)*larPrec;
 
}

/*------------------------------------------------*/
/*!
**  \brief determine les bornes du recalage
**   en fonction du recalage choisi
**
**  \param  transfo : le type transformation choisi
**  \param  min_param : bornes min des parametres
**  \param  max_param : bornes max des parametres
**  \param  prec_param : precision des parametres
**  \param  param : parametres de depart du recalage
*/
/*------------------------------------------------*/
void init_bornes_matching_skew(double *min_param, double *max_param, double *prec_param, double *param, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision, double dx, double dy, double dz)
{
 double lar, larPrec=1.0;
 int k;

 lar=0.1; //facteurs de skew

 lar=lar*get_facteur_intervalle_recherche(FieldOfResearch);
 larPrec=get_facteur_precision(matchPrecision);

 for (k=0;k<3;k++)
 {min_param[k]=param[k]-lar;max_param[k]=param[k]+lar;} 

 //definition de la precision des parametres en fonction de la precision des images
 prec_param[0]=dx/(2.0*BRAIN_RADIUS)*larPrec;
 prec_param[1]=dy/(2.0*BRAIN_RADIUS)*larPrec;
 prec_param[2]=dz/(2.0*BRAIN_RADIUS)*larPrec;

}

/*------------------------------------------------*/
/*!
**  \brief determine les bornes du recalage
**   en fonction du recalage choisi
**
**  \param  transfo : le type transformation choisi
**  \param  min_param : bornes min des parametres
**  \param  max_param : bornes max des parametres
**  \param  prec_param : precision des parametres
**  \param  FieldOfResearch : Sert �initialiser max_param et min_param
**                            Il sert donc a definir l'intervalle de recherche de la transformation de recalage
**                            Il peut etre mis �3 valeurs SMALL, MEDIUM et LARGE
**                            On s'arrange pour que LARGE corresponde aux variations de transformation maximum
**                            dans les cas connus (typiquement +/- 30 degres en rotation et +/- 5cm en translation
**                            MEDIUM et SMALL correspondent alors respectivement a ces intervalles divises par 2 et par 4
**  \param  matchPrecision : Sert �inititaliser prec_param
**                           Il sert donc a definir a partir quelle variation minimale des parametres
**                           on considere que le recalage ne varira plus de fa�n significative
**                           Il peut etre mis �NORMAL (precision d'1/2 voxel), PRECIS (1/20 voxel)
**                           TRES_PRECIS (1/200 voxel), PRECISION_MAX (1/20000 voxel)
**  \param  param : parametres de depart du recalage
*/
/*------------------------------------------------*/
void init_bornes_matching(TRANSFO3D transfo, double *min_param, double *max_param, double *prec_param, double *param, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision, grphic3d *imref, grphic3d *imreca)
{
 double dminRef, dminReca, res;

 dminRef=MINI(imref->dz, MINI(imref->dx, imref->dy));
 dminReca=MINI(imreca->dz, MINI(imreca->dx, imreca->dy));
 res=MAXI(dminRef, dminReca);  
 
 switch (transfo)
 {
  case RIGID3D :    init_bornes_matching_rot(min_param, max_param, prec_param, param,FieldOfResearch,matchPrecision,res,res,res);
                    init_bornes_matching_trans(&(min_param[3]), &(max_param[3]), &(prec_param[3]), &(param[3]),FieldOfResearch,matchPrecision,res,res,res);
                    init_bornes_matching_centre(&(min_param[6]), &(max_param[6]), &(prec_param[6]), &(param[6]),FieldOfResearch,matchPrecision);
                    break;
  case RIGIDZOOM3D :init_bornes_matching_rot(min_param, max_param, prec_param, param,FieldOfResearch,matchPrecision,res,res,res);
                    init_bornes_matching_zoom(&(min_param[3]), &(max_param[3]), &(prec_param[3]), &(param[3]),FieldOfResearch,matchPrecision,res,res,res);
                    init_bornes_matching_trans(&(min_param[6]), &(max_param[6]), &(prec_param[6]), &(param[6]),FieldOfResearch,matchPrecision,res,res,res);
                    init_bornes_matching_centre(&(min_param[9]), &(max_param[9]), &(prec_param[9]), &(param[9]),FieldOfResearch,matchPrecision);
                    break;
  case AFFINE3D :   init_bornes_matching_affine(min_param, max_param, prec_param, param,FieldOfResearch,matchPrecision,res,res,res);
                    init_bornes_matching_trans(&(min_param[9]), &(max_param[9]), &(prec_param[9]), &(param[9]),FieldOfResearch,matchPrecision,res,res,res);
                    init_bornes_matching_centre(&(min_param[12]), &(max_param[12]), &(prec_param[12]), &(param[12]),FieldOfResearch,matchPrecision);
                    break;
  case RIGIDGLOBALZOOM3D : init_bornes_matching_rot(min_param, max_param, prec_param, param,FieldOfResearch,matchPrecision,res,res,res);
                    init_bornes_matching_global_zoom(&(min_param[3]), &(max_param[3]), &(prec_param[3]), &(param[3]),FieldOfResearch,matchPrecision,res,res,res);
                    init_bornes_matching_trans(&(min_param[4]), &(max_param[4]), &(prec_param[4]), &(param[4]),FieldOfResearch,matchPrecision,res,res,res);
                    init_bornes_matching_centre(&(min_param[7]), &(max_param[7]), &(prec_param[7]), &(param[7]),FieldOfResearch,matchPrecision);
                    break;
  case AFFINEDECOUPLE3D :init_bornes_matching_rot(min_param, max_param, prec_param, param,FieldOfResearch,matchPrecision,res,res,res);
                    init_bornes_matching_trans(&(min_param[3]), &(max_param[3]), &(prec_param[3]), &(param[3]),FieldOfResearch,matchPrecision,res,res,res);
                    init_bornes_matching_zoom(&(min_param[6]), &(max_param[6]), &(prec_param[6]), &(param[6]),FieldOfResearch,matchPrecision,res,res,res);
                    init_bornes_matching_skew(&(min_param[9]), &(max_param[9]), &(prec_param[9]), &(param[9]),FieldOfResearch,matchPrecision,res,res,res);
                    init_bornes_matching_centre(&(min_param[12]), &(max_param[12]), &(prec_param[12]), &(param[12]),FieldOfResearch,matchPrecision);
                    break;
  default :         init_bornes_matching_rot(min_param, max_param, prec_param, param,FieldOfResearch,matchPrecision,res,res,res);
                    init_bornes_matching_trans(&(min_param[3]), &(max_param[3]), &(prec_param[3]), &(param[3]),FieldOfResearch,matchPrecision,res,res,res);
                    init_bornes_matching_centre(&(min_param[6]), &(max_param[6]), &(prec_param[6]), &(param[6]),FieldOfResearch,matchPrecision);
                    break;
 }
 
}


/*------------------------------------------------*/
/*!
**  \brief Determine les parametres suivant lesquels on optimise pas
**
**  \param  transfo : le type transformation choisi
**  \param  min_param : bornes min des parametres
**  \param  max_param : bornes max des parametres
**  \param  prec_param : precision des parametres
**  \param  param : parametres de depart du recalage
*/
/*------------------------------------------------*/
void init_bornes_matching_hybride(TRANSFO3D transfo, double *min_param, double *max_param, double *prec_param, double *param)
{
 int i;
 switch (transfo)
 {
  case RIGID3D :    for (i=0;i<3;i++)
                                        if (_mtch_3d_hybride_param[i]==0)
                                                {min_param[i]=max_param[i]=param[i];}
                                                
                                        for (i=3;i<9;i++)
                                        if (_mtch_3d_hybride_param[i+3]==0)
                                                {min_param[i]=max_param[i]=param[i];}
                                      break;
    
    case RIGIDGLOBALZOOM3D :
                                        for (i=0;i<3;i++)
                                        if (_mtch_3d_hybride_param[i]==0)
                                                {min_param[i]=max_param[i]=param[i];}
                                                
                                        if (_mtch_3d_hybride_param[15]==0)
                                                {min_param[3]=max_param[3]=param[3];}
                                        
                                        for (i=4;i<10;i++)
                                        if (_mtch_3d_hybride_param[i+2]==0)
                                                {min_param[i]=max_param[i]=param[i];}
                                        
                                        break;      
                                                
  case RIGIDZOOM3D : for (i=0;i<12;i++)
                                        if (_mtch_3d_hybride_param[i]==0)
                                                {min_param[i]=max_param[i]=param[i];}
                    break;
  default :         for (i=0;i<12;i++)
                                        if (_mtch_3d_hybride_param[i]==0)
                                                {min_param[i]=max_param[i]=param[i];}
                    break;
 }


 
}

int imx_matching_preprocessing_3d_p(grphic3d *imref,grphic3d *imreca, grphic3d **imtref, grphic3d **imtreca, grphic3d **imtres, transf3d **preProcessingRefTrf, transf3d **preProcessingRecaTrf, dist_func_t  distance)
{
 transf3d *resulttrfReca = NULL;
 grphic3d *im_temp=NULL;

 if (imref->dx==0) imref->dx=1.0; if (imref->dy==0) imref->dy=1.0; if (imref->dz==0) imref->dz=1.0;
 if (imreca->dx==0) imreca->dx=1.0; if (imreca->dy==0) imreca->dy=1.0; if (imreca->dz==0) imreca->dz=1.0;
 //on s'assure que les champs des images sont a jour
 imx_inimaxminpixel_3d_p(imref);
 imx_inimaxminpixel_3d_p(imreca);

 //recadrage eventuel de l'image ref et de l'image reca au cas ou celles-ci seraient inutilement grandes
 //!!!valable uniquement pour le recalage de cerveaux
 //dans le cadre general: a faire un cadre plus modulaire de preprocessing!
 (*preProcessingRefTrf) = imx_matching_recadrage_cerveau_3d_p(imref, imtref);
 resulttrfReca = imx_matching_recadrage_cerveau_3d_p(imreca, imtreca);
 if (resulttrfReca)
 {
  (*preProcessingRecaTrf) = cr_copy_transf3d(resulttrfReca);
  imx_inverser_transformation_3d_p(resulttrfReca, (*preProcessingRecaTrf));
 }
 else (*preProcessingRecaTrf)=NULL;

 if ((distance==erreur_quad_3d)||(distance==erreur_quad_robust_geman_3d)||(distance==erreur_quad_robust2_3d))
 {
  im_temp=cr_grphic3d(*imtreca);
  //normalisation de imreca par rapport a imref: resultat dans imtreca
  imx_norm_seuil_meanecty_3d_p(*imtreca,*imtref,im_temp);
  if (*imtreca==imreca) {*imtreca=cr_grphic3d(im_temp);imx_copie_3d_p(im_temp,*imtreca); }
  else {imx_copie_3d_p(im_temp,*imtreca);}
  free_grphic3d(im_temp);
 }

 imx_inimaxminpixel_3d_p(*imtref);
 imx_inimaxminpixel_3d_p(*imtreca);

 //on ne recale que les images positives
 if ((*imtref)->min_pixel<0) imx_add_coe_3d_p(*imtref, -(*imtref)->min_pixel*(*imtref)->rcoeff, *imtref);
 if ((*imtreca)->min_pixel<0) imx_add_coe_3d_p(*imtreca, -(*imtreca)->min_pixel*(*imtreca)->rcoeff, *imtreca);

#ifdef WIN32
 //on s'impose une limite de taille pour les images a recaler
 imx_reduire_image_recalage(*imtref);
 imx_reduire_image_recalage(*imtreca);
#endif

 //creation de l'image ou placer les resultats de transformations
 if (imtres) *imtres=cr_grphic3d(*imtref);

 if (resulttrfReca) free_transf3d(resulttrfReca);

 return 0;
}

transf3d * imx_matching_recadrage_cerveau_3d_p(grphic3d *imdeb, grphic3d **imres)
{
  transf3d *resulttrf = NULL;
  double cgx   =0.0;
  double cgy   =0.0;
  double cgz   =0.0;
  double taille_boite_utile = 120.0; // rayon util de l'image qu'on veut considerer (en mm)
  int width,height,depth;
  int wdthDeb, wdthEnd, hghtDeb, hghtEnd, dpthDeb, dpthEnd;
  int i, j, k;
  TYPEMRI3D ***mri;
  char logString[1024];

  if (!imdeb && !imres) return NULL;

  width = imdeb->width;
  height= imdeb->height;
  depth = imdeb->depth;
  mri = imdeb->mri;

  imx_find_cg_with_seuil_3d_p(imdeb,&cgx,&cgy,&cgz);

  //on verifie si l'image est plus grande
  //qu'une boite de taille_boite_utile mm autour du centre de gravite
  wdthDeb=(int)floor(cgx-taille_boite_utile/imdeb->dx+0.5); wdthEnd=(int)floor(cgx+taille_boite_utile/imdeb->dx+0.5);
  hghtDeb=(int)floor(cgy-taille_boite_utile/imdeb->dy+0.5); hghtEnd=(int)floor(cgy+taille_boite_utile/imdeb->dy+0.5);
  dpthDeb=(int)floor(cgz-taille_boite_utile/imdeb->dz+0.5); dpthEnd=(int)floor(cgz+taille_boite_utile/imdeb->dz+0.5);
  if ((wdthDeb<=0)&&(wdthEnd>=width)&&(hghtDeb<=0)&&(hghtEnd>=height)&&(dpthDeb<=0)&&(dpthEnd>=depth))
  {
   (*imres)=imdeb;
   return NULL;
  }

  // la partie utile de l'image est petite par
  // rapport a la taille de l'image -> crop
  sprintf(logString, "PreProcessing !: %d<X<%d; %d<Y<%d; %d<Z<%d\n",
          wdthDeb, wdthEnd, hghtDeb, hghtEnd, dpthDeb, dpthEnd);
  aff_log(logString);

  wdthDeb=MAXI(0, wdthDeb); wdthEnd=MINI(width, wdthEnd);
  hghtDeb=MAXI(0, hghtDeb); hghtEnd=MINI(height, hghtEnd);
  dpthDeb=MAXI(0, dpthDeb); dpthEnd=MINI(depth, dpthEnd);

  (*imres)=cr_grphic3d_modif(wdthEnd-wdthDeb, hghtEnd-hghtDeb, dpthEnd-dpthDeb, 0.0, 1.0, 0);
  imx_copie_param_3d_p(imdeb, *imres);
  (*imres)->width=wdthEnd-wdthDeb; (*imres)->height=hghtEnd-hghtDeb, (*imres)->depth=dpthEnd-dpthDeb;

  for (i=wdthDeb;i<wdthEnd;i++)
  {
   for (j=hghtDeb;j<hghtEnd;j++)
   {
    for (k=dpthDeb;k<dpthEnd;k++)
    {
     (*imres)->mri[i-wdthDeb][j-hghtDeb][k-dpthDeb]=mri[i][j][k];
    }
   }
  }

  resulttrf = cr_transf3d_p(imdeb,NULL);
  resulttrf->param = CALLOC(15,double);
  resulttrf->nb_param = 9;
  resulttrf->typetrans = RIGID3D;
  resulttrf->param[3]=-wdthDeb*imdeb->dx;
  resulttrf->param[4]=-hghtDeb*imdeb->dy;
  resulttrf->param[5]=-dpthDeb*imdeb->dz;

  return resulttrf;
}

void imx_matching_recadrage_cerveau_3d(int im_deb, int im_res)
{
 grphic3d *imdeb, *imres, *imtres;

 imdeb=ptr_img_3d(im_deb);
 imres=ptr_img_3d(im_res);

 imx_matching_recadrage_cerveau_3d_p(imdeb, &imtres);

 imx_copie_3d_p(imtres, imres);
}

void matching_recadrage_cerveau_3d()
{
 int im_deb,im_res;

 /*question sur les images a recaler*/
 im_deb=GET_PLACE3D(TEXT0031);
 im_res=GET_PLACE3D(TEXT0006);

 imx_matching_recadrage_cerveau_3d(im_deb, im_res);

 show_picture_3d(im_res);
}

transf3d *imx_trouver_transf_initiale_matching_3d_p(grphic3d *imref, grphic3d *imreca, transf3d *inittrf)
{
 transf3d *transfoini=NULL;
 int k;
 double *res;

 transfoini=cr_transf3d_p(imref,NULL);

 if (inittrf!=NULL)
 {
  //verifie la compatibilite de la transfo avec l'image ref
  if ((inittrf->width!=imref->width)||(inittrf->height!=imref->height)||(inittrf->depth!=imref->depth)||
      (inittrf->dx!=imref->dx)||(inittrf->dy!=imref->dy)||(inittrf->dz!=imref->dz))
  {
   fprintf(stderr,"Pb avec la transfo initiale!!: initialisation par centres de gravite\n");
  }
  else
  {
   copie_transf_3d(inittrf, transfoini);
   aff_log("\n INITIALISATION AVEC UNE TRANSFORMATION\n");
   return transfoini;   
  } 
 }

 // si pas d'initialisation par fichier, alignement centre de gravite
 aff_log("\n ALIGNEMENT DES CENTRES DE GRAVITE \n");
 res=gravite_to_field_3d(imref,imreca,NULL);
 transfoini->param=CALLOC(9, double);
 transfoini->nb_param=9;
 transfoini->typetrans=RIGID3D;
 for (k=0;k<6;k++) transfoini->param[k]=res[k];
 gravite_to_rigid_3d(transfoini->param);
 free(res);

 return transfoini;
}

/*! \ingroup TransfoFonction @{ */
/*******************************************************************************
**     rigid_to_field_champ_3d(nb_param,param,points_choisis,champ,champ,imref,imreca)
*/
/*!     remplit un champ dense 3D pour une transformation rigide
**  applique autour du centre (xc,yc) et avec une translation (tx,ty)
**      Les parametres dans le tableau param, sont ordonnes ainsi
**  tetax tetay tetaz tx ty tz xc yc zc
**  \param nb_param : inutilise
**  \param points_choisis  : points transformes
**  \param param  : tableau de param
**  \param champ : le champ (E/S)
**  \param imref : utile pour creer un champ d'une transfo anisotrope (espace imreca -> espace imref)
**  \param imref : utile pour creer un champ d'une transfo anisotrope (espace imreca -> espace imref)
**  \note : si imref et imreca == NULL : la champ est genere tel quel a partir de param.
**  \retval : 1
*******************************************************************************/
int rigid_to_field_champ_3d(int nb_param,  const double *param, field3d* points_choisis, field3d *champ,grphic3d *imref, grphic3d *imreca)
{
  int i,j,k,wdth,hght,dpth;
  double x,y,z;
  double tetax,tetay,tetaz,tx,ty,tz,xc,yc,zc;
  double cxx,sxx,cyy,syy,czz,szz;
  double a11,a12,a13,a21,a22,a23,a31,a32,a33;
  double tpx11,tpx21,tpx31,tpy12,tpy22,tpy32;
  vector3d ***data;
  transf3d* transfOriginal   = NULL;
  transf3d* transfAnisotrope = NULL;
  vector3d ***pts_choisis;
  
 if (points_choisis==NULL)
 {
  return rigid_to_field_3d(nb_param,  param, champ, imref, imreca);
 }

  transfOriginal = cr_transf3d(0,0,0,NULL);
  transfOriginal->param = CALLOC(15,double);
  transfOriginal->typetrans = RIGID3D;
  transfOriginal->nb_param = nb_param;
  for (i=0;i<nb_param;i++)
    transfOriginal->param[i] = param[i];


  if (imref && imreca)
    transfAnisotrope = cr_transf_anisotrope_3d(transfOriginal,imref,imreca);

  wdth=points_choisis->width; hght=points_choisis->height; dpth=points_choisis->depth;
  data=champ->raw;
  pts_choisis=points_choisis->raw;
  
  if (transfAnisotrope && transfAnisotrope->typetrans==AFFINE3D)
  {
    a11=transfAnisotrope->param[0];a12=transfAnisotrope->param[1];a13=transfAnisotrope->param[2];
    a21=transfAnisotrope->param[3];a22=transfAnisotrope->param[4];a23=transfAnisotrope->param[5];
    a31=transfAnisotrope->param[6];a32=transfAnisotrope->param[7];a33=transfAnisotrope->param[8];
    tx=transfAnisotrope->param[9];ty=transfAnisotrope->param[10];tz=transfAnisotrope->param[11];
    xc=transfAnisotrope->param[12];yc=transfAnisotrope->param[13];zc=transfAnisotrope->param[14];
  }
  else
  {
    tetax=param[0];tetay=param[1];tetaz=param[2];
    tx=param[3];ty=param[4];tz=param[5];
    xc=param[6];yc=param[7];zc=param[8];
    cxx=cos(tetax);sxx=sin(tetax);cyy=cos(tetay);syy=sin(tetay);czz=cos(tetaz);szz=sin(tetaz);
    a11=cyy*czz;a12=-cyy*szz;a13=-syy;
    a21=cxx*szz-sxx*syy*czz;
    a22=cxx*czz+sxx*syy*szz;
    a23=-sxx*cyy;
    a31=sxx*szz+cxx*syy*czz;
    a32=sxx*czz-cxx*syy*szz;
    a33=cxx*cyy;
  }

  a11=a11-1.0;a22=a22-1.0;a33=a33-1.0;
  for (i=0;i<wdth;i++)
   {
    for (j=0;j<hght;j++)
    {
     for (k=0;k<dpth;k++)
      {
       x=(double)pts_choisis[i][j][k].x-xc;tpx11=a11*x;tpx21=a21*x;tpx31=a31*x;
       y=(double)pts_choisis[i][j][k].y-yc;tpy12=a12*y;tpy22=a22*y;tpy32=a32*y;
       z=(double)pts_choisis[i][j][k].z-zc;
       data[i][j][k].x=(float)(tpx11+tpy12+a13*z+tx);
       data[i][j][k].y=(float)(tpx21+tpy22+a23*z+ty);
       data[i][j][k].z=(float)(tpx31+tpy32+a33*z+tz);
      }
    }
   }
  if (transfOriginal) free_transf3d(transfOriginal);
  if (transfAnisotrope) free_transf3d(transfAnisotrope);  

  return(1); 
}
/*! @} */

/*! \ingroup TransfoFonction @{ */
/*******************************************************************************
**     rigidz_to_field_champ_3d(nb_param,param,points_choisis,champ,champ,imref,imreca)
*/
/*!     remplit un champ dense 3D pour une transformation rigide + zoom
**  applique autour du centre (xc,yc) et avec une translation (tx,ty)
**      Les parametres dans le tableau param, sont ordonnes ainsi
**  tetax tetay tetaz tx ty tz xc yc zc
**  \param nb_param : inutilise
**  \param points_choisis  : points transformes
**  \param param  : tableau de param
**  \param champ : le champ (E/S)
**  \param imref : utile pour creer un champ d'une transfo anisotrope (espace imreca -> espace imref)
**  \param imref : utile pour creer un champ d'une transfo anisotrope (espace imreca -> espace imref)
**  \note : si imref et imreca == NULL : la champ est genere tel quel a partir de param.
**  \retval : 1
*******************************************************************************/
int rigidz_to_field_champ_3d(int nb_param,  const double *param, field3d* points_choisis, field3d *champ,grphic3d *imref, grphic3d *imreca)
{
  int i,j,k,wdth,hght,dpth;
  double x,y,z;
  double tetax,tetay,tetaz,tx,ty,tz,xc,yc,zc,sx,sy,sz;
  double cxx,sxx,cyy,syy,czz,szz;
  double a11,a12,a13,a21,a22,a23,a31,a32,a33;
  double tpx11,tpx21,tpx31,tpy12,tpy22,tpy32;
  vector3d ***data;
  transf3d* transfOriginal   = NULL;
  transf3d* transfAnisotrope = NULL;
  vector3d ***pts_choisis;
  
 if (points_choisis==NULL)
 {
  return rigid_to_field_3d(nb_param,  param, champ, imref, imreca);
 }

  transfOriginal = cr_transf3d(0,0,0,NULL);
  transfOriginal->param = CALLOC(15,double);
  transfOriginal->typetrans = RIGIDZOOM3D;
  transfOriginal->nb_param = nb_param;
  for (i=0;i<nb_param;i++)
    transfOriginal->param[i] = param[i];


  if (imref && imreca)
    transfAnisotrope = cr_transf_anisotrope_3d(transfOriginal,imref,imreca);

  wdth=champ->width;hght=champ->height;dpth=champ->depth;
  data=champ->raw;
  pts_choisis=points_choisis->raw;

  if (transfAnisotrope && transfAnisotrope->typetrans==AFFINE3D)
  {
    a11=transfAnisotrope->param[0];a12=transfAnisotrope->param[1];a13=transfAnisotrope->param[2];
    a21=transfAnisotrope->param[3];a22=transfAnisotrope->param[4];a23=transfAnisotrope->param[5];
    a31=transfAnisotrope->param[6];a32=transfAnisotrope->param[7];a33=transfAnisotrope->param[8];
    tx=transfAnisotrope->param[9];ty=transfAnisotrope->param[10];tz=transfAnisotrope->param[11];
    xc=transfAnisotrope->param[12];yc=transfAnisotrope->param[13];zc=transfAnisotrope->param[14];
  }
  else
  {
    tetax=param[0];tetay=param[1];tetaz=param[2];
    sx=param[3];sy=param[4];sz=param[5];
    tx=param[6];ty=param[7];tz=param[8];
    xc=param[9];yc=param[10];zc=param[11];

    cxx=cos(tetax);sxx=sin(tetax);cyy=cos(tetay);syy=sin(tetay);czz=cos(tetaz);szz=sin(tetaz);
    a11=sx*cyy*czz;a12=-sy*cyy*szz;a13=-sz*syy;
    a21=sx*(cxx*szz-sxx*syy*czz);
    a22=sy*(cxx*czz+sxx*syy*szz);
    a23=-sz*sxx*cyy;
    a31=sx*(sxx*szz+cxx*syy*czz);
    a32=sy*(sxx*czz-cxx*syy*szz);
    a33=sz*cxx*cyy;
  }

  a11=a11-1.0;a22=a22-1.0;a33=a33-1.0;
  for (i=0;i<wdth;i++)
   {
    for (j=0;j<hght;j++)
    {
     for (k=0;k<dpth;k++)
      {
       x=(double)pts_choisis[i][j][k].x-xc;tpx11=a11*x;tpx21=a21*x;tpx31=a31*x;
       y=(double)pts_choisis[i][j][k].y-yc;tpy12=a12*y;tpy22=a22*y;tpy32=a32*y;
       z=(double)pts_choisis[i][j][k].z-zc;
       data[i][j][k].x=(float)(tpx11+tpy12+a13*z+tx);
       data[i][j][k].y=(float)(tpx21+tpy22+a23*z+ty);
       data[i][j][k].z=(float)(tpx31+tpy32+a33*z+tz);
      }
    }
   }

  return(1);
}
/*! @} */

/*! \ingroup TransfoFonction @{ */
/*******************************************************************************
**     affine_to_field_champ_3d(nb_param,param,points_choisis,champ,champ,imref,imreca)
*/
/*!     remplit un champ dense 3D pour une transformation affine
**  applique autour du centre (xc,yc) et avec une translation (tx,ty)
**      Les parametres dans le tableau param, sont ordonnes ainsi
**  tetax tetay tetaz tx ty tz xc yc zc
**  \param nb_param : inutilise
**  \param points_choisis  : points transformes
**  \param param  : tableau de param
**  \param champ : le champ (E/S)
**  \param imref : utile pour creer un champ d'une transfo anisotrope (espace imreca -> espace imref)
**  \param imref : utile pour creer un champ d'une transfo anisotrope (espace imreca -> espace imref)
**  \note : si imref et imreca == NULL : la champ est genere tel quel a partir de param.
**  \retval : 1
*******************************************************************************/
int affine_to_field_champ_3d(int nb_param,  const double *param, field3d* points_choisis, field3d *champ, grphic3d *imref, grphic3d *imreca)
{
  int i,j,k,wdth,hght,dpth;

  double x,y,z;
  double /*tetax,tetay,tetaz,*/tx,ty,tz,xc,yc,zc/*,sx,sy,sz*/;
  double a11,a12,a13,a21,a22,a23,a31,a32,a33;
  double tpx11,tpx21,tpx31,tpy12,tpy22,tpy32;
  vector3d ***data;
  transf3d* transfOriginal   = NULL;
  transf3d* transfAnisotrope = NULL;
  vector3d ***pts_choisis;
  
 if (points_choisis==NULL)
 {
  return rigid_to_field_3d(nb_param,  param, champ, imref, imreca);
 }

  transfOriginal = cr_transf3d(0,0,0,NULL);
  transfOriginal->param = CALLOC(15,double);
  transfOriginal->typetrans =AFFINE3D;
  transfOriginal->nb_param = nb_param;
  for (i=0;i<nb_param;i++)
    transfOriginal->param[i] = param[i];


  if (imref && imreca)
    transfAnisotrope = cr_transf_anisotrope_3d(transfOriginal,imref,imreca);

  wdth=champ->width;hght=champ->height;dpth=champ->depth;
  data=champ->raw;
  pts_choisis=points_choisis->raw;
  
  if (transfAnisotrope && transfAnisotrope->typetrans==AFFINE3D)
  {
    a11=transfAnisotrope->param[0];a12=transfAnisotrope->param[1];a13=transfAnisotrope->param[2];
    a21=transfAnisotrope->param[3];a22=transfAnisotrope->param[4];a23=transfAnisotrope->param[5];
    a31=transfAnisotrope->param[6];a32=transfAnisotrope->param[7];a33=transfAnisotrope->param[8];
    tx=transfAnisotrope->param[9];ty=transfAnisotrope->param[10];tz=transfAnisotrope->param[11];
    xc=transfAnisotrope->param[12];yc=transfAnisotrope->param[13];zc=transfAnisotrope->param[14];
  }
  else
  {
    a11=param[0];a12=param[1];a13=param[2];
    a21=param[3];a22=param[4];a23=param[5];
    a31=param[6];a32=param[7];a33=param[8];
    tx=param[9];ty=param[10];tz=param[11];
    xc=param[12];yc=param[13];zc=param[14];
  }

  a11=a11-1.0;a22=a22-1.0;a33=a33-1.0;
  for (i=0;i<wdth;i++)
   {
    for (j=0;j<hght;j++)
    {
     for (k=0;k<dpth;k++)
      {
       x=(double)pts_choisis[i][j][k].x-xc;tpx11=a11*x;tpx21=a21*x;tpx31=a31*x;
       y=(double)pts_choisis[i][j][k].y-yc;tpy12=a12*y;tpy22=a22*y;tpy32=a32*y;
       z=(double)pts_choisis[i][j][k].z-zc;
       data[i][j][k].x=(float)(tpx11+tpy12+a13*z+tx);
       data[i][j][k].y=(float)(tpx21+tpy22+a23*z+ty);
       data[i][j][k].z=(float)(tpx31+tpy32+a33*z+tz);
      }
    }
   }

  return(1);
}
/*! @} */

int rigid_global_zoom_to_field_champ_3d(int nb_param,  const double *param, field3d* points_choisis, field3d *champ, grphic3d *imref, grphic3d *imreca)
{
 double *param_zoom;
 int nb_param_zoom=12;
 int i, res;

 //on met tout comme il faut en zoom anisotrope
 param_zoom=CALLOC(nb_param_zoom,double);
 //3 premiers param=angles d'euler
 for (i=0;i<3;i++) param_zoom[i]=param[i];
 //troisieme=scale global
 for (i=3;i<6;i++) param_zoom[i]=param[3];
 //ensuite translations et centre
 for (i=6;i<12;i++) param_zoom[i]=param[i-2];
 //on utilise ce qui est deja fait
 res=rigidz_to_field_champ_3d(nb_param_zoom, param_zoom, points_choisis, champ, imref, imreca);

 FREE(param_zoom);

 return res;
}

int affine_decouple_to_field_champ_3d(int nb_param,  const double *param, field3d* points_choisis, field3d *champ,grphic3d *imref, grphic3d *imreca)
{
 double matrice_affine[15];
 int i, res=0;

 //on utilise ce qui existe deja
 for (i=0;i<nb_param;i++) matrice_affine[i]=param[i];
 affine_decouple_to_affine_3d(matrice_affine);
 res=affine_to_field_champ_3d(nb_param, matrice_affine, points_choisis, champ, imref, imreca);

 return res;
}


/*******************************************************************************
**     imx_calc_image_resultat_3d_p(imref, imreca, imres, champ, the_transf, nb_param, param)
*/
/*!     calcule l'image resultant d'une transformation calculee a partir d'un recalage
**  \param imref: l'image de reference
**  \param imreca: l'image recalee
**  \param imres: l'image resultat
**  \param nb_param : nombre de parametre de la transformation
**  \param param  : tableau de parametres de la transformation a appliquer
**  \retval : 0 en cas de succes, 1 en cas d'echec
*******************************************************************************/
int imx_calc_image_resultat_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, field3d *champ, transf3d *transfres)
{
 grphic3d *imtmp=NULL;
 if (!imres) return 1;

 transf_to_field_3d_noalloc(transfres, champ,  imref,  imreca);
 if (imreca==imres)
 {
  imtmp=cr_grphic3d(imref);
  if (!imtmp) { fprintf(stderr, "erreur d'allocation memoire dans imx_calc_image_resultat_3d_p\n"); return 1; }
  imx_copie_param_3d_p(imref,imtmp);
  inter_qsinc3_3d(imreca,champ,imtmp);
  imx_copie_3d_p(imtmp, imres);
  free_grphic3d(imtmp);
 }
 else
 {
  imx_copie_param_3d_p(imref,imres);
  inter_qsinc3_3d(imreca,champ,imres);
 }
 //inter_Bspline_3d_3(imreca,champ,imres);
 imx_inimaxminpixel_3d_p(imres);
 return 0;

}

int imx_calc_transfo_sav_image_res_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, transf3d *transfres, int save_type)
{
 transf3d *tmp_transf=NULL;
 field3d *champ1=NULL;
 int err=0;
 
 champ1=cr_field3d(imref->width,imref->height,imref->depth);
 if (!champ1) { fprintf (stderr,"erreur allocation memoire dans  imx_calc_transfo_sav_image_res_3d_p\n"); return 1; }

 err=imx_calc_image_resultat_3d_p(imref, imreca, imres, champ1, transfres);
 if (err) goto end_func;
 
 if (save_type==0)
 {
  tmp_transf=field_to_transf_3d(champ1,imref,imreca);
  copie_transf_3d(tmp_transf, transfres);
  free_transf3d(tmp_transf);
 }

end_func:
 
 free_field3d(champ1);

 return err;
}
 
/*******************************************************************************
**     imx_aff_temps_calcul(tdebut, tfin)
*/
/*!     affiche le temps de calcul d'un recalage
**  \param tdebut, tfin: les temps de debut et fin du timer
**  \retval : 0 en cas de succes
*******************************************************************************/
int imx_aff_temps_calcul(int tdebut, int tfin)
{
 int h,m,s,ttotal;
 ttotal=tfin-tdebut;
 h=ttotal/3600;
 m=(ttotal-h*3600)/60;
 s=(ttotal-h*3600-m*60);
 aff_log("TEMPS TOTAL: %dh %dm %ds \n",h,m,s);

 return 0;
}

/*******************************************************************************
**     imx_query_transformation_fct(transfo)
*/
/*!     determine la fonction de transformation correspondant a un type de transformation donne
**  \param transfo: le type de transfo voulu
**  \retval : la fonction permettant le passage des parametres au champs
*******************************************************************************/
transf_func_t imx_query_transformation_fct(TRANSFO3D transfo)
{
 transf_func_t transformation_fct=NULL;
 if (!transfo_type_est_parametree(transfo))
 {
  fprintf (stderr, "type de transformation inadapte dans imx_query_transformation_fct");
  return NULL;
 }
 switch (transfo)
 {
  case RIGID3D: transformation_fct=rigid_to_field_3d; break;
  case RIGIDZOOM3D: transformation_fct=rigidz_to_field_3d; break;
  case AFFINE3D: transformation_fct=affine_to_field_3d; break;
  case BSPLINE3D: transformation_fct=base_to_field_3d; break;
  case RIGIDGLOBALZOOM3D: transformation_fct=rigid_global_zoom_to_field_3d; break;
  case AFFINEDECOUPLE3D: transformation_fct=affine_decouple_to_field_3d; break;
  default: break;
 }
 return transformation_fct;
}

/*******************************************************************************
**     transfo_type_est_parametree(transfo)
*/
/*!     determine si le type de transformation correspond a une transformation
**      utilisant des parametres ou un champ
**  \param transfo: le type de transfo voulu
**  \retval : vrai si le type de transformation utilise des parametres
*******************************************************************************/
bool transfo_type_est_parametree(TRANSFO3D transfo)
{
 if ((transfo==RIGID3D)||
    (transfo==RIGIDZOOM3D)||
    (transfo==AFFINE3D)||
    (transfo==BSPLINE3D)||
    (transfo==RIGIDGLOBALZOOM3D)||
    (transfo==AFFINEDECOUPLE3D))
 {
  return TRUE;
 }
 else return FALSE;
}

/*******************************************************************************
**     imx_copie_dim_param_to_transf3d_p(imsrc, transfo)
*/
/*!     copie les infos de dimensions et de resolution d'une image vers une transformation
**  \param imsrc: l'image source
**  \param transfo: la transfo a modifier
**  \retval : 0 en cas de succes
*******************************************************************************/
int imx_copie_dim_param_to_transf3d_p(grphic3d *imsrc, transf3d *transfo)
{
 transfo->width=imsrc->width;
 transfo->height=imsrc->height;
 transfo->depth=imsrc->depth;
 transfo->dx=imsrc->dx;
 transfo->dy=imsrc->dy;
 transfo->dz=imsrc->dz;

 return 0;
}

/*******************************************************************************
**     imx_matching_simple_3d_p(imref, imreca, imres, dist, inter, minimisation, transfres, initrf, FieldOfResearch, contrainte, champ_ini, champ_fin)
*/
/*!     recalage de base sans pretraitement et qui minimise selon le meme type de transformation que inittrf
**  \param imref: l'image de reference
**  \param imreca: l'image recalee
**  \param imres: l'image resultat
**  \param reste : voir les parametres passes aux differentes methodes de matching
**  \retval : l'energie minimale trouvee
*******************************************************************************/
double imx_matching_simple_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, dist_func_t dist, InterpolationFct inter, min_func_t minimisation, transf3d *transfres, transf3d *inittrf, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision, contrainte_func_t contrainte, field3d *champ_ini, field3d *champ_fin)
{
 double min_param[15],max_param[15],prec_param[15];
 double *param=NULL;
 int nb_param;
 TRANSFO3D transfoType;
 double energieMini=0.0;
 transf_func_t transf=NULL;
 ERREUR_RECALAGE err=NO_ERR;
 ptr_minimisation_param min_parameters=CALLOC(1, minimisation_param);

 //toutes les images passees en parametres doivent etre non nulles
 if ((!imref)||(!imreca)||(!imres))
 {
  fprintf(stderr,"pb dans les images passees en parametre de imx_matching_simple_3d_p");
 }

 //on recherche un transformation du meme type que celui de inittrf
 //donc il doit etre non nul et doit etre defini par un ensemble de parametres
 if (!inittrf)
 {
  fprintf (stderr, "tranfo initiale nulle dans imx_matching_simple_3d_p");
  return energieMini;
 }
 transfoType=inittrf->typetrans;
 if (!transfo_type_est_parametree(transfoType))
 {
  fprintf (stderr, "type de tranfo initiale non autorisee dans imx_matching_simple_3d_p");
  return energieMini;
 }
 nb_param=inittrf->nb_param;
 transf=imx_query_transformation_fct(transfoType);

 //creation de la transformation temporaire
 copie_transf_3d(inittrf, transfres);

 //positionnement des parametres de depart
 param=transfres->param;

 //initialisation du recalage
 init_bornes_matching(transfoType, min_param, max_param, prec_param, param, FieldOfResearch, matchPrecision, imref, imreca);

if (_mtch_3d_hybride_param)
     init_bornes_matching_hybride(transfoType, min_param, max_param, prec_param, param);


 //minimisation

//-----------------------------
 min_parameters->imref=imref;
 min_parameters->imreca=imreca;
 min_parameters->imres=imres;
 min_parameters->champ_ini=champ_ini;
 min_parameters->champ_fin=champ_fin;
 min_parameters->inter=inter;
 min_parameters->dist=dist;
 min_parameters->contrainte=contrainte;
 min_parameters->min_param=min_param;
 min_parameters->max_param=max_param;
 min_parameters->prec_param=prec_param;
 min_parameters->err=&err;
 min_parameters->transfres=transfres;
 min_parameters->dist_param=cr_distParam_from_minParam(min_parameters, minimisation);

// energieMini=minimisation(imref,imreca,imres,champ_ini,champ_fin,transf,inter,dist,contrainte,nb_param,min_param,max_param,prec_param,param,donnees_VP,&err);
 energieMini=minimisation(min_parameters);
//----------------------------------

 //liberation des variables temporaires de minimisation

 free_dist_param(min_parameters->dist_param);
 FREE(min_parameters);

 return energieMini;
}

/*******************************************************************************
**     imx_matching_multires_simple_3d_p(imref, imreca, dist, inter, minimisation, transfres, initrf, contrainte, champ_ini, start_resol, end_resol)
*/
/*!     recalage en multiresolution de base sans pretraitement et qui minimise selon le meme type de transformation que inittrf
**  \param imref: l'image de reference
**  \param imreca: l'image recalee
**  \param start_resol, end_resol: resolutions max et min du recalage
**  \param reste : voir les parametres passes aux differentes methodes de matching
**  \retval : l'energie minimale trouvee
*******************************************************************************/
double imx_matching_multires_simple_3d_p(grphic3d *imref, grphic3d *imreca, dist_func_t dist, InterpolationFct inter, min_func_t minimisation, transf3d *transfres, transf3d *inittrf, eMatchPrecision matchPrecision, contrainte_func_t contrainte, field3d *champ_ini, int start_resol,int end_resol)
{
 int i;
 double energieMin=0.0;
 grphic3d *imtref=NULL, *imtreca=NULL, *imtres=NULL;
 float nResol;
 eResearchInterv FieldOfResearch;
 char logString[64];
 field3d *champ1;

 //toutes les images passees en parametres doivent etre non nulles
 if ((!imref)||(!imreca))
 {
  fprintf(stderr,"pb dans les images passees en parametre de imx_matching_simple_3d_p");
 }

 //on recherche un transformation du meme type que celui de inittrf
 //donc il doit etre non nul et doit etre defini par un ensemble de parametres
 if (!inittrf)
 {
  fprintf (stderr, "tranfo initiale nulle dans imx_matching_simple_3d_p");
  return energieMin;
 }

 if (end_resol>start_resol)
 {
  PUT_WARN("ERROR in imx_matching_multires_simple_3d_p\n"); return 0;
 }

 //creations des images sous echantillonnees
 imtreca=cr_grphic3d(imreca);
 imtref = cr_grphic3d(imref);

 //recalage en multiresolution sans changement de degre de liberte
 for (i=start_resol;i>=end_resol;i--)
 {
    // calcul des images sous echantillonees
    nResol = (float)pow(2,i);
    //nResol = (float)(i+1);

    sprintf(logString,"\n\n GRILLE : %d\n",(int)floor(nResol));
    aff_log(logString);

    //espace de recherche diminuant avec la resolution
    switch (i)
    {
      case 0 : FieldOfResearch = SMALL; break;
      case 1 : FieldOfResearch = MEDIUM; break;
      case 2 : FieldOfResearch = LARGE; break;
      default: FieldOfResearch = LARGE; break;
    }

    pyramide_gaussienne_recalage_3d_p(imref, imreca, imtref, imtreca, nResol);

    imtres=cr_grphic3d(imtref);

    //mis a jour des parametres aux nouvelles dimensions de l'image ref
    imx_copie_dim_param_to_transf3d_p(imtref, inittrf);

    //recalage
    {
     champ1=cr_field3d(imtref->width,imtref->height,imtref->depth);

     energieMin=imx_matching_simple_3d_p(imtref, imtreca, imtres, dist, inter, minimisation, transfres, inittrf, FieldOfResearch, matchPrecision, contrainte, champ_ini, champ1);

     imx_aff_param(transfres->typetrans, transfres->param);
     free_field3d(champ1);
    }

    if (imtres) free_grphic3d(imtres); imtres=NULL;

    //on se sert de la transformation resultante pour initialiser la transfo de depart
    copie_transf_3d(transfres,inittrf);
 }

 //liberation memoire des images temporaires
 free_grphic3d(imtref);
 free_grphic3d(imtreca);

 return energieMin;
}
/*******************************************************************************
**     determination_minima_grille_angulaire(energies, indices, nb_angles)
*/
/*!    determine les minima locaux du tableau d'energies calculees
**
**     utilisation:
**     \param       energies: energies calculees
**     \param       indices: tableau d'indices ou seront places les indices correpondant a des minima locaux d'energie
**     \param       nb_angles : nb de points ou l'energie   a ete calculee
**     \retval      nb de minima locaux trouves
*******************************************************************************/
int determination_minima_grille_angulaire(double *energies, int *indices, int nb_angles)
{
 int i,j,k;
 int i_param;
 int nb_minima_locaux=0;
 int i_voisin_bas, i_voisin_haut, j_voisin_bas, j_voisin_haut, k_voisin_bas, k_voisin_haut;

 for (i=0;i<nb_angles;i++)
 {
  for (j=0;j<nb_angles;j++)
  {
   for (k=0;k<nb_angles;k++)
   {
    i_param=i*nb_angles*nb_angles+j*nb_angles+k;
    i_voisin_bas=MAXI(i-1,0); j_voisin_bas=MAXI(j-1,0); k_voisin_bas=MAXI(k-1,0);
    i_voisin_haut=MINI(i+1,nb_angles-1); j_voisin_haut=MINI(j+1,nb_angles-1); k_voisin_haut=MINI(k+1,nb_angles-1);
    if ((energies[i_param]<=energies[i_voisin_bas*nb_angles*nb_angles+j*nb_angles+k])&&
        (energies[i_param]<=energies[i*nb_angles*nb_angles+j_voisin_bas*nb_angles+k])&&
        (energies[i_param]<=energies[i*nb_angles*nb_angles+j*nb_angles+k_voisin_bas])&&
        (energies[i_param]<=energies[i_voisin_haut*nb_angles*nb_angles+j*nb_angles+k])&&
        (energies[i_param]<=energies[i*nb_angles*nb_angles+j_voisin_haut*nb_angles+k])&&
        (energies[i_param]<=energies[i*nb_angles*nb_angles+j*nb_angles+k_voisin_haut]))
     {
      indices[nb_minima_locaux]=i_param;
      nb_minima_locaux++;
     }
   }
  }
 }

 return nb_minima_locaux;
}

/*******************************************************************************
**     evaluer_parametres_grille_fine(param_calc_fin, param_calc_coarse, nb_angles_fin, nb_angles_grossier, transfoType)
*/
/*!    evalue les parametres correspondant a la grille angulaire fine a partir de ceux
**     correspodant a la grille angulaire grossiere par interpolation trilineaire
**
**     utilisation:
**     \param       param_calc_fin: parametres pour la grille fine
**     \param       param_calc_coarse: parametres pour la grille grossiere
**     \param       nb_angles_fin, nb_angles_grossier : nb de divisions angulaires en grilles fine et grossiere
**     \param       transfoType: type de transfo utilisee dans le recalage courant
**     \retval      0 si tout s'est bien passe
*******************************************************************************/
int evaluer_parametres_grille_fine(double **param_calc_fin, double **param_calc_coarse, int nb_angles_fin, int nb_angles_grossier, TRANSFO3D transfoType)
{
 int i,j,k,l;
 int pos_param_rot=0, pos_param_scale, pos_param_trans, pos_param_centre;
 double i_rel, j_rel, k_rel;
 int i_coarse, j_coarse, k_coarse;
 double i_coarsed, j_coarsed, k_coarsed;
 int i_coarse_min, i_coarse_max, j_coarse_min, j_coarse_max, k_coarse_min, k_coarse_max;
 int i_param, i_param_coarse;
 double pas_angulaire_fin[3]={0.0}, min_angulaire_fin[3]={0.0};
 double *scale_tab=NULL;
 int nb_pts_coarse=nb_angles_grossier*nb_angles_grossier*nb_angles_grossier;
 double scale_median=1.0;

 if (transfoType==RIGIDGLOBALZOOM3D) { pos_param_scale=3; pos_param_trans=4; pos_param_centre=7; }
 else if (transfoType==RIGID3D) { pos_param_scale=0; pos_param_trans=3; pos_param_centre=6; }
 else
 {
  fprintf (stderr,"mauvais type de transformation initiale dans evaluer_parametres_grille_fine\n");
  return 1;
 }

 if (nb_angles_fin<=nb_angles_grossier)
 {
  fprintf (stderr,"le grille fine est plus grossiere que la grille grossiere!\n");
  return 2;
 }

 //determination du scale median
 if (transfoType==RIGIDGLOBALZOOM3D)
 {
  scale_tab=CALLOC(nb_pts_coarse,double);
  for (i=0;i<nb_pts_coarse;i++) scale_tab[i]=param_calc_coarse[i][pos_param_scale];
  scale_median=quick_select(nb_pts_coarse/2,scale_tab,nb_pts_coarse);
  aff_log ("scale median:%f\n",scale_median);
  FREE(scale_tab);
 }

 //calcul des parametres angulaires
 for (i=0;i<3;i++)
 {
  pas_angulaire_fin[i]=(param_calc_coarse[nb_angles_grossier*nb_angles_grossier*nb_angles_grossier-1][pos_param_rot+i]-param_calc_coarse[0][pos_param_rot+i])/(double)(nb_angles_fin-1);
  min_angulaire_fin[i]=param_calc_coarse[0][pos_param_rot+i];
 }

 //mise en place des parametres
 for (i=0;i<nb_angles_fin;i++)
 {
  for (j=0;j<nb_angles_fin;j++)
  {
   for (k=0;k<nb_angles_fin;k++)
   {
    i_param=i*nb_angles_fin*nb_angles_fin+j*nb_angles_fin+k;
    param_calc_fin[i_param][pos_param_rot]=i*pas_angulaire_fin[0]+min_angulaire_fin[0];
    param_calc_fin[i_param][pos_param_rot+1]=j*pas_angulaire_fin[1]+min_angulaire_fin[1];
    param_calc_fin[i_param][pos_param_rot+2]=k*pas_angulaire_fin[2]+min_angulaire_fin[2];

    //parametres en translation
    //il sont interpoles a partir des resultats trouves pour la grille grossiere par interpolation trilineaire
    i_rel=modf((double)(i*(nb_angles_grossier-1))/(double)(nb_angles_fin-1),&i_coarsed);
    j_rel=modf((double)(j*(nb_angles_grossier-1))/(double)(nb_angles_fin-1),&j_coarsed);
    k_rel=modf((double)(k*(nb_angles_grossier-1))/(double)(nb_angles_fin-1),&k_coarsed);

    i_coarse=(int)floor(i_coarsed); j_coarse=(int)floor(j_coarsed); k_coarse=(int)floor(k_coarsed);

    i_param_coarse=i_coarse*nb_angles_grossier*nb_angles_grossier+j_coarse*nb_angles_grossier+k_coarse;
    i_coarse_min=i_coarse; j_coarse_min=j_coarse; k_coarse_min=k_coarse;
    i_coarse_max=MINI(i_coarse_min+1,nb_angles_grossier-1); j_coarse_max=MINI(j_coarse_min+1,nb_angles_grossier-1); k_coarse_max=MINI(k_coarse_min+1,nb_angles_grossier-1);

    //determination des parametres interpoles en translation
    for (l=0;l<3;l++)
    {
     param_calc_fin[i_param][pos_param_trans+l]=(1.0-i_rel)*(1.0-j_rel)*(1.0-k_rel)*param_calc_coarse[i_coarse_min*nb_angles_grossier*nb_angles_grossier+j_coarse_min*nb_angles_grossier+k_coarse_min][pos_param_trans+l] +
                                                i_rel*(1.0-j_rel)*(1.0-k_rel)*param_calc_coarse[i_coarse_max*nb_angles_grossier*nb_angles_grossier+j_coarse_min*nb_angles_grossier+k_coarse_min][pos_param_trans+l] +
                                                (1.0-i_rel)*j_rel*(1.0-k_rel)*param_calc_coarse[i_coarse_min*nb_angles_grossier*nb_angles_grossier+j_coarse_max*nb_angles_grossier+k_coarse_min][pos_param_trans+l] +
                                                (1.0-i_rel)*(1.0-j_rel)*k_rel*param_calc_coarse[i_coarse_min*nb_angles_grossier*nb_angles_grossier+j_coarse_min*nb_angles_grossier+k_coarse_max][pos_param_trans+l] +
                                                i_rel*j_rel*(1.0-k_rel)*param_calc_coarse[i_coarse_max*nb_angles_grossier*nb_angles_grossier+j_coarse_max*nb_angles_grossier+k_coarse_min][pos_param_trans+l] +
                                                i_rel*(1.0-j_rel)*k_rel*param_calc_coarse[i_coarse_max*nb_angles_grossier*nb_angles_grossier+j_coarse_min*nb_angles_grossier+k_coarse_max][pos_param_trans+l] +
                                                (1.0-i_rel)*j_rel*k_rel*param_calc_coarse[i_coarse_min*nb_angles_grossier*nb_angles_grossier+j_coarse_max*nb_angles_grossier+k_coarse_max][pos_param_trans+l] +
                                                i_rel*j_rel*k_rel*param_calc_coarse[i_coarse_max*nb_angles_grossier*nb_angles_grossier+j_coarse_max*nb_angles_grossier+k_coarse_max][pos_param_trans+l];
    }
    //positionnement du centre
    for (l=0;l<3;l++) param_calc_fin[i_param][pos_param_centre+l]=param_calc_coarse[0][pos_param_centre+l];

    //position du zoom global
    if (transfoType==RIGIDGLOBALZOOM3D)
    {
     param_calc_fin[i_param][pos_param_scale]=scale_median;
    }

   }
  }
 }

 return 0;
}

/*******************************************************************************
**     initialiser_parametres_grossier(param_calc_coarse, nb_angles_grossier, max_param, min_param, param, transfoType)
*/
/*!    initialise les parametres correspondant a la grille grossiere a partir des parametres en entree du recalage
**
**     utilisation:
**     \param       param_calc_coarse: parametres pour la grille grossiere
**     \param       nb_angles_grossier : nb de divisions angulaires en grille grossiere
**     \param       max_param, min_param, param: parametres en entree du recalage
**     \param       transfoType: type de transfo utilisee dans le recalage courant
**     \retval      0 si tout s'est bien passe
*******************************************************************************/
int initialiser_parametres_grossier(double **param_calc_coarse, int nb_angles_grossier, double *max_param, double *min_param, double *param, TRANSFO3D transfoType)
{
 int i,j,k,l;
 int pos_param_rot=0, pos_param_scale, pos_param_trans, pos_param_centre;
 int i_param;
 double pas_angulaire_grossier[3]={0.0};

 if (transfoType==RIGIDGLOBALZOOM3D) { pos_param_scale=3; pos_param_trans=4; pos_param_centre=7; }
 else if (transfoType==RIGID3D) { pos_param_scale=0; pos_param_trans=3; pos_param_centre=6; }
 else
 {
  fprintf (stderr,"mauvais type de transformation initiale dans initialiser_parametres_grossier\n");
  return 1;
 }

 for (i=0;i<3;i++)
 {
  pas_angulaire_grossier[i]=(max_param[pos_param_rot+i]-min_param[pos_param_rot+i])/(double)(nb_angles_grossier-1);
 }

 for (i=0;i<nb_angles_grossier;i++)
 {
  for (j=0;j<nb_angles_grossier;j++)
  {
   for (k=0;k<nb_angles_grossier;k++)
   {
    i_param=i*nb_angles_grossier*nb_angles_grossier+j*nb_angles_grossier+k;

    //positionnement des parametres en translation
    for (l=0;l<3;l++) param_calc_coarse[i_param][pos_param_trans+l]=param[pos_param_trans+l];

    //positionnement des parametres de rotation
    param_calc_coarse[i_param][pos_param_rot]=(double)i*pas_angulaire_grossier[0]+min_param[pos_param_rot];
    param_calc_coarse[i_param][pos_param_rot+1]=(double)j*pas_angulaire_grossier[1]+min_param[pos_param_rot+1];
    param_calc_coarse[i_param][pos_param_rot+2]=(double)k*pas_angulaire_grossier[2]+min_param[pos_param_rot+2];
   
    //positionnement du zoom global
    if (transfoType==RIGIDGLOBALZOOM3D)
    {
     param_calc_coarse[i_param][pos_param_scale]=1.0;
    }

    //positionnement du centre
    for (l=0;l<3;l++) param_calc_coarse[i_param][pos_param_centre+l]=param[pos_param_centre+l];
   }
  }
 }

 return 0;
}

/*******************************************************************************
**     min_multistart_multires_8mm_3d(imref, imreca, imres, dist, inter, minimisation, transfres, inittrf, FieldOfResearch, contrainte, champ_ini, champ_fin, nb_minima, minima)
*/
/*!    minimisations multiples realisees a la resolution la plus basse:
**     - il y a d'abord des minimisations a rotations fixes a partir d'une grille grossiere en rotation
**     - puis des minimisations rigides a partir d'une grille fine en rotation
**     - puis un renvoi des nb_minima meilleures solutions trouvees
**
**     utilisation:
**     \param       imref, imreca, imres, dist, inter, minimisation, transfres, inittrf, FieldOfResearch, contrainte, champ_ini, champ_fin : param de matching
**     \param       nb_minima : nombre d'optima trouves
**     \param       minima : tableau ou seront stockes les vecteurs parametres correspondant a ces optima
*******************************************************************************/
double imx_matching_multistart_multires_8mm_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, dist_func_t dist, InterpolationFct inter, min_func_t minimisation, transf3d *transfres, transf3d *inittrf, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision, contrainte_func_t contrainte, field3d *champ_ini, field3d *champ_fin, int *nb_minima, double **minima)
{
 int i,j;
 int nb_angles_grossier=3, nb_angles_fin=7;
 double energie_min=DBL_MAX;
 ERREUR_RECALAGE err=NO_ERR;

 //variables utilisees en grille grossiere
 double **param_calc_coarse=NULL;
 double *energie_coarse=NULL;
 int nb_pts_coarse=nb_angles_grossier*nb_angles_grossier*nb_angles_grossier;

 //variable utilisees en grille fine
 double **param_calc_fin=NULL;
 double *energie_fin=NULL;
 int nb_pts_fin=nb_angles_fin*nb_angles_fin*nb_angles_fin;

 int nb_minima_locaux=0;
 int *indices_minima_locaux=NULL;
 int pos_param_rot=0, pos_param_scale, pos_param_trans, pos_param_centre;

 TRANSFO3D transfoType;
 transf_func_t transf=NULL;
 int nb_param=0;
 double max_param[15]={0.0}, min_param[15]={0.0}, prec_param[15]={0.0};

 double *param=NULL;
 ptr_minimisation_param min_parameters=CALLOC(1, minimisation_param);

 //toutes les images passees en parametres doivent etre non nulles
 if ((!imref)||(!imreca)||(!imres))
 {
  RAISE_ERR_TEXT(&err, ERREUR_ALLOC_MEMOIRE, "pb dans les images passees en parametre de imx_matching_multistart_multires_8mm_3d_p\n");
 }

 //inittrf doit etre non nul et de type rigide ou rigide+zoom global
 if (!inittrf)
 {
  RAISE_ERR_TEXT(&err, ERREUR_ALLOC_MEMOIRE, "tranfo initiale nulle dans imx_matching_multistart_multires_8mm_3d_p\n");
 }
 transfoType=inittrf->typetrans;

 if ((transfoType!=RIGID3D)&&(transfoType!=RIGIDGLOBALZOOM3D))
 {
  RAISE_ERR_TEXT(&err, TRANSFO_INADAPTEE, "mauvais type de transformation initiale dans imx_matching_multistart_multires_8mm_3d_p\n");
 }

 //--------------------- INITS ---------------------------//

 if (transfoType==RIGIDGLOBALZOOM3D)
 { pos_param_scale=3; pos_param_trans=4; pos_param_centre=7; }
 else
 { pos_param_scale=0; pos_param_trans=3; pos_param_centre=6; }

 nb_param=inittrf->nb_param;
// for (i=0;i<nb_param;i++) param[i]=inittrf->param[i];
 copie_transf_3d(inittrf, transfres);
 param=transfres->param;

 //allocations memoire
 param_calc_coarse=CALLOC(nb_pts_coarse, double*);
 energie_coarse=CALLOC(nb_pts_coarse, double);
 param_calc_fin=CALLOC(nb_pts_fin, double*);
 energie_fin=CALLOC(nb_pts_fin, double);

 if ((!param_calc_coarse)||(!energie_coarse)||(!param_calc_fin)||(!energie_fin))
 {
  RAISE_ERR_TEXT(&err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans imx_matching_multistart_multires_8mm_3d_p\n");
 }

 for (i=0;i<nb_pts_coarse;i++)
 {
  param_calc_coarse[i]=CALLOC(nb_param, double);
  if (!param_calc_coarse[i])
  {
   RAISE_ERR_TEXT(&err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans imx_matching_multistart_multires_8mm_3d_p\n");
  }
 }
 for (i=0;i<nb_pts_fin;i++)
 {
  param_calc_fin[i]=CALLOC(nb_param, double);
  if (!param_calc_fin[i])
  {
   RAISE_ERR_TEXT(&err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans imx_matching_multistart_multires_8mm_3d_p\n");
  }
 }

 indices_minima_locaux=CALLOC(nb_pts_fin, int);

 //initialisation des recalages
 transf=imx_query_transformation_fct(transfoType);

//-----------------------------
 min_parameters->imref=imref;
 min_parameters->imreca=imreca;
 min_parameters->imres=imres;
 min_parameters->champ_ini=champ_ini;
 min_parameters->champ_fin=champ_fin;
 min_parameters->inter=inter;
 min_parameters->dist=dist;
 min_parameters->contrainte=contrainte;
 min_parameters->min_param=min_param;
 min_parameters->max_param=max_param;
 min_parameters->prec_param=prec_param;
 min_parameters->err=&err;
 min_parameters->transfres=transfres;
 min_parameters->dist_param=cr_distParam_from_minParam(min_parameters, minimisation);

 //--------------------- MINIMISATIONS GRILLE GROSSIERE ---------------------------//

 init_bornes_matching(transfoType, min_param, max_param, prec_param, param, FieldOfResearch, matchPrecision, imref, imreca);
 
 initialiser_parametres_grossier(param_calc_coarse, nb_angles_grossier, max_param, min_param, param, transfoType);

 //toutes les minimisations se font a parametres de rotation fixes
 init_bornes_matching(transfoType, min_param, max_param, prec_param, param, SMALL, matchPrecision, imref, imreca);
 prec_param[pos_param_rot]=0.0; prec_param[pos_param_rot+1]=0.0; prec_param[pos_param_rot+2]=0.0;

 //minimisations sur chaque points de la grille en ne faisant varier que la translation
 //et (eventuellement) le zoom global
 for (i=0;i<nb_pts_coarse;i++)
 {
  aff_log ("minimisation no: %d sur %d\n", i+1, nb_pts_coarse);

 min_parameters->transfres->param=param_calc_coarse[i];
//  energie_coarse[i]=minimisation(imref, imreca, imres, champ_ini, champ_fin, transf, inter, dist, contrainte, nb_param, min_param, max_param, prec_param, param_calc_coarse[i], donnees_VP, &err);
 energie_coarse[i]=minimisation(min_parameters);

 min_parameters->transfres->param=param;

  if (erreur_critique(&err)) goto end_func;
 }

 //--------------------- DETERMINATION MINIMA LOCAUX ---------------------------//

 //evaluation des parametres de la grille fine par interpolation lineaire
 evaluer_parametres_grille_fine(param_calc_fin, param_calc_coarse, nb_angles_fin, nb_angles_grossier, transfoType);

 //calcul des energies sur la grille fine a partir des parametres interpoles
 for (i=0;i<nb_pts_fin;i++)
 {
  min_parameters->transfres->param=param_calc_fin[i];
//  energie_fin[i]=fonc_transfo_3d(imref, imreca, imres, transf, inter, dist, contrainte, NULL, champ_fin, nb_param, param_calc_fin[i], donnees_VP, points_calcules, NULL,&err);
  energie_fin[i]=energie_transfo_3d(min_parameters);
  min_parameters->transfres->param=param;
  if (erreur_critique(&err)) goto end_func;
 }

 //determination des minima locaux a partir des energies trouvees pour la grille fine calculees ci-dessus
 nb_minima_locaux=determination_minima_grille_angulaire(energie_fin, indices_minima_locaux, nb_angles_fin);

 //--------------------- MINIMISATIONS SUR MINIMA LOCAUX ---------------------------//

 //on va minimiser sur tous les parametres a partir de tous les minima locaux trouves
 for (i=0;i<nb_minima_locaux;i++)
 {
  aff_log ("recalage minima locaux no: %d sur %d\n", i+1,nb_minima_locaux);
  init_bornes_matching(transfoType, min_param, max_param, prec_param, param_calc_fin[indices_minima_locaux[i]], SMALL, matchPrecision, imref, imreca);

 min_parameters->transfres->param=param_calc_fin[indices_minima_locaux[i]];
//  energie_fin[indices_minima_locaux[i]]=minimisation(imref, imreca, imres, NULL, champ_fin, transf, inter, dist, NULL, nb_param, min_param, max_param, prec_param, param_calc_fin[indices_minima_locaux[i]], donnees_VP, &err);
 energie_fin[indices_minima_locaux[i]]=minimisation(min_parameters);

 min_parameters->transfres->param=param;

  if (erreur_critique(&err)) goto end_func;
 }
 aff_log ("\n");

 //--------------------- RELEASE ---------------------------//

 //on selectionne les meilleurs points
 if (nb_minima_locaux<(*nb_minima)) (*nb_minima)=nb_minima_locaux;

 tri_double_index(energie_fin, nb_minima_locaux, indices_minima_locaux);

 for (i=0;i<(*nb_minima);i++)
 {
  for (j=0;j<nb_param;j++)
  {
   minima[i][j]=param_calc_fin[indices_minima_locaux[i]][j];
  }
 }

 //on met l'optimum trouve dans transfres
 for (i=0;i<nb_param;i++) transfres->param[i]=param_calc_fin[indices_minima_locaux[0]][i];
 energie_min=energie_fin[indices_minima_locaux[0]];

end_func:

 //liberations memoire
 free_dist_param(min_parameters->dist_param);
 FREE(min_parameters);

 if (indices_minima_locaux) FREE(indices_minima_locaux);
 if (energie_fin) FREE(energie_fin);
 for (i=0;i<nb_pts_fin;i++) { if (param_calc_fin[i]) FREE(param_calc_fin[i]); }
 if (param_calc_fin) FREE(param_calc_fin);
 if (energie_coarse) FREE(energie_coarse);
 for (i=0;i<nb_pts_coarse;i++) { if (param_calc_coarse[i]) FREE(param_calc_coarse[i]); }
 if (param_calc_coarse) FREE(param_calc_coarse);

 aff_log ("Energie min: %f\n",energie_min);
 aff_log("imx_matching_multistart_multires_8mm_3d_p\n");
 return energie_min;
}

/*******************************************************************************
**     min_multistart_multires_4mm_3d(imref, imreca, imres, dist, inter, minimisation, transfres, inittrf, FieldOfResearch, contrainte, champ_ini, champ_fin, nb_minima, minima)
*/
/*!    minimisation multistart a partir d'optima locaux obtenus a la resolution inferieure passes en entree
**     et de ces optimas auxquels on a rajoute des perturbations en rotation
**
**     utilisation:
**     \param       imref, imreca, imres, dist, inter, minimisation, transfres, inittrf, FieldOfResearch, contrainte, champ_ini, champ_fin : param de matching
**     \param       nb_minima : nombre d'optima
**     \param       minima : tableau des vecteurs parametres correspondant a ces optima
*******************************************************************************/
double imx_matching_multistart_multires_4mm_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, dist_func_t dist, InterpolationFct inter, min_func_t minimisation, transf3d *transfres, transf3d *inittrf, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision, contrainte_func_t contrainte, field3d *champ_ini, field3d *champ_fin, int nb_minima, double **minima)
{
 int i,j,k;
 double **param_opt=NULL, *energies=NULL;
 int *indexes=NULL;
 int pos_param_rot=0, pos_param_scale=3;
 int nb_pts=0, nb_perturbations=0;
 double energie_min=DBL_MAX;
 transf_func_t transf=NULL;
 double max_param[15]={0.0}, min_param[15]={0.0}, prec_param[15]={0.0};
 TRANSFO3D transfoType;
 int nb_param=0;
 ERREUR_RECALAGE err=NO_ERR;

 double *param=NULL;

 ptr_minimisation_param min_parameters=CALLOC(1, minimisation_param);

 //perturbations en rotation
 //correspond a la precision esperee a la resolution inferieure
 //l'image passee en entree est supposee isotrope
 double drot=imref->dx*2.0/(2.0*BRAIN_RADIUS);
 double dscale=imref->dx*2.0/(2.0*BRAIN_RADIUS);

 //toutes les images passees en parametres doivent etre non nulles
 if ((!imref)||(!imreca)||(!imres))
 {
  RAISE_ERR_TEXT(&err, ERREUR_ALLOC_MEMOIRE, "pb dans les images passees en parametre de imx_matching_multistart_multires_4mm_3d_p\n");
 }

 //inittrf doit etre non nul et de type rigide ou rigide+zoom global
 if (!inittrf)
 {
  RAISE_ERR_TEXT(&err, ERREUR_ALLOC_MEMOIRE, "tranfo initiale nulle dans imx_matching_multistart_multires_4mm_3d_p\n");
 }
 transfoType=inittrf->typetrans;

 if ((transfoType!=RIGID3D)&&(transfoType!=RIGIDGLOBALZOOM3D))
 {
  RAISE_ERR_TEXT(&err, TRANSFO_INADAPTEE, "mauvais type de transformation initiale dans imx_matching_multistart_multires_4mm_3d_p\n");
 }

 if (transfoType==RIGIDGLOBALZOOM3D) nb_perturbations=10;
 else nb_perturbations=7;

 nb_pts=nb_perturbations*nb_minima;

 transf=imx_query_transformation_fct(transfoType);

 nb_param=inittrf->nb_param;
// for (i=0;i<nb_param;i++) param[i]=inittrf->param[i];

 //allocations memoire
 param_opt=CALLOC(nb_pts,double *);
 indexes=CALLOC(nb_pts,int);
 energies=CALLOC(nb_pts,double);
 if ((!param_opt)||(!indexes)||(!energies))
 {
  RAISE_ERR_TEXT(&err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans imx_matching_multistart_multires_4mm_3d_p\n");
 }

 for (i=0;i<nb_pts;i++)
 {
  param_opt[i]=CALLOC(nb_param, double);
  if (!param_opt[i])
  {
   RAISE_ERR_TEXT(&err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans imx_matching_multistart_multires_4mm_3d_p\n");
  }
 }

 //---------------- CREATION DES PARAMETRES PERTURBES ---------------//

 //on va creer des un ensembles de vecteurs de parametres correspondant aux optima locaux passes en entree
 //et a ces optima perturbes de drot
 for (i=0;i<nb_pts;i++) indexes[i]=i;
 //mise en place des params
 for (i=0;i<nb_minima;i++)
 {
  for (j=0;j<nb_perturbations;j++)
  {
   for (k=0;k<nb_param;k++) param_opt[i*nb_perturbations+j][k]=minima[i][k];
  }
 }
 //rajout des perturbations
 for (i=0;i<nb_minima;i++)
 {
  param_opt[i*nb_perturbations+1][pos_param_rot]=minima[i][pos_param_rot]-drot;
  param_opt[i*nb_perturbations+2][pos_param_rot]=minima[i][pos_param_rot]+drot;
  param_opt[i*nb_perturbations+3][pos_param_rot+1]=minima[i][pos_param_rot+1]-drot;
  param_opt[i*nb_perturbations+4][pos_param_rot+1]=minima[i][pos_param_rot+1]+drot;
  param_opt[i*nb_perturbations+5][pos_param_rot+2]=minima[i][pos_param_rot+2]-drot;
  param_opt[i*nb_perturbations+6][pos_param_rot+2]=minima[i][pos_param_rot+2]+drot;
  if (transfoType==RIGIDGLOBALZOOM3D)
  {
   param_opt[i*nb_perturbations+7][pos_param_scale]=minima[i][pos_param_scale];
   param_opt[i*nb_perturbations+8][pos_param_scale]=minima[i][pos_param_scale]-dscale;
   param_opt[i*nb_perturbations+9][pos_param_scale]=minima[i][pos_param_scale]+dscale;
  }
 }

 //---------------- MINIMISATIONS SUR LES PARAMETRES PERTURBES ---------------//
  copie_transf_3d(inittrf, transfres);
  param=transfres->param;

//-----------------------------
 min_parameters->imref=imref;
 min_parameters->imreca=imreca;
 min_parameters->imres=imres;
 min_parameters->champ_ini=champ_ini;
 min_parameters->champ_fin=champ_fin;
 min_parameters->inter=inter;
 min_parameters->dist=dist;
 min_parameters->contrainte=contrainte;
 min_parameters->min_param=min_param;
 min_parameters->max_param=max_param;
 min_parameters->prec_param=prec_param;
 min_parameters->err=&err;
 min_parameters->transfres=transfres;
 min_parameters->dist_param=cr_distParam_from_minParam(min_parameters, minimisation);

 //minimisations a partir de chaque vecteur parametre passe en entree
 for (i=0;i<nb_pts;i++)
 {
   init_bornes_matching(transfoType, min_param, max_param, prec_param, param_opt[i], SMALL, matchPrecision, imref, imreca);

   min_parameters->transfres->param=param_opt[i];
  // energies[i]=minimisation(imref, imreca, imres, NULL, champ_fin, transf, inter, dist, NULL, nb_param, min_param, max_param, prec_param, param_opt[i], donnees_VP, &err);
   energies[i]=minimisation(min_parameters);
 }
 min_parameters->transfres->param=param;

 //tri des energies trouvees
 tri_double_index(energies, nb_pts, indexes);

 //on met l'optimum trouve dans transfres
 for (i=0;i<nb_param;i++) transfres->param[i]=param_opt[indexes[0]][i];
 energie_min=energies[indexes[0]];

 printf ("Energie min: %f\n",energie_min);

end_func:

 //liberations memoire

 free_dist_param(min_parameters->dist_param);
 FREE(min_parameters);

 if (energies) FREE(energies);
 if (indexes) FREE(indexes);
 for (i=0;i<nb_pts;i++) { if (param_opt[i]) FREE(param_opt[i]); }
 if (param_opt) FREE(param_opt);

 return energie_min;
}

eMatchPrecision imx_query_matching_precision()
{
 eMatchPrecision matchPrecision;
 int matchPrecInt=0;
 char *query[5];

    query[0]="Normal";
    query[1]="Precis (1e-1)";
    query[2]="Tres precis (1e-2)";
    query[3]="Precision max (1e-4)";
    query[4]=NULL;

 matchPrecInt=GETV_QCM("Precision requise pour le matching",(char **)query);

 switch (matchPrecInt)
 {
  case 0: matchPrecision=NORMAL; break;
  case 1: matchPrecision=PRECIS; break;
  case 2: matchPrecision=TRES_PRECIS; break;
  case 3: matchPrecision=PRECISION_MAX; break;
  default: matchPrecision=NORMAL; break;
 }

 return matchPrecision;
}

ptr_distance_param cr_distParam_from_minParam(ptr_minimisation_param min_parameters, min_func_t minimisation)
{
 ptr_distance_param dist_param=NULL;
 VP_histos *donnees_VP=NULL;
 double ***weights=NULL;
 //recuperation des donnees utiles dans les parametres de minimisation
 dist_func_t distance=min_parameters->dist;
 grphic3d *imref=min_parameters->imref;
 grphic3d *imreca=min_parameters->imreca;
 grphic3d *imres=min_parameters->imres;
 field3d *champ_fin=min_parameters->champ_fin; 
 ERREUR_RECALAGE *err=min_parameters->err;
 int i, j ,k;

 dist_param=CALLOC(1, distance_param);

  //  Initialisation de _sigma_robust dans le cas de l'ancien recalage robuste
  if ((distance==erreur_quad_robust_geman_3d)
    ||(distance==erreur_woods_robust_3d))
  {
   calc_robust_stat_3d(imref,imreca,0,0);
   aff_log("\n");
   aff_log("Valeur de rejet du robuste _sigma_robust=%f \n",_sigma_robust);
  }

  // initialisation pour l'interpolation volume partiel et les distances basees sur l'histogramme conjoint
  if (dist_utilisant_VP(distance))
  {

   dist_param->imref=imref;
   dist_param->imreca=imreca;

   if ((minimisation==min_desc_grad_3d)
     ||(minimisation==min_grad_conj_3d)
     ||(minimisation==min_desc_grad_mod_3d)
     ||(minimisation==min_quasi_newton_mod_3d))
   {
    donnees_VP = cr_VP_histos(imreca,imref,1);
   } else {
    donnees_VP = cr_VP_histos(imreca,imref,0);
   }
    printf("VP: %d\n", donnees_VP->max1);

   dist_param->donnees_VP=donnees_VP;

  }
  else
  {

   dist_param->imref=imref;
   dist_param->imreca=imres;
   dist_param->donnees_VP=NULL;

  }

  dist_param->champ=champ_fin;
  dist_param->err=err;

  //inits specifiques aux distances basees sur un echantillonage stochastique
  if (distance==erreur_IM_VP_stochastique_3d)
  {
   dist_param->points_calcules=choix_points_aleatoires(imref->width, imref->height, imref->depth, err);
  }
  else
  {
   dist_param->points_calcules=NULL;
  }

 //creation du tableau des poids robustes
  if ((distance==erreur_quad_robust2_3d)
    ||(distance==erreur_woods_robust2_3d)
    ||(distance==erreur_CR_robust_3d))
  {
   weights=CALLOC(imref->width,double**);
//   if (!weights) {
//    RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans min_simplex_robuste_3d\n");
//   }
   for (i=0;i<(int)imref->width;i++)
   {
    weights[i]=CALLOC(imref->height,double*);
//     if (!weights[i]) {
//      RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans min_simplex_robuste_3d\n");
//     }
    for (j=0;j<(int)imref->height;j++)
    {
     weights[i][j]=CALLOC(imref->depth,double);
//      if (!weights[i][j]) {
//       RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans min_simplex_robuste_3d\n");
//      }
    }
   }

  //init
  for (i=0;i<(int)imref->width;i++)
  {
   for (j=0;j<(int)imref->height;j++)
   {
    for (k=0;k<(int)imref->depth;k++)
    {
     weights[i][j][k]=1.0;
    }
   }
  }

  dist_param->weights=weights;
 }
 else
 {
  dist_param->weights=NULL;
 }

 dist_param->nbclRef=imx_calculer_nombre_de_classes_3d_p(imref, &(dist_param->maxRef));
 dist_param->nbclReca=imx_calculer_nombre_de_classes_3d_p(imreca, &(dist_param->maxReca));

 return dist_param;
}

int free_dist_param(ptr_distance_param dist_param)
{
 grphic3d *imref=dist_param->imref;
 double ***weights=dist_param->weights;
 field3d *points_calcules=dist_param->points_calcules;
 int i, j;

 if (dist_param->donnees_VP)
 {
  free_VP_histos(&(dist_param->donnees_VP));
 }

 if (weights)
 {
  for (i=0;i<(int)imref->width;i++)
  {
    for (j=0;j<(int)imref->height;j++) { if (weights[i][j]) FREE(weights[i][j]); }
    if (weights[i]) FREE(weights[i]);
  }
  FREE(weights);
 }

 if (points_calcules)
 {
  free_field3d(points_calcules);
 }

 FREE(dist_param);

 return 0;
}



/*******************************************************************************
**        matching_rigide_ICP_3d                                                 
**                                                                        
**                             
*******************************************************************************/
/*! 
     recalage rigide de deux images 3D binaire (sans zoom)        
  
*/
void matching_rigide_ICP_3d()
{
 int im_ref,im_reca,im_res,i;
 char *nomfichres=NULL;
 char *nomfichtrf=NULL;
 int dist_type,inter_type,min_type,save_type;
 eResearchInterv FieldOfResearch;
 char *quest[3];

 //demande des images
 im_reca=GET_PLACE3D(TEXT0030);
 im_ref=GET_PLACE3D(TEXT0031);
 im_res=GET_PLACE3D(TEXT0006);

 //distance
 for(i=0;i<3;i++)
    quest[i]=CALLOC(80,char);

  strcpy(quest[0],"ICP ");
  strcpy(quest[1],"ICP sym");
  strcpy(quest[2],"\0");

 dist_type=GETV_QCM("Distance",(char **)quest);
 
 for(i=0;i<3;i++)
    free(quest[i]);

 //interpolation lin�aire
 inter_type=1; //imx_query_interpolation_type_3D(dist_type);
 //minimisation
 min_type=imx_query_minimisation_type_3D(dist_type);
 //fichier d'enregistrement du resultat
 nomfichres=quest_save_result_3d(&save_type);
 //intervalle de recherche
 FieldOfResearch = MEDIUM;

//matchPrecision=NORMAL

 //matching
 imx_matching_rigide_ICP_3d(im_ref,im_reca,im_res,dist_type,inter_type,min_type,save_type,nomfichres,nomfichtrf,FieldOfResearch,PRECIS);

 if (nomfichres) free(nomfichres);
 if (nomfichtrf) free(nomfichtrf);

 show_picture_3d(im_res);


}

/*******************************************************************************
**        imx_matching_rigide_ICP_3d                                             
*/                                                                       
/*!       recalage rigide sans zoom de 2 images 3D                                   
** 
**  \param im_ref : numero de l'image de reference, 
**  \param im_reca  : numero de l'image a recaler,
**  \param im_res : numero de l'image resultat
**  les autres parametres sont expliques dans ::imx_matching_rigide_3d_p    
**  \retval 1
**  
*******************************************************************************/
int imx_matching_rigide_ICP_3d(int im_ref, int im_reca, int im_res, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres, char *nomfichtrf,eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision)
{
  grphic3d *imref,*imreca,*imres;
  transf3d *inittrf=NULL;
  transf3d *transfres=NULL;
 
  imref=ptr_img_3d(im_ref);
  imreca=ptr_img_3d(im_reca);
  imres=ptr_img_3d(im_res);
   
  if (nomfichtrf !=NULL)
     inittrf=load_transf_3d(nomfichtrf);
  transfres = cr_transf3d_p(imref,NULL); 

  init_log(NULL);
  aff_log("\n RECALAGE RIGIDE SANS ZOOM ");
  aff_log("\n");
  aff_log("recalage de %s sur %s \n",(imreca->patient).name,(imref->patient).name);

  imx_inimaxminpixel_3d_p(imref);
  imx_inimaxminpixel_3d_p(imreca);

  imx_matching_rigide_ICP_3d_p(imref,imreca,imres,dist_type,inter_type,min_type,save_type,transfres,inittrf,FieldOfResearch,matchPrecision);

  /*enregistrement du champ resultat*/
  if (nomfichres!=NULL)
  { 
    save_transf_3d(transfres,nomfichres);
  }

  end_log();

  if (inittrf)
    free_transf3d(inittrf);

  free_transf3d(transfres);
 
  return(1);
}

/*! \ingroup Recalage  @{ */
/*******************************************************************************
**        imx_matching_rigide_ICP_3d_p                                            
*/
/*!                                                                        
**       \brief recalage rigide entre 2 images 3D                                    
**
**  \param imref : images de reference, 
**  \param imreca  : image a recaler,
**  \param imres : images resultat (E/S)
**  \param dist_type : numero de la fonction de distance 
**                          \li 0 : Quadratic
**                          \li 1 : Quadratic Robust
**                          \li 2 : woods
**                          \li 3 : woods robust
**                          \li 4 : IM Old
**                          \li 5 : IM
**                          \li 6 : IM NORM STUDHOLME
**                          \li 7 : IM NORM MAES1
**                          \li 8 : IM NORM MAES2
**                          \li 9 : Entropie conjointe
**                          \li 10 : Correlation ratio
**  \param inter_type : numero de la fonction d'interpolation des images
**                          \li 0 : nearest
**                          \li 1 : linear
**                          \li 2 : sin card
**                          \li 3 : quick sin card2
**                          \li 4 : quick sin card3 
**                          \li 5 : bspline 2
**                          \li 6 : bspline 3
**                          \li 7 : bspline 4
**                          \li 8 : bspline 5
**  \param min_type : numero de la fonction de minimisation utilisee aux plus faibles
    echelles de la transformation 
**                          \li 0 : ICM
**                          \li 1 : simplex
**                          \li 2 : Descente gradient
**                          \li 3 : Gradient conjugue
**                          \li 4 : Descente gradient modifiee
**                          \li 5 : Quasi Newton modifiee
**  \param save_type : Type de sauvegarde
**                          \li 0 : champ
**                          \li 1 : parametres
**  \param nomfichres :  nom du fichier de sauvegarde des transformations (*.trf)
**  \param nomfichtrf : initialisation du champ de depart a la valeur d'un fichier (*.trf)\n
**                      si NULL -> alignement centre de gravite
**  \param FieldOfResearch :  voir init_bornes_matching
**  \param matchPrecision  :  voir init_bornes_matching
**  \retval 1
**
*******************************************************************************/
int imx_matching_rigide_ICP_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int dist_type, int inter_type, int min_type, int save_type, transf3d *transfres, transf3d *inittrf,eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision)
{
 grphic3d *imtref,*imtreca,*imtres,*imtmp;
 field3d *champ1;
 InterpolationFct interpol;
 dist_func_t  distance;
 min_func_t minimisation;
 double  mini;
 int    tdebut,tfin;
 transf3d *preProcessingRefTrf=NULL, *preProcessingRecaTrf=NULL, *transfoini=NULL;
 
#ifndef WIN32
 setvbuf(stdout, (char *)NULL, _IONBF, 0); //setbuf(stdout,NULL);
#endif
 tdebut=time(NULL);
 
 //---------mise en place des parametres du recalage-------/
  
  if (dist_type == 1)
    distance=erreur_ICP_sym_3d;
    else
    distance=erreur_ICP_3d;
  
    interpol=imx_choose_interpolation_fct(inter_type);
     
  minimisation=imx_choose_minimisation_fct(min_type,dist_type);  

 // ------------------PREPROCESS----------------------
 {
 // imx_matching_preprocessing_3d_p(imref, imreca, &imtref, &imtreca, &imtres, &preProcessingRefTrf, &preProcessingRecaTrf, distance);
 }

  imtreca=cr_grphic3d(imreca);imx_copie_3d_p(imreca,imtreca);
  imtref=cr_grphic3d(imref);imx_copie_3d_p(imref,imtref);
  imtres=cr_grphic3d(imref);imx_copie_3d_p(imref,imtres);
    imtreca->mask=cr_grphic3d(imreca);imx_copie_3d_p(imreca,imtreca->mask);
    imtref->mask=cr_grphic3d(imref);imx_copie_3d_p(imref,imtref->mask); 
    imtmp=cr_grphic3d(imreca);imx_copie_3d_p(imreca,imtmp);
    imtres->mask=cr_grphic3d(imtres);imx_copie_3d_p(imtres,imtres->mask); 
    

  /* calcul du mask de imref avec l'image compl�mentaire */
        imx_inv_3d_p(imtref,imtref->mask);
        imx_inv_3d_p(imreca,imtmp);
        
 // ------------------INITIALISATION PAR TRANSFO DE DEPART---------------------- 
 {
  //on veut une transformation de depart compatible avec le recalage courant
  if ((inittrf)&&(inittrf->typetrans==RIGID3D))
  {
   transfoini=imx_trouver_transf_initiale_matching_3d_p(imtref, imtreca, inittrf);
  }
  else transfoini=imx_trouver_transf_initiale_matching_3d_p(imtref, imtreca, NULL);
 }

/* calcul des cartes de distance */
    imx_chamfer_distance_3d_p(imreca,imtreca);
    imx_chamfer_distance_3d_p(imtmp,imtreca->mask);

    imx_inimaxminpixel_3d_p(imtreca);
    imx_inimaxminpixel_3d_p(imtreca->mask);

  imx_mul_coe_3d_p(imtreca,10.0,imtreca);
    imx_mul_coe_3d_p(imtreca->mask,10.0,imtreca->mask);

  imx_add_coe_3d_p(imtreca,1.0,imtreca);
  imx_add_coe_3d_p(imtreca->mask,1.0,imtreca->mask);

    
    imtres->max_pixel=imtreca->max_pixel;
    imtres->mask->max_pixel=imtreca->mask->max_pixel;
    imtres->rcoeff=imtreca->rcoeff;
    imtres->mask->rcoeff=imtreca->mask->rcoeff;




 //------------------RECALAGE RIGIDE----------------------------------------/
 aff_log("\n TRANSFO RIGIDE \n");
 {
  champ1=cr_field3d(imtref->width,imtref->height,imtref->depth);

  mini=imx_matching_simple_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, FieldOfResearch, matchPrecision, NULL, NULL, champ1);

  free_field3d(champ1);
 }

 //on remet le tout dans le repere des images originales
 imx_copie_dim_param_to_transf3d_p(imref, transfres);
// if (preProcessingRecaTrf) comb_transf_3d(preProcessingRecaTrf,transfres,transfres,0);
 //if (preProcessingRefTrf) comb_transf_3d(transfres,preProcessingRefTrf,transfres,0);

 //if (preProcessingRecaTrf) free_transf3d(preProcessingRecaTrf);
 //if (preProcessingRefTrf) free_transf3d(preProcessingRefTrf);
 if (imtres) free_grphic3d(imtres);
 if (imtref!=imref) free_grphic3d(imtref);
 if (imtreca!=imreca) free_grphic3d(imtreca);
 if (imtmp)  free_grphic3d(imtmp);

 // ------------------CALCUL DE L'IMAGE A AFFICHER ET DE LA TRANSFO A SAUVEGARDER----------------------
 imx_aff_param(transfres->typetrans, transfres->param);
 imx_calc_transfo_sav_image_res_3d_p(imref, imreca, imres, transfres, save_type);

 //liberation memoire des variables dynamiques 
 if (transfoini) free_transf3d(transfoini);
 
 tfin=time(NULL);
 //affichage du temps de calcul total 
 {
  imx_aff_temps_calcul(tdebut, tfin);
 }

 return(1);  
}



/*******************************************************************************
**        matching_rigide_zoom_ICP_3d                                                 
**                                                                        
**       recalage rigide + zoom de deux images binaires 3D (avec zoom)                                
*******************************************************************************/
void matching_rigide_zoom_ICP_3d()
{
int im_ref,im_reca,im_res,i;
 char *nomfichres=NULL;
 char *nomfichtrf=NULL,*quest[3];
 int dist_type,inter_type,min_type,save_type;
 eResearchInterv FieldOfResearch;

 //demande des images
 im_reca=GET_PLACE3D(TEXT0030);
 im_ref=GET_PLACE3D(TEXT0031);
 im_res=GET_PLACE3D(TEXT0006);

 //distance
 for(i=0;i<3;i++)
    quest[i]=CALLOC(80,char);

  strcpy(quest[0],"ICP ");
  strcpy(quest[1],"ICP sym");
  strcpy(quest[2],"\0");

 dist_type=GETV_QCM("Distance",(char **)quest);
 
 for(i=0;i<3;i++)
    free(quest[i]);
        
    //interpolation lin�aire
 inter_type=1; //imx_query_interpolation_type_3D(dist_type);
 //minimisation
 min_type=imx_query_minimisation_type_3D(dist_type);
 //fichier d'enregistrement du resultat
 nomfichres=quest_save_result_3d(&save_type);
 //intervalle de recherche
 FieldOfResearch = MEDIUM;

//matchPrecision=NORMAL

 //matching
 imx_matching_rigide_zoom_ICP_3d(im_ref,im_reca,im_res,dist_type,inter_type,min_type,save_type,nomfichres,nomfichtrf,FieldOfResearch,PRECIS);

 if (nomfichres) free(nomfichres);
 if (nomfichtrf) free(nomfichtrf);

 show_picture_3d(im_res);
}


/*******************************************************************************
**        imx_matching_rigide_zoom_ICP_3d                                             
**                                                                       
**       recalage rigide avec zoom de 2 images 3D                                   
*******************************************************************************/
int imx_matching_rigide_zoom_ICP_3d(int im_ref, int im_reca, int im_res, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres, char *nomfichtrf, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision)
{
 grphic3d *imref,*imreca,*imres;
 transf3d *transfres=NULL;
 transf3d *inittrf=NULL;

 imref=ptr_img_3d(im_ref);
 imreca=ptr_img_3d(im_reca);
 imres=ptr_img_3d(im_res);

 if (nomfichtrf !=NULL)
     inittrf=load_transf_3d(nomfichtrf);
 transfres = cr_transf3d_p(imref, NULL);

 imx_inimaxminpixel_3d_p(imref);
 imx_inimaxminpixel_3d_p(imreca);

 init_log(NULL);
 aff_log("\n RECALAGE RIGIDE AVEC ZOOM ");
 aff_log("\n");
 aff_log("recalage de %s sur %s \n",(imreca->patient).name,(imref->patient).name);

 imx_matching_rigide_zoom_ICP_3d_p(imref,imreca,imres,dist_type,inter_type,min_type,save_type,transfres,inittrf,FieldOfResearch,matchPrecision);

 /*enregistrement du champ resultat*/
 if (nomfichres!=NULL)
    save_transf_3d(transfres,nomfichres);

 end_log();

 free_transf3d(transfres);
 
  return(1);
}


/*! \ingroup     Recalage  @{ */
/*******************************************************************************
**        imx_matching_rigide_zoom_ICP_3d_p
**                                         
**       \param :pour les parametres voir imx_matching_rigide_3d_p
**                                                              
**       recalage rigide avec zoom de 2 images 3D                                    
*******************************************************************************/
int imx_matching_rigide_zoom_ICP_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int dist_type, int inter_type, int min_type, int save_type, transf3d *transfres, transf3d *inittrf,eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision)
{
 int wdth,hght,dpth;
 grphic3d *imtref,*imtreca,*imtres,*imtmp;
 field3d *champ1;
 InterpolationFct interpol;
 dist_func_t distance;
 min_func_t minimisation;
 double  mini;
 int    tdebut,tfin;
 transf3d *preProcessingRefTrf=NULL, *preProcessingRecaTrf=NULL, *transfoini=NULL;
 
 wdth=imref->width;hght=imref->height;dpth=imref->depth;
#ifndef WIN32
 setvbuf(stdout, (char *)NULL, _IONBF, 0); //setbuf(stdout,NULL);
#endif
 tdebut=time(NULL);

 //---------mise en place des parametres du recalage-------/

  if (dist_type == 1)
    distance=erreur_ICP_sym_3d;
    else
    distance=erreur_ICP_3d;
    
  interpol=imx_choose_interpolation_fct(inter_type);

  minimisation=imx_choose_minimisation_fct(min_type,dist_type);

 // ------------------PREPROCESS----------------------
 {
// imx_matching_preprocessing_3d_p(imref, imreca, &imtref, &imtreca, &imtres, &preProcessingRefTrf, &preProcessingRecaTrf, distance);
 }

  imtreca=cr_grphic3d(imreca);imx_copie_3d_p(imreca,imtreca);
  imtref=cr_grphic3d(imref);imx_copie_3d_p(imref,imtref);
  imtres=cr_grphic3d(imref);imx_copie_3d_p(imref,imtres);
    imtreca->mask=cr_grphic3d(imreca);imx_copie_3d_p(imreca,imtreca->mask);
    imtref->mask=cr_grphic3d(imref);imx_copie_3d_p(imref,imtref->mask); 
    imtmp=cr_grphic3d(imreca);imx_copie_3d_p(imreca,imtmp);
    imtres->mask=cr_grphic3d(imtres);imx_copie_3d_p(imtres,imtres->mask); 
    
    
 // ------------------INITIALISATION PAR TRANSFO DE DEPART----------------------
 {
  //on veut une transformation de depart compatible avec le recalage courant
  if ((inittrf)&&(inittrf->typetrans==RIGID3D))
  {
   transfoini=imx_trouver_transf_initiale_matching_3d_p(imtref, imtreca, inittrf);
  }
  else transfoini=imx_trouver_transf_initiale_matching_3d_p(imtref, imtreca, NULL);
 }

 // ------------------MISE A ZERO DES ELEMENTS SUR LESQUELS ON OPTIMISE PAS----------------------
    if (_mtch_3d_hybride_param)
        {int i;
            for (i=0;i<3;i++)
                if (_mtch_3d_hybride_param[i]==0)
                {transfoini->param[i]=0.0;}
            for (i=3;i<6;i++)
                if (_mtch_3d_hybride_param[i+3]==0)
                {transfoini->param[i]=0.0;}
            
            for (i=6;i<9;i++)
                if (_mtch_3d_hybride_param[i+3]==0)
                {transfoini->param[i]=_mtch_3d_hybride_param[i+6];}
            
                
        }

 /* calcul du mask de imref avec l'image compl�mentaire */
        imx_inv_3d_p(imtref,imtref->mask);
        imx_inv_3d_p(imreca,imtmp);
    
    
/* calcul des cartes de distance */
    imx_chamfer_distance_3d_p(imreca,imtreca);
    imx_chamfer_distance_3d_p(imtmp,imtreca->mask);

    imx_inimaxminpixel_3d_p(imtreca);
    imx_inimaxminpixel_3d_p(imtreca->mask);

  imx_mul_coe_3d_p(imtreca,10.0,imtreca);
    imx_mul_coe_3d_p(imtreca->mask,10.0,imtreca->mask);

  imx_add_coe_3d_p(imtreca,1.0,imtreca);
  imx_add_coe_3d_p(imtreca->mask,1.0,imtreca->mask);

    
    imtres->max_pixel=imtreca->max_pixel;
    imtres->mask->max_pixel=imtreca->mask->max_pixel;
    imtres->rcoeff=imtreca->rcoeff;
    imtres->mask->rcoeff=imtreca->mask->rcoeff;


    

 //------------------RECALAGE ----------------------------------------/
 {
 wdth=imtref->width;hght=imtref->height;dpth=imtref->depth;
 champ1=cr_field3d(wdth,hght,dpth);

 //recalage rigide
 aff_log("\n TRANSFO RIGIDE \n");
 {
  mini=imx_matching_simple_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, FieldOfResearch, matchPrecision, NULL, NULL, champ1);
 }

 //construction de la transfo initiale a injecter en debut du recalage rigide+zoom
 imx_aff_param(RIGID3D, transfres->param);
 copie_transf_3d(transfres, transfoini);
 
if (_mtch_3d_hybride_param)
    {if (_mtch_3d_hybride_param[15]!=0)
        {
        transfoini->param=REALLOC(transfoini->param,transfoini->nb_param+1,double);
        transfoini->nb_param=rigid_to_rigid_global_zoom_3d(transfoini->param);
        transfoini->typetrans=RIGIDGLOBALZOOM3D; 
        aff_log("\n TRANSFO RIGIDE+ZOOM GLOBAL \n");
        }
    else
        {transfoini->param=REALLOC(transfoini->param,transfoini->nb_param+3,double);
        transfoini->nb_param=rigid_to_rigidz_3d(transfoini->param);
        transfoini->typetrans=RIGIDZOOM3D;
        aff_log("\n TRANSFO RIGIDE+ZOOM \n");
        }
    }
    else
        {transfoini->param=REALLOC(transfoini->param,transfoini->nb_param+3,double);
        transfoini->nb_param=rigid_to_rigidz_3d(transfoini->param);
        transfoini->typetrans=RIGIDZOOM3D;
        aff_log("\n TRANSFO RIGIDE+ZOOM \n");
        }
        
 //recalage rigide+zoom
 {
  mini=imx_matching_simple_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, FieldOfResearch, matchPrecision, NULL, NULL, champ1);
 }

 free_field3d(champ1);
 }

 //on remet le tout dans le repere des images originales
 imx_copie_dim_param_to_transf3d_p(imref, transfres);
 //if (preProcessingRecaTrf) comb_transf_3d(preProcessingRecaTrf,transfres,transfres,0);
 //if (preProcessingRefTrf) comb_transf_3d(transfres,preProcessingRefTrf,transfres,0);

 //if (preProcessingRecaTrf) free_transf3d(preProcessingRecaTrf);
 //if (preProcessingRefTrf) free_transf3d(preProcessingRefTrf);
 if (imtres) free_grphic3d(imtres);
 if (imtref!=imref) free_grphic3d(imtref);
 if (imtreca!=imreca) free_grphic3d(imtreca);
 if (imtmp) free_grphic3d(imtmp);
 
 // ------------------CALCUL DE L'IMAGE A AFFICHER ET DE LA TRANSFO A SAUVEGARDER----------------------
 imx_aff_param(transfres->typetrans, transfres->param);
 imx_calc_transfo_sav_image_res_3d_p(imref, imreca, imres, transfres, save_type);

 //liberation memoire des variables dynamiques
 if (transfoini) free_transf3d(transfoini);

 tfin=time(NULL);
 //affichage du temps de calcul total
 {
  imx_aff_temps_calcul(tdebut, tfin);
 }
 end_log();

 return(1);  
}

/*******************************************************************************
**        matching_hybride_ICP_3d                                                 
*/
/*!                                                                        
**      \brief recalage rigide + zoom suivant un nombre restreint de parametres
**                  entre de deux images binaires 3D                                
*******************************************************************************/
void matching_hybride_ICP_3d()
{
int im_ref,im_reca,im_res,e,i;
 char *nomfichres=NULL;
 char *nomfichtrf=NULL,*quest[3];
 int dist_type,inter_type,min_type,save_type;
 eResearchInterv FieldOfResearch;
char    s[250];

_mtch_3d_hybride_param = (double *) malloc (16*sizeof(double));

 //demande des images
 im_reca=GET_PLACE3D(TEXT0030);
 im_ref=GET_PLACE3D(TEXT0031);
 im_res=GET_PLACE3D(TEXT0006);

  //distance
 for(i=0;i<3;i++)
    quest[i]=CALLOC(80,char);

  strcpy(quest[0],"ICP ");
  strcpy(quest[1],"ICP sym");
  strcpy(quest[2],"\0");

 dist_type=GETV_QCM("Distance",(char **)quest);
 
 for(i=0;i<3;i++)
    free(quest[i]);
        
 //interpolation lin�aire
 inter_type=1; //imx_query_interpolation_type_3D(dist_type);
 //minimisation
 min_type=imx_query_minimisation_type_3D(dist_type);
 //fichier d'enregistrement du resultat
 nomfichres=quest_save_result_3d(&save_type);
 //intervalle de recherche
 FieldOfResearch = MEDIUM;

// Questions sur les parametres a optimiser
sprintf(s,"Optimiser teta_x ?\n");
    _mtch_3d_hybride_param[0]=GET_EXIST(s,&e);

sprintf(s,"Optimiser teta_y ?\n");
    _mtch_3d_hybride_param[1]=GET_EXIST(s,&e);

sprintf(s,"Optimiser teta_z ?\n");
    _mtch_3d_hybride_param[2]=GET_EXIST(s,&e);

sprintf(s,"Optimiser zoom global ?\n");
    _mtch_3d_hybride_param[15]=GET_EXIST(s,&e);

if (_mtch_3d_hybride_param[15]==0)
{sprintf(s,"Optimiser zoom_x ?\n");
    _mtch_3d_hybride_param[3]=GET_EXIST(s,&e);}

if (_mtch_3d_hybride_param[15]==0)
{sprintf(s,"Optimiser zoom_y ?\n");
    _mtch_3d_hybride_param[4]=GET_EXIST(s,&e);}

if (_mtch_3d_hybride_param[15]==0)
{sprintf(s,"Optimiser zoom_z ?\n");
    _mtch_3d_hybride_param[5]=GET_EXIST(s,&e);}

sprintf(s,"Optimiser t_x ?\n");
    _mtch_3d_hybride_param[6]=GET_EXIST(s,&e);

sprintf(s,"Optimiser t_y ?\n");
    _mtch_3d_hybride_param[7]=GET_EXIST(s,&e);

sprintf(s,"Optimiser t_z ?\n");
    _mtch_3d_hybride_param[8]=GET_EXIST(s,&e);

_mtch_3d_hybride_param[9]=1;
_mtch_3d_hybride_param[10]=1;
_mtch_3d_hybride_param[11]=1;


if ((_mtch_3d_hybride_param[1]!=0)&&(_mtch_3d_hybride_param[6]==0)&&(_mtch_3d_hybride_param[9]==1))
{
_mtch_3d_hybride_param[12]=GET_DOUBLE("xc",0,&e);
_mtch_3d_hybride_param[9]=0;
}

if ((_mtch_3d_hybride_param[2]!=0)&&(_mtch_3d_hybride_param[6]==0)&&(_mtch_3d_hybride_param[9]==1))
{
_mtch_3d_hybride_param[12]=GET_DOUBLE("xc",0,&e);
_mtch_3d_hybride_param[9]=0;
}


if ((_mtch_3d_hybride_param[0]!=0)&&(_mtch_3d_hybride_param[7]==0)&&(_mtch_3d_hybride_param[10]==1))
{
_mtch_3d_hybride_param[13]=GET_DOUBLE("yc",0,&e);
_mtch_3d_hybride_param[10]=0;
}

if ((_mtch_3d_hybride_param[2]!=0)&&(_mtch_3d_hybride_param[7]==0)&&(_mtch_3d_hybride_param[10]==1))
{
_mtch_3d_hybride_param[13]=GET_DOUBLE("yc",0,&e);
_mtch_3d_hybride_param[10]=0;
}

if ((_mtch_3d_hybride_param[0]!=0)&&(_mtch_3d_hybride_param[8]==0)&&(_mtch_3d_hybride_param[11]==1))
{
_mtch_3d_hybride_param[14]=GET_DOUBLE("zc",0,&e);
_mtch_3d_hybride_param[11]=0;
}



if ((_mtch_3d_hybride_param[1]!=0)&&(_mtch_3d_hybride_param[8]==0)&&(_mtch_3d_hybride_param[11]==1))
{
_mtch_3d_hybride_param[14]=GET_DOUBLE("zc",0,&e);
_mtch_3d_hybride_param[11]=0;
}


 //matching
 imx_matching_rigide_zoom_ICP_3d(im_ref,im_reca,im_res,dist_type,inter_type,min_type,save_type,nomfichres,nomfichtrf,FieldOfResearch,PRECIS);

 if (nomfichres) free(nomfichres);
 if (nomfichtrf) free(nomfichtrf);

free(_mtch_3d_hybride_param);

 show_picture_3d(im_res);}




/*******************************************************************************
**        matching_affine_ICP_3d                                                 
*/
/*!                                                                        
**      \brief recalage affine de deux images binaires 3D                                
*******************************************************************************/
void matching_affine_ICP_3d()
{
int im_ref,im_reca,im_res,i;
 char *nomfichres=NULL;
 char *nomfichtrf=NULL,*quest[3];
 int dist_type,inter_type,min_type,save_type;
 eResearchInterv FieldOfResearch;

 //demande des images
 im_reca=GET_PLACE3D(TEXT0030);
 im_ref=GET_PLACE3D(TEXT0031);
 im_res=GET_PLACE3D(TEXT0006);

  //distance
 for(i=0;i<3;i++)
    quest[i]=CALLOC(80,char);

  strcpy(quest[0],"ICP ");
  strcpy(quest[1],"ICP sym");
  strcpy(quest[2],"\0");

 dist_type=GETV_QCM("Distance",(char **)quest);
 
 for(i=0;i<3;i++)
    free(quest[i]);
 //interpolation lin�aire
 inter_type=1; //imx_query_interpolation_type_3D(dist_type);
 //minimisation
 min_type=imx_query_minimisation_type_3D(dist_type);
 //fichier d'enregistrement du resultat
 nomfichres=quest_save_result_3d(&save_type);
 //intervalle de recherche
 FieldOfResearch = MEDIUM;

//matchPrecision=NORMAL

 //matching
 imx_matching_affine_ICP_3d(im_ref,im_reca,im_res,dist_type,inter_type,min_type,save_type,nomfichres,nomfichtrf,FieldOfResearch,PRECIS);

 if (nomfichres) free(nomfichres);
 if (nomfichtrf) free(nomfichtrf);

 show_picture_3d(im_res);}



/*******************************************************************************
**        imx_matching_affine_ICP_3d                                             
*/
/*!     \brief  recalage affine de 2 images 3D                                    
**
**      \param  imref : numero de l'image de ref
**      \param  imreca : numero de l'image a recaler
**      \param  imres : numero de l'image ou l'on stocke le resultat du recalage de imreca
    (peut etre NULL)
**      \param  dist_type : numero de la fonction de distance (0 pour distance quadratique)
                            \li 0 : QUAD
                            \li 1 : QUAD ROBUST
                            \li 2 : WOODS
                            \li 3 : WOODS ROBUST
                            \li 4 : IM
                            \li 
                            
**      \param  inter_type : numero de la fonction d'interpolation des images (1 pour lineaire)
                            \li 0 : NEAREST
                            \li 1 : LINEAR
                            \li 2 : SINC
                            \li 3 : QSINC2
                            \li 4 : QSINC3
                            \li 5 : BSPLINE2
                            \li 6 : BSPLINE3 
                            \li 7 : BSPLINE4
                            \li 8 : BSPLINE5


**      \param  min_type : numero de la fonction de minimisation utilisee aux plus faibles
    echelles de la transformation (5 pour quasi-newton)
                            \li 0 : ICM
                            \li 1 : SIMPLEX 
                            \li 2 : DESCGRAD 
                            \li 3 : GRADCONJ
                            \li 4 : DESCGRADMOD
                            \li 5 : QUASINEWTONMOD

**      \param save_type : 
**                          \li 0 : champ
**                          \li 1 : parametres
**      \param  nomfichres :  nom du fichier de sauvegarde de la transformation
**      \retval 1
*******************************************************************************/
int imx_matching_affine_ICP_3d(int im_ref, int im_reca, int im_res, int dist_type, int inter_type, int min_type, int save_type, char *nomfichres, char *nomfichtrf, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision)
{
 grphic3d *imref,*imreca,*imres;
 transf3d *inittrf=NULL;
 transf3d *transfres=NULL;
 imref=ptr_img_3d(im_ref);
 imreca=ptr_img_3d(im_reca);
 imres=ptr_img_3d(im_res);

 if (nomfichtrf !=NULL)
     inittrf=load_transf_3d(nomfichtrf);  
 transfres=cr_transf3d_p(imref,NULL);

 imx_inimaxminpixel_3d_p(imref);
 imx_inimaxminpixel_3d_p(imreca);

 init_log(NULL);
 aff_log("\n RECALAGE AFFINE ");
 aff_log("\n");
 aff_log("recalage de %s sur %s \n",(imreca->patient).name,(imref->patient).name);

 imx_matching_affine_ICP_3d_p(imref,imreca,imres,dist_type,inter_type,min_type,save_type,transfres,inittrf,FieldOfResearch,matchPrecision);

 /*enregistrement du champ resultat*/
 if (nomfichres!=NULL)
 { 
   save_transf_3d(transfres,nomfichres);
 } 

 end_log();
/*   if (inittrf)
     free_transf3d(inittrf);
*/
  free_transf3d(transfres);
  return(1);
}

/*! \ingroup     Recalage @{ */
/*******************************************************************************
**        imx_matching_affine_ICP_3d_p                                            
**                                                                        
*/
/*!     \brief  recalage affine de 2 images 3D                                    
**
**      \param  imref : image de reference
**      \param  imreca : image a recaler
**      \param  imres : image resultat (E/S)
**      \param  dist_type : numero de la fonction de distance 
                            \li 0 : QUAD
                            \li 1 : QUAD ROBUST
                            \li 2 : WOODS
                            \li 3 : WOODS ROBUST
                            \li 4 : IM
                            \li 
                            
**      \param inter_type : numero de la fonction d'interpolation des images 
                            \li 0 : NEAREST
                            \li 1 : LINEAR
                            \li 2 : SINC
                            \li 3 : QSINC2
                            \li 4 : QSINC3
                            \li 5 : BSPLINE2
                            \li 6 : BSPLINE3 
                            \li 7 : BSPLINE4
                            \li 8 : BSPLINE5

**      \param  min_type : numero de la fonction de minimisation utilisee aux plus faibles
    echelles de la transformation 
                            \li 0 : ICM
                            \li 1 : SIMPLEX 
                            \li 2 : DESCGRAD 
                            \li 3 : GRADCONJ
                            \li 4 : DESCGRADMOD
                            \li 5 : QUASINEWTONMOD

**      \param save_type : 
**                          \li 0 : champ
**                          \li 1 : parametres
**      \param  nomfichres :   nom du fichier de sauvegarde de la transfo (*.trf)
**  \param FieldOfResearch :  voir init_bornes_matching
**  \param matchPrecision  :  voir init_bornes_matching
**      \retval 1
*******************************************************************************/
int imx_matching_affine_ICP_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres, int dist_type, int inter_type, int min_type, int save_type, transf3d *transfres,transf3d *inittrf, eResearchInterv FieldOfResearch, eMatchPrecision matchPrecision)
{
 int    wdth,hght,dpth,k;
 grphic3d *imtref,*imtreca,*imtres,*imtmp;
 field3d *champ1;
 InterpolationFct interpol;
 dist_func_t distance;
 min_func_t minimisation;
 double  mini,sous;
 int    tdebut,tfin;
 transf3d *preProcessingRefTrf=NULL, *preProcessingRecaTrf=NULL, *transfoini=NULL;

 wdth=imref->width;hght=imref->height;dpth=imref->depth;
#ifndef WIN32
 setvbuf(stdout, (char *)NULL, _IONBF, 0); //setbuf(stdout,NULL);
#endif
 tdebut=time(NULL);

 //---------mise en place des parametres du recalage-------/

  if (dist_type == 1)
    distance=erreur_ICP_sym_3d;
    else
    distance=erreur_ICP_3d;

  interpol=imx_choose_interpolation_fct(inter_type);

  minimisation=imx_choose_minimisation_fct(min_type,dist_type);

 // ------------------PREPROCESS----------------------
 {
  //imx_matching_preprocessing_3d_p(imref, imreca, &imtref, &imtreca, &imtres, &preProcessingRefTrf, &preProcessingRecaTrf, distance);
 }

    imtreca=cr_grphic3d(imreca);imx_copie_3d_p(imreca,imtreca);
  imtref=cr_grphic3d(imref);imx_copie_3d_p(imref,imtref);
  imtres=cr_grphic3d(imref);imx_copie_3d_p(imref,imtres);
    imtreca->mask=cr_grphic3d(imreca);imx_copie_3d_p(imreca,imtreca->mask);
    imtref->mask=cr_grphic3d(imref);imx_copie_3d_p(imref,imtref->mask); 
    imtmp=cr_grphic3d(imreca);imx_copie_3d_p(imreca,imtmp);
    imtres->mask=cr_grphic3d(imtres);imx_copie_3d_p(imtres,imtres->mask); 
        
 // ------------------INITIALISATION PAR TRANSFO DE DEPART----------------------
 {
  //on veut une transformation de depart compatible avec le recalage courant
  if ((inittrf)&&(inittrf->typetrans==RIGID3D))
  {
   transfoini=imx_trouver_transf_initiale_matching_3d_p(imtref, imtreca, inittrf);
  }
  else transfoini=imx_trouver_transf_initiale_matching_3d_p(imtref, imtreca, NULL);
 }
for (k=0;k<6;k++) transfoini->param[k]=0.0;

    /* calcul du mask de imref avec l'image compl�mentaire */
    imx_inv_3d_p(imtref,imtref->mask);
    imx_inv_3d_p(imreca,imtmp);
    
    
    /* calcul des cartes de distance */
    imx_chamfer_distance_3d_p(imreca,imtreca);
    imx_chamfer_distance_3d_p(imtmp,imtreca->mask);

    imx_inimaxminpixel_3d_p(imtreca);
    imx_inimaxminpixel_3d_p(imtreca->mask);

  imx_mul_coe_3d_p(imtreca,10.0,imtreca);
    imx_mul_coe_3d_p(imtreca->mask,10.0,imtreca->mask);

  imx_add_coe_3d_p(imtreca,1.0,imtreca);
  imx_add_coe_3d_p(imtreca->mask,1.0,imtreca->mask);

    
    imtres->max_pixel=imtreca->max_pixel;
    imtres->mask->max_pixel=imtreca->mask->max_pixel;
    imtres->rcoeff=imtreca->rcoeff;
    imtres->mask->rcoeff=imtreca->mask->rcoeff;




 //------------------RECALAGE ----------------------------------------/
 {
 wdth=imtref->width;hght=imtref->height;dpth=imtref->depth;
 champ1=cr_field3d(wdth,hght,dpth);

 //recalage rigide
 aff_log("\n TRANSFO RIGIDE \n");
 {
  mini=imx_matching_simple_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, FieldOfResearch, matchPrecision, NULL, NULL, champ1);
 }

 //construction de la transfo initiale a injecter en debut du recalage rigide+zoom
 imx_aff_param(RIGID3D, transfres->param);
 copie_transf_3d(transfres, transfoini);
 transfoini->param=REALLOC(transfoini->param,transfoini->nb_param+3,double);
 transfoini->nb_param=rigid_to_rigidz_3d(transfoini->param);
 transfoini->typetrans=RIGIDZOOM3D;

 //recalage rigide+zoom
 aff_log("\n TRANSFO RIGIDE+ZOOM \n");
 {
  mini=imx_matching_simple_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, FieldOfResearch, matchPrecision, NULL, NULL, champ1);
 }

 //construction de la transfo initiale a injecter en debut du recalage affine
 imx_aff_param(RIGIDZOOM3D, transfres->param);
 copie_transf_3d(transfres, transfoini);
 transfoini->param=REALLOC(transfoini->param,transfoini->nb_param+3,double);
 transfoini->nb_param=rigidz_to_affine_decouple_3d(transfoini->param);
 transfoini->typetrans=AFFINEDECOUPLE3D;

 //recalage affine
 aff_log("\n TRANSFO AFFINE \n");
 {
  mini=imx_matching_simple_3d_p(imtref, imtreca, imtres, distance, interpol, minimisation, transfres, transfoini, FieldOfResearch, matchPrecision, NULL, NULL, champ1);
 }

 free_field3d(champ1);
 }

 //on remet le tout dans le repere des images originales
 transfres->nb_param=affine_decouple_to_affine_3d(transfres->param);
 transfres->typetrans=AFFINE3D;
 imx_copie_dim_param_to_transf3d_p(imref, transfres);
 //if (preProcessingRecaTrf) comb_transf_3d(preProcessingRecaTrf,transfres,transfres,0);
 //if (preProcessingRefTrf) comb_transf_3d(transfres,preProcessingRefTrf,transfres,0);

 //if (preProcessingRecaTrf) free_transf3d(preProcessingRecaTrf);
 //if (preProcessingRefTrf) free_transf3d(preProcessingRefTrf);
 if (imtres) free_grphic3d(imtres);
 if (imtref!=imref) free_grphic3d(imtref);
 if (imtreca!=imreca) free_grphic3d(imtreca);
 if (imtmp) free_grphic3d(imtmp);
 
 // ------------------CALCUL DE L'IMAGE A AFFICHER ET DE LA TRANSFO A SAUVEGARDER----------------------
 imx_aff_param(transfres->typetrans, transfres->param);
 imx_calc_transfo_sav_image_res_3d_p(imref, imreca, imres, transfres, save_type);

 //liberation memoire des variables dynamiques
 if (transfoini) free_transf3d(transfoini);

 tfin=time(NULL);
 //affichage du temps de calcul total
 {
  imx_aff_temps_calcul(tdebut, tfin);
 }
 end_log();

 return(1);  
}


