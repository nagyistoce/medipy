/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*******************************************************************************/
/*!	   \file:		gradient_distance_3d.c
***
***		project:	Imagix 2.01
***
***
***		\brief description:    
***
***
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***
********************************************************************************/

#include <config.h>
#include <stdio.h>
#include <stdlib.h> 	/* Needed when you use GET_FLOAT	*/
#include <math.h> 	/* Needed when you use math function	*/
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

#include "math/imx_bspline.h"	/* pour les fonctions d'interpolation */

#include "noyau/io/imx_log.h"
#include "math/oper_3d.h"
#include "traitement/trai_3d.h"

#include "recalage/gradient_distance_3d.h"


/*******************************************************************************
********************************************************************************
************ GRADIENT DES TRANSFORMATIONS ERREUR QUADRATIQUE *******************
********************************************************************************
*******************************************************************************/
/*! \ingroup GradientFonction @{ */
/*******************************************************************************
**     gradient_rigid_quad_3d(I2,Im,nb_param,param,gradI2,grad)
*/
/*!
**     Calcul du gradient de la transfo rigide
**     \param I2 : image final = imref
**     \param Im : image transforme de I1 (Imreca) a l'iteration k
**	   \param nb_param  : inutilise
**	   \param param : tableau de parametres (tetax, tetay etc....)
**     \param gradI2= gradient de l'image I2
**     \param grad : vecteur 1D contenant le gradient de l'energie
**	   \retval 1
*******************************************************************************/
int	gradient_rigid_quad_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP)
{
 double	xc,yc,zc,x,y,z,gx,gy,gz,tetax,tetay,tetaz,diff;
 double cxx,sxx,cyy,syy,czz,szz;
 double **R,**DRX,**DRY,**DRZ,**RT,**RTDRX,**RTDRY,**RTDRZ;
 int 	wdth,hght,dpth,i,j,k;
 vector3d ***data;

 imx_gradient_3d_p(Im,gradI2,2,4);

 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 data=gradI2->raw;

 tetax=param[0];tetay=param[1];tetaz=param[2];
 xc=param[6];yc=param[7];zc=param[8];
 cxx=cos(tetax);sxx=sin(tetax);cyy=cos(tetay);syy=sin(tetay);czz=cos(tetaz);szz=sin(tetaz);
 /*calcul de la matrice de rotation R*/
 R=alloc_dmatrix(3,3);
 R[0][0]=cyy*czz;		R[0][1]=-cyy*szz;		R[0][2]=-syy;
 R[1][0]=cxx*szz-sxx*syy*czz;	R[1][1]=cxx*czz+sxx*syy*szz;	R[1][2]=-sxx*cyy;
 R[2][0]=sxx*szz+cxx*syy*czz;	R[2][1]=sxx*czz-cxx*syy*szz;	R[2][2]=cxx*cyy;
 /*calcul des matrices DRX,DRY et DRZ derive de R par rapport a tetx,tetay et tetaz*/
 DRX=alloc_dmatrix(3,3);
 DRX[0][0]=0.0;			DRX[0][1]=0.0;			DRX[0][2]=-0.0;
 DRX[1][0]=-sxx*szz-cxx*syy*czz;DRX[1][1]=-sxx*czz+cxx*syy*szz;	DRX[1][2]=-cxx*cyy;
 DRX[2][0]=cxx*szz-sxx*syy*czz;	DRX[2][1]=cxx*czz+sxx*syy*szz;	DRX[2][2]=-sxx*cyy;
 DRY=alloc_dmatrix(3,3);

 DRY[0][0]=-syy*czz;		DRY[0][1]=syy*szz;		DRY[0][2]=-cyy;
 DRY[1][0]=-sxx*cyy*czz;	DRY[1][1]=sxx*cyy*szz;		DRY[1][2]=sxx*syy;
 DRY[2][0]=cxx*cyy*czz;		DRY[2][1]=-cxx*cyy*szz;		DRY[2][2]=-cxx*syy;
 DRZ=alloc_dmatrix(3,3);
 DRZ[0][0]=-cyy*szz;		DRZ[0][1]=-cyy*czz;		DRZ[0][2]=0.0;
 DRZ[1][0]=cxx*czz+sxx*syy*szz;	DRZ[1][1]=-cxx*szz+sxx*syy*czz;	DRZ[1][2]=0.0;
 DRZ[2][0]=sxx*czz-cxx*syy*szz;	DRZ[2][1]=-sxx*szz-cxx*syy*czz;	DRZ[2][2]=0.0;
 /*calcul de RT: transpose= inverse de R*/
 RT=matrix_transp(R,3,3);free_dmatrix(R,3,3);
 /*calcul de RTDRX, RTDRY, RTDRZ: RTDRX=RT*DRX*/
 RTDRX=multiply_matrix(RT,3,3,DRX,3,3);free_dmatrix(DRX,3,3);
 RTDRY=multiply_matrix(RT,3,3,DRY,3,3);free_dmatrix(DRY,3,3);
 RTDRZ=multiply_matrix(RT,3,3,DRZ,3,3);free_dmatrix(DRZ,3,3);

 /*mise a zero du vecteur grad*/
 for (i=0;i<9;i++) grad[i]=0.0;

 /*calcul du gradient*/
 for (i=0;i<wdth;i++)
 {
  x=(double)i-xc;
  for (j=0;j<hght;j++)
  {
   y=(double)j-yc;
    for (k=0;k<dpth;k++)
    {
     z=(double)k-zc;

     gx=data[i][j][k].x;gy=data[i][j][k].y;gz=data[i][j][k].z;
     diff=2.0*(double)(Im->mri[i][j][k]-I2->mri[i][j][k]);
     if (diff!=0.0)
     {
      if (gx!=0.0)
      {
       grad[0]+=diff*gx*(RTDRX[0][0]*x+RTDRX[0][1]*y+RTDRX[0][2]*z);
       grad[1]+=diff*gx*(RTDRY[0][0]*x+RTDRY[0][1]*y+RTDRY[0][2]*z);
       grad[2]+=diff*gx*(RTDRZ[0][0]*x+RTDRZ[0][1]*y+RTDRZ[0][2]*z);
      }
       if (gy!=0.0)
      {
       grad[0]+=diff*gy*(RTDRX[1][0]*x+RTDRX[1][1]*y+RTDRX[1][2]*z);
       grad[1]+=diff*gy*(RTDRY[1][0]*x+RTDRY[1][1]*y+RTDRY[1][2]*z);
       grad[2]+=diff*gy*(RTDRZ[1][0]*x+RTDRZ[1][1]*y+RTDRZ[1][2]*z);
      }
      if (gz!=0.0)
      {
       grad[0]+=diff*gz*(RTDRX[2][0]*x+RTDRX[2][1]*y+RTDRX[2][2]*z);
       grad[1]+=diff*gz*(RTDRY[2][0]*x+RTDRY[2][1]*y+RTDRY[2][2]*z);
       grad[2]+=diff*gz*(RTDRZ[2][0]*x+RTDRZ[2][1]*y+RTDRZ[2][2]*z);
      }

      grad[3]+=diff*(gx*RT[0][0]+gy*RT[1][0]+gz*RT[2][0]);
      grad[4]+=diff*(gx*RT[0][1]+gy*RT[1][1]+gz*RT[2][1]);
      grad[5]+=diff*(gx*RT[0][2]+gy*RT[1][2]+gz*RT[2][2]);
     }
    }
  }
 }

 for (i=0;i<6;i++) grad[i]=grad[i]/(double)(wdth*hght*dpth);

 free_dmatrix(RTDRX,3,3);free_dmatrix(RTDRY,3,3);free_dmatrix(RTDRZ,3,3);free_dmatrix(RT,3,3);
 return(1);
}
/*! @} */


/*! \ingroup GradientFonction @{ */
/*******************************************************************************
**     gradient_rigid_quad2_3d(I2,Im,nb_param,param,gradI2,grad)
**
**     Calcul du gradient de la transfo rigide
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2
**     grad: vecteur 1D contenant le gradient de l'energie
*******************************************************************************/
int	gradient_rigid_quad2_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP)
{
 double	xc,yc,zc,x,y,z,gx,gy,gz,tetax,tetay,tetaz,diff;
 double cxx,sxx,cyy,syy,czz,szz;
 double **R,**DRX,**DRY,**DRZ,**RT,**RTDRX,**RTDRY,**RTDRZ;
 int 	wdth,hght,dpth,i,j,k;
 vector3d ***data;

 imx_gradient_3d_p(Im,gradI2,2,4);

 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 data=gradI2->raw;

 tetax=param[0];tetay=param[1];tetaz=param[2];
 xc=param[6];yc=param[7];zc=param[8];
 cxx=cos(tetax);sxx=sin(tetax);cyy=cos(tetay);syy=sin(tetay);czz=cos(tetaz);szz=sin(tetaz);
 /*calcul de la matrice de rotation R*/
 R=alloc_dmatrix(3,3);
 R[0][0]=cyy*czz;		R[0][1]=-cyy*szz;		R[0][2]=-syy;
 R[1][0]=cxx*szz-sxx*syy*czz;	R[1][1]=cxx*czz+sxx*syy*szz;	R[1][2]=-sxx*cyy;
 R[2][0]=sxx*szz+cxx*syy*czz;	R[2][1]=sxx*czz-cxx*syy*szz;	R[2][2]=cxx*cyy;
 /*calcul des matrices DRX,DRY et DRZ derive de R par rapport a tetx,tetay et tetaz*/
 DRX=alloc_dmatrix(3,3);
 DRX[0][0]=0.0;			DRX[0][1]=0.0;			DRX[0][2]=-0.0;
 DRX[1][0]=-sxx*szz-cxx*syy*czz;DRX[1][1]=-sxx*czz+cxx*syy*szz;	DRX[1][2]=-cxx*cyy;
 DRX[2][0]=cxx*szz-sxx*syy*czz;	DRX[2][1]=cxx*czz+sxx*syy*szz;	DRX[2][2]=-sxx*cyy;
 DRY=alloc_dmatrix(3,3);
 DRY[0][0]=-syy*czz;		DRY[0][1]=syy*szz;		DRY[0][2]=-cyy;
 DRY[1][0]=-sxx*cyy*czz;	DRY[1][1]=sxx*cyy*szz;		DRY[1][2]=sxx*syy;
 DRY[2][0]=cxx*cyy*czz;		DRY[2][1]=-cxx*cyy*szz;		DRY[2][2]=-cxx*syy;
 DRZ=alloc_dmatrix(3,3);
 DRZ[0][0]=-cyy*szz;		DRZ[0][1]=-cyy*czz;		DRZ[0][2]=0.0;
 DRZ[1][0]=cxx*czz+sxx*syy*szz;	DRZ[1][1]=-cxx*szz+sxx*syy*czz;	DRZ[1][2]=0.0;
 DRZ[2][0]=sxx*czz-cxx*syy*szz;	DRZ[2][1]=-sxx*szz-cxx*syy*czz;	DRZ[2][2]=0.0;
 /*calcul de RT: transpose= inverse de R*/
 RT=matrix_transp(R,3,3);free_dmatrix(R,3,3);
 /*calcul de RTDRX, RTDRY, RTDRZ: RTDRX=RT*DRX*/
 RTDRX=multiply_matrix(RT,3,3,DRX,3,3);
 RTDRY=multiply_matrix(RT,3,3,DRY,3,3);
 RTDRZ=multiply_matrix(RT,3,3,DRZ,3,3);

 /*mise a zero du vecteur grad*/
 for (i=0;i<9;i++) grad[i]=0.0;

 /*calcul du gradient*/
 for (i=0;i<wdth;i++)
 {
  x=(double)i-xc;
  for (j=0;j<hght;j++)
  {
   y=(double)j-yc;
    for (k=0;k<dpth;k++)
    {
     z=(double)k-zc;
     gx=data[i][j][k].x;gy=data[i][j][k].y;gz=data[i][j][k].z;
     diff=2.0*(double)(Im->mri[i][j][k]-I2->mri[i][j][k]);
     if (diff!=0.0)
     {
      if (gx!=0.0)
      {
       grad[0]+=diff*gx*(DRX[0][0]*x+DRX[0][1]*y+DRX[0][2]*z);
       grad[1]+=diff*gx*(DRY[0][0]*x+DRY[0][1]*y+DRY[0][2]*z);
       grad[2]+=diff*gx*(DRZ[0][0]*x+DRZ[0][1]*y+DRZ[0][2]*z);
      }
       if (gy!=0.0)
      {
       grad[0]+=diff*gy*(DRX[1][0]*x+DRX[1][1]*y+DRX[1][2]*z);
       grad[1]+=diff*gy*(DRY[1][0]*x+DRY[1][1]*y+DRY[1][2]*z);
       grad[2]+=diff*gy*(DRZ[1][0]*x+DRZ[1][1]*y+DRZ[1][2]*z);
      }
      if (gz!=0.0)
      {
       grad[0]+=diff*gz*(DRX[2][0]*x+DRX[2][1]*y+DRX[2][2]*z);
       grad[1]+=diff*gz*(DRY[2][0]*x+DRY[2][1]*y+DRY[2][2]*z);
       grad[2]+=diff*gz*(DRZ[2][0]*x+DRZ[2][1]*y+DRZ[2][2]*z);
      }

      grad[3]+=diff*gx;
      grad[4]+=diff*gy;
      grad[5]+=diff*gz;
     }
    }
  }
 }

 for (i=0;i<6;i++) grad[i]=grad[i]/(double)(wdth*hght*dpth);
 free_dmatrix(DRX,3,3);free_dmatrix(DRY,3,3);free_dmatrix(DRZ,3,3);
 free_dmatrix(RTDRX,3,3);free_dmatrix(RTDRY,3,3);free_dmatrix(RTDRZ,3,3);free_dmatrix(RT,3,3);
 return(1);
}
/*! @} */


/*! \ingroup GradientFonction @{ */
/*******************************************************************************
**     gradient_rigidz_quad_3d(I2,Im,nb_param,param,gradI2,grad)
*/
/*!
**     Calcul du gradient de la transfo rigide+zoom
**     \param I2 : image final = imref
**     \param Im : image transforme de I1 (Imreca) a l'iteration k
**	   \param nb_param  : inutilise
**	   \param param : tableau de parametres (tetax, tetay etc....)
**     \param gradI2= gradient de l'image I2
**     \param grad : vecteur 1D contenant le gradient de l'energie
**	   \retval 1
******************************************************************************/
int	gradient_rigidz_quad_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP)
{
 double	xc,yc,zc,x,y,z,gx,gy,gz,sx,sy,sz,tetax,tetay,tetaz,diff;
 double cxx,sxx,cyy,syy,czz,szz;
 double **R,**SI,**S,**RT,**SIRT,**DRX,**DRY,**DRZ,**DRXS,**DRYS,**DRZS;
 double **SIRTDRXS,**SIRTDRYS,**SIRTDRZS;
 int 	wdth,hght,dpth,i,j,k;
 vector3d ***data;

 imx_gradient_3d_p(Im,gradI2,2,4);

 wdth=I2->width;hght=I2->height;dpth=I2->depth;

 data=gradI2->raw;

 tetax=param[0];tetay=param[1];tetaz=param[2];
 sx=param[3];sy=param[4];sz=param[5];
 xc=param[9];yc=param[10];zc=param[11];
 cxx=cos(tetax);sxx=sin(tetax);cyy=cos(tetay);syy=sin(tetay);czz=cos(tetaz);szz=sin(tetaz);
 /*calcul de la matrice de rotation R*/
 R=alloc_dmatrix(3,3);
 R[0][0]=cyy*czz;		R[0][1]=-cyy*szz;		R[0][2]=-syy;
 R[1][0]=cxx*szz-sxx*syy*czz;	R[1][1]=cxx*czz+sxx*syy*szz;	R[1][2]=-sxx*cyy;
 R[2][0]=sxx*szz+cxx*syy*czz;	R[2][1]=sxx*czz-cxx*syy*szz;	R[2][2]=cxx*cyy;
 /*calcul des matrices S et SI*/
 S=alloc_dmatrix(3,3);SI=alloc_dmatrix(3,3);
 for (i=0;i<3;i++) for (j=0;j<3;j++) S[i][j]=SI[i][j]=0.0;
 S[0][0]=sx;S[1][1]=sy;S[2][2]=sz;
 for (i=0;i<3;i++) SI[i][i]=1.0/S[i][i];
 /*calcul des matrices DRX,DRY et DRZ derive de R par rapport a tetx,tetay et tetaz*/
 DRX=alloc_dmatrix(3,3);
 DRX[0][0]=0.0;			DRX[0][1]=0.0;			DRX[0][2]=-0.0;
 DRX[1][0]=-sxx*szz-cxx*syy*czz;DRX[1][1]=-sxx*czz+cxx*syy*szz;	DRX[1][2]=-cxx*cyy;
 DRX[2][0]=cxx*szz-sxx*syy*czz;	DRX[2][1]=cxx*czz+sxx*syy*szz;	DRX[2][2]=-sxx*cyy;
 DRY=alloc_dmatrix(3,3);
 DRY[0][0]=-syy*czz;		DRY[0][1]=syy*szz;		DRY[0][2]=-cyy;
 DRY[1][0]=-sxx*cyy*czz;	DRY[1][1]=sxx*cyy*szz;		DRY[1][2]=sxx*syy;
 DRY[2][0]=cxx*cyy*czz;		DRY[2][1]=-cxx*cyy*szz;		DRY[2][2]=-cxx*syy;
 DRZ=alloc_dmatrix(3,3);
 DRZ[0][0]=-cyy*szz;		DRZ[0][1]=-cyy*czz;		DRZ[0][2]=0.0;
 DRZ[1][0]=cxx*czz+sxx*syy*szz;	DRZ[1][1]=-cxx*szz+sxx*syy*czz;	DRZ[1][2]=0.0;
 DRZ[2][0]=sxx*czz-cxx*syy*szz;	DRZ[2][1]=-sxx*szz-cxx*syy*czz;	DRZ[2][2]=0.0;
 /*calcul de RT: transpose= inverse de R*/
 RT=matrix_transp(R,3,3);free_dmatrix(R,3,3);
 /*calcul de SIRT=SI*RT*/
 SIRT=multiply_matrix(SI,3,3,RT,3,3);free_dmatrix(RT,3,3);free_dmatrix(SI,3,3);
 /*calcul de DRXS, DRYS, DRZS: DRXS=DRX*S*/
 DRXS=multiply_matrix(DRX,3,3,S,3,3);free_dmatrix(DRX,3,3);
 DRYS=multiply_matrix(DRY,3,3,S,3,3);free_dmatrix(DRY,3,3);
 DRZS=multiply_matrix(DRZ,3,3,S,3,3);free_dmatrix(DRZ,3,3);free_dmatrix(S,3,3);
 /*calcul de SIRTDRXS, SIRTDRYS, SIRTDRZS: SIRTDRXS=SI*RT*DRX*S*/
 SIRTDRXS=multiply_matrix(SIRT,3,3,DRXS,3,3);free_dmatrix(DRXS,3,3);
 SIRTDRYS=multiply_matrix(SIRT,3,3,DRYS,3,3);free_dmatrix(DRYS,3,3);
 SIRTDRZS=multiply_matrix(SIRT,3,3,DRZS,3,3);free_dmatrix(DRZS,3,3);



 /*mise a zero du vecteur grad*/
 for (i=0;i<12;i++) grad[i]=0.0;

 /*calcul du gradient*/
 for (i=0;i<wdth;i++)
 {
  x=(double)i-xc;
  for (j=0;j<hght;j++)
  {
   y=(double)j-yc;
    for (k=0;k<dpth;k++)
    {
     z=(double)k-zc;
     gx=data[i][j][k].x;gy=data[i][j][k].y;gz=data[i][j][k].z;
     diff=2.0*(double)(Im->mri[i][j][k]-I2->mri[i][j][k]);
     if (diff!=0.0)
     {
      if (gx!=0.0)
      {
       grad[0]+=diff*gx*(SIRTDRXS[0][0]*x+SIRTDRXS[0][1]*y+SIRTDRXS[0][2]*z);
       grad[1]+=diff*gx*(SIRTDRYS[0][0]*x+SIRTDRYS[0][1]*y+SIRTDRYS[0][2]*z);
       grad[2]+=diff*gx*(SIRTDRZS[0][0]*x+SIRTDRZS[0][1]*y+SIRTDRZS[0][2]*z);
      }
       if (gy!=0.0)
      {
       grad[0]+=diff*gy*(SIRTDRXS[1][0]*x+SIRTDRXS[1][1]*y+SIRTDRXS[1][2]*z);
       grad[1]+=diff*gy*(SIRTDRYS[1][0]*x+SIRTDRYS[1][1]*y+SIRTDRYS[1][2]*z);
       grad[2]+=diff*gy*(SIRTDRZS[1][0]*x+SIRTDRZS[1][1]*y+SIRTDRZS[1][2]*z);
      }
      if (gz!=0.0)
      {
       grad[0]+=diff*gz*(SIRTDRXS[2][0]*x+SIRTDRXS[2][1]*y+SIRTDRXS[2][2]*z);
       grad[1]+=diff*gz*(SIRTDRYS[2][0]*x+SIRTDRYS[2][1]*y+SIRTDRYS[2][2]*z);
       grad[2]+=diff*gz*(SIRTDRZS[2][0]*x+SIRTDRZS[2][1]*y+SIRTDRZS[2][2]*z);
      }

      grad[3]+=diff*gx/sx;grad[4]+=diff*gy/sy;grad[5]+=diff*gz/sz;
      grad[6]+=diff*(gx*SIRT[0][0]+gy*SIRT[1][0]+gz*SIRT[2][0]);
      grad[7]+=diff*(gx*SIRT[0][1]+gy*SIRT[1][1]+gz*SIRT[2][1]);
      grad[8]+=diff*(gx*SIRT[0][2]+gy*SIRT[1][2]+gz*SIRT[2][2]);
     }
    }
  }
 }

 for (i=0;i<9;i++) grad[i]=grad[i]/(double)(wdth*hght*dpth);

 free_dmatrix(SIRTDRXS,3,3);free_dmatrix(SIRTDRYS,3,3);free_dmatrix(SIRTDRZS,3,3);free_dmatrix(SIRT,3,3);
 return(1);
}
/*! @} */


/*! \ingroup GradientFonction @{ */
/*******************************************************************************
**     gradient_rigidz_quad2_3d(I2,Im,nb_param,param,gradI2,grad)
**
**     Calcul du gradient de la transfo rigide+zoom
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2
**     grad: vecteur 1D contenant le gradient de l'energie
*******************************************************************************/
int	gradient_rigidz_quad2_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP)
{
 double	xc,yc,zc,x,y,z,gx,gy,gz,sx,sy,sz,tetax,tetay,tetaz,diff;
 double cxx,sxx,cyy,syy,czz,szz;
 double **R,**SI,**S,**RT,**SIRT,**DRX,**DRY,**DRZ,**DRXS,**DRYS,**DRZS;
 double **SIRTDRXS,**SIRTDRYS,**SIRTDRZS;
 int 	wdth,hght,dpth,i,j,k;
 vector3d ***data;

 imx_gradient_3d_p(Im,gradI2,2,4);

 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 data=gradI2->raw;

 tetax=param[0];tetay=param[1];tetaz=param[2];
 sx=param[3];sy=param[4];sz=param[5];
 xc=param[9];yc=param[10];zc=param[11];
 cxx=cos(tetax);sxx=sin(tetax);cyy=cos(tetay);syy=sin(tetay);czz=cos(tetaz);szz=sin(tetaz);
 /*calcul de la matrice de rotation R*/
 R=alloc_dmatrix(3,3);
 R[0][0]=cyy*czz;		R[0][1]=-cyy*szz;		R[0][2]=-syy;
 R[1][0]=cxx*szz-sxx*syy*czz;	R[1][1]=cxx*czz+sxx*syy*szz;	R[1][2]=-sxx*cyy;
 R[2][0]=sxx*szz+cxx*syy*czz;	R[2][1]=sxx*czz-cxx*syy*szz;	R[2][2]=cxx*cyy;
 /*calcul des matrices S et SI*/
 S=alloc_dmatrix(3,3);SI=alloc_dmatrix(3,3);
 for (i=0;i<3;i++) for (j=0;j<3;j++) S[i][j]=SI[i][j]=0.0;
 S[0][0]=sx;S[1][1]=sy;S[2][2]=sz;
 for (i=0;i<3;i++) SI[i][i]=1.0/S[i][i];
 /*calcul des matrices DRX,DRY et DRZ derive de R par rapport a tetx,tetay et tetaz*/
 DRX=alloc_dmatrix(3,3);
 DRX[0][0]=0.0;			DRX[0][1]=0.0;			DRX[0][2]=-0.0;
 DRX[1][0]=-sxx*szz-cxx*syy*czz;DRX[1][1]=-sxx*czz+cxx*syy*szz;	DRX[1][2]=-cxx*cyy;
 DRX[2][0]=cxx*szz-sxx*syy*czz;	DRX[2][1]=cxx*czz+sxx*syy*szz;	DRX[2][2]=-sxx*cyy;
 DRY=alloc_dmatrix(3,3);
 DRY[0][0]=-syy*czz;		DRY[0][1]=syy*szz;		DRY[0][2]=-cyy;
 DRY[1][0]=-sxx*cyy*czz;	DRY[1][1]=sxx*cyy*szz;		DRY[1][2]=sxx*syy;
 DRY[2][0]=cxx*cyy*czz;		DRY[2][1]=-cxx*cyy*szz;		DRY[2][2]=-cxx*syy;
 DRZ=alloc_dmatrix(3,3);
 DRZ[0][0]=-cyy*szz;		DRZ[0][1]=-cyy*czz;		DRZ[0][2]=0.0;
 DRZ[1][0]=cxx*czz+sxx*syy*szz;	DRZ[1][1]=-cxx*szz+sxx*syy*czz;	DRZ[1][2]=0.0;
 DRZ[2][0]=sxx*czz-cxx*syy*szz;	DRZ[2][1]=-sxx*szz-cxx*syy*czz;	DRZ[2][2]=0.0;
 /*calcul de RT: transpose= inverse de R*/
 RT=matrix_transp(R,3,3);free_dmatrix(R,3,3);
 /*calcul de SIRT=SI*RT*/
 SIRT=multiply_matrix(SI,3,3,RT,3,3);free_dmatrix(RT,3,3);free_dmatrix(SI,3,3);
 /*calcul de DRXS, DRYS, DRZS: DRXS=DRX*S*/
 DRXS=multiply_matrix(DRX,3,3,S,3,3);free_dmatrix(DRX,3,3);
 DRYS=multiply_matrix(DRY,3,3,S,3,3);free_dmatrix(DRY,3,3);
 DRZS=multiply_matrix(DRZ,3,3,S,3,3);free_dmatrix(DRZ,3,3);free_dmatrix(S,3,3);
 /*calcul de SIRTDRXS, SIRTDRYS, SIRTDRZS: SIRTDRXS=SI*RT*DRX*S*/
 SIRTDRXS=multiply_matrix(SIRT,3,3,DRXS,3,3);
 SIRTDRYS=multiply_matrix(SIRT,3,3,DRYS,3,3);
 SIRTDRZS=multiply_matrix(SIRT,3,3,DRZS,3,3);



 /*mise a zero du vecteur grad*/
 for (i=0;i<12;i++) grad[i]=0.0;

 /*calcul du gradient*/
 for (i=0;i<wdth;i++)

 {
  x=(double)i-xc;
  for (j=0;j<hght;j++)
  {
   y=(double)j-yc;
    for (k=0;k<dpth;k++)
    {
     z=(double)k-zc;
     gx=data[i][j][k].x;gy=data[i][j][k].y;gz=data[i][j][k].z;
     diff=2.0*(double)(Im->mri[i][j][k]-I2->mri[i][j][k]);
     if (diff!=0.0)
     {
      if (gx!=0.0)
      {
       grad[0]+=diff*gx*(DRXS[0][0]*x+DRXS[0][1]*y+DRXS[0][2]*z);
       grad[1]+=diff*gx*(DRYS[0][0]*x+DRYS[0][1]*y+DRYS[0][2]*z);
       grad[2]+=diff*gx*(DRZS[0][0]*x+DRZS[0][1]*y+DRZS[0][2]*z);
      }
       if (gy!=0.0)
      {
       grad[0]+=diff*gy*(DRXS[1][0]*x+DRXS[1][1]*y+DRXS[1][2]*z);
       grad[1]+=diff*gy*(DRYS[1][0]*x+DRYS[1][1]*y+DRYS[1][2]*z);
       grad[2]+=diff*gy*(DRZS[1][0]*x+DRZS[1][1]*y+DRZS[1][2]*z);
      }
      if (gz!=0.0)
      {
       grad[0]+=diff*gz*(DRXS[2][0]*x+DRXS[2][1]*y+DRXS[2][2]*z);
       grad[1]+=diff*gz*(DRYS[2][0]*x+DRYS[2][1]*y+DRYS[2][2]*z);
       grad[2]+=diff*gz*(DRZS[2][0]*x+DRZS[2][1]*y+DRZS[2][2]*z);
      }

      grad[3]+=diff*gx/sx;grad[4]+=diff*gy/sy;grad[5]+=diff*gz/sz;
      grad[6]+=diff*gx;
      grad[7]+=diff*gy;
      grad[8]+=diff*gz;
     }
    }
  }
 }

 for (i=0;i<9;i++) grad[i]=grad[i]/(double)(wdth*hght*dpth);

 free_dmatrix(DRXS,3,3);free_dmatrix(DRYS,3,3);free_dmatrix(DRZS,3,3);
 free_dmatrix(SIRTDRXS,3,3);free_dmatrix(SIRTDRYS,3,3);free_dmatrix(SIRTDRZS,3,3);free_dmatrix(SIRT,3,3);
 return(1);
}
/*! @} */



/*! \ingroup GradientFonction @{ */
/*******************************************************************************
**     gradient_affine_quad_3d(I2,Im,nb_param,param,gradI2,grad)
*/
/*!
**     Calcul du gradient de la transfo affine
**     \param I2 : image final = imref
**     \param Im : image transforme de I1 (Imreca) a l'iteration k
**	   \param nb_param  : inutilise
**	   \param param : tableau de parametres (tetax, tetay etc....)
**     \param gradI2= gradient de l'image I2
**     \param grad : vecteur 1D contenant le gradient de l'energie
**	   \retval 1
******************************************************************************/
int	gradient_affine_quad_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP)
{
 double	xc,yc,zc,x,y,z,gx,gy,gz,diff,a0,a1,a2;
 int 	wdth,hght,dpth,i,j,k;
 double **A,**AI;
 vector3d ***data;

 imx_gradient_3d_p(Im,gradI2,2,4);

 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 data=gradI2->raw;

 /*remplissage de A*/
 A=alloc_dmatrix(3,3);
 A[0][0]=param[0];A[0][1]=param[1];A[0][1]=param[2];
 A[1][0]=param[3];A[1][1]=param[4];A[1][1]=param[5];
 A[2][0]=param[6];A[2][1]=param[7];A[2][1]=param[8];
 xc=param[12];yc=param[13];zc=param[14];

 /*calcul de AI inverse de A*/
 AI=matrix_inversion(A,3);free_dmatrix(A,3,3);

 /*mise a zero du vecteur grad*/
 for (i=0;i<15;i++) grad[i]=0.0;

 /*calcul du gradient*/
 for (i=0;i<wdth;i++)
 {
  x=(double)i-xc;
  for (j=0;j<hght;j++)
  {
   y=(double)j-yc;
   for (k=0;k<dpth;k++)
   {
    z=(double)k-zc;

    gx=data[i][j][k].x;gy=data[i][j][k].y;gz=data[i][j][k].z;
    diff=2.0*(double)(Im->mri[i][j][k]-I2->mri[i][j][k]);
    a0=gx*AI[0][0]+gy*AI[1][0]+gz*AI[2][0];
    a1=gx*AI[0][1]+gy*AI[1][1]+gz*AI[2][1];
    a2=gx*AI[0][2]+gy*AI[1][2]+gz*AI[2][2];
    grad[0]+=diff*a0*x;grad[1]+=diff*a1*x;grad[2]+=diff*a2*x;
    grad[3]+=diff*a0*y;grad[4]+=diff*a1*y;grad[5]+=diff*a2*y;
    grad[6]+=diff*a0*z;grad[7]+=diff*a1*z;grad[8]+=diff*a2*z;
    grad[9]+=diff*a0;grad[10]+=diff*a1;grad[11]+=diff*a2;
   }
  }
 }

 for (i=0;i<12;i++) grad[i]=grad[i]/(double)(wdth*hght*dpth);

 return(1);
}
/*! @} */


/*! \ingroup GradientFonction @{ */
/*******************************************************************************
**     gradient_affine_quad2_3d(I2,Im,nb_param,param,gradI2,grad)
**
**     Calcul du gradient de la transfo affine
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2
**     grad: vecteur 1D contenant le gradient de l'energie
*******************************************************************************/
int	gradient_affine_quad2_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP)
{
 double	xc,yc,zc,x,y,z,gx,gy,gz,diff,a0,a1,a2;
 int 	wdth,hght,dpth,i,j,k;
 double **A,**AI;
 vector3d ***data;

 imx_gradient_3d_p(Im,gradI2,2,4);

 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 data=gradI2->raw;

 /*remplissage de A*/
 A=alloc_dmatrix(3,3);
 A[0][0]=param[0];A[0][1]=param[1];A[0][1]=param[2];
 A[1][0]=param[3];A[1][1]=param[4];A[1][1]=param[5];
 A[2][0]=param[6];A[2][1]=param[7];A[2][1]=param[8];
 xc=param[12];yc=param[13];zc=param[14];

 /*calcul de AI inverse de A*/
 AI=matrix_inversion(A,3);free_dmatrix(A,3,3);

 /*mise a zero du vecteur grad*/
 for (i=0;i<15;i++) grad[i]=0.0;

 /*calcul du gradient*/
 for (i=0;i<wdth;i++)
 {
  x=(double)i-xc;
  for (j=0;j<hght;j++)
  {
   y=(double)j-yc;
   for (k=0;k<dpth;k++)
   {
    z=(double)k-zc;

    gx=data[i][j][k].x;gy=data[i][j][k].y;gz=data[i][j][k].z;
    diff=2.0*(double)(Im->mri[i][j][k]-I2->mri[i][j][k]);
    a0=gx;
    a1=gy;
    a2=gz;
    grad[0]+=diff*a0*x;grad[1]+=diff*a1*x;grad[2]+=diff*a2*x;
    grad[3]+=diff*a0*y;grad[4]+=diff*a1*y;grad[5]+=diff*a2*y;
    grad[6]+=diff*a0*z;grad[7]+=diff*a1*z;grad[8]+=diff*a2*z;
    grad[9]+=diff*a0;grad[10]+=diff*a1;grad[11]+=diff*a2;
   }
  }
 }

 for (i=0;i<12;i++) grad[i]=grad[i]/(double)(wdth*hght*dpth);

 return(1);
}
/*! @} */


/*! \ingroup GradientFonction @{ */
/*******************************************************************************
**     gradient_base_quad_3d(I2,Im,nb_param,param,gradI2,grad)
*/
/*!
**     Calcul du gradient de la transfo sur la base de fonctions
**     \param I2 : image final = imref
**     \param Im : image transforme de I1 (Imreca) a l'iteration k
**	   \param nb_param  : inutilise
**	   \param param : tableau de parametres (tetax, tetay etc....)
**     \param gradI2= gradient de l'image I2
**     \param grad : vecteur 1D contenant le gradient de l'energie
**	   \retval 1
******************************************************************************/
int gradient_base_quad_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP)
{
 int 	i,j,k,l,wdth,hght,dpth;
 int 	i0,i1,j0,j1,k0,k1;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double aa1,aa2,aa3,bb1,bb2,bb3,cc1,cc2,cc3,d23,d31,d12;
 double A1,A2,A3,B1,B2,B3,C1,C2,C3;
 double ax,ay,az,gx,gy,gz,diff,det=0.,detmin=0.,detm=0.,temp;
 /*double *Jtab;*/
 vector3d ***data;

 /*Jtab=CALLOC(nb_param,double);detm=0.0;*/

 imx_gradient_3d_p(Im,gradI2,2,4);

 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 data=gradI2->raw;
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;



 /*mise a zero du gradient*/
 for (i=0;i<nb_param;i++) grad[i]=0.0;

 detmin=1.0;
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    gx=data[i][j][k].x;gy=data[i][j][k].y;gz=data[i][j][k].z;
    diff=2.0*(double)(Im->mri[i][j][k]-I2->mri[i][j][k]);

    /*initialisation de la jacobienne a zero en ce point*/


    aa1=1.0;aa2=aa3=bb1=0.0;bb2=1.0;bb3=cc1=cc2=0.0;cc3=1.0;

    for (l=0;l<nb_param/3;l++)
    {
     i0=x0[l];i1=x1[l];j0=y0[l];j1=y1[l];k0=z0[l];k1=z1[l];
     if (i>=i0 && j>=j0 && k>=k0 && i<i1 && j<j1 && k<k1)
     {
      ax=param[3*l];ay=param[3*l+1];az=param[3*l+2];
      aa1=aa1+ax*dfx[i-i0]*fy[j-j0]*fz[k-k0];
      bb1=bb1+ax*fx[i-i0]*dfy[j-j0]*fz[k-k0];
      cc1=cc1+ax*fx[i-i0]*fy[j-j0]*dfz[k-k0];
      aa2=aa2+ay*dfx[i-i0]*fy[j-j0]*fz[k-k0];
      bb2=bb2+ay*fx[i-i0]*dfy[j-j0]*fz[k-k0];
      cc2=cc2+ay*fx[i-i0]*fy[j-j0]*dfz[k-k0];
      aa3=aa3+az*dfx[i-i0]*fy[j-j0]*fz[k-k0];
      bb3=bb3+az*fx[i-i0]*dfy[j-j0]*fz[k-k0];
      cc3=cc3+az*fx[i-i0]*fy[j-j0]*dfz[k-k0];
     }
    }
    /*on inverse la matrice*/
    d23=bb2*cc3-bb3*cc2;d31=bb3*cc1-bb1*cc3;d12=bb1*cc2-bb2*cc1;
    det=aa1*d23+aa2*d31+aa3*d12;

    if (det<detmin) detmin=det;
    if (det<detm) detm=det;
    if (det!=0.0)
     {
      /*calcul des coeffs de la matrice inverse*/
      A1=d23;A2=aa3*cc2-aa2*cc3;A3=aa2*bb3-aa3*bb2;
      B1=d31;B2=aa1*cc3-aa3*cc1;B3=aa3*bb1-aa1*bb3;
      C1=d12;C2=aa2*cc1-aa1*cc2;C3=aa1*bb2-aa2*bb1;

      for (l=0;l<nb_param/3;l++)
      {
       i0=x0[l];i1=x1[l];j0=y0[l];j1=y1[l];k0=z0[l];k1=z1[l];
       if (i>=i0 && j>=j0 && k>=k0 && i<i1 && j<j1 && k<k1)
       {
	temp=diff*fx[i-i0]*fy[j-j0]*fz[k-k0]/det;
        grad[3*l]=grad[3*l]+temp*(A1*gx+A2*gy+A3*gz);
        grad[3*l+1]=grad[3*l+1]+temp*(B1*gx+B2*gy+B3*gz);
        grad[3*l+2]=grad[3*l+2]+temp*(C1*gx+C2*gy+C3*gz);
       }
      }
     }
  }

  for (i=0;i<nb_param;i++) grad[i]=grad[i]/(double)(wdth*hght*dpth);

  return(1);
}
/*! @} */


/*! \ingroup GradientFonction @{ */
/*******************************************************************************
**     gradient_base_quad2_3d(I2,Im,nb_param,param,gradI2,grad)
*/
/*!
**     Calcul du gradient de la transfo sur la base de fonctions
**     \param I2 : image final = imref
**     \param Im : image transforme de I1 (Imreca) a l'iteration k
**	   \param nb_param  : nbre d'element de param
**	   \param param : tableau de parametres (tetax, tetay etc....)
**     \param gradI2= gradient de l'image I2
**     \param grad : vecteur 1D contenant le gradient de l'energie
**	   \retval 1
******************************************************************************/
int gradient_base_quad2_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP)
{
 int 	i,j,k,l,wdth,hght,dpth;
 int 	i0,i1,j0,j1,k0,k1;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 /*double aa1,aa2,aa3,bb1,bb2,bb3,cc1,cc2,cc3,d23,d31,d12;
 double A1,A2,A3,B1,B2,B3,C1,C2,C3;*/
 double /*ax,ay,az,*/gx,gy,gz,diff,/*det,detmin,detm,*/temp;
 /*double *Jtab;*/
 vector3d ***data;

 /*Jtab=CALLOC(nb_param,double);detm=0.0;*/

 imx_gradient_3d_p(Im,gradI2,2,4);

 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 data=gradI2->raw;
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;



 /*mise a zero du gradient*/
 for (i=0;i<nb_param;i++) grad[i]=0.0;

 /*on rempli le vecteur gradient*/
 for (l=0;l<nb_param/3;l++)
  {
   i0=x0[l];i1=x1[l];j0=y0[l];j1=y1[l];k0=z0[l];k1=z1[l];
   for (i=i0;i<i1;i++)
    for (j=j0;j<j1;j++)
     for (k=k0;k<k1;k++)
     {
      gx=data[i][j][k].x;gy=data[i][j][k].y;gz=data[i][j][k].z;
      diff=2.0*(double)(Im->mri[i][j][k]-I2->mri[i][j][k]);
      temp=diff*fx[i-i0]*fy[j-j0]*fz[k-k0];
      grad[3*l]=grad[3*l]+temp*gx;
      grad[3*l+1]=grad[3*l+1]+temp*gy;
      grad[3*l+2]=grad[3*l+2]+temp*gz;
    }
  }


  for (i=0;i<nb_param;i++) grad[i]=grad[i]/(double)(wdth*hght*dpth);

  return(1);
}

/*! @}*/


/*******************************************************************************

********************************************************************************
***** GRADIENT DES TRANSFORMATIONS INFORMATION MUTUELLE ET ASSOCIES ************
********************************************************************************
*******************************************************************************/


/*! \ingroup GradientFonction @{ */
/*******************************************************************************
**     gradient_rigid_entropie_conjointe_3d(I2,Im,nb_param,param,gradI2,grad,donnees_VP)
*/
/*!
**     Calcul du gradient de la transfo rigide pour l'entropie conjointe volume partiel
**     \param I2 : image final = imref
**     \param Im : image transforme de I1 (Imreca) a l'iteration k
**	   \param nb_param  : inutilise
**	   \param param : tableau de parametres (tetax, tetay etc....)
**     \param gradI2= gradient de l'image I2
**     \param grad : vecteur 1D contenant le gradient de l'energie selon les parametres
**		respectifs: tetax, tetay, tetaz, tx, ty, et tz
**     \param donnees_VP: les donnees des histos de volume partiel
**	   \retval 1
*******************************************************************************/
int	gradient_rigid_entropie_conjointe_3d(grphic3d *Imref, grphic3d *Imreca, int nb_param, double *param, field3d *gradImreca, double *grad, VP_histos *donnees_VP)
{
 return calcul_gradient_VP_3d(ENTROPIE_CONJOINTE, RIGID3D, Imreca, Imref, nb_param, param, grad, donnees_VP);
}
/*! @} */

/*! \ingroup GradientFonction @{ */
/*******************************************************************************
**     gradient_rigid_IM_3d(I2,Im,nb_param,param,gradI2,grad,donnees_VP)
*/
/*!
**     Calcul du gradient de la transfo rigide pour l'info mutuelle volume partiel
**     \param I2 : image final = imref
**     \param Im : image transforme de I1 (Imreca) a l'iteration k
**	   \param nb_param  : inutilise
**	   \param param : tableau de parametres (tetax, tetay etc....)
**     \param gradI2= gradient de l'image I2
**     \param grad : vecteur 1D contenant le gradient de l'energie selon les parametres
**		respectifs: tetax, tetay, tetaz, tx, ty, et tz
**     \param donnees_VP: les donnees des histos de volume partiel
**	   \retval 1
*******************************************************************************/
int	gradient_rigid_IM_3d(grphic3d *Imref, grphic3d *Imreca, int nb_param, double *param, field3d *gradImreca, double *grad, VP_histos *donnees_VP)
{
 return calcul_gradient_VP_3d(IM, RIGID3D, Imreca, Imref, nb_param, param, grad, donnees_VP);
}
/*! @} */

/*! \ingroup GradientFonction @{ */
/*******************************************************************************
**     gradient_rigid_IMNS_3d(I2,Im,nb_param,param,gradI2,grad,donnees_VP)
*/
/*!
**     Calcul du gradient de la transfo rigide pour l'info mutuelle normalisee de Studholme volume partiel
**     \param I2 : image final = imref
**     \param Im : image transforme de I1 (Imreca) a l'iteration k
**	   \param nb_param  : inutilise
**	   \param param : tableau de parametres (tetax, tetay etc....)
**     \param gradI2= gradient de l'image I2
**     \param grad : vecteur 1D contenant le gradient de l'energie selon les parametres
**		respectifs: tetax, tetay, tetaz, tx, ty, et tz
**     \param donnees_VP: les donnees des histos de volume partiel
**	   \retval 1
*******************************************************************************/
int	gradient_rigid_IMNS_3d(grphic3d *Imref, grphic3d *Imreca, int nb_param, double *param, field3d *gradImreca, double *grad, VP_histos *donnees_VP)
{
 return calcul_gradient_VP_3d(IMNS, RIGID3D, Imreca, Imref, nb_param, param, grad, donnees_VP);
}
/*! @} */

/*! \ingroup GradientFonction @{ */
/*******************************************************************************
**     gradient_rigid_IMNM1_3d(I2,Im,nb_param,param,gradI2,grad,donnees_VP)
*/
/*!
**     Calcul du gradient de la transfo rigide pour l'info mutuelle normalisee de Maes 1 volume partiel
**     \param I2 : image final = imref
**     \param Im : image transforme de I1 (Imreca) a l'iteration k
**	   \param nb_param  : inutilise
**	   \param param : tableau de parametres (tetax, tetay etc....)
**     \param gradI2= gradient de l'image I2
**     \param grad : vecteur 1D contenant le gradient de l'energie selon les parametres
**		respectifs: tetax, tetay, tetaz, tx, ty, et tz
**     \param donnees_VP: les donnees des histos de volume partiel
**	   \retval 1
*******************************************************************************/
int	gradient_rigid_IMNM1_3d(grphic3d *Imref, grphic3d *Imreca, int nb_param, double *param, field3d *gradImreca, double *grad, VP_histos *donnees_VP)
{
 return calcul_gradient_VP_3d(IMNM1, RIGID3D, Imreca, Imref, nb_param, param, grad, donnees_VP);
}
/*! @} */

/*! \ingroup GradientFonction @{ */
/*******************************************************************************
**     gradient_rigid_IMNM2_3d(I2,Im,nb_param,param,gradI2,grad,donnees_VP)
*/
/*!
**     Calcul du gradient de la transfo rigide pour l'info mutuelle normalisee de Maes 2 volume partiel
**     \param I2 : image final = imref
**     \param Im : image transforme de I1 (Imreca) a l'iteration k
**	   \param nb_param  : inutilise
**	   \param param : tableau de parametres (tetax, tetay etc....)
**     \param gradI2= gradient de l'image I2
**     \param grad : vecteur 1D contenant le gradient de l'energie selon les parametres
**		respectifs: tetax, tetay, tetaz, tx, ty, et tz
**     \param donnees_VP: les donnees des histos de volume partiel
**	   \retval 1
*******************************************************************************/
int	gradient_rigid_IMNM2_3d(grphic3d *Imref, grphic3d *Imreca, int nb_param, double *param, field3d *gradImreca, double *grad, VP_histos *donnees_VP)
{
 return calcul_gradient_VP_3d(IMNM2, RIGID3D, Imreca, Imref, nb_param, param, grad, donnees_VP);
}
/*! @} */

/*! \ingroup GradientFonction @{ */
/*******************************************************************************
**     gradient_rigidz_entropie_conjointe_3d(I2,Im,nb_param,param,gradI2,grad,donnees_VP)
*/
/*!
**     Calcul du gradient de la transfo rigide pour l'entropie conjointe volume partiel
**     \param I2 : image final = imref
**     \param Im : image transforme de I1 (Imreca) a l'iteration k
**	   \param nb_param  : inutilise
**	   \param param : tableau de parametres (tetax, tetay etc....)
**     \param gradI2= gradient de l'image I2
**     \param grad : vecteur 1D contenant le gradient de l'energie selon les parametres
**		respectifs: tetax, tetay, tetaz, tx, ty, et tz
**     \param donnees_VP: les donnees des histos de volume partiel
**	   \retval 1
*******************************************************************************/
int	gradient_rigidz_entropie_conjointe_3d(grphic3d *Imref, grphic3d *Imreca, int nb_param, double *param, field3d *gradImreca, double *grad, VP_histos *donnees_VP)
{
 return calcul_gradient_VP_3d(ENTROPIE_CONJOINTE, RIGIDZOOM3D, Imreca, Imref, nb_param, param, grad, donnees_VP);
}
/*! @} */

/*! \ingroup GradientFonction @{ */
/*******************************************************************************
**     gradient_rigidz_IM_3d(I2,Im,nb_param,param,gradI2,grad,donnees_VP)
*/
/*!
**     Calcul du gradient de la transfo rigide pour l'info mutuelle volume partiel
**     \param I2 : image final = imref
**     \param Im : image transforme de I1 (Imreca) a l'iteration k
**	   \param nb_param  : inutilise
**	   \param param : tableau de parametres (tetax, tetay etc....)
**     \param gradI2= gradient de l'image I2
**     \param grad : vecteur 1D contenant le gradient de l'energie selon les parametres
**		respectifs: tetax, tetay, tetaz, tx, ty, et tz
**     \param donnees_VP: les donnees des histos de volume partiel
**	   \retval 1
*******************************************************************************/
int	gradient_rigidz_IM_3d(grphic3d *Imref, grphic3d *Imreca, int nb_param, double *param, field3d *gradImreca, double *grad, VP_histos *donnees_VP)
{
 return calcul_gradient_VP_3d(IM, RIGIDZOOM3D, Imreca, Imref, nb_param, param, grad, donnees_VP);
}
/*! @} */

/*! \ingroup GradientFonction @{ */
/*******************************************************************************
**     gradient_rigidz_IMNS_3d(I2,Im,nb_param,param,gradI2,grad,donnees_VP)
*/
/*!
**     Calcul du gradient de la transfo rigide pour l'info mutuelle normalisee de Studholme volume partiel
**     \param I2 : image final = imref
**     \param Im : image transforme de I1 (Imreca) a l'iteration k
**	   \param nb_param  : inutilise
**	   \param param : tableau de parametres (tetax, tetay etc....)
**     \param gradI2= gradient de l'image I2
**     \param grad : vecteur 1D contenant le gradient de l'energie selon les parametres
**		respectifs: tetax, tetay, tetaz, tx, ty, et tz
**     \param donnees_VP: les donnees des histos de volume partiel
**	   \retval 1
*******************************************************************************/
int	gradient_rigidz_IMNS_3d(grphic3d *Imref, grphic3d *Imreca, int nb_param, double *param, field3d *gradImreca, double *grad, VP_histos *donnees_VP)
{
 return calcul_gradient_VP_3d(IMNS, RIGIDZOOM3D, Imreca, Imref, nb_param, param, grad, donnees_VP);
}
/*! @} */

/*! \ingroup GradientFonction @{ */
/*******************************************************************************
**     gradient_rigidz_IMNM1_3d(I2,Im,nb_param,param,gradI2,grad,donnees_VP)
*/
/*!
**     Calcul du gradient de la transfo rigide pour l'info mutuelle normalisee de Maes 1 volume partiel
**     \param I2 : image final = imref
**     \param Im : image transforme de I1 (Imreca) a l'iteration k
**	   \param nb_param  : inutilise
**	   \param param : tableau de parametres (tetax, tetay etc....)
**     \param gradI2= gradient de l'image I2
**     \param grad : vecteur 1D contenant le gradient de l'energie selon les parametres
**		respectifs: tetax, tetay, tetaz, tx, ty, et tz
**     \param donnees_VP: les donnees des histos de volume partiel
**	   \retval 1

*******************************************************************************/
int	gradient_rigidz_IMNM1_3d(grphic3d *Imref, grphic3d *Imreca, int nb_param, double *param, field3d *gradImreca, double *grad, VP_histos *donnees_VP)
{
 return calcul_gradient_VP_3d(IMNM1, RIGIDZOOM3D, Imreca, Imref, nb_param, param, grad, donnees_VP);
}
/*! @} */

/*! \ingroup GradientFonction @{ */
/*******************************************************************************
**     gradient_rigidz_IMNM2_3d(I2,Im,nb_param,param,gradI2,grad,donnees_VP)
*/
/*!
**     Calcul du gradient de la transfo rigide pour l'info mutuelle normalisee de Maes 2 volume partiel
**     \param I2 : image final = imref
**     \param Im : image transforme de I1 (Imreca) a l'iteration k
**	   \param nb_param  : inutilise
**	   \param param : tableau de parametres (tetax, tetay etc....)
**     \param gradI2= gradient de l'image I2
**     \param grad : vecteur 1D contenant le gradient de l'energie selon les parametres
**		respectifs: tetax, tetay, tetaz, tx, ty, et tz
**     \param donnees_VP: les donnees des histos de volume partiel
**	   \retval 1
*******************************************************************************/
int	gradient_rigidz_IMNM2_3d(grphic3d *Imref, grphic3d *Imreca, int nb_param, double *param, field3d *gradImreca, double *grad, VP_histos *donnees_VP)
{
 return calcul_gradient_VP_3d(IMNM2, RIGIDZOOM3D, Imreca, Imref, nb_param, param, grad, donnees_VP);
}
/*! @} */


/*! \ingroup GradientFonction @{ */
/*******************************************************************************
**     gradient_affine_entropie_conjointe_3d(I2,Im,nb_param,param,gradI2,grad,donnees_VP)
*/
/*!
**     Calcul du gradient de la transfo rigide pour l'entropie conjointe volume partiel
**     \param I2 : image final = imref
**     \param Im : image transforme de I1 (Imreca) a l'iteration k
**	   \param nb_param  : inutilise
**	   \param param : tableau de parametres (tetax, tetay etc....)
**     \param gradI2= gradient de l'image I2
**     \param grad : vecteur 1D contenant le gradient de l'energie selon les parametres
**		respectifs: tetax, tetay, tetaz, tx, ty, et tz
**     \param donnees_VP: les donnees des histos de volume partiel
**	   \retval 1
*******************************************************************************/
int	gradient_affine_entropie_conjointe_3d(grphic3d *Imref, grphic3d *Imreca, int nb_param, double *param, field3d *gradImreca, double *grad, VP_histos *donnees_VP)
{
 return calcul_gradient_VP_3d(ENTROPIE_CONJOINTE, AFFINE3D, Imreca, Imref, nb_param, param, grad, donnees_VP);
}
/*! @} */

/*! \ingroup GradientFonction @{ */
/*******************************************************************************
**     gradient_affine_IM_3d(I2,Im,nb_param,param,gradI2,grad,donnees_VP)
*/
/*!
**     Calcul du gradient de la transfo rigide pour l'info mutuelle volume partiel
**     \param I2 : image final = imref
**     \param Im : image transforme de I1 (Imreca) a l'iteration k
**	   \param nb_param  : inutilise
**	   \param param : tableau de parametres (tetax, tetay etc....)
**     \param gradI2= gradient de l'image I2
**     \param grad : vecteur 1D contenant le gradient de l'energie selon les parametres
**		respectifs: tetax, tetay, tetaz, tx, ty, et tz
**     \param donnees_VP: les donnees des histos de volume partiel
**	   \retval 1
*******************************************************************************/
int	gradient_affine_IM_3d(grphic3d *Imref, grphic3d *Imreca, int nb_param, double *param, field3d *gradImreca, double *grad, VP_histos *donnees_VP)
{
 return calcul_gradient_VP_3d(IM, AFFINE3D, Imreca, Imref, nb_param, param, grad, donnees_VP);
}
/*! @} */

/*! \ingroup GradientFonction @{ */
/*******************************************************************************
**     gradient_affine_IMNS_3d(I2,Im,nb_param,param,gradI2,grad,donnees_VP)
*/
/*!
**     Calcul du gradient de la transfo rigide pour l'info mutuelle normalisee de Studholme volume partiel
**     \param I2 : image final = imref
**     \param Im : image transforme de I1 (Imreca) a l'iteration k
**	   \param nb_param  : inutilise
**	   \param param : tableau de parametres (tetax, tetay etc....)
**     \param gradI2= gradient de l'image I2
**     \param grad : vecteur 1D contenant le gradient de l'energie selon les parametres
**		respectifs: tetax, tetay, tetaz, tx, ty, et tz
**     \param donnees_VP: les donnees des histos de volume partiel
**	   \retval 1
*******************************************************************************/
int	gradient_affine_IMNS_3d(grphic3d *Imref, grphic3d *Imreca, int nb_param, double *param, field3d *gradImreca, double *grad, VP_histos *donnees_VP)
{
 return calcul_gradient_VP_3d(IMNS, AFFINE3D, Imreca, Imref, nb_param, param, grad, donnees_VP);
}
/*! @} */

/*! \ingroup GradientFonction @{ */
/*******************************************************************************
**     gradient_affine_IMNM1_3d(I2,Im,nb_param,param,gradI2,grad,donnees_VP)
*/
/*!
**     Calcul du gradient de la transfo rigide pour l'info mutuelle normalisee de Maes 1 volume partiel
**     \param I2 : image final = imref
**     \param Im : image transforme de I1 (Imreca) a l'iteration k
**	   \param nb_param  : inutilise
**	   \param param : tableau de parametres (tetax, tetay etc....)
**     \param gradI2= gradient de l'image I2
**     \param grad : vecteur 1D contenant le gradient de l'energie selon les parametres
**		respectifs: tetax, tetay, tetaz, tx, ty, et tz
**     \param donnees_VP: les donnees des histos de volume partiel
**	   \retval 1
*******************************************************************************/
int	gradient_affine_IMNM1_3d(grphic3d *Imref, grphic3d *Imreca, int nb_param, double *param, field3d *gradImreca, double *grad, VP_histos *donnees_VP)
{
 return calcul_gradient_VP_3d(IMNM1, AFFINE3D, Imreca, Imref, nb_param, param, grad, donnees_VP);
}
/*! @} */

/*! \ingroup GradientFonction @{ */
/*******************************************************************************
**     gradient_affine_IMNM2_3d(I2,Im,nb_param,param,gradI2,grad,donnees_VP)
*/
/*!
**     Calcul du gradient de la transfo rigide pour l'info mutuelle normalisee de Maes 2 volume partiel
**     \param I2 : image final = imref
**     \param Im : image transforme de I1 (Imreca) a l'iteration k
**	   \param nb_param  : inutilise
**	   \param param : tableau de parametres (tetax, tetay etc....)
**     \param gradI2= gradient de l'image I2
**     \param grad : vecteur 1D contenant le gradient de l'energie selon les parametres
**		respectifs: tetax, tetay, tetaz, tx, ty, et tz
**     \param donnees_VP: les donnees des histos de volume partiel
**	   \retval 1
*******************************************************************************/
int	gradient_affine_IMNM2_3d(grphic3d *Imref, grphic3d *Imreca, int nb_param, double *param, field3d *gradImreca, double *grad, VP_histos *donnees_VP)
{
 return calcul_gradient_VP_3d(IMNM2, AFFINE3D, Imreca, Imref, nb_param, param, grad, donnees_VP);
}
/*! @} */

/*******************************************************************************
**     calcul_transformation_derivee_rigid_3d(param, DTRANS)
*/
/*!
**     calcul des derivees de la matrice de transformation rigide selon
**     les parametres de la transformation
**     \param param : les parametres courant de la transformation
**     \param DTRANS : matrices derivees, DTRANS[i] correspond a la matrice 3x4
**            derivee de la matrice de transformation selon le parametre no i
*******************************************************************************/
void calcul_transformation_derivee_rigid_3d(double *param, double *** DTRANS)
{
 double	xc,yc,zc,tetax,tetay,tetaz,tx,ty,tz;
 double cxx,sxx,cyy,syy,czz,szz;

 //recuperation des parametres de la transformation
 tetax=param[0];tetay=param[1];tetaz=param[2];
 tx=param[3];ty=param[4];tz=param[5];
 xc=param[6];yc=param[7];zc=param[8];

 cxx=cos(tetax);sxx=sin(tetax);cyy=cos(tetay);syy=sin(tetay);czz=cos(tetaz);szz=sin(tetaz);

 //calcul des matrices derivees de la matrice de transformation
 //Derivee selon tetax
 DTRANS[0][0][0]=0.0;                  DTRANS[0][0][1]=0.0;                  DTRANS[0][0][2]=0.0;      DTRANS[0][0][3]=0.0;
 DTRANS[0][1][0]=-sxx*szz-cxx*syy*czz; DTRANS[0][1][1]=-sxx*czz+cxx*syy*szz; DTRANS[0][1][2]=-cxx*cyy; DTRANS[0][1][3]=-(-sxx*szz-cxx*syy*czz)*xc-(-sxx*czz+cxx*syy*szz)*yc+cxx*cyy*zc;
 DTRANS[0][2][0]=cxx*szz-sxx*syy*czz;  DTRANS[0][2][1]=cxx*czz+sxx*syy*szz;  DTRANS[0][2][2]=-sxx*cyy; DTRANS[0][2][3]=-(cxx*szz-sxx*syy*czz)*xc-(cxx*czz+sxx*syy*szz)*yc+sxx*cyy*zc;

 //Derivee selon tetay
 DTRANS[1][0][0]=-syy*czz;     DTRANS[1][0][1]=syy*szz;      DTRANS[1][0][2]=-cyy;     DTRANS[1][0][3]=-(-syy*czz)*xc-(syy*szz)*yc+cyy*zc;
 DTRANS[1][1][0]=-sxx*cyy*czz; DTRANS[1][1][1]=sxx*cyy*szz;  DTRANS[1][1][2]=sxx*syy;  DTRANS[1][1][3]=-(-sxx*cyy*czz)*xc-(sxx*cyy*szz)*yc-sxx*syy*zc;
 DTRANS[1][2][0]=cxx*cyy*czz;  DTRANS[1][2][1]=-cxx*cyy*szz; DTRANS[1][2][2]=-cxx*syy; DTRANS[1][2][3]=-(cxx*cyy*czz)*xc-(-cxx*cyy*szz)*yc+cxx*syy*zc;

 //Derivee selon tetaz
 DTRANS[2][0][0]=-cyy*szz;            DTRANS[2][0][1]=-cyy*czz;             DTRANS[2][0][2]=0.0; DTRANS[2][0][3]=-(-cyy*szz)*xc-(-cyy*czz)*yc;
 DTRANS[2][1][0]=cxx*czz+sxx*syy*szz; DTRANS[2][1][1]=-cxx*szz+sxx*syy*czz; DTRANS[2][1][2]=0.0; DTRANS[2][1][3]=-(cxx*czz+sxx*syy*szz)*xc-(-cxx*szz+sxx*syy*czz)*yc;
 DTRANS[2][2][0]=sxx*czz-cxx*syy*szz; DTRANS[2][2][1]=-sxx*szz-cxx*syy*czz; DTRANS[2][2][2]=0.0; DTRANS[2][2][3]=-(sxx*czz-cxx*syy*szz)*xc-(-sxx*szz-cxx*syy*czz)*yc;

 //Derivee selon tx
 DTRANS[3][0][0]=0.0; DTRANS[3][0][1]=0.0; DTRANS[3][0][2]=0.0; DTRANS[3][0][3]=1.0;
 DTRANS[3][1][0]=0.0; DTRANS[3][1][1]=0.0; DTRANS[3][1][2]=0.0; DTRANS[3][1][3]=0.0;
 DTRANS[3][2][0]=0.0; DTRANS[3][2][1]=0.0; DTRANS[3][2][2]=0.0; DTRANS[3][2][3]=0.0;

 //Derivee selon ty
 DTRANS[4][0][0]=0.0; DTRANS[4][0][1]=0.0; DTRANS[4][0][2]=0.0; DTRANS[4][0][3]=0.0;
 DTRANS[4][1][0]=0.0; DTRANS[4][1][1]=0.0; DTRANS[4][1][2]=0.0; DTRANS[4][1][3]=1.0;
 DTRANS[4][2][0]=0.0; DTRANS[4][2][1]=0.0; DTRANS[4][2][2]=0.0; DTRANS[4][2][3]=0.0;

 //Derivee selon tz
 DTRANS[5][0][0]=0.0; DTRANS[5][0][1]=0.0; DTRANS[5][0][2]=0.0; DTRANS[5][0][3]=0.0;
 DTRANS[5][1][0]=0.0; DTRANS[5][1][1]=0.0; DTRANS[5][1][2]=0.0; DTRANS[5][1][3]=0.0;
 DTRANS[5][2][0]=0.0; DTRANS[5][2][1]=0.0; DTRANS[5][2][2]=0.0; DTRANS[5][2][3]=1.0;

 return;
}

/*******************************************************************************
**     calcul_transformation_derivee_rigid_zoom_3d(param, DTRANS)
*/
/*!
**     calcul des derivees de la matrice de transformation rigide+zoom selon
**     les parametres de la transformation
**     \param param : les parametres courant de la transformation
**     \param DTRANS : matrices derivees, DTRANS[i] correspond a la matrice 3x4
**            derivee de la matrice de transformation selon le parametre no i
*******************************************************************************/
void calcul_transformation_derivee_rigid_zoom_3d(double *param, double *** DTRANS)
{
 double	xc,yc,zc,tetax,tetay,tetaz,tx,ty,tz,sx,sy,sz;
 double cxx,sxx,cyy,syy,czz,szz;

 //recuperation des parametres de la transformation
 tetax=param[0];tetay=param[1];tetaz=param[2];
 sx=param[3];sy=param[4];sz=param[5];
 tx=param[6];ty=param[7];tz=param[8];
 xc=param[9];yc=param[10];zc=param[11];

 cxx=cos(tetax);sxx=sin(tetax);cyy=cos(tetay);syy=sin(tetay);czz=cos(tetaz);szz=sin(tetaz);

 //calcul des matrices derivees de la matrice de transformation
 //Derivee selon tetax
 DTRANS[0][0][0]=0.0;                       DTRANS[0][0][1]=0.0;                       DTRANS[0][0][2]=0.0;           DTRANS[0][0][3]=0.0;
 DTRANS[0][1][0]=sx*(-sxx*szz-cxx*syy*czz); DTRANS[0][1][1]=sy*(-sxx*czz+cxx*syy*szz); DTRANS[0][1][2]=sz*(-cxx*cyy); DTRANS[0][1][3]=-sx*(-sxx*szz-cxx*syy*czz)*xc-sy*(-sxx*czz+cxx*syy*szz)*yc+sz*cxx*cyy*zc;
 DTRANS[0][2][0]=sx*(cxx*szz-sxx*syy*czz);  DTRANS[0][2][1]=sy*(cxx*czz+sxx*syy*szz);  DTRANS[0][2][2]=sz*(-sxx*cyy); DTRANS[0][2][3]=-sx*(cxx*szz-sxx*syy*czz)*xc-sy*(cxx*czz+sxx*syy*szz)*yc+sz*sxx*cyy*zc;

 //Derivee selon tetay
 DTRANS[1][0][0]=sx*(-syy*czz);     DTRANS[1][0][1]=sy*(syy*szz);      DTRANS[1][0][2]=sz*(-cyy);     DTRANS[1][0][3]=-sx*(-syy*czz)*xc-sy*(syy*szz)*yc+sz*cyy*zc;
 DTRANS[1][1][0]=sx*(-sxx*cyy*czz); DTRANS[1][1][1]=sy*(sxx*cyy*szz);  DTRANS[1][1][2]=sz*(sxx*syy);  DTRANS[1][1][3]=-sx*(-sxx*cyy*czz)*xc-sy*(sxx*cyy*szz)*yc-sz*sxx*syy*zc;
 DTRANS[1][2][0]=sx*(cxx*cyy*czz);  DTRANS[1][2][1]=sy*(-cxx*cyy*szz); DTRANS[1][2][2]=sz*(-cxx*syy); DTRANS[1][2][3]=-sx*(cxx*cyy*czz)*xc-sy*(-cxx*cyy*szz)*yc+sz*cxx*syy*zc;

 //Derivee selon tetaz
 DTRANS[2][0][0]=sx*(-cyy*szz);            DTRANS[2][0][1]=sy*(-cyy*czz);             DTRANS[2][0][2]=0.0; DTRANS[2][0][3]=-sx*(-cyy*szz)*xc-sy*(-cyy*czz)*yc;
 DTRANS[2][1][0]=sx*(cxx*czz+sxx*syy*szz); DTRANS[2][1][1]=sy*(-cxx*szz+sxx*syy*czz); DTRANS[2][1][2]=0.0; DTRANS[2][1][3]=-sx*(cxx*czz+sxx*syy*szz)*xc-sy*(-cxx*szz+sxx*syy*czz)*yc;
 DTRANS[2][2][0]=sx*(sxx*czz-cxx*syy*szz); DTRANS[2][2][1]=sy*(-sxx*szz-cxx*syy*czz); DTRANS[2][2][2]=0.0; DTRANS[2][2][3]=-sx*(sxx*czz-cxx*syy*szz)*xc-sy*(-sxx*szz-cxx*syy*czz)*yc;

 //Derivee selon sx
 DTRANS[3][0][0]=cyy*czz;             DTRANS[3][0][1]=0.0; DTRANS[3][0][2]=0.0; DTRANS[3][0][3]=-cyy*czz*xc;
 DTRANS[3][1][0]=cxx*szz-sxx*syy*czz; DTRANS[3][1][1]=0.0; DTRANS[3][1][2]=0.0; DTRANS[3][1][3]=-(cxx*szz-sxx*syy*czz)*xc;
 DTRANS[3][2][0]=sxx*szz+cxx*syy*czz; DTRANS[3][2][1]=0.0; DTRANS[3][2][2]=0.0; DTRANS[3][2][3]=-(sxx*szz+cxx*syy*czz)*xc;

 //Derivee selon sy
 DTRANS[4][0][0]=0.0; DTRANS[4][0][1]=-cyy*szz;            DTRANS[4][0][2]=0.0; DTRANS[4][0][3]=cyy*szz*yc;
 DTRANS[4][1][0]=0.0; DTRANS[4][1][1]=cxx*czz+sxx*syy*szz; DTRANS[4][1][2]=0.0; DTRANS[4][1][3]=-(cxx*czz+sxx*syy*szz)*yc;
 DTRANS[4][2][0]=0.0; DTRANS[4][2][1]=sxx*czz-cxx*syy*szz; DTRANS[4][2][2]=0.0; DTRANS[4][2][3]=-(sxx*czz-cxx*syy*szz)*yc;

 //Derivee selon sz
 DTRANS[5][0][0]=0.0; DTRANS[5][0][1]=0.0; DTRANS[5][0][2]=-syy;     DTRANS[5][0][3]=syy*zc;
 DTRANS[5][1][0]=0.0; DTRANS[5][1][1]=0.0; DTRANS[5][1][2]=-sxx*cyy; DTRANS[5][1][3]=sxx*cyy*zc;
 DTRANS[5][2][0]=0.0; DTRANS[5][2][1]=0.0; DTRANS[5][2][2]=cxx*cyy;  DTRANS[5][2][3]=-cxx*cyy*zc;

 //Derivee selon tx
 DTRANS[6][0][0]=0.0; DTRANS[6][0][1]=0.0; DTRANS[6][0][2]=0.0; DTRANS[6][0][3]=1.0;
 DTRANS[6][1][0]=0.0; DTRANS[6][1][1]=0.0; DTRANS[6][1][2]=0.0; DTRANS[6][1][3]=0.0;
 DTRANS[6][2][0]=0.0; DTRANS[6][2][1]=0.0; DTRANS[6][2][2]=0.0; DTRANS[6][2][3]=0.0;

 //Derivee selon ty
 DTRANS[7][0][0]=0.0; DTRANS[7][0][1]=0.0; DTRANS[7][0][2]=0.0; DTRANS[7][0][3]=0.0;
 DTRANS[7][1][0]=0.0; DTRANS[7][1][1]=0.0; DTRANS[7][1][2]=0.0; DTRANS[7][1][3]=1.0;
 DTRANS[7][2][0]=0.0; DTRANS[7][2][1]=0.0; DTRANS[7][2][2]=0.0; DTRANS[7][2][3]=0.0;

 //Derivee selon tz
 DTRANS[8][0][0]=0.0; DTRANS[8][0][1]=0.0; DTRANS[8][0][2]=0.0; DTRANS[8][0][3]=0.0;
 DTRANS[8][1][0]=0.0; DTRANS[8][1][1]=0.0; DTRANS[8][1][2]=0.0; DTRANS[8][1][3]=0.0;
 DTRANS[8][2][0]=0.0; DTRANS[8][2][1]=0.0; DTRANS[8][2][2]=0.0; DTRANS[8][2][3]=1.0;

 return;
}

/*******************************************************************************
**     calcul_transformation_derivee_affine_3d(param, DTRANS)
*/
/*!
**     calcul des derivees de la matrice de transformation affine selon
**     les parametres de la transformation
**     \param param : les parametres courant de la transformation
**     \param DTRANS : matrices derivees, DTRANS[i] correspond a la matrice 3x4
**            derivee de la matrice de transformation selon le parametre no i
*******************************************************************************/
void calcul_transformation_derivee_affine_3d(double *param, double *** DTRANS)
{
 double	xc,yc,zc;

 xc=param[12];yc=param[13];zc=param[14];

 DTRANS[0][0][0]=1.0; DTRANS[0][0][1]=0.0; DTRANS[0][0][2]=0.0; DTRANS[0][0][3]=-xc;
 DTRANS[0][1][0]=0.0; DTRANS[0][1][1]=0.0; DTRANS[0][1][2]=0.0; DTRANS[0][1][3]=0.0;
 DTRANS[0][2][0]=0.0; DTRANS[0][2][1]=0.0; DTRANS[0][2][2]=0.0; DTRANS[0][2][3]=0.0;

 DTRANS[1][0][0]=0.0; DTRANS[1][0][1]=1.0; DTRANS[1][0][2]=0.0; DTRANS[1][0][3]=-yc;
 DTRANS[1][1][0]=0.0; DTRANS[1][1][1]=0.0; DTRANS[1][1][2]=0.0; DTRANS[1][1][3]=0.0;
 DTRANS[1][2][0]=0.0; DTRANS[1][2][1]=0.0; DTRANS[1][2][2]=0.0; DTRANS[1][2][3]=0.0;

 DTRANS[2][0][0]=0.0; DTRANS[2][0][1]=0.0; DTRANS[2][0][2]=1.0; DTRANS[2][0][3]=-zc;
 DTRANS[2][1][0]=0.0; DTRANS[2][1][1]=0.0; DTRANS[2][1][2]=0.0; DTRANS[2][1][3]=0.0;
 DTRANS[2][2][0]=0.0; DTRANS[2][2][1]=0.0; DTRANS[2][2][2]=0.0; DTRANS[2][2][3]=0.0;

 DTRANS[3][0][0]=0.0; DTRANS[3][0][1]=0.0; DTRANS[3][0][2]=0.0; DTRANS[3][0][3]=0.0;
 DTRANS[3][1][0]=1.0; DTRANS[3][1][1]=0.0; DTRANS[3][1][2]=0.0; DTRANS[3][1][3]=-xc;
 DTRANS[3][2][0]=0.0; DTRANS[3][2][1]=0.0; DTRANS[3][2][2]=0.0; DTRANS[3][2][3]=0.0;

 DTRANS[4][0][0]=0.0; DTRANS[4][0][1]=0.0; DTRANS[4][0][2]=0.0; DTRANS[4][0][3]=0.0;
 DTRANS[4][1][0]=0.0; DTRANS[4][1][1]=1.0; DTRANS[4][1][2]=0.0; DTRANS[4][1][3]=-yc;
 DTRANS[4][2][0]=0.0; DTRANS[4][2][1]=0.0; DTRANS[4][2][2]=0.0; DTRANS[4][2][3]=0.0;

 DTRANS[5][0][0]=0.0; DTRANS[5][0][1]=0.0; DTRANS[5][0][2]=0.0; DTRANS[5][0][3]=0.0;
 DTRANS[5][1][0]=0.0; DTRANS[5][1][1]=0.0; DTRANS[5][1][2]=1.0; DTRANS[5][1][3]=-zc;
 DTRANS[5][2][0]=0.0; DTRANS[5][2][1]=0.0; DTRANS[5][2][2]=0.0; DTRANS[5][2][3]=0.0;

 DTRANS[6][0][0]=0.0; DTRANS[6][0][1]=0.0; DTRANS[6][0][2]=0.0; DTRANS[6][0][3]=0.0;
 DTRANS[6][1][0]=0.0; DTRANS[6][1][1]=0.0; DTRANS[6][1][2]=0.0; DTRANS[6][1][3]=0.0;
 DTRANS[6][2][0]=1.0; DTRANS[6][2][1]=0.0; DTRANS[6][2][2]=0.0; DTRANS[6][2][3]=-xc;

 DTRANS[7][0][0]=0.0; DTRANS[7][0][1]=0.0; DTRANS[7][0][2]=0.0; DTRANS[7][0][3]=0.0;
 DTRANS[7][1][0]=0.0; DTRANS[7][1][1]=0.0; DTRANS[7][1][2]=0.0; DTRANS[7][1][3]=0.0;
 DTRANS[7][2][0]=0.0; DTRANS[7][2][1]=1.0; DTRANS[7][2][2]=0.0; DTRANS[7][2][3]=-yc;

 DTRANS[8][0][0]=0.0; DTRANS[8][0][1]=0.0; DTRANS[8][0][2]=0.0; DTRANS[8][0][3]=0.0;
 DTRANS[8][1][0]=0.0; DTRANS[8][1][1]=0.0; DTRANS[8][1][2]=0.0; DTRANS[8][1][3]=0.0;
 DTRANS[8][2][0]=0.0; DTRANS[8][2][1]=0.0; DTRANS[8][2][2]=1.0; DTRANS[8][2][3]=-zc;

 DTRANS[9][0][0]=0.0; DTRANS[9][0][1]=0.0; DTRANS[9][0][2]=0.0; DTRANS[9][0][3]=1.0;
 DTRANS[9][1][0]=0.0; DTRANS[9][1][1]=0.0; DTRANS[9][1][2]=0.0; DTRANS[9][1][3]=0.0;
 DTRANS[9][2][0]=0.0; DTRANS[9][2][1]=0.0; DTRANS[9][2][2]=0.0; DTRANS[9][2][3]=0.0;

 DTRANS[10][0][0]=0.0; DTRANS[10][0][1]=0.0; DTRANS[10][0][2]=0.0; DTRANS[10][0][3]=0.0;
 DTRANS[10][1][0]=0.0; DTRANS[10][1][1]=0.0; DTRANS[10][1][2]=0.0; DTRANS[10][1][3]=1.0;
 DTRANS[10][2][0]=0.0; DTRANS[10][2][1]=0.0; DTRANS[10][2][2]=0.0; DTRANS[10][2][3]=0.0;

 DTRANS[11][0][0]=0.0; DTRANS[11][0][1]=0.0; DTRANS[11][0][2]=0.0; DTRANS[11][0][3]=0.0;
 DTRANS[11][1][0]=0.0; DTRANS[11][1][1]=0.0; DTRANS[11][1][2]=0.0; DTRANS[11][1][3]=0.0;
 DTRANS[11][2][0]=0.0; DTRANS[11][2][1]=0.0; DTRANS[11][2][2]=0.0; DTRANS[11][2][3]=1.0;

 return;
}

/*******************************************************************************
**     gradient_entropie_conjointe_3d(nbVoxelsOverlap, donnees_VP, dEnergie)
*/
/*!
**     calcul des derivees de l'entropie conjointe selon Tij ou T est la matrice
**     3x4 de la transformation
**     \param nbVoxelsOverlap : le nombre de voxels dans la zone de recouvrement
**            des images
**     \param donnees_VP : les donnees des histos de volume partiel
**     \param dEnergie : la matrice allouee des derivees de l'energie, dEnergie[i][j]
**            represente la derivee de l'energie selon Tij
**     \retval 1
*******************************************************************************/
int gradient_entropie_conjointe_3d(double nbVoxelsOverlap, VP_histos *donnees_VP, double **dEnergie)
{
 int i,j;
 int nbcl1, nbcl2;
 double ** histogramme_conjoint=NULL;
 double *histogramme_marginal1=NULL, *histogramme_marginal2=NULL;
 double entropie_conjointe=0.0;

/////// INITIALISATIONS ///////
 nbcl1 = donnees_VP->nb_cl1;
 nbcl2 = donnees_VP->nb_cl2;
 histogramme_conjoint = donnees_VP->histo_conjoint;
 histogramme_marginal1 = donnees_VP->histo_marginal1;
 histogramme_marginal2 = donnees_VP->histo_marginal2;
 entropie_conjointe = donnees_VP->entropie_conjointe;

 //calcul de l'histogramme gradient d'energie
 for (i=1;i<nbcl1;i++)
  for (j=1;j<nbcl2;j++)
   if (histogramme_conjoint[i][j]!=0)
    dEnergie[i][j]=-(log(histogramme_conjoint[i][j]/nbVoxelsOverlap)+entropie_conjointe);

 return(1);
}

/*******************************************************************************
**     gradient_entropie_conjointe_3d(nbVoxelsOverlap, donnees_VP, dEnergie)
*/
/*!
**     calcul des derivees de l'info mutuelle selon Tij ou T est la matrice
**     3x4 de la transformation
**     \param nbVoxelsOverlap : le nombre de voxels dans la zone de recouvrement
**            des images
**     \param donnees_VP : les donnees des histos de volume partiel
**     \param dEnergie : la matrice allouee des derivees de l'energie, dEnergie[i][j]
**            represente la derivee de l'energie selon Tij
**     \retval 1
*******************************************************************************/
int gradient_IM_3d(double nbVoxelsOverlap, VP_histos *donnees_VP, double **dEnergie)
{
  int i,j;
 int nbcl1, nbcl2;
 double ** histogramme_conjoint=NULL;
 double *histogramme_marginal1=NULL, *histogramme_marginal2=NULL;
 double info_mutuelle=0.0;

/////// INITIALISATIONS ///////
 nbcl1 = donnees_VP->nb_cl1;
 nbcl2 = donnees_VP->nb_cl2;
 histogramme_conjoint = donnees_VP->histo_conjoint;
 histogramme_marginal1 = donnees_VP->histo_marginal1;
 histogramme_marginal2 = donnees_VP->histo_marginal2;
 info_mutuelle = donnees_VP->entropie1 + donnees_VP->entropie2 - donnees_VP->entropie_conjointe;

 //calcul de l'histogramme gradient d'energie
 for (i=1;i<nbcl1;i++)
  for (j=1;j<nbcl2;j++)
   if (histogramme_conjoint[i][j]!=0)
    dEnergie[i][j]=-(log(nbVoxelsOverlap*histogramme_conjoint[i][j]/(histogramme_marginal1[i]*histogramme_marginal2[j]))-info_mutuelle);

 return(1);
}

/*******************************************************************************
**     gradient_entropie_conjointe_3d(nbVoxelsOverlap, donnees_VP, dEnergie)
*/
/*!
**     calcul des derivees de l'info mutuelle normalisee de Studholme selon Tij
**     ou T est la matrice 3x4 de la transformation
**     \param nbVoxelsOverlap : le nombre de voxels dans la zone de recouvrement
**            des images
**     \param donnees_VP : les donnees des histos de volume partiel
**     \param dEnergie : la matrice allouee des derivees de l'energie, dEnergie[i][j]
**            represente la derivee de l'energie selon Tij
**     \retval 1
*******************************************************************************/
int gradient_IMNS_3d(double nbVoxelsOverlap, VP_histos *donnees_VP, double **dEnergie)
{
  int i,j;
 int nbcl1, nbcl2;
 double ** histogramme_conjoint=NULL;
 double *histogramme_marginal1=NULL, *histogramme_marginal2=NULL;
 double info_mutuelle_norm_studholme=0.0;

/////// INITIALISATIONS ///////
 nbcl1 = donnees_VP->nb_cl1;
 nbcl2 = donnees_VP->nb_cl2;
 histogramme_conjoint = donnees_VP->histo_conjoint;
 histogramme_marginal1 = donnees_VP->histo_marginal1;
 histogramme_marginal2 = donnees_VP->histo_marginal2;
 info_mutuelle_norm_studholme = (donnees_VP->entropie1 + donnees_VP->entropie2)/(donnees_VP->entropie_conjointe);

 //calcul de l'histogramme gradient d'energie
 for (i=1;i<nbcl1;i++)
  for (j=1;j<nbcl2;j++)
   if (histogramme_conjoint[i][j]!=0)
    dEnergie[i][j]=-(info_mutuelle_norm_studholme*log(histogramme_conjoint[i][j]/nbVoxelsOverlap)-log(histogramme_marginal1[i]*histogramme_marginal2[j]/(nbVoxelsOverlap*nbVoxelsOverlap)))/(donnees_VP->entropie_conjointe);

 return(1);
}

/*******************************************************************************
**     gradient_entropie_conjointe_3d(nbVoxelsOverlap, donnees_VP, dEnergie)
*/
/*!
**     calcul des derivees de l'info mutuelle normalisee de Maes 1 selon Tij
**     ou T est la matrice 3x4 de la transformation
**     \param nbVoxelsOverlap : le nombre de voxels dans la zone de recouvrement
**            des images
**     \param donnees_VP : les donnees des histos de volume partiel
**     \param dEnergie : la matrice allouee des derivees de l'energie, dEnergie[i][j]
**            represente la derivee de l'energie selon Tij
**     \retval 1
*******************************************************************************/
int gradient_IMNM1_3d(double nbVoxelsOverlap, VP_histos *donnees_VP, double **dEnergie)
{
 int i,j;
 int nbcl1, nbcl2;
 double ** histogramme_conjoint=NULL;
 double *histogramme_marginal1=NULL, *histogramme_marginal2=NULL;
 double info_mutuelle_norm_maes1=0.0,entropie1=0.0,entropie2=0.0,entropie_conjointe=0.0;

/////// INITIALISATIONS ///////
 nbcl1 = donnees_VP->nb_cl1;
 nbcl2 = donnees_VP->nb_cl2;
 histogramme_conjoint = donnees_VP->histo_conjoint;
 histogramme_marginal1 = donnees_VP->histo_marginal1;
 histogramme_marginal2 = donnees_VP->histo_marginal2;
 entropie_conjointe = donnees_VP->entropie_conjointe;
 entropie1 = donnees_VP->entropie1;
 entropie2 = donnees_VP->entropie2;
 info_mutuelle_norm_maes1 = (entropie1+entropie2-entropie_conjointe)/(entropie1+entropie2);

 //calcul de l'histogramme gradient d'energie
 for (i=1;i<nbcl1;i++)
  for (j=1;j<nbcl2;j++)
   if (histogramme_conjoint[i][j]!=0)
    {
     dEnergie[i][j]=(log(nbVoxelsOverlap*histogramme_conjoint[i][j]/(histogramme_marginal1[i]*histogramme_marginal2[j]))+log(histogramme_marginal1[i]*histogramme_marginal2[j]/(nbVoxelsOverlap*nbVoxelsOverlap))*info_mutuelle_norm_maes1);
     dEnergie[i][j] = -dEnergie[i][j]/(histogramme_marginal1[i]+histogramme_marginal2[j]);
    }

 return(1);
}

/*******************************************************************************
**     gradient_entropie_conjointe_3d(nbVoxelsOverlap, donnees_VP, dEnergie)
*/
/*!
**     calcul des derivees de l'info mutuelle normalisee de Maes 2 selon Tij
**     ou T est la matrice 3x4 de la transformation
**     \param nbVoxelsOverlap : le nombre de voxels dans la zone de recouvrement
**            des images
**     \param donnees_VP : les donnees des histos de volume partiel
**     \param dEnergie : la matrice allouee des derivees de l'energie, dEnergie[i][j]
**            represente la derivee de l'energie selon Tij
**     \retval 1
*******************************************************************************/
int gradient_IMNM2_3d(double nbVoxelsOverlap, VP_histos *donnees_VP, double **dEnergie)
{
  int i,j;
 int nbcl1, nbcl2;
 double ** histogramme_conjoint=NULL;
 double *histogramme_marginal1=NULL, *histogramme_marginal2=NULL;
 double info_mutuelle_norm_maes2=0.0,entropie1=0.0,entropie2=0.0,entropie_conjointe=0.0;

/////// INITIALISATIONS ///////
 nbcl1 = donnees_VP->nb_cl1;
 nbcl2 = donnees_VP->nb_cl2;
 histogramme_conjoint = donnees_VP->histo_conjoint;
 histogramme_marginal1 = donnees_VP->histo_marginal1;
 histogramme_marginal2 = donnees_VP->histo_marginal2;
 entropie_conjointe = donnees_VP->entropie_conjointe;
 entropie1 = donnees_VP->entropie1;
 entropie2 = donnees_VP->entropie2;
 info_mutuelle_norm_maes2 = 2.0*entropie_conjointe-entropie1-entropie2;

 //calcul de l'histogramme gradient d'energie
 for (i=1;i<nbcl1;i++)
  for (j=1;j<nbcl2;j++)
   if (histogramme_conjoint[i][j]!=0)
    dEnergie[i][j]=(log(histogramme_marginal1[i]*histogramme_marginal2[j]/(histogramme_conjoint[i][j]*histogramme_conjoint[i][j]))-info_mutuelle_norm_maes2);

 return(1);
}

/*******************************************************************************
**     calcul_gradient_VP_3d(im_type, transfo, nb_param, param, grad, donnees_VP)
*/
/*!
**     calcul du gradient de l'energie basee sur la volume partiel selon
**     les parametres de la transformation
**     les histogrammes de gradient doivent avoir ete calcules au prealable
**     par la methode 'calcul_histogrammes_gradient_VP'
**     \param im_type : le type d'energie utilisee
**     \param transfo : le type de transformation utilisee
**     \param nb_param : le nombre de parametres
**     \param param : le vecteur des parametres courants
**     \param grad : le vecteur de gradient alloue
**     \param donnees_VP : les donnees des histos de volume partiel
**     \retval 1
*******************************************************************************/
int calcul_gradient_VP_3d(IM3DTYPE im_type, TRANSFO3D transfo, grphic3d *Imreca, grphic3d *Imref, int nb_param, double *param, double *grad, VP_histos *donnees_VP)
{
 int h,i,j,k,l;
 int nbcl1, nbcl2;
 double dH,normGrad;
 double **** histogrammes_gradient_base_IM=NULL;
 double ** histogramme_conjoint=NULL;
 double *histogramme_marginal1=NULL, *histogramme_marginal2=NULL;
 double nbVoxelsOverlap=0.0;
 double *** DTRANS;
 double **dEnergie=NULL;
 int nb_param_util;
 field3d *champ;
 transf_func_t the_transf;
 int wdth,hght,dpth;
 ptr_distance_param dist_param=CALLOC(1, distance_param);

//////// INITIALISATIONS ////////
 switch (transfo)
 {
  case RIGID3D : nb_param_util=6;break;
  case RIGIDZOOM3D : nb_param_util=9;break;
  case AFFINE3D : nb_param_util=12;break;
  default : nb_param_util=6;break;
 }

 switch (transfo)
 {
  case RIGID3D: the_transf=rigid_to_field_3d; break;
  case RIGIDZOOM3D: the_transf=rigidz_to_field_3d; break;
  case AFFINE3D: the_transf=affine_to_field_3d; break;
  default: the_transf=rigid_to_field_3d; break;
 }

 wdth=Imref->width;hght=Imref->height;dpth=Imref->depth;
 champ=cr_field3d(wdth,hght,dpth);

 the_transf(nb_param,param,champ,Imref,Imreca);
 dist_param->imreca=Imreca; dist_param->imref=Imref;
 dist_param->champ=champ; dist_param->donnees_VP=donnees_VP;
// erreur_IM_VP_3d(Imreca, Imref, champ, donnees_VP, NULL, NULL, NULL);
 erreur_IM_VP_3d(dist_param);

 calcul_histogrammes_gradient_VP(Imreca, Imref, champ, donnees_VP, NULL);

 nbcl1 = donnees_VP->nb_cl1;
 nbcl2 = donnees_VP->nb_cl2;
 histogrammes_gradient_base_IM = donnees_VP->histos_gradient;
 histogramme_conjoint = donnees_VP->histo_conjoint;
 histogramme_marginal1 = donnees_VP->histo_marginal1;
 histogramme_marginal2 = donnees_VP->histo_marginal2;

 //on calcule le nb de voxel dans la zone de recouvrement
 for (i=1;i<nbcl1;i++)
  for (j=1;j<nbcl2;j++)
   nbVoxelsOverlap += histogramme_conjoint[i][j];

 //allocation memoire et initialisation
 dEnergie = alloc_dmatrix(nbcl1, nbcl2);
 for (i=1;i<nbcl1;i++)
  for (j=1;j<nbcl2;j++)
   dEnergie[i][j] = 0.0;

 //allocation et initialisation des matrices derivees de la matrice de transformation selon
 //les parametres respectifs: tetax, tetay, tetaz, tx, ty, et tz
 DTRANS = CALLOC(nb_param_util,double **);
 for (i=0;i<nb_param_util;i++)
 {
  DTRANS[i]=CALLOC(3, double *);
  for (j=0;j<3;j++)
  {
   DTRANS[i][j]=CALLOC(4, double);
  }
 }


////CALCUL DU GRADIENT////

 //calcul du gradient de l'energie selon Tij avec T matrice 3x4 de la transformation
 //les derivee sont stockees dans la matrice dEnergie
 switch (im_type)
 {
  case ENTROPIE_CONJOINTE: gradient_entropie_conjointe_3d(nbVoxelsOverlap, donnees_VP, dEnergie); break;
  case IM: gradient_IM_3d(nbVoxelsOverlap, donnees_VP, dEnergie); break;
  case IMNS: gradient_IMNS_3d(nbVoxelsOverlap, donnees_VP, dEnergie); break;
  case IMNM1: gradient_IMNM1_3d(nbVoxelsOverlap, donnees_VP, dEnergie); break;
  case IMNM2: gradient_IMNM2_3d(nbVoxelsOverlap, donnees_VP, dEnergie); break;
  default: gradient_entropie_conjointe_3d(nbVoxelsOverlap, donnees_VP, dEnergie); break;
 }

 //calcul des matrices derivees selon les parametres de la transformation
 //elles sont placees dans DTRANS avec DTRANS[i]=derivee selon le parametre no i
 switch (transfo)
 {
  case RIGID3D : calcul_transformation_derivee_rigid_3d(param, DTRANS);break;
  case RIGIDZOOM3D : calcul_transformation_derivee_rigid_zoom_3d(param, DTRANS);break;
  case AFFINE3D : calcul_transformation_derivee_affine_3d(param, DTRANS);break;
  default : calcul_transformation_derivee_rigid_3d(param, DTRANS);break;
 }

 //mise a zero du vecteur grad
 for (i=0;i<nb_param;i++) grad[i]=0.0;

 //calcul du gradient
 for (h=0;h<nb_param_util;h++)
 {
  //calcul du gradient selon la direction no h
  for (i=1;i<nbcl1;i++)
  {
   for (j=1;j<nbcl2;j++)
   {
     dH = 0.0;

     for (k=0;k<3;k++)
      for (l=0;l<4;l++)
       dH += DTRANS[h][k][l]*histogrammes_gradient_base_IM[k][l][i][j];

     grad[h] += dEnergie[i][j]*dH;
	 }
  }
  grad[h]=grad[h]/nbVoxelsOverlap;
 }

 //normalisation du gradient
 normGrad = 0.0;
 for (i=0;i<nb_param_util;i++)
  normGrad += grad[i]*grad[i];

 normGrad = sqrt(normGrad);

 if (normGrad>0)
  for (i=0;i<nb_param_util;i++)
   grad[i] = grad[i]/normGrad;

//////// LIBERATIONS MEMOIRE ////////
 free_dmatrix(dEnergie,nbcl1,nbcl2);

 for (i=0;i<nb_param_util;i++)
 {
  for (j=0;j<3;j++)
   FREE(DTRANS[i][j]);

  FREE(DTRANS[i]);
 }
 free(DTRANS);

 free_field3d(champ);

 FREE(dist_param);

 return(0);
}


/*------------------------------------------------*/
/*!
**	\brief determine la fonction de gradient
**
**	\param  the_transf : le type de transformation applique
**	\param  dist : le type de distance choisie
**	\retval gradient_func_t : un pointeur sur la fonction gradient
*/
/*------------------------------------------------*/
//gradient_func_t imx_choose_gradient_fct(transf_func_t the_transf, dist_func_t dist)
gradient_func_t imx_choose_gradient_fct(TRANSFO3D transfo_type, dist_func_t dist)
{
 gradient_func_t gradient=NULL;

 if (dist==erreur_IM_VP_3d) //erreur information mutuelle volume partiel
 {
  if (transfo_type==RIGID3D) {gradient=gradient_rigid_IM_3d;}
  else if (transfo_type==RIGIDZOOM3D) {gradient=gradient_rigidz_IM_3d;}
  else if (transfo_type==AFFINE3D) {gradient=gradient_affine_IM_3d;}

 }
 else if (dist==erreur_IMNS_VP_3d) //erreur information mutuelle normalisee de studholme volume partiel
 {
  if (transfo_type==RIGID3D) {gradient=gradient_rigid_IMNS_3d;}
  else if (transfo_type==RIGIDZOOM3D) {gradient=gradient_rigidz_IMNS_3d;}
  else if (transfo_type==AFFINE3D) {gradient=gradient_affine_IMNS_3d;}
 }
 else if (dist==erreur_IMNM1_VP_3d) //erreur information mutuelle normalisee de Maes 1 volume partiel
 {
  if (transfo_type==RIGID3D) {gradient=gradient_rigid_IMNM1_3d;}
  else if (transfo_type==RIGIDZOOM3D) {gradient=gradient_rigidz_IMNM1_3d;}
  else if (transfo_type==AFFINE3D) {gradient=gradient_affine_IMNM1_3d;}
 }
 else if (dist==erreur_IMNM2_VP_3d) /*erreur information mutuelle normalisee de Maes 2 volume partiel*/
 {
  if (transfo_type==RIGID3D) {gradient=gradient_rigid_IMNM2_3d;}
  else if (transfo_type==RIGIDZOOM3D) {gradient=gradient_rigidz_IMNM2_3d;}
  else if (transfo_type==AFFINE3D) {gradient=gradient_affine_IMNM2_3d;}
 }
 else if (dist==erreur_entropie_conjointe_VP_3d) //erreur entropie conjointe volume partiel
 {
  if (transfo_type==RIGID3D) {gradient=gradient_rigid_entropie_conjointe_3d;}
  else if (transfo_type==RIGIDZOOM3D) {gradient=gradient_rigidz_entropie_conjointe_3d;}
  else if (transfo_type==AFFINE3D) {gradient=gradient_affine_entropie_conjointe_3d;}
 }
 else if (dist==erreur_quad_3d) //on utilise l'erreur quadratique
 {
 if (transfo_type==RIGID3D) {gradient=gradient_rigid_quad_3d;}
 else if (transfo_type==RIGIDZOOM3D) {gradient=gradient_rigidz_quad_3d;}
 else if (transfo_type==AFFINE3D) {gradient=gradient_affine_quad_3d;}
 else if (transfo_type==BSPLINE3D) {gradient=gradient_base_quad_3d;}
 }

 return gradient;
}


