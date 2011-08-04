/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/**-----------------------------------------------------------------------
*** 
*** file:  operation.c
***
*** project: Imagix 1.01
***   
***
*** description: Operation mathematique et logique
*** 
*** 
*** Copyright (c) 1993, ULP-IPB Strasbourg.
*** All rights are reserved.
***
***     Last user action by: Mr. ARMSPACH on July 23th 1993
***
***---------------------------------------------------------------------*/
#include <config.h>
#include <stdio.h> 
#include <stdlib.h>  /* Needed when you use GET_FLOAT */
#include <math.h>  /* Needed when you use math function */
#include <string.h> 

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_2d.h"
#include "noyau/imx_lang.h"
#include "math/operation.h"


/*******************************************
** --  op_math() ----------------------------
**
** Op�rations math�matiques
********************************************/

int imx_add (int im_1, int im_2, int im_res);
int imx_sub (int im_1, int im_2, int im_res);
int imx_mul (int im_1, int im_2, int im_res);
int imx_div (int im_1, int im_2, int im_res);
int imx_abs (int im_1, int im_res);
int imx_add_coe (int im_1, float coefficient, int im_res);
int imx_sub_coe (int im_1, float coefficient, int im_res);
int imx_mul_coe (int im_1, float coefficient, int im_res);
int imx_div_coe (int im_1, float coefficient, int im_res);
int imx_and (int im_1, int im_2, int im_res);
int imx_or (int im_1, int im_2, int im_res);
int imx_eg (int im_1, int im_2, int im_res);
int imx_inv (int im_1, int im_res);
int imx_add_p (grphic *im1, grphic *im2, grphic *imres);
int imx_brukermax (float amax, int *maxpixel, float *icomp, float *rcoeff);
int imx_sub_p (grphic *im1, grphic *im2, grphic *imres);
int imx_mul_p (grphic *im1, grphic *im2, grphic *imres);
int imx_div_p (grphic *im1, grphic *im2, grphic *imres);
int imx_mul_coe_p (grphic *im1, float coefficient, grphic *imres);
int imx_inimaxpixel_p (grphic *im1);
int imx_iniminpixel_p (grphic *im1);
int imx_inimaxminpixel_p (grphic *im1);
int imx_iniparaimg_p (grphic *im1);
int imx_copie_param_p (grphic *im1, grphic *imres);

void op_math(void)
{ 
  char *quest[6];
  int answer_nr;
  int i;
  int im_1=0,im_2=0,im_res=0;

 /* Question concerning type of operation  */
 /* and image     */ 
  
  for(i=0;i<6;i++)
    quest[i]=CALLOC(80,char);
  strcpy(quest[0],TEXT0075);
  strcpy(quest[1],TEXT0076);
  strcpy(quest[2],TEXT0077);
  strcpy(quest[3],TEXT0078);
  strcpy(quest[4],TEXT0079);
  strcpy(quest[5],"\0");
 /* Put a question box at position 150,200. */
  answer_nr=GETV_QCM(TEXT0043,(const char **)quest);
 /*     Recherche numero image      */
  if (answer_nr!=4) {
    im_1=GET_PLACE("Image 1:");
    im_2=GET_PLACE("Operation with image 2:");
    im_res=GET_PLACE(TEXT0006);
    }
  else {
    im_1=GET_PLACE("Image 1:");
    im_res=GET_PLACE(TEXT0006);
    }

  switch ((int)answer_nr) {
  
 case 0:  /*  Addition */
  
  imx_add(im_1,im_2,im_res); 
  show_picture(im_res);
  break;

 case 1:  /*  Soustraction */

  imx_sub(im_1,im_2,im_res);
  show_picture(im_res);
  break;
 
 case 2:  /*  Multiplication */

  imx_mul(im_1,im_2,im_res);
  show_picture(im_res);
  break;
 
 case 3:  /*  Division */

  imx_div(im_1,im_2,im_res);
  show_picture(im_res);
  break;
 
 case 4:  /*  Absolu */

  imx_abs(im_1,im_res);
  show_picture(im_res);
  break;
 
 default:
  fprintf(stderr,"Warning: in Mathematical Operation\n");
  break;

  }


}
/*******************************************
** -- op_math_coe() ----------------------------
**
**    Operations mathematiques avec coeeficient
********************************************/
void op_math_coe(void)
{ 
  char *quest[5];
  int answer_nr;
  int i,err;
  int im_1,im_res;
  float coefficient;

 /* Question concerning type of operation  */
 /* and image     */ 
  
  for(i=0;i<5;i++)
    quest[i]=CALLOC(80,char);
  strcpy(quest[0],TEXT0080);
  strcpy(quest[1],TEXT0081);
  strcpy(quest[2],TEXT0082);
  strcpy(quest[3],TEXT0083);
  strcpy(quest[4],"\0");
 /* Put a question box at position 150,200. */
  answer_nr=GETV_QCM(TEXT0043,(const char **)quest);
 /*     Image       */
  im_1=GET_PLACE("Image 1:");
l1:coefficient= GET_FLOAT("Coefficient:", 0, &err);
  im_res=GET_PLACE(TEXT0006);

  switch ((int)answer_nr) {
  
 case 0:  /*  Addition */
  
  imx_add_coe(im_1,coefficient,im_res); 
  show_picture(im_res);
  break;

 case 1:  /*  Soustraction */

  imx_sub_coe(im_1,coefficient,im_res);
  show_picture(im_res);
  break;
 
 case 2:  /*  Multiplication */

  imx_mul_coe(im_1,coefficient,im_res);
  show_picture(im_res);
  break;
 
 case 3:  /*  Division */

  if (coefficient==0) goto l1;
  imx_div_coe(im_1,coefficient,im_res);
  show_picture(im_res);
  break;
 
 default:
  fprintf(stderr,"Warning: in Mathematical Operation with coef\n");
  break;

  }

}
/*******************************************
** --  op_log() ----------------------------
**
**    Operations logiques sue les images
********************************************/
void op_log(void)
{ 
  char *quest[5];
  int answer_nr;
  int i;
  int im_1=0,im_2=0,im_res=0;

 /* Question concerning type of operation  */
 /* and image     */ 
  
  for(i=0;i<5;i++)
    quest[i]=CALLOC(80,char);
  strcpy(quest[0],"And");
  strcpy(quest[1],"Or");
  strcpy(quest[2],"Egalite");
  strcpy(quest[3],"Inversion");
  strcpy(quest[4],"\0");
 /* Put a question box at position 150,200. */
  answer_nr=GETV_QCM(TEXT0043,(const char **)quest);
 /*     Image       */
  im_1=GET_PLACE("Image 1:");
  if (answer_nr!=3)
    im_2=GET_PLACE("op�ration with image 2:");
  im_res=GET_PLACE(TEXT0006);

  switch ((int)answer_nr) {
  
 case 0:  /*  ET logique */
  
  imx_and(im_1,im_2,im_res); 
  show_picture(im_res);
  break;

 case 1:  /*  OU logique */

  imx_or(im_1,im_2,im_res);
  show_picture(im_res);
  break;
 
 case 2:  /*  Egalite */

  imx_eg(im_1,im_2,im_res);
  show_picture(im_res);
  break;
 
 case 3:  /*  Inversion */

  imx_inv(im_1,im_res);
  show_picture(im_res);
  break;
 
 default:
  fprintf(stderr,"Warning: in Mathematical Operation\n");
  break;

  }

}

/*******************************************
** --  imx_add() ----------------------------
**     Additive function
**     im_1,im_2 : im_1 operation with im_2
**     im_res : Result image
********************************************/
int imx_add(int im_1, int im_2, int im_res)
{ 
  grphic *im1,*im2,*imres;

  im1=ptr_img(im_1);
  im2=ptr_img(im_2);
  imres=ptr_img(im_res);

  imx_add_p(im1,im2,imres);

  return(1);

}

/*******************************************
** --  imx_add_p() ----------------------------
**     Additive function
**     im1,im2 : im1 operation with im2
**     imres : Result image
********************************************/
int imx_add_p(grphic *im1, grphic *im2, grphic *imres)
{ 
  int i,j,wdth,hght;
  int err;
  float amp1,amp2,ampres;
  float amax;
  float rcoeff;
  int maxpixel;
  float icomp;


  hght=MINI(im1->height,im2->height);
  wdth=MINI(im1->width,im2->width);
/*    Recherche de amax   */
  amax=0;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++) {
      amp1=im1->rcoeff*im1->mri[i][j];
      amp2=im2->rcoeff*im2->mri[i][j];
      ampres=amp1+amp2;
      amax=MAXI(ampres,amax);
      }

/*   Calcul de maxpixel et icomp   */
  err=imx_brukermax(amax,&maxpixel,&icomp,&rcoeff);

/*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      {
      amp1=im1->rcoeff*im1->mri[i][j];
      amp2=im2->rcoeff*im2->mri[i][j];
      ampres=amp1+amp2;
      imres->mri[i][j]=ampres/rcoeff;
      }

  imres->max_pixel=maxpixel;
  imres->min_pixel=0;
  imres->cutoff_max=maxpixel;
  imres->cutoff_min=0;
  imres->icomp=icomp;
  imres->rcoeff=rcoeff;
  imres->height=hght;
  imres->width=wdth;
  return(1);

}

/*******************************************
** --  imx_sub() ----------------------------
**     Soustractive function
**     im_1,im_2 : im_1 operation with im_2
**     im_res : Result image
********************************************/
int imx_sub(int im_1, int im_2, int im_res)
{ 
  grphic *im1,*im2,*imres;

  im1=ptr_img(im_1);
  im2=ptr_img(im_2);
  imres=ptr_img(im_res);

  imx_sub_p(im1,im2,imres);

  if (MESG_DEBUG) printf(" fin sub\n");
  return(1);

}

/*******************************************
** --  imx_sub_p() ----------------------------
**     Soustractive function
**     im1,im2 : im1 operation with im2
**     imres : Result image
********************************************/
int imx_sub_p(grphic *im1, grphic *im2, grphic *imres)
{ 
  int i,j,wdth,hght;
  int err;
  float amp1,amp2,ampres;
  float amax,amin;
  float rcoeff;
  int maxpixel;
  float icomp;

  hght=MINI(im1->height,im2->height);
  wdth=MINI(im1->width,im2->width);

/*    Recherche de amax   */
  amax=0; amin=LONGMAX;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      {
      amp1=im1->rcoeff*im1->mri[i][j];
      amp2=im2->rcoeff*im2->mri[i][j];
      ampres=amp1-amp2;
      if (ampres>amax) amax=ampres;
      amin=MINI(ampres,amin);
      }
/*   Calcul de maxpixel et icomp   */
  amax=MAXI(abs(amax),abs(amin));
  err=imx_brukermax(amax,&maxpixel,&icomp,&rcoeff);
if (MESG_DEBUG) printf(" sub hgt wdt maxpixel1 rcoeff %d;%d;%d;%f\n",hght,wdth,maxpixel,rcoeff);

/*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      {
      amp1=im1->rcoeff*im1->mri[i][j];
      amp2=im2->rcoeff*im2->mri[i][j];
      ampres=amp1-amp2;
      imres->mri[i][j]=ampres/rcoeff;
      }

  imres->max_pixel=maxpixel;
  imres->min_pixel=0;
  imres->cutoff_max=maxpixel;
  imres->cutoff_min=0;
  imres->icomp=icomp;
  imres->rcoeff=rcoeff;
  imres->height=hght;
  imres->width=wdth;

  return(1);

}

/*******************************************
** --  imx_subabs() ----------------------------
**     Soustractive absolute function
**     im_1,im_2 : im_1 operation with im_2
**     im_res : Result image
********************************************/
int imx_subabs(int im_1, int im_2, int im_res)
{ 
  grphic *im1,*im2,*imres;
  int i,j,wdth,hght;
  int err;
  float amp1,amp2,ampres;
  float amax;
  float rcoeff;
  int maxpixel;
  float icomp;

  im1=ptr_img(im_1);
  im2=ptr_img(im_2);
  imres=ptr_img(im_res);

  hght=MINI(im1->height,im2->height);
  wdth=MINI(im1->width,im2->width);


/*    Recherche de amax   */
  amax=0;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      {
      amp1=im1->rcoeff*im1->mri[i][j];
      amp2=im2->rcoeff*im2->mri[i][j];
      ampres=fabs( (double) amp1-amp2);
      if (ampres>amax)
  amax=ampres;
      }
/*   Calcul de maxpixel et icomp   */
  err=imx_brukermax(amax,&maxpixel,&icomp,&rcoeff);
if (MESG_DEBUG) printf(" sub hgt wdt maxpixel1 rcoeff %d;%d;%d;%f\n",hght,wdth,maxpixel,rcoeff);

/*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      {
      amp1=im1->rcoeff*im1->mri[i][j];
      amp2=im2->rcoeff*im2->mri[i][j];
      ampres=fabs( (double) amp1-amp2);
      imres->mri[i][j]=ampres/rcoeff;
      }

  imres->max_pixel=maxpixel;
  imres->min_pixel=0;
  imres->cutoff_max=maxpixel;
  imres->cutoff_min=0;
  imres->icomp=icomp;
  imres->rcoeff=rcoeff;
  imres->height=hght;
  imres->width=wdth;
  if (MESG_DEBUG) printf(" fin sub\n");
  return(1);

}


/*******************************************
** --  imx_mul() ----------------------------
**     multiplicative function
**     im_1,im_2 : im_1 operation with im_2
**     im_res : Result image
********************************************/
int imx_mul(int im_1, int im_2, int im_res)
{ 
  grphic *im1,*im2,*imres;
  

  im1=ptr_img(im_1);
  im2=ptr_img(im_2);
  imres=ptr_img(im_res);

  imx_mul_p(im1,im2,imres);

  return(1);

}


/*******************************************
** --  imx_mul_p() ----------------------------
**     multiplicative function
**     im1,im2 : im1 operation with im2
**     imres : Result image
********************************************/

int imx_mul_p(grphic *im1, grphic *im2, grphic *imres)
{
  int i,j,wdth,hght;
  int err;
  float amp1,amp2,ampres;
  float amax;
  float rcoeff;
  int maxpixel;
  float icomp;

  hght=MINI(im1->height,im2->height);
  wdth=MINI(im1->width,im2->width);
/*    Recherche de amax   */
  amax=0;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      {
      amp1=im1->rcoeff*im1->mri[i][j];
      amp2=im2->rcoeff*im2->mri[i][j];
      ampres=amp1*amp2;
      if (ampres>amax)
  amax=ampres;
      }
/*   Calcul de maxpixel et icomp   */
  err=imx_brukermax(amax,&maxpixel,&icomp,&rcoeff);

/*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      {
      amp1=im1->rcoeff*im1->mri[i][j];
      amp2=im2->rcoeff*im2->mri[i][j];
      ampres=amp1*amp2;
      imres->mri[i][j]=ampres/rcoeff;
      }

  imres->max_pixel=maxpixel;
  imres->min_pixel=0;
  imres->cutoff_max=maxpixel;
  imres->cutoff_min=0;
  imres->icomp=icomp;
  imres->rcoeff=rcoeff;
  imres->height=hght;
  imres->width=wdth;
  return(1);


}

/*******************************************
** --  imx_div() ----------------------------
**     Divide function
**     im_1,im_2 : im_1 operation with im_2
**     im_res : Result image
********************************************/
int imx_div(int im_1, int im_2, int im_res)
{ 
  grphic *im1,*im2,*imres;

  im1=ptr_img(im_1);
  im2=ptr_img(im_2);
  imres=ptr_img(im_res);

  imx_div_p(im1,im2,imres);

  return(1);
}

/*******************************************
** --  imx_div_p() ----------------------------
**     Divide function
**     im1,im2 : im1 operation with im2
**     imres : Result image
********************************************/
int imx_div_p(grphic *im1, grphic *im2, grphic *imres)
{ 
  int i,j,wdth,hght;
  int err;
  float amp1,amp2,ampres;
  float amax;
  float rcoeff;
  int maxpixel;
  float icomp;

  hght=MINI(im1->height,im2->height);
  wdth=MINI(im1->width,im2->width);
/*    Recherche de amax   */
  amax=0;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      {
      amp1=im1->rcoeff*im1->mri[i][j];
      amp2=im2->rcoeff*im2->mri[i][j];
      if (amp2==0)
 ampres=0;
      else
 ampres=amp1/amp2;
      if (ampres>amax)
  amax=ampres;
      }
/*   Calcul de maxpixel et icomp   */
  err=imx_brukermax(amax,&maxpixel,&icomp,&rcoeff);

/*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      {
      amp1=im1->rcoeff*im1->mri[i][j];
      amp2=im2->rcoeff*im2->mri[i][j];
      if (amp2==0)
 ampres=0;
      else
 ampres=amp1/amp2;
      imres->mri[i][j]=ampres/rcoeff;
      }

  imres->max_pixel=maxpixel;
  imres->min_pixel=0;
  imres->cutoff_max=maxpixel;
  imres->cutoff_min=0;
  imres->icomp=icomp;
  imres->rcoeff=rcoeff;
  imres->height=hght;
  imres->width=wdth;
  return(1);
}

/*******************************************
** --  imx_abs() ----------------------------
**     Divide function
**     im_1,im_2 : im_1 operation with im_2
**     im_res : Result image
********************************************/
int imx_abs(int im_1, int im_res)
{ 
  grphic *im1,*imres;
  int i,j,wdth,hght;
  int err;
  float amp1,amax;
  float rcoeff;
  int maxpixel;
  float icomp;

  im1=ptr_img(im_1);
  imres=ptr_img(im_res);

  hght=im1->height;
  wdth=im1->width;
/*    Recherche de amax   */
  amax=0;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      {
      amp1=fabs( (double) im1->rcoeff*im1->mri[i][j]);
      if (amp1>amax)
  amax=amp1;
      }
/*   Calcul de maxpixel et icomp   */
  err=imx_brukermax(amax,&maxpixel,&icomp,&rcoeff);

/*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      {
      amp1=fabs( (double) im1->rcoeff*im1->mri[i][j]);
      imres->mri[i][j]=amp1/rcoeff;
      }

  imres->max_pixel=maxpixel;
  imres->min_pixel=0;
  imres->cutoff_max=maxpixel;
  imres->cutoff_min=0;
  imres->icomp=icomp;
  imres->rcoeff=rcoeff;
  imres->height=hght;
  imres->width=wdth;
  return(1);

}

/*******************************************
** --  imx_add_coe() ----------------------------
**     Additive function with coeff
**     im_1,coef : im_1 operation with coef
**     im_res : Result image
********************************************/
int imx_add_coe(int im_1, float coefficient, int im_res)
{ 
  grphic *im1,*imres;
  int i,j,wdth,hght;
  int err;
  float amp1,ampres;
  float amax;
  float rcoeff;
  int maxpixel;
  float icomp;

  im1=ptr_img(im_1);
  imres=ptr_img(im_res);

  hght=im1->height;
  wdth=im1->width;
/*    Recherche de amax   */
  amax=im1->rcoeff*im1->max_pixel+coefficient;
/*   Calcul de maxpixel et icomp   */
  err=imx_brukermax(amax,&maxpixel,&icomp,&rcoeff);

/*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      {
      amp1=im1->rcoeff*im1->mri[i][j];
      ampres=amp1+coefficient;
      imres->mri[i][j]=ampres/rcoeff;
      }

  imres->max_pixel=maxpixel;
  imres->min_pixel=0;
  imres->cutoff_max=maxpixel;
  imres->cutoff_min=0;
  imres->icomp=icomp;
  imres->rcoeff=rcoeff;
  imres->height=hght;
  imres->width=wdth;
  return(1);

}
/*******************************************
** --  imx_mul_coe() ----------------------------
**     Multiplicative function with coeff
**     im_1,coef : im_1 operation with coef
**     im_res : Result image
********************************************/
int imx_mul_coe(int im_1, float coefficient, int im_res)
{ 
  grphic *im1,*imres;
  im1=ptr_img(im_1);
  imres=ptr_img(im_res);
  
  imx_mul_coe_p(im1,coefficient,imres);
  return(1);

}
/********************************************/
int imx_mul_coe_p(grphic *im1, float coefficient, grphic *imres)
{
  int i,j,wdth,hght;
  int err;
  float amp1,ampres;
  float amax;
  float rcoeff;
  int maxpixel;
  float icomp;

  hght=im1->height;
  wdth=im1->width;
/*    Recherche de amax   */
  amax=im1->rcoeff*im1->max_pixel*coefficient;
/*   Calcul de maxpixel et icomp   */
  err=imx_brukermax(amax,&maxpixel,&icomp,&rcoeff);

/*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      {
      amp1=im1->rcoeff*im1->mri[i][j];
      ampres=amp1*coefficient;
      imres->mri[i][j]=ampres/rcoeff;
      }

  imres->max_pixel=maxpixel;
  imres->min_pixel=0;
  imres->cutoff_max=maxpixel;
  imres->cutoff_min=0;
  imres->icomp=icomp;
  imres->rcoeff=rcoeff;
  imres->height=hght;
  imres->width=wdth;

  return(1);
}


/*******************************************
** --  imx_sub_coe() ----------------------------
**     Soustractive function with coeff
**     im_1,coef : im_1 operation with coef
**     im_res : Result image
********************************************/
int imx_sub_coe(int im_1, float coefficient, int im_res)
{ 
  grphic *im1,*imres;
  int i,j,wdth,hght;
  int err;
  float amp1,ampres;
  float amax;
  float rcoeff;
  int maxpixel;
  float icomp;

  im1=ptr_img(im_1);
  imres=ptr_img(im_res);

  hght=im1->height;
  wdth=im1->width;
/*    Recherche de amax   */
  amax=im1->rcoeff*im1->max_pixel-coefficient;
/*   Calcul de maxpixel et icomp   */
  err=imx_brukermax(amax,&maxpixel,&icomp,&rcoeff);

/*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      {
      amp1=im1->rcoeff*im1->mri[i][j];
      ampres=amp1-coefficient;
      imres->mri[i][j]=ampres/rcoeff;
      }

  imres->max_pixel=maxpixel;
  imres->min_pixel=0;
  imres->cutoff_max=maxpixel;
  imres->cutoff_min=0;
  imres->icomp=icomp;
  imres->rcoeff=rcoeff;
  imres->height=hght;
  imres->width=wdth;
  return(1);

}

/*******************************************
** --  imx_div_coe() ----------------------------
**     Additive function with coeff
**     im_1,coef : im_1 operation with coef
**     im_res : Result image
********************************************/
int imx_div_coe(int im_1, float coefficient, int im_res)
{ 
  grphic *im1,*imres;
  int i,j,wdth,hght;
  int err;
  float amp1,ampres;
  float amax;
  float rcoeff;
  int maxpixel;
  float icomp;

  im1=ptr_img(im_1);
  imres=ptr_img(im_res);

  hght=im1->height;
  wdth=im1->width;
/*    Recherche de amax   */
  if (coefficient!=0)
    amax=im1->rcoeff*im1->max_pixel/coefficient;
  else
    return(0);
/*   Calcul de maxpixel et icomp   */
  err=imx_brukermax(amax,&maxpixel,&icomp,&rcoeff);

/*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      {
      amp1=im1->rcoeff*im1->mri[i][j];
      ampres=amp1/coefficient;
      imres->mri[i][j]=ampres/rcoeff;
      }

  imres->max_pixel=maxpixel;
  imres->min_pixel=0;
  imres->cutoff_max=maxpixel;
  imres->cutoff_min=0;
  imres->icomp=icomp;
  imres->rcoeff=rcoeff;
  imres->height=hght;
  imres->width=wdth;
  return(1);

}


/*******************************************
** --  imx_and() ----------------------------
**     Logical AND function
**     im_1,im_2 : im_1 operation with im_2
**     im_res : Result image
********************************************/
int imx_and(int im_1, int im_2, int im_res)
{ 
  grphic *im1,*im2,*imres;
  int i,j,wdth,hght;
  int err;
  float amp1,ampres;
  float amax;
  float rcoeff;
  int maxpixel;
  float icomp;

  im1=ptr_img(im_1);
  im2=ptr_img(im_2);
  imres=ptr_img(im_res);

  hght=MINI(im1->height,im2->height);
  wdth=MINI(im1->width,im2->width);
/*    Recherche de amax   */
  amax=0.;
  ampres=0.;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      {
      amp1=(im1->rcoeff)*(im1->mri[i][j]);
      if (im1->mri[i][j]&&im2->mri[i][j])ampres=amp1;
      if (ampres>amax)
  amax=ampres;
      }
/*   Calcul de maxpixel et icomp   */
  err=imx_brukermax(amax,&maxpixel,&icomp,&rcoeff);
/*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      {
      amp1=im1->rcoeff*im1->mri[i][j];
      ampres=0;
      if (im1->mri[i][j]&&im2->mri[i][j])ampres=amp1;
      imres->mri[i][j]=ampres/rcoeff;
      }

  imres->max_pixel=maxpixel;
  imres->min_pixel=0;
  imres->cutoff_max=maxpixel;
  imres->cutoff_min=0;
  imres->icomp=icomp;
  imres->rcoeff=rcoeff;
  imres->height=hght;
  imres->width=wdth;
  return(1);

}

/*******************************************
** --  imx_or() ----------------------------
**     Logical OR function
**     im_1,im_2 : im_1 operation with im_2
**     im_res : Result image
********************************************/
int imx_or(int im_1, int im_2, int im_res)
{ 
  grphic *im1,*im2,*imres;
  int i,j,wdth,hght;
  int err;
  float amp1,amp2,ampres;
  float amax;
  float rcoeff;
  int maxpixel;
  float icomp;

  im1=ptr_img(im_1);
  im2=ptr_img(im_2);
  imres=ptr_img(im_res);

  hght=MINI(im1->height,im2->height);
  wdth=MINI(im1->width,im2->width);
/*    Recherche de amax   */
  amax=0;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      {
      amp1=im1->rcoeff*im1->mri[i][j];
      amp2=im2->rcoeff*im2->mri[i][j];
      ampres=amp1;
      if (im1->mri[i][j]==0)ampres=amp2;
      if (ampres>amax)
  amax=ampres;
      }
/*   Calcul de maxpixel et icomp   */
  err=imx_brukermax(amax,&maxpixel,&icomp,&rcoeff);

/*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      {
      amp1=im1->rcoeff*im1->mri[i][j];
      amp2=im2->rcoeff*im2->mri[i][j];
      ampres=amp1;
      if (im1->mri[i][j]==0)ampres=amp2;
      imres->mri[i][j]=ampres/rcoeff;
      }

  imres->max_pixel=maxpixel;
  imres->min_pixel=0;
  imres->cutoff_max=maxpixel;
  imres->cutoff_min=0;
  imres->icomp=icomp;
  imres->rcoeff=rcoeff;
  imres->height=hght;
  imres->width=wdth;
  return(1);

}

/*******************************************
** --  imx_eg() ----------------------------
**     Logical EGALITE function
**     im_1,im_2 : im_1 operation with im_2
**     im_res : Result image
********************************************/
int imx_eg(int im_1, int im_2, int im_res)
{ 
  grphic *im1,*im2,*imres;
  int i,j,wdth,hght;
  int err;
  float amp1,amp2,ampres;
  float amax;
  float rcoeff;
  int maxpixel;
  float icomp;

  im1=ptr_img(im_1);
  im2=ptr_img(im_2);
  imres=ptr_img(im_res);

  hght=MINI(im1->height,im2->height);
  wdth=MINI(im1->width,im2->width);
/*    Recherche de amax   */
  amax=0.;
  ampres=0.;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      {
      amp1=im1->rcoeff*im1->mri[i][j];
      amp2=im2->rcoeff*im2->mri[i][j];
      if (amp1==amp2)ampres=amp1;
      if (ampres>amax)
  amax=ampres;
      }
/*   Calcul de maxpixel et icomp   */
  err=imx_brukermax(amax,&maxpixel,&icomp,&rcoeff);

/*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      {
      amp1=im1->rcoeff*im1->mri[i][j];
      amp2=im2->rcoeff*im2->mri[i][j];
      ampres=0.;
      if (amp1==amp2)ampres=amp1;
      imres->mri[i][j]=ampres/rcoeff;
      }

  imres->max_pixel=maxpixel;
  imres->min_pixel=0;
  imres->cutoff_max=maxpixel;
  imres->cutoff_min=0;
  imres->icomp=icomp;
  imres->rcoeff=rcoeff;
  imres->height=hght;
  imres->width=wdth;
  return(1);

}

/*******************************************
** --  imx_inv() ----------------------------
**     Logical INVERSION function
**     Soustraction de l'image de son maximum
**     im_1 : Inversion de im_1
**     im_res : Result image
********************************************/
int imx_inv(int im_1, int im_res)
{ 
  grphic *im1,*imres;
  int i,j,wdth,hght;
  int err;
  float ampres;
  float amax;
  float rcoeff;
  int maxpixel;
  float icomp;

  im1=ptr_img(im_1);
  imres=ptr_img(im_res);

  hght=im1->height;
  wdth=im1->width;
/*    Recherche de amax   */
  amax=0;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      {
      ampres=im1->rcoeff*(im1->max_pixel-im1->mri[i][j]);
      if (ampres>amax)
  amax=ampres;
      }
/*   Calcul de maxpixel et icomp   */
  err=imx_brukermax(amax,&maxpixel,&icomp,&rcoeff);

/*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      {
      ampres=im1->rcoeff*(im1->max_pixel-im1->mri[i][j]);
      imres->mri[i][j]=ampres/rcoeff;
      }

  imres->max_pixel=maxpixel;
  imres->min_pixel=0;
  imres->cutoff_max=maxpixel;
  imres->cutoff_min=0;
  imres->icomp=icomp;
  imres->rcoeff=rcoeff;
  imres->height=hght;
  imres->width=wdth;
  return(1);

}

