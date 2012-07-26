/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*!*************************************************************************
***	
***	\file:		oper_3d.c
***
***	project:	Imagix 1.01
***			
***
***	\brief description:	Operation mathematique et logique sur image 3D
***	
***	
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***     Last user action by: Mr. ARMSPACH on Dec 15th 1996
***
**************************************************************************/
#include <config.h>
#include <stdio.h>
#include <stdlib.h> 	/* Needed when you use GET_FLOAT	*/
#include <math.h> 	/* Needed when you use math function	*/
#include <string.h> 	/* Needed for strcpy	*/
#include <float.h> 

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/gui/imx_picture3d.h"
#include "noyau/imx_3d.h"
#include "noyau/imx_lang.h"
#include "noyau/io/imx_file.h"
#include "math/oper_3d.h"


/*******************************************
 ** --  op_math_3d() ----------------------------
 **
 ** Op�ations math�atiques
 ********************************************/

#ifndef __DOXYGEN_SHOULD_SKIP_THIS__
int imx_add_3d (int im_1, int im_2, int im_res);
int imx_sub_3d (int im_1, int im_2, int im_res);
int imx_mul_3d (int im_1, int im_2, int im_res);
int imx_div_3d (int im_1, int im_2, int im_res);
int imx_sqr_3d (int im_1, int im_res);
int imx_abs_3d (int im_1, int im_res);

void imx_add_coe_3d (int im_1, float coefficient, int im_res);
void imx_sub_coe_3d (int im_1, float coefficient, int im_res);
void imx_mul_coe_3d (int im_1, float coefficient, int im_res);
void imx_div_coe_3d (int im_1, float coefficient, int im_res);

int imx_and_3d (int im_1, int im_2, int im_res);
int imx_or_3d (int im_1, int im_2, int im_res);
int imx_eg_3d (int im_1, int im_2, int im_res);
int imx_inv_3d (int im_1, int im_res);

int imx_add_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
int imx_sub_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
int imx_subabs_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
int imx_mul_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
int imx_div_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
int imx_sqr_3d_p (grphic3d *im1, grphic3d *imres);
int imx_abs_3d_p (grphic3d *im1, grphic3d *imres);

int imx_copie_param_3d_p (grphic3d *im1, grphic3d *im2);
int imx_mul_coe_3d_p (grphic3d *im1, float coefficient, grphic3d *imres);
int imx_div_coe_3d_p (grphic3d *im1, float coefficient, grphic3d *imres);
int imx_add_coe_3d_p (grphic3d *im1, float coefficient, grphic3d *imres);
int imx_sub_coe_3d_p (grphic3d *im1, float coefficient, grphic3d *imres);

int imx_and_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
int imx_or_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
int imx_inv_3d_p (grphic3d *im1, grphic3d *imres);
int imx_eg_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);

int imx_inimaxpixel_3d_p (grphic3d *im1);
int imx_inimaxminpixel_3d_p (grphic3d *im1);
int imx_iniminpixel_3d_p (grphic3d *im1);
int imx_iniparaimg_3d_p (grphic3d *im1);
#endif /*__DOXYGEN_SHOULD_SKIP_THIS__*/

/*****************************************/
/*! \brief fonction du menu pour l'appel des operations math

permet l'addition, la soustraction, la multiplication etc... de 2 images
******************************************/
void	op_math_3d(void)
{ 
  char *quest[11];
  int answer_nr;
  int i;
  int im_1,im_2,im_res;

  /* Question concerning type of operation  */
  /* and image					*/ 
  
  for(i=0;i<11;i++)
    quest[i]=CALLOC(80,char);
    
  strcpy(quest[0],TEXT0075);
  strcpy(quest[1],TEXT0076);
  strcpy(quest[2],TEXT0077);
  strcpy(quest[3],TEXT0078);
  strcpy(quest[4],TEXT0079);
  strcpy(quest[5],"Racine carree");
  strcpy(quest[6],"Log");
  strcpy(quest[7],"Exponent");
  strcpy(quest[8],"Max");
  strcpy(quest[9],"Min");
  strcpy(quest[10],"\0");
  /* Put a question box at position 150,200. */
  answer_nr=GETV_QCM(TEXT0043,(const char **)quest);

  for(i=0;i<11;i++)
    FREE(quest[i]);

  /*     Recherche numero image      */
  if ((answer_nr<4)||(answer_nr>7)) {
    im_1=GET_PLACE3D("Image 1:");
    im_2=GET_PLACE3D("Operation with image 2:");
    im_res=GET_PLACE3D(TEXT0006);
  }
  else {
    im_1=GET_PLACE3D("Image 1:");
    im_2=0;
    im_res=GET_PLACE3D(TEXT0006);
  }

  switch ((int)answer_nr) {
  
  case 0: 	/*  Addition */
		
    imx_add_3d(im_1,im_2,im_res);	
    show_picture_3d(im_res);
    break;

  case 1:		/*  Soustraction */

    imx_sub_3d(im_1,im_2,im_res);
    show_picture_3d(im_res);
    break;
	
  case 2:		/*  Multiplication */

    imx_mul_3d(im_1,im_2,im_res);
    show_picture_3d(im_res);
    break;
	
  case 3:		/*  Division */

    imx_div_3d(im_1,im_2,im_res);
    show_picture_3d(im_res);
    break;
	
  case 4:		/*  Absolu */

    imx_abs_3d(im_1,im_res);
    show_picture_3d(im_res);
    break;
  case 5:		/*  racine carree */

    imx_sqr_3d(im_1,im_res);
    show_picture_3d(im_res);
    break;
	
  case 6:		/*  Logarithme */

    imx_log_3d(im_1,im_res);
    show_picture_3d(im_res);
    break;
  
  case 7:		/*  Exponent */

    imx_exponent_3d(im_1,im_res);
    show_picture_3d(im_res);
    break;  
  
  case 8:		/*  Max */

    imx_max_3d(im_1,im_2,im_res);
    show_picture_3d(im_res);
    break;
  
  case 9:		/*  Min */

    imx_min_3d(im_1,im_2,im_res);
    show_picture_3d(im_res);
    break;
  
  
  default:
    fprintf(stderr,"Warning: in Mathematical Operation\n");
    break;

  }

  return;
}

/*******************************************
 ** --  op_math_coe_3d() ----------------------------
 **
 **    Operations mathematiques avec coeeficient
 ********************************************/
void	op_math_coe_3d(void)
{ 
  char *quest[7];
  int answer_nr;
  int i,err=0;
  int ii,jj,kk;
  int im_1,im_res;
  float coefficient;
  grphic3d *im1;


   

  /* Question concerning type of operation  */
  /* and image					*/ 
  
	quest[0]="Addition avec coeff";
	quest[1]="Soustraction avec coeff";
	quest[2]="Multiplication par coeff";
	quest[3]="Division par coeff";
	quest[4]="Soustraction de la plus petite valeure de l image a l image";
	quest[5]="\0";
	quest[6]=NULL;

  /* Put a question box at position 150,200. */
  answer_nr=GETV_QCM("Votre reponse",(const char **)quest);

  /*     Image   				*/
  im_1 = GET_PLACE3D("Image 1:");
  if ((int)answer_nr!=4)
	 coefficient = GET_FLOAT("Coefficient:", 0, &err);
	 
  im_res=GET_PLACE3D(TEXT0006);

  switch ((int)answer_nr) {
  
  case 0: 	/*  Addition */
		
    imx_add_coe_3d(im_1,coefficient,im_res);	
    show_picture_3d(im_res);
    break;

  case 1:		/*  Soustraction */

    imx_sub_coe_3d(im_1,coefficient,im_res);
    show_picture_3d(im_res);
    break;
	
  case 2:		/*  Multiplication */

    imx_mul_coe_3d(im_1,coefficient,im_res);
    show_picture_3d(im_res);
    break;
	
  case 3:		/*  Division */

    if (coefficient==0) 
    {
	PUT_ERROR("Impossible de diviser par 0 !");
	return;
    }

    imx_div_coe_3d(im_1,coefficient,im_res);
    show_picture_3d(im_res);
    break;
	
  case 4:		/*  Translation de maniere a ce que la plus petite valeur de l image soit mise a 0 */
	coefficient = FLT_MAX;
	im1=ptr_img_3d(im_1);
	for (ii = 0; ii < im1->width ; ii++)
    for (jj = 0; jj < im1->height; jj++)
    for (kk = 0; kk < im1->depth ; kk++)
		if (coefficient > (im1->mri[ii][jj][kk] * im1->rcoeff))
			{
			coefficient = im1->mri[ii][jj][kk] * im1->rcoeff;
			}
	imx_sub_coe_3d(im_1,coefficient,im_res);
    show_picture_3d(im_res);
    break;

  default:
    fprintf(stderr,"Warning: in Mathematical Operation with coef\n");
    break;

  }

  return;
}

/*******************************************
 ** --  op_log_3d() ----------------------------
 */
/*!
**    Operations logiques sue les images
**
**
********************************************/
void	op_log_3d(void)
{ 
  char *quest[6];
  int answer_nr;
  int i;
  int im_1=0,im_2=0,im_res=0;

  /* Question concerning type of operation  */
  /* and image					*/ 
  
  for(i=0;i<6;i++)
    quest[i]=CALLOC(80,char);
  strcpy(quest[0],"And");
  strcpy(quest[1],"Or");
  strcpy(quest[2],"Egalite");
  strcpy(quest[3],"Inversion");
  strcpy(quest[4],"And (2 images)");
  strcpy(quest[5],"\0");
  answer_nr=GETV_QCM(TEXT0043,(const char **)quest);

  for(i=0;i<6;i++)
    FREE(quest[i]);

  /*     Image   				*/
  im_1=GET_PLACE3D("Image 1:");
  im_2=0;
  if (answer_nr!=3)
    im_2=GET_PLACE3D("Image 2:");
  if (answer_nr!=4)
    im_res=GET_PLACE3D(TEXT0006);

  switch ((int)answer_nr) {
 
  case 0: 	/*  ET logique */
		
    imx_and_3d(im_1,im_2,im_res);	
    show_picture_3d(im_res);
    break;

  case 1:		/*  OU logique */

    imx_or_3d(im_1,im_2,im_res);
    show_picture_3d(im_res);
    break;
	
  case 2:		/*  Egalite */

    imx_eg_3d(im_1,im_2,im_res);
    show_picture_3d(im_res);
    break;
	
  case 3:		/*  Inversion */

    imx_inv_3d(im_1,im_res);
    show_picture_3d(im_res);
    break;

  case 4:		/*  Et logique sur 2 images */

    imx_and2img_3d(im_1,im_2);
    show_picture_3d(im_1);
    show_picture_3d(im_2);
    break;
	
  default:
    PUT_WARN("Warning: in Mathematical Operation\n");
    break;

  }

  return;
}

/*******************************************
 ** --  imx_add_3d() ----------------------------
 */
/*!    \brief Addition function
**     \param im_1,im_2 : im_1 operation with im_2
**     \param im_res : Result image
**	   \retval 1
********************************************/
int	imx_add_3d(int im_1, int im_2, int im_res)
{
  grphic3d *im1,*im2,*imres;
   
  im1=ptr_img_3d(im_1);
  im2=ptr_img_3d(im_2);
  imres=ptr_img_3d(im_res);

  imx_add_3d_p(im1,im2,imres);
  
  return(1);
}        

/*******************************************
 ** --  imx_sub_3d() ----------------------------
 */
/*!    \brief Soustractive function
**     \param im_1,im_2 : im_1 operation with im_2
**     \param im_res : Result image
**	   \retval 1
********************************************/
int	imx_sub_3d(int im_1, int im_2, int im_res)
{
  grphic3d *im1,*im2,*imres;
   
  im1=ptr_img_3d(im_1);
  im2=ptr_img_3d(im_2);
  imres=ptr_img_3d(im_res);

  imx_sub_3d_p(im1,im2,imres);
  
  return(1);
}        

/*******************************************
 ** --  imx_subabs_3d() ----------------------------
 **     Absolute soustractive function
 **     im_1,im_2 : im_1 operation with im_2
 **     im_res : Result image
 ********************************************/
int	imx_subabs_3d(int im_1, int im_2, int im_res)
{
  grphic3d *im1,*im2,*imres;
   
  im1=ptr_img_3d(im_1);
  im2=ptr_img_3d(im_2);
  imres=ptr_img_3d(im_res);

  imx_subabs_3d_p(im1,im2,imres);
  
  return(1);
}    


/*******************************************
 ** --  imx_mul_3d() ----------------------------
 */
/*!    \brief multiplicative function
**     \param im_1,im_2 : im_1 operation with im_2
**     \param im_res : Result image 
** 	   \retval 1
********************************************/
int	imx_mul_3d(int im_1, int im_2, int im_res)
{
  grphic3d *im1,*im2,*imres;
   
  im1=ptr_img_3d(im_1);
  im2=ptr_img_3d(im_2);
  imres=ptr_img_3d(im_res);

  imx_mul_3d_p(im1,im2,imres);
  
  return(1);
}        

/*******************************************
 ** --  imx_div_3d() ----------------------------
 */
/*!    \brief Divide function
**     \param im_1,im_2 : im_1 operation with im_2
**     \param im_res : Result image
**	   \retval 1
********************************************/
int	imx_div_3d(int im_1, int im_2, int im_res)
{
  grphic3d *im1,*im2,*imres;
   
  im1=ptr_img_3d(im_1);
  im2=ptr_img_3d(im_2);
  imres=ptr_img_3d(im_res);

  imx_div_3d_p(im1,im2,imres);
  
  return(1);
}        

/*******************************************
 ** --  imx_abs_3d() ----------------------------
 */
/*!    \brief Absolute function
**     \param im_1,im_2 : im_1 operation with im_2
**     \param im_res : Result image
**	   \retval 1
********************************************/
int	imx_abs_3d(int im_1, int im_res)
{
  grphic3d *im1,*imres;
   
  im1=ptr_img_3d(im_1);
  imres=ptr_img_3d(im_res);

  imx_abs_3d_p(im1,imres);
  
  return(1);
}       
 
/*******************************************
 ** --  imx_sqr_3d() ----------------------------
 */
/*!    \brief Root square function
**     \param im_1,im_2 : im_1 operation with im_2
**     \param im_res : Result image
**	   \retval 1
********************************************/
int	imx_sqr_3d(int im_1, int im_res)
{
  grphic3d *im1,*imres;
   
  im1=ptr_img_3d(im_1);
  imres=ptr_img_3d(im_res);

  imx_sqr_3d_p(im1,imres);
  
  return(1);
}  

/*******************************************
 ** --  imx_log_3d() ----------------------------
 */
/*!    \brief Logarithm function
**     \param im_1,im_2 : im_1 operation with im_2
**     \param im_res : Result image
**	   \retval 1
********************************************/
int	imx_log_3d(int im_1, int im_res)
{
  grphic3d *im1,*imres;
   
  im1=ptr_img_3d(im_1);
  imres=ptr_img_3d(im_res);

  imx_log_3d_p(im1,imres);
  
  return(1);
}
/*********************************************/  
int	imx_exponent_3d(int im_1, int im_res)
{
  grphic3d *im1,*imres;
   
  im1=ptr_img_3d(im_1);
  imres=ptr_img_3d(im_res);

  imx_exponent_3d_p(im1,imres);
  
  return(1);
}  
/*******************************************
 ** --  imx_max_3d() ----------------------------
 */
/*!    \brief Max function
**     \param im_1,im_2 : im_1 operation with im_2
**     \param im_res : Result image
**	   \retval 1
********************************************/
int	imx_max_3d(int im_1, int im_2, int im_res)
{
  grphic3d *im1,*im2,*imres;
   
  im1=ptr_img_3d(im_1);
  im2=ptr_img_3d(im_2);
  imres=ptr_img_3d(im_res);

  imx_max_3d_p(im1,im2,imres);
  
  return(1);
}        


/*******************************************
 ** --  imx_min_3d() ----------------------------
 */
/*!    \brief Min function
**     \param im_1,im_2 : im_1 operation with im_2
**     \param im_res : Result image
**	   \retval 1
********************************************/
int	imx_min_3d(int im_1, int im_2, int im_res)
{
  grphic3d *im1,*im2,*imres;
   
  im1=ptr_img_3d(im_1);
  im2=ptr_img_3d(im_2);
  imres=ptr_img_3d(im_res);

  imx_min_3d_p(im1,im2,imres);
  
  return(1);
}        


/*******************************************
 ** --  imx_add_3d_p() ----------------------------
 */
/*!    \ingroup ToolBoxImage3d 
**
**     \brief Additive function
**
**     \param im_1,im_2 : im_1 operation with im_2
**     \param im_res : Result image (E/S)
**	   \retval 1
********************************************/
int	imx_add_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{ 
  int i,j,k,wdth,hght,dpth;
  int err;
  float amp1,amp2,ampres;
  float amax,amin,rcoeff;
  float im1rcoeff,im2rcoeff;

  im1rcoeff=im1->rcoeff;
  im2rcoeff=im2->rcoeff;

  hght=MINI(im1->height,im2->height);
  wdth=MINI(im1->width,im2->width);
  dpth=MINI(im1->depth,im2->depth);

  imx_copie_param_3d_p(im1,imres);
  /*    Recherche de amax   */
  amax=0; amin=LONGMAX;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=im1rcoeff*im1->mri[i][j][k];
	  amp2=im2rcoeff*im2->mri[i][j][k];
	  ampres=amp1+amp2;
	  if (ampres>amax)
	    amax=ampres;
	  amin=MINI(ampres,amin);
	}
  /*   Calcul de maxpixel et icomp   */
  err=imx_brukermax_3d(amax,amin,imres); /* Attention modif im2 si =imres */
  rcoeff=imres->rcoeff;

  /*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=im1rcoeff*im1->mri[i][j][k];
	  amp2=im2rcoeff*im2->mri[i][j][k];
	  ampres=amp1+amp2;
	  imres->mri[i][j][k]=(TYPEMRI3D)floor(ampres/rcoeff+0.5);
	}

  return(1);

}

/*******************************************
 ** --  imx_sub_3d_p() ----------------------------
 */
/*!    \ingroup ToolBoxImage3d 
**
**     \brief Soustractive function function
**     \param im_1,im_2 : im_1 operation with im_2
**     \param im_res : Result image (E/S)
**	   \retval 1
********************************************/
int	imx_sub_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{ 
  int i,j,k,wdth,hght,dpth;
  int err;
  float amp1,amp2,ampres;
  float amax,amin;
  float rcoeff;
  float im1rcoeff,im2rcoeff;

  im1rcoeff=im1->rcoeff;
  im2rcoeff=im2->rcoeff;
  hght=MINI(im1->height,im2->height);
  wdth=MINI(im1->width,im2->width);
  dpth=MINI(im1->depth,im2->depth);

  imx_copie_param_3d_p(im1,imres);

  /*    Recherche de amax   */
  amax=0; amin=LONGMAX;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=im1rcoeff*im1->mri[i][j][k];
	  amp2=im2rcoeff*im2->mri[i][j][k];
	  ampres=amp1-amp2;
	  if (ampres>amax)
	    amax=ampres;
	  amin=MINI(ampres,amin);
	}
  /*   Calcul de maxpixel et icomp   */
  err=imx_brukermax_3d(amax,amin,imres);
  rcoeff=imres->rcoeff;

  /*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=im1rcoeff*im1->mri[i][j][k];
	  amp2=im2rcoeff*im2->mri[i][j][k];
	  ampres=amp1-amp2;
	  imres->mri[i][j][k]=(TYPEMRI3D)floor(ampres/rcoeff+0.5);
	}

  if (MESG3D_DEBUG) printf(" fin sub\n");
  return(1);

}
/*******************************************
 ** --  imx_subabs_3d_p() ----------------------------
 **     Absolute soustractive function
 **     im_1,im_2 : im_1 operation with im_2
 **     im_res : Result image
 ********************************************/
int	imx_subabs_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{ 
  int i,j,k,wdth,hght,dpth;
  int err;
  float amp1,amp2,ampres;
  float amax,amin;
  float rcoeff;
  float im1rcoeff,im2rcoeff;

  im1rcoeff=im1->rcoeff;
  im2rcoeff=im2->rcoeff;
  hght=MINI(im1->height,im2->height);
  wdth=MINI(im1->width,im2->width);
  dpth=MINI(im1->depth,im2->depth);

  copie_image_param_adapt_buffer(im1, imres);

  /*    Recherche de amax   */
  amax=0; amin=LONGMAX;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=im1rcoeff*im1->mri[i][j][k];
	  amp2=im2rcoeff*im2->mri[i][j][k];
	  ampres=(float)fabs ( (double) amp1-amp2);
	  if (ampres>amax)
	    amax=ampres;
	  amin=MINI(amin,ampres);
	}
  /*   Calcul de maxpixel et icomp   */
  err=imx_brukermax_3d(amax,amin,imres);
  rcoeff=imres->rcoeff;

  /*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=im1rcoeff*im1->mri[i][j][k];
	  amp2=im2rcoeff*im2->mri[i][j][k];
	  ampres=(float)fabs ( (double) amp1-amp2);
	  imres->mri[i][j][k]=(TYPEMRI3D)floor(ampres/rcoeff+0.5);
	}

  if (MESG3D_DEBUG) printf(" fin sub\n");
  return(1);

}

/*******************************************
 ** --  imx_mul_3d_p() ----------------------------
 */
/*!  \ingroup ToolBoxImage3d 
**
**   \brief   multiplicative function
**   \param  im_1,im_2 : im_1 operation with im_2
**   \param  im_res : Result image (E/S)
**   \retval 1
********************************************/
int	imx_mul_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{ 
  int i,j,k,wdth,hght,dpth;
  int err;
  float amp1,amp2,ampres;
  float amax,amin;
  float rcoeff;
  float im1rcoeff,im2rcoeff;

  hght=MINI(im1->height,im2->height);
  wdth=MINI(im1->width,im2->width);
  dpth=MINI(im1->depth,im2->depth);
  im1rcoeff=im1->rcoeff;
  im2rcoeff=im2->rcoeff;

  imx_copie_param_3d_p(im1,imres);

  /*    Recherche de amax   */
  amax=0; amin=LONGMAX;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=im1rcoeff*im1->mri[i][j][k];
	  amp2=im2rcoeff*im2->mri[i][j][k];
	  ampres=amp1*amp2;
	  if (ampres>amax)
	    amax=ampres;
	  amin=MINI(amin,ampres);
	}
  /*   Calcul de maxpixel et icomp   */
  err=imx_brukermax_3d(amax,amin,imres);
  rcoeff=imres->rcoeff;

  /*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=im1rcoeff*im1->mri[i][j][k];
	  amp2=im2rcoeff*im2->mri[i][j][k];
	  ampres=amp1*amp2;
	  imres->mri[i][j][k]=(TYPEMRI3D)floor(ampres/rcoeff+0.5);
	}

  return(1);

}

/*******************************************
 ** --  imx_div_3d_p() ----------------------------
 */
/*!    \ingroup ToolBoxImage3d 
**
**     \brief Divide function
**     \param im_1,im_2 : im_1 operation with im_2
**     \param im_res : Result image (E/S)
**	   \retval 1
********************************************/
int	imx_div_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{ 
  int i,j,k,wdth,hght,dpth;
  int err;
  float amp1,amp2,ampres;
  float amax,amin;
  float rcoeff;
  float im1rcoeff,im2rcoeff;

  hght=MINI(im1->height,im2->height);
  wdth=MINI(im1->width,im2->width);
  dpth=MINI(im1->depth,im2->depth);
  im1rcoeff=im1->rcoeff;
  im2rcoeff=im2->rcoeff;

  imx_copie_param_3d_p(im1,imres);

  /*    Recherche de amax   */
  amax=0; amin=LONGMAX;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=im1rcoeff*im1->mri[i][j][k];
	  amp2=im2rcoeff*im2->mri[i][j][k];
	  if (amp2==0)
	    ampres=0;
	  else
	    ampres=amp1/amp2;
	  if (ampres>amax)
	    amax=ampres;
	  amin=MINI(amin,ampres);
	}
  /*   Calcul de maxpixel et icomp   */
  err=imx_brukermax_3d(amax,amin,imres);
  rcoeff=imres->rcoeff;

  /*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=im1rcoeff*im1->mri[i][j][k];
	  amp2=im2rcoeff*im2->mri[i][j][k];
	  if (amp2==0)
	    ampres=0;
	  else
	    ampres=amp1/amp2;
	  imres->mri[i][j][k]=(TYPEMRI3D)floor(ampres/rcoeff+0.5);
	}

  return(1);

}

/*******************************************
 ** --  imx_abs_3d_p() ----------------------------
 */
/*!    \ingroup ToolBoxImage3d 
**
**     \brief Absolute function
**     \param im_1 : image
**     \param im_res : Result image (E/S)
**	   \retval 1
********************************************/
int	imx_abs_3d_p(grphic3d *im1, grphic3d *imres)
{ 
  int i,j,k,wdth,hght,dpth;
  int err;
  float amp1,amax,amin;
  double rcoeff,im1rcoeff;

  hght=im1->height;
  wdth=im1->width;
  dpth=im1->depth;
  im1rcoeff= (double) im1->rcoeff;

  imx_copie_param_3d_p(im1,imres);

  /*    Recherche de amax   */
  amax=0; amin=LONGMAX;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=(float)fabs( im1rcoeff * (double) im1->mri[i][j][k]);
	  if (amp1>amax)
	    amax=amp1;
	  amin=MINI(amin,amp1);
	}
  /*   Calcul de maxpixel et icomp   */
  err=imx_brukermax_3d(amax,amin,imres);
  rcoeff=(double) imres->rcoeff;

  /*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=(float)fabs( im1rcoeff* (double) im1->mri[i][j][k]);
	  imres->mri[i][j][k]=(TYPEMRI3D)floor(amp1/rcoeff+0.5);
	}

  return(1);

}

/*******************************************
 ** --  imx_sqr_3d_p() ----------------------------
 */
/*!    \ingroup ToolBoxImage3d 
**
**     \brief Square function
**     \param im_1 : image
**     \param im_res : Result image (E/S)
**	   \retval 1
********************************************/
int	imx_sqr_3d_p(grphic3d *im1, grphic3d *imres)
{ 
  int i,j,k,wdth,hght,dpth;
  int err;
  float amp1,amax,amin;
  double rcoeff,im1rcoeff;

  hght=im1->height;
  wdth=im1->width;
  dpth=im1->depth;
  im1rcoeff= (double) im1->rcoeff;

  imx_copie_param_3d_p(im1,imres);

  /*    Recherche de amax   */
  amax=0; amin=LONGMAX;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=(float)(1.0*sqrt( im1rcoeff * (double) im1->mri[i][j][k]>0?im1rcoeff * (double) im1->mri[i][j][k]:0));
	  if (amp1>amax)
	    amax=amp1;
	  amin=MINI(amin,amp1);
	}
  /*   Calcul de maxpixel et icomp   */
  err=imx_brukermax_3d(amax,amin,imres);
  rcoeff=(double) imres->rcoeff;

  /*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=(float)(1.0*sqrt( im1rcoeff * (double) im1->mri[i][j][k]>0?im1rcoeff * (double) im1->mri[i][j][k]:0));
	  imres->mri[i][j][k]=(TYPEMRI3D)floor(amp1/rcoeff+0.5);
	}

  return(1);

}


/*******************************************
 ** --  imx_log_3d_p() ----------------------------
 */
/*!    \ingroup ToolBoxImage3d 
**
**     \brief Log function
**     \param im_1 : image
**     \param im_res : Result image (E/S)
**	   \retval 1
********************************************/
int imx_log_3d_p(grphic3d *im1, grphic3d *imres)
{ 
  int i,j,k,wdth,hght,dpth;
  int err;
  float amp1,amax,amin;
  double rcoeff,im1rcoeff;

  hght=im1->height;
  wdth=im1->width;
  dpth=im1->depth;
  im1rcoeff= (double) im1->rcoeff;

  imx_copie_param_3d_p(im1,imres);

  /*    Recherche de amax   */
  amax=-1.0*LONGMAX; amin=LONGMAX;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
      	if (im1->mri[i][j][k]>0)
	  {
	    amp1=(float)(1.0*log( im1rcoeff * (double) im1->mri[i][j][k]));
	    if (amp1>amax)
	      amax=amp1;
	    if (amp1<amin)
	      amin=amp1;
	  }
			
			
			
			
  /*   Calcul de maxpixel et icomp   */
  err=imx_brukermax_3d(amax,amin,imres);
  rcoeff=(double) imres->rcoeff;

  /*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	if (im1->mri[i][j][k]>0)
	  {
	    amp1=(float)(1.0*log( im1rcoeff * (double) im1->mri[i][j][k]));
	    imres->mri[i][j][k]=(TYPEMRI3D)floor(amp1/rcoeff+0.5);
	  }
	else
	  imres->mri[i][j][k]=0;

  return(1);

}
/********************** Exponent **************/
int imx_exponent_3d_p(grphic3d *im1, grphic3d *imres)
{ 
  int i,j,k,wdth,hght,dpth;
  int err;
  float amp1,amax,amin;
  double rcoeff,im1rcoeff;

  hght=im1->height;
  wdth=im1->width;
  dpth=im1->depth;
  im1rcoeff= (double) im1->rcoeff;

  imx_copie_param_3d_p(im1,imres);

  /*    Recherche de amax   */
  amax=-1.0*LONGMAX; amin=LONGMAX;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
      //	if (im1->mri[i][j][k]>0)
	  {
	    amp1=(float)(1.0*exp( im1rcoeff * (double) im1->mri[i][j][k]));
	    if (amp1>amax)
	      amax=amp1;
	    if (amp1<amin)
	      amin=amp1;
	  }
			
			
			
			
  /*   Calcul de maxpixel et icomp   */
  err=imx_brukermax_3d(amax,amin,imres);
  rcoeff=(double) imres->rcoeff;

  /*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	//if (im1->mri[i][j][k]>0)
	   {
	    amp1=(float)(1.0*exp( im1rcoeff * (double) im1->mri[i][j][k]));
	    imres->mri[i][j][k]=(TYPEMRI3D)floor(amp1/rcoeff+0.5);
	   }
	//else
	 // imres->mri[i][j][k]=0;

  return(1);

}


/*******************************************
 ** --  imx_max_3d_p() ----------------------------
 */
/*!    \ingroup ToolBoxImage3d 
**
**     \brief Max function
**
**     \param im_1,im_2 : im_1 operation with im_2
**     \param im_res : Result image (E/S)
**	   \retval 1
********************************************/
int	imx_max_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{ 
  int i,j,k,wdth,hght,dpth;
  int err;
  float amp1,amp2,ampres;
  float amax,amin,rcoeff;
  float im1rcoeff,im2rcoeff;

  im1rcoeff=im1->rcoeff;
  im2rcoeff=im2->rcoeff;

  hght=MINI(im1->height,im2->height);
  wdth=MINI(im1->width,im2->width);
  dpth=MINI(im1->depth,im2->depth);

  imx_copie_param_3d_p(im1,imres);
  /*    Recherche de amax   */
  amax=0; amin=LONGMAX;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=im1rcoeff*im1->mri[i][j][k];
	  amp2=im2rcoeff*im2->mri[i][j][k];
	  ampres=MAXI(amp1,amp2);
	  if (ampres>amax)
	    amax=ampres;
	  amin=MINI(ampres,amin);
	}
  /*   Calcul de maxpixel et icomp   */
  err=imx_brukermax_3d(amax,amin,imres); /* Attention modif im2 si =imres */
  rcoeff=imres->rcoeff;

  /*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=im1rcoeff*im1->mri[i][j][k];
	  amp2=im2rcoeff*im2->mri[i][j][k];
	  ampres=MAXI(amp1,amp2);
	  imres->mri[i][j][k]=(TYPEMRI3D)floor(ampres/rcoeff+0.5);
	}

  return(1);

}


/*******************************************
 ** --  imx_min_3d_p() ----------------------------
 */
/*!    \ingroup ToolBoxImage3d 
**
**     \brief Max function
**
**     \param im_1,im_2 : im_1 operation with im_2
**     \param im_res : Result image (E/S)
**	   \retval 1
********************************************/
int	imx_min_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{ 
  int i,j,k,wdth,hght,dpth;
  int err;
  float amp1,amp2,ampres;
  float amax,amin,rcoeff;
  float im1rcoeff,im2rcoeff;

  im1rcoeff=im1->rcoeff;
  im2rcoeff=im2->rcoeff;

  hght=MINI(im1->height,im2->height);
  wdth=MINI(im1->width,im2->width);
  dpth=MINI(im1->depth,im2->depth);

  imx_copie_param_3d_p(im1,imres);
  /*    Recherche de amax   */
  amax=0; amin=LONGMAX;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=im1rcoeff*im1->mri[i][j][k];
	  amp2=im2rcoeff*im2->mri[i][j][k];
	  ampres=MINI(amp1,amp2);
	  if (ampres>amax)
	    amax=ampres;
	  amin=MINI(ampres,amin);
	}
  /*   Calcul de maxpixel et icomp   */
  err=imx_brukermax_3d(amax,amin,imres); /* Attention modif im2 si =imres */
  rcoeff=imres->rcoeff;

  /*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=im1rcoeff*im1->mri[i][j][k];
	  amp2=im2rcoeff*im2->mri[i][j][k];
	  ampres=MINI(amp1,amp2);
	  imres->mri[i][j][k]=(TYPEMRI3D)floor(ampres/rcoeff+0.5);
	}

  return(1);

}
/*******************************************
 ** --  imx_mul_coe_3d() ----------------------------
 **     Multiplicative function with coeff
 **     im_1,coef : im_1 operation with coef
 **     im_res : Result image
 ********************************************/

void	imx_mul_coe_3d(int im_1, float coefficient, int im_res)
{ 
  grphic3d *im1,*imres;

  im1=ptr_img_3d(im_1);
  imres=ptr_img_3d(im_res);
  
  imx_mul_coe_3d_p(im1,coefficient,imres);
}

/*******************************************
 ** --  imx_div_coe_3d() ----------------------------
 **     Divide function with coeff
 **     im_1,coef : im_1 operation with coef
 **     im_res : Result image
 ********************************************/

void	imx_div_coe_3d(int im_1, float coefficient, int im_res)
{ 
  grphic3d *im1,*imres;

  im1=ptr_img_3d(im_1);
  imres=ptr_img_3d(im_res);
  
  imx_div_coe_3d_p(im1,coefficient,imres);
}

/*******************************************
 ** --  imx_add_coe_3d() ----------------------------
 **     Additive function with coeff
 **     im_1,coef : im_1 operation with coef
 **     im_res : Result image
 ********************************************/

void	imx_add_coe_3d(int im_1, float coefficient, int im_res)
{ 
  grphic3d *im1,*imres;

  im1=ptr_img_3d(im_1);
  imres=ptr_img_3d(im_res);
  
  imx_add_coe_3d_p(im1,coefficient,imres);
}

/*******************************************
 ** --  imx_sub_coe_3d() ----------------------------
 **     Substractuve function with coeff
 **     im_1,coef : im_1 operation with coef
 **     im_res : Result image
 ********************************************/

void	imx_sub_coe_3d(int im_1, float coefficient, int im_res)
{ 
  grphic3d *im1,*imres;

  im1=ptr_img_3d(im_1);
  imres=ptr_img_3d(im_res);
  
  imx_sub_coe_3d_p(im1,coefficient,imres);
}

/*******************************************
 ** --  imx_mul_coe_3d_p() ----------------------------
 **     Multiplicative function with coeff
 **     im_1,coef : im_1 operation with coef
 **     im_res : Result image
 ********************************************/

int	imx_mul_coe_3d_p(grphic3d *im1, float coefficient, grphic3d *imres)
{
  int i,j,k,wdth,hght,dpth;
  int err;
  float amp1,ampres;
  float amax,amin;
  float rcoeff;
  float im1rcoeff;
  
  hght=im1->height;
  wdth=im1->width;
  dpth=im1->depth;
  im1rcoeff=im1->rcoeff;
  
  /*    Recherche de amax   */
  if (coefficient > 0 ) {
    amax=im1rcoeff*im1->max_pixel*coefficient;
    amin=im1rcoeff*im1->min_pixel*coefficient;
  }
  else  {	
    amin=im1rcoeff*im1->max_pixel*coefficient;
    amax=im1rcoeff*im1->min_pixel*coefficient;
  }
  /*   Calcul de maxpixel et icomp   */
  imx_copie_param_3d_p(im1,imres);
  err=imx_brukermax_3d(amax,amin,imres);
  rcoeff=imres->rcoeff;
  
  /*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=im1rcoeff*im1->mri[i][j][k];
	  ampres=amp1*coefficient;
	  imres->mri[i][j][k]=(TYPEMRI3D)floor(ampres/rcoeff+0.5);
	}

  return(1);
}

/*******************************************
 ** --  imx_div_coe_3d_p() ----------------------------
 **     Additive function with coeff
 **     im_1,coef : im_1 operation with coef
 **     im_res : Result image
 ********************************************/
int	imx_div_coe_3d_p(grphic3d *im1, float coefficient, grphic3d *imres)
{ 
  int i,j,k,wdth,hght,dpth;
  int err;
  float amp1,ampres;
  float amax,amin;
  float rcoeff,im1rcoeff;

  hght=im1->height;
  wdth=im1->width;
  dpth=im1->depth;
  im1rcoeff=im1->rcoeff;


  /*    Recherche de amax   */
  if (coefficient!=0) {
    amax=im1rcoeff*im1->max_pixel/coefficient;
    amin=im1rcoeff*im1->min_pixel/coefficient;
  }
  else
    return(0);
  /*   Calcul de maxpixel et icomp   */
  imx_copie_param_3d_p(im1,imres);
  err=imx_brukermax_3d(amax,amin,imres);
  rcoeff=imres->rcoeff;

  /*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=im1rcoeff*im1->mri[i][j][k];
	  ampres=amp1/coefficient;
	  imres->mri[i][j][k]=(TYPEMRI3D)floor(ampres/rcoeff+0.5);
	}

  return(1);

}

/*******************************************
 ** --  imx_add_coe_3d_p() ----------------------------
 **     Additive function with coeff
 **     im_1,coef : im_1 operation with coef
 **     im_res : Result image
 ********************************************/
int	imx_add_coe_3d_p(grphic3d *im1, float coefficient, grphic3d *imres)
{ 
  int i,j,k,wdth,hght,dpth;
  int err;
  float amp1,ampres;
  float amax,amin;
  float rcoeff,im1rcoeff;

  hght=im1->height;
  wdth=im1->width;
  dpth=im1->depth;
  im1rcoeff=im1->rcoeff;

  /*    Recherche de amax   */
  amax=im1rcoeff*im1->max_pixel+coefficient;
  amin=im1rcoeff*im1->min_pixel+coefficient;
  /*   Calcul de maxpixel et icomp   */
  imx_copie_param_3d_p(im1,imres);
  err=imx_brukermax_3d(amax,amin,imres);
  rcoeff=imres->rcoeff;

  /*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=im1rcoeff*im1->mri[i][j][k];
	  ampres=amp1+coefficient;
	  imres->mri[i][j][k]=(TYPEMRI3D)floor(ampres/rcoeff+.5);
	}

  imx_inimaxminpixel_3d_p(imres);

  return(1);

}

/*******************************************
 ** --  imx_sub_coe_3d_p() ----------------------------
 **     Soustractive function with coeff
 **     im_1,coef : im_1 operation with coef
 **     im_res : Result image
 ********************************************/
int	imx_sub_coe_3d_p(grphic3d *im1, float coefficient, grphic3d *imres)
{ 
  int i,j,k,wdth,hght,dpth;
  int err;
  float amp1,ampres;
  float amax,amin;
  float rcoeff,im1rcoeff;

  hght=im1->height;
  wdth=im1->width;
  dpth=im1->depth;
  im1rcoeff=im1->rcoeff;

  /*    Recherche de amax   */
  amax=im1rcoeff*im1->max_pixel-coefficient;
  amin=im1rcoeff*im1->min_pixel-coefficient;
  /*   Calcul de maxpixel et icomp   */
  imx_copie_param_3d_p(im1,imres);
  err=imx_brukermax_3d(amax,amin,imres);
  rcoeff=imres->rcoeff;

  /*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=im1rcoeff*im1->mri[i][j][k];
	  ampres=amp1-coefficient;
	  imres->mri[i][j][k]=(TYPEMRI3D)floor(ampres/rcoeff+0.5);
	}

  return(1);

}

/*******************************************
 ** --  imx_and_3d() ----------------------------
 */
/*!    \brief Logical AND function
**     \param im_1,im_2 : im_1 operation with im_2
**     \param im_res : Result image 
**	   \retval 1
********************************************/
int	imx_and_3d(int im_1, int im_2, int im_res)
{ 
  grphic3d *im1,*im2,*imres;
  

  im1=ptr_img_3d(im_1);
  im2=ptr_img_3d(im_2);
  imres=ptr_img_3d(im_res);

  imx_and_3d_p(im1,im2,imres);
  
  return(1);
}

/*******************************************
 ** --  imx_or_3d() ----------------------------
 */
/*!    \brief Logical OR function
**     \param im_1,im_2 : im_1 operation with im_2
**     \param im_res : Result image 
**	   \retval 1
********************************************/
int	imx_or_3d(int im_1, int im_2, int im_res)
{ 
  grphic3d *im1,*im2,*imres;
  

  im1=ptr_img_3d(im_1);
  im2=ptr_img_3d(im_2);
  imres=ptr_img_3d(im_res);

  imx_or_3d_p(im1,im2,imres);
  
  return(1);
}

/*******************************************
 ** --  imx_inv_3d() ----------------------------
 */
/*!    \brief Logical INVERSION function
**
**     Soustraction de l'image de son maximum
**     \param im_1 :numero de l'image
**     \param im_res :numero de l'image resultat
**	   \retval 1
********************************************/
int	imx_inv_3d(int im_1, int im_res)
{ 
  grphic3d *im1,*imres;
  

  im1=ptr_img_3d(im_1);
  imres=ptr_img_3d(im_res);

  imx_inv_3d_p(im1,imres);
  return(1);

}

/*******************************************
 ** --  imx_eg_3d() ----------------------------
 */
/*!    \brief Logical EGALITE function
**     \param im_1,im_2 :numero des images
**     \param im_res :numero de l'image resultat
**	   \retval 1
********************************************/
int	imx_eg_3d(int im_1, int im_2, int im_res)
{ 
  grphic3d *im1,*im2,*imres;
  

  im1=ptr_img_3d(im_1);
  im2=ptr_img_3d(im_2);
  imres=ptr_img_3d(im_res);

  imx_eg_3d_p(im1,im2,imres);
  
  return(1);
}

/*******************************************
 ** --  imx_and2img_3d() ----------------------------
 */
/*!    \brief Logical And applied on 2 images
**     \param im_1,im_2 :numero des images
**	   \retval 1
********************************************/
int	imx_and2img_3d(int im_1, int im_2)
{ 
  grphic3d *im1,*im2;
  

  im1=ptr_img_3d(im_1);
  im2=ptr_img_3d(im_2);

  imx_and2img_3d_p(im1,im2);
  
  return(1);
}

/*******************************************
 ** --  imx_and_3d_p() ----------------------------
 */
/*! \ingroup ToolBoxImage3d
**	   \brief Logical AND function
**     \param im1,im2 : images
**     \param imres :image resultat (E/S)
**	   \retval 1
********************************************/
int	imx_and_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{
  
  int i,j,k,wdth,hght,dpth;
  int err;
  float amp1,ampres;
  float amax,amin,im1rcoeff,im2rcoeff;

  hght=MINI(im1->height,im2->height);
  wdth=MINI(im1->width,im2->width);
  dpth=MINI(im1->depth,im2->depth);
  im1rcoeff=im1->rcoeff; 
  im2rcoeff=im2->rcoeff; 

  imx_copie_param_3d_p(im1,imres);

  /*    Recherche de l'amplitude max reelle dans l'image resultat */
  amax=0.; amin=LONGMAX;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
        {
	  amp1=im1rcoeff*im1->mri[i][j][k];
	  if (im1->mri[i][j][k]&&im2->mri[i][j][k]) ampres=amp1;
	  else ampres=0.;
	  amax=MAXI(amax,ampres);
	  amin=MINI(amin,ampres);
        }
        
  /*   Calcul de maxpixel et icomp en fonction de amax  */
  err=imx_brukermax_3d(amax,amin,imres);
  
  /*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
        {
	  amp1=im1rcoeff*im1->mri[i][j][k];
	  if (im1->mri[i][j][k]&&im2->mri[i][j][k]) ampres=amp1;
	  else ampres=0;
	  imres->mri[i][j][k]=(TYPEMRI3D)floor(ampres/(imres->rcoeff));
        }
  
  return(1);
}

/*******************************************
 ** --  imx_or_3d_p() ----------------------------
 */
/*! \ingroup ToolBoxImage3d
**	   \brief Logical OR function
**     \param im1,im2 : images
**     \param imres :image resultat (E/S)
**	   \retval 1
********************************************/
int	imx_or_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{ 
  int i,j,k,wdth,hght,dpth;
  int err;
  float amp1,amp2,ampres;
  float amax,amin;
  float rcoeff,im1rcoeff,im2rcoeff;

  im1rcoeff=im1->rcoeff;
  im2rcoeff=im2->rcoeff;
  hght=MINI(im1->height,im2->height);
  wdth=MINI(im1->width,im2->width);
  dpth=MINI(im1->depth,im2->depth);

  imx_copie_param_3d_p(im1,imres);

  /*    Recherche de amax   */
  amax=0; amin=LONGMAX;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=im1rcoeff*im1->mri[i][j][k];
	  amp2=im2rcoeff*im2->mri[i][j][k];
	  ampres=amp1;
	  if (im1->mri[i][j][k]==0)ampres=amp2;
	  if (ampres>amax)
	    amax=ampres;
	  amin=MINI(amin,ampres);
	}
  /*   Calcul de maxpixel et icomp   */
  err=imx_brukermax_3d(amax,amin,imres);
  rcoeff=imres->rcoeff;

  /*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=im1rcoeff*im1->mri[i][j][k];
	  amp2=im2rcoeff*im2->mri[i][j][k];
	  ampres=amp1;
	  if (im1->mri[i][j][k]==0)ampres=amp2;
	  imres->mri[i][j][k]=(TYPEMRI3D)floor(ampres/rcoeff);
	}

  return(1);

}

/*******************************************
 ** --  imx_eg_3d_p() ----------------------------
 */
/*! \ingroup ToolBoxImage3d
**    \brief Logical EGALITE function
**     \param im1,im2 : images
**     \param imres :image resultat (E/S)
**	   \retval 1
********************************************/
int	imx_eg_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{ 
  int i,j,k,wdth,hght,dpth;
  int err;
  float amp1,amp2,ampres;
  float amax,amin;
  float rcoeff,im1rcoeff,im2rcoeff;

  hght=MINI(im1->height,im2->height);
  wdth=MINI(im1->width,im2->width);
  dpth=MINI(im1->depth,im2->depth);
  im1rcoeff=im1->rcoeff;
  im2rcoeff=im2->rcoeff;

  imx_copie_param_3d_p(im1,imres);

  /*    Recherche de amax   */
  amax=0.;amin=LONGMAX;
  ampres=0.;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=im1rcoeff*im1->mri[i][j][k];
	  amp2=im2rcoeff*im2->mri[i][j][k];
	  if (amp1==amp2)ampres=amp1;
	  if (ampres>amax)
	    amax=ampres;
	  amin=MINI(ampres,amin);
	}
  /*   Calcul de maxpixel et icomp   */
  err=imx_brukermax_3d(amax,amin,imres);
  rcoeff=imres->rcoeff;

  /*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=im1rcoeff*im1->mri[i][j][k];
	  amp2=im2rcoeff*im2->mri[i][j][k];
	  ampres=0.;
	  if (amp1==amp2)ampres=amp1;
	  imres->mri[i][j][k]=(TYPEMRI3D)floor(ampres/rcoeff);
	}

  return(1);

}


/*******************************************
 ** --  imx_and2img_3d_p() ----------------------------
 */
/*! \ingroup ToolBoxImage3d
**
**     Logical And between 2 imgs
**    
**     \param im_1 : Image to And with im_2
**     \param im_2 : Image to And with im_1
**	   \retval 1
********************************************/
int    imx_and2img_3d_p(grphic3d *im1, grphic3d *im2)
{
  int i,j,k;
  int w,h,d;

  if ( (im1->width != im2->width) ||
       (im1->height!= im2->height) ||
       (im1->depth != im2->depth) )
    {
      PUT_WARN("\nAttention les 2 images n'ont pas la meme taille\n");
    }

  w = MINI(im1->width,im2->width);
  h = MINI(im1->height,im2->height);
  d = MINI(im1->depth,im2->depth);

  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      for (k=0;k<d;k++)
	{
	  if ( !(im1->mri[i][j][k] && im2->mri[i][j][k]))
	    {
	      im1->mri[i][j][k]=0;
	      im2->mri[i][j][k]=0;
	    }
	} 

  return(1);
}



/*******************************************
 ** --  imx_inv_3d_p() ----------------------------
 */
/*! \ingroup ToolBoxImage3d
**
**     Logical INVERSION function
**    
**     \param im_1 : Image to invert
**     \param im_res : Result image (E/S)
**	   \retval 1
********************************************/
int    imx_inv_3d_p(grphic3d *im1, grphic3d *imres)
{
  int i,j,k;
  int w,h,d;
  imx_copie_param_3d_p(im1,imres); 
   
  w = im1->width;
  h = im1->height;
  d = im1->depth;

  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      for (k=0;k<d;k++)
	{
	  if ((im1->mri[i][j][k])==0)
	    {imres->mri[i][j][k]=1;}
	  else
	    {imres->mri[i][j][k]=0;}
	} 
  imres->icomp=0;
  imres->rcoeff=1;
  imres->max_pixel=1;
  imres->cutoff_max=1;
  return(1);
}

