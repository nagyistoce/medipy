/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*!**********************************************************************
***	
***	\file		norm_3d.c
***
***	project:	Imagix 2.01
***			
***
***	\brief description:	normalisation d'images
***	
***	
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***     Last user action by: Mr. ARMSPACH on Dec 15th 1996
***
*************************************************************************/
#include <config.h>
#include <string.h>

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h" 
#include "noyau/gui/imx_picture3d.h"
#include "noyau/imx_3d.h"
#include "noyau/imx_lang.h"
#include "math/oper_3d.h"            // -> imx_mul_coe_3d_p()
#include "noyau/io/imx_export_file.h"    // -> PUTI()
#include "traitement/norm_3d.h"
#include "detection/detection_anomalies_3d.h"

/***********************************************************
*** Normalisation moyenne/ecart-type avec des seuils d�finis pas l'utilisateur
**********************************************************/	
void norm_seuil_coeff_3d(void)
{
	int im_1,im_res, err;
	float mean,std;
		
	im_1=GET_PLACE3D(TEXT0231);
	im_res=GET_PLACE3D(TEXT0006);
	mean=GET_FLOAT("Moyenne ?", 100, &err);
	std=GET_FLOAT("Ecart type ?", 50, &err);
		
	
	imx_norm_seuil_coeff_3d(im_1,im_res,mean,std);
	
	show_picture_3d(im_res);
}


/***********************************************************
*** Normalisation moyenne/ecart-type avec des seuils d�finis pas l'utilisateur
**********************************************************/	
void imx_norm_seuil_coeff_3d(int im_1,int im_res,float mean,float std)
{
	grphic3d *im,*imres;
  
	im=ptr_img_3d(im_1);	
	imres=ptr_img_3d(im_res);
 	imx_norm_seuil_coeff_3d_p(im,imres, mean, std);
}

/***********************************************************
*** Normalisation moyenne/ecart-type avec des seuils d�finis pas l'utilisateur
**********************************************************/	
void imx_norm_seuil_coeff_3d_p(grphic3d * im,grphic3d * imres,float mean,float std)
{
int i,j,k;
int height,width,depth,nb_tot;
double moy1,var1,sigma1;
double rcoeff1,rcoeff;
double aux1,aux2;
float max,min;
int err;

height=im->height;
width=im->width;
depth=im->depth;
nb_tot=0;
rcoeff1=im->rcoeff;
moy1=0;
var1=0;
aux1=0;
aux2=0;
max=0;

for (i=0;i<height;i++)
   for (j=0;j<width;j++)
      for (k=0;k<depth;k++)
		if  (im->mri[i][j][k]>0)
		{
		moy1 += im->mri[i][j][k];
		nb_tot++;
		}
moy1= rcoeff1*moy1/(double)nb_tot;
	
	
for (i=0;i<height;i++)
   for (j=0;j<width;j++)
      for (k=0;k<depth;k++)
		if  (im->mri[i][j][k]>0)
      	{
		var1 +=(rcoeff1*im->mri[i][j][k]-moy1)*(rcoeff1*im->mri[i][j][k]-moy1);
		}

var1= var1/(double)nb_tot;
sigma1=sqrt(var1);
	
aux1=rcoeff1*std/sigma1;
aux2=mean-moy1*std/sigma1;

imx_copie_param_3d_p(im,imres);
max=(float)(aux1*(im->max_pixel)+aux2);
min=(float)(aux1*(im->min_pixel)+aux2);


/*   Calcul de imres->max_pixel et imres->icomp imres->rcoeff a partir de max*/
err= imx_brukermax_3d(max,min,imres);
rcoeff= imres->rcoeff;

for (i=0;i<width;i++)
   for (j=0;j<height;j++)
      for (k=0;k<depth;k++)
      	if  (im->mri[i][j][k]>0)
		{
      	imres->mri[i][j][k] = (TYPEMRI3D)floor((aux1*(im->mri[i][j][k])+aux2)/rcoeff);
      	}
		else
		imres->mri[i][j][k]=0;
		
imx_inimaxminpixel_3d_p(imres);
}




/*    norm_roi_3d()   ------------------------------------------------------
***	NORMALISATION D'UNE IMAGE 3D PAR RAPPORT A UNE AUTRE A L'AIDE D'UNE ROI
***
***
**********************************************************/	
void norm_roi_3d(void)
{
	int im_1,im_2,im_res;
	int type_norm;

	im_1=GET_PLACE3D(TEXT0231);
	im_2=GET_PLACE3D(TEXT0031);
	im_res=GET_PLACE3D(TEXT0006);

  type_norm=imx_query_normalisation_type_noReca_3d();
	
	imx_norm_roi_3d(im_1,im_2,im_res,type_norm);
	
	show_picture_3d(im_res);
}

/***    imx_norm_roi_3d()   ------------------------------------------------------
***	NORMALISATION D'UNE IMAGE 3D PAR RAPPORT A UNE AUTRE A L'AIDE D'UNE ROI
***	im_1	:pointeur de l'image 1
***	im_2	:pointeur de l'image 2
***     roi 	: region d'interet
***	im_res	:pointeur de l'image resultat
***	type_norm :type de normalisation
***		0 : moyenne
***		1 : moyenne +ecty
***		2 : regression lineaire
***		3 : rapport des moyennes
***
***
**********************************************************/	
void imx_norm_roi_3d(int im_1, int im_2, int im_res, int type_norm)
{
	grphic3d *im1,*im2,*imres;
  norm_func_t normalisation=NULL;

	im1=ptr_img_3d(im_1);	
	im2=ptr_img_3d(im_2);
	imres=ptr_img_3d(im_res);

  if (!ptr_has_mask_3d(im_1) || !ptr_has_mask_3d(im_2))
  {
    PUT_WARN("!! NO MASK !!\n");
    return ;
  }

  normalisation=imx_choose_normalisation_roi_noReca_fct(type_norm);
  normalisation(im1,im2,imres);
	
}

/*    norm_roi_3d()   ------------------------------------------------------
***	NORMALISATION D'UNE IMAGE 3D PAR RAPPORT A UNE AUTRE A L'AIDE D'UNE ROI
***
***
**********************************************************/
void norm_roi_reca_3d(void)
{
	int im_1,im_2,im_res;
	int type_norm;

	im_1=GET_PLACE3D(TEXT0231);
	im_2=GET_PLACE3D(TEXT0031);
	im_res=GET_PLACE3D(TEXT0006);

  type_norm=imx_query_normalisation_type_reca_3d();

	imx_norm_roi_reca_3d(im_1,im_2,im_res,type_norm);

	show_picture_3d(im_res);
}

/***    imx_norm_roi_3d()   ------------------------------------------------------
***	NORMALISATION D'UNE IMAGE 3D PAR RAPPORT A UNE AUTRE A L'AIDE D'UNE ROI
***	im_1	:pointeur de l'image 1
***	im_2	:pointeur de l'image 2
***     roi 	: region d'interet
***	im_res	:pointeur de l'image resultat
***	type_norm :type de normalisation
***		0 : moyenne
***		1 : moyenne +ecty
***		2 : regression lineaire
***		3 : rapport des moyennes
***
***
**********************************************************/
void imx_norm_roi_reca_3d(int im_1, int im_2, int im_res, int type_norm)
{
	grphic3d *im1,*im2,*imres;
  norm_func_t normalisation=NULL;

	im1=ptr_img_3d(im_1);
	im2=ptr_img_3d(im_2);
	imres=ptr_img_3d(im_res);

  if (!ptr_has_mask_3d(im_1) || !ptr_has_mask_3d(im_2))
  {
    PUT_WARN("!! NO MASK !!\n");
    return ;
  }

  normalisation=imx_choose_normalisation_roi_reca_fct(type_norm);
  normalisation(im1,im2,imres);

}

/*    norm_im_3d()   --------------------------------------
***	NORMALISATION D'UNE IMAGE PAR RAPPORT A UNE AUTRE NON SEUILLEE
***
***
***
**********************************************************/	
void norm_im_3d(void)
{
	int im_1,im_2,im_res;
	int i,type_norm;
	char *quest[6];	

	for(i=0;i<6;i++) quest[i]=CALLOC(80,char);
	strcpy(quest[0],TEXT0229);
	strcpy(quest[1],TEXT0230);
	strcpy(quest[2],TEXT0233);
	strcpy(quest[3],TEXT0234);
	strcpy(quest[4],TEXT1234);
	strcpy(quest[5],"\0");

	im_1=GET_PLACE3D(TEXT0231);
	im_2=GET_PLACE3D(TEXT0031);
	im_res=GET_PLACE3D(TEXT0006);
	type_norm=GETV_QCM(TEXT0232,(char **)quest);

	imx_norm_im_3d(im_1,im_2,im_res,type_norm);
	
	for(i=0;i<6;i++)
           FREE(quest[i]);

	show_picture_3d(im_res);

}

/****   imx_norm_im_3d    ****************************************
****
***  NORMALISATION D'UNE IMAGE PAR RAPPORT A UNE AUTRE NON SEUILLEE
***
***	im_1	:pointeur de l'image a normaliser
***	im_2	:pointeur de l'image de reference
***	im_res	:pointeur de l'image resultat
***	type_norm :type de normalisation
***		0 : moyenne
***		1 : moyenne +ecty
***		2 : regression lineaire
***		3 : rapport des moyennes
***   4 : rapport min/max
***
*************************************************************/
void imx_norm_im_3d(int im_1, int im_2, int im_res, int type_norm)
{
	grphic3d *im1,*im2,*imres;
	
	im1=ptr_img_3d(im_1);	
	im2=ptr_img_3d(im_2);
	imres=ptr_img_3d(im_res);
	
        switch (type_norm) {
	case 0:  /* Moyenne */
	       imx_norm_im_mean_3d_p(im1,im2,imres);
	       break;
	case 1:  /* Moyenne + ecty */
	       imx_norm_im_meanecty_3d_p(im1,im2,imres);
	       break;
	case 2:  /* Regression lineaire */
	       imx_norm_im_rl_3d_p(im1,im2,imres);
	       break;
	case 3:  /* rapport des moyennes */
	       imx_norm_im_3d_p(im1,im2,imres);
	       break;
	case 4:  /* rapport des min/max */
	       imx_norm_im_minmax_3d_p(im1,im2,imres);
	       break;
	default:  /* defaut moyenne +ecty */
	       imx_norm_im_meanecty_3d_p(im1,im2,imres);
	       break;
          }
	
}


/*     norm_seuil_3d()   -----------------------------------------
***	NORMALISATION D'UNE IMAGE 3D PAR RAPPORT A UNE IMAGE 
***	SEUILLE
***	im_1	:pointeur de l'image a normaliser
***	im_2	:pointeur de l'image de reference
***	im_res	:pointeur de l'image resultat
***
**********************************************************/	
void norm_seuil_3d(void)
{
	int im_1,im_2,im_res;
	int i,type_norm;
	char *quest[6];	

	for(i=0;i<6;i++) quest[i]=CALLOC(80,char);
	strcpy(quest[0],TEXT0229);
	strcpy(quest[1],TEXT0230);
	strcpy(quest[2],TEXT0233);
	strcpy(quest[3],TEXT0234);
	strcpy(quest[4],TEXT1234);
	strcpy(quest[5],"\0");
	
	im_1=GET_PLACE3D(TEXT0231);
	im_2=GET_PLACE3D(TEXT0031);
	im_res=GET_PLACE3D(TEXT0006);
	type_norm=GETV_QCM(TEXT0232,(char **)quest);
	
	imx_norm_seuil_3d(im_1,im_2,im_res,type_norm);
        
	for(i=0;i<6;i++)
           FREE(quest[i]);
	
	show_picture_3d(im_res);
}

/***********************************************************
***	imx_norm_seuil_3d()   
***	NORMALISATION D'UNE IMAGE 3D PAR RAPPORT A UNE IMAGE 
***	SEUILLE
***	im_1	:pointeur de l'image a normaliser
***	im_2	:pointeur de l'image de reference
***	im_res	:pointeur de l'image resultat
***	type_norm :type de normalisation
***		0 : moyenne
***		1 : moyenne +ecty
***		2 : regression lineaire
***		3 : rapport des moyennes
***
***
**********************************************************/	
void imx_norm_seuil_3d(int im_1, int im_2, int im_res, int type_norm)
{	
	grphic3d *im1,*im2,*imres;
	
	im1=ptr_img_3d(im_1);	
	im2=ptr_img_3d(im_2);
	imres=ptr_img_3d(im_res);
	
        switch (type_norm) {
	case 0:  /* Moyenne */
	       imx_norm_seuil_mean_3d_p(im1,im2,imres);
	       break;
	case 1:  /* Moyenne + ecty */
	       imx_norm_seuil_meanecty_3d_p(im1,im2,imres);
	       break;
	case 2:  /* Regression lineaire */
	       imx_norm_seuil_rl_3d_p(im1,im2,imres);
	       break;
	case 3:  /* rapport des moyennes */
	       imx_norm_seuil_3d_p(im1,im2,imres);
	       break;
	case 4:  /* rapport des moyennes */
	       imx_norm_seuil_minmax_3d_p(im1,im2,imres);
	       break;
	default:  /* defaut moyenne +ecty */
	       imx_norm_seuil_meanecty_3d_p(im1,im2,imres);
	       break;
          }
	
}

/***    imx_norm_roi_3d_p()   ------------------------------------------------------
***	NORMALISATION D'UNE IMAGE 3D PAR RAPPORT A UNE AUTRE A L'AIDE DES MASQUES
***
***   Normalisation en utilisant le rapport des moyennes
***  
***************************************************************/
void imx_norm_roi_3d_p(grphic3d *im1, grphic3d *im2,  grphic3d *imres)
{
  int taille1=0, taille2=0;
  double moy1, moy2;
  grphic3d *mask1,*mask2;

  mask1 = im1->mask;
  mask2 = im2->mask;
  if (!mask1 || !mask2)
  { fprintf (stderr, "les images n'ont pas de masque dans imx_norm_roi_3d_p\n"); return ;}

  //calcul des tailles des masques et des moyennes des images sur les masques
  imx_inimaxminpixel_3d_p(im1);
  imx_inimaxminpixel_3d_p(im2);
  taille1=imx_calc_moy_im_mask_3d_p(im1,&moy1);
  taille2=imx_calc_moy_im_mask_3d_p(im2,&moy2);

  if ((taille1<=0)||(taille2<=0))
  { fprintf (stderr, "impossible de normaliser pour cause de taile de masque nulle dans imx_norm_roi_3d_p\n"); return; }

  //normalisation lineaire de l'image1
  imx_norm_lin_im_3d_p(im1, imres, moy2/moy1, 0.0);

  //Mettre le meme rcoeff pour les deux images
  imx_norm_rcoeff_3d_p(im2,imres);
}

/***    imx_norm_roi_3d_p()   ------------------------------------------------------
***	NORMALISATION D'UNE IMAGE 3D PAR RAPPORT A UNE AUTRE A L'AIDE DES MASQUES
***
***   Normalisation en utilisant le rapport min/max
***
***************************************************************/
void imx_norm_roi_minmax_3d_p(grphic3d *im1, grphic3d *im2,  grphic3d *imres)
{
  int max1,min1;
  int max2,min2;
  double coeffA,coeffB;
  int total1,total2;
  grphic3d *mask1,*mask2;
  
  mask1 = im1->mask;
  mask2 = im2->mask;
  if (!mask1 || !mask2)
  { fprintf (stderr, "les images n'ont pas de masque dans imx_norm_roi_minmax_3d_p\n"); return ;}

  //calcul des tailles des masques et des moyennes des images sur les masques
  imx_inimaxminpixel_3d_p(im1);
  imx_inimaxminpixel_3d_p(im2);
  total1=imx_calc_taille_mask_3d_p(im1);
  total2=imx_calc_taille_mask_3d_p(im2);
  max1=im1->max_pixel; min1=im1->min_pixel;
  max2=im2->max_pixel; min2=im2->min_pixel;

  if ((total1<=0)||(total2<=0))
  { fprintf (stderr, "impossible de normaliser pour cause de taile de masque nulle dans imx_norm_roi_minmax_3d_p\n"); return; }

  imx_copie_param_3d_p(im1,imres);

  coeffA=((max2-min2)*im2->rcoeff) / ((max1-min1)*im1->rcoeff);
  coeffB=min2*im2->rcoeff - min1*im1->rcoeff*coeffA;

  //normalisation lineaire de l'image1
  imx_norm_lin_im_3d_p(im1, imres, coeffA, coeffB);

  //Mettre le meme rcoeff pour les deux images
  imx_norm_rcoeff_3d_p(im2,imres);
}

/***    imx_norm_roi_meanecty_3d_p()   ------------------------------------------------------
***	NORMALISATION D'UNE IMAGE 3D PAR RAPPORT A UNE AUTRE A L'AIDE DES MASQUES
***
***    normalisation moyenne et ecart type
***
****************************************************************/
void imx_norm_roi_meanecty_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{
  double moy1=0.0,moy2=0.0,sygma1=0.0,sygma2=0.0;
  double alpha,beta;
  int taille1=0,taille2=0;
  grphic3d *mask1,*mask2;
  
  mask1 = im1->mask;
  mask2 = im2->mask;
  if (!mask1 || !mask2)
  { fprintf (stderr, "les images n'ont pas de masque dans imx_norm_roi_meanecty_3d_p\n"); return ;}

  //calcul des tailles des masques et des moyennes des images sur les masques
  imx_inimaxminpixel_3d_p(im1);
  imx_inimaxminpixel_3d_p(im2);
  taille1=imx_calc_moy_im_mask_3d_p(im1,&moy1);
  taille2=imx_calc_moy_im_mask_3d_p(im2,&moy2);

  if ((taille1<=0)||(taille2<=0))
  { fprintf (stderr, "impossible de normaliser pour cause de taile de masque nulle dans imx_norm_roi_meanecty_3d_p\n"); return; }

  //calcul des ecarts types des images
  imx_calc_sigma_im_mask_3d_p(im1, taille1, moy1, &sygma1);
  imx_calc_sigma_im_mask_3d_p(im2, taille2, moy2, &sygma2);

  alpha=sygma2/sygma1;
  beta=moy2-moy1*sygma2/sygma1;

  //normalisation lineaire de l'image1
  imx_norm_lin_im_3d_p(im1, imres, alpha, beta);

  //Mettre le meme rcoeff pour les deux images 
  imx_norm_rcoeff_3d_p(im2,imres);
  return;	
}

/***    imx_norm_roi_mean_3d_p()   ------------------------------------------------------
***	NORMALISATION D'UNE IMAGE 3D PAR RAPPORT A UNE AUTRE A L'AIDE DES MASQUES
***
***    normalisation moyenne
***
****************************************************************/
void imx_norm_roi_mean_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{
  double moy1=0.0,moy2=0.0;
  double ecartMoy=0.0;
  int taille1=0, taille2=0;
  grphic3d *mask1,*mask2;

  mask1 = im1->mask;
  mask2 = im2->mask;
  if (!mask1 || !mask2)
  { fprintf (stderr, "les images n'ont pas de masque dans imx_norm_roi_mean_3d_p\n"); return ;}

  //calcul des tailles des masques et des moyennes des images sur les masques
  imx_inimaxminpixel_3d_p(im1);
  imx_inimaxminpixel_3d_p(im2);
  taille1=imx_calc_moy_im_mask_3d_p(im1,&moy1);
  taille2=imx_calc_moy_im_mask_3d_p(im2,&moy2);

  if ((taille1<=0)||(taille2<=0))
  { fprintf (stderr, "impossible de normaliser pour cause de taile de masque nulle dans imx_norm_roi_mean_3d_p\n"); return; }

  ecartMoy=moy2-moy1;

  //normalisation lineaire de l'image1
  imx_norm_lin_im_3d_p(im1, imres, 1.0, ecartMoy);

  // Mettre le meme rcoeff pour les deux images 
  imx_norm_rcoeff_3d_p(im2,imres);
  return;	
}

/***    imx_norm_roi_rl_3d_p()   ------------------------------------------------------
***	NORMALISATION D'UNE IMAGE 3D PAR RAPPORT A UNE AUTRE A L'AIDE DES MASQUES
***
***    normalisation en utilisant la regression lineaire
***
****************************************************************/
void imx_norm_roi_rl_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{
  int taille1=0,taille2=0,taille_covar=0;
  double moy1=0.0,moy2=0.0,covar=0.0,sygma1=0.0;
  double alpha,beta;
  grphic3d *mask1,*mask2;

  mask1 = im1->mask;
  mask2 = im2->mask;
  if (!mask1 || !mask2)
  { fprintf (stderr, "les images n'ont pas de masque dans imx_norm_roi_rl_3d_p\n"); return ;}

  //calcul des tailles des masques et des moyennes des images sur les masques
  imx_inimaxminpixel_3d_p(im1);
  imx_inimaxminpixel_3d_p(im2);
  taille1=imx_calc_moy_im_mask_3d_p(im1,&moy1);
  taille2=imx_calc_moy_im_mask_3d_p(im2,&moy2);

  if ((taille1<=0)||(taille2<=0))
  { fprintf (stderr, "impossible de normaliser pour cause de taille de masque nulle dans imx_norm_roi_rl_3d_p\n"); return; }

  //calcul des ecarts types et covariance des images
  imx_calc_sigma_im_mask_3d_p(im1, taille1, moy1, &sygma1);
  taille_covar=imx_calc_covar_im_mask_3d_p(im1, im2, moy1, moy2, &covar);

  if (taille_covar<=0)
  { fprintf (stderr, "impossible de normaliser pour cause de taille d'overlap des masques nulle dans imx_norm_roi_rl_3d_p\n"); return; }

  alpha=covar/(sygma1*sygma1);
  beta=moy2-(moy1*covar)/(sygma1*sygma1);

  //normalisation lineaire de l'image1
  imx_norm_lin_im_3d_p(im1, imres, alpha, beta);

  // Mettre le meme rcoeff pour les deux images
  imx_norm_rcoeff_3d_p(im2,imres);
  return;
}

/***    imx_norm_roi_rl_3d_p()   ------------------------------------------------------
***	NORMALISATION D'UNE IMAGE 3D PAR RAPPORT A UNE AUTRE A L'AIDE DES MASQUES
***
***    normalisation en utilisant la regression lineaire mais a 1 parametre
***
****************************************************************/
void imx_norm_roi_rl_1_param_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{
  int taille1=0,taille2=0,taille_covar=0;
  double moy1=0.0,moy2=0.0,covar=0.0,sygma1=0.0;
  double alpha;
  grphic3d *mask1,*mask2;

  mask1 = im1->mask;
  mask2 = im2->mask;
  if (!mask1 || !mask2)
  { fprintf (stderr, "les images n'ont pas de masque dans imx_norm_roi_rl_3d_p\n"); return ;}

  //calcul des tailles des masques et des moyennes des images sur les masques
  imx_inimaxminpixel_3d_p(im1);
  imx_inimaxminpixel_3d_p(im2);
  taille1=imx_calc_moy_im_mask_3d_p(im1,&moy1);
  taille2=imx_calc_moy_im_mask_3d_p(im2,&moy2);

  if ((taille1<=0)||(taille2<=0))
  { fprintf (stderr, "impossible de normaliser pour cause de taille de masque nulle dans imx_norm_roi_rl_3d_p\n"); return; }

  //calcul des ecarts types et covariance des images
  imx_calc_sigma_im_mask_3d_p(im1, taille1, moy1, &sygma1);
  taille_covar=imx_calc_covar_im_mask_3d_p(im1, im2, moy1, moy2, &covar);

  if (taille_covar<=0)
  { fprintf (stderr, "impossible de normaliser pour cause de taille d'overlap des masques nulle dans imx_norm_roi_rl_3d_p\n"); return; }

  alpha=covar/(sygma1*sygma1);

  //normalisation lineaire de l'image1
  imx_norm_lin_im_3d_p(im1, imres, alpha, 0.0);

  // Mettre le meme rcoeff pour les deux images
  imx_norm_rcoeff_3d_p(im2,imres);
  return;
}

int imx_calc_taille_mask_3d_p(grphic3d *im)
{
 TDimension width, height, depth;
 unsigned int i,j,k;
 
 int taille=0;
 grphic3d *mask=NULL;

 height=im->height; width=im->width; depth=im->depth;
 mask=im->mask;
 if (!mask) { fprintf (stderr, "l'image n'a pas de masque dans imx_calc_taille_mask_3d_p\n"); return 0; }

 for(i=0;i<width;i++)
  for(j=0;j<height;j++)
   for(k=0;k<depth;k++)
   {
    if (mask->mri[i][j][k]!=0)
    {
       taille++;
    }
   }

 return taille;
}

int imx_calc_moy_im_mask_3d_p(grphic3d *im, double *moyenne)
{
 TDimension width, height, depth;
 unsigned int i,j,k; 
 int taille=0;
 grphic3d *mask=NULL;
 double moy=0.0;

 height=im->height; width=im->width; depth=im->depth;
 mask=im->mask;
 if (!mask) { fprintf (stderr, "l'image n'a pas de masque dans imx_calc_moy_im_mask_3d_p\n"); *moyenne=0.0; return 0; }

 for(i=0;i<width;i++)
  for(j=0;j<height;j++)
   for(k=0;k<depth;k++)
   {
    if (mask->mri[i][j][k]!=0)
    {
       moy += im->mri[i][j][k];
       taille++;
    }
   }

 if (taille>0) moy=im->rcoeff*moy/taille;
 *moyenne=moy;

 return taille;
}

int imx_calc_sigma_im_mask_3d_p(grphic3d *im, int taille, double moyenne, double *sigma)
{
 TDimension width, height, depth;
 unsigned int i,j,k;
 grphic3d *mask=NULL;
 double var=0.0, rcoeff;
 double val=0.0;

 height=im->height; width=im->width; depth=im->depth;
 mask=im->mask;
 if (!mask) { fprintf (stderr, "l'image n'a pas de masque dans imx_calc_sigma_im_mask_3d_p\n"); *sigma=0.0; return 1; }
 rcoeff=im->rcoeff;

 for(i=0;i<width;i++)
  for(j=0;j<height;j++)
   for(k=0;k<depth;k++)
   {
    if (mask->mri[i][j][k]!=0)
    {
      val=rcoeff*im->mri[i][j][k]-moyenne;
      var +=val*val/taille;
    }
   }

 var=sqrt(var);
 *sigma=var;

 return 0;
}

int imx_calc_covar_im_mask_3d_p(grphic3d *im1, grphic3d *im2, double moyenne1, double moyenne2, double *covariance)
{
 TDimension width, height, depth;
 unsigned int i,j,k;
 grphic3d *mask1=NULL, *mask2=NULL;
 double covar=0.0;
 int taille_covar=0;
 double rcoeff1,rcoeff2;

 mask1 = im1->mask;
 mask2 = im2->mask;
 if (!mask1 || !mask2)
 { fprintf (stderr, "les images n'ont pas de masque dans imx_calc_covar_im_mask_3d_p\n"); *covariance=0.0; return 0;}

 rcoeff1=im1->rcoeff; rcoeff2=im2->rcoeff;
 height=im1->height; width=im1->width; depth=im1->depth;

 for(i=0;i<width;i++)
  for(j=0;j<height;j++)
   for(k=0;k<depth;k++)
   {
    if (mask2->mri[i][j][k]!=0 && mask1->mri[i][j][k]!=0)
    {
      covar +=(rcoeff1*im1->mri[i][j][k]-moyenne1)*(rcoeff2*im2->mri[i][j][k]-moyenne2);
      taille_covar++;
    }
   }

 if (taille_covar>0) covar/=taille_covar;
 *covariance=covar;

 return taille_covar;
}

int imx_norm_lin_im_3d_p(grphic3d *im1, grphic3d *imres, double alpha, double beta)
{
 int i,j,k;
 float max, min;
 int width, height, depth;
 double rcoeff1,rcoeff;
 int err;
 double alpha2=0.0, beta2=0.0;
 TYPEMRI3D val1=0;
 TYPEMRI3D ***im1MRI=NULL, ***imresMRI=NULL;

 //Calcul de imres->maxpixel et imres->icomp et imres->rcoeff a partir de max et min
 rcoeff1=im1->rcoeff;
 imx_copie_param_3d_p(im1,imres);
 max=(float)(alpha*(rcoeff1*im1->max_pixel)+beta);
 min=(float)(alpha*(rcoeff1*im1->min_pixel)+beta);
 err=imx_brukermax_3d(max,min,imres);
 rcoeff=imres->rcoeff;

 //normalisation
 height=(int)im1->height; width=(int)im1->width; depth=(int)im1->depth;
 alpha2=alpha*rcoeff1/rcoeff; beta2=beta/rcoeff+0.5;
 im1MRI=im1->mri; imresMRI=imres->mri;
 for (i=0;i<width;i++)
  for (j=0;j<height;j++)
   for (k=0;k<depth;k++)
   {
    val1=im1MRI[i][j][k];
    if (val1==0) imresMRI[i][j][k]=0;
    else imresMRI[i][j][k]=(TYPEMRI3D)floor(alpha2*val1+beta2);
   }

 return err;
}

/****   imx_norm_im_3d_p    ****************************************
****
***  NORMALISATION D'UNE IMAGE PAR RAPPORT A UNE AUTRE (meme min et max)
***
**************************************************************/
void imx_norm_im_minmax_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{
  int i,j,k,err,hght,wdth,dpth;
  int max1,min1;
  int max2,min2;
  float max,min,rcoeff;
  double coeffA,coeffB;

  imx_copie_param_3d_p(im1,imres);

  dpth=im1->depth;
  hght=im1->height;
  wdth=im1->width;
  max1=0;
  min1=0;
  for(i=0;i<wdth;i++)
  for(j=0;j<hght;j++)
  for(k=0;k<dpth;k++)
  {
    max1 = MAXI(im1->mri[i][j][k],max1);
    min1 = MINI(im1->mri[i][j][k],min1);
  }

  dpth=im2->depth;
  hght=im2->height;
  wdth=im2->width;
  max2=0;
  min2=0;
  for(i=0;i<wdth;i++)
  for(j=0;j<hght;j++)
  for(k=0;k<dpth;k++)
  {
    max2 = MAXI(im2->mri[i][j][k],max2);
    min2 = MINI(im2->mri[i][j][k],min2);
  }
  
  coeffA=((max2-min2)*im2->rcoeff) / ((max1-min1)*im1->rcoeff);
  coeffB=min2*im2->rcoeff - min1*im1->rcoeff*coeffA;
  
  imx_copie_param_3d_p(im1,imres);
  max=(float)(max2*im2->rcoeff);
  min=(float)(min2*im2->rcoeff);
  err=imx_brukermax_3d(max,min,imres);
  rcoeff=imres->rcoeff;

  dpth=im1->depth;
  hght=im1->height;
  wdth=im1->width;
  for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
  for (k=0;k<dpth;k++)
      imres->mri[i][j][k]=(TYPEMRI3D)floor(((im1->mri[i][j][k]*im1->rcoeff)*coeffA+coeffB)/rcoeff);
}


/****   imx_norm_im_3d_p    ****************************************
****
***  NORMALISATION D'UNE IMAGE PAR RAPPORT A UNE AUTRE NON SEUILLEE
***
**************************************************************/
void imx_norm_im_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{
	int i,j,k,err,hght,wdth,dpth;
	double total1,total2;
	double aux1,aux2;	
	float max,min,rcoeff;
	
	imx_copie_param_3d_p(im1,imres);

  	total1=0.;
	total2=0.;
	
	dpth=im1->depth;
  	hght=im1->height;
    wdth=im1->width;
    	
    for(i=0;i<wdth;i++)
  	  for(j=0;j<hght;j++)
	    for(k=0;k<dpth;k++)
		total1=total1+im1->mri[i][j][k];
	
	
    for(i=0;i<wdth;i++)
  	  for(j=0;j<hght;j++)
	    for(k=0;k<dpth;k++)
		total2=total2+im2->mri[i][j][k];

        imx_copie_param_3d_p(im1,imres);
        aux1=total1*im1->rcoeff;
        aux2=total2*im2->rcoeff;
        max=(float)(aux2*im1->rcoeff*im1->max_pixel/aux1);
        min=(float)(aux2*im1->rcoeff*im1->min_pixel/aux1);
        err=imx_brukermax_3d(max,min,imres);
        rcoeff=imres->rcoeff;

	if (total1!=0)
          {
	  for (i=0;i<wdth;i++)
            for (j=0;j<hght;j++)
               for (k=0;k<dpth;k++)
      	          imres->mri[i][j][k]=(TYPEMRI3D)floor(aux2*im1->rcoeff*im1->mri[i][j][k]/(aux1*rcoeff));
           }	
	else
	  PUT_ERROR("Impossible to normalize");


}

/****   imx_norm_im_meanecty_3d_p    ****************************************
****
***  NORMALISATION D'UNE IMAGE PAR RAPPORT A UNE AUTRE NON SEUILLEE
***
**************************************************************/
void imx_norm_im_meanecty_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{
int i,j,k;
int height,width,depth;
double moy1,moy2,var1,var2,sygma1,sygma2;
double rcoeff1,rcoeff2,rcoeff;
double aux1,aux2;
float taille;
float max,min;
int err;

	
height=im1->height;
width=im1->width;
depth=im1->depth;
taille=(float)(height*width*depth);
rcoeff1=im1->rcoeff;
rcoeff2=im2->rcoeff;
moy1=0;
moy2=0;
var1=0;
var2=0;
aux1=0;
aux2=0;
max=0;

for (i=0;i<height;i++)
   for (j=0;j<width;j++)
      for (k=0;k<depth;k++)
	{
	moy1 += im1->mri[i][j][k];
	moy2 += im2->mri[i][j][k];
	}
moy1= rcoeff1*moy1/taille;
moy2= rcoeff2*moy2/taille;
	
	
for (i=0;i<height;i++)
   for (j=0;j<width;j++)
      for (k=0;k<depth;k++)
      	{
var1 +=(rcoeff1*im1->mri[i][j][k]-moy1)*(rcoeff1*im1->mri[i][j][k]-moy1);
var2 +=(rcoeff2*im2->mri[i][j][k]-moy2)*(rcoeff2*im2->mri[i][j][k]-moy2);
	}

var1= var1/taille;
var2=var2/taille;
sygma1=sqrt(var1);
sygma2=sqrt(var2);
	
aux1=rcoeff1*sygma2/sygma1;
aux2=moy2-moy1*sygma2/sygma1;

imx_copie_param_3d_p(im1,imres);
max=(float)(aux1*(im1->max_pixel)+aux2);
min=(float)(aux1*(im1->min_pixel)+aux2);


/*   Calcul de imres->max_pixel et imres->icomp imres->rcoeff a partir de max*/
err= imx_brukermax_3d(max,min,imres);
rcoeff= imres->rcoeff;

for (i=0;i<width;i++)
   for (j=0;j<height;j++)
      for (k=0;k<depth;k++)
      {
      	imres->mri[i][j][k] = (TYPEMRI3D)floor((aux1*(im1->mri[i][j][k])+aux2)/rcoeff);
	 
      }
imres->min_pixel=0;
imres->width=width;
imres->height=height;
imres->depth=depth;

/* Mettre le meme rcoeff pour les deux images */
imx_norm_rcoeff_3d_p(im2,imres);	

return ;	
}

/****   imx_norm_im_mean_3d_p    ****************************************
****
***  NORMALISATION D'UNE IMAGE PAR RAPPORT A UNE AUTRE NON SEUILLEE
***  Normalisation de la moyenne
***
**************************************************************/
void imx_norm_im_mean_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{
int i,j,k;
int height,width,depth;
double moy1,moy2;
double rcoeff1,rcoeff2,rcoeff;
double aux1,aux2;
float taille;
float max,min;
int err;

	
height=im1->height;
width=im1->width;
depth=im1->depth;
taille=(float)(height*width*depth);
rcoeff1=im1->rcoeff;
rcoeff2=im2->rcoeff;
moy1=0;
moy2=0;
aux1=0;
aux2=0;
max=0;

for (i=0;i<height;i++)
   for (j=0;j<width;j++)
      for (k=0;k<depth;k++)
	{
	moy1 += im1->mri[i][j][k];
	moy2 += im2->mri[i][j][k];
	}
moy1= rcoeff1*moy1/taille;
moy2= rcoeff2*moy2/taille;
	
	
aux1=rcoeff1;
aux2=moy2-moy1;

imx_copie_param_3d_p(im1,imres);
max=(float)(aux1*(im1->max_pixel)+aux2);
min=(float)(aux1*(im1->min_pixel)+aux2);


/*   Calcul de imres->max_pixel et imres->icomp imres->rcoeff a partir de max*/
err= imx_brukermax_3d(max,min,imres);
rcoeff= imres->rcoeff;

for (i=0;i<width;i++)
   for (j=0;j<height;j++)
      for (k=0;k<depth;k++)
      {
      	imres->mri[i][j][k] = (TYPEMRI3D)floor((aux1*(im1->mri[i][j][k])+aux2)/rcoeff);
	 
      }
imres->min_pixel=0;
imres->width=width;
imres->height=height;
imres->depth=depth;

/* Mettre le meme rcoeff pour les deux images */
imx_norm_rcoeff_3d_p(im2,imres);	

return ;	
}

/***    imx_norm_im_rl_3d_p()   ------------------------------------------------------
***	NORMALISATION D'UNE IMAGE 3D PAR RAPPORT A UNE AUTRE
***
***    normalisation en utilisant la regression lineaire
***
****************************************************************/
void imx_norm_im_rl_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{
int i,j,k;
int height,width,depth;
float taille;
double moy1,moy2,var1,var2,covar,sygma1,sygma2,rho;
double rcoeff1,rcoeff2,rcoeff;
double aux1,aux2;
float max, min;
int err;

moy1=0.0;
moy2=0.0;
var1=0.0;
var2=0.0;
covar=0.0;

height=im1->height;
width=im1->width;
depth=im1->depth;
taille=(float)(height*width*depth);
rcoeff1=im1->rcoeff;
rcoeff2=im2->rcoeff;

/*cacul des moyennes des images*/
for(i=0;i<width;i++)
    for(j=0;j<height;j++)
	for(k=0;k<depth;k++)
	{
		moy1 += im1->mri[i][j][k];
		moy2 += im2->mri[i][j][k];
	} 

moy1 = rcoeff1*moy1/taille;
moy2 = rcoeff2*moy2/taille;

/*calcul des variances des images*/
for(i=0;i<width;i++)
   for(j=0;j<height;j++)
	for(k=0;k<depth;k++)
	{
	  var1 +=(rcoeff1*im1->mri[i][j][k]-moy1)*(rcoeff1*im1->mri[i][j][k]-moy1);
	  var2 +=(rcoeff2*im2->mri[i][j][k]-moy2)*(rcoeff2*im2->mri[i][j][k]-moy2);
	  covar +=(rcoeff1*im1->mri[i][j][k]-moy1)*(rcoeff2*im2->mri[i][j][k]-moy2);
	}
	

var1 = var1/taille;
var2 = var2/taille;
covar = covar/taille;
sygma1 = sqrt(var1);
sygma2 = sqrt(var2);
rho = covar/(sygma1*sygma2);

aux1=(rho*sygma2)/sygma1;
aux2=moy2-((moy1*rho*sygma2)/sygma1);


/*  Calcul de imres->maxpixel et imres->icomp et imres->rcoeff a partir de max*/
imx_copie_param_3d_p(im1,imres);
max=(float)(aux1*(rcoeff1*im1->max_pixel)+aux2);
min=(float)(aux1*(rcoeff1*im1->min_pixel)+aux2);
err=imx_brukermax_3d(max,min,imres);
rcoeff=imres->rcoeff;


for (i=0;i<width;i++)
   for (j=0;j<height;j++)
	for(k=0;k<depth;k++)
	    imres->mri[i][j][k]=(TYPEMRI3D)floor((aux1*rcoeff1*im1->mri[i][j][k]+aux2)/rcoeff);

imres->min_pixel=0;
imres->width=width;
imres->height=height;
imres->depth=depth;

return;
}

/************************************************************
*** imx_norm_seuil_3d_p()
***	NORMALISATION D'UNE IMAGE 3D PAR RAPPORT A UNE IMAGE 
***	SEUILLE
***************************************************************/
void imx_norm_seuil_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{
	int i,j,k,hght,wdth,dpth,err;
	double total1,total2,coef,rcoeff;
    double aux1,aux2;
	float min,max;
    int m,l;
	
	total1=0.;
	total2=0.;
	m=0;
	l=0;
	
	dpth=im1->depth;
  	hght=im1->height;
    	wdth=im1->width;
    	
    for(i=0;i<wdth;i++)
  	  for(j=0;j<hght;j++)
	    for(k=0;k<dpth;k++)    		{
    		 	if (im1->mri[i][j][k]!=0)
			{
				m=m+1;
				total1=total1+im1->mri[i][j][k];
			}
		}		
	
    for(i=0;i<wdth;i++)
  	  for(j=0;j<hght;j++)
	    for(k=0;k<dpth;k++)
    		{
    			if (im2->mri[i][j][k]!=0)
			{
				l=l+1;
				total2=total2+im2->mri[i][j][k];
			}
		}		

        imx_copie_param_3d_p(im1,imres);
        aux1=total1*im1->rcoeff/m;
        aux2=total2*im2->rcoeff/l;
        max=(float)(aux2*im1->rcoeff*im1->max_pixel/aux1);
        min=(float)(aux2*im1->rcoeff*im1->min_pixel/aux1);
        err=imx_brukermax_3d(max,min,imres);
        rcoeff=imres->rcoeff;

	if (total1!=0)
          {
	  for (i=0;i<wdth;i++)
            for (j=0;j<hght;j++)
               for (k=0;k<dpth;k++)
      	          imres->mri[i][j][k]=(TYPEMRI3D)floor(aux2*im1->rcoeff*im1->mri[i][j][k]/(aux1*rcoeff));
           }	
	else
	  PUT_ERROR("Impossible to normalize");
	


	
	if (total1!=0)
	{
		coef=(total2*im2->rcoeff*m)/(total1*im1->rcoeff*l);
		imx_mul_coe_3d_p(im1,(float)coef,imres);
		return;
	}
	else
	{
		PUTI("Erreur : Veuillez verifier l'existance\n    de l'image ",1,"_\n");
		return;
	}
		
}


/************************************************************
*** imx_norm_seuil_3d_p()
***	NORMALISATION D'UNE IMAGE 3D PAR RAPPORT A UNE IMAGE
***	SEUILLE
***************************************************************/
void imx_norm_seuil_minmax_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{
  int i,j,k,err,hght,wdth,dpth;
  int max1,min1;
  int max2,min2;
  float max,min,rcoeff;
  double coeffA,coeffB;
  int total1,total2;
  imx_copie_param_3d_p(im1,imres);

  dpth=im1->depth;
  hght=im1->height;
  wdth=im1->width;
  max1=0;
  min1=0;
  total1=0;

  for(i=0;i<wdth;i++)
  for(j=0;j<hght;j++)
  for(k=0;k<dpth;k++)
  {
    if (im1->mri[i][j][k])
    {
      max1 = MAXI(im1->mri[i][j][k],max1);
      min1 = MINI(im1->mri[i][j][k],min1);
      total1++;
    }
  }

  dpth=im2->depth;
  hght=im2->height;
  wdth=im2->width;
  max2=0;
  min2=0;
  total2=0;
  
  for(i=0;i<wdth;i++)
  for(j=0;j<hght;j++)
  for(k=0;k<dpth;k++)
  {
    if (im1->mri[i][j][k])
    {
      max2 = MAXI(im2->mri[i][j][k],max2);
      min2 = MINI(im2->mri[i][j][k],min2);
      total2++;
    }
  }

  if (total1 && total2)
  {
    coeffA=((max2-min2)*im2->rcoeff) / ((max1-min1)*im1->rcoeff);
    coeffB=min2*im2->rcoeff - min1*im1->rcoeff*coeffA;

    imx_copie_param_3d_p(im1,imres);
    max=(float)(max2*im2->rcoeff);
    min=(float)(min2*im2->rcoeff);
    err=imx_brukermax_3d(max,min,imres);
    rcoeff=imres->rcoeff;

    dpth=im1->depth;
    hght=im1->height;
    wdth=im1->width;
    for (i=0;i<wdth;i++)
    for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
      if (im1->mri[i][j][k])
        imres->mri[i][j][k]=(TYPEMRI3D)floor(((im1->mri[i][j][k]*im1->rcoeff)*coeffA+coeffB)/rcoeff);    
  }
  else
    PUT_ERROR("Impossible to normalize\n");

}
/****   imx_norm_seuil_3d_p    ****************************************
****
***  NORMALISATION D'UNE IMAGE PAR RAPPORT A UNE AUTRE NON SEUILLEE
***
**************************************************************/
void imx_norm_seuiltoto_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{
	int i,j,k,err,hght1,wdth1,dpth1,hght2,wdth2,dpth2;
	double total1,total2;
	double aux1,aux2;	
	float max,min,rcoeff;
	/*  A faire */
	imx_copie_param_3d_p(im1,imres);

  	total1=0.;
	total2=0.;
	
	dpth1=im1->depth;
  	hght1=im1->height;
    wdth1=im1->width;
    	
    dpth2=im2->depth;
  	hght2=im2->height;
    wdth2=im2->width;
    		
    for(i=0;i<wdth1;i++)
  	  for(j=0;j<hght1;j++)
	    for(k=0;k<dpth1;k++)
		total1=total1+im1->mri[i][j][k];
	
	
    for(i=0;i<wdth2;i++)
  	  for(j=0;j<hght2;j++)
	    for(k=0;k<dpth2;k++)
		total2=total2+im2->mri[i][j][k];

    aux1=total1*im1->rcoeff/(dpth1*wdth1*hght1);
    aux2=total2*im2->rcoeff/(dpth2*wdth2*hght2);
    imx_copie_param_3d_p(im1,imres);
    max=(float)(aux2*im1->max_pixel/aux1);
    min=(float)(aux2*im1->min_pixel/aux1);
    err=imx_brukermax_3d(max,min,imres);
    rcoeff=imres->rcoeff;

	if (total1!=0)
          {
	  for (i=0;i<wdth1;i++)
            for (j=0;j<hght1;j++)
               for (k=0;k<dpth1;k++)
      	          imres->mri[i][j][k]=(TYPEMRI3D)floor(aux2*im1->mri[i][j][k]/(aux1*rcoeff));
           }	
	else
	  PUT_ERROR("Impossible to normalize");
	
}

/******************************************************************
*** imx_norm_seuil_meanecty_3d_p 
**/
/*!  NORMALISATION D'UNE IMAGE (im1) PAR RAPPORT A UNE AUTRE (im2) SEUILLEE
***
***	\param im1 : image source
***	\param im2 : image ref
***	\param imres : image resultat (E/S)
**************************************************************/
void imx_norm_seuil_meanecty_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{
int i,j,k,l,m;
int height,width,depth;
double moy1,moy2,var1,var2,sygma1,sygma2;
double rcoeff1,rcoeff2,rcoeff;
double aux1,aux2;
float max,min;
int err;



m=0;
l=0;
moy1=0;
moy2=0;
var1=0;
var2=0;
height=im1->height;
width=im1->width;
depth=im1->depth;
rcoeff1=im1->rcoeff;
rcoeff2=im2->rcoeff;

/*cacul des moyennes des images*/	
for(i=0;i<width;i++)
  for(j=0;j<height;j++)
     for(k=0;k<depth;k++)
     {
      if (im1->mri[i][j][k])
      {
        l++;
        moy1 += im1->mri[i][j][k];
      }
     }		
height=im2->height;
width=im2->width;
depth=im2->depth;
for(i=0;i<width;i++)
  for(j=0;j<height;j++)
     for(k=0;k<depth;k++)
     {
      if (im2->mri[i][j][k])
      {
        m++;
        moy2 += im2->mri[i][j][k];
      }
     }

moy1= rcoeff1*moy1/l;
moy2= rcoeff2*moy2/m;

/*calcul des variances des images*/
height=im1->height;
width=im1->width;
depth=im1->depth;
for(i=0;i<width;i++)
  for(j=0;j<height;j++)
     for(k=0;k<depth;k++)
      {
        if (im1->mri[i][j][k])
          var1 +=(rcoeff1*im1->mri[i][j][k]-moy1)*(rcoeff1*im1->mri[i][j][k]-moy1);
      }

height=im2->height;
width=im2->width;
depth=im2->depth;
for(i=0;i<width;i++)
  for(j=0;j<height;j++)
     for(k=0;k<depth;k++)
      {
        if (im2->mri[i][j][k])
          var2 +=(rcoeff2*im2->mri[i][j][k]-moy2)*(rcoeff2*im2->mri[i][j][k]-moy2);
      }

var1= var1/l;
var2=var2/m;
sygma1=sqrt(var1);
sygma2=sqrt(var2);
	
aux1=rcoeff1*sygma2/sygma1;
aux2=moy2-moy1*sygma2/sygma1;

imx_copie_param_3d_p(im1,imres);
max=(float)(aux1*(im1->max_pixel)+aux2);
min=(float)(aux1*(im1->min_pixel)+aux2);

/*   Calcul de imres->maxpixel et imres->icomp et imres->rcoeff a partir de max*/
err=imx_brukermax_3d(max,min,imres);
rcoeff=imres->rcoeff;

  
height=im1->height;
width=im1->width;
depth=im1->depth;
for (i=0;i<width;i++)
   for (j=0;j<height;j++)
      for(k=0;k<depth;k++)
       if (im1->mri[i][j][k]!=0)
        {imres->mri[i][j][k]=(TYPEMRI3D)floor((aux1*im1->mri[i][j][k]+aux2)/rcoeff);}
       else
        {imres->mri[i][j][k]=0;}
	
/* Modif JP 2/2001 imres->min_pixel=0;*/
imres->width=width;
imres->height=height;
imres->depth=depth;

/* Mettre le meme rcoeff pour les deux images */
imx_norm_rcoeff_3d_p(im2,imres);
return;	
}

/****   imx_norm_seuil_mean_3d_p    ****************************************
****
***  NORMALISATION D'UNE IMAGE PAR RAPPORT A UNE AUTRE NON SEUILLEE
***
**************************************************************/
void imx_norm_seuil_mean_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{
int i,j,k,l,m;
int height,width,depth;
double moy1,moy2;
double rcoeff1,rcoeff2,rcoeff;
double aux1,aux2;
float max,min;
int err;


m=0;
l=0;
moy1=0;
moy2=0;
height=im1->height;
width=im1->width;
depth=im1->depth;
rcoeff1=im1->rcoeff;
rcoeff2=im2->rcoeff;

/*cacul des moyennes des images*/	
for(i=0;i<width;i++)
  for(j=0;j<height;j++)
     for(k=0;k<depth;k++)
    	{
    	if (im1->mri[i][j][k]!=0)
		{
		l++;
		moy1 += im1->mri[i][j][k];
		}
	if (im2->mri[i][j][k]!=0)
		{
		m++;
		moy2 += im2->mri[i][j][k];
		}
	}		
	
moy1= rcoeff1*moy1/l;
moy2= rcoeff2*moy2/m;

aux1=rcoeff1;
aux2=moy2-moy1;

imx_copie_param_3d_p(im1,imres);
max=(float)(aux1*(im1->max_pixel)+aux2);
min=(float)(aux1*(im1->min_pixel)+aux2);

/*   Calcul de imres->maxpixel et imres->icomp et imres->rcoeff a partir de max*/
  err=imx_brukermax_3d(max,min,imres);
  rcoeff=imres->rcoeff;

  
for (i=0;i<width;i++)
   for (j=0;j<height;j++)
      for(k=0;k<depth;k++)
       if (im1->mri[i][j][k]!=0)
        {imres->mri[i][j][k]=(TYPEMRI3D)floor((aux1*im1->mri[i][j][k]+aux2)/rcoeff);}
       else
        {imres->mri[i][j][k]=0;}
	
imres->min_pixel=0;
imres->width=width;
imres->height=height;
imres->depth=depth;

/* Mettre le meme rcoeff pour les deux images */
imx_norm_rcoeff_3d_p(im2,imres);
return;	
}

/***    imx_norm_im_rl_3d_p()   ------------------------------------------------------
***	NORMALISATION D'UNE IMAGE 3D PAR RAPPORT A UNE AUTRE
***
***    normalisation en utilisant la regression lineaire
***
****************************************************************/
void imx_norm_seuil_rl_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{
int i,j,k,l,m,lm;
int height,width,depth;
double moy1,moy2,var1,var2,covar,sygma1,sygma2,rho;
double rcoeff1,rcoeff2,rcoeff;
double aux1,aux2;
float max, min;
int err;

l=0;
m=0;
lm=0;
moy1=0.0;
moy2=0.0;
var1=0.0;
var2=0.0;
covar=0.0;

height=im1->height;
width=im1->width;
depth=im1->depth;
rcoeff1=im1->rcoeff;
rcoeff2=im2->rcoeff;

/*cacul des moyennes des images*/
for(i=0;i<width;i++)
  for(j=0;j<height;j++)
	for(k=0;k<depth;k++)
	{
            if (im1->mri[i][j][k]!=0)
	    {
		l++;
		moy1 += im1->mri[i][j][k];
	    }
	    if (im2->mri[i][j][k]!=0)
	    {
	        m++;
		moy2 += im2->mri[i][j][k];
	    }
	} 

moy1 = rcoeff1*moy1/l;
moy2 = rcoeff2*moy2/m;

/*calcul des variances des images*/
for(i=0;i<width;i++)
  for(j=0;j<height;j++)
	for(k=0;k<depth;k++)
	{	
            if (im1->mri[i][j][k]!=0)
		var1 +=(rcoeff1*im1->mri[i][j][k]-moy1)*(rcoeff1*im1->mri[i][j][k]-moy1);
	    if (im2->mri[i][j][k]!=0)
	       var2 +=(rcoeff2*im2->mri[i][j][k]-moy2)*(rcoeff2*im2->mri[i][j][k]-moy2);
	    if (im1->mri[i][j][k]!=0 && im2->mri[i][j][k]!=0)
	    {
	       lm++;
	       covar +=(rcoeff1*im1->mri[i][j][k]-moy1)*(rcoeff2*im2->mri[i][j][k]-moy2);
	    }
	}

var1 = var1/l;
var2 = var2/m;
covar = covar/lm;
sygma1 = sqrt(var1);
sygma2 = sqrt(var2);
rho = covar/(sygma1*sygma2);

aux1=(rho*sygma2)/sygma1;
aux2=moy2-((moy1*rho*sygma2)/sygma1);


/*  Calcul de imres->maxpixel et imres->icomp et imres->rcoeff a partir de max*/
imx_copie_param_3d_p(im1,imres);
max=(float)(aux1*(rcoeff1*im1->max_pixel)+aux2);
min=(float)(aux1*(rcoeff1*im1->min_pixel)+aux2);
err=imx_brukermax_3d(max,min,imres);
rcoeff=imres->rcoeff;



for (i=0;i<width;i++)
  for (j=0;j<height;j++)
	for(k=0;k<depth;k++)
	if (im1->mri[i][j][k]!=0)
	    {imres->mri[i][j][k]=(TYPEMRI3D)floor((aux1*rcoeff1*im1->mri[i][j][k]+aux2)/rcoeff);}
	else
	    {imres->mri[i][j][k]=0;}


imres->min_pixel=0;
imres->width=width;
imres->height=height;
imres->depth=depth;

return;
}

/************************ norm_local_3d() ***********************************/ 
/*	Normalisation d'une image 3D par rapport �une autre image seuill� */ 
/*	Normalisation des moyennes et �arts types                          */ 
/*	im_1	:image �normaliser                                         */ 
/*	im_2	:image r��ence                                            */
/*	im_res	:image resultat                                             */
/*                                                                          */
/****************************************************************************/
void norm_local_3d(void)
{
	int im_1,im_2,im_res;
	int N,err=0;

	im_1=GET_PLACE3D(TEXT0231);
	im_2=GET_PLACE3D(TEXT0031);
	im_res=GET_PLACE3D(TEXT0006);
	
	while((N= GET_INT(TEXT0187, 0, &err))<=0)
          PUT_WNDMSG(TEXT0189);
        CLEAR_WNDMSG();

	imx_norm_local_3d(im_1,im_2,im_res,N);
	
	show_picture_3d(im_res);
}

/****************************************************************************/
void imx_norm_local_3d(int im_1, int im_2, int im_res, int N)
{	
	grphic3d *im1,*im2,*imres;
	
	im1=ptr_img_3d(im_1);	
	im2=ptr_img_3d(im_2);
	imres=ptr_img_3d(im_res);
	
    imx_norm_local_3d_p(im1,im2,imres,N);
}



/****************************************************************************/
int imx_norm_local_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres, int N)
{

float *mri_temp;
int width,height,depth,width2,height2,depth2;
int i,j,k,f,g,h,w,l,m;
int inf1,inf2,inf3,sup1,sup2,sup3;
int err;
int N3;
float sum1,sum2,sum1_square,sum2_square;
float moy1,moy2,var1,var2,sygma1,sygma2;
float rcoeff1,rcoeff2,rcoeff;
float aux1,aux2;
float norm;
float max=0,min=LONGMAX;

width=im1->width;
height=im1->height;
depth=im1->depth;
w=N/2;
width2=width-w;
height2=height-w;
depth2=depth-w;
rcoeff1=im1->rcoeff;
rcoeff2=im2->rcoeff;
N3=N*N*N;

imx_copie_param_3d_p(im1,imres);

if (( mri_temp = (float *)malloc(sizeof(double)*width*height*depth))==NULL)
{
  printf("failled to allocate memory\n");
  return (-1);
}

for (i=0;i<width;i++)
 for (j=0;j<height;j++)
  for (k=0;k<depth;k++)
    mri_temp[i*height*depth+j*depth+k]=0.0;

for(i=w;i<width2;i++)
 for(j=w;j<height2;j++)
  for (k=w;k<depth2;k++)
  {
    if (im1->mri[i][j][k]==0) norm=0;
    else{
    l=0;
    m=0;
    moy1=0.0;
    moy2=0.0;
    var1=0.0;
    var2=0.0;
    sum1=0.0;
    sum2=0.0;
    sum1_square=0.0;
    sum2_square=0.0;
    inf1=i-w;
    sup1=i+w;
    inf2=j-w;
    sup2=j+w;
    inf3=k-w;
    sup3=k+w;
    for(f=inf1;f<=sup1;f++)
     for(g=inf2;g<=sup2;g++)
      for(h=inf3;h<=sup3;h++)
      {
        if (im1->mri[f][g][h]!=0)
	{
	  l++;
	  sum1 += rcoeff1*im1->mri[f][g][h];
          sum1_square += (rcoeff1*im1->mri[f][g][h]*rcoeff1*im1->mri[f][g][h]);
	}
	if (im2->mri[f][g][h]!=0)
	{
	  m++;
          sum2 += rcoeff2*im2->mri[f][g][h];
	  sum2_square += (rcoeff2*im2->mri[f][g][h]*rcoeff2*im2->mri[f][g][h]);
	}
      }	
    
    moy1=sum1/l;
    moy2=sum2/m;
    var1 = sum1_square/l-moy1*moy1;
    var2 = sum2_square/m-moy2*moy2;
    sygma1 = (float)sqrt(var1);
    sygma2 = (float)sqrt(var2);
    
    aux1=sygma2/sygma1;
    aux2=moy2-((moy1*sygma2)/sygma1);
    norm=aux1*rcoeff1*im1->mri[i][j][k]+aux2;
    }

    if((mri_temp[i*height2*depth2+j*depth2+k]=norm)>max) max=norm;
    min=MINI(min,norm);
    
  }

/* Calcul de imres->maxpixel et imres->icomp et imres->rcoeff a partir de max*/
err=imx_brukermax_3d(max,min,imres);
rcoeff=imres->rcoeff;

for (i=w;i<width2;i++)
 for (j=w;j<height2;j++)
  for (k=w;k<depth2;k++)
  imres->mri[i][j][k]=(TYPEMRI3D)floor(mri_temp[i*height2*depth2+j*depth2+k]/rcoeff);


free(mri_temp);

return(1);
}




/************************ norm_coronal_3d() *********************************/ 
/*	Normalisation coupe par coupe ( coupes coronles ) d'une             */
/*             image 3D par rapport �une autre image seuill�              */
/*									    */
/*		Normalisation des moyennes et �arts types                  */ 
/*									    */
/*			im_1	:image �normaliser                         */ 
/*			im_2	:image r��ence                            */
/*			im_res	:image resultat                             */
/*                                                                          */
/****************************************************************************/
void norm_coronal_3d(void)
{
	int im_1,im_2,im_res;

	im_1=GET_PLACE3D(TEXT0231);
	im_2=GET_PLACE3D(TEXT0031);
	im_res=GET_PLACE3D(TEXT0006);
	
	imx_norm_coronal_3d(im_1,im_2,im_res);
	
	show_picture_3d(im_res);
}

/****************************************************************************/

void imx_norm_coronal_3d(int im_1, int im_2, int im_res)
{	
	grphic3d *im1,*im2,*imres;
	
	im1=ptr_img_3d(im_1);	
	im2=ptr_img_3d(im_2);
	imres=ptr_img_3d(im_res);
	
	imx_norm_coronal_3d_p(im1,im2,imres); 
}



/****************************************************************************/

void imx_norm_coronal_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{
int i,j,k,l,m;
int height,width,depth;
float *mri_temp;
float moy1,moy2,var1,var2,sygma1,sygma2;
float rcoeff1,rcoeff2,rcoeff;
float max=0,min=LONGMAX;
float aux1,aux2,norm;
int err;


height=im1->height;
width=im1->width;
depth=im1->depth;
rcoeff1=im1->rcoeff;
rcoeff2=im2->rcoeff;


if (( mri_temp = (float *)malloc(sizeof(double)*width*height*depth))==NULL)
{
  printf("failled to allocate memory\n");
  return ;
}

imx_copie_param_3d_p(im1,imres);

for (i=0;i<width;i++)
 for (j=0;j<height;j++)
  for (k=0;k<depth;k++)
  {
    mri_temp[i*height*depth+j*depth+k]=0.0;
    imres->mri[i][j][k]=0;
  }

for (k=0;k<depth;k++)   /* une coupe coronale */
{
 m=0;
 l=0;
 moy1=0;
 moy2=0;
 var1=0;
 var2=0;

 /*cacul des moyennes des images*/
 for(i=0;i<width;i++)
   for(j=0;j<height;j++)
   {
    if (im1->mri[i][j][k]!=0)
    {
     l++;
     moy1 += im1->mri[i][j][k];
    }
    if (im2->mri[i][j][k]!=0)
    {
     m++;
     moy2 += im2->mri[i][j][k];
    }
   }		
	
 moy1= rcoeff1*moy1/l;
 moy2= rcoeff2*moy2/m;

 /*calcul des variances des images*/
 for(i=0;i<width;i++)
   for(j=0;j<height;j++)
   {
   if (im1->mri[i][j][k]!=0)
    var1 +=(rcoeff1*im1->mri[i][j][k]-moy1)*(rcoeff1*im1->mri[i][j][k]-moy1);
   if (im2->mri[i][j][k]!=0)
    var2 +=(rcoeff2*im2->mri[i][j][k]-moy2)*(rcoeff2*im2->mri[i][j][k]-moy2);
   }

 var1= var1/l;
 var2=var2/m;
 sygma1=(float)sqrt(var1);
 sygma2=(float)sqrt(var2);
	
 for(i=0;i<width;i++)
   for(j=0;j<height;j++)
   if (im1->mri[i][j][k]==0) norm=0;
   else
   { 
    aux1=(sygma2)/sygma1;
    aux2=moy2-((moy1*sygma2)/sygma1);
    norm=aux1*rcoeff1*im1->mri[i][j][k]+aux2;
    if((mri_temp[i*height*depth+j*depth+k]=norm)>max) max=norm;
    min=MINI(min,norm);
   }

}

/* Calcul de imres->maxpixel et imres->icomp et imres->rcoeff a partir de max*/
err=imx_brukermax_3d(max,min,imres);
rcoeff=imres->rcoeff;

for (i=0;i<width;i++)
 for (j=0;j<height;j++)
  for(k=0;k<depth;k++)
  if (im1->mri[i][j][k]!=0)
   imres->mri[i][j][k]=(TYPEMRI3D)floor(mri_temp[i*height*depth+j*depth+k]/rcoeff);
  else
   imres->mri[i][j][k]=0;

imres->min_pixel=0;

free(mri_temp);

return;	
}


/************************ norm_sagittal_3d() ********************************/ 
/*	Normalisation coupe par coupe ( coupes sagittales ) d'une image 3D  */
/*                 par rapport �une autre image seuill�                   */
/*									    */
/*		Normalisation des moyennes et �arts types                  */ 
/*									    */
/*			im_1	:image �normaliser                         */ 
/*			im_2	:image r��ence                            */
/*			im_res	:image resultat                             */
/*                                                                          */
/****************************************************************************/
void norm_sagittal_3d(void)
{
	int im_1,im_2,im_res;

	im_1=GET_PLACE3D(TEXT0231);
	im_2=GET_PLACE3D(TEXT0031);
	im_res=GET_PLACE3D(TEXT0006);
	
	imx_norm_sagittal_3d(im_1,im_2,im_res);
	
	show_picture_3d(im_res);
}

/****************************************************************************/

void imx_norm_sagittal_3d(int im_1, int im_2, int im_res)
{	
	grphic3d *im1,*im2,*imres;
	
	im1=ptr_img_3d(im_1);	
	im2=ptr_img_3d(im_2);
	imres=ptr_img_3d(im_res);
	
	imx_norm_sagittal_3d_p(im1,im2,imres); 
}



/****************************************************************************/

void imx_norm_sagittal_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{
int i,j,k,l,m;
int height,width,depth;
float *mri_temp;
float moy1,moy2,var1,var2,sygma1,sygma2;
float rcoeff1,rcoeff2,rcoeff;
float max=0,min=LONGMAX;
float aux1,aux2,norm;
int err;


height=im1->height;
width=im1->width;
depth=im1->depth;
rcoeff1=im1->rcoeff;
rcoeff2=im2->rcoeff;


if (( mri_temp = (float *)malloc(sizeof(float)*width*height*depth))==NULL)
{
  printf("failled to allocate memory\n");
  return ;
}

imx_copie_param_3d_p(im1,imres);

for (i=0;i<width;i++)
 for (j=0;j<height;j++)
  for (k=0;k<depth;k++)
  {
    mri_temp[i*height*depth+j*depth+k]=0.0;
    imres->mri[i][j][k]=0;
  }
 
for(i=0;i<width;i++)     /* une coupe sagittale */
{
 m=0;
 l=0;
 moy1=0;
 moy2=0;
 var1=0;
 var2=0;

 /*cacul des moyennes des images*/
 for(j=0;j<height;j++)
   for(k=0;k<depth;k++)	
   {
    if (im1->mri[i][j][k]!=0)
    {
     l++;
     moy1 += im1->mri[i][j][k];
    }
    if (im2->mri[i][j][k]!=0)
    {
     m++;
     moy2 += im2->mri[i][j][k];
    }
   }		
	
 moy1= rcoeff1*moy1/l;
 moy2= rcoeff2*moy2/m;

 /*calcul des variances des images*/
 for(j=0;j<height;j++)
   for(k=0;k<depth;k++)	
   {
   if (im1->mri[i][j][k]!=0)
    var1 +=(rcoeff1*im1->mri[i][j][k]-moy1)*(rcoeff1*im1->mri[i][j][k]-moy1);
   if (im2->mri[i][j][k]!=0)
    var2 +=(rcoeff2*im2->mri[i][j][k]-moy2)*(rcoeff2*im2->mri[i][j][k]-moy2);
   }

 var1= var1/l;
 var2=var2/m;
 sygma1=(float)sqrt(var1);
 sygma2=(float)sqrt(var2);
	
 for(j=0;j<height;j++)
   for(k=0;k<depth;k++)	
   if (im1->mri[i][j][k]==0) norm=0;
   else
   { 
    aux1=(sygma2)/sygma1;
    aux2=moy2-((moy1*sygma2)/sygma1);
    norm=aux1*rcoeff1*im1->mri[i][j][k]+aux2;
    if((mri_temp[i*height*depth+j*depth+k]=norm)>max) max=norm;
    min=MINI(min,norm);
  }

}

/* Calcul de imres->maxpixel et imres->icomp et imres->rcoeff a partir de max*/
err=imx_brukermax_3d(max,min,imres);
rcoeff=imres->rcoeff;

for (i=0;i<width;i++)
 for (j=0;j<height;j++)
  for(k=0;k<depth;k++)
  if (im1->mri[i][j][k]!=0)
   imres->mri[i][j][k]=(TYPEMRI3D)floor(mri_temp[i*height*depth+j*depth+k]/rcoeff);
  else
   imres->mri[i][j][k]=0;


free(mri_temp);

return;	
}


/************************ norm_transverse_3d() ******************************/ 
/*	Normalisation coupe par coupe ( coupes transverses ) d'une image 3D */
/*                 par rapport �une autre image seuill�                   */
/*									    */
/*		Normalisation des moyennes et �arts types                  */ 
/*									    */
/*			im_1	:image �normaliser                         */ 
/*			im_2	:image r��ence                            */
/*			im_res	:image resultat                             */
/*                                                                          */
/****************************************************************************/
void norm_transverse_3d(void)
{
	int im_1,im_2,im_res;

	im_1=GET_PLACE3D(TEXT0231);
	im_2=GET_PLACE3D(TEXT0031);
	im_res=GET_PLACE3D(TEXT0006);
	
	imx_norm_transverse_3d(im_1,im_2,im_res);
	
	show_picture_3d(im_res);
}

/****************************************************************************/

void imx_norm_transverse_3d(int im_1, int im_2, int im_res)
{	
	grphic3d *im1,*im2,*imres;
	
	im1=ptr_img_3d(im_1);	
	im2=ptr_img_3d(im_2);
	imres=ptr_img_3d(im_res);
	
	imx_norm_transverse_3d_p(im1,im2,imres); 
}



/****************************************************************************/

void imx_norm_transverse_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{
int i,j,k,l,m;
int height,width,depth;
float *mri_temp;
float moy1,moy2,var1,var2,sygma1,sygma2;
float rcoeff1,rcoeff2,rcoeff;
float max=0,min=LONGMAX;
float aux1,aux2,norm;
int err;


height=im1->height;
width=im1->width;
depth=im1->depth;
rcoeff1=im1->rcoeff;
rcoeff2=im2->rcoeff;


if (( mri_temp = (float *)malloc(sizeof(float)*width*height*depth))==NULL)
{
  printf("failled to allocate memory\n");
  return ;
}

imx_copie_param_3d_p(im1,imres);

for (i=0;i<width;i++)
 for (j=0;j<height;j++)
  for (k=0;k<depth;k++)
  {
    mri_temp[i*height*depth+j*depth+k]=0.0;
    imres->mri[i][j][k]=0;
  }

for(j=0;j<height;j++)	 /* une coupe transverse */
{
 m=0;
 l=0;
 moy1=0;
 moy2=0;
 var1=0;
 var2=0;

 /*cacul des moyennes des images*/
 for(i=0;i<width;i++)
   for (k=0;k<depth;k++)
   {
    if (im1->mri[i][j][k]!=0)
    {
     l++;
     moy1 += im1->mri[i][j][k];
    }
    if (im2->mri[i][j][k]!=0)
    {
     m++;
     moy2 += im2->mri[i][j][k];
    }
   }		
	
 moy1= rcoeff1*moy1/l;
 moy2= rcoeff2*moy2/m;

 /*calcul des variances des images*/
 for(i=0;i<width;i++)
   for (k=0;k<depth;k++)
   {
   if (im1->mri[i][j][k]!=0)
    var1 +=(rcoeff1*im1->mri[i][j][k]-moy1)*(rcoeff1*im1->mri[i][j][k]-moy1);
   if (im2->mri[i][j][k]!=0)
    var2 +=(rcoeff2*im2->mri[i][j][k]-moy2)*(rcoeff2*im2->mri[i][j][k]-moy2);
   }

 var1= var1/l;
 var2=var2/m;
 sygma1=(float)sqrt(var1);
 sygma2=(float)sqrt(var2);
	
 for(i=0;i<width;i++)
   for (k=0;k<depth;k++)
   if (im1->mri[i][j][k]==0) norm=0;
   else
   { 
    aux1=(sygma2)/sygma1;
    aux2=moy2-((moy1*sygma2)/sygma1);
    norm=aux1*rcoeff1*im1->mri[i][j][k]+aux2;
    if((mri_temp[i*height*depth+j*depth+k]=norm)>max) max=norm;
    min=MINI(min,norm);
   }

    
}

/* Calcul de imres->maxpixel et imres->icomp et imres->rcoeff a partir de max*/
err=imx_brukermax_3d(max,min,imres);
rcoeff=imres->rcoeff;

for (i=0;i<width;i++)
 for (j=0;j<height;j++)
  for(k=0;k<depth;k++)
  {
  if (im1->mri[i][j][k]!=0)
   imres->mri[i][j][k]=(TYPEMRI3D)floor(mri_temp[i*height*depth+j*depth+k]/rcoeff);
  else
   imres->mri[i][j][k]=0;
  }

/* imres->min_pixel=0;*/

free(mri_temp);

return;	
}

/************************ norm_tarik_3d() ***********************************/ 
/*	Normalisation d'une image 3D par rapport �une autrei		    */ 
/*	Normalisation O.V.N.I.			                            */ 
/*                                                                          */
/****************************************************************************/
void norm_tarik_3d(void)
{
	int im_1,im_2,im_res;
	int N,err=0;

	im_1=GET_PLACE3D(TEXT0231);
	im_2=GET_PLACE3D(TEXT0031);
	im_res=GET_PLACE3D(TEXT0006);
	
	while((N= GET_INT(TEXT0187, 0, &err))<=0)
          PUT_WNDMSG(TEXT0189);
        CLEAR_WNDMSG();

	imx_norm_tarik_3d(im_1,im_2,im_res,N);
	
	show_picture_3d(im_res);
}

/****************************************************************************/
void imx_norm_tarik_3d(int im_1, int im_2, int im_res, int N)
{	
	grphic3d *im1,*im2,*imres;
	
	im1=ptr_img_3d(im_1);	
	im2=ptr_img_3d(im_2);
	imres=ptr_img_3d(im_res);
	
        imx_norm_tarik_3d_p(im1,im2,imres,N);
}


/****************************************************************************/
int imx_norm_tarik_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres, int N)
{

float *mri_temp;
int width,height,depth,width2,height2,depth2;
int i,j,k,f,g,h,w,l,m;
int inf1,inf2,inf3,sup1,sup2,sup3;
int err;
int N3;
float moy1,moy2,sygma1,sygma2;
float rcoeff1,rcoeff2,rcoeff;
float aux1,aux2;
float norm;
float max=0,min=LONGMAX;


width=im1->width;
height=im1->height;
depth=im1->depth;
w=N/2;
width2=width-w;
height2=height-w;
depth2=depth-w;
rcoeff1=im1->rcoeff;
rcoeff2=im2->rcoeff;
N3=N*N*N;

imx_copie_param_3d_p(im1,imres);

if (( mri_temp = (float *)malloc(sizeof(double)*width*height*depth))==NULL)
{
  printf("failled to allocate memory\n");
  return (-1);
}

for (i=0;i<width;i++)
 for (j=0;j<height;j++)
  for (k=0;k<depth;k++)
  {
    mri_temp[i*height*depth+j*depth+k]=0.0;
    imres->mri[i][j][k]=0;
  }

for(i=w;i<width2;i++)
 for(j=w;j<height2;j++)
  for (k=w;k<depth2;k++)
  {
    if (im1->mri[i][j][k]==0) norm=0;
    else{
    l=0;
    m=0;
    moy1=0.0;
    moy2=0.0;
    inf1=i-w;
    sup1=i+w;
    inf2=j-w;
    sup2=j+w;
    inf3=k-w;
    sup3=k+w;
    for(f=inf1;f<=sup1;f++)
     for(g=inf2;g<=sup2;g++)
      for(h=inf3;h<=sup3;h++)
      {
	if (im1->mri[f][g][h]!=0)
	{
	  l++;
	  moy1 += im1->mri[f][g][h];  
	}
	if (im2->mri[f][g][h]!=0)
	{
	  m++;
	  moy2 += im2->mri[f][g][h]; 
        }
      }	
    moy1=moy1/l;
    moy2=moy2/m;
    
    sygma1 = (float)fabs(im1->mri[i][j][k]-moy1);
    sygma2 = (float)fabs(im2->mri[i][j][k]-moy2);
    
    aux1=(rcoeff2*sygma2)/(rcoeff1*sygma1);
    aux2=rcoeff2*moy2-rcoeff1*moy1*aux1;
    norm=aux1*rcoeff1*im1->mri[i][j][k]+aux2;
    /*
    printf("moy1=%g moy2=%g sygma1=%g sygma2=%g\n",moy1,moy2,sygma1,sygma2);
    printf("I1=%g I2=%g norm=%g\n",rcoeff1*im1->mri[i][j][k],rcoeff2*im2->mri[i][j][k],norm);
    */
    }
    
    if((mri_temp[i*height2*depth2+j*depth2+k]=norm)>max) max=norm;
    min=MINI(min,norm);
    
  }

/* Calcul de imres->maxpixel et imres->icomp et imres->rcoeff a partir de max*/
err=imx_brukermax_3d(max,min,imres);
rcoeff=imres->rcoeff;

for (i=w;i<width2;i++)
 for (j=w;j<height2;j++)
  for (k=w;k<depth2;k++)
  imres->mri[i][j][k]=(int)floor(mri_temp[i*height2*depth2+j*depth2+k]/rcoeff);

imres->min_pixel=0;

free(mri_temp);

return(1);
}


/*******************************************
** --  imx_norm_rcoeff_3d_p() ---- Musse Juin 1998 ------
*/
/*! 
** normalise la valeur du rcoeff des deux images
** pour qu'elles aient le meme rcoeff au final
** \remark l'image modifiee est celle dont le rcoeff est le
** plus petit (pour des problemes de depassement de dynamique)
**	\param im1 : une image 3d (E/S)
**  \param im2 : une image 3d (E/S)
********************************************/
int imx_norm_rcoeff_3d_p(grphic3d *im1, grphic3d *im2)
{
float  rmod,rref,imod,iref;
grphic3d *immod,*imref;
int    wdth,hght,dpth,i,j,k;


/*choix de l'image a modifier et de l'image de reference*/
/*en fonction de la valeur des rcoeffs des images*/
if (im1->rcoeff<im2->rcoeff)
  {	
    immod=im1;rmod=im1->rcoeff;imod=im1->icomp;
    imref=im2;rref=im2->rcoeff;iref=im2->icomp;
  }
else
  {	
    immod=im2;rmod=im2->rcoeff;imod=im2->icomp;
    imref=im1;rref=im1->rcoeff;iref=im1->icomp;
  }


/*transformation de l'image immod*/
wdth=immod->width;hght=immod->height;dpth=immod->depth;
for (i=0;i<wdth;i++)
 for (j=0;j<hght;j++)
  for (k=0;k<dpth;k++)
   {
    immod->mri[i][j][k]=(int)floor((((float)immod->mri[i][j][k])*rmod/rref));
   }
/*mise a jour des valeurs de icomp et rcoeff*/
immod->rcoeff=rref;
immod->icomp=iref; 
/*mise a jour de valeurs de max et min pixel*/
immod->max_pixel=(long)((float)immod->max_pixel*rmod/rref);
immod->min_pixel=(long)((float)immod->min_pixel*rmod/rref);
immod->cutoff_max= immod->max_pixel;
immod->cutoff_min= immod->min_pixel;
return(1);
}

/*******************************************************************************
**        normalise_mixture_gaussienne                                                
**                                                                                                   
*******************************************************************************/

void normalise_mixture_gaussienne(void)
{
#ifdef HAS_ATLAS_PERFUSION	
	int im_1,im_2,im_res;
	int nb_gaussienne,e;
	double facteur;
	grphic3d *im1,*im2,*imres;
	
	im_1 = GET_PLACE3D(TEXT0231);
	im_2 = GET_PLACE3D(TEXT0031);
	im_res=GET_PLACE3D(TEXT0006);
	
	nb_gaussienne= GET_INT("nb de gaussiennes", 4, &e);
	facteur = GET_DOUBLE("facteur", 0.05, &e);
	
	im1=ptr_img_3d(im_1);	
	im2=ptr_img_3d(im_2);
	imres=ptr_img_3d(im_res);
	
	CPP_normalisation_GM(im1,im2,imres,nb_gaussienne,facteur);
	
	show_picture_3d(im_res);
#endif

}

int imx_query_normalisation_type_noReca_3d()
{
 int norm_type=0, i;
 char *quest[5];

 for (i=0;i<5;i++) quest[i]=CALLOC(80, char);

 strcpy(quest[0],TEXT0229);
 strcpy(quest[1],TEXT0230);
 strcpy(quest[2],TEXT0234);
 strcpy(quest[3],TEXT1234);
 strcpy(quest[4],"\0");

 norm_type=GETV_QCM(TEXT0232,(char **)quest);

 for (i=0;i<5;i++) FREE(quest[i]);

 return norm_type;
}

int imx_query_normalisation_type_reca_3d()
{
 int norm_type=0, i;
 char *quest[4];

 for (i=0;i<4;i++) quest[i]=CALLOC(80, char);

 strcpy(quest[0],TEXT0233);
 strcpy(quest[1],TEXT1235);
 strcpy(quest[2],TEXT1236);
 strcpy(quest[3],"\0");

 norm_type=GETV_QCM(TEXT0232,(char **)quest);

 for (i=0;i<4;i++) FREE(quest[i]);

 return norm_type;
}

norm_func_t imx_choose_normalisation_roi_noReca_fct(int norm_type)
{
 norm_func_t res_fct=NULL;
 switch (norm_type)
 {
	case 0:  // Moyenne 
	       res_fct=imx_norm_roi_mean_3d_p;
	       break;
	case 1:  // Moyenne + ecty 
	       res_fct=imx_norm_roi_meanecty_3d_p;
	       break;
	case 2:  // rapport des moyennes 
	       res_fct=imx_norm_roi_3d_p;
	       break;
	case 3:  // rapport des minmax 
	       res_fct=imx_norm_roi_minmax_3d_p;
	       break;
	default:  // defaut moyenne +ecty 
	       res_fct=imx_norm_roi_meanecty_3d_p;
	       break;
 }

 return res_fct;
}

norm_func_t imx_choose_normalisation_roi_reca_fct(int norm_type)
{
 norm_func_t res_fct=NULL;
 switch (norm_type)
 {
	case 0:  // regression lineaire a 2 parametres
	       res_fct=imx_norm_roi_rl_3d_p;
	       break;
	case 1:  // regression lineaire a 1 parametre
	       res_fct=imx_norm_roi_rl_1_param_3d_p;
	       break;
	case 2:  // regression lineaire total least squares
	       res_fct=imx_regression_TLS_roi_3d_p;
	       break;
	default:  // regression lineaire total least squares
	       res_fct=imx_regression_TLS_roi_3d_p;
	       break;
 }

 return res_fct;
}
