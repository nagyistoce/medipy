/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*!************************************************************************
***	
***	\file:		segm_seuillage_3d.c
***
***	project:	Imagix 2.01
***			
***
***	\brief description:    Fichier source pour segmentation 3D
***	
***	
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	Copyright (c) 1997, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***
*************************************************************************/

#include <config.h>
#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_3d.h" 
#include "noyau/gui/imx_picture3d.h"
#include "noyau/imx_lang.h" 
#include "segmentation/segm_seuillage_3d.h" 
#include "traitement/trai_3d.h" 
#include "math/oper_3d.h" 




/*******************************************
** --  threh_cor_3d() ----------------------------
**
** extraction d'une partie de l'image par valeur
********************************************/
void    threh_cor_3d(void)
{
    int im_1,im_res;
    int /*unused variable i,j,*/err=0;
    float seuil;
    /*unused variable grphic *im1;*/

    im_1=GET_PLACE3D(TEXT0044);
    im_res=GET_PLACE3D("Resultat du seuillage:");
    seuil= GET_FLOAT("Valeur de seuil ", 0, &err);
    if(err!=0) return;


    imx_threh_cor_3d(im_1,im_res,seuil,1.);
    show_picture_3d(im_res);

}

/*******************************************
** --  imx_threh_cor_3d() ----------------------------
**
** Seuillage d'une image de correlation. Garde tout ce qui est
** < bas et > haut et >-bas et <-haut
** bas : seuil bas
** haut: seuil haut
** im_1 : image seuille
** im_res : resultat du seuillage
********************************************/
int    imx_threh_cor_3d(int im_1, int im_res, float bas, float haut)
{
    grphic3d *im1,*imres;
    unsigned int i,j,k;

    im1=ptr_img_3d(im_1);
    imres=ptr_img_3d(im_res);
    bas=bas/im1->rcoeff;
    haut=haut/im1->rcoeff;
    if(PROG_DEBUG) printf("BAS=%f HAUT=%f\n",bas,haut);
    for (i=0;i<im1->width;i++)
      for (j=0;j<im1->height;j++)
        for (k=0;k<im1->depth;k++)
        {
            if (((im1->mri[i][j][k]>=bas)&&(im1->mri[i][j][k]<=haut))
                ||((im1->mri[i][j][k]>=-haut)&&(im1->mri[i][j][k]<=-bas)))
            {
                imres->mri[i][j][k]=im1->mri[i][j][k];
            }
            else
            {
                imres->mri[i][j][k]=0;
            }
        }

    imx_copie_param_3d_p(im1,imres);

    return(1);
}

/*******************************************
** --  threhcluster_cor_3d() ----------------------------
**
** extraction d'une partie de l'image par valeur
** On donne 2 valeur. La plus haute sert de point 
** de depart a un cluster dont la valuer la plus faible
** est le criter d'arret
********************************************/
void    threhcluster_cor_3d(void)
{
    int im_1,im_res;
    int err=0;
    float seuilbas,seuilhaut;
 
    im_1=GET_PLACE3D(TEXT0044);
    im_res=GET_PLACE3D(TEXT0051);
    seuilbas= GET_FLOAT(TEXT0236, 0, &err);
    seuilhaut= GET_FLOAT(TEXT0237, 0, &err);
    if(err!=0) return;


    imx_threhcluster_cor_3d(im_1,im_res,seuilbas,seuilhaut);
    show_picture_3d(im_res);

}

/*******************************************
** --  imx_threhcluster_cor_3d() ----------------------------
**
** Seuillage d'une image de correlation.
** Ne garde que des cluster dont le seed point est superieur
** a haut et le critere d'arret est superieur a bas.
**
** bas : seuil bas
** haut: seuil haut
** im_1 : image seuille
** im_res : resultat du seuillage
********************************************/
int    imx_threhcluster_cor_3d(int im_1, int im_res, float bas, float haut)
{
  grphic3d *im1,*imres;

  im1=ptr_img_3d(im_1);
  imres=ptr_img_3d(im_res);

  imx_threhcluster_cor_3d_p(im1,imres,bas,haut);

  return(1);
}

/*******************************************
** --  imx_threhcluster_cor_3d_p() ----------------------------
**
** Seuillage d'une image de correlation.
** Ne garde que des cluster dont le seed point est superieur
** a haut et le critere d'arret est superieur a bas.
**
** bas : seuil bas
** haut: seuil haut
** im_1 : image seuille
** im_res : resultat du seuillage
********************************************/
int    imx_threhcluster_cor_3d_p(grphic3d *im1,grphic3d *imres, float bas, float haut)
{
  grphic3d *imtemp;
  int i,j,k;
  int wdth,hght,dpth;
  int ntab,typevois;
  long valeur;
  float tab[1];

  imtemp=cr_grphic3d(im1);
  wdth=im1->width; 
  hght=im1->height; 
  dpth=im1->depth; 
  
  bas=bas/im1->rcoeff;
  haut=haut/im1->rcoeff;
  if(PROG_DEBUG) printf("BAS=%f HAUT=%f\n",bas,haut);

  tab[0]=bas;
  ntab=1;
  valeur=1;
  typevois=0; /* (0:croix 1:cube) */
  imx_copie_param_3d_p(im1,imres);
  imtemp->rcoeff=1;
  imtemp->icomp=0;
  imtemp->max_pixel=1;
  imtemp->min_pixel=0;

  for (i=0;i<wdth;i++)
    for (j=0;j<hght;j++)
      for (k=0;k<dpth;k++)
        {
        if (im1->mri[i][j][k] > haut && imtemp->mri[i][j][k] == 0) {
		  growing_3d_p(i,j,k,im1,imtemp,fct_grow_con_pos_3d
						,tab,ntab,valeur,typevois);
	  	  }
		}      
      
  imx_mul_3d_p(im1,imtemp,imres);

  return(1);
}

