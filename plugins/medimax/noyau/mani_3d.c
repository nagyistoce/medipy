/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/**-----------------------------------------------------------------------
**/	
/*!	\file		mani_3d.c 
***
***
***	project:	Imagix 1.01
***			
***
*** \brief	description:    Manipulation d'images 3d
***	
***	
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***
***---------------------------------------------------------------------*/
#include <config.h>
#include <stdio.h> 
#include <stdlib.h> 	/* Needed when you use GET_FLOAT	*/
#include <math.h> 	/* Needed when you use math function	*/
#include <string.h>
 
#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_lang.h" 
#include "noyau/imx_3d.h"
#include "outils/imx_misc.h"

#include "recalage/chps_3d.h"
#include "mani_3d.h"
#include "traitement/trai_3d.h"
#include "math/oper_3d.h"
#include "noyau/gui/imx_picture3d.h"
#include "math/imx_matrix.h"
#include "noyau/io/imx_export_file.h"
#include "noyau/zone_3d.h"
#include "recalage/transformations_3d.h"

/***********************************************
**
**---------  imx_shift_bou_3d()  ---------------------
**
**   Shifter une image avec boucle
**
**   shft_x ,  shift suivant x
**   shft_y , shift suivant y
**   shft_z , sifht suivant z
**   fonction du menu
**
**   G.ENGER juin 1996, ARMSPACH 1997
**
************************************************/
int imx_shift_bou_3d(void)
{
  int im_1,im_res;
  int shft_x,shft_y,shft_z,err=0;

  grphic3d *im1;
  grphic3d *imres;

  im_1=GET_PLACE3D(TEXT0135);
  im_res=GET_PLACE3D(TEXT0006);
  shft_x= GET_INT(TEXT0136, 0, &err);
  shft_y= GET_INT(TEXT0137, 0, &err);
  shft_z= GET_INT(TEXT0165, 0, &err);
  
  im1=ptr_img_3d(im_1);
  imres=ptr_img_3d(im_res);

  imx_shift_bou_3d_p(im1,imres,shft_x,shft_y,shft_z);

  show_picture_3d(im_res);
return(1);
}

/*************************************************
**
** --------  imx_shift_bou_3d_p(im1,imres,shft_x,shft_y,shft_z)  -------------------
**
**   Shifter une image avec boucle
**
**   shft_x ,  shift suivant x
**   shft_y , shift suivant y
**   shft_z ,  shift suivant Z
**   fonction avec pointeurs
**
**   G.ENGER juin 1996 JP ARMSPACH 1997
**
*****************************************************/
int imx_shift_bou_3d_p(grphic3d *im1, grphic3d *imres, int shft_x, int shft_y, int shft_z)
{
  int *v;
  int i,j,k;
  int wdth,hght,dpth;

  wdth=im1->width;            /* im1 et im2 memes dimensions */
  hght=im1->height;
  dpth=im1->depth;

/*  Afin que le sp ne fait pas un core il faut oligatoirement un 
    shift positif.  Remodifier 5/2001                                               */
  if (shft_x < 0) shft_x=wdth+shft_x;
  else shft_x=shft_x;
  if (shft_y < 0) shft_y=hght+shft_y;
  else shft_y=shft_y;
  if (shft_z < 0) shft_z=dpth+shft_z;
  else shft_z=shft_z;

  v=alloc_ivector(wdth);

  for (j=0;j<hght;j++)
  for (k=0;k<dpth;k++)
	{
	for (i=0;i<wdth;i++)
		{
		if((i+shft_x)>=wdth)
			v[i+shft_x-wdth]=im1->mri[i][j][k];
		else
			v[i+shft_x]=im1->mri[i][j][k];
		}

	for (i=0;i<wdth;i++)
		imres->mri[i][j][k]=v[i];
	}

  free_ivector(v,wdth);
  v=alloc_ivector(hght);

  for (i=0;i<wdth;i++)
  for (k=0;k<dpth;k++)
	{
	for (j=0;j<hght;j++)
		{	
		if((j+shft_y)>=hght)
			v[j+shft_y-hght]=imres->mri[i][j][k];
		else
			v[j+shft_y]=imres->mri[i][j][k];
		}

	for (j=0;j<hght;j++)
		imres->mri[i][j][k]=v[j];
	}

  free_ivector(v,hght);
  v=alloc_ivector(dpth);

  for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
	{
	for (k=0;k<dpth;k++)
		{	
		if((k+shft_z)>=dpth)
			v[k+shft_z-dpth]=imres->mri[i][j][k];
		else
			v[k+shft_z]=imres->mri[i][j][k];
		}

	for (k=0;k<dpth;k++)
		imres->mri[i][j][k]=v[k];
	}

  free_ivector(v,dpth);

  imx_copie_param_3d_p(im1,imres);

return(1);
}

/*******************************************
** --  imx_copie_3d()  ---------------------
**
**  Operation de copie des images
**  im_deb : Image de depart
**  im_res : Image resultat
*******************************************/
int    imx_copie_3d(int im_deb, int im_res)
{
   grphic3d *imdeb,*imres;

  imdeb=ptr_img_3d(im_deb);
  imres=ptr_img_3d(im_res);

  imx_copie_3d_p(imdeb,imres);


  return(1);

}

/*********************************************
**  ---  imx_copie_3d_p()
*/
/*!  \brief Operation de copie des images 3d
 
	\param imdeb : pointeur de l'image de depart
	\param imres : pointeur de l'image resultat
	\retval int 1, imres est modifiee
*/
/**********************************************/
int   imx_copie_3d_p(grphic3d *imdeb, grphic3d *imres)
{   
   int i,j,k;
   int w,h,d;
   int err=0;
   TYPEMRI3D ***imdebMRI=NULL, ***imresMRI=NULL;

   //dimensions de l'image de depart
   w = (int)imdeb->width; h = (int)imdeb->height; d = (int)imdeb->depth;

   //on s'assure que le buffer est suffisant et on recopie tous les parametres
   err=resize_grphic3d_buffer(imres, (TDimension)w, (TDimension)h, (TDimension)d);
   if (err) return err;
   imx_copie_param_3d_p(imdeb,imres);

   //recopie du buffer de donnees
   imdebMRI=imdeb->mri; imresMRI=imres->mri;
   for(i=0;i<w;i++)
     for(j=0;j<h;j++)
       for(k=0;k<d;k++)
       {
        imresMRI[i][j][k]=imdebMRI[i][j][k];
       }

  return err;
}

/*******************************************
** --  imx_copie_positif_3d()  ---------------------
**
**  Operation de copie des images en ne gardant que la valeur positive
**  im_deb : Image de depart
**  im_res : Image resultat
*******************************************/
int    imx_copie_positif_3d(int im_deb, int im_res)
{
   grphic3d *imdeb,*imres;

  imdeb=ptr_img_3d(im_deb);
  imres=ptr_img_3d(im_res);

  imx_copie_positif_3d_p(imdeb,imres);

  return(1);

}

/*********************************************
**  ---  imx_copie_positif_3d_p()
**
**  Operation de copie des images 3d en ne gardant que la valeur positive
** 
**  imdeb : pointeur de l'image de depart
**  imres : pointeur de l'image resultat
**********************************************/
int   imx_copie_positif_3d_p(grphic3d *imdeb, grphic3d *imres)
{   
   int i,j,k;
   int wdth,hght,dpth;
   
   wdth=imdeb->width;hght=imdeb->height;dpth=imdeb->depth;
   imx_copie_param_3d_p(imdeb,imres);
    
   for(i=0;i<wdth;i++)
     for(j=0;j<hght;j++)
       for(k=0;k<dpth;k++)
	     {
	     if (imdeb->mri[i][j][k] >= 0) 
		   imres->mri[i][j][k]=imdeb->mri[i][j][k];
		 else
		   imres->mri[i][j][k]=0; 
         }
  imres->min_pixel=0;
  imres->cutoff_min=0;

  return(1);

}

/*******************************************
** --  imx_copie_positif_min 3d()  ---------------------
**
**  Operation de copie des images en  remplacant valeur negative
**  par valeur min d'un voisinage
**  im_deb : Image de depart
**  im_res : Image resultat
**  wnd : fenetre d'observation
*******************************************/
int    imx_copie_positif_min_3d(int im_deb, int im_res, int wnd)
{
   grphic3d *imdeb,*imres;

  imdeb=ptr_img_3d(im_deb);
  imres=ptr_img_3d(im_res);

  imx_copie_positif_min_3d_p(imdeb,imres,wnd);

  return(1);

}

/*********************************************
**  ---  imx_copie_positif_min_3d_p()
**
**  Operation de copie des images en  remplacant valeur negative
**  par valeur min d'un voisinage
**  imdeb : Image de depart
**  imres : Image resultat
**  wnd : fenetre d'observation
**********************************************/
int   imx_copie_positif_min_3d_p(grphic3d *imdeb, grphic3d *imres, int wnd)
{   
 int i,j,k,l,m,n;
 int err=0;
 int width,height,depth;
 int i1, i2, j1, j2, k1, k2;
 TYPEMRI3D ***imdebMRI=NULL, ***imtempMRI=NULL;
 TYPEMRI3D val=0;
 grphic3d *imtemp;
   
 if ((!imdeb)||(!imres))
  { fprintf(stderr,"les deux images doivent etre non nulles dans imx_median_3d_p\n"); return 1; }

 // creation image temporaire
 imtemp = cr_grphic3d(imdeb);

 //on s'assure que le buffer de imres est suffisant et on recopie tous les parametres
 err=resize_grphic3d_buffer(imres, (TDimension)imdeb->width,
                                    (TDimension)imdeb->height, (TDimension)imdeb->depth);
 if (err) return err;
 imx_copie_param_3d_p(imdeb,imres);

 width=(int)imdeb->width; height=(int)imdeb->height; depth=(int)imdeb->depth;

 //on filtre toute l'image et on gere les effets de bords en ne tenant compte que du voisinage dans l'image
  wnd=wnd/2; imdebMRI=imdeb->mri; imtempMRI=imtemp->mri;
  for(i=0;i<width;i++)
    {
    i1=MAXI(i-wnd, 0); i2=MINI(width-1, i+wnd);
    for(j=0;j<height;j++)
      {
      j1=MAXI(j-wnd, 0); j2=MINI(height-1, j+wnd);
      for(k=0;k<depth;k++)
        {
        k1=MAXI(k-wnd, 0); k2=MINI(depth-1, k+wnd);
        if (imdebMRI[i][j][k] >= 0) 
 	   imtemp->mri[i][j][k]=imdebMRI[i][j][k];
        else {
	  val=(TYPEMRI3D)imdeb->max_pixel;
          for(l=i1;l<=i2;l++)
            for(m=j1;m<=j2;m++)
              for(n=k1;n<=k2;n++)
                if (imdebMRI[l][m][n] >= 0) val=MINI(val,imdebMRI[l][m][n]);
 	  if (val == imdeb->max_pixel)
	    imtemp->mri[i][j][k]=imdebMRI[i][j][k];
	  else
	    imtempMRI[i][j][k]=(TYPEMRI3D)val; 
          }
        }
      }
    }
  
 imx_copie_3d_p(imtemp,imres);
 imx_inimaxminpixel_3d_p(imres);
 free_grphic3d(imtemp);
 
 return(0);
}

/*********************************************
**  ---  imx_maxrcoeff_multicoupe_p()
**
**  Operation de copie des images 3d
** 
**  imres : pointeur de l'image resultat   3d
**  nbslice :  nombre de slice de l'image 3d
**  max : tableau des max des slices
**  read_format_3d : format de transfer de l'image 
**********************************************/
int   imx_maxrcoeff_multicoupe_p(grphic3d *imres, int nbslice, double *tabrcoeff3d, float max, int read_format_3d)
{   
  int i,j,slice;
  int wdth,hght,dpth;
  float rcoeff,min=0;
   
/*MISE EN PLACE DES DIMENSIONS DES IMAGES EN FONCTION DES FORMATS DE LECTURE*/
  wdth=imres->width;
  hght=imres->height;        
  dpth=imres->depth;        

  switch (read_format_3d) 
    {
    case 1:  /*on echange k et j*/
           wdth=wdth;
           hght=dpth;
	   dpth=hght;
           break;
    case 2:  /*on echange k et i*/
 	   wdth=wdth;
 	   hght=hght;
	   dpth=dpth;
 	   break;
   case 3:  /*on echange i et j*/
           wdth=hght;
           hght=wdth;
	   dpth=hght;
          break;
   case 4:  /*on echange rien*/
           wdth=wdth;
           hght=hght;
	   dpth=dpth;
          break;
   case 5:  /*on echange rien*/
           wdth=wdth;
           hght=hght;
	   dpth=dpth;
          break;
    case 6:  /*on echange k et j*/
           wdth=wdth;
           hght=dpth;
	   dpth=hght;
           break;
   }

  imx_brukermax_3d(max,min,imres);

  for (slice=0;slice<nbslice;slice++)  {
    rcoeff=(float)(tabrcoeff3d[slice]/imres->rcoeff);
    for(i=0;i<wdth;i++)
      for(j=0;j<hght;j++) {
        switch (read_format_3d) {
          case 1:  /* Coupe transversale */
            imres->mri[i][slice][j]=(TYPEMRI3D)floor(rcoeff*(float)imres->mri[i][slice][j]);
	    break;
          case 2:   /* Format IPB paravision*/
            imres->mri[slice][j][i]=(TYPEMRI3D)floor(rcoeff*(float)imres->mri[slice][j][i]);
	    break;
          case 3:  /*  Format Roufach 1 paravision*/
            imres->mri[j][i][slice]=(TYPEMRI3D)floor(rcoeff*(float)imres->mri[j][i][slice]);
	    break;
          case 4:  /*  Format Roufach 2 paravision*/
            imres->mri[i][j][slice]=(TYPEMRI3D)floor(rcoeff*(float)imres->mri[i][j][slice]);
	    break;
          case 5:  /*  Format IPB perso*/
            imres->mri[i][j][slice]=(TYPEMRI3D)floor(rcoeff*(float)imres->mri[i][j][slice]);
	    break;
          case 6:  /* Coupe transversale */
            imres->mri[i][nbslice-slice-1][j]=(TYPEMRI3D)floor(rcoeff*(float)imres->mri[i][nbslice-slice-1][j]);
	    break;
	  
	  }
        }
      }
 
  return(1);
}

/*******************************************
** --  copie_byte_3d() ----------------------------
**
** Op�ations de copie d'une image forcee a 1
********************************************/
void	copie_byte_3d(void)
{ 
  int im_1,im_res;

 /* Question concerning type of operation  */
 /* and image					*/ 
  
 /*     Image   				*/
  im_1=GET_PLACE3D(TEXT0028);
  im_res=GET_PLACE3D(TEXT0023);
		
  imx_copie_byte_3d(im_1,im_res);	
  show_picture_3d(im_res);

  return;

}
/*******************************************
** --  copie_byte_maxmin_3d() ----------------------------
**
** Op�ations de copie d'une image forcee a 1
********************************************/
void	copie_byte_maxmin_3d(void)
{ 
    int im_1,im_res;
    
    /* Question concerning type of operation  */
    /* and image					*/ 
    
    /*     Image   				*/
    im_1=GET_PLACE3D(TEXT0028);
    im_res=GET_PLACE3D(TEXT0023);
    
    imx_copie_byte_maxmin_3d(im_1,im_res);	
    show_picture_3d(im_res);
    
    return;
}


/*******************************************
** --  imx_histogram_equalization_3d_p() ----------------------------
**
**   histogramme egalisation
**
**  im1:	image de depart
**  imres:	image resultat
**  roi_: 	masque de traitement (mettre a 1 pour totalite image)
**  Nbpas :	nombre de pas de l'histogramme
**  y :		histogramme
********************************************/
int	imx_histogram_equalization_3d_p(grphic3d *im1, grphic3d *imres, grphic3d *roi, float *hist, long int Nbpas)
{ 
  int 	i,j,k,wdth,hght,dpth;
  long *histeq;
  double total,sum,max;

  wdth=im1->width;
  hght=im1->height;
  dpth=im1->depth;

/*   Creation du tableu histeq   */
  histeq=CALLOC(Nbpas,long);

  total=0.0;max=im1->max_pixel;
  for (i=0;i<Nbpas;i++) total=total+(double)hist[i];
  for (i=0;i<Nbpas;i++)
    {
    sum=0.0;
    for (j=0;j<=i;j++)
       sum=sum+(double)hist[j];
       histeq[i]= (long) (max*sum/total + (double) .5);
    }

  for (i=0;i<wdth;i++)
    for (j=0;j<hght;j++)
      for (k=0;k<dpth;k++)
        {
        if (roi->mri[i][j][k]!=0)
          imres->mri[i][j][k]= (TYPEMRI3D)floor(histeq[im1->mri[i][j][k]*255/im1->max_pixel]);
        else
	      imres->mri[i][j][k]=im1->mri[i][j][k];
        }

  FREE(histeq);

  return(1);
}


/********************************************
** imx_copie_f1_3d() 
*/
/*!
** \brief Op�ations de copie d'une image forcee a 1
** \param im_deb : image de depart
** \param im_res : image d'arrive (E/S)
** \retval 1
********************************************/
int	imx_copie_f1_3d(int im_deb, int im_res)
{ 
  grphic3d *imdeb,*imres;
  int i,j,k;
  int hght,wdth,dpth;

  imdeb=ptr_img_3d(im_deb);
  imres=ptr_img_3d(im_res);
  hght=imdeb->height;
  wdth=imdeb->width;
  dpth=imdeb->depth;

  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
       for(k=0;k<dpth;k++)
          {
           if (imdeb->mri[i][j][k]==0)
             imres->mri[i][j][k]=0;
           else
	     imres->mri[i][j][k]=1;
           }

  imx_copie_param_3d(im_deb,im_res);
  imres->max_pixel=1;
  imres->min_pixel=0;
  imres->cutoff_max=1;
  imres->cutoff_min=0;
  imres->icomp=0;
  imres->rcoeff=1;

  return(1);
}

/*******************************************
** --  imx_copie_f1_3d_p() ----------------------------
**
** Op�ations de copie d'une image forcee a 1
**  im_deb : image de depart
**  im_res : image d'arrive
********************************************/
int	imx_copie_f1_3d_p(grphic3d *imdeb,grphic3d *imres)
{ 
  int i,j,k;
  int hght,wdth,dpth;

  hght=imdeb->height;
  wdth=imdeb->width;
  dpth=imdeb->depth;

  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
       for(k=0;k<dpth;k++)
          {
           if (imdeb->mri[i][j][k]==0)
             imres->mri[i][j][k]=0;
           else
	     imres->mri[i][j][k]=1;
           }

  imx_copie_param_3d_p(imdeb,imres);
  imres->max_pixel=1;
  imres->min_pixel=0;
  imres->cutoff_max=1;
  imres->cutoff_min=0;
  imres->icomp=0;
  imres->rcoeff=1;

  return(1);
}

/*******************************************
** --  imx_copie_byte_3d() ----------------------------
**
** Op�ations de copie d'une image forcee a 1
**  im_deb : image de depart
**  im_res : image d'arrive
********************************************/
int	imx_copie_byte_3d(int im_deb, int im_res)
{ 
  grphic3d *imdeb,*imres;
  int i,j,k;
  int hght,wdth,dpth;
  float val,rcoeff;

  imdeb=ptr_img_3d(im_deb);
  imres=ptr_img_3d(im_res);
  hght=imdeb->height;
  wdth=imdeb->width;
  dpth=imdeb->depth;
  rcoeff=imdeb->rcoeff; 

  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++) 
      for(k=0;k<dpth;k++) {
      val=rcoeff*imdeb->mri[i][j][k];
      imres->mri[i][j][k] = (TYPEMRI3D)floor( val );
      }

  imx_copie_param_3d(im_deb,im_res);
  imx_inimaxminpixel_3d(im_res);
  imres->icomp=0;
  imres->rcoeff=1;

  return(1);
}

/*******************************************
** --  imx_copie_byte_maxmin_3d() ----------------------------
**
** Op�ations de copie d'une image forcee a 1
**  im_deb : image de depart
**  im_res : image d'arrive
********************************************/
int	imx_copie_byte_maxmin_3d(int im_deb, int im_res)
{ 
  grphic3d *imdeb,*imres;
  int i,j,k;
  int hght,wdth,dpth;
  float val,rcoeff;

  imdeb=ptr_img_3d(im_deb);
  imres=ptr_img_3d(im_res);
  hght=imdeb->height;
  wdth=imdeb->width;
  dpth=imdeb->depth;
  rcoeff=imdeb->rcoeff; 

  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)  
      for(k=0;k<dpth;k++) { 
      val=rcoeff*imdeb->mri[i][j][k];
      imres->mri[i][j][k] =  (TYPEMRI3D)floor( val );
      }

  imx_copie_param_3d(im_deb,im_res);
  imres->max_pixel=100;
  imres->min_pixel=-100;
  imres->cutoff_max=100;
  imres->cutoff_min=-100;
  imres->icomp=0;
  imres->rcoeff=1;

  return(1);
}


/****************************************************************************/
/*  									    */
/* -----cut_head_3d()------------------------------------------------------*/
/*! crop une image sur une hauteur definie par l'utilisateur          */
/****************************************************************************/

void cut_head_3d(void)
{
   int im_deb, im_res;
   double hd_hght; 
   int e=0;

   im_deb=GET_PLACE3D("Image Originale ?" );
   im_res=GET_PLACE3D("Image Resultat ?");
   hd_hght= GET_FLOAT("Length in cm ?", 0, &e);
   
   imx_cut_head_3d(im_deb,im_res,hd_hght);

  }


/****************************************************************************/
/*  									    */
/* -----imx_cut_head_3d()------------------------------------------------------*/
/*! crop une image sur une hauteur definie par l'utilisateur 
	\param im_deb : numero de l'image source         
	\parma im_res : numero de l'image resultat
	\param hd_hght : hauteur en cm a concerver								*/	
/****************************************************************************/

void imx_cut_head_3d(int im_deb, int im_res, double hd_hght)
{
   grphic3d *imdeb,*imres;


   imdeb=ptr_img_3d(im_deb);
   imres=ptr_img_3d(im_res);
   imx_cut_head_3d_p(imdeb,imres,hd_hght);
   
   show_picture_3d(im_res);


  }

/****************************************************************************/
/*  									    */
/* -----imx_cut_head_3d_p()------------------------------------------------------*/
/*! crop une image sur une hauteur definie par l'utilisateur 
	\param im_deb :image source         
	\parma im_res : image resultat (E/S)
	\param hd_hght : hauteur en cm a concerver								*/	
/****************************************************************************/

void imx_cut_head_3d_p(grphic3d *imdeb, grphic3d *imres, double hd_hght)
{
   int i,j,k;
   int width,height,depth;
   int nb_slices;
   int starty,endy;
   int flag=0;


   width=imdeb->width;
   height=imdeb->height;
   depth=imdeb->depth;


   nb_slices=(int)(0.5+(10.0*hd_hght/imdeb->dy));

   printf("hd_hght=%f dy=%f slices kept=%d\n",hd_hght,imdeb->dy,nb_slices);

   j=0;

   while(flag==0 && j<height)
        {
	 for(i=0;i<width;i++)
	    for(k=0;k<depth;k++)
	       if(imdeb->mri[i][j][k]!=0)
		 flag=1;
         j++;
	};  

    starty=j-1;
    if(starty+nb_slices<height)
      endy=j+nb_slices;
      
	 else
	    endy=height;
     
  printf("keeping %f cm between slices %d - %d \n",hd_hght,starty,endy); 


   for(j=0;j<starty;j++)
      for(i=0;i<width;i++)
         for(k=0;k<depth;k++)
            imres->mri[i][j][k]=0;

   for(j=starty;j<endy;j++)
      for(i=0;i<width;i++)
         for(k=0;k<depth;k++)
            imres->mri[i][j][k]=imdeb->mri[i][j][k];

   for(j=endy;j<height;j++)
      for(i=0;i<width;i++)
         for(k=0;k<depth;k++)
            imres->mri[i][j][k]=0;
 
  imx_copie_param_3d_p(imdeb,imres);
  imx_inimaxpixel_3d_p(imres);
  }





/*******************************************
** --  modify_pixsize_3d() ---- JPA ----
**
**  Permet de modifier la taille du pixel 
**  de l' image selectionnne en 3D
**
********************************************/
void  modify_pixsize_3d(void)
{
  int im_choix,err=0;
  float  valeur;
  grphic3d *imchoix;
  char str[80];
  
  im_choix=GET_WND3D(NULL);
  imchoix=ptr_img_3d(im_choix);
  sprintf(str,TEXT0108" ( %f )\n",imchoix->dx);
  valeur= GET_FLOAT(str, 0, &err);
  if (!err)
   {
   imchoix->dx=valeur;
   sprintf(str,TEXT0109" ( %f )\n",imchoix->dy);
   valeur= GET_FLOAT(str, 0, &err);
   if (!err)
    { 
    imchoix->dy=valeur;
    sprintf(str,TEXT0110" ( %f )\n",imchoix->dz);
    valeur= GET_FLOAT(str, 0, &err);
    if (!err)
     imchoix->dz=valeur;
    }
   }
  show_picture_3d(im_choix);
}

void imx_modify_norrcoeff_3d_p(grphic3d* imchoix,int valeur)
{
  int i,j,k,wdth,hght,dpth;
  double coeff,result;


  wdth=imchoix->width;
  hght=imchoix->height;
  dpth=imchoix->depth;

   coeff=(double) imchoix->rcoeff/pow((double) 2,(double) valeur);
   for (i=0;i<wdth;i++)
     for (j=0;j<hght;j++)
       for (k=0;k<dpth;k++)
         {
         result=(double) coeff*imchoix->mri[i][j][k];
         imchoix->mri[i][j][k]=(TYPEMRI3D)floor(result);
         }

   imchoix->icomp=(float) valeur;
   imchoix->rcoeff= (float)pow((double) 2,(double) imchoix->icomp );

}



void imx_modify_norrcoeff_3d(int pos, int coeff)
{
	grphic3d* im;
	im = ptr_img_3d(pos);

	imx_modify_norrcoeff_3d_p(im,coeff);

    imx_inimaxminpixel_3d(pos);
    show_picture_3d(pos);
}

/*******************************************
** --  modify_norrcoef_3d() ---- Armspach July 1999 ------
*/
/*!  Permet de modifier la valeur de rcoeff 
**  pour l' image selectionnner et ainsi normaliser
**  plusieurs image avec le meme rcoeff
**
********************************************/

void modify_norrcoeff_3d(void)
{
  	int err=0,valeur;
	int pos;
    char str[80];
	
	pos = GET_WND3D(NULL);
  	sprintf(str,TEXT0107);
  	valeur= GET_INT(str, 0, &err);

  if (!err)
	imx_modify_norrcoeff_3d(pos, valeur);
}

/*******************************************
** --  rot_droite_3d() ----------------------------
**
** Rotation �droite de l'image 3D
**  im_deb : image de depart
**  im_res : image d'arrive
********************************************/
int    rot_droite_3d(int im_deb, int im_res)
{   
  grphic3d *imdeb,*imres,*imtmp;
  int i,j,k;
  int hght,wdth,dpth;

  imdeb=ptr_img_3d(im_deb);
  // image temporaire pour l'image resultat
  imtmp=ptr_img_3d((int)NULL);

  hght=imdeb->height;
  wdth=imdeb->width;
  dpth=imdeb->depth;
    
  // formule toute simple qui marche meme si c'est pas un cube
  for (i=0 ; i<wdth ; i++)
    for (j=0 ; j<hght ; j++)
      for (k=0 ; k<dpth ; k++){
	imtmp->mri[hght-j-1][i][k]=imdeb->mri[i][j][k];
      }
  // recopie des parametres de l'img de depart dans l'img tmp
  imx_copie_param_3d_p(imdeb,imtmp);
  imres=ptr_img_3d(im_res);
  // mise a jour des dimensions de l'image temporaire
  // pas la peine de le faire si l'image de depart est un cube
  if ( (wdth != hght)||(wdth != dpth) ) {
     imtmp->height=wdth; 
     imtmp->width=hght;
     imtmp->depth=dpth;
    }
  imtmp->dy=imdeb->dx; 
  imtmp->dx=imdeb->dy;

  imx_copie_param_3d_p(imtmp,imres);
  imx_copie_3d_p(imtmp,imres);

  show_picture_3d(im_res);
  return(1);
}

/*******************************************
** --  rot_gauche_3d() ----------------------------
**
** Rotation �gauche de l'image 3D
**  im_deb : image de depart
**  im_res : image d'arrive
********************************************/
int    rot_gauche_3d(int im_deb, int im_res)
{
  grphic3d *imdeb,*imres,*imtmp;
  int i,j,k;
  int hght,wdth,dpth;

  imdeb=ptr_img_3d(im_deb);
  // image temporaire pour l'image resultat
  imtmp=ptr_img_3d((int)NULL);

  hght=imdeb->height;
  wdth=imdeb->width;
  dpth=imdeb->depth;
    
  // formule toute simple qui marche meme si c'est pas un cube
  for (i=0 ; i<wdth ; i++)
    for (j=0 ; j<hght ; j++)
      for (k=0 ; k<dpth ; k++){
	imtmp->mri[j][wdth-i-1][k]=imdeb->mri[i][j][k];//marat 03/07/2002
      }
  // recopie des parametres de l'img temp dans l'img res
  imx_copie_param_3d_p(imdeb,imtmp);
  imres=ptr_img_3d(im_res);

  // mise a jour des dimensions de l'image temporaire
  // pas la peine de le faire si l'image de depart est un cube
   if ( (wdth != hght)||(wdth != dpth) ) {
    imtmp->height=wdth;
    imtmp->width=hght;
    imtmp->depth=dpth;
    }
  imtmp->dy=imdeb->dx; 
  imtmp->dx=imdeb->dy;

  imx_copie_param_3d_p(imtmp,imres);
  imx_copie_3d_p(imtmp,imres);

  show_picture_3d(im_res);
  return(1);
}

/*******************************************
** --  rot_bas_3d() ----------------------------
**
** Rotation vers le bas de l'image 3D
**  im_deb : image de depart
**  im_res : image d'arrive
********************************************/
int    rot_bas_3d(int im_deb, int im_res)
{
  grphic3d *imdeb,*imres,*imtmp;
  int i,j,k;
  int hght,wdth,dpth;

  imdeb=ptr_img_3d(im_deb);
  // image temporaire pour l'image resultat
  imtmp=ptr_img_3d((int)NULL);

  hght=imdeb->height;
  wdth=imdeb->width;
  dpth=imdeb->depth;
    
  // formule toute simple qui marche meme si c'est pas un cube
  for (i=0 ; i<wdth ; i++)
    for (j=0 ; j<hght ; j++)
      for (k=0 ; k<dpth ; k++){
//	imtmp->mri[i][dpth-k-1][hght-j-1]=imdeb->mri[i][j][k];
	imtmp->mri[i][dpth-1-k][j]=imdeb->mri[i][j][k];
      }
  // recopie des parametres de l'img temp dans l'img res
  imx_copie_param_3d_p(imdeb,imtmp);
  imres=ptr_img_3d(im_res);

  // mise a jour des dimensions de l'image temporaire
  // pas la peine de le faire si l'image de depart est un cube
   if ( (wdth != hght)||(wdth != dpth) ) {
    imtmp->height=dpth;
    imtmp->width=wdth;
    imtmp->depth=hght;
    }
  imtmp->dy=imdeb->dz; 
  imtmp->dz=imdeb->dy;

  imx_copie_param_3d_p(imtmp,imres);
  imx_copie_3d_p(imtmp,imres);

  show_picture_3d(im_res);
  return(1);
}

/*******************************************
** --  rot_haut_3d() ----------------------------
**
** Rotation vers le haut de l'image 3D
**  im_deb : image de depart
**  im_res : image d'arrive
********************************************/
int    rot_haut_3d(int im_deb, int im_res)
{
  grphic3d *imdeb,*imres,*imtmp;
  int i,j,k;
  int hght,wdth,dpth;

  imdeb=ptr_img_3d(im_deb);
  // image temporaire pour l'image resultat
  imtmp=ptr_img_3d((int)NULL);

  hght=imdeb->height;
  wdth=imdeb->width;
  dpth=imdeb->depth;
    
  // formule toute simple qui marche meme si c'est pas un cube
  for (i=0 ; i<wdth ; i++)
    for (j=0 ; j<hght ; j++)
      for (k=0 ; k<dpth ; k++){
//	imtmp->mri[i][k][j]=imdeb->mri[i][j][k];
	imtmp->mri[i][k][hght-1-j]=imdeb->mri[i][j][k];
      }
  // recopie des parametres de l'img temp dans l'img res
  imx_copie_param_3d_p(imdeb,imtmp);
  imres=ptr_img_3d(im_res);

  // mise a jour des dimensions de l'image temporaire
  // pas la peine de le faire si l'image de depart est un cube
   if ( (wdth != hght)||(wdth != dpth) ) {
    imtmp->height=dpth;
    imtmp->width=wdth;
    imtmp->depth=hght;
    }
  imtmp->dy=imdeb->dz; 
  imtmp->dz=imdeb->dy;

  imx_copie_param_3d_p(imtmp,imres);
  imx_copie_3d_p(imtmp,imres);

  show_picture_3d(im_res);
  return(1);

}


/*******************************************
** --  rot_pro_dr_3d() ----------------------------
**
** Rotation vers la gauche suivant un axe vertical
** de l'image 3D
**  im_deb : image de depart
**  im_res : image d'arrive
********************************************/
int    rot_pro_dr_3d(int im_deb, int im_res)
{
  grphic3d *imdeb,*imres,*imtmp;
  int i,j,k;
  int hght,wdth,dpth;

  imdeb=ptr_img_3d(im_deb);
  // image temporaire pour l'image resultat
  imtmp=ptr_img_3d((int)NULL);

  hght=imdeb->height;
  wdth=imdeb->width;
  dpth=imdeb->depth;
    
  // formule toute simple qui marche meme si c'est pas un cube
  for (i=0 ; i<wdth ; i++)
    for (j=0 ; j<hght ; j++)
      for (k=0 ; k<dpth ; k++){
//	imtmp->mri[dpth-k-1][j][wdth-i-1]=imdeb->mri[i][j][k];
	imtmp->mri[dpth-k-1][j][i]=imdeb->mri[i][j][k];
      }
  // recopie des parametres de l'img temp dans l'img res
  imx_copie_param_3d_p(imdeb,imtmp);
  imres=ptr_img_3d(im_res);

  // mise a jour des dimensions de l'image temporaire
  // pas la peine de le faire si l'image de depart est un cube
   if ( (wdth != hght)||(wdth != dpth) ) {
    imtmp->height=hght;
    imtmp->width=dpth;
    imtmp->depth=wdth;
    }
  imtmp->dx=imdeb->dz; 
  imtmp->dz=imdeb->dx;

  imx_copie_param_3d_p(imtmp,imres);
  imx_copie_3d_p(imtmp,imres); 

  show_picture_3d(im_res);
  return(1);
}

/*******************************************
** --  rot_pro_gau_3d() ----------------------------
**
** Rotation vers la gauche suivant un axe vertical
** de l'image 3D
**  im_deb : image de depart
**  im_res : image d'arrive
********************************************/
int    rot_pro_gau_3d(int im_deb, int im_res)
{
  grphic3d *imdeb,*imres,*imtmp;
  int i,j,k;
  int hght,wdth,dpth;

  imdeb=ptr_img_3d(im_deb);
  // image temporaire pour l'image resultat
  imtmp=ptr_img_3d((int)NULL);

  hght=imdeb->height;
  wdth=imdeb->width;
  dpth=imdeb->depth;
    
  // formule toute simple qui marche meme si c'est pas un cube
  for (i=0 ; i<wdth ; i++)
    for (j=0 ; j<hght ; j++)
      for (k=0 ; k<dpth ; k++){
//	imtmp->mri[k][j][i]=imdeb->mri[i][j][k];
	imtmp->mri[k][j][wdth-i-1]=imdeb->mri[i][j][k];
      }
  // recopie des parametres de l'img temp dans l'img res
  imx_copie_param_3d_p(imdeb,imtmp);
  imres=ptr_img_3d(im_res);

  // mise a jour des dimensions de l'image temporaire
  // pas la peine de le faire si l'image de depart est un cube
   if ( (wdth != hght)||(wdth != dpth) ) {
    imtmp->height=hght;
    imtmp->width=dpth;
    imtmp->depth=wdth;
    }
  imtmp->dx=imdeb->dz; 
  imtmp->dz=imdeb->dx;

  imx_copie_param_3d_p(imtmp,imres);
  imx_copie_3d_p(imtmp,imres);

  show_picture_3d(im_res);
  return(1);
}


void repositionner_image_3d()
{
  int im_deb,im_res;
  int vue_sagittale, pos_nez, pos_sommet_crane;
  char *query_vue_sagittale[5];
	char *query_nez[6];
	char *query_sommet_crane1[4];
	char *query_sommet_crane2[4];

	query_vue_sagittale[0]="vue 1";
	query_vue_sagittale[1]="vue 2";
	query_vue_sagittale[2]="vue 3";
	query_vue_sagittale[3]="\0";
	query_vue_sagittale[4]=NULL;
    
	query_nez[0]="gauche de l'image";
	query_nez[1]="haut de l'image";
	query_nez[2]="droite de l'image";
	query_nez[3]="bas de l'image";
	query_nez[4]="\0";
	query_nez[5]=NULL;

	query_sommet_crane1[0]="gauche de l'image";
	query_sommet_crane1[1]="droite de l'image";
	query_sommet_crane1[2]="\0";
	query_sommet_crane1[3]=NULL;

	query_sommet_crane2[0]="haut de l'image";
	query_sommet_crane2[1]="bas de l'image";
	query_sommet_crane2[2]="\0";
	query_sommet_crane2[3]=NULL;
  
  im_deb=GET_PLACE3D(TEXT0001);
  im_res=GET_PLACE3D(TEXT0006);

  vue_sagittale=GETV_QCM("vue sagittale",(char **)query_vue_sagittale);
  vue_sagittale++;
  
  pos_nez=GETV_QCM("position du nez",(char **)query_nez);
  if ((pos_nez==0)||(pos_nez==2))
  {
    pos_sommet_crane=GETV_QCM("position du sommet du crane",(char **)query_sommet_crane2);
    pos_sommet_crane=2*pos_sommet_crane+1;
  }
  else
  {
    pos_sommet_crane=GETV_QCM("position du sommet du crane",(char **)query_sommet_crane1);
    pos_sommet_crane=2*pos_sommet_crane;
  }
  
  imx_repositionner_image_3d(im_deb, im_res, vue_sagittale, pos_nez, pos_sommet_crane);

}

int imx_repositionner_image_3d(int im_deb, int im_res, int vue_sagittale, int pos_nez, int pos_sommet_crane)
{
  int i;
  
  //on ramene la vue sagittale dans la bonne vue
  if (vue_sagittale==1)
  {
   rot_pro_dr_3d(im_deb, im_res);
   im_deb=im_res;
  }
  else if (vue_sagittale==3)
  {
    rot_gauche_3d(im_deb, im_res);
    pos_nez--; if (pos_nez<0) pos_nez=3;
    pos_sommet_crane--; if (pos_sommet_crane<0) pos_sommet_crane=3;
    im_deb=im_res;
  }
  
  //on se ramene dans la bonne configuration de l'image sagittale
  if ((pos_sommet_crane+1==pos_nez)||((pos_sommet_crane==3)&&(pos_nez==0)))
  {
    rot_gauche_3d(im_deb, im_res);
    rot_gauche_3d(im_res, im_res);
    if (pos_nez==1) pos_nez=3;
    else if (pos_nez==3) pos_nez=1;
    if (pos_sommet_crane==1) pos_sommet_crane=3;
    else if (pos_sommet_crane==3) pos_sommet_crane=1;
    im_deb=im_res;
  }
  
  //on tourne selon l'axe x autant de fois que necessaire
  for (i=0;i<pos_nez;i++)
  {
    rot_bas_3d(im_deb, im_res);
    im_deb=im_res;
  }

  if (im_deb!=im_res) imx_copie_3d(im_deb, im_res);

  show_picture_3d(im_res);

  return 0;
}


/****************************************************************************/
/*  									    */
/* -----set_size_3d()------------------------------------------------*/
/*  									    */
/****************************************************************************/
void
set_size_3d(void)
{
	int im_1;
	int nw,nh,nd;
	im_1=GET_PLACE3D("Image 1 ?" );
	nw= GET_INT("new size w", 0, NULL);
	nh= GET_INT("new size h", 0, NULL);
	nd= GET_INT("new size d", 0, NULL);

	imx_set_size_3d(im_1,nw,nh,nd);
   
	show_picture_3d(im_1);

}

/****************************************************************************/
/*  									    */
/* -----imx_set_size_3d()------------------------------------------*/
/*  									    */
/****************************************************************************/
void imx_set_size_3d(int im_1,int nw,int nh,int nd)
{
	grphic3d *im1;

	im1=ptr_img_3d(im_1);
	imx_set_size_3d_p(im1,nw,nh,nd);
}

/****************************************************************************/
/*  									    */
/* -----imx_set_size_3d_p()------------------------------------------*/
/*  									    */
/****************************************************************************/
void imx_set_size_3d_p(grphic3d *im1,int nw,int nh,int nd)
{
	int i,j,k;
	int w1,h1,d1;

	w1=im1->width;
	h1=im1->height;
	d1=im1->depth;


	im1->width=nw;
	im1->height=nh;
	im1->depth=nd;

	for(i=w1;i<nw;i++)
		for(j=0;j<nh;j++)
			for(k=0;k<nd;k++)
				im1->mri[i][j][k]=0;

	for(j=h1;j<nh;j++)
		for(i=0;i<nw;i++)
			for(k=0;k<nd;k++)
				im1->mri[i][j][k]=0;

  	for(k=d1;k<nd;k++)
		for(i=0;i<nw;i++)
			for(j=0;j<nh;j++)
				im1->mri[i][j][k]=0;

	imx_inimaxminpixel_3d_p(im1);

}


/****************************************************************************/
/*  									    */
/* -----set_size_cg_3d()------------------------------------------------*/
/*  									    */
/****************************************************************************/
void
set_size_cg_3d(void)
{
	int im_1;
	int nw,nh,nd;
	
	im_1=GET_PLACE3D("Image 1 ?" );
	nw= GET_INT("new size w", 0, NULL);
	nh= GET_INT("new size h", 0, NULL);
	nd= GET_INT("new size d", 0, NULL);

	imx_set_size_cg_3d(im_1,nw,nh,nd);
   	show_picture_3d(im_1);

}

/****************************************************************************/
/*  									    */
/* -----imx_set_size_cg_3d_p()------------------------------------------*/
/*  									    */
/****************************************************************************/
void imx_set_size_cg_3d(int im_1,int nw,int nh,int nd)
{
	grphic3d *im1;

	im1=ptr_img_3d(im_1);

	imx_set_size_cg_3d_p(im1,nw,nh,nd);

}

/*******************************************************************************/
/*  									    	*/
/* -----imx_set_size_cg_3d_p()------------------------------------------	*/
/*  									   	*/
/*   Toujours verifier avant d'utiliser cette fonction que nw, nh, nd sont soient*/
/*   plus grand ou plus petit que im1->width,im1->hght, im1->dpth 		*/	
/*   Ne gere pas un melange. 							*/
/******************************************************************************/
void imx_set_size_cg_3d_p(grphic3d *im1,unsigned int nw,unsigned int nh,unsigned int nd)
{
    unsigned int i,j,k;
    unsigned int ii,jj,kk;
    unsigned int w1,h1,d1;
    unsigned int dw1,fw1,dh1,fh1,dd1,fd1;
    int flag=0;
    grphic3d *imtemp,*masktemp=NULL;
    grphic3d *maskim1=NULL;

    /* Controle */
    w1=im1->width;   
    h1=im1->height;   
    d1=im1->depth;   
    flag= (nw-w1>=0) ? ((nh-h1>=0)?(((nd-d1)>=0?1:0)):0):((nh-h1<0)?(((nd-d1)<0?1:0)):0);
    flag= (nw>im1->mri_width)?0:(nh>im1->mri_height)?0:(nd>im1->mri_depth)?0:flag;
    if (flag==0) { 
      PUT_ERROR("Nouvelle Taille incoherente");
      return;
      }  

    if( ptr_has_mask_3d_p(im1) ) maskim1=ptr_mask_3d_p(im1);
    imtemp=cr_grphic3d(NULL);
    if( maskim1 != NULL ) masktemp=ptr_mask_3d_p(imtemp);
    w1=im1->width;   
    h1=im1->height;   
    d1=im1->depth;   
    
    if(nw>=w1||nh>=h1||nd>=d1) { /* agrandissement */
      dw1=(nw-w1)/2; fw1=nw-dw1;
      dh1=(nh-h1)/2; fh1=nh-dh1;
      dd1=(nd-d1)/2; fd1=nd-dd1;
      for(i=dw1,ii=0;i<fw1;i++,ii++)
        for(j=dh1,jj=0;j<fh1;j++,jj++)
           for(k=dd1,kk=0;k<fd1;k++,kk++) {
     	      imtemp->mri[i][j][k]=im1->mri[ii][jj][kk]; 
              if( maskim1 != NULL ) masktemp->mri[i][j][k]=maskim1->mri[ii][jj][kk];
	      }
      }	

    else if(nw<=w1||nh<=h1||nd<=d1) { /* reduction */
      dw1=(w1-nw)/2; fw1=w1-dw1;
      dh1=(h1-nh)/2; fh1=h1-dh1;
      dd1=(d1-nd)/2; fd1=d1-dd1;
      for(i=dw1,ii=0;i<fw1;i++,ii++)
	for(j=dh1,jj=0;j<fh1;j++,jj++)
	  for(k=dd1,kk=0;k<fd1;k++,kk++) { 
	     imtemp->mri[ii][jj][kk]=im1->mri[i][j][k]; 
             if( maskim1 != NULL ) masktemp->mri[ii][jj][kk]=maskim1->mri[i][j][k];
             }
      }	
   
   else { printf("PROBLEME dans set_size_cg_3d\n"); }

   for(i=0;i<nw;i++)
     for(j=0;j<nh;j++)
       for(k=0;k<nd;k++)
	 {
	 im1->mri[i][j][k]=imtemp->mri[i][j][k];
         if( maskim1 != NULL ) maskim1->mri[i][j][k]=masktemp->mri[i][j][k];
         }
	  
    im1->width=nw;
    im1->height=nh;
    im1->depth=nd;
    im1->x3dc=im1->width/2;
    im1->y3dc=im1->height/2;
    im1->z3dc=im1->depth/2;
    imx_inimaxminpixel_3d_p(im1);
    if (maskim1)
    {
      maskim1->width=nw;
      maskim1->height=nh;
      maskim1->depth=nd;
      maskim1->x3dc=im1->width/2;
      maskim1->y3dc=im1->height/2;
      maskim1->z3dc=im1->depth/2;
      imx_inimaxminpixel_3d_p(maskim1);
    }
	
    free_grphic3d(imtemp);
}

/****************************************************************************/
/*  									    */
/* -----set_size_cg_trf_3d()------------------------------------------------*/
/*  									    */
/****************************************************************************/
void
set_size_cg_trf_3d(void)
{
	int im_1;
	int nw,nh,nd;
	char *nomfichres;
	char str[1024];
	int err = 0;
  
	im_1=GET_PLACE3D("Image 1 ?" );
	nw= GET_INT("new size w", 0, NULL);
	nh= GET_INT("new size h", 0, NULL);
	nd= GET_INT("new size d", 0, NULL);
	/*question sur l'enregistrement du champ resultat*/
	nomfichres=CALLOC(1024,char);
	sprintf(str,"Save transforamtion ? in %s",_CurrentPath);
	strcpy(nomfichres,SAVE_FILE(str,NULL,&err));
	if ( err )
	  return; 
	put_file_extension(nomfichres,".trf",nomfichres);

	imx_set_size_cg_trf_3d(im_1,nw,nh,nd,nomfichres);
	show_picture_3d(im_1);
	if (nomfichres)
		free(nomfichres);
}


/*******************************************************************************/
/*  									    	*/
/* -----imx_set_size_cg_trf_3d_p()------------------------------------------	*/
/*  									   	*/
/*   Toujours verifier avant d'utiliser cette fonction que nw, nh, nd sont soient*/
/*   plus grand ou plus petit que im1->width,im1->hght, im1->dpth 		*/
/*   Ne gere pas un melange. 							*/
/******************************************************************************/
void imx_set_size_cg_trf_3d_p(grphic3d *im1,unsigned int nw,unsigned int nh,unsigned int nd,transf3d *transfo)
{
    unsigned int i,j,k;
    unsigned int w1,h1,d1;
    unsigned int dw1,fw1,dh1,fh1,dd1,fd1;
    int flag=0;
    grphic3d *imtemp,*masktemp=NULL;
    grphic3d *maskim1=NULL;
    double *param;

    /* Controle */
    w1=im1->width;
    h1=im1->height;
    d1=im1->depth;
    flag= (nw-w1>=0) ? ((nh-h1>=0)?(((nd-d1)>=0?1:0)):0):((nh-h1<0)?(((nd-d1)<0?1:0)):0);
    flag= (nw>im1->mri_width)?0:(nh>im1->mri_height)?0:(nd>im1->mri_depth)?0:flag;
    if (flag==0) {
      PUT_ERROR("Nouvelle Taille incoherente");
      return;
      }

    if( ptr_has_mask_3d_p(im1) ) maskim1=ptr_mask_3d_p(im1);
    imtemp=cr_grphic3d(NULL);
    if( maskim1 != NULL ) masktemp=ptr_mask_3d_p(imtemp);
    w1=im1->width;
    h1=im1->height;
    d1=im1->depth;

    transfo->dx = im1->dx;
    transfo->dy = im1->dy;
    transfo->dz = im1->dz;
    param = transfo->param;

    if(nw>=w1||nh>=h1||nd>=d1) { /* agrandissement */
      dw1=(nw-w1)/2; fw1=nw-dw1;
      dh1=(nh-h1)/2; fh1=nh-dh1;
      dd1=(nd-d1)/2; fd1=nd-dd1;
    }
    else if(nw<=w1||nh<=h1||nd<=d1) { /* reduction */
      dw1=(w1-nw)/2; fw1=w1-dw1;
      dh1=(h1-nh)/2; fh1=h1-dh1;
      dd1=(d1-nd)/2; fd1=d1-dd1;
    }
    else { printf("PROBLEME dans set_size_cg_3d\n"); }

    param[0]=0.0 ;param[1]=0.0 ;param[2]=0.0 ; //rotation
    param[3]=-(dw1*im1->dx) ; param[4]=-(dh1*im1->dy); param[5]=-(dd1*im1->dz); //translation
    param[6]=0.0 ;param[7]=0.0 ;param[8]=0.0 ; //centre

    imx_apply_transf_3d_p(imtemp, im1, transfo, inter_linear_3d);
    if( maskim1 != NULL )
      imx_apply_transf_3d_p(masktemp, maskim1, transfo, inter_linear_3d);

    for(i=0;i<nw;i++)
     for(j=0;j<nh;j++)
       for(k=0;k<nd;k++)
    {
      im1->mri[i][j][k]=imtemp->mri[i][j][k];
      if( maskim1 != NULL ) maskim1->mri[i][j][k]=masktemp->mri[i][j][k];
    }

    im1->width=nw;
    im1->height=nh;
    im1->depth=nd;
    im1->x3dc=im1->width/2;
    im1->y3dc=im1->height/2;
    im1->z3dc=im1->depth/2;
    imx_inimaxminpixel_3d_p(im1);
    if (maskim1)
    {
      maskim1->width=nw;
      maskim1->height=nh;
      maskim1->depth=nd;
      maskim1->x3dc=im1->width/2;
      maskim1->y3dc=im1->height/2;
      maskim1->z3dc=im1->depth/2;
      imx_inimaxminpixel_3d_p(maskim1);
    }

    free_grphic3d(imtemp);
}


/****************************************************************************/
/*  									    */
/* -----imx_set_size_cg_3d_p()------------------------------------------*/
/*  									    */
/****************************************************************************/
void imx_set_size_cg_trf_3d(int im_1,int nw,int nh,int nd,char *nomfichres)
{
  grphic3d *im1;
  transf3d *transfo;
  
  im1=ptr_img_3d(im_1);
  transfo = cr_transf3d(nw,nh,nd,NULL);
  transfo->param = MALLOC(9,double);
  transfo->typetrans = RIGID3D;
  transfo->nb_param = 9;

  imx_set_size_cg_trf_3d_p(im1,nw,nh,nd,transfo);

  /*enregistrement du champ resultat*/
  if (nomfichres!=NULL)
  {
    save_transf_3d(transfo,nomfichres);
  }
  free_transf3d(transfo);
}




/****************************************************************************/
/*  									    */
/* -----set_x3dcy3dcz3dc_3d()------------------------------------------------*/
/*  									    */
/****************************************************************************/
void
set_x3dcy3dcz3dc_3d(void)
{
int im_1;
grphic3d *im1;
	
  im_1=GET_PLACE3D("Image ?" );

  im1=ptr_img_3d(im_1);
  imx_set_x3dcy3dcz3dc_3d_3d_p(im1);
  show_picture_3d(im_1);
}

/*******************************************************************************/
/*  									    	*/
/* -----imx_set_x3dcy3dcz3dc_3d_p()------------------------------------------	*/
/*  									   	*/
/*                                                                       */
/*    		                                                         */	
/*   							                        */
/******************************************************************************/
void imx_set_x3dcy3dcz3dc_3d_3d_p(grphic3d *im1)
{
	  
    im1->x3dc=im1->width/2;
    im1->y3dc=im1->height/2;
    im1->z3dc=im1->depth/2;

}

/****************************************************************************/
/*  									    */
/* -----equalize_slices_3d()------------------------------------------------*/
/*  									    */
/****************************************************************************/

void  equalize_slices_3d(void)
{
   int im_1, im_2;

   im_1=GET_PLACE3D("Image 1 ?" );
   im_2=GET_PLACE3D("Image 2 ?");
   
   imx_equalize_slices_3d(im_1,im_2);
   
   show_picture_3d(im_1);
   show_picture_3d(im_2);

  }

/****************************************************************************/
/*  									    */
/* -----imx_equalize_slices_3d_p()------------------------------------------*/
/*  									    */
/* Attention : Ne plus utiliser avec une image temporaire		    */ 
/*  	       cree avec cr_grphic3d car la taille n'est pas controle       */
/*  									    */
/****************************************************************************/

void imx_equalize_slices_3d_p(grphic3d *im1, grphic3d *im2)
{
   int i,j,k;
   int w1,h1,d1;
   int w2,h2,d2;
   int nw,nh,nd;


   w1=im1->width;   
   h1=im1->height;   
   d1=im1->depth;   
   w2=im2->width;   
   h2=im2->height;   
   d2=im2->depth;   

   nw=MAXI(w1,w2);
   nh=MAXI(h1,h2);
   nd=MAXI(d1,d2);

   im1->width=nw;
   im2->width=nw;
   im1->height=nh;
   im2->height=nh;
   im1->depth=nd;
   im2->depth=nd;



        for(i=w1;i<nw;i++)
           for(j=0;j<nh;j++)
              for(k=0;k<nd;k++) 
                 im1->mri[i][j][k]=0; 
  
        for(i=w2;i<nw;i++)
           for(j=0;j<nh;j++)
              for(k=0;k<nd;k++) 
                 im2->mri[i][j][k]=0; 
  
  
        for(j=h1;j<nh;j++)
           for(i=0;i<nw;i++)
              for(k=0;k<nd;k++) 
                 im1->mri[i][j][k]=0; 
  
        for(j=h2;j<nh;j++)
           for(i=0;i<nw;i++)
              for(k=0;k<nd;k++) 
                 im2->mri[i][j][k]=0; 
  
        for(k=d1;k<nd;k++) 
           for(i=0;i<nw;i++)
              for(j=0;j<nh;j++)
                 im1->mri[i][j][k]=0; 
  
        for(k=d2;k<nd;k++) 
           for(i=0;i<nw;i++)
              for(j=0;j<nh;j++)
                 im2->mri[i][j][k]=0; 

  imx_inimaxminpixel_3d_p(im1);
  imx_inimaxminpixel_3d_p(im2);

  
}

/****************************************************************************/
/*  									    */
/* -----imx_equalize_slices_3d()------------------------------------------*/
/*  									    */
/* Attention : Ne plus utiliser avec une image temporaire		    */ 
/*  	       cree avec cr_grphic3d car la taille n'est pas controle       */
/*  									    */
/****************************************************************************/

void imx_equalize_slices_3d(int im_1, int im_2)
{
   grphic3d *im1,*im2;


   im1=ptr_img_3d(im_1);
   im2=ptr_img_3d(im_2);
   
   imx_equalize_slices_3d_p(im1, im2);
}

/****************************************************************************/
/*  									    */
/* -----equalize_dxdydz_3d()------------------------------------------------*/
/*  									    */
/****************************************************************************/

void  equalize_dxdydz_3d(void)
{
   int im_1, im_2;
   grphic3d *im1,*im2;

   im_1=GET_PLACE3D("Adapt dx dy dz from this image ?" );
   im_2=GET_PLACE3D(TEXT0031);
   
   im1=ptr_img_3d(im_1);
   im2=ptr_img_3d(im_2);
   
   im1->dx=im2->dx;
   im1->dy=im2->dy;
   im1->dz=im2->dz;
   
  }
/****************************************************************************/
/*  									    */
/* -----imx_equalize_cutoff_3d()------------------------------------------------*/
/*  									    */
/****************************************************************************/

void imx_equalize_cutoff_3d(void)
{
   int im_1, im_2;
   grphic3d *im1,*im2;

   im_1=GET_PLACE3D("Equalize cuttoff:src  image" );
   im_2=GET_PLACE3D("Equalize cuttoff:dest image");
   
   im1=ptr_img_3d(im_1);
   im2=ptr_img_3d(im_2);
   
   im2->cutoff_max=(TYPEMRI3D)floor(im1->cutoff_max*im1->rcoeff/im2->rcoeff);
   im2->cutoff_min=(TYPEMRI3D)floor(im1->cutoff_min*im1->rcoeff/im2->rcoeff);

   show_picture_3d(im_2);   
}


/****************************************************************************/
/*  									    								*/
/* -----imx_adapte_taille_3d()----------------------------------------------*/
/*  									    								*/
/*  Copie de im_1 dans imres, imres est cree s'il est NULL, avec			*/
/*  adaptation de la taille de imres a wdth,hght,dpth.						*/
/*																			*/
/* Attention : Ne plus utiliser avec une image temporaire cree avec    		*/ 
/*  	        cr_grphic3d pour imres car la taille n'est pas controle     */
/*  									    								*/
/*  									    								*/
/*   im_1 : image a copier dans imres                                       */
/*   imres : Image resultat. Si NULL cr�tion d'une image temporaire       	*/
/*   wdth,hght,dpth : nouvelle taille de imres                                       */
/*  									    								*/
/*  									    								*/
/****************************************************************************/
grphic3d *imx_adapte_taille_3d(int im_1, grphic3d *imres, int wdth, int hght, int dpth)
{
   grphic3d *im1;

   im1=ptr_img_3d(im_1);

  return(imx_adapte_taille_3d_p(im1, imres, wdth, hght, dpth));
}


/*****************************************************/
/*!
**  Copie de im_1 dans imres,  avec  		 
**  adaptation de la taille de imres a wdth,hght,dpth.					 
**	\param im1 : image source
**	\param imres : image resultat (E/S)
**	\param wdth,hght,dpth : dimension
**	\retval : un pointeur sur l'image resultat
*****************************************************/
grphic3d *imx_adapte_taille_3d_p(grphic3d *im1, grphic3d *imres, int wdth, int hght, int dpth)
{
   int i,j,k;
   int w1,h1,d1;
   int nw,nh,nd;

   w1=im1->width;   
   h1=im1->height;   
   d1=im1->depth;   

   if (wdth < w1 || hght <h1 || dpth <d1 ) {
    PUT_ERROR("Cannot change the size in imx_adapte_taille_3d");
    return((grphic3d *)-1);
    }

   if ( imres==(grphic3d *)NULL ) {
   /*  Astuce pas tres propre. Mais ca marche. Bien verifier w1, h1,d1 */
     im1->width=wdth; im1->height=hght; im1->depth=dpth;
     imres=cr_grphic3d(im1);
     im1->width=w1; im1->height=h1; im1->depth=d1;
     }

   imx_copie_3d_p(im1,imres);
   imres->width=nw=wdth;
   imres->height=nh=hght;
   imres->depth=nd=dpth;

   for(i=w1;i<nw;i++)
     for(j=0;j<nh;j++)
       for(k=0;k<nd;k++) 
           imres->mri[i][j][k]=0; 
  
    for(j=h1;j<nh;j++)
      for(i=0;i<nw;i++)
         for(k=0;k<nd;k++) 
            imres->mri[i][j][k]=0; 
  
     for(k=d1;k<nd;k++) 
       for(i=0;i<nw;i++)
         for(j=0;j<nh;j++)
            imres->mri[i][j][k]=0; 
  

  imx_inimaxminpixel_3d_p(imres);

  return((grphic3d *)imres);

}




/*
** -- Rot_Dr_3d() -------------------------------------------------------------
**
** >> Last user action by: F. BOULEAU on July 23rd, 1999
**
*/

/******************************************************************************
** -- Miroir_3d() -------------------------------------------------------------
**
** >> Last user action by: B. MARAT on July  2002
**
******************************************************************************/
int Miroir_ttaxe_3d(int im_deb,int im_res,int axe)
{
  grphic3d *imdeb, *imres;

  imdeb = ptr_img_3d(im_deb);
  imres = ptr_img_3d(im_res);
  
  Miroir_ttaxe_3d_p(imdeb,imres,axe);
  
	show_picture_3d(im_res);
	return(axe);
}
/* Fonction effectuant un miroir sur l'axe "axe" (axe valant 1,2 ou 3, pour x,y et z respectivement),
 * avec des pointeurs en argument, sans affichage.
 * Hennequin 07/04/2005
 */
int Miroir_ttaxe_3d_p(grphic3d *imdeb,grphic3d *imres,int axe)
{
	grphic3d *imtmp;
	int i,j,k;
	int hght=imdeb->height,wdth=imdeb->width,dpth=imdeb->depth;

 	imtmp=cr_grphic3d(imdeb);

	for (i=0 ; i<wdth ; i++)
		for (j=0 ; j<hght ; j++)
			for (k=0 ; k<dpth ; k++)
				switch(axe)
				{
					case	1	:
						imtmp->mri[wdth-i-1][j][k]=imdeb->mri[i][j][k];
						break;
					case	2	:
						imtmp->mri[i][hght-j-1][k]=imdeb->mri[i][j][k];
						break;
					case	3	:
						imtmp->mri[i][j][dpth-k-1]=imdeb->mri[i][j][k];
						break;
					default	:	return 0;
				}
	
	imx_copie_param_3d_p(imtmp,imres);
	imx_copie_3d_p(imtmp,imres);
  	free_grphic3d(imtmp);
	
	return(axe);
}

int Miroir_3d_x(int im_deb,int im_res)	{return Miroir_ttaxe_3d(im_deb,im_res,1);}
int Miroir_3d_y(int im_deb,int im_res)	{return Miroir_ttaxe_3d(im_deb,im_res,2);}
int Miroir_3d_z(int im_deb,int im_res)	{return Miroir_ttaxe_3d(im_deb,im_res,3);}

void Miroir_3d(void)
{
	char *quest[5];
	int i,im_deb,im_res,retour=0;
	for (i=0;i<5;i++) 
		if ( (quest[i]=(char *)calloc(20,sizeof(char)))==NULL )
	 		{ printf("Erreur d'allocation : quest : Miroir_3d\n") ; 	PUT_ERROR("Erreur d'allocation");}

	im_deb=GET_PLACE3D(TEXT0001);
	im_res=GET_PLACE3D(TEXT0006);
	strcpy(quest[0],"l'Axe x") ;
	strcpy(quest[1],"l'Axe y") ;
	strcpy(quest[2],"l'Axe z") ;
	strcpy(quest[3],"Quitter") ;
	strcpy(quest[4],"\0") ;
	switch(GETV_QCM("Miroir selon :",(char **)quest) )
	{
		case	0	:		retour=Miroir_3d_x(im_deb,im_res);	break;
		case	1	:		retour=Miroir_3d_y(im_deb,im_res);	break;
		case	2	:		retour=Miroir_3d_z(im_deb,im_res);	break;
		case	3	:		retour=4;
		default		:	   retour=0;
	}
/*	
for(i=0;i<5;i++)
	free(quest[i]);
*/
if(retour!=0)
	printf("Action miroir correctement effectu� axe:%d\n",retour);
else
	printf("Erreur pendant l'action de miroir\n");
return;
}

/*
** -- Rot_Dr_3d() -------------------------------------------------------------
**
** >> Last user action by: F. BOULEAU on July 23rd, 1999
**
*/
void Rot_Dr_3d(void)
{
            int im_deb,im_res;

            im_deb=GET_PLACE3D(TEXT0001);
            im_res=GET_PLACE3D(TEXT0006);

            rot_droite_3d(im_deb,im_res);
        }

/******************************************************************************
** -- Rot_Ga_3d() -------------------------------------------------------------
**
** >> Last user action by: F. BOULEAU on July 23rd, 1999
**
******************************************************************************/

void Rot_Ga_3d(void)
{
            int im_deb,im_res;

            im_deb=GET_PLACE3D(TEXT0001);
            im_res=GET_PLACE3D(TEXT0006);
            rot_gauche_3d(im_deb,im_res);
        }

/******************************************************************************
** -- Rot_Ba_3d() -------------------------------------------------------------
**
** >> Last user action by: F. BOULEAU on July 23rd, 1999
**
******************************************************************************/

void Rot_Ba_3d(void)
{
            int im_deb,im_res;

            im_deb=GET_PLACE3D(TEXT0001);
            im_res=GET_PLACE3D(TEXT0006);
            rot_bas_3d(im_deb,im_res);
	}
	
/******************************************************************************
** -- Rot_Ha_3d() -------------------------------------------------------------
**
** >> Last user action by: F. BOULEAU on July 23rd, 1999
**
******************************************************************************/

void Rot_Ha_3d(void)
{
            int im_deb,im_res;

            im_deb=GET_PLACE3D(TEXT0001);
            im_res=GET_PLACE3D(TEXT0006);
            rot_haut_3d(im_deb,im_res);
        }

/******************************************************************************
** -- Rot_Pro_Dr_3d() ---------------------------------------------------------
**
** >> Last user action by: F. BOULEAU on July 23rd, 1999
**
******************************************************************************/

void Rot_Pro_Dr_3d(void)
{
            int im_deb,im_res;

            im_deb=GET_PLACE3D(TEXT0001);
            im_res=GET_PLACE3D(TEXT0006);
            rot_pro_dr_3d(im_deb,im_res);
        }

/******************************************************************************
** -- Rot_Pro_Gau_3d() ---------------------------------------------------------
**
** >> Last user action by: F. BOULEAU on July 23rd, 1999
**
******************************************************************************/

void Rot_Pro_Gau_3d(void)
{
            int im_deb,im_res;

            im_deb=GET_PLACE3D(TEXT0001);
            im_res=GET_PLACE3D(TEXT0006);
            rot_pro_gau_3d(im_deb,im_res);
        }

/******************************************************************************
** -- Rot_3d() ----------------------------------------------------------------
**
** >> Last user action by: F. BOULEAU on July 23rd, 1999
**
******************************************************************************/
/*
void Rot_3d(void)
{
            int im_deb,im_res;

            im_deb=GET_PLACE3D(TEXT0001);
            im_res=GET_PLACE3D(TEXT0006);
            rot_3d(im_deb,im_res);
        }*/

/******************************************************************************
** -- Move_3d_Inc() -----------------------------------------------------------
**
** >> Last user action by: F. BOULEAU on July 23rd, 1999
**
******************************************************************************/

void Move_3d_Inc(void)
{ 
            int im_deb,im_res;

            im_deb=GET_PLACE3D(TEXT0001);
            im_res=GET_PLACE3D(TEXT0006);
            decalage_3d(im_deb,im_res);
        }



/******************************************************************************
** -- Copie_3d() --------------------------------------------------------------
**
** >> Last user action by: F. BOULEAU on July 23rd, 1999
**
******************************************************************************/
void Copie_3d(void)
{
            int im_deb,im_res;

            im_deb=GET_PLACE3D(TEXT0001);
            im_res=GET_PLACE3D(TEXT0006);
            imx_copie_3d(im_deb,im_res);
            show_picture_3d(im_res);
        }

/******************************************************************************
** -- Copie_positif_3d() --------------------------------------------------------------
**
** >> Last user action by: F. BOULEAU on July 23rd, 1999
**
******************************************************************************/
void Copie_positif_3d(void)
{
            int im_deb,im_res;

            im_deb=GET_PLACE3D(TEXT0001);
            im_res=GET_PLACE3D(TEXT0006);
            imx_copie_positif_3d(im_deb,im_res);
            show_picture_3d(im_res);
        }

/******************************************************************************
** -- Copie_positif_min_3d() --------------------------------------------------------------
**
** >> Last user action by: JP ARMSAPCH 2004
**
******************************************************************************/
void Copie_positif_min_3d(void)
{
            int im_deb,im_res,wnd=0;

            im_deb=GET_PLACE3D(TEXT0001);
            im_res=GET_PLACE3D(TEXT0006);
	    wnd=imx_query_filter_dimension();
	      
            imx_copie_positif_min_3d(im_deb,im_res,wnd);
            show_picture_3d(im_res);
        }

/******************************************************************************
**  Copie_F1_3d()
** >> Last user action by: F. BOULEAU on July 23rd, 1999
*/
/*!
**
**	\brief copie une image et force le resultat a 1
**
******************************************************************************/
void Copie_F1_3d(void)
{ 
            int im_deb,im_res;

            im_deb=GET_PLACE3D(TEXT0001);
            im_res=GET_PLACE3D(TEXT0006);

            imx_copie_f1_3d(im_deb,im_res);
            show_picture_3d(im_res);
}

/******************************************************************************
**  Copie_imgmask_3d_p()
** >> Last user action by: JPA on Fevrier 2004
*/
/*!
**
**	\brief copie une image et le masque dans une autre image
**  \return 0 si succes
**
******************************************************************************/
int imx_Copie_imgmask_3d_p(grphic3d *imdeb,grphic3d *imres)
{
    if (!imdeb || !imres)
      return 1;

    imx_copie_3d_p(imdeb,imres);
    if (ptr_has_mask_3d_p(imdeb) )
        imx_copie_3d_p( ptr_mask_3d_p(imdeb),ptr_mask_3d_p(imres) );
    return 0;
}



/******************************************************************************
**  imx_Copie_imgmask_3d()
** >> Last user action by: JPA on Fevrier 2004
*/
/*!
**
**	\brief copie une image et le masque dans une autre image
**  \return 0 si succes
**
******************************************************************************/
int imx_Copie_imgmask_3d(int im_deb, int im_res)
{
    grphic3d *imdeb,*imres;

    imdeb = ptr_img_3d(im_deb);
    imres = ptr_img_3d(im_res);
    return imx_Copie_imgmask_3d_p(imdeb,imres);
}



/******************************************************************************
**  Copie_imgmask_3d()
** >> Last user action by: JPA on Fevrier 2004
*/
/*!
**
**	\brief copie une image et le masque dans une autre image
**
******************************************************************************/
void Copie_imgmask_3d(void)
{ 
            int im_deb,im_res;

            im_deb=GET_PLACE3D(TEXT0001);
            im_res=GET_PLACE3D(TEXT0006);
            imx_Copie_imgmask_3d(im_deb,im_res);
            show_picture_3d(im_res);
}

/******************************************************************************
** -- cent_3d() -----------------------------------------------------------
*/
/*!
**  \brief centre une image pour faire en sorte que le milieu de la composasnte
**  quelle contient (entour� de zero) soit au centre de l'image
**
******************************************************************************/
void cent_3d(void)
{ 
    int im_deb,im_res;
    long decal_x,decal_y,decal_z;

    im_deb=GET_PLACE3D(TEXT0001);
    im_res=GET_PLACE3D(TEXT0006);
    imx_cent_3d(im_deb,im_res,&decal_x,&decal_y,&decal_z);

    PUTI("Decalage X=",decal_x," _");
    PUTI("Y=",decal_y," _");
    PUTI("Z=",decal_z," ");

    show_picture_3d(im_res);
}

/******************************************************************************
** -- imx_cent_3d() -----------------------------------------------------------
*/
/*! centre une image pour faire en sorte que le milieu de la composasnte
**  quelle contient (entour� de zero) soit au centre de l'image
**  \param im_deb : numero de l'image source
**  \param im_res : numero de l'image dest
**  \retval 1 
******************************************************************************/
int imx_cent_3d(int im_deb, int im_res, long int *decal_x, long int *decal_y, long int *decal_z )
{ 
  grphic3d *imdeb,*imres;

  imdeb=ptr_img_3d(im_deb);
  imres=ptr_img_3d(im_res);
  imx_cent_3d_p(imdeb,imres,decal_x,decal_y,decal_z);
  return(1);	
}

/******************************************************************************
** -- imx_cent_3d_p(imdeb,imres) -----------------------------------------------------------
*/
/*!  centre une image pour faire en sorte que le milieu de la composasnte
**  quelle contient (entour� de zero) soit au centre de l'image
**  \param imdeb : image source
**  \param imres : image dest (E/S)
**  \retval 1
******************************************************************************/
int imx_cent_3d_p(grphic3d *imdeb, grphic3d *imres, long int *decal_x, long int *decal_y, long int *decal_z )
{
    int i,j,k,deci,decj,deck,imin,jmin,kmin,imax,jmax,kmax;
    int wdth,hght,dpth;

    wdth=imdeb->width;hght=imdeb->height;dpth=imdeb->depth;
    
    /*calcul du decalage �realiser*/
    imin=wdth;jmin=hght;kmin=dpth;
    imax=jmax=kmax=0;
    for(i=0;i<wdth;i++)
     for (j=0;j<hght;j++)
      for (k=0;k<dpth;k++)
       if ((imdeb->mri[i][j][k])!=0)
       {
        if (i<imin) imin=i;
	if (j<jmin) jmin=j;
	if (k<kmin) kmin=k;
	if (i>imax) imax=i;
	if (j>jmax) jmax=j;
	if (k>kmax) kmax=k;	
       }
    deci=(wdth-1-imax-imin)/2; 
    decj=(hght-1-jmax-jmin)/2; 
    deck=(dpth-1-kmax-kmin)/2; 
    
    imx_decalage_3d_p(imdeb,imres,-deci,-decj,-deck); 

    *decal_x=deci;
    *decal_y=decj;
    *decal_z=deck;

    imx_copie_param_3d_p(imdeb,imres);    
    return(1);
 }

/*******************************************
** --  bary_3d() ----------------------------
**
** Met l'image au barycentre
********************************************/
void    bary_3d(void)
{
    int im_1,im_res;
    long bary_x,bary_y,bary_z;
    long tx,ty,tz;
	grphic3d *im1;

    /* Question concerning type of operation  */
    /* and image                    */

    /*     Image                   */
    im_1=GET_PLACE3D(TEXT0021);
    im_res=GET_PLACE3D(TEXT0023);

    imx_bary_3d(im_1,im_res,&bary_x,&bary_y,&bary_z);

    im1=ptr_img_3d(im_1);
    PUTI("Barycentre X=",bary_x," _");
    PUTI("Y=",bary_y," _");
    PUTI("Z=",bary_z," _");
    tx=im1->width/2-bary_x;
    ty=im1->height/2-bary_y;
    tz=im1->depth/2-bary_z;
    PUTI("Decalage DX=",tx," _");
    PUTI("DY=",ty," _");
    PUTI("DZ=",tz," ");

    show_picture_3d(im_res);

}

/*******************************************
** --  imx_bary_3d() ----------------------------
**
** Met le barycentre de l'image au centre de la fenetre
**  im_deb : image de depart
**  im_res : image d'arrive
********************************************/
int    imx_bary_3d_p(grphic3d *imdeb, grphic3d * imres, long int *bary_x, long int *bary_y, long int *bary_z)
{
    int tx,ty,tz;

    imx_calcul_barycentre_de_img_3d_p(imdeb,bary_x,bary_y,bary_z);

    tx=imdeb->width/2-*bary_x;
    ty=imdeb->height/2-*bary_y;
    tz=imdeb->depth/2-*bary_z;
	
    imx_shift_bou_3d_p(imdeb,imres,tx,ty,tz);

    imx_copie_param_3d_p(imdeb,imres);

    return(1);
}

/*******************************************
** --  imx_bary_3d() ----------------------------
**
** Met le barycentre de l'image au centre de la fenetre
**  im_deb : image de depart
**  im_res : image d'arrive
********************************************/
int    imx_bary_3d(int im_deb, int im_res, long int *bary_x, long int *bary_y, long int *bary_z)
{
    grphic3d *imdeb,*imres;

    imdeb=ptr_img_3d(im_deb);
    imres=ptr_img_3d(im_res);

    imx_bary_3d_p(imdeb,imres,bary_x,bary_y,bary_z);
	
    return(1);
}

/*******************************************
** --  clean_board_3d() ----------------------------
**
** Met a zero le bord de l'image
********************************************/
void    clean_board_3d(void)
{
    int im_1,im_res,e=0;
    long nbpoints;

    /* Question concerning type of operation  */
    /* and image                    */

    /*     Image                   */
    im_1=GET_PLACE3D(TEXT0021);
    im_res=GET_PLACE3D(TEXT0023);
    nbpoints= GET_INT(TEXT0111, 0, &e);
    if (e) return;

    imx_clean_board_3d(im_1,im_res,nbpoints);

    show_picture_3d(im_res);

}

/*******************************************
** --  imx_clean_board_3d()  ---------------------
**
**  Met a zero le bord de l'images
**  im_deb : Image de depart
**  im_res : Image resultat
**  nbpts : nombre points mis a zero
*******************************************/
int    imx_clean_board_3d(int im_deb, int im_res, int nbpts)
{
   grphic3d *imdeb,*imres;

  imdeb=ptr_img_3d(im_deb);
  imres=ptr_img_3d(im_res);

  imx_clean_board_3d_p(imdeb,imres,nbpts);

  return(1);
}

/*******************************************
** --  imx_clean_board_3d_p()  ---------------------
**
**  Met a zero le bord de l'images
**  im_deb : Image de depart
**  im_res : Image resultat
**  nbpts : nombre points mis a zero
*******************************************/
int    imx_clean_board_3d_p(grphic3d *imdeb, grphic3d *imres, int nbpts)
{
  int i,j,k;
  int wdth,hght,dpth;

  wdth=imdeb->width;hght=imdeb->height;dpth=imdeb->depth;

  imx_copie_3d_p(imdeb,imres);

  for(i=0;i<wdth;i++)
	for (j=0;j<nbpts;j++)
	  for (k=0;k<dpth;k++)
	    imres->mri[i][j][k]=0;
        
  for(i=0;i<wdth;i++)
	for (j=hght-nbpts;j<hght;j++)
	  for (k=0;k<dpth;k++)
	    imres->mri[i][j][k]=0;
        
  for(i=0;i<nbpts;i++)
	for (j=0;j<hght;j++)
	  for (k=0;k<dpth;k++)
	    imres->mri[i][j][k]=0;
        
  for(i=wdth-nbpts;i<wdth;i++)
	for (j=0;j<hght;j++)
	  for (k=0;k<dpth;k++)
	    imres->mri[i][j][k]=0;

  for(i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	  for (k=0;k<nbpts;k++)
	    imres->mri[i][j][k]=0;
        
  for(i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	  for (k=dpth-nbpts;k<dpth;k++)
	    imres->mri[i][j][k]=0;
  
  return(1);
}


/******************************************************************************/
/*      energy_inter_images_3d() ----------------------------------------------
****************************************************************************/

void energie_inter_images_3d(void)
{
  int im_1,im_2;

  im_1=GET_PLACE3D("Image 1 ?");
  im_2=GET_PLACE3D("Image 2 ?");

  imx_energie_inter_images_3d(im_1,im_2);

 }

/****************************************************************************/

void imx_energie_inter_images_3d(int im_1, int im_2)
{
  grphic3d *im1,*im2;
  double energie;

  im1=ptr_img_3d(im_1);
  im2=ptr_img_3d(im_2);

  energie=imx_energie_inter_images_3d_p(im1,im2);

  PUTF("Energie ",(float)energie," ");
 }

/****************************************************************************/

 double imx_energie_inter_images_3d_p(grphic3d *im1, grphic3d *im2)
{
  double energie=0.0;
  int i,j,k;
  int val;
  int width,height,depth;

  width=im1->width;
  height=im1->height;
  depth=im1->depth;

  for(i=0;i<width;i++)
     for(j=0;j<height;j++)
        for(k=0;k<depth;k++)
	   {
	    val=(int)floor((im1->rcoeff)*(im1->mri[i][j][k])-(im2->rcoeff)*(im2->mri[i][j][k]));
	    energie+=(double)(val*val);
	   }

      energie=energie/((double)(width*height*depth));
      energie=sqrt(energie);

  return(energie);
 }

/******************************************************************************/
/*      energy_intra_image_3d() ----------------------------------------------
****************************************************************************/

void energie_intra_image_3d(void)
{
  int im_1;

  im_1=GET_PLACE3D("Image 1 ?");

  imx_energie_intra_image_3d(im_1);

 }

/****************************************************************************/

void imx_energie_intra_image_3d(int im_1)
{
  grphic3d *im1;
  double energie;

  im1=ptr_img_3d(im_1);

  energie=imx_energie_intra_image_3d_p(im1);

  PUTF("Energie ",(float)energie," ");
 }

/****************************************************************************/

 double imx_energie_intra_image_3d_p(grphic3d *im1)
{
  double energie=0.0;
  int i,j,k;
  int val;
  int width,height,depth;

  width=im1->width;
  height=im1->height;
  depth=im1->depth;

  for(i=0;i<width;i++)
     for(j=0;j<height;j++)
	for(k=0;k<depth;k++)
	   {
	    val=(int)floor((im1->rcoeff)*(im1->mri[i][j][k]));
	    energie+=(double)(val*val);
	   }

      energie=energie/((double)(width*height*depth));
      energie=sqrt(energie);

  return(energie);
 }

/********************************************************
**
**	imx_barycentre_binaire_3d_p ()
**  Calcul le barycentre d'une image binaire
**
**  img : pointeur image
**  seuil : 0 si img->mri<seuil
	    1 si img->mri>=seuil
**  x,y,z : (int) barycentre
**
******************************************************/
void imx_barycentre_binaire_3d_p(grphic3d *img, int seuil, float *x, float *y, float *z)
{
  int i,j,k;
  double xm=0.0,ym=0.0,zm=0.0,weight=0.0;
  int wdth,hght,dpth;

  wdth=img->width;
  hght=img->height;
  dpth=img->depth;

  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	if(img->mri[i][j][k]>=seuil)
             {
	      xm+=(double)i;
	      ym+=(double)j;
	      zm+=(double)k;
              weight+=1.0;
             }

  if(weight!=0.0) {
     *x=(float)(xm/weight);
     *y=(float)(ym/weight);
     *z=(float)(zm/weight);
    }
  else {
    PUT_WARN("The image is empty.\n Cannot calculate barycenter (x,y,z) of roi \n");
    *x=0.;
    *y=0.;
    *z=0.;
    }

 }

/*******************************************
**
** calcul le barycentre de l'image.
**
*******************************************/
void imx_calcul_barycentre_de_img_3d_p(grphic3d *im1, long int *x, long int *y, long int *z)
{
  int i,j,k;
  int w,h,d;
  double xm=0.0,ym=0.0,zm=0.0,weight=0.0;

  w = im1->width;
  h = im1->height;
  d = im1->depth;

  for(k=0;k<d;k++)
     for(j=0;j<h;j++)
        for(i=0;i<w;i++)
           {	
	    xm+=(double)((i+1))*(double)im1->mri[i][j][k];
	    ym+=(double)(j+1)*(double)im1->mri[i][j][k];
	    zm+=(double)(k+1)*(double)im1->mri[i][j][k];
            weight+=(double)im1->mri[i][j][k];
           }
 
    if(weight!=0.0) 
      {
       *x=(long int)(0.5+xm/weight);
       *y=(long int)(0.5+ym/weight);
       *z=(long int)(0.5+zm/weight);
      }
   
	  else
        {
        PUT_WARN("The image is empty.\n Cannot calculate barycenter (x,y,z) of img \n");
        *x=0;
        *y=0;
        *z=0;
        }
}

/***************************************************
**
**  imx_barycentre_pondere_3d_p ()
**  Calcul le barycentre d'une image en utilisant les valeur de mri
**
**  img : pointeur image
**  seuil : 0 si img->mri<seuil
	    1 si img->mri>=seuil
**  x,y,z : (int) barycentre
**
******************************************************/
void imx_barycentre_pondere_3d_p(grphic3d *img, int seuil, float *x, float *y, float *z)
{
  int i,j,k;
  double xm=0.0,ym=0.0,zm=0.0,weight=0.0;
  int wdth,hght,dpth;

  wdth=img->width;
  hght=img->height;
  dpth=img->depth;

  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	    if(img->mri[i][j][k]>=seuil)
          {
	      xm+=(double)(i+1)*img->mri[i][j][k];
	      ym+=(double)(j+1)*img->mri[i][j][k];
	      zm+=(double)(k+1)*img->mri[i][j][k];
          weight+=img->mri[i][j][k];
          }

  if(weight!=0.0) {
     *x=(float)(xm/weight);
     *y=(float)(ym/weight);
     *z=(float)(zm/weight);
    }
  else {
    PUT_WARN("The image is empty.\n Cannot calculate barycenter (x,y,z) of img \n");
    *x=0.;
    *y=0.;
    *z=0.;
    }

 }

/*******************************************
** --  imx_decalage_3d  ------------------------
**
**  Decalage de l'image 3d
**  im_deb : image de depart
**  im_res : image resultat
********************************************/
int   imx_decalage_3d(int im_deb, int im_res, int deci, int decj, int deck)
{
    grphic3d *imdeb,*imres;

    imdeb=ptr_img_3d(im_deb);
    imres=ptr_img_3d(im_res);

    imx_decalage_3d_p(imdeb,imres,deci,decj,deck);

    return(1);

}

/*******************************************
** --  imx_decalage_3d_p -------------------
*/
/*!  Decalage de l'image 3d
**  \param imdeb : pointeur sur image de depart
**  \param imres : pointeur sur image resultat (E/S)
**  \param deci : decalage vers la droite si deci>0 , vers la gauche sinon
**  \param decj : decalage vers le bas si decj>0 , vers le haut sinon
**  \param deck : ...
**  \retval 1
********************************************/
int    imx_decalage_3d_p(grphic3d *imdeb, grphic3d *imres, int deci, int decj, int deck)
{
    int i,j,k,ii,jj,kk;
    int hght,wdth,dpth;

    hght=imdeb->height;
    wdth=imdeb->width;
    dpth=imdeb->depth;
     
   if(deci>=0 && decj>=0 && deck>=0)
    {
        for (i=0;i<wdth;i++)
        {
            ii=i+deci;
            for (j=0;j<hght;j++)
            {
                jj=j+decj;
                for(k=0;k<dpth;k++)
                {
                    kk=k+deck;
                    if (ii<0||jj<0||kk<0||ii>=wdth||jj>=hght||kk>=dpth)
                        imres->mri[i][j][k]=0;
                    else
                        imres->mri[i][j][k]=imdeb->mri[ii][jj][kk];
                }
            }
        }
    }
    if(deci<0 && decj>=0 && deck>=0)
    {
        for(i=wdth-1;i>=0;i--)
        {
            ii=i+deci;
            for(j=0;j<hght;j++)
            {
                jj=j+decj;
                for(k=0;k<dpth;k++)
                {
                    kk=k+deck;
                    if (ii<0||jj<0||kk<0||ii>=wdth||jj>=hght||kk>=dpth)
                        imres->mri[i][j][k]=0;
                    else
                        imres->mri[i][j][k]=imdeb->mri[ii][jj][kk];
                }
            }
        }
    }
    if(deci>=0 && decj<0 && deck>=0)
    {
        for(i=0;i<wdth;i++)
        {
            ii=i+deci;
            for(j=hght-1;j>=0;j--)
            {
                jj=j+decj;
                for(k=0;k<dpth;k++)
                {
                    kk=k+deck;
                    if (ii<0||jj<0||kk<0||ii>=wdth||jj>=hght||kk>=dpth)
                        imres->mri[i][j][k]=0;
                    else
                        imres->mri[i][j][k]=imdeb->mri[ii][jj][kk];
                }
            }
        }
    }
    if(deci>=0 && decj>=0 && deck<0)
    {
        for(i=0;i<wdth;i++)
        {
            ii=i+deci;
            for(j=0;j<hght;j++)
            {
                jj=j+decj;
                for(k=dpth-1;k>=0;k--)
                {
                    kk=k+deck;
                    if (ii<0||jj<0||kk<0||ii>=wdth||jj>=hght||kk>=dpth)
                        imres->mri[i][j][k]=0;
                    else
                        imres->mri[i][j][k]=imdeb->mri[ii][jj][kk];
                }
            }
        }
    }
    if(deci<0 && decj<0 && deck>=0)
    {
        for(i=wdth-1;i>=0;i--)
        {
            ii=i+deci;
            for(j=hght-1;j>=0;j--)
            {
                jj=j+decj;
                for(k=0;k<dpth;k++)
                {
                    kk=k+deck;
                    if (ii<0||jj<0||kk<0||ii>=wdth||jj>=hght||kk>=dpth)
                        imres->mri[i][j][k]=0;
                    else
                        imres->mri[i][j][k]=imdeb->mri[ii][jj][kk];
                }
            }
        }
    }
    if(deci<0 && decj>=0 && deck<0)
    {
        for(i=wdth-1;i>=0;i--)
        {
            ii=i+deci;
            for(j=0;j<hght;j++)
            {
                jj=j+decj;
                for(k=dpth-1;k>=0;k--)
                {
                    kk=k+deck;
                    if (ii<0||jj<0||kk<0||ii>=wdth||jj>=hght||kk>=dpth)
                        imres->mri[i][j][k]=0;
                    else
                        imres->mri[i][j][k]=imdeb->mri[ii][jj][kk];
                }
            }
        }
    }
    if(deci>=0 && decj<0 && deck<0)
    {
        for(i=0;i<wdth;i++)
        {
            ii=i+deci;
            for(j=hght-1;j>=0;j--)
            {
                jj=j+decj;
                for(k=dpth-1;k>=0;k--)
                {
                    kk=k+deck;
                    if (ii<0||jj<0||kk<0||ii>=wdth||jj>=hght||kk>=dpth)
                        imres->mri[i][j][k]=0;
                    else
                        imres->mri[i][j][k]=imdeb->mri[ii][jj][kk];
                }
            }
        }
    }
    if(deci<0 && decj<0 && deck<0)
    {
        for(i=wdth-1;i>=0;i--)
        {
            ii=i+deci;
            for(j=hght-1;j>=0;j--)
            {
                jj=j+decj;
                for(k=dpth-1;k>=0;k--)
                {
                    kk=k+deck;
                    if (ii<0||jj<0||kk<0||ii>=wdth||jj>=hght||kk>=dpth)
                        imres->mri[i][j][k]=0;
                    else
                        imres->mri[i][j][k]=imdeb->mri[ii][jj][kk];
                }
            }
        }
    }

    return(1);

}


/*******************************************
**    --  decalage_3d()  ---------------------
**
**  Decalage de l'image 3d
**
********************************************/
int       decalage_3d(int im_dep, int im_res)
{
    char *quest[8];
    int answer_nr,i;
    int deci=0,decj=0,deck=0;

    /* Direction du decalage    */

    for(i=0;i<8;i++)
        quest[i]=CALLOC(80,char);
    strcpy(quest[0],TEXT0070);
    strcpy(quest[1],TEXT0071);
    strcpy(quest[2],TEXT0072);
    strcpy(quest[3],TEXT0073);
    strcpy(quest[4],TEXT0123);
    strcpy(quest[5],TEXT0124);
    strcpy(quest[6],TEXT0074);
    strcpy(quest[7],"\0");
    answer_nr=GETV_QCM(TEXT0043,quest);

    imx_copie_3d(im_dep,im_res);

    while(answer_nr!=6)
    {
        switch (answer_nr) {

        case 0:
            deci=-1;
            decj=0;
            deck=0;
            break;

        case 1:
            deci=1;
            decj=0;
            deck=0;
            break;

        case 2:
            deci=0;
            decj=-1;
            deck=0;
            break;

        case 3:
            deci=0;
            decj=1;
            deck=0;
            break;

        case 4:
            deci=0;
            decj=0;
            deck=-1;
            break;

        case 5:
            deci=0;
            decj=0;
            deck=1;
            break;
        }
        imx_decalage_3d(im_res,im_res,deci,decj,deck);
        show_picture_3d(im_res);
        answer_nr=GETV_QCM(TEXT0043,quest);
    }

    for(i=0;i<6;i++)
        FREE(quest[i]);

    return(1);

}


/*----------------------------------------------------------------------*/
/*   rescal3d								*/
/*! \brief  Rescale the grey_levels in the image to 0.. 255			*/
/*! \param im1resc : prt sur l'image (E/S)							*/
/*----------------------------------------------------------------------*/

void rescal3d(grphic3d *im1resc)
{   
    int	i, j,k;
    double scale;
    int w,h,d;

    scale = 255.0 / (double)(im1resc->max_pixel - im1resc->min_pixel);
      
    w = im1resc->width;
    h = im1resc->height;
    d = im1resc->depth;

    for (i=0; i<w; i++)
      for (j=0; j<h; j++) 
	for (k=0; k<d; k++) 
	im1resc->mri[i][j][k]= (TYPEMRI3D)floor(labs ((long)(scale * (im1resc->mri[i][j][k] - im1resc->min_pixel))));
 }



