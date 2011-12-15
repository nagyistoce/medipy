/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*******************************************************************************/	
/*!	   \file:		mtch_robust_3d.c
***
***		project:	Imagix 2.01
***			
***
***		\brief description:    Fichier pourle matching robust des images 3D
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

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_3d.h"
#include "recalage/mtch_robust_3d.h"

/*************************************************************************************/
/*! calcul de variables globales
	 - _moy_imref_3d, _moy_imreca_3d (moyenne des images)
	 - _ect_imref_3d, _ect_imreca_3d (ecart type) 
	 - _sigma_robust (sigam)
	 \param imref, imreca : images
	 \param energy_function : ne sert pas
	 \param robust_estimator : ne sert pas
	 \retval : fonction ne retourne rien, mais modifie les variables globales _moy_imref_3d,_moy_imreca_3d,_ect_imref_3d,_ect_imreca_3d,_sigma_robust
*/
void calc_robust_stat_3d(grphic3d *imref, grphic3d *imreca, int energy_function, int robust_estimator)
{
  double n_moy,n_sigma;
  int retour;
  double ect_imref_3d, ect_imreca_3d;

   _moy_imref_3d=calc_moy_image_3d(imref);
   _moy_imreca_3d=calc_moy_image_3d(imreca);
    
   ect_imref_3d=calc_ecart_type_image_3d(imref);
   ect_imreca_3d=calc_ecart_type_image_3d(imreca);
  
   retour=define_noise_param_3d(imref,imreca,energy_function,&n_moy,&n_sigma);
   
   if (retour==1) /* Il y a du bruit dans l'image */
     _sigma_robust=n_moy+10.0*n_sigma;
   else if (retour==2) /* Image segmente, on considere au pif
                              JPA Juillet 1999 */
      _sigma_robust=n_moy;

   if(RECA_DEBUG)
     {
     printf("IMREF  : moyenne = %f sigma = %f\n",_moy_imref_3d,ect_imref_3d);
     printf("IMRECA : moyenne = %f sigma = %f\n",_moy_imreca_3d,ect_imreca_3d);
     printf("NOISE  : moyenne = %f sigma = %f\n",n_moy,n_sigma);
     printf("ROBUST : sigma = %f\n",_sigma_robust);
     }
 }

/**************************************************/
/*! calcul la moyenne des valeurs de mri[][][] d'une image 3d
	\param image : image 3d
	\retval la moyenne
**********************************/ 
 double calc_moy_image_3d(grphic3d *image)
{
  int i,j,k,hght,wdth,dpth;
  double moy=0.0,size;

  wdth=image->width;
  hght=image->height;
  dpth=image->depth;
  size=(double)(hght*wdth*dpth);

  for(i=0;i<wdth;i++)
     for(j=0;j<hght;j++)
        for(k=0;k<dpth;k++)
           moy+=(double)image->mri[i][j][k]/size;


  return(moy);

 }

/**************************************************/
/*! calcul de l'ecart type des valeurs de mri[][][] d'une image 3d
	\param image : image 3d
	\retval l'ecart type
**********************************/ 
 double calc_ecart_type_image_3d(grphic3d *image)
{
  int i,j,k,hght,wdth,dpth;
  double moy=0.0,moy2=0.0,ect=0.0,size;

  wdth=image->width;
  hght=image->height;
  dpth=image->depth;
  size=(double)(hght*wdth*dpth);

  for(i=0;i<wdth;i++)
     for(j=0;j<hght;j++)
        for(k=0;k<dpth;k++) {
           moy+=(double)image->mri[i][j][k]/size;
           moy2+=(double)image->mri[i][j][k]*(double)image->mri[i][j][k]/size;
           }

  ect=sqrt(moy2-moy*moy);

  return(ect);

 }
/*******************************   
**  define_noise_param_3d  
*/
/*!  Recherche de la valeur moyenne et de l'ecty de la difference
**     dans le bruit 
**
**	\param  im1,im2 : images 
**	\param  energy_function : ne sert pas
**	\param  moy : moyenne (E/S) 
**	\param  sigma : ecart type
**
**	\retval 1 si bruit dans l'image \n
** 			2 si nmoys et sigma = moy et ecty de l'image
**	
**
**************************************************************/
int  define_noise_param_3d(grphic3d *im1, grphic3d *im2, int energy_function, double *moy, double *sigma)
{
  int i,j,k,w;
  int wdth,hght,dpth;
  double m=0.0,m2=0.0;
  double x;
  int pixels=0;

  w=(im1->width)/10;
  
  /* Recherche de la valeur moyenne et de l'ecty de la difference
     dans le bruit */
  for(i=0;i<w;i++)
     for(j=0;j<w;j++)
        for(k=0;k<w;k++)
               {
                x=(double)(im1->mri[i][j][k]-im2->mri[i][j][k]);
                m+=(double)(x/(double)(w*w*w));
	        m2+=(double)((x*x)/(double)(w*w*w));
               }

  
  if (m!=0) {
    *sigma=sqrt(m2-m*m);
    *moy=m;
    return (1);
    }
  else /* S'il n'y a pas de bruit dans l'image  on calcule la valeur moyenne
          et l'ecty de l'image  */
    {
    wdth=im1->width; hght=im1->height; dpth=im1->depth;
    m=0.0;m2=0.0;pixels=0;
    for(i=0;i<wdth;i++)
      for(j=0;j<hght;j++)
	    for(k=0;k<dpth;k++)
	      {
	      if (im1->mri[i][j][k]!=0)
	        {
	        m+=(double) im1->mri[i][j][k];
	        m2+= (double) im1->mri[i][j][k]*(double) im1->mri[i][j][k];
	        pixels++;
	        }
              } 
    m/=pixels;
    m2/=pixels;
    *sigma=sqrt(m2-m*m);
    *moy=m;
    return (2);
    }
}

