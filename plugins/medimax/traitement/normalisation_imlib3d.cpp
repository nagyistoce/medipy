/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include <config.h>

#ifdef HAS_IMLIB3D

#include<ImLib3D/CppTools.hpp>
#include<ImLib3D/Image3Dlinear.hpp>
#include<ImLib3D/Modules/SplineInterpolation.hpp>
#include<ImLib3D/TestPatterns.hpp>
#include<ImLib3D/Modules/Normalization.hpp>

#include "interface_imlib3d/ImLib3DImagixImageConversion.hpp"
#include "interface_imlib3d/LibImxCPP.hpp"
#include "interface_imlib3d/cpp_hooks.hpp"
#include "atlas/perfusion/MMInterface/NormalizationGMM.hpp"

extern "C"
{
#include	"noyau/imagix.h"
#include	"noyau/gtkimagix.h"
#include	"noyau/imx_3d.h"
#include    "noyau/imx_lang.h"
#include    "traitement/norm_3d.h"
#include "noyau/io/file_3d.h"

}
#endif



#ifdef HAS_IMLIB3D

/* ---------------------------------------------------------------- */
extern "C"
{
  void CPP_normalisation_histo(grphic3d *im1,grphic3d *im2,grphic3d *imres)
  {	unsigned int i,j,k;
  grphic3d *mask1,*mask2;
 
  // remplissage du masque
  mask1=ptr_mask_3d_p(im1);
  mask2=ptr_mask_3d_p(im2);
	
  for (i=0;i<im1->width;i++)
    for (j=0;j<im1->height;j++)
      for (k=0;k<im1->depth;k++)
	if (im1->mri[i][j][k]!=0)
	  mask1->mri[i][j][k]=1;
  mask1->width = im1->width;
  mask1->height = im1->height;
  mask1->depth = im1->depth;
	
  for (i=0;i<im2->width;i++)
    for (j=0;j<im2->height;j++)
      for (k=0;k<im2->depth;k++)
	if (im2->mri[i][j][k]!=0)
	  mask2->mri[i][j][k]=1;
  mask2->width = im2->width;
  mask2->height = im2->height;
  mask2->depth = im2->depth;
  // convert source arguments
  Image3Df im3D1,im3D2,im3Dres;
  ImagixToImLib3D(*im1,im3D1);
  ImagixToImLib3D(*im2,im3D2);
  IP3D::NormalizeJointHistogram(im3D1,im3D2,im3Dres,100,200);
  // convert result argument
  ImLib3DToImagix(im3Dres,imres);
  //on applique le masque
  for (i=0;i<imres->width;i++)
    for (j=0;j<imres->height;j++)
      for (k=0;k<imres->depth;k++)
	if (im1->mri[i][j][k]==0)
	  imres->mri[i][j][k]=0;
			
		
  //copie des dx,dy,dz
  imres->dx=im1->dx;
  imres->dy=im1->dy;
  imres->dz=im1->dz;
	
  }
}

/* ---------------------------------------------------------------- */

extern "C"
{
  void CPP_normalisation_GM(grphic3d *im1,grphic3d *im2,grphic3d *imres,int nb_gaussienne, double facteur)
  {	
 
#ifdef HAS_ATLAS_PERFUSION

    unsigned int i,j,k;
    grphic3d *mask1,*mask2;
 		double max;
		grphic3d *imtmp1,*imtmp2;
		
imx_inimaxminpixel_3d_p(im2);
imx_inimaxminpixel_3d_p(im1);

imtmp2=cr_grphic3d(im2);
imx_copie_param_3d_p(im2,imtmp2);

imtmp1=cr_grphic3d(im1);
imx_copie_param_3d_p(im1,imtmp1);

for (i=0;i<im1->width;i++)
	for (j=0;j<im1->height;j++)
	for (k=0;k<im1->depth;k++)
		{
		imtmp1->mri[i][j][k]=im1->mri[i][j][k];
		}

for (i=0;i<im2->width;i++)
	for (j=0;j<im2->height;j++)
	for (k=0;k<im2->depth;k++)
		{
		imtmp2->mri[i][j][k]=im2->mri[i][j][k];
		}
 

max=im2->max_pixel*im2->rcoeff;

if (max>1000)
	{
	printf("Rescaling des intensit�\n");
	for (i=0;i<im2->width;i++)
	for (j=0;j<im2->height;j++)
	for (k=0;k<im2->depth;k++)
		{
		im2->mri[i][j][k]=(int)(im2->mri[i][j][k]*im2->rcoeff*1000.0/max);
		}
	im2->max_pixel=1000;
	im2->rcoeff=1;
	}

imx_norm_seuil_meanecty_3d_p(im1,im2,im1);

/*	for (i=0;i<im2->width;i++)
	for (j=0;j<im2->height;j++)
	for (k=0;k<im2->depth;k++)
		{
		im1->mri[i][j][k]=(int)(1.0*im1->mri[i][j][k]*im1->rcoeff*im2->max_pixel/im1->max_pixel/im2->rcoeff);
		}
rcoeff1=im1->rcoeff;
im1->rcoeff=im2->rcoeff;
*/
	
				/*save_mri_ipb_3d_p("/home/noblet/im1.ipb",im1);
				save_mri_ipb_3d_p("/home/noblet/im2.ipb",im2);*/
				
    // remplissage du masque
    mask1=ptr_mask_3d_p(im1);
    mask2=ptr_mask_3d_p(im2);
	
    for (i=0;i<im1->width;i++)
    for (j=0;j<im1->height;j++)
		for (k=0;k<im1->depth;k++)
	  	if (im1->mri[i][j][k]!=0)
	    mask1->mri[i][j][k]=1;
    
		mask1->width = im1->width;
    mask1->height = im1->height;
    mask1->depth = im1->depth;
	
    for (i=0;i<im2->width;i++)
    for (j=0;j<im2->height;j++)
		for (k=0;k<im2->depth;k++)
	  if (im2->mri[i][j][k]!=0)
	    mask2->mri[i][j][k]=1;
    mask2->width = im2->width;
    mask2->height = im2->height;
    mask2->depth = im2->depth;
	

	
    // convert source arguments
    Image3Df im3D1,im3D2,im3Dres;
    ImagixToImLib3D(*im1,im3D1);
    ImagixToImLib3D(*im2,im3D2);
    IP3D::NormalizeGMM(im3D1,im3D2,im3Dres,nb_gaussienne,facteur);
    // convert result argument
    ImLib3DToImagix(im3Dres,imres);
    //on applique le masque
    for (i=0;i<imres->width;i++)
      for (j=0;j<imres->height;j++)
	for (k=0;k<imres->depth;k++)
	  if (im1->mri[i][j][k]==0)
	    imres->mri[i][j][k]=0;
			
    //copie des dx,dy,dz
    imres->dx=im1->dx;
    imres->dy=im1->dy;
    imres->dz=im1->dz;

imx_copie_param_3d_p(imtmp1,im1);
imx_copie_param_3d_p(imtmp2,im2);

if (max>1000)
	{
	for (i=0;i<imres->width;i++)
	for (j=0;j<imres->height;j++)
	for (k=0;k<imres->depth;k++)
		{
		imres->mri[i][j][k]=(int)(imres->mri[i][j][k]*max*imres->rcoeff/1000.0/im2->rcoeff);
		}
imres->rcoeff=im2->rcoeff;
	}

for (i=0;i<im2->width;i++)
	for (j=0;j<im2->height;j++)
	for (k=0;k<im2->depth;k++)
		{
		im2->mri[i][j][k]=imtmp2->mri[i][j][k];
		}
for (i=0;i<im1->width;i++)
	for (j=0;j<im1->height;j++)
	for (k=0;k<im1->depth;k++)
		{
		im1->mri[i][j][k]=imtmp1->mri[i][j][k];
		}
free_grphic3d(imtmp1);
free_grphic3d(imtmp2);

#endif		

  }
}

#endif /*!HAS_IMLIB3D*/
