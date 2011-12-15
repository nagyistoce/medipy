/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/******************************************************************************
 **/	
/*!	\file:		multiresolution_3d.c
***
***			
***
***	\brief description:    Recalage rigide avec multiresolution
***	
***	
***	Copyright (c) 2003, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***
********************************************************************************/

#include <config.h>
#include <string.h>

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/io/imx_log.h"
#include "noyau/imx_3d.h"
#include "noyau/gui/imx_picture3d.h"
#include "noyau/mani_3d.h"
#include "recalage/chps_3d.h"
#include "recalage/mtch_3d.h"
#include "math/imx_matrix.h"
#include "traitement/trai_3d.h"
#include "recalage/multiresolution_3d.h"

#undef		ERROR
#define		ERROR				TRUE

#define MAXF 200L		/* Maximum size of the filter */

#define SPLINE			"Spline"			/* Spline filter (l2-norm) */
#define SPLINE_L2		"Spline L2"			/* Spline filter (L2-norm) */
#define SPLINE_CENT		"Centered Spline"	/* Centered Spline filter (l2-norm) */
#define SPLINE_CENT_L2	"Centered Spline L2"/* Centered Spline filter (L2-norm) */

#define SIGMA_GAUSSIEN_RESAMPLING 0.425
#define TAILLE_GAUSSIENNE_RESAMPLING 3

#define MAX_TAILLE_IMAGE_RECALAGE 256

/******************************************************
 **
 */
//
/*!	\imx_ReSampleDimension_3d_p d'une image pour lui appliquer une dimension de voxel
//  passee en argument.
//
//  /param imsrc : image de depart
//  /param imdst : image resultat , si NULL une nouvelle image est allouee
//  /param dx : dx de l'image resultat
//  /param dy : dy de l'image resultat
//  /param dz : dz de l'image resultat
//  /param inter_type : type d'interpolation
//		0: NEAREST
1: LINEAR
2: SINC
3: QSINC2
4: QSINC3
5: BSPLINE2
6: BSPLINE3
7: BSPLINE4
8: BSPLINE5
9: LABEL
//  \param bFiltrage : not implemented
//  \return 0 si succes
//
// WARNING n'utiliser que pour une "reduction" !!!
//
//
*/
int imx_ReSampleDimension_3d_p(grphic3d *imsrc, grphic3d *imdst, float dx, float dy, float dz, int inter_type, int bFiltrage)
{
  grphic3d *retimg3d=NULL;
  transf3d *transf3d;
  InterpolationFct inter_fct;
  grphic3d *imsrc_filtered=NULL;
  float nZoomOut;
  int err;

  if (imdst == imsrc)
    retimg3d = cr_grphic3d(imsrc);
  else
    retimg3d = imdst;
  nZoomOut = dx/imsrc->dx;
  nZoomOut = (nZoomOut>dy/imsrc->dy) ? nZoomOut : dy/imsrc->dy;
  nZoomOut = (nZoomOut>dz/imsrc->dz) ? nZoomOut : dz/imsrc->dz;

  printf("ZoomOut dans filtrage = %f\n",nZoomOut);   
        
#if 0
  if (bFiltrage)
    imsrc_filtered=LowPassButerworth(imsrc,NULL,1.0/nZoomOut, 2);
  else
#endif /*0*/
    {
      imsrc_filtered=cr_grphic3d(imsrc);
      imx_copie_3d_p(imsrc,imsrc_filtered);
    }
  transf3d = cr_transf3d(imsrc->width,imsrc->height,imsrc->depth,NULL);
  transf3d->typetrans=RIGIDZOOM3D;
  transf3d->nb_param=12;
  transf3d->param=(double*)malloc(12*sizeof(double));
  transf3d->param[0]=0;transf3d->param[1]=0;transf3d->param[2]=0;
  transf3d->param[3]=imsrc->dx/dx;transf3d->param[4]=imsrc->dy/dy;transf3d->param[5]=imsrc->dz;
  transf3d->param[6]=0;transf3d->param[7]=0;transf3d->param[8]=0;
  transf3d->param[9]=imsrc->width/2;transf3d->param[10]=imsrc->height/2;transf3d->param[11]=imsrc->depth/2;

  inter_fct = imx_choose_interpolation_fct(inter_type);

  if (inter_fct == inter_VP_3d) // inter_VP_3d ne convient pas (pas une vrai interpolation)
    inter_fct = inter_linear_3d;

  err = imx_apply_transf_3d_p(retimg3d, imsrc_filtered, transf3d, inter_fct);
  retimg3d->width = (int)floor(retimg3d->width*retimg3d->dx/dx+0.5);
  retimg3d->height= (int)floor(retimg3d->height*retimg3d->dy/dy+0.5);
  retimg3d->depth = (int)floor(retimg3d->depth*retimg3d->dz/dz+0.5);
  retimg3d->dx    = dx;
  retimg3d->dy    = dy;
  retimg3d->dz    = dz;

  free_transf3d(transf3d);
  free_grphic3d(imsrc_filtered);
  if (imdst == imsrc)
    {
      imx_copie_3d_p(retimg3d,imdst);
      free_grphic3d(retimg3d);
    } 

  return 0;
}


/******************************************************
 **
 */
//
/*!	\briefResample d'une image d'un facteur 1/nZoomOut
//  ex : nZoomOut = 2 -> reduction par 2 de l'image
//
//  \param imsrc : image de depart
//  \param imdst : image resultat , si NULL une nouvelle image est allouee
//  \param inter_type : type d'interpolation 
//		0: NEAREST
1: LINEAR
2: SINC
3: QSINC2
4: QSINC3
5: BSPLINE2
6: BSPLINE3
7: BSPLINE4
8: BSPLINE5
9: LABEL
//  \param nZoomOut : facteur d'echelle
//  \param bFiltrage : not implemented
//  \return 0 si succes
//
*/
int imx_ReSample_3d_p(grphic3d *imsrc, grphic3d *imdst, int inter_type, float nZoomOut)
{
  grphic3d *retimg3d=NULL;
  transf3d *transf3d;
  InterpolationFct inter_fct;
  int err;
  TDimension newWidth, newHeight, newDepth;

  if (imdst == imsrc)
    retimg3d = cr_grphic3d(imsrc);
  else
    retimg3d = imdst;

  newWidth  = (TDimension)floor(imsrc->width/nZoomOut+0.5);
  newHeight = (TDimension)floor(imsrc->height/nZoomOut+0.5);
  newDepth  = (TDimension)floor(imsrc->depth/nZoomOut+0.5);
  
  //creation de la transformation a appliquer
  transf3d = cr_transf3d(newWidth,newHeight,newDepth,NULL);
  transf3d->dx=imsrc->dx*nZoomOut; transf3d->dy=imsrc->dy*nZoomOut; transf3d->dz=imsrc->dz*nZoomOut;
  transf3d->typetrans=RIGID3D;
  transf3d->nb_param=9;
  transf3d->param=MALLOC(9,double);
  transf3d->param[0]=0.0;transf3d->param[1]=0.0;transf3d->param[2]=0.0;
  transf3d->param[3]=0.0;transf3d->param[4]=0.0;transf3d->param[5]=0.0;
  transf3d->param[6]=0.0;transf3d->param[7]=0.0;transf3d->param[8]=0.0;

  inter_fct = imx_choose_interpolation_fct(inter_type);

  if (inter_fct == inter_VP_3d) // inter_VP_3d ne convient pas (pas une vrai interpolation)
    inter_fct = inter_linear_3d;

  //construction de la nouvelle image
  err = imx_apply_transf_3d_p(retimg3d, imsrc, transf3d, inter_fct);

  free_transf3d(transf3d);
  if (imdst == imsrc)
    {
      imx_copie_3d_p(retimg3d,imdst);
      free_grphic3d(retimg3d);
    }
  return 0;
}

int imx_ReSample_3d(int im_src, int im_dst, int inter_type, float nZoomOut)
{
  grphic3d *imsrc,*imdst;
  imsrc = ptr_img_3d(im_src);
  imdst = ptr_img_3d(im_dst);

  imx_ReSample_3d_p(imsrc, imdst, inter_type, nZoomOut);
  Refresh_3d();
	
  return 1;
	
}

void ReSample_3d()
{
  int im_src, im_dst;
  float zoom;
  int err=0, inter_type;
  im_src = GET_WND3D("image depart");
  im_dst = GET_WND3D("image resultat");
  inter_type = imx_query_interpolation_type_3D(0);
  zoom=GET_FLOAT("facteur de reduction?", 3.0, &err);
	
  if (err)
    return;
		
  imx_ReSample_3d(im_src,im_dst,inter_type,zoom);

  show_picture_3d(im_dst);

}

//--------------------------------------------------------------------------------------------------------//
//                             reduit une image avec un filtrage de Butterworth                           //
//--------------------------------------------------------------------------------------------------------//
/*
  int imx_ReSample_Butterworth_3d_p(grphic3d *imsrc, grphic3d *imdst, double nZoomOut)
  {
  grphic3d *imsrc_filtered=NULL;

  if (nZoomOut<=0.0) { fprintf(stderr,"facteur de zoom negatif ou nul dans imx_ReSample_Butterworth_3d_p\n"); return 1; }

  imsrc_filtered=cr_grphic3d(imsrc);
  LowPassButterworth(imsrc,imsrc_filtered,1.0/nZoomOut, 2);

  //subsampling en prenant un point sur deux
  imx_ReSample_3d_p(imsrc_filtered, imdst, 0, nZoomOut);

  free_grphic3d(imsrc_filtered);
  return 0;
  }

  int imx_ReSample_Butterworth_3d(int im_src, int im_dst, double nZoomOut)
  {
  grphic3d *imsrc,*imdst;
  imsrc = ptr_img_3d(im_src);
  imdst = ptr_img_3d(im_dst);

  imx_ReSample_Butterworth_3d_p(imsrc, imdst, nZoomOut);
  Refresh_3d();

  return 0;
  }

  void ReSample_Butterworth_3d()
  {
  int im_src, im_dst;
  float zoom;
  int err;

  im_src = GET_WND3D("image depart");
  im_dst = GET_WND3D("image resultat");
  zoom=GET_FLOAT(float,"facteur de reduction?","2",&err);

  imx_ReSample_Butterworth_3d(im_src,im_dst,zoom);
  }
*/
//--------------------------------------------------------------------------------------------------------//
//                             reduit une image avec un filtrage gaussien                                 //
//--------------------------------------------------------------------------------------------------------//

int imx_ReSample_gaussien_3d_p(grphic3d *imsrc, grphic3d *imdst, float nZoomOut)
{
  double sigma_gaussien = SIGMA_GAUSSIEN_RESAMPLING;
  grphic3d *imsrc_filtered=NULL;
  int taille=((int)floor(nZoomOut)/2)*2+1;
  //int taille=TAILLE_GAUSSIENNE_RESAMPLING;

  if (nZoomOut<=0.0) { fprintf(stderr,"facteur de zoom negatif ou nul dans imx_ReSample_gaussien_3d_p\n"); return 1; }

  imsrc_filtered=cr_grphic3d(imsrc);
  imx_filtre_gaussien_anisotrope_3d_p(imsrc, imsrc_filtered, taille, sigma_gaussien*nZoomOut, sigma_gaussien*nZoomOut, sigma_gaussien*nZoomOut);

  //subsampling en prenant un point sur deux
  imx_ReSample_3d_p(imsrc_filtered, imdst, 0, nZoomOut);

  free_grphic3d(imsrc_filtered);
  return 0;
}

int imx_ReSample_gaussien_3d(int im_src, int im_dst, float nZoomOut)
{
  grphic3d *imsrc,*imdst;
  imsrc = ptr_img_3d(im_src);
  imdst = ptr_img_3d(im_dst);

  imx_ReSample_gaussien_3d_p(imsrc, imdst, nZoomOut);
  Refresh_3d();

  return 0;
}

void ReSample_gaussien_3d()
{
  int im_src, im_dst;
  float zoom;
  int err;

  im_src = GET_WND3D("image depart");
  im_dst = GET_WND3D("image resultat");
  zoom=GET_FLOAT("facteur de reduction?", 2,&err);

  imx_ReSample_gaussien_3d(im_src,im_dst,zoom);
}

int imx_filtrage_gaussien_relatif_3d_p(grphic3d *imref, grphic3d *imsrc, grphic3d *imdst)
{
  double sigma_gaussien = SIGMA_GAUSSIEN_RESAMPLING;
  int taille=TAILLE_GAUSSIENNE_RESAMPLING;
  double fac_recax=0.0, fac_recay=0.0, fac_recaz=0.0;

  //on ne va filtrer l'image reca que dans les directions ou elle a une resolution plus fine que l'image ref filtree
  fac_recax=imref->dx/imsrc->dx; if (fac_recax<=1.0) fac_recax=0.0;
  fac_recay=imref->dy/imsrc->dy; if (fac_recay<=1.0) fac_recay=0.0;
  fac_recaz=imref->dz/imsrc->dz; if (fac_recaz<=1.0) fac_recaz=0.0;
  taille=((int)floor(MAXI(MAXI(fac_recax, fac_recay), fac_recaz))/2)*2+1;
  if ((fac_recax==0.0)&&(fac_recay==0.0)&&(fac_recaz==0.0)) imx_copie_3d_p(imsrc,imdst);
  else { imx_copie_param_3d_p(imsrc, imdst); imx_filtre_gaussien_anisotrope_3d_p(imsrc, imdst, taille, sigma_gaussien*fac_recax, sigma_gaussien*fac_recay, sigma_gaussien*fac_recaz); }

  return 0;
}

int imx_filtrage_gaussien_relatif_3d(int im_ref, int im_src, int im_dst)
{
  grphic3d *imref=NULL, *imsrc=NULL, *imdst=NULL;
  imsrc = ptr_img_3d(im_src);
  imdst = ptr_img_3d(im_dst);
  imref = ptr_img_3d(im_ref);

  imx_filtrage_gaussien_relatif_3d_p(imref, imsrc, imdst);
  Refresh_3d();

  return 0;
}

void filtrage_gaussien_relatif_3d()
{
  int im_ref, im_src, im_dst;

  im_ref = GET_WND3D("image de reference");
  im_src = GET_WND3D("image depart");
  im_dst = GET_WND3D("image resultat");

  imx_filtrage_gaussien_relatif_3d(im_ref, im_src, im_dst);
}

/*******************************************************************************
 **        pyramide_gaussienne_recalage_3d_p(imref, imreca, imrefres, imrecares, subsamplingfactor)
 */
/*!       redimensionne imref en le filtrant par une gaussienne de taille taille
**        et en en prenant un voxel tous les taille voxels
**        imreca est egalement filtree par le meme filtre
**
**	\param imref, imrefres: l'image de reference et le resultat de sa transformation
**	\param imreca, imrecares  : l'image a recaler et le resultat de sa transformation
**	\param taille : taille du filtre gaussien applique (!! elle doit etre impaire !!)
**	\retval 0
**
*******************************************************************************/
int pyramide_gaussienne_recalage_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imrefres, grphic3d *imrecares, double subsamplingfactor)
{
  if (subsamplingfactor>1.0) { imx_ReSample_gaussien_3d_p(imref, imrefres, (float)subsamplingfactor); imx_filtrage_gaussien_relatif_3d_p(imrefres,imreca,imrecares); }
  else { imx_copie_3d_p(imref,imrefres); imx_copie_3d_p(imreca,imrecares); }

  return 0;
}

/*******************************************************************************
 **        imx_reduire_image_recalage(im)
 */
/*!       reduit une image trop grosse a la taille maximale admise pour le recalage
**
**	\param im: l'image a reduire
**	\retval 0 en cas de succes
**
*******************************************************************************/
int imx_reduire_image_recalage(grphic3d *im)
{
  int new_width, new_height, new_depth;
  int width, height, depth;
  double fac_recax=0.0, fac_recay=0.0, fac_recaz=0.0;
  double sigma_gaussien = SIGMA_GAUSSIEN_RESAMPLING;
  int taille=TAILLE_GAUSSIENNE_RESAMPLING;
  grphic3d *im_tmp=NULL;
  transf3d *transfo=NULL;

  width=im->width; height=im->height; depth=im->depth;
  new_width=MINI(width, MAX_TAILLE_IMAGE_RECALAGE);
  new_height=MINI(height, MAX_TAILLE_IMAGE_RECALAGE);
  new_depth=MINI(depth, MAX_TAILLE_IMAGE_RECALAGE);

  if ((new_width==width)&&(new_height==height)&&(new_depth==depth)) return 0;

  printf ("l'image est trop volumineuse: on la reduit!!\n");
  im_tmp=cr_grphic3d_modif(new_width, new_height, new_depth, im->icomp, im->rcoeff, im->max_pixel);

  if (width>MAX_TAILLE_IMAGE_RECALAGE) fac_recax=(double)width/(double)MAX_TAILLE_IMAGE_RECALAGE; else fac_recax=0.0;
  if (height>MAX_TAILLE_IMAGE_RECALAGE) fac_recay=(double)height/(double)MAX_TAILLE_IMAGE_RECALAGE; else fac_recay=0.0;
  if (depth>MAX_TAILLE_IMAGE_RECALAGE) fac_recaz=(double)depth/(double)MAX_TAILLE_IMAGE_RECALAGE; else fac_recaz=0.0;

  //filtrage gaussien de l'image de depart
  taille=((int)floor(MAXI(MAXI(fac_recax, fac_recay), fac_recaz))/2)*2+1;
  if ((fac_recax>0.0)||(fac_recay>0.0)||(fac_recaz>0.0))
    { imx_filtre_gaussien_anisotrope_3d_p(im, im, taille, sigma_gaussien*fac_recax, sigma_gaussien*fac_recay, sigma_gaussien*fac_recaz); }

  //creation de la transformation zoom
  transfo=cr_transf3d(new_width, new_height, new_depth, NULL);
  if (width>MAX_TAILLE_IMAGE_RECALAGE) transfo->dx=im->dx*(float)width/(float)MAX_TAILLE_IMAGE_RECALAGE; else transfo->dx=im->dx;
  if (height>MAX_TAILLE_IMAGE_RECALAGE) transfo->dy=im->dy*(float)height/(float)MAX_TAILLE_IMAGE_RECALAGE; else transfo->dy=im->dy;
  if (depth>MAX_TAILLE_IMAGE_RECALAGE) transfo->dz=im->dz*(float)depth/(float)MAX_TAILLE_IMAGE_RECALAGE; else transfo->dz=im->dz;
  transfo->typetrans=RIGID3D;
  transfo->param=CALLOC(9, double);
  imx_apply_transf_3d_p(im_tmp, im, transfo, inter_nearest_3d);
  free_transf3d(transfo);

  //recopie de l'image transformee
  imx_inimaxminpixel_3d_p(im_tmp);
  resize_grphic3d_buffer(im, new_width, new_height, new_depth);
  imx_copie_3d_p(im_tmp, im);

  free_grphic3d(im_tmp);

  return 0;
}

//--------------------------------------------------------------------------------------------------------//
//    reduit une image par une puissance de 2 par la methode de reduction BSpline+moindres carres d'Unser       //
//--------------------------------------------------------------------------------------------------------//

/* ----------------------------------------------------------------------------
   Function: 	PyramidFilterSplinel2

   Purpose:  	Initializes down- and up-sampling filter arrays for
   least squares splines of order 0 to 3.  (little l_2 norm)
   g : reduce filter
   h : expand filter

   Author:		Michael Unser, NIH, BEIP, May 1992

   ---------------------------------------------------------------------------- */
void PyramidFilterSplinel2(double g[],long *ng,double *h,long *nh,long Order)
{

  switch (Order) {

  case 0L :
    *ng = 1L;
    *nh = 1L;
    break;

  case 1L :
    g[0]  =  0.707107;
    g[1]  =  0.292893;
    g[2]  = -0.12132;
    g[3]  = -0.0502525;
    g[4]  =  0.0208153;
    g[5]  =  0.00862197;
    g[6]  = -0.00357134;
    g[7]  = -0.0014793;
    g[8]  =  0.000612745;
    *ng = 9L;
    h[0]  = 1.;
    h[1]  = 0.5;
    *nh = 2L;
    break;
  case 2L :
    g[0]  =  0.617317;
    g[1]  =  0.310754;
    g[2]  = -0.0949641;
    g[3]  = -0.0858654;
    g[4]  =  0.0529153;
    g[5]  =  0.0362437;
    g[6]  = -0.0240408;
    g[7]  = -0.0160987;
    g[8]  =  0.0107498;
    g[9]  =  0.00718418;
    g[10] = -0.00480004;
    g[11] = -0.00320734;
    g[12] =  0.00214306;
    g[13] =  0.00143195;
    g[14] = -0.0009568;
    g[15] = -0.000639312;
    *ng = 16L;
    h[0]  =  1.;
    h[1]  =  0.585786;
    h[2]  =  0;
    h[3]  = -0.100505;
    h[4]  =  0;
    h[5]  =  0.0172439;
    h[6]  =  0;
    h[7]  = -0.00295859;
    h[8]  =  0;
    h[9]  =  0.000507614;
    *nh = 10L;
    break;
  case 3L :
    g[0]  =  0.596797;
    g[1]  =  0.313287;
    g[2]  = -0.0827691;
    g[3]  = -0.0921993;
    g[4]  =  0.0540288;
    g[5]  =  0.0436996;
    g[6]  = -0.0302508;
    g[7]  = -0.0225552;
    g[8]  =  0.0162251;
    g[9]  =  0.0118738;
    g[10] = -0.00861788;
    g[11] = -0.00627964;
    g[12] =  0.00456713;
    g[13] =  0.00332464;
    g[14] = -0.00241916;
    g[15] = -0.00176059;
    g[16] =  0.00128128;
    g[17] =  0.000932349;
    g[18] = -0.000678643;
    g[19] = -0.000493682;
    *ng = 20L;
    h[0]  =  1.;
    h[1]  =  0.600481;
    h[2]  =  0;
    h[3]  = -0.127405;
    h[4]  =  0;
    h[5]  =  0.034138;
    h[6]  =  0;
    h[7]  = -0.00914725;
    h[8]  =  0;
    h[9]  =  0.002451;
    h[10] =  0;
    h[11] = -0.000656743;
    *nh = 12L;
    break;
  default :
    *ng = -1L;
    *nh = -1L;
    break;
  }

}

/* ----------------------------------------------------------------------------
   Function: 	PyramidFilterSplineL2

   Purpose:  	Initializes down- and up-sampling filter arrays for
   L2 spline pyramid of order 0 to 5.
   g : reduce filter
   h : expand filter

   Author:		Michael Unser, NIH, BEIP, May 1992

   ---------------------------------------------------------------------------- */
void PyramidFilterSplineL2(double g[],long *ng,double h[],long *nh,long Order)
{

  switch (Order) {
  case 0L :
    *ng = 1L;
    *nh = 1L;
    break;

  case 1L :
    g[0]  =  0.683013;
    g[1]  =  0.316987;
    g[2]  = -0.116025;
    g[3]  = -0.0849365;
    g[4]  =  0.0310889;
    g[5]  =  0.0227587;
    g[6]  = -0.00833025;
    g[7]  = -0.00609817;
    g[8]  =  0.00223208;
    g[9]  =  0.001634;
    g[10] = -0.000598085;
    g[11] = -0.000437829;
    g[12] =  0.000160256;
    g[13] =  0.000117316;
    *ng = 14L;
    h[0]  =  1.;
    h[1]  =  0.5;
    *nh = 2L;
    break;

  case 3L :
    g[0]  =  0.594902;
    g[1]  =  0.31431;
    g[2]  = -0.0816632;
    g[3]  = -0.0942586;
    g[4]  =  0.0541374;
    g[5]  =  0.0454105;
    g[6]  = -0.0307778;
    g[7]  = -0.0236728;
    g[8]  =  0.0166858;
    g[9]  =  0.0125975;
    g[10] = -0.00895838;
    g[11] = -0.00673388;
    g[12] =  0.00479847;
    g[13] =  0.00360339;
    g[14] = -0.00256892;
    g[15] = -0.00192868;
    g[16] =  0.00137514;
    g[17] =  0.00103237;
    g[18] = -0.000736093;
    g[19] = -0.000552606;
    g[20] =  0.000394017;
    g[21] =  0.000295799;
    g[22] = -0.00021091;
    g[23] = -0.000158335;
    g[24] =  0.000112896;
    *ng = 25L;
    h[0]  =  1.;
    h[1]  =  0.600481;
    h[2]  =  0.0;
    h[3]  = -0.127405;
    h[4]  =  0;
    h[5]  =  0.034138;
    h[6]  =  0;
    h[7]  = -0.00914725;
    h[8]  =  0;
    h[9]  =  0.002451;
    h[10] =  0;
    h[11] = -0.000656743;
    *nh = 12L;
    break;

  case 5L :
    g[0]  =  0.564388;
    g[1]  =  0.316168;
    g[2]  = -0.0597634;
    g[3]  = -0.0998708;
    g[4]  =  0.0484525;
    g[5]  =  0.0539099;
    g[6]  = -0.0355614;
    g[7]  = -0.033052;
    g[8]  =  0.0246347;
    g[9]  =  0.0212024;
    g[10] = -0.0166097;
    g[11] = -0.0138474;
    g[12] =  0.0110719;
    g[13] =  0.00911006;
    g[14] = -0.00734567;
    g[15] = -0.0060115;
    g[16] =  0.00486404;
    g[17] =  0.00397176;
    g[18] = -0.00321822;
    g[19] = -0.00262545;
    g[20] =  0.00212859;
    g[21] =  0.00173587;
    g[22] = -0.0014077;
    g[23] = -0.0011478;
    g[24] =  0.000930899;
    g[25] =  0.000758982;
    g[26] = -0.000615582;
    g[27] = -0.000501884;
    g[28] =  0.000407066;
    g[29] =  0.000331877;
    g[30] = -0.00026918;
    g[31] = -0.000219459;
    g[32] =  0.000178;
    g[33] =  0.00014512;
    g[34] = -0.000117706;
    *ng = 35L;
    h[0]  =  1.;
    h[1]  =  0.619879;
    h[2]  =  0.0;
    h[3]  = -0.167965;
    h[4]  =  0;
    h[5]  =  0.0686374;
    h[6]  =  0;
    h[7]  = -0.0293948;
    h[8]  =  0.0;
    h[9]  =  0.0126498;
    h[10] =  0;
    h[11] = -0.00544641;
    h[12] =  0.0;
    h[13] =  0.00234508;
    h[14] =  0;
    h[15] = -0.00100973;
    h[16] =  0.0;
    h[17] =  0.000434766;
    h[18] =  0;
    h[19] = -0.000187199;
    *nh = 20L;
    break;

  default :
    *ng = -1L;
    *nh = -1L;
    aff_log( "Spline filters only defined for n=0,1,3,5");
    break;
  }
}

/* ----------------------------------------------------------------------------
   Function:	PyramidFilterCentered

   Purpose:	Initializes down- and up-sampling filter arrays for
   least squares CENTERED splines of order 0 to 4.  (little l_2 norm)
   g : reduce filter
   h : expand filter

   Note:		filter arrays should be defined as
   float g[20],h[20] filter arrays
   short *ng,*nh;	number of taps
   short Order;	order of the spline

   Author:		Patrick Brigger, NIH, BEIP	May 1996
   Daniel Sage, EPFL, Biomedical Imaging Group, November 1999

   ---------------------------------------------------------------------------- */
void PyramidFilterCentered(double g[],long *ng,double h[],long *nh,long Order)
{
  switch (Order) {
  case 0 :
    g[0] = 1;
    *ng=1;
    h[0] = 2;
    *nh=1;
    break;

  case 1 :
    g[0]  =  1.;
    g[1]  =  0.333333;
    g[2]  = -0.333333;
    g[3]  = -0.111111;
    g[4]  =  0.111111;
    g[5]  =  0.037037;
    g[6]  = -0.037037;
    g[7]  = -0.0123457;
    g[8]  =  0.0123457;
    g[9]  =  0.00411523;
    g[10] = -0.00411523;
    g[11] = -0.00137174;
    g[12] =  0.00137174;
    g[13] =  0.000457247;
    g[14] = -0.000457247;
    g[15] = -0.000152416;
    g[16] =  0.000152416;
    g[17] =  0.0000508053;
    g[18] = -0.0000508053;
    g[19] = -0.0000169351;
    g[20] =  0.0000169351;
    *ng = 21;
    h[0] =  1;
    h[1] =  0.5;
    *nh = 2;
    break;

  case 2 :
    g[0]  =  0.738417;
    g[1]  =  0.307916;
    g[2]  = -0.171064;
    g[3]  = -0.0799199;
    g[4]  =  0.0735791;
    g[5]  =  0.03108;
    g[6]  = -0.0307862;
    g[7]  = -0.0128561;
    g[8]  =  0.0128425;
    g[9]  =  0.00535611;
    g[10] = -0.00535548;
    g[11] = -0.00223325;
    g[12] =  0.00223322;
    g[13] =  0.000931242;
    g[14] = -0.00093124;
    g[15] = -0.000388322;
    g[16] =  0.000388322;
    g[17] =  0.000161928;
    g[18] = -0.000161928;
    g[19] = -0.0000675233;
    g[20] =  0.0000675233;
    *ng = 21;
    h[0]  =  1.20711;
    h[1]  =  0.585786;
    h[2]  = -0.12132;
    h[3]  = -0.100505;
    h[4]  =  0.0208153;
    h[5]  =  0.0172439;
    h[6]  = -0.00357134;
    h[7]  = -0.00295859;
    h[8]  =  0.000612745;
    h[9]  =  0.000507614;
    h[10] = -0.00010513;
    *nh = 11;
    break;

  case 3 :
    g[0]  =  0.708792;
    g[1]  =  0.328616;
    g[2]  = -0.165157;
    g[3]  = -0.114448;
    g[4]  =  0.0944036;
    g[5]  =  0.0543881;
    g[6]  = -0.05193;
    g[7]  = -0.0284868;
    g[8]  =  0.0281854;
    g[9]  =  0.0152877;
    g[10] = -0.0152508;
    g[11] = -0.00825077;
    g[12] =  0.00824629;
    g[13] =  0.00445865;
    g[14] = -0.0044582;
    g[15] = -0.00241009;
    g[16] =  0.00241022;
    g[17] =  0.00130278;
    g[18] = -0.00130313;
    g[19] = -0.000704109;
    g[20] =  0.000704784;
    *ng = 21;
    h[0]  =  1.13726;
    h[1]  =  0.625601;
    h[2]  = -0.0870191;
    h[3]  = -0.159256;
    h[4]  =  0.0233167;
    h[5]  =  0.0426725;
    h[6]  = -0.00624769;
    h[7]  = -0.0114341;
    h[8]  =  0.00167406;
    h[9]  =  0.00306375;
    h[10] = -0.000448564;
    h[11] = -0.000820929;
    h[12] =  0.000120192;
    h[13] =  0.000219967;
    h[14] = -0.0000322054;
    h[15] = -0.00005894;
    *nh = 16;
    break;

  case 4 :
    g[0]  =  0.673072;
    g[1]  =  0.331218;
    g[2]  = -0.139359;
    g[3]  = -0.12051;
    g[4]  =  0.086389;
    g[5]  =  0.0611801;
    g[6]  = -0.0542989;
    g[7]  = -0.034777;
    g[8]  =  0.033388;
    g[9]  =  0.0206275;
    g[10] = -0.0203475;
    g[11] = -0.0124183;
    g[12] =  0.0123625;
    g[13] =  0.00751369;
    g[14] = -0.00750374;
    g[15] = -0.00455348;
    g[16] =  0.00455363;
    g[17] =  0.00276047;
    g[18] = -0.00276406;
    g[19] = -0.00167279;
    g[20] =  0.00167938;
    *ng = 21;
    h[0]  =  1.14324;
    h[1]  =  0.643609;
    h[2]  = -0.0937888;
    h[3]  = -0.194993;
    h[4]  =  0.030127;
    h[5]  =  0.0699433;
    h[6]  = -0.0108345;
    h[7]  = -0.0252663;
    h[8]  =  0.00391424;
    h[9]  =  0.00912967;
    h[10] = -0.00141437;
    h[11] = -0.00329892;
    h[12] =  0.000511068;
    h[13] =  0.00119204;
    h[14] = -0.00018467;
    h[15] = -0.000430732;
    h[16] =  0.0000667289;
    h[17] =  0.000155641;
    h[18] = -0.0000241119;
    h[19] = -0.0000562395;
    *nh = 20;
    break;

  default :
    g[0] = 1.;
    *ng = 1;
    h[0] = 2.;
    *nh = 1;
    aff_log( "Spline filters only defined for n=0,1,2,3,4");
    break;
  }
}

/* ----------------------------------------------------------------------------
   Function:	PyramidFilterCenteredL2

   Purpose:	Initializes the symmetric down- and up-sampling filter arrays for
   L2 spline pyramid of order 0 to 5 when the downsampled grid is centered.
   These filters have then to be followed by a Haar filter.
   g: reduce filter
   h: expand filter

   Note:		filter arrays should be defined as
   float g[35],h[35]	filter arrays
   short *ng,*nh	number of taps
   short Order	order of the spline

   Author:		Patrick Brigger, NIH, BEIP,	April 1996
   Daniel Sage, EPFL, Biomedical Imaging Group, November 1999

   ---------------------------------------------------------------------------- */
void PyramidFilterCenteredL2(double g[],long *ng,double h[],long *nh,long Order)
{
  switch (Order) {
  case 0 :
    g[0] = 1.;
    *ng = 1;
    h[0] = 2.;
    *nh = 1;
    break;

  case 1 :
    g[0]  =  0.820272;
    g[1]  =  0.316987;
    g[2]  = -0.203044;
    g[3]  = -0.0849365;
    g[4]  =  0.0544056;
    g[5]  =  0.0227587;
    g[6]  = -0.0145779;
    g[7]  = -0.00609817;
    g[8]  =  0.00390615;
    g[9]  =  0.001634;
    g[10] = -0.00104665;
    g[11] = -0.000437829;
    g[12] =  0.000280449;
    g[13] =  0.000117316;
    g[14] = -0.000075146;
    g[15] = -0.0000314347;
    g[16] =  0.0000201353;
    *ng = 17;
    h[0] =  1.20096;
    h[1] =  0.473076;
    h[2] = -0.0932667;
    h[3] =  0.0249907;
    h[4] = -0.00669625;
    h[5] =  0.00179425;
    h[6] = -0.000480769;
    h[7] =  0.000128822;
    h[8] = -0.0000345177;
    *nh = 9;
    break;

  case 2 :
    g[0]  =  0.727973;
    g[1]  =  0.314545;
    g[2]  = -0.167695;
    g[3]  = -0.0893693;
    g[4]  =  0.0768426;
    g[5]  =  0.0354175;
    g[6]  = -0.0331015;
    g[7]  = -0.0151496;
    g[8]  =  0.0142588;
    g[9]  =  0.00651781;
    g[10] = -0.00613959;
    g[11] = -0.00280621;
    g[12] =  0.00264356;
    g[13] =  0.00120827;
    g[14] = -0.00113825;
    g[15] = -0.000520253;
    g[16] =  0.000490105;
    g[17] =  0.000224007;
    g[18] = -0.000211028;
    g[19] = -0.0000964507;
    g[20] =  0.0000908666;
    *ng = 21;
    h[0]  =  1.20711;
    h[1]  =  0.585786;
    h[2]  = -0.12132;
    h[3]  = -0.100505;
    h[4]  =  0.0208153;
    h[5]  =  0.0172439;
    h[6]  = -0.00357134;
    h[7]  = -0.00295859;
    h[8]  =  0.000612745;
    h[9]  =  0.000507614;
    h[10] = -0.00010513;
    *nh = 11;
    break;

  case 3 :
    g[0]  =  0.70222;
    g[1]  =  0.328033;
    g[2]  = -0.159368;
    g[3]  = -0.113142;
    g[4]  =  0.0902447;
    g[5]  =  0.0530861;
    g[6]  = -0.0492084;
    g[7]  = -0.0274987;
    g[8]  =  0.0264529;
    g[9]  =  0.0146073;
    g[10] = -0.0141736;
    g[11] = -0.0078052;
    g[12] =  0.00758856;
    g[13] =  0.00417626;
    g[14] = -0.00406225;
    g[15] = -0.00223523;
    g[16] =  0.00217454;
    g[17] =  0.00119638;
    g[18] = -0.00116412;
    g[19] = -0.000640258;
    g[20] =  0.000623379;
    *ng = 21;
    h[0]  =  1.15089;
    h[1]  =  0.623278;
    h[2]  = -0.0961988;
    h[3]  = -0.155743;
    h[4]  =  0.0259827;
    h[5]  =  0.041346;
    h[6]  = -0.0067263;
    h[7]  = -0.0112084;
    h[8]  =  0.00187221;
    h[9]  =  0.00296581;
    h[10] = -0.000481593;
    h[11] = -0.000805427;
    h[12] =  0.000134792;
    h[13] =  0.000212736;
    h[14] = -0.00003447;
    *nh = 15;
    break;

  case 4:
    g[0]  =  0.672101;
    g[1]  =  0.331667;
    g[2]  = -0.138779;
    g[3]  = -0.121385;
    g[4]  =  0.0864024;
    g[5]  =  0.0618776;
    g[6]  = -0.0545165;
    g[7]  = -0.0352403;
    g[8]  =  0.0335951;
    g[9]  =  0.0209537;
    g[10] = -0.0205211;
    g[11] = -0.0126439;
    g[12] =  0.0124959;
    g[13] =  0.0076682;
    g[14] = -0.00760135;
    g[15] = -0.00465835;
    g[16] =  0.00462238;
    g[17] =  0.00283148;
    g[18] = -0.00281055;
    g[19] = -0.00172137;
    g[20] =  0.00170884;
    *ng = 21;
    h[0]  =  1.14324;
    h[1]  =  0.643609;
    h[2]  = -0.0937888;
    h[3]  = -0.194993;
    h[4]  =  0.030127;
    h[5]  =  0.0699433;
    h[6]  = -0.0108345;
    h[7]  = -0.0252663;
    h[8]  =  0.00391424;
    h[9]  =  0.00912967;
    h[10] = -0.00141437;
    h[11] = -0.00329892;
    h[12] =  0.000511068;
    h[13] =  0.00119204;
    h[14] = -0.00018467;
    h[15] = -0.000430732;
    h[16] =  0.0000667289;
    h[17] =  0.000155641;
    h[18] = -0.0000241119;
    h[19] = -0.0000562396;
    *nh = 20;
    break;

  default :
    g[0] = 1.;
    *ng = 1;
    h[0] = 2.;
    *nh = 1;
    aff_log( "Spline filters only defined for n=0,1,2,3,4");
    break;
  }
}


/* ----------------------------------------------------------------------------
   Function:	PyramidFilterCenteredL2Derivate

   Purpose:	Initializes the symmetric down- and up-sampling filter arrays for
   L2 DERIVATIVE spline pyramid of order 0 to 5 when the downsampled
   grid is centered.
   These filters have then to be followed by a Derivative Haar filter.
   g : reduce filter
   h : expand filter
   Note:		filter arrays should be defined as
   float g[35],h[35]	filter arrays
   short *ng,*nh	number of taps
   short Order	order of the spline

   Author:		Patrick Brigger, NIH, BEIP,	April 1996
   Daniel Sage, EPFL, Biomedical Imaging Group, November 1999

   ---------------------------------------------------------------------------- */
void PyramidFilterCenteredL2Derivate(double g[],long *ng,double h[],long *nh,long Order)
{
  switch (Order) {
  case 0 :
    g[0] = 1.;
    *ng=1;
    h[0] = 2.;
    *nh=1;
    break;

  case 1 :
    g[0]  =  0.820272;
    g[1]  =  0.316987;
    g[2]  = -0.203044;
    g[3]  = -0.0849365;
    g[4]  =  0.0544056;
    g[5]  =  0.0227587;
    g[6]  = -0.0145779;
    g[7]  = -0.00609817;
    g[8]  =  0.00390615;
    g[9]  =  0.001634;
    g[10] = -0.00104665;
    g[11] = -0.000437829;
    g[12] =  0.000280449;
    g[13] =  0.000117316;
    g[14] = -0.000075146;
    g[15] = -0.0000314347;
    g[16] =  0.0000201353;
    *ng = 17;
    h[0]  =  1.20096;
    h[1]  =  1.20096;
    h[2]  = -0.254809;
    h[3]  =  0.068276;
    h[4]  = -0.0182945;
    h[5]  =  0.004902;
    h[6]  = -0.00131349;
    h[7]  =  0.000351947;
    h[8]  = -0.000094304;
    h[9]  =  0.0000252687;
    *nh = 10;
    break;

  case 2 :
    g[0]  =  0.727973;
    g[1]  =  0.314545;
    g[2]  = -0.167695;
    g[3]  = -0.0893693;
    g[4]  =  0.0768426;
    g[5]  =  0.0354175;
    g[6]  = -0.0331015;
    g[7]  = -0.0151496;
    g[8]  =  0.0142588;
    g[9]  =  0.00651781;
    g[10] = -0.00613959;
    g[11] = -0.00280621;
    g[12] =  0.00264356;
    g[13] =  0.00120827;
    g[14] = -0.00113825;
    g[15] = -0.000520253;
    g[16] =  0.000490105;
    g[17] =  0.000224007;
    g[18] = -0.000211028;
    g[19] = -0.0000964507;
    g[20] =  0.0000908666;
    *ng = 21;
    h[0]  =  1.20711;
    h[1]  =  0.585786;
    h[2]  = -0.12132;
    h[3]  = -0.100505;
    h[4]  =  0.0208153;
    h[5]  =  0.0172439;
    h[6]  = -0.00357134;
    h[7]  = -0.00295859;
    h[8]  =  0.000612745;
    h[9]  =  0.000507614;
    h[10] = -0.00010513;
    *nh = 11;
    break;

  case 3 :
    g[0]  =  0.70222;
    g[1]  =  0.328033;
    g[2]  = -0.159368;
    g[3]  = -0.113142;
    g[4]  =  0.0902447;
    g[5]  =  0.0530861;
    g[6]  = -0.0492084;
    g[7]  = -0.0274987;
    g[8]  =  0.0264529;
    g[9]  =  0.0146073;
    g[10] = -0.0141736;
    g[11] = -0.0078052;
    g[12] =  0.00758856;
    g[13] =  0.00417626;
    g[14] = -0.00406225;
    g[15] = -0.00223523;
    g[16] =  0.00217454;
    g[17] =  0.00119638;
    g[18] = -0.00116412;
    g[19] = -0.000640258;
    g[20] =  0.000623379;
    *ng = 21;
    h[0]  =  1.15089;
    h[1]  =  0.623278;
    h[2]  = -0.0961988;
    h[3]  = -0.155743;
    h[4]  =  0.0259827;
    h[5]  =  0.041346;
    h[6]  = -0.0067263;
    h[7]  = -0.0112084;
    h[8]  =  0.00187221;
    h[9]  =  0.00296581;
    h[10] = -0.000481593;
    h[11] = -0.000805427;
    h[12] =  0.000134792;
    h[13] =  0.000212736;
    h[14] = -0.00003447;
    *nh = 15;
    break;

  case 4:
    g[0]  =  0.672101;
    g[1]  =  0.331667;
    g[2]  = -0.138779;
    g[3]  = -0.121385;
    g[4]  =  0.0864024;
    g[5]  =  0.0618776;
    g[6]  = -0.0545165;
    g[7]  = -0.0352403;
    g[8]  =  0.0335951;
    g[9]  =  0.0209537;
    g[10] = -0.0205211;
    g[11] = -0.0126439;
    g[12] =  0.0124959;
    g[13] =  0.0076682;
    g[14] = -0.00760135;
    g[15] = -0.00465835;
    g[16] =  0.00462238;
    g[17] =  0.00283148;
    g[18] = -0.00281055;
    g[19] = -0.00172137;
    g[20] =  0.00170884;
    *ng = 21;
    h[0]  =  1.14324;
    h[1]  =  0.643609;
    h[2]  = -0.0937888;
    h[3]  = -0.194993;
    h[4]  =  0.030127;
    h[5]  =  0.0699433;
    h[6]  = -0.0108345;
    h[7]  = -0.0252663;
    h[8]  =  0.00391424;
    h[9]  =  0.00912967;
    h[10] = -0.00141437;
    h[11] = -0.00329892;
    h[12] =  0.000511068;
    h[13] =  0.00119204;
    h[14] = -0.00018467;
    h[15] = -0.000430732;
    h[16] =  0.0000667289;
    h[17] =  0.000155641;
    h[18] = -0.0000241119;
    h[19] = -0.0000562396;
    *nh = 20;
    break;

  default :
    g[0]  = 1.;
    *ng=1;
    h[0]  = 2.;
    *nh=1;
    aff_log( "Spline filters only defined for n=0,1,2,3,4");
    break;
  }
}

/* ----------------------------------------------------------------------------

Function:
GetPyramidFilter

Purpose:
Get the coefficients of the filter (reduce and expand filter)
Return the coefficients in g[ng] and in h[nh]

Convention:
g[ng] for the reduce filter
h[nh] for the expansion filter

Parameters:
Filter is the name of the filter

Order is the order for the filters based on splines
For the "Spline" filter, Order is 0, 1, 2 or 3
For the "Spline L2" filter, Order is 0, 1, 3 or 5
For the "Centered Spline" filter, Order is 0, 1, 2, 3 or 4
For the "Centered Spline L2" filter, Order is 0, 1, 2, 3 or 4

IsCentered is a return value indicates if the filter is a centered filter
TRUE if it is a centered filter
FALSE if it is not a centered filter

---------------------------------------------------------------------------- */
int GetPyramidFilter(
		     char *Filter,
		     long Order,
		     double g[], long *ng,
		     double h[], long *nh,
		     short *IsCentered)
{

  ng[0] = -1L;
  nh[0] = -1L;
  *IsCentered = FALSE;

  if ( !strcmp(Filter, "Spline"))	{
    PyramidFilterSplinel2(g, ng, h, nh, Order);
    *IsCentered = FALSE;
  }

  if ( !strcmp(Filter, "Spline L2")) {
    PyramidFilterSplineL2(g, ng, h, nh, Order);
    *IsCentered = FALSE;
  }

  if ( !strcmp(Filter, "Centered Spline")) {
    PyramidFilterCentered(g, ng, h, nh, Order);
    *IsCentered = TRUE;
  }

  if ( !strcmp(Filter, "Centered Spline L2"))	{
    PyramidFilterCenteredL2(g, ng, h, nh, Order);
    *IsCentered = TRUE;
  }

  if ( ng[0] == -1L && nh[0] == -1L) {
    aff_log( "This familly filters is unknown");
    return(ERROR);
  }
  return( !ERROR);

}

/* ----------------------------------------------------------------------------

Function:	ReduceCentered_1D

Purpose:	Reduces an image by a factor of two
The reduced image grid is between the finer grid

Parameters:
In[NxIn] is the input signal	(NxIn should be greater than 2 and even)
Out[NxIn/2] is the output signal
g[ng] is an array that contains the coefficients of the filter

Author:
Michael Unser, NIH, BEIP, June 1994
Patrick Brigger, NIH, BEIP,	May 1996, modified
Daniel Sage, EPFL, Biomedical Imaging Group, April 1999

---------------------------------------------------------------------------- */
void ReduceCentered_1D(	double In[], long NxIn,
			double Out[],
			double g[], long ng)
{
  double 	*y_tmp;
  long 	k, i, i1, i2;
  long	kk, kn, nred, n;

  nred = NxIn/2L;
  n = nred*2L;
  kn = 2L*n;

  /* --- Allocated memory for a temporary buffer --- */
  y_tmp = (double *)malloc( (size_t)(n*(long)sizeof(double)));
  if ( y_tmp == (double *)NULL) {
    aff_log("Out of memory in reduce_centered!");
    return;
  }

  /* --- Apply the symmetric filter to all the coefficients --- */
  for (k=0L; k<n; k++) {
    y_tmp[k] = In[k]*g[0 ];
    for (i=1L; i<ng; i++) {
      i1 = k-i;
      i2 = k+i;
      if (i1 < 0L) {
	i1 = (2L*n-1-i1) % kn;
	if (i1 >= n)
	  i1 = kn-i1-1L;
      }
      if (i2 > (n-1L)) {
	i2 = i2 % kn;
	if (i2 >= n) i2 = kn-i2-1L;
      }
      y_tmp[k] = y_tmp[k] + g[i]*(In[i1]+In[i2]);
    }
  }

  /* --- Now apply the Haar and perform downsampling --- */
  for(kk=0L; kk<nred; kk++) {
    k = 2L*kk;
    Out[kk] = (y_tmp[k] + y_tmp[k+1])/2.;
  }

  /* --- Free allocated memory --- */
  free(y_tmp);

}

/* ----------------------------------------------------------------------------

Function:
ReduceStandard_1D

Purpose:
Basic function to reduce a 1D signal

Parameters:
In[NxIn] is the input signal	(NxIn should be greater than 2 and even)
Out[NxIn/2] is the output signal
g[ng] is an array that contains the coefficients of the filter

Author:
Michael Unser, NIH, BEIP, June 1994
Daniel Sage, EPFL, Biomedical Imaging Group, April 1999

---------------------------------------------------------------------------- */
void ReduceStandard_1D(
		       double In[], long NxIn,
		       double Out[],
		       double g[], long ng)
{
  long k, i, i1, i2;
  long kk, kn, nred, n;

  nred = NxIn/2L;
  n  = nred*2L;
  kn = n-1L;			/* kn=n-2; DS Modified */

  if (ng<2L) {     	/* length filter < 2 */
    for (kk=0L; kk<nred; kk++) {
      k  = 2L*kk;
      i2 = k+1L;
      if (i2 > n-1L)
	i2 = kn-i2;
      Out[kk] = (In[k]+In[i2])/2.;
    }
  }

  else {
    for (kk=0L; kk<nred; kk++) {
      k = 2L*kk;
      Out[kk] = In[k]*g[0];
      for (i=1L; i<ng; i++) {
	i1 = k-i;
	i2 = k+i;
	if (i1<0L) {
	  i1 = (-i1) % kn;
	  if (i1 > n-1)
	    i1=kn-i1;
	}
	if (i2 > n-1L) {
	  i2 = i2 % kn;
	  if (i2 > n-1L)
	    i2=kn-i2;
	}
	Out[kk] = Out[kk] + g[i]*(In[i1]+In[i2]);
      }
    }
  }

}

/* ----------------------------------------------------------------------------

Function:
Reduce_1D

Purpose:
Router function to call ReduceStandard_1D or ReduceCentered_1D

---------------------------------------------------------------------------- */
void Reduce_1D(
	       double In[], long NxIn,
	       double Out[],
	       double g[], long ng,
	       short IsCentered
	       )
{

  if (IsCentered)
    ReduceCentered_1D( In, NxIn, Out, g, ng);
  else
    ReduceStandard_1D( In, NxIn, Out, g, ng);
}

/* ----------------------------------------------------------------------------

Function:
imx_ReSample_BSpline_3d_p

Purpose:
Reduces an image by a factor of power of two in corresponding dimensions.

Note:
Expects the output array (Out) to be allocated.

Parameters:
Input image:  	In[NxIn*NyIn]
Output image: 	Out[NxIn/2*NyIn/2]
nb of resolution to reduce: nbResol

----------------------------------------------------------------------------
                                                                              */

int imx_ReSample_BSpline_3d_p(grphic3d *imsrc, grphic3d *imdst, int nbResolX, int nbResolY, int nbResolZ)
{
  double ***tmp_res;
  double 	*InBuffer;		/* Input buffer to 1D process */
  double	*OutBuffer;		/* Output buffer to 1D process */
  unsigned int i,j,k;
  int l;
  TDimension wdth, hght, dpth;
  double  g[MAXF];					/* Coefficients of the reduce filter */
  long  	ng;							/* Number of coefficients of the reduce filter */
  double	h[MAXF];					/* Coefficients of the expansion filter */
  long 	nh;							/* Number of coefficients of the expansion filter */
  short	IsCentered;					/* Equal TRUE if the filter is a centered spline, FALSE otherwise */
/*   grphic3d *imtemp=NULL; */
  
  /* Get the filter coefficients for the Spline (order = 3) filter*/
  printf("Filter: " SPLINE_CENT_L2 "\n");
  if (GetPyramidFilter( SPLINE_CENT_L2, 3, g, &ng, h, &nh, &IsCentered) == ERROR) {
    printf("Unable to load the filter coeffiients");
    return 1;
  }

  /* --- Define dimension of the output --- */
  wdth=imsrc->width; hght=imsrc->height; dpth=imsrc->depth;

  /* --- Allocate a temporary image --- */
  // image temporaire en double
  tmp_res=alloc_dmatrix_3d(wdth, hght, dpth);
  if (tmp_res == (double ***)NULL) {
    aff_log("Unable to allocate memory");
    return(ERROR);
  }

  InBuffer=CALLOC(MAXI(wdth, MAXI(hght, dpth)), double);
  OutBuffer=CALLOC(MAXI(wdth, MAXI(hght, dpth)) / 2,double);

  if (OutBuffer == (double *)NULL) {
    free_dmatrix_3d(tmp_res);
    FREE(InBuffer);
    aff_log("Unable to allocate memory");
	    return(ERROR);
  }
  
  if (InBuffer == (double *)NULL) {
    free_dmatrix_3d(tmp_res);
    aff_log("Unable to allocate memory");
    return(ERROR);
  }



  for (i=0;i<wdth;i++)
    for (j=0;j<hght;j++)
      for (k=0;k<dpth;k++)
	tmp_res[i][j][k]=(double)imsrc->mri[i][j][k];

  for (l=0;l<nbResolX;l++)
    {
      /* --- X processing --- */
      if (wdth>1)
	{
	  for (j=0;j<hght;j++)
	    {
	      for (k=0;k<dpth;k++)
		{
		  for (i=0;i<wdth;i++) InBuffer[i]=tmp_res[i][j][k];
		  Reduce_1D(InBuffer, wdth, OutBuffer, g, ng, IsCentered);
		  for (i=0;i<wdth/2;i++) tmp_res[i][j][k]=OutBuffer[i];
		}
	    }
	}
      wdth=wdth/2;
    }

  for (l=0;l<nbResolY;l++)
    {
      /* --- Y processing --- */
      if (hght>1)
	{
	  for (i=0;i<wdth;i++)
	    {
	      for (k=0;k<dpth;k++)
		{
		  for (j=0;j<hght;j++) InBuffer[j] = tmp_res[i][j][k];
		  Reduce_1D(InBuffer, hght, OutBuffer, g, ng, IsCentered);
		  for (j=0;j<hght/2;j++) tmp_res[i][j][k] = OutBuffer[j];
		}
	    }
	}
      hght=hght/2;
    }

  for (l=0;l<nbResolZ;l++)
    {

      /* --- Z processing --- */
      if (dpth>1)
	{
	  for (i=0;i<wdth;i++)
	    {
	      for (j=0;j<hght;j++)
		{
		  for (k=0;k<dpth;k++) InBuffer[k] = tmp_res[i][j][k];
		  Reduce_1D(InBuffer, dpth, OutBuffer, g, ng, IsCentered);
		  for (k=0;k<dpth/2;k++) tmp_res[i][j][k] = OutBuffer[k];
		}
	    }
	}
      dpth=dpth/2;
    }
  
  //creation de l'image finale   
/*   if (imsrc == imdst) */
/*     { */
/*       imtemp = cr_grphic3d_modif(wdth, hght,dpth, 1.0, 1.0, 0); */
/*       if (!imtemp) */
/* 	{ fprintf (stderr, "erreur d'allocation memoire dans imx_compute_zoom_BSpline_3d_p\n"); goto end_func; } */
/*       imtemp->dx=(float)(imsrc->width/wdth)*imsrc->dx; */
/*       imtemp->dy=(float)(imsrc->height/hght)*imsrc->dy; */
/*       imtemp->dz=(float)(imsrc->depth/dpth)*imsrc->dz; */
/*       normaliser_3d(imsrc, imtemp, tmp_res, wdth, hght,dpth); */
/*       imx_copie_3d_p(imtemp, imdst); */
/*     } */
/*   else */
/*     { */
      imx_copie_param_3d_p(imsrc, imdst);
      imdst->width=wdth; imdst->height=hght; imdst->depth=dpth;
      imdst->dx=(float)(imsrc->width/wdth)*imsrc->dx;
      imdst->dy=(float)(imsrc->height/hght)*imsrc->dy;
      imdst->dz=(float)(imsrc->depth/dpth)*imsrc->dz;
      //normalisation
      normaliser_3d(imsrc, imdst, tmp_res, wdth, hght,dpth);
/*     } */

  imx_inimaxminpixel_3d_p(imdst);

/*  end_func: */
 
  /* --- Free the temporary image --- */
  free_dmatrix_3d(tmp_res);

  return (!ERROR);
}

/* ----------------------------------------------------------------------------

Function:
imx_ReSample_BSpline_3d

Purpose:
Reduces an image by a factor of power of two in corresponding dimensions.

Note:
Expects the output array (Out) to be allocated.

Parameters:
Input image:  	im_dst
Output image: 	im_dst
nbResol:			nb of resolution to reduce

---------------------------------------------------------------------------- */
int imx_ReSample_BSpline_3d(int im_src, int im_dst, int nbResolX, int nbResolY, int nbResolZ)
{
  grphic3d *imsrc,*imdst;
  imsrc = ptr_img_3d(im_src);
  imdst = ptr_img_3d(im_dst);

  imx_ReSample_BSpline_3d_p(imsrc, imdst, nbResolX, nbResolY, nbResolZ);
  Refresh_3d();

  return 0;
}

/* ----------------------------------------------------------------------------

Function:
ReSample_BSpline_3d

Purpose:
Reduces an image by a factor of power of two in corresponding dimensions.

Note:
Expects the output array (Out) to be allocated.

Parameters:

---------------------------------------------------------------------------- */
void ReSample_BSpline_3d()
{
  int im_src, im_dst;
  int zoomX,zoomY,zoomZ;
  int err=0;

  im_src = GET_WND3D("image depart");
  im_dst = GET_WND3D("image resultat");
  zoomX=GET_INT("puissance de 2 pour reduction X", 1, &err);
  zoomY=GET_INT("puissance de 2 pour reduction Y", 1, &err);
  zoomZ=GET_INT("puissance de 2 pour reduction Z", 1, &err);
  imx_ReSample_BSpline_3d(im_src,im_dst,zoomX,zoomY,zoomZ);
}

