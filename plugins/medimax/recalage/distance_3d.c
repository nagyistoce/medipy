/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*******************************************************************************/   
/*!    \file:       distance_3d.c
***
***     project:    Imagix 2.01
***         
***
***     \brief description:    Calcul de la distance interimage
*** 
*** 
*** Copyright (c) 1993, ULP-IPB Strasbourg.
*** All rights are reserved.
***
***
********************************************************************************/

#include <config.h>
#include <float.h>
#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_3d.h"
#include "math/imx_matrix.h"
#include "outils/imx_sort.h"
#include "recalage/chps_3d.h"
#include "recalage/distance_3d.h"
#include "recalage/mtch_3d.h"
#include "noyau/io/imx_log.h"
#include "noyau/mani_3d.h"
#include "recalage/erreur_recalage.h"

#define NB_CLASSES_WOODS_3D 250
//#define   NB_CLASSES_IM_3D 250
#define NB_CLASSES_IM_3D 256
#define NB_CLASSES_CR_3D 251
#define EPS 1e-6

/*******************************************************************************
********************************************************************************
******************************** DISTANCES *************************************
********************************************************************************
*******************************************************************************/ 
/*! \ingroup     DistanceFonction  @{ */
/*******************************************************************************
**     erreur_quad_3d(im1,im2)                                        
*/                                                                    
/*!    erreur quadratique moyenne entre deux images
        \param im1, im2 : les 2 images 
        \retval l'erreur quadratic
*******************************************************************************/
double erreur_quad_3d(ptr_distance_param dist_param)
{
 unsigned int i,j,k;
 TDimension wdth,hght,dpth;
 double r=DBL_MAX,diff;//,m1,m2;
 grphic3d *im1=NULL, *im2=NULL;
 TYPEMRI3D ***im1MRI=NULL, ***im2MRI=NULL;
 ERREUR_RECALAGE *err=NULL;

 im1=dist_param->imreca; im2=dist_param->imref;
 im1MRI=im1->mri; im2MRI=im2->mri;
 err=dist_param->err;

 wdth=im1->width;hght=im1->height;dpth=im1->depth;
 if ((wdth!=im2->width)||(hght!=im2->height)||(dpth!=im2->depth))
 {
  RAISE_ERR_TEXT(err, PB_TAILLES_IMAGES, "tailles d'images differentes dans erreur_quad_3d\n");
 }

 r=0.0;
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
    {
     diff=(double)(im1MRI[i][j][k]-im2MRI[i][j][k]);
     r=r+diff*diff;
    }
 r=sqrt(r);

  if (err!=NULL) *err=NO_ERR;
end_func:

 return(r);
}
/*! @} */

/*! \ingroup     DistanceFonction  @{ */
/*******************************************************************************
**     erreur_quad_robust_geman_3d(im1,im2)                                        
*/                                                                    
/*!    erreur quadratique moyenne robuste entre deux images.
**    Utilisation estimateur robuste Geman Mc Clure
        \param im1, im2 : les 2 images 
        \retval double : l'erreur quadratic
*******************************************************************************/
double erreur_quad_robust_geman_3d(ptr_distance_param dist_param)
 {
 unsigned int i,j,k;
 TDimension wdth,hght,dpth;
 double r=DBL_MAX,diff2,diff,w,weight;
 double sigma,sigma2;
 grphic3d *im1=NULL, *im2=NULL;
 TYPEMRI3D ***im1MRI=NULL, ***im2MRI=NULL;
 ERREUR_RECALAGE *err=NULL;

 im1=dist_param->imreca; im2=dist_param->imref;
 im1MRI=im1->mri; im2MRI=im2->mri;
 err=dist_param->err;

 wdth=im1->width;hght=im1->height;dpth=im1->depth;
 if ((wdth!=im2->width)||(hght!=im2->height)||(dpth!=im2->depth))
 {
  RAISE_ERR_TEXT(err, PB_TAILLES_IMAGES, "tailles d'images differentes dans erreur_quad_robust_geman_3d\n");
 }

 r=0.0;weight=0.0;
 sigma=_sigma_robust; sigma2=sigma*sigma;
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
    {
     diff=(double)(im1MRI[i][j][k]-im2MRI[i][j][k]);
     diff2=diff*diff;
     w=2.0*sigma*sigma/((sigma2+diff2)*(sigma2+diff2));
     weight+=w;
     r=r+diff2*w;
    }
 r/=weight;
 r=sqrt(r);

  if (err!=NULL) *err=NO_ERR;
end_func:

 return(r);
}
/*! @} */


/*******************************************************************************
**     old_erreur_woods_3d(im1,im2)                                         
**                                                                   
**    erreur de woods entre deux images (ancienne version olivier 
**   a supprimer apres accord 
*******************************************************************************/
double old_erreur_woods_3d(ptr_distance_param dist_param)
{
 unsigned int i,j,k;
 TDimension wdth,hght,dpth;
 int nbcl,maxpix,cl,nb;
 double r=DBL_MAX;
 double *moy=NULL,*var=NULL,*surf=NULL;
 grphic3d *im1=NULL, *im2=NULL;
 TYPEMRI3D ***im1MRI=NULL, ***im2MRI=NULL;
 ERREUR_RECALAGE *err=NULL;

 im1=dist_param->imreca; im2=dist_param->imref;
 im1MRI=im1->mri; im2MRI=im2->mri;
 err=dist_param->err;

 wdth=im1->width;hght=im1->height;dpth=im1->depth;
 if ((wdth!=im2->width)||(hght!=im2->height)||(dpth!=im2->depth))
 {
  RAISE_ERR_TEXT(err, PB_TAILLES_IMAGES, "tailles d'images differentes dans old_erreur_woods_3d\n");
 }

 //calcul du nombre de classes a utiliser
 nbcl=NB_CLASSES_WOODS_3D;
 maxpix=0;
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
    if (im1MRI[i][j][k]>maxpix) maxpix=im1MRI[i][j][k];

 if (nbcl>maxpix/2) nbcl=maxpix/2;

 if (nbcl<2)
 {
  RAISE_ERR_TEXT(err, NB_BINS_NUL, "nb de classes nulles dans old_erreur_woods_3d\n");
 }

 //allocations memoire
 moy=CALLOC(nbcl,double);
 var=CALLOC(nbcl,double);
 surf=CALLOC(nbcl,double);

 if ((!moy)||(!var)||(!surf))
 {
  RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans old_erreur_woods_3d\n");
 }

 //calcul des caracteristiques
 for (i=0;i<(unsigned int)nbcl;i++) {moy[i]=0.0;var[i]=0.0;surf[i]=0.0;}

 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    if (im1MRI[i][j][k]<=0) cl=0;
    else cl=(int)((double)(nbcl-1)*(double)(im1MRI[i][j][k])/(double)(maxpix+1))+1;
    surf[cl]+=1.0;
    moy[cl]+=(double)im2MRI[i][j][k];
   }

 for (cl=0;cl<nbcl;cl++)
  if (surf[cl]!=0.0)
   moy[cl]=moy[cl]/surf[cl];

 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
    {
     if (im1MRI[i][j][k]<=0) cl=0;
     else cl=(int)((double)(nbcl-1)*(double)(im1MRI[i][j][k])/(double)(maxpix+1))+1;
     if (cl<0) cl=0;
     if (cl>nbcl-1) cl=nbcl-1;
     var[cl]+=((double)im2MRI[i][j][k]-moy[cl])*((double)im2MRI[i][j][k]-moy[cl]);
    }

 r=0.0;nb=0;
 for (cl=0;cl<nbcl;cl++)
  if (surf[cl]!=0.0)
   {var[cl]=sqrt(var[cl])/surf[cl];r=r+var[cl];nb++;}

 r=r/(double)(nb);
  if (err!=NULL) *err=NO_ERR;
end_func:

 if (moy) FREE(moy);
 if (var) FREE(var);
 if (surf) FREE(surf);
 return(r);
}

/*! \ingroup DistanceFonction @{ */
/*******************************************************************************
**     erreur_woods_3d(im1,im2)                                         
*/                                                                    
/*!    erreur de woods entre deux images       
        \param im1, im2 : les 2 images 
        \retval double : l'erreur de woods
*******************************************************************************/
double erreur_woods_3d(ptr_distance_param dist_param)
{
 unsigned int i,j,k;
 TDimension wdth,hght,dpth;
 int nbcl,maxpix,cl;
 double r=DBL_MAX,valtemp;
 double *moy=NULL,*var=NULL,*surf=NULL;
 long all_points=0;
 grphic3d *im1=NULL, *im2=NULL;
 TYPEMRI3D ***im1MRI=NULL, ***im2MRI=NULL;
 ERREUR_RECALAGE *err=NULL;

 im1=dist_param->imreca; im2=dist_param->imref;
 im1MRI=im1->mri; im2MRI=im2->mri;
 err=dist_param->err;

 wdth=im1->width;hght=im1->height;dpth=im1->depth;
 if ((wdth!=im2->width)||(hght!=im2->height)||(dpth!=im2->depth))
 {
  RAISE_ERR_TEXT(err, PB_TAILLES_IMAGES, "tailles d'images differentes dans erreur_woods_3d\n");
 }

 //calcul du nombre de classes
 nbcl=NB_CLASSES_WOODS_3D;
 maxpix=im1->max_pixel;
 if (nbcl>maxpix) nbcl=maxpix+1;

 if (nbcl<2)
 {
  RAISE_ERR_TEXT(err, NB_BINS_NUL, "nb de classes nulles dans erreur_woods_3d\n");
 }

 //allocations memoire
 moy=CALLOC(nbcl,double);
 var=CALLOC(nbcl,double);
 surf=CALLOC(nbcl,double);

 if ((!moy)||(!var)||(!surf))
 {
  RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans erreur_woods_3d\n");
 }

//calcul des caracteristiques

 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    if (im1MRI[i][j][k]<=0) cl=0;
    else cl=(int)floor((double)(nbcl-1)*(double)(im1MRI[i][j][k])/(double)(maxpix)+0.5);

    if (cl > nbcl-1) cl=nbcl-1;
    surf[cl]+=1.0;
    valtemp=(double)im2MRI[i][j][k];
    moy[cl]+=valtemp;
    var[cl]+=valtemp*valtemp;
    all_points+=1;
   }

 for (cl=0;cl<nbcl;cl++)
  if (surf[cl]!=0.0) {
   moy[cl]=moy[cl]/surf[cl];
   var[cl]=var[cl]/surf[cl];
   var[cl]=var[cl]-moy[cl]*moy[cl];
   }

 r=0.0;
 for (cl=0;cl<nbcl;cl++)
  if (surf[cl]!=0.0)
   {r=r+sqrt(var[cl])*surf[cl]/all_points;}

  if (err!=NULL) *err=NO_ERR;
end_func:

 if (moy) FREE(moy);
 if (var) FREE(var);
 if (surf) FREE(surf);
 return(r);
}

/*! @} */

/*! \ingroup DistanceFonction  @{ */
/*******************************************************************************
**     erreur_woods_robust_3d(im1,im2)                                         
*/                                                                    
/*!    erreur de woods robuste entre deux images (cf. article de C. Nikou)     
        \param im1, im2 : les 2 images 
        \retval double : l'erreur de woods
*******************************************************************************/
double erreur_woods_robust_3d(ptr_distance_param dist_param)
{
 unsigned int i,j,k;
 TDimension wdth,hght,dpth;
 int nbcl=NB_CLASSES_WOODS_3D,maxpix,cl;
 double r=DBL_MAX;
 double diff,w,sigma,coeff;
 double *moy=NULL,*var=NULL,*surf=NULL,*weight=NULL;
 double *region[NB_CLASSES_WOODS_3D]={NULL};
 long all_points=0;
 grphic3d *im1=NULL, *im2=NULL;
 TYPEMRI3D ***im1MRI=NULL, ***im2MRI=NULL;
 ERREUR_RECALAGE *err=NULL;

 im1=dist_param->imreca; im2=dist_param->imref;
 im1MRI=im1->mri; im2MRI=im2->mri;
 err=dist_param->err;

 wdth=im1->width;hght=im1->height;dpth=im1->depth;
 if ((wdth!=im2->width)||(hght!=im2->height)||(dpth!=im2->depth))
 {
  RAISE_ERR_TEXT(err, PB_TAILLES_IMAGES, "tailles d'images differentes dans erreur_woods_robust_3d\n");
 }

 //calcul du nombre de classes
 maxpix=im1->max_pixel;
 if (nbcl>maxpix) nbcl=maxpix+1;

 if (nbcl<2)
 {
  RAISE_ERR_TEXT(err, NB_BINS_NUL, "nb de classes nulles dans erreur_woods_robust_3d\n");
 }

 //allocations memoire
 moy=CALLOC(nbcl,double);
 var=CALLOC(nbcl,double);
 surf=CALLOC(nbcl,double);
 weight=CALLOC(nbcl,double);

 if ((!moy)||(!var)||(!surf)||(!weight))
 {
  RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans erreur_woods_robust_3d\n");
 }

 for (i=0;i<(unsigned int)nbcl;i++)
 {
  region[i]=alloc_dvector(1);
  if (!region[i])
  {
   RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans erreur_woods_robust_3d\n");
  }
 }

 //calcul des caracteristiques

 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    if (im1MRI[i][j][k]<=0) cl=0;
    else cl=(int)floor((double)(nbcl-1)*(double)(im1MRI[i][j][k])/(double)(maxpix)+0.5);
    if (cl > nbcl-1) cl=nbcl-1;
    region[cl][(int) surf[cl]]=(double) im2MRI[i][j][k];
    surf[cl]+=1.0;
    all_points+=1;
    region[cl]=realloc_dvector(region[cl],(int) (surf[cl]+1.0));
   }

 sigma=10.;
 coeff=255./im2->max_pixel;

 for (cl=0;cl<nbcl;cl++)
  if (surf[cl]!=0.0) {
   moy[cl]=quick_select((int) floor(surf[cl]/2.0),region[cl],(int)floor(surf[cl]));
   for (j=0;j<(unsigned int)floor(surf[cl]);j++) {
     diff=fabs((region[cl][j]-moy[cl]))*coeff;
     w=2.0*diff*sigma*sigma/((sigma*sigma+diff*diff)*(sigma*sigma+diff*diff));
     weight[cl]+=w;
     var[cl]+=diff*diff*w;
     }
   }

 r=0.0;
 for (cl=0;cl<nbcl;cl++)
  if (surf[cl]!=0.0 && weight[cl]!=0.0)
   {r+=(surf[cl]/all_points)*(var[cl]/weight[cl]);}

  if (err!=NULL) *err=NO_ERR;
end_func:

 for (i=0;i<(unsigned int)nbcl;i++) if (region[i]) free_dvector(region[i],(int)floor(surf[i]));
 if (moy) FREE(moy);
 if (var) FREE(var);
 if (surf) FREE(surf);
 if (weight) FREE(weight);

 return(r);
}

/*! @} */

/*******************************************************************************
************** DISTANCES RELATIVES A L'HISTOGRAMME CONJOINT ********************
*******************************************************************************/

/*! \ingroup DistanceFonction  @{ */
/*******************************************************************************
**     erreur_IM_3d(im1,im2)                                         
*/                                                                    
/*!    erreur information mutuelle        
        \param im1, im2 : les 2 images 
        \retval double : l'erreur d'information mut.
*******************************************************************************/
double erreur_IM_3d(ptr_distance_param dist_param)
{
 unsigned int i,j,k;
 TDimension wdth,hght,dpth;
 int nbcl1=NB_CLASSES_IM_3D,nbcl2=NB_CLASSES_IM_3D,max1,max2,cl1,cl2;
 double **histo=NULL,*histo1=NULL,*histo2=NULL;
 double r=DBL_MAX,s;
 grphic3d *im1=NULL, *im2=NULL;
 TYPEMRI3D ***im1MRI=NULL, ***im2MRI=NULL;
 ERREUR_RECALAGE *err=NULL;

 im1=dist_param->imreca; im2=dist_param->imref;
 im1MRI=im1->mri; im2MRI=im2->mri;
 err=dist_param->err;

 wdth=im1->width;hght=im1->height;dpth=im1->depth;
 if ((wdth!=im2->width)||(hght!=im2->height)||(dpth!=im2->depth))
 {
  RAISE_ERR_TEXT(err, PB_TAILLES_IMAGES, "tailles d'images differentes dans erreur_IM_3d\n");
 }

 //calcul des maxpixel dans les deux images
 max1=im1->mri[0][0][0];max2=im2->mri[0][0][0];
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    if (im1MRI[i][j][k]>max1) max1=im1MRI[i][j][k];
    if (im2MRI[i][j][k]>max2) max2=im2MRI[i][j][k];
   }

 //calcul du nombre de classes
 if (nbcl1>max1) nbcl1=(max1+1);
 if (nbcl2>max2) nbcl2=(max2+1);

 if ((nbcl1<2)||(nbcl2<2))
 {
  RAISE_ERR_TEXT(err, NB_BINS_NUL, "nb de classes nulles dans erreur_IM_3d\n");
 }

 //allocations memoire
 histo1=CALLOC(nbcl1,double);histo2=CALLOC(nbcl2,double);
 histo=alloc_dmatrix(nbcl1,nbcl2);

 if ((!histo1)||(!histo2)||(!histo))
 {
  RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans erreur_IM_3d\n");
 }

 for (i=0;i<(unsigned int)nbcl1;i++)
  for (j=0;j<(unsigned int)nbcl2;j++)
   histo[i][j]=0.0;

 //calcul de l'histogramme conjoint et des histogrammes seuls
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    if (im1MRI[i][j][k]<=0) cl1=0;
    else cl1=(int)floor((double)(nbcl1-1)*(double)(im1MRI[i][j][k])/(double)(max1+1))+1;
    if (im2MRI[i][j][k]<=0) cl2=0;
    else cl2=(int)floor((double)(nbcl2-1)*(double)(im2MRI[i][j][k])/(double)(max2+1))+1;
    histo[cl1][cl2]=histo[cl1][cl2]+1.0;
    histo1[cl1]=histo1[cl1]+1.0;
    histo2[cl2]=histo2[cl2]+1.0;
   }

 s=(double)(wdth*hght*dpth);
 for (i=0;i<(unsigned int)nbcl1;i++) histo1[i]=histo1[i]/s;
 for (j=0;j<(unsigned int)nbcl2;j++) histo2[j]=histo2[j]/s;

 for (i=0;i<(unsigned int)nbcl1;i++)
  for (j=0;j<(unsigned int)nbcl2;j++)
  {
   histo[i][j]=histo[i][j]/s;
  }

 //calcul de la distance
 r=0.0;
 for (i=0;i<(unsigned int)nbcl1;i++)
  for (j=0;j<(unsigned int)nbcl2;j++)
  if (histo[i][j]!=0.0)
  r=r+histo[i][j]*log(histo[i][j]/histo1[i]/histo2[j]);

  r=100.0/r;

  if (err!=NULL) *err=NO_ERR;

end_func:

 if (histo) free_dmatrix(histo,nbcl1,nbcl2);
 if (histo1) FREE(histo1);
 if (histo) FREE(histo2);

 return(r);
}

/*! @} */

/*! \ingroup DistanceFonction  @{ */
/*******************************************************************************
**     erreur_entropie_conjointe_simple_3d(im1,im2)                                         
*/                                                                    
/*!    erreur entropie conjointe par sommation de voxels      
        \param im1, im2 : les 2 images 
        \retval double : l'erreur d'entropie conjointe.
*******************************************************************************/
double erreur_entropie_conjointe_simple_3d(ptr_distance_param dist_param)
{
 double entropie_conjointe=0.0;

 entropie_conjointe = calcul_info_mutuelle_3d(ENTROPIE_CONJOINTE, dist_param);

 return entropie_conjointe;
}

/*! @} */

/*! \ingroup DistanceFonction  @{ */
/*******************************************************************************
**     erreur_IM_simple_3d(im1,im2)                                         
*/                                                                    
/*!    erreur information mutuelle la plus simple a calculer      
        \param im1, im2 : les 2 images 
        \retval double : l'erreur d'information mut.
*******************************************************************************/
double erreur_IM_simple_3d(ptr_distance_param dist_param)
{
 double info_mutuelle=0.0;

 info_mutuelle = calcul_info_mutuelle_3d(IM, dist_param);

 return info_mutuelle;
}

/*! @} */

/*! \ingroup DistanceFonction  @{ */
/*******************************************************************************
**     erreur_IM_Normalisee_Studholme_simple_3d(im1,im2)                                         
*/                                                                    
/*!    erreur information mutuelle la plus simple a calculer      
        \param im1, im2 : les 2 images 
        \retval double : l'erreur d'information mut.
*******************************************************************************/
double erreur_IM_Normalisee_Studholme_simple_3d(ptr_distance_param dist_param)
{
 double info_mutuelle_normalisee=0.0;

 info_mutuelle_normalisee = calcul_info_mutuelle_3d(IMNS, dist_param);

 return info_mutuelle_normalisee;
}

/*! @} */

/*! \ingroup DistanceFonction  @{ */
/*******************************************************************************
**     erreur_IM_Normalisee_Maes1_simple_3d(im1,im2)                                         
*/                                                                    
/*!    erreur information mutuelle la plus simple a calculer      
        \param im1, im2 : les 2 images 
        \retval double : l'erreur d'information mut.
*******************************************************************************/
double erreur_IM_Normalisee_Maes1_simple_3d(ptr_distance_param dist_param)
{
 double info_mutuelle_normalisee=0.0;

 info_mutuelle_normalisee = calcul_info_mutuelle_3d(IMNM1, dist_param);

 return info_mutuelle_normalisee;
}

/*! @} */

/*! \ingroup DistanceFonction  @{ */
/*******************************************************************************
**     erreur_IM_Normalisee_Maes2_simple_3d(im1,im2)                                         
*/                                                                    
/*!    erreur information mutuelle la plus simple a calculer      
        \param im1, im2 : les 2 images 
        \retval double : l'erreur d'information mut.
*******************************************************************************/
double erreur_IM_Normalisee_Maes2_simple_3d(ptr_distance_param dist_param)
{
 double info_mutuelle_normalisee=0.0;

 info_mutuelle_normalisee = calcul_info_mutuelle_3d(IMNM2, dist_param);

 return info_mutuelle_normalisee;
}

/*! @} */


/*! \ingroup DistanceFonction  @{ */
/*******************************************************************************
**     erreur_CR_3d(im1,im2)                                         
*/                                                                    
/*!    erreur rapport de correlation de Ayache        
        \param im1, im2 : les 2 images 
        \retval double : l'erreur dde rapport de correlation
*******************************************************************************/
double erreur_CR_3d(ptr_distance_param dist_param)
{
 int i,j,k;
 int wdth,hght,dpth;
 int nbcl1,nbcl2,max1,max2,/*cl1,*/cl2, nbcl, max;
 double cl1;
 double nbVoxelsCommuns=0.0, moyenneTot=0.0, varianceTot=0.0;
 double *nbVoxelsParClasse=NULL, *moyennesParClasse=NULL, *variancesParClasse=NULL;
 double cr=DBL_MAX;
 grphic3d *im1=NULL, *im2=NULL;
 TYPEMRI3D ***im1MRI=NULL, ***im2MRI=NULL;
 ERREUR_RECALAGE *err=NULL;

 im1=dist_param->imreca; im2=dist_param->imref;
 im1MRI=im1->mri; im2MRI=im2->mri;
 err=dist_param->err;

 wdth=(int)im1->width;hght=(int)im1->height;dpth=(int)im1->depth;
 if ((wdth!=(int)im2->width)||(hght!=(int)im2->height)||(dpth!=(int)im2->depth))
 {
  RAISE_ERR_TEXT(err, PB_TAILLES_IMAGES, "tailles d'images differentes dans erreur_CR_3d\n");
 }

 //calcul des maxpixel dans les deux images
 //et du nombre de classes de valeurs de voxels qu'on va utiliser
// nbcl1=imx_calculer_nombre_de_classes_3d_p(im1, &max1);
// nbcl2=imx_calculer_nombre_de_classes_3d_p(im2, &max2);

 nbcl1=dist_param->nbclReca; nbcl2=dist_param->nbclRef;
 max1=dist_param->maxReca; max2=dist_param->maxRef;

 if ((nbcl1<2)||(nbcl2<2))
 {
  RAISE_ERR_TEXT(err, NB_BINS_NUL, "nb de classes nulles dans erreur_CR_3d\n");
 }

 //allocations memoires et initialisation
 nbVoxelsParClasse = CALLOC(nbcl2,double);
 moyennesParClasse = CALLOC(nbcl2,double);
 variancesParClasse = CALLOC(nbcl2,double);

 if ((!nbVoxelsParClasse)||(!moyennesParClasse)||(!variancesParClasse))
 {
  RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans erreur_CR_3d\n");
 }

 nbcl=nbcl2-1; max=max2+1;
 //calcul du rapport de correlation
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    //on fait une classe a part pour les intensites nulles qui ne sont porteuses d'aucune information
    if (im1MRI[i][j][k]<=0) cl1=0.0;
    //else cl1=(int)floor((double)(nbcl1-1)*(double)(im1->mri[i][j][k])/(double)(max1+1))+1;
    else cl1= (double)im1MRI[i][j][k];
    if (im2MRI[i][j][k]<=0) cl2=0;
    else cl2=MINI((nbcl*im2MRI[i][j][k])/max+1, nbcl);

      //on ne fait le calcul que sur la zone de recouvrement des 2 images
      if ((cl1>EPS)&&(cl2>0))
      {
       nbVoxelsCommuns+=1.0;
       moyenneTot += cl1;
       varianceTot += (cl1*cl1);
       nbVoxelsParClasse[cl2] += 1.0;
       moyennesParClasse[cl2] += cl1;
       variancesParClasse[cl2] += (cl1*cl1);
      }

   }

 //si il n'y a pas de recouvrement des 2 images on renvoie la valeur max
 if (nbVoxelsCommuns<EPS) return 1.0;

 moyenneTot = moyenneTot/nbVoxelsCommuns;
 varianceTot = (varianceTot/nbVoxelsCommuns) - (moyenneTot*moyenneTot);

 //si la variance totale est 0 ce qui ne se produit que lorsque tous les pixels de imref
 //ont uniformement pour classe 1 ou 0 pour les differentes classes de imreca dans la zone de recouvrement
 //on ne peut en tirer aucune information et on renvoie la valeur max
 if (varianceTot<EPS) return 1.0;

 //on ne fait le calcul que sur la zone de recouvrement des 2 images
 for (i=1;i<nbcl2;i++)
 {
  if (nbVoxelsParClasse[i]>0)
  {
   moyennesParClasse[i] = moyennesParClasse[i]/nbVoxelsParClasse[i];
   variancesParClasse[i] = (variancesParClasse[i]/nbVoxelsParClasse[i]) - (moyennesParClasse[i]*moyennesParClasse[i]);
  }
 }

 //calcul du critere
 cr=0.0;
 for (i=1;i<nbcl2;i++)
  cr += nbVoxelsParClasse[i]*variancesParClasse[i];

 cr = cr/(nbVoxelsCommuns*varianceTot);

  if (err!=NULL) *err=NO_ERR;

end_func:

 //liberation memoire
 if (nbVoxelsParClasse) FREE(nbVoxelsParClasse);
 if (moyennesParClasse) FREE(moyennesParClasse);
 if (variancesParClasse) FREE(variancesParClasse);

 return cr;

}

/*! @} */

/*! \ingroup DistanceFonction  @{ */
/*******************************************************************************
**     erreur_Woods_vraie_3d(im1,im2,champ,donnees_VP)
*/
/*!    l'erreur critere de Woods, la vraie... 
        \param im1, im2 : les 2 images
        \retval double : l'erreur dde rapport de correlation
*******************************************************************************/
double erreur_Woods_vraie_3d(ptr_distance_param dist_param)
{
 unsigned int i,j,k;
 TDimension wdth,hght,dpth;
 int nbcl1,nbcl2,max1,max2,cl1,cl2;
 int nbVoxelsCommuns=0;
 double *nbVoxelsParClasse=NULL, *moyennesParClasse=NULL, *variancesParClasse=NULL;
 double cr=DBL_MAX;
 grphic3d *im1=NULL, *im2=NULL;
 TYPEMRI3D ***im1MRI=NULL, ***im2MRI=NULL;
 ERREUR_RECALAGE *err=NULL;

 im1=dist_param->imreca; im2=dist_param->imref;
 im1MRI=im1->mri; im2MRI=im2->mri;
 err=dist_param->err;

 wdth=im1->width;hght=im1->height;dpth=im1->depth;
 if ((wdth!=im2->width)||(hght!=im2->height)||(dpth!=im2->depth))
 {
  RAISE_ERR_TEXT(err, PB_TAILLES_IMAGES, "tailles d'images differentes dans erreur_Woods_vraie_3d\n");
 }

 //calcul des maxpixel dans les deux images
 //et du nombre de classes de valeurs de voxels qu'on va utiliser
// nbcl1=imx_calculer_nombre_de_classes_3d_p(im1, &max1);
// nbcl2=imx_calculer_nombre_de_classes_3d_p(im2, &max2);

 nbcl1=dist_param->nbclReca; nbcl2=dist_param->nbclRef;
 max1=dist_param->maxReca; max2=dist_param->maxRef;

 if ((nbcl1<2)||(nbcl2<2))
 {
  RAISE_ERR_TEXT(err, NB_BINS_NUL, "nb de classes nulles dans erreur_Woods_vraie_3d\n");
 }

 //allocations memoires et initialisation
 nbVoxelsParClasse = CALLOC(nbcl2,double);
 moyennesParClasse = CALLOC(nbcl2,double);
 variancesParClasse = CALLOC(nbcl2,double);

 if ((!nbVoxelsParClasse)||(!moyennesParClasse)||(!variancesParClasse))
 {
  RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans erreur_Woods_vraie_3d\n");
 }

 //calcul du rapport de correlation
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    //on fait une classe a part pour les intensites nulles qui ne sont porteuses d'aucune information
    if (im1->mri[i][j][k]<=0) cl1=0;
    //else cl1=((nbcl1-1)*im1->mri[i][j][k])/(max1+1)+1;
    else cl1=im1->mri[i][j][k];
    if (im2->mri[i][j][k]<=0) cl2=0;
    else cl2=MINI(((nbcl2-1)*im2->mri[i][j][k])/(max2+1)+1, (nbcl2-1));

    //on ne fait les calculs que sur la zone de recouvrement des images
    if (cl1>0 && cl2>0)
    {
     //on rajoute 1.0 a la valeur de cl1 :
     //ne change pas l'info de l'image et permet de contourner le pb des moyennes nulles
     nbVoxelsParClasse[cl2] += 1.0;
     moyennesParClasse[cl2] += ((double)cl1+1.0);
     variancesParClasse[cl2] += ((double)cl1+1.0)*((double)cl1+1.0);
     nbVoxelsCommuns++;
    }

   }

 //si il n'y a pas de recouvrement des 2 images on renvoie la valeur max
 if (nbVoxelsCommuns==0) return 1.0;

 //on ne fait les calculs que sur la zone de recouvrement
 for (i=1;i<(unsigned int)nbcl2;i++)
 {
  if (nbVoxelsParClasse[i]>0)
  {
   moyennesParClasse[i] = moyennesParClasse[i]/nbVoxelsParClasse[i];
   variancesParClasse[i] = (variancesParClasse[i]/nbVoxelsParClasse[i]) - (moyennesParClasse[i]*moyennesParClasse[i]);
  }
 }

 //calcul du critere
 cr=0.0;
 for (i=1;i<(unsigned int)nbcl2;i++)
  if (nbVoxelsParClasse[i]>0)
   cr += nbVoxelsParClasse[i]*sqrt(variancesParClasse[i]);
   //cr += nbVoxelsParClasse[i]*sqrt(variancesParClasse[i])/moyennesParClasse[i];

 cr = cr/(double)nbVoxelsCommuns;

 if (err!=NULL) *err=NO_ERR;

end_func:

 //liberation memoire
 if (nbVoxelsParClasse) FREE(nbVoxelsParClasse);
 if (moyennesParClasse) FREE(moyennesParClasse);
 if (variancesParClasse) FREE(variancesParClasse);

 return cr;

}

/*! @} */

/*! \ingroup DistanceFonction  @{ */
/*******************************************************************************
**     erreur_woods_robust2_3d(im1,im2,champ,donnees_VP,points_calcules,weights)
*/
/*!    le critere de Woods modifiee par estimation robuste
        \param im1, im2 : les 2 images
        \retval double : l'erreur 
*******************************************************************************/
double erreur_woods_robust2_3d(ptr_distance_param dist_param)
{
 unsigned int i,j,k;
 TDimension wdth,hght,dpth;
 int nbcl1=0,nbcl2=0,max1,max2,cl1,cl2;
 int nbVoxelsCommuns=0;
 double *moyennesParClasse=NULL, *variancesParClasse=NULL;
 int *nbVoxelsParClasse=NULL;
 double cr=DBL_MAX;
 double **region=NULL;
 double diff, w;
 double **weights_par_region=NULL;
 double *region_weights=NULL;
 double nbVoxelsTot=0.0;
 int nb_voxels_courant=0;
 TYPEMRI3D val1=0, val2=0;
 grphic3d *im1=NULL, *im2=NULL;
 TYPEMRI3D ***im1MRI=NULL, ***im2MRI=NULL;
 double ***weights=NULL;
 ERREUR_RECALAGE *err=NULL;

 im1=dist_param->imreca; im2=dist_param->imref;
 im1MRI=im1->mri; im2MRI=im2->mri;
 err=dist_param->err;
 weights=dist_param->weights;

 wdth=im1->width;hght=im1->height;dpth=im1->depth;
 if ((wdth!=im2->width)||(hght!=im2->height)||(dpth!=im2->depth))
 {
  RAISE_ERR_TEXT(err, PB_TAILLES_IMAGES, "tailles d'images differentes dans erreur_woods_robust2_3d\n");
 }

 //calcul des maxpixel dans les deux images
 //et du nombre de classes de valeurs de voxels qu'on va utiliser
// nbcl1=imx_calculer_nombre_de_classes_3d_p(im1, &max1);
// nbcl2=imx_calculer_nombre_de_classes_3d_p(im2, &max2);

 nbcl1=dist_param->nbclReca; nbcl2=dist_param->nbclRef;
 max1=dist_param->maxReca; max2=dist_param->maxRef;

 if ((nbcl1<2)||(nbcl2<2))
 {
  RAISE_ERR_TEXT(err, NB_BINS_NUL, "nb de classes nulles dans erreur_woods_robust2_3d\n");
 }

 region=CALLOC(nbcl2,double *);
 for (i=0;i<(unsigned int)nbcl2;i++) region[i]=alloc_dvector(1);

 if (weights!=NULL)
 {
  weights_par_region=CALLOC(nbcl2,double *);
  if (!weights_par_region) RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans erreur_woods_robust2_3d\n");
  for (i=0;i<(unsigned int)nbcl2;i++)
  {
   weights_par_region[i]=alloc_dvector(1);
   if (!weights_par_region[i]) RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans erreur_woods_robust2_3d\n");
  }
  region_weights=CALLOC(nbcl2,double);
  if (!region_weights) RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans erreur_woods_robust2_3d\n");
 }

 //allocations memoires et initialisation
 nbVoxelsParClasse = CALLOC(nbcl2,int);
 moyennesParClasse = CALLOC(nbcl2,double);
 variancesParClasse = CALLOC(nbcl2,double);

 if ((!nbVoxelsParClasse)||(!moyennesParClasse)||(!variancesParClasse))
 {
  RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans erreur_woods_robust2_3d\n");
 }

 //calcul du rapport de correlation
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    //on fait une classe a part pour les intensites nulles qui ne sont porteuses d'aucune information
    val1=im1MRI[i][j][k]; val2=im2MRI[i][j][k];
    if (val2<=0) cl2=0;
    else cl2=MINI(((nbcl2-1)*val2)/(max2+1)+1, (nbcl2-1));
    if (val1<=0) cl1=0;
    else cl1=val1;
    if ((cl1>0)&&(cl2>0))
    {
     nb_voxels_courant=nbVoxelsParClasse[cl2];
     region[cl2][nb_voxels_courant]=(double) cl1;

     if (weights_par_region) weights_par_region[cl2][nb_voxels_courant]=(double) weights[i][j][k];

     nbVoxelsParClasse[cl2]++;
     nbVoxelsCommuns++;
     region[cl2]=realloc_dvector(region[cl2],(nb_voxels_courant+1));

     if (weights_par_region) weights_par_region[cl2]=realloc_dvector(weights_par_region[cl2],(nb_voxels_courant+1));
    }
   }

 //si il n'y a pas de recouvrement des 2 images on renvoie la valeur max
 if (nbVoxelsCommuns==0) return 1.0;

 //calcul des moyennes et des variances par classe
 //on ne fait le calcul que sur la zone de recouvrement des images

 //si il n'y a pas de poids, on calcule une variance classique
 if (weights==NULL)
 {
  for (i=1;i<(unsigned int)nbcl2;i++)
  {
   if (nbVoxelsParClasse[i]>0)
   {
    moyennesParClasse[i]=quick_select(nbVoxelsParClasse[i]/2,region[i],nbVoxelsParClasse[i]);
    for (j=0;j<(unsigned int)nbVoxelsParClasse[i];j++)
    {
     diff=region[i][j]-moyennesParClasse[i];
     w=diff*diff;
     variancesParClasse[i]+=w;
    }
    nbVoxelsTot+=nbVoxelsParClasse[i];
   }
  }

  //normalisation des variances
  for (i=1;i<(unsigned int)nbcl2;i++)
  {
   if (nbVoxelsParClasse[i]>0) variancesParClasse[i]=variancesParClasse[i]/(double)nbVoxelsParClasse[i];
  }
 }
 //si il y a de poids, on calcule une variance ponderee
 else
 {
  //calcul des poids totaux associes a chaque region
  for (i=1;i<(unsigned int)nbcl2;i++)
  {
   for (j=0;j<(unsigned int)nbVoxelsParClasse[i];j++)
   {
    region_weights[i]+=weights_par_region[i][j];
   }
  }

  //calcul des variances non normalisees par region
  for (i=1;i<(unsigned int)nbcl2;i++)
  {
   if (nbVoxelsParClasse[i]>0)
   {
    moyennesParClasse[i]=quick_select(nbVoxelsParClasse[i]/2,region[i],nbVoxelsParClasse[i]);
    for (j=0;j<(unsigned int)nbVoxelsParClasse[i];j++)
    {
     diff=region[i][j]-moyennesParClasse[i];
     w=weights_par_region[i][j]*diff*diff;
     variancesParClasse[i]+=w;
    }
    nbVoxelsTot+=nbVoxelsParClasse[i];
   }
  }

  //normalisation des variances par region
  for (i=1;i<(unsigned int)nbcl2;i++)
  {
   if (region_weights[i]>0) variancesParClasse[i]=variancesParClasse[i]/region_weights[i];
  }
 }

 //calcul du critere
 cr=0.0;
 for (i=1;i<(unsigned int)nbcl2;i++)
  if (nbVoxelsParClasse[i]>0)
   cr += (double)nbVoxelsParClasse[i]*sqrt(variancesParClasse[i]);

 cr=cr/nbVoxelsTot;

  if (err!=NULL) *err=NO_ERR;

end_func:

 //liberation memoire
 if (weights_par_region!=NULL)
 {
  for (i=0;i<(unsigned int)nbcl2;i++) { if (weights_par_region[i]) free_dvector(weights_par_region[i],(nbVoxelsParClasse[i]+1)); }
  if (weights_par_region) FREE(weights_par_region);
  if (region_weights) FREE(region_weights);
 }
 for (i=0;i<(unsigned int)nbcl2;i++) { if (region[i]) free_dvector(region[i],(nbVoxelsParClasse[i]+1)); }
 if (region) FREE(region);
 if (nbVoxelsParClasse) FREE(nbVoxelsParClasse);
 if (moyennesParClasse) FREE(moyennesParClasse);
 if (variancesParClasse) FREE(variancesParClasse);

 return cr;
}

/*! @} */

/*! \ingroup DistanceFonction  @{ */
/*******************************************************************************
**     erreur_quad_robust2_3d(im1,im2,champ,donnees_VP,points_calcules,weights)
*/
/*!    l'erreur quadratique ponderee pour estimation robuste
        \param im1, im2 : les 2 images
        \retval double : l'erreur
*******************************************************************************/
double erreur_quad_robust2_3d(ptr_distance_param dist_param)
{
 unsigned int i,j,k;
 TDimension wdth,hght,dpth,nb;
 double r=DBL_MAX,diff;
 double wTot=0.0;
 TYPEMRI3D m1=0,m2=0;
 grphic3d *im1=NULL, *im2=NULL;
 TYPEMRI3D ***im1MRI=NULL, ***im2MRI=NULL;
 double ***weights=NULL;
 ERREUR_RECALAGE *err=NULL;

 im1=dist_param->imreca; im2=dist_param->imref;
 im1MRI=im1->mri; im2MRI=im2->mri;
 err=dist_param->err;
 weights=dist_param->weights;

 wdth=im1->width;hght=im1->height;dpth=im1->depth;
 if ((wdth!=im2->width)||(hght!=im2->height)||(dpth!=im2->depth))
 {
  RAISE_ERR_TEXT(err, PB_TAILLES_IMAGES, "tailles d'images differentes dans erreur_quad_robust2_3d\n");
 }

 r=0.0;nb=0;
 //si il n'y a pas de poids, on est dans le cadre des moindres carres classique
 if (weights==NULL)
 {
  for (i=0;i<wdth;i++)
   for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
    {
     m1=im1MRI[i][j][k];m2=im2MRI[i][j][k];
     if ((m1>0)&&(m2>0))
     {
      diff=(double)(m1-m2);
      r+=diff*diff;
      wTot+=1.0;
     }
    }
 }
 //si il y a des poids, on obtient un critere des moindres carres ponderes
 else
 {
  for (i=0;i<wdth;i++)
   for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
    {
     m1=im1MRI[i][j][k];m2=im2MRI[i][j][k];
     if ((m1>0)&&(m2>0))
     {
      diff=(double)(m1-m2);
      r+=diff*diff*weights[i][j][k];
      wTot+=weights[i][j][k];
     }
    }
 }
 r=sqrt(r/wTot);

  if (err!=NULL) *err=NO_ERR;

end_func:

 return(r);
}

/*! @} */

/*! \ingroup DistanceFonction  @{ */
/*******************************************************************************
**     erreur_CR_robust_3d(im1,im2,champ,donnees_VP,points_calcules,weights)
*/
/*!    le critere de Woods modifiee par estimation robuste
        \param im1, im2 : les 2 images
        \retval double : l'erreur
*******************************************************************************/
double erreur_CR_robust_3d(ptr_distance_param dist_param)
{
 unsigned int i,j,k;
 TDimension wdth,hght,dpth;
 int nbcl1=0,nbcl2=0,max1,max2,cl1,cl2;
 int nbVoxelsCommuns=0;
 double *moyennesParClasse=NULL, *variancesParClasse=NULL;
 int *nbVoxelsParClasse=NULL;
 double cr=DBL_MAX;
 double **region=NULL;
 double diff, w;
 double **weights_par_region=NULL;
 double moyTot=0.0, varTot=0.0;
 double *region_weights=NULL;
 grphic3d *im1=NULL, *im2=NULL;
 TYPEMRI3D ***im1MRI=NULL, ***im2MRI=NULL;
 double ***weights=NULL;
 ERREUR_RECALAGE *err=NULL;

 im1=dist_param->imreca; im2=dist_param->imref;
 im1MRI=im1->mri; im2MRI=im2->mri;
 err=dist_param->err;
 weights=dist_param->weights;

 wdth=im1->width;hght=im1->height;dpth=im1->depth;
 if ((wdth!=im2->width)||(hght!=im2->height)||(dpth!=im2->depth))
 {
  RAISE_ERR_TEXT(err, PB_TAILLES_IMAGES, "tailles d'images differentes dans erreur_CR_robust_3d\n");
 }

 //calcul des maxpixel dans les deux images
 //et du nombre de classes de valeurs de voxels qu'on va utiliser
// nbcl1=imx_calculer_nombre_de_classes_3d_p(im1, &max1);
// nbcl2=imx_calculer_nombre_de_classes_3d_p(im2, &max2);
 nbcl1=dist_param->nbclReca; nbcl2=dist_param->nbclRef;
 max1=dist_param->maxReca; max2=dist_param->maxRef;

 if ((nbcl1<2)||(nbcl2<2))
 {
  RAISE_ERR_TEXT(err, NB_BINS_NUL, "nb de classes nulles dans erreur_CR_robust_3d\n");
 }

 region=CALLOC(nbcl2,double *);
 for (i=0;i<(unsigned int)nbcl2;i++) region[i]=alloc_dvector(1);

 //allocation du tableau des poids affectes aux differents voxels dans chaque region
 if (weights!=NULL)
 {
  weights_par_region=CALLOC(nbcl2,double *);
  if (!weights_par_region) RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans erreur_CR_robust_3d\n");
  for (i=0;i<(unsigned int)nbcl2;i++)
  {
   weights_par_region[i]=alloc_dvector(1);
   if (!weights_par_region[i]) RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans erreur_CR_robust_3d\n");
  }
  region_weights=CALLOC(nbcl2,double);
  if (!region_weights) RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans erreur_CR_robust_3d\n");
 }

 //allocations memoires et initialisation
 nbVoxelsParClasse = CALLOC(nbcl2,int);
 moyennesParClasse = CALLOC(nbcl2,double);
 variancesParClasse = CALLOC(nbcl2,double);
 if ((!nbVoxelsParClasse)||(!moyennesParClasse)||(!variancesParClasse))
 {
  RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans erreur_CR_robust_3d\n");
 }

 //calcul du rapport de correlation
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    //on fait une classe a part pour les intensites nulles qui ne sont porteuses d'aucune information
    if (im2MRI[i][j][k]<=0) cl2=0;
    else cl2=MINI(((nbcl2-1)*im2MRI[i][j][k])/(max2+1)+1, (nbcl2-1));
    if (im1MRI[i][j][k]<=0) cl1=0;
    else cl1=im1MRI[i][j][k];
    if ((cl1>0)&&(cl2>0))
    {
     region[cl2][nbVoxelsParClasse[cl2]]=(double) cl1;

     if (weights_par_region) weights_par_region[cl2][nbVoxelsParClasse[cl2]]=(double) weights[i][j][k];

     nbVoxelsParClasse[cl2]++;
     nbVoxelsCommuns++;
     region[cl2]=realloc_dvector(region[cl2],(nbVoxelsParClasse[cl2]+1));

     if (weights_par_region) weights_par_region[cl2]=realloc_dvector(weights_par_region[cl2],(nbVoxelsParClasse[cl2]+1));

     moyTot+=(double)cl1;
     varTot+=(double)(cl1*cl1);
    }
   }

 //si il n'y a pas de recouvrement des 2 images on renvoie la valeur max
 if (nbVoxelsCommuns==0) return 1.0;

 //on ne fait le calcul que sur la zone de recouvrement des images

 //calcul des variances non normalisees par classe
 //si il n'y a pas de poids, on calcule une variance classique
 if (weights==NULL)
 {
  for (i=1;i<(unsigned int)nbcl2;i++)
  {
   if (nbVoxelsParClasse[i]>0)
   {
    moyennesParClasse[i]=quick_select(nbVoxelsParClasse[i]/2,region[i],nbVoxelsParClasse[i]);
    for (j=0;j<(unsigned int)nbVoxelsParClasse[i];j++)
    {
     diff=region[i][j]-moyennesParClasse[i];
     w=diff*diff;
     variancesParClasse[i]+=w;
    }
   }
  }

  //normalisation des variances
  for (i=1;i<(unsigned int)nbcl2;i++)
  {
   if (nbVoxelsParClasse[i]>0) variancesParClasse[i]=variancesParClasse[i]/(double)nbVoxelsParClasse[i];
   else variancesParClasse[i]=0.0;
  }
 }
 //si il y a de poids, on calcule une variance ponderee
 else
 {

  //calcul des variances non normalisees par region
  for (i=1;i<(unsigned int)nbcl2;i++)
  {
   if (nbVoxelsParClasse[i]>0)
   {
    moyennesParClasse[i]=quick_select(nbVoxelsParClasse[i]/2,region[i],nbVoxelsParClasse[i]);
    for (j=0;j<(unsigned int)nbVoxelsParClasse[i];j++)
    {
     //calcul des poids totaux associes a chaque region
     region_weights[i]+=weights_par_region[i][j];
     diff=region[i][j]-moyennesParClasse[i];
     w=weights_par_region[i][j]*diff*diff;
     variancesParClasse[i]+=w;
    }
   }
  }

  //normalisation des variances par region
  for (i=1;i<(unsigned int)nbcl2;i++)
  {
   if (region_weights[i]>0) variancesParClasse[i]=variancesParClasse[i]/region_weights[i];
   else variancesParClasse[i]=0;
//   else printf ("i: %d, variancesParClasse: %f, region_weights: %f\n", i, variancesParClasse[i], region_weights[i]);
  }
 }

 //calcul de la variance totale
 moyTot = moyTot/nbVoxelsCommuns;
 varTot = (varTot/nbVoxelsCommuns) - (moyTot*moyTot);

 //calcul du critere
 cr=0.0;
 for (i=1;i<(unsigned int)nbcl2;i++)
  if (nbVoxelsParClasse[i]>0)
   cr += (double)nbVoxelsParClasse[i]*variancesParClasse[i];

 cr=cr/((double)nbVoxelsCommuns*varTot);

  if (err!=NULL) *err=NO_ERR;

end_func:

 //liberation memoire
 if (weights_par_region!=NULL)
 {
  for (i=0;i<(unsigned int)nbcl2;i++) {if (weights_par_region[i]) free_dvector(weights_par_region[i],(nbVoxelsParClasse[i]+1)); }
  if (weights_par_region) FREE(weights_par_region);
  if (region_weights) FREE(region_weights);
 }
 for (i=0;i<(unsigned int)nbcl2;i++) { if (region[i]) free_dvector(region[i],(nbVoxelsParClasse[i]+1)); }
 if (region) FREE(region);
 if (nbVoxelsParClasse) FREE(nbVoxelsParClasse);
 if (moyennesParClasse) FREE(moyennesParClasse);
 if (variancesParClasse) FREE(variancesParClasse);

 return cr;
}

/*! @} */

/*******************************************************************************
********** DISTANCES RELATIVES A L'INTERPOLATION VOLUME PARTIEL ****************
*******************************************************************************/

/*! \ingroup DistanceFonction  @{ */
/*******************************************************************************
**     erreur_entropie_conjointe_VP_3d(Imreca, Imref, champ, donnees_VP)
*/
/*!    erreur entropie conjointe volume partiel
        \param Imreca, Imref : les 2 images
        \param champ : le champ de deformation
        \param donnees_VP : structure contenant les donnees pour le volume partiel
        \retval double : l'erreur d'information mut.
*******************************************************************************/
double erreur_entropie_conjointe_VP_3d(ptr_distance_param dist_param)
{
 double entropie_conjointe=0.0;

 entropie_conjointe = calcul_info_mutuelle_VP_3d(ENTROPIE_CONJOINTE, dist_param);

 return entropie_conjointe;
}

/*! @} */

/*! \ingroup DistanceFonction  @{ */
/*******************************************************************************
**     erreur_IM_VP_3d(Imreca, Imref, champ, donnees_VP)
*/
/*!    erreur information mutuelle volume partiel
        \param Imreca, Imref : les 2 images
        \param champ : le champ de deformation
        \param donnees_VP : structure contenant les donnees pour le volume partiel
        \retval double : l'erreur d'information mut.
*******************************************************************************/
double erreur_IM_VP_3d(ptr_distance_param dist_param)
{
 double info_mutuelle=0.0;

 info_mutuelle = calcul_info_mutuelle_VP_3d(IM, dist_param);

 return info_mutuelle;
}

/*! @} */

/*! \ingroup DistanceFonction  @{ */
/*******************************************************************************
**     erreur_IMNS_VP_3d(Imreca, Imref, champ, donnees_VP)
*/
/*!    erreur information mutuelle normalisee de studholme volume partiel
        \param Imreca, Imref : les 2 images
        \param champ : le champ de deformation
        \param donnees_VP : structure contenant les donnees pour le volume partiel
        \retval double : l'erreur d'information mut.
*******************************************************************************/
double erreur_IMNS_VP_3d(ptr_distance_param dist_param)
{
 double info_mutuelle_normalisee=0.0;

 info_mutuelle_normalisee = calcul_info_mutuelle_VP_3d(IMNS, dist_param);

 return info_mutuelle_normalisee;
}

/*! @} */

/*! \ingroup DistanceFonction  @{ */
/*******************************************************************************
**     erreur_IMNM1_VP_3d(Imreca, Imref, champ, donnees_VP)
*/
/*!    erreur information mutuelle normalisee de Maes 1 volume partiel
        \param Imreca, Imref : les 2 images
        \param champ : le champ de deformation
        \param donnees_VP : structure contenant les donnees pour le volume partiel
        \retval double : l'erreur d'information mut.
*******************************************************************************/
double erreur_IMNM1_VP_3d(ptr_distance_param dist_param)
{
 double info_mutuelle_normalisee=0.0;

 info_mutuelle_normalisee = calcul_info_mutuelle_VP_3d(IMNM1, dist_param);

 return info_mutuelle_normalisee;
}

/*! @} */

/*! \ingroup DistanceFonction  @{ */
/*******************************************************************************
**     erreur_IMNM2_VP_3d(Imreca, Imref, champ, donnees_VP)
*/
/*!    erreur information mutuelle normalisee de Maes 2 volume partiel
        \param Imreca, Imref : les 2 images
        \param champ : le champ de deformation
        \param donnees_VP : structure contenant les donnees pour le volume partiel
        \retval double : l'erreur d'information mut.
*******************************************************************************/
double erreur_IMNM2_VP_3d(ptr_distance_param dist_param)
{
 double info_mutuelle_normalisee=0.0;

 info_mutuelle_normalisee = calcul_info_mutuelle_VP_3d(IMNM2, dist_param);

 return info_mutuelle_normalisee;
}

/*! @} */

/*! \ingroup DistanceFonction  @{ */
/*******************************************************************************
**        erreur_IM_VP_stochastique_3d(Imreca, Imref, points_choisis, champ)
**
**        calcul de l'information mutuelle, methode stochastique
**
**        Imreca, Imref: images a comparer
**        points_choisis: points choisis aleatoirement
**        champ: champs de la transformation appliquee
********************************************************************************/
double erreur_IM_VP_stochastique_3d(ptr_distance_param dist_param)
{
 double info_mutuelle=0.0;

// info_mutuelle = calcul_info_mutuelle_VP_stochastique_3d(IM, dist_param->imreca, dist_param->imref, dist_param->champ, dist_param->donnees_VP, dist_param->points_calcules, dist_param->err);
 info_mutuelle = calcul_info_mutuelle_VP_stochastique_3d(IM, dist_param);

 return info_mutuelle;
}

/*! @} */

/*! \ingroup DistanceFonction  @{ */
/*******************************************************************************
**        erreur_IM_BSpline_3d(Imreca, Imref, points_choisis, champ)
**
**        calcul de l'information mutuelle, methode de construction de l'histogramme par des BSplines
**
**        Imreca, Imref: images a comparer
**        points_choisis: points choisis aleatoirement
**        champ: champs de la transformation appliquee
********************************************************************************/
double erreur_IM_BSpline_3d(ptr_distance_param dist_param)
{
 double info_mutuelle=0.0;

// info_mutuelle=calcul_info_mutuelle_BSpline_3d(IM, dist_param->imreca, dist_param->imref, dist_param->champ, dist_param->donnees_VP, dist_param->err);
 info_mutuelle=calcul_info_mutuelle_BSpline_3d(IM, dist_param);

 return info_mutuelle;
}
/*! @} */

/*! \ingroup DistanceFonction  @{ */
/*******************************************************************************
**        erreur_IM_Region_3d(Imreca, Imref, points_choisis, champ)
**
**        calcul de l'information mutuelle, methode combinant des informations statistiques et
**        spatiales
**
**        Imreca, Imref: images a comparer
**        points_choisis: points choisis aleatoirement
**        champ: champs de la transformation appliquee
********************************************************************************/
double   erreur_IM_Region_3d(ptr_distance_param dist_param)
{
 double info_mutuelle=0.0;

 info_mutuelle=calcul_info_mutuelle_region_3d(IM, dist_param);

 return info_mutuelle;
}
/*! @} */

/*! \ingroup DistanceFonction  @{ */
/*******************************************************************************
**        erreur_coeff_correlation_d(Imreca, Imref, points_choisis, champ)
**
**        calcul du coefficient de correlation
**
**        Imreca, Imref: images a comparer
**        points_choisis: points choisis aleatoirement
**        champ: champs de la transformation appliquee
********************************************************************************/
double erreur_coeff_correlation_d(ptr_distance_param dist_param)
{
 int i,j,k;
 int wdth,hght,dpth;
 double moy1=0.0, moy2=0.0, var1=0.0, var2=0.0, corr=0.0;
 grphic3d *im1=NULL, *im2=NULL;
 TYPEMRI3D ***im1MRI=NULL, ***im2MRI=NULL;
 TYPEMRI3D val1, val2;
 ERREUR_RECALAGE *err=NULL;
 double nb_overlap=0.0;

 im1=dist_param->imreca; im2=dist_param->imref;
 im1MRI=im1->mri; im2MRI=im2->mri;
 err=dist_param->err;

 wdth=(int)im1->width;hght=(int)im1->height;dpth=(int)im1->depth;
 if ((wdth!=(int)im2->width)||(hght!=(int)im2->height)||(dpth!=(int)im2->depth))
 {
  RAISE_ERR_TEXT(err, PB_TAILLES_IMAGES, "tailles d'images differentes dans erreur_coeff_correlation_d\n");
 }

 //calcul des moyennes et zone d'overlap
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    val1=im1MRI[i][j][k]; val2=im2MRI[i][j][k];
    if ((val1>0)&&(val2>0))
    {
     nb_overlap+=1.0;
     moy1+=(double)val1;
     moy2+=(double)val2;
    }
   }

 if (nb_overlap==0.0) return 0.0;
 moy1/=nb_overlap; moy2/=nb_overlap;

 //calcul des variances et du coefficient de correlation
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    val1=im1MRI[i][j][k]; val2=im2MRI[i][j][k];
    if ((val1>0)&&(val2>0))
    {
     corr+=((double)val1-moy1)*((double)val2-moy2);
     var1+=((double)(val1*val1))/nb_overlap;
     var2+=((double)(val2*val2))/nb_overlap;
    }
   }

  if ((var1<=0.0)||(var2<=0.0)) return 0.0;
  var1-=(moy1*moy1); var2-=(moy2*moy2);
  corr/=(nb_overlap*var1*var2);

  if (err!=NULL) *err=NO_ERR;

end_func:

 //on retourne -corr car on doit maximiser le coefficient
 return (-corr);
}
/*! @} */

/*******************************************************************************
**     calcul_info_mutuelle_3d(im_type, im1, im2, champ, donnees_VP)
*/
/*!    calcul de la distance se basant sur l'information mutuelle simple
    \param im_type: le type de distance utilisee
        \param im1, im2 : les 2 images
        \param champ : le champ de deformation
        \param donnees_VP : structure contenant les donnees pour le volume partiel
        \retval double : l'erreur d'information mut.
*******************************************************************************/
double calcul_info_mutuelle_3d(IM3DTYPE im_type, ptr_distance_param dist_param)
{
 int i,j;
 int nbcl1,nbcl2,max1,max2;
 double **histo_conjoint=NULL,*histo1=NULL,*histo2=NULL;
 double info_mutuelle=DBL_MAX, entropie_conjointe=0.0, entropie1=0.0, entropie2=0.0;
 double nbVoxelsOverlap=0.0, nbVoxelsHorsOverlap=0.0;
 grphic3d *im1=dist_param->imreca;
 grphic3d *im2=dist_param->imref;
 ERREUR_RECALAGE *err=dist_param->err;

 //calcul des maxpixel dans les deux images
 //et du nombre de classes de valeurs de voxels qu'on va utiliser
// nbcl1=imx_calculer_nombre_de_classes_3d_p(im1, &max1);
// nbcl2=imx_calculer_nombre_de_classes_3d_p(im2, &max2);
 nbcl1=dist_param->nbclReca; nbcl2=dist_param->nbclRef;
 max1=dist_param->maxReca; max2=dist_param->maxRef;

 if ((nbcl1<2)||(nbcl2<2))
 {
  RAISE_ERR_TEXT(err, NB_BINS_NUL, "nb de classes nulles dans calcul_info_mutuelle_3d\n");
 }

 //allocation memoire et initialisation des histos
 histo1=CALLOC(nbcl1,double);histo2=CALLOC(nbcl2,double);
 histo_conjoint=alloc_dmatrix(nbcl1,nbcl2);
 if ((!histo1)||(!histo2)||(!histo_conjoint))
 {
  RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans calcul_info_mutuelle_3d\n");
 }
 for (i=0;i<nbcl1;i++)
  for (j=0;j<nbcl2;j++)
   histo_conjoint[i][j]=0.0;

 //appel a la fonction de calcul des histogrammes appropriee
 calcul_histo_conjoint_simple(im1, im2, histo_conjoint, nbcl1, nbcl2, max1, max2, err);
 if ((err)&&(*err!=NO_ERR)) goto end_func;

 //on calcule le nb de voxel dans la zone de recouvrement
 for (i=1;i<nbcl1;i++)
  for (j=1;j<nbcl2;j++)
   nbVoxelsOverlap += histo_conjoint[i][j];

 //si aucun recouvrement on renvoie le max de l'energie
 if (nbVoxelsOverlap==0.0)
 {
  info_mutuelle=max_erreur_info_mutuelle(im_type, nbcl1, nbcl2);
  RAISE_ERR_TEXT(err, OVERLAP_NUL, "zone d'overlap nulle dans calcul_info_mutuelle_3d\n");
 }
 
 nbVoxelsHorsOverlap = (double)(im1->width*im1->height*im1->depth)-nbVoxelsOverlap;

 info_mutuelle=calc_info_mutuelle_histos(im_type, histo_conjoint, histo1, histo2, nbcl1, nbcl2, nbVoxelsOverlap, nbVoxelsHorsOverlap, &entropie_conjointe, &entropie1, &entropie2);

 if (err) *err=NO_ERR;

end_func:

 //liberation memoire
 if (histo_conjoint) free_dmatrix(histo_conjoint,nbcl1,nbcl2);
 if (histo1) FREE(histo1);
 if (histo2) FREE(histo2);

 return info_mutuelle;
}

/*******************************************************************************
**     calcul_info_mutuelle_VP_3d(im_type, Imreca, Imref, champ, donnees_VP)
*/
/*!    calcul de la distance se basant sur l'information calculee avec le volume partiel
    \param im_type: le type de distance utilisee
        \param Imreca, Imref : les 2 images
        \param champ : le champ de deformation
        \param donnees_VP : structure contenant les donnees pour le volume partiel
        \retval double : l'erreur d'information mut.
*******************************************************************************/
//double calcul_info_mutuelle_VP_3d(IM3DTYPE im_type, grphic3d *Imreca, grphic3d *Imref, field3d *champ, VP_histos * donnees_VP, ERREUR_RECALAGE *err)
double calcul_info_mutuelle_VP_3d(IM3DTYPE im_type, ptr_distance_param dist_param)
{
 int nbcl1,nbcl2;
 int i,j;
 double **histo_conjoint=NULL,*histo1=NULL,*histo2=NULL;
 double info_mutuelle=DBL_MAX, entropie_conjointe=DBL_MAX, entropie1=DBL_MAX, entropie2=DBL_MAX;
 double nbVoxelsOverlap=0.0, nbVoxelsHorsOverlap=0.0;
 grphic3d *Imref=dist_param->imref;
 VP_histos * donnees_VP=dist_param->donnees_VP;
 ERREUR_RECALAGE *err=dist_param->err;

 if (donnees_VP==NULL)
 {
   RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "la structure devant contenant l'histogramme conjoint pour le volume partiel n'est pas allouee dans calcul_info_mutuelle_VP_3d\n");
 }
  
 histo_conjoint=donnees_VP->histo_conjoint;
 histo1=donnees_VP->histo_marginal1;
 histo2=donnees_VP->histo_marginal2;

 //appel a la fonction de calcul des histogrammes appropriee
// calcul_histogrammes_VP(Imreca, Imref, champ, donnees_VP, err);
 calcul_histogrammes_VP(dist_param);
 if ((err)&&(*err!=NO_ERR)) goto end_func;

// nbcl1 = donnees_VP->nb_cl1;
// nbcl2 = donnees_VP->nb_cl2;
 nbcl1=dist_param->nbclReca; nbcl2=dist_param->nbclRef;
// max1=dist_param->maxReca; max2=dist_param->maxRef;
 
 //on calcule le nb de voxel dans la zone de recouvrement
 for (i=1;i<nbcl1;i++)
  for (j=1;j<nbcl2;j++)
   nbVoxelsOverlap += histo_conjoint[i][j];  

 //si aucun recouvrement on renvoie le max de l'energie
 if (nbVoxelsOverlap==0.0)
 {
  RAISE_ERR_TEXT(err, OVERLAP_NUL, "zone d'overlap nulle dans calcul_info_mutuelle_VP_3d\n");
 }

 nbVoxelsHorsOverlap = (double)(Imref->width*Imref->height*Imref->depth)-nbVoxelsOverlap;

 info_mutuelle=calc_info_mutuelle_histos(im_type, histo_conjoint, histo1, histo2, nbcl1, nbcl2, nbVoxelsOverlap, nbVoxelsHorsOverlap, &entropie_conjointe, &entropie1, &entropie2);

 if (err) *err=NO_ERR;

end_func:

  //update des donnees
 donnees_VP->entropie_conjointe = entropie_conjointe;
 donnees_VP->entropie1 = entropie1;
 donnees_VP->entropie2 = entropie2;

 return info_mutuelle;
}

/*******************************************************************************
**     calcul_IM(entropie_conjointe, entropie1, entropie2)
*/
/*!    calcul de l'information mutuelle telle que proposee par Collignonn et al
       a partir des entropies des images et de l'entropie conjointe
    \param entropie_conjointe: l'entropie conjointe
        \param entropie1, entropie2 : les entropies des images
        \retval double : l'information mutuelle
*******************************************************************************/
double calcul_IM(double entropie_conjointe, double entropie1, double entropie2)
{
 double info_mutuelle;
  
 info_mutuelle = entropie1 + entropie2 - entropie_conjointe;
  
 return info_mutuelle;
}

/*******************************************************************************
**     calcul_IM(entropie_conjointe, entropie1, entropie2)
*/
/*!    calcul de l'information mutuelle normalisee de Studholme
       a partir des entropies des images et de l'entropie conjointe
    \param entropie_conjointe: l'entropie conjointe
        \param entropie1, entropie2 : les entropies des images
        \retval double : l'information mutuelle
*******************************************************************************/
double calcul_IMNS(double entropie_conjointe, double entropie1, double entropie2)
{
 double info_mutuelle_normalisee;

 if (entropie_conjointe==0.0) return 0.0;
 
 info_mutuelle_normalisee = (entropie1 + entropie2)/entropie_conjointe;
 
 return info_mutuelle_normalisee; 
}

/*******************************************************************************
**     calcul_IM(entropie_conjointe, entropie1, entropie2)
*/
/*!    calcul de l'information mutuelle normalisee de Maes 1
       a partir des entropies des images et de l'entropie conjointe
    \param entropie_conjointe: l'entropie conjointe
        \param entropie1, entropie2 : les entropies des images
        \retval double : l'information mutuelle
*******************************************************************************/
double calcul_IMNM1(double entropie_conjointe, double entropie1, double entropie2)
{
 double info_mutuelle_normalisee;

 if ((entropie1 + entropie2)==0.0) return 0.0;
 
 info_mutuelle_normalisee = (entropie1 + entropie2 - entropie_conjointe)/(entropie1 + entropie2);

 return info_mutuelle_normalisee; 
}

/*******************************************************************************
**     calcul_IM(entropie_conjointe, entropie1, entropie2)
*/
/*!    calcul de l'information mutuelle normalisee de Maes 2
       a partir des entropies des images et de l'entropie conjointe
    \param entropie_conjointe: l'entropie conjointe
        \param entropie1, entropie2 : les entropies des images
        \retval double : l'information mutuelle
*******************************************************************************/  
double calcul_IMNM2(double entropie_conjointe, double entropie1, double entropie2)
{
 double info_mutuelle_normalisee;

 info_mutuelle_normalisee = 2.0*entropie_conjointe - entropie1 - entropie2;

 return info_mutuelle_normalisee; 
}

/*******************************************************************************
**    calcul_histo_conjoint_simple(im1,im2,histo_conjoint,nbcl1,nbcl2,max1,max2)
*/
/*!   calcul de l'histogramme conjoint en comptant le nombre de 
      pixels dans chaque case et divisant par le nombre total de pixels
        \param Imreca, im2 : les 2 images
        \param histo_joint : le tableau de l'histo alloue
        \param nbcl1, nbcl2 : le nombre de classes d'histo pour chaque image
        \param max1, max2 : les valeurs max des pixels dans chaque image
*******************************************************************************/
void calcul_histo_conjoint_simple(grphic3d *im1, grphic3d *im2, double **histo_conjoint, int nbcl1, int nbcl2, int max1, int max2, ERREUR_RECALAGE *err)
{
 int i,j,k;
 int wdth,hght,dpth;
 int cl1,cl2;
 TYPEMRI3D ***im1Mri=NULL, ***im2MRI=NULL;
 
 //calcul des dimensions des images
 //les 2 images sont censees avoir les meme dimensions
 wdth=(int)im1->width;hght=(int)im1->height;dpth=(int)im1->depth;
 if ((wdth!=(int)im2->width)||(hght!=(int)im2->height)||(dpth!=(int)im2->depth))
 {
  RAISE_ERR_TEXT(err, PB_TAILLES_IMAGES, "tailles d'images differentes dans calcul_histo_conjoint_simple\n");
 }

 //calcul de l'histogramme conjoint
 im1Mri=im1->mri; im2MRI=im2->mri;
 nbcl1=nbcl1-1; nbcl2=nbcl2-1;
 max1=max1+1; max2=max2+1;
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    if (im1Mri[i][j][k]<=0) cl1=0;
    else cl1=(nbcl1*im1Mri[i][j][k])/max1+1;
    if (im2MRI[i][j][k]<=0) cl2=0;
    else cl2=(nbcl2*im2MRI[i][j][k])/max2+1;
    cl1=MINI(cl1,nbcl1);
    cl2=MINI(cl2,nbcl2);
    histo_conjoint[cl1][cl2]+=1.0; 
   }

  if (err!=NULL) *err=NO_ERR;

end_func:
  
  return;
}

void calcul_histo_conjoint_lineaire(grphic3d *im1, grphic3d *im2, double **histo_conjoint, int nbcl1, int nbcl2, int max1, int max2, ERREUR_RECALAGE *err)
{
 unsigned int i,j,k;
 TDimension wdth,hght,dpth;
 int cl11, cl12, cl21, cl22;
 double cl1, cl2, dcl1, dcl2;

 //calcul des dimensions des images
 //les 2 images sont censees avoir les meme dimensions
 wdth=im1->width;hght=im1->height;dpth=im1->depth;
 if ((wdth!=im2->width)||(hght!=im2->height)||(dpth!=im2->depth))
 {
  RAISE_ERR_TEXT(err, PB_TAILLES_IMAGES, "tailles d'images differentes dans calcul_histo_conjoint_lineaire\n");
 }

 //calcul de l'histogramme conjoint
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    cl1 = (double)((nbcl1-1)*im1->mri[i][j][k])/(double)(max1+1);
    cl2 = (double)((nbcl2-1)*im2->mri[i][j][k])/(double)(max2+1);
    cl11 = (int)floor(cl1); cl12 = (int)floor(cl1)+1; cl21 = (int)floor(cl2); cl22 = (int)floor(cl2)+1;
    dcl1 = (int)cl12-cl1; dcl2=(double)cl22-cl2;
    histo_conjoint[cl11][cl21]+=  dcl1*dcl2;
    histo_conjoint[cl12][cl21]+=  (1.0-dcl1)*dcl2;
    histo_conjoint[cl11][cl22]+=  dcl1*(1.0-dcl2);
    histo_conjoint[cl12][cl22]+=  (1.0-dcl1)*(1.0-dcl2);   
   }

  if (err!=NULL) *err=NO_ERR;

end_func:

  return;
}

/*******************************************************************************
**    calcul_histos_marginaux(histo_conjoint,histo1,histo2,nbcl1,nbcl2)
*/
/*!   calcul de l'histogramme d'une image
        \param im : l'image
        \param histo : le tableau de l'histo alloue
        \param nbcl : le nombre de classes d'histo pour l'image
        \param max : la valeur max des pixels dans l'image
*******************************************************************************/
void calcul_histos_marginaux(double **histo_conjoint, double *histo1, double *histo2, int nbcl1, int nbcl2, double nbVoxelsHorsOverlap)
{
 int i,j;
 
 /*on met tous les voxels hors de la zone de recouvrement dans la classe nulle*/
 histo1[0]=nbVoxelsHorsOverlap;
 histo2[0]=nbVoxelsHorsOverlap;
  
 /*les histos marginaux sont deductibles de l'histo conjoint en sommant ses valeurs par ligne et par colonne 
 dans la zone de recouvrement*/
 for (i=1;i<nbcl1;i++)
 {
  histo1[i]=0.0;
  for (j=1;j<nbcl2;j++)
   histo1[i] += histo_conjoint[i][j];
 } 
 
 for (i=1;i<nbcl2;i++)
 {
  histo2[i]=0.0;
  for (j=1;j<nbcl1;j++)
   histo2[i] += histo_conjoint[j][i];
 } 
 
  return;
}

/*******************************************************************************
**    calcul_entropie_conjointe(histo_conjoint,nbcl1,nbcl2)
*/
/*!   calcul de l'entropie conjointe a partir de l'histo conjoint.
      On assume que l'histogramme est normalise par le nombre de pixels.
        \param histo_joint : l'histo conjoint
        \param nbcl1, nbcl2 : le nombre de classes d'histo pour chaque image
        \retval double : l'entropie conjointe calculee
*******************************************************************************/
double calcul_entropie_conjointe(double **histo_conjoint, int nbcl1, int nbcl2, double nbVoxelsOverlap)
{
 int i,j;
 double entropie_conjointe=0.0;
   
 //calcul de l'entropie
 for (i=1;i<nbcl1;i++) 
  for (j=1;j<nbcl2;j++)
   if (histo_conjoint[i][j]!=0.0)
   {
    if (histo_conjoint[i][j]>0.0) entropie_conjointe -= histo_conjoint[i][j]*log(histo_conjoint[i][j]);
    else fprintf(stderr,"valeur negative dans l'histo conjoint!!");
   }
   
 entropie_conjointe = entropie_conjointe/nbVoxelsOverlap + log(nbVoxelsOverlap);
 
 return entropie_conjointe;
}


/*******************************************************************************
**    calcul_entropie(histo,nbcl,nbPixels)
*/
/*!   calcul de l'entropie a partir de l'histogramme
        \param histo : l'histogramme
        \param nbcl : le nombre de classe d'histo pour l'image
        \retval double : l'entropie calculee
*******************************************************************************/
double calcul_entropie(double *histo, int nbcl, double nbVoxelsOverlap)
{
 int i;
 double entropie=0.0;
  
 //calcul de l'entropie
 for (i=1;i<nbcl;i++)
  if (histo[i]!=0.0) 
   entropie -= histo[i]*log(histo[i]);
   
 entropie = entropie/nbVoxelsOverlap + log(nbVoxelsOverlap);
 
 return entropie;
}

/*******************************************************************************
**    calcul_histogrammes_VP(Imreca,Imref,champ,donnees_VP)
*/
/*!
**    fonction de calcul de l'histogramme conjoint pour le volume partiel
**  \param imreca : image a recaler
**  \param champ : le champ
**  \param donnees_VP : structure ou mettre les donnees calculees
**  \retval int : 1 en cas de succes, 0 dans le cas contraire
*******************************************************************************/
int calcul_histogrammes_VP(ptr_distance_param dist_param)
{
 int wdth,hght,dpth;
 int wdthReca,hghtReca,dpthReca;
 int i,j,k,i2,j2,k2,i3,j3,k3;
 double x,y,z,xrel,yrel,zrel;
 double w[8]={0.0};
 TYPEMRI3D v[8]={0};
 int max1,max2,nbcl1,nbcl2;
 vector3d ***data=NULL;
 double ** histo_conjoint=NULL;
 int pixelRef;
 grphic3d *imrecaQuantif=NULL;
 grphic3d *Imreca=dist_param->imreca;
 grphic3d *Imref=dist_param->imref;
 field3d *champ=dist_param->champ;
 VP_histos * donnees_VP=dist_param->donnees_VP;
 ERREUR_RECALAGE *err=dist_param->err;

 int nbcl=0, max=0;

 wdthReca=Imreca->width;
 hghtReca=Imreca->height;
 dpthReca=Imreca->depth;
 
///////// Initialisations /////////
  if (donnees_VP==NULL) {
   RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "la structure devant contenant l'histogramme conjoint pour le volume partiel n'est pas allouee dans calcul_histogrammes_VP\n");
  }
  
  histo_conjoint = donnees_VP->histo_conjoint;

  wdth=Imref->width;hght=Imref->height;dpth=Imref->depth;

 if (champ==NULL) {
  RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "champ nul dans calcul_histogrammes_VP\n");
 }

 if (champ->width!=wdth || champ->height!=hght || champ->depth!=dpth)
 {
  RAISE_ERR_TEXT(err, PB_TAILLE_CHAMP, "le champ n'est pas de la bonne taille dans calcul_histogrammes_VP\n");
 }

 data = champ->raw;

 //determination du nombre de bins de chaque histogramme
// nbcl1=donnees_VP->nb_cl1;
// nbcl2=donnees_VP->nb_cl2;
// max1=donnees_VP->max1;
// max2=donnees_VP->max2;
 nbcl1=dist_param->nbclReca; nbcl2=dist_param->nbclRef;
 max1=dist_param->maxReca; max2=dist_param->maxRef;

 //quantification de image reca non transformee
 imrecaQuantif=cr_grphic3d(Imreca);
 if (!imrecaQuantif) {
  RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation de imrecaQuantif dans calcul_histogrammes_VP\n");
 }
 imx_copie_param_3d_p(Imreca,imrecaQuantif);

{
 nbcl=nbcl1-1; max=max1+1;

 for (i=0;i<wdthReca;i++)
  for (j=0;j<hghtReca;j++)
   for (k=0;k<dpthReca;k++)
   {
    if (Imreca->mri[i][j][k]<=0) imrecaQuantif->mri[i][j][k]=0;
    else imrecaQuantif->mri[i][j][k]=MINI((nbcl*Imreca->mri[i][j][k])/max+1, nbcl);
   }   
}

 //initialisation de l'histogramme conjoint
 for (i=0;i<nbcl1;i++)
  for (j=0;j<nbcl2;j++)
   histo_conjoint[i][j]=0.0;

///////// calcul des histogrammes ///////////

 nbcl=nbcl2-1; max=max2+1;

 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    //calcul du voxel transforme
    if (data!=NULL) {
     x=(double)i+(double)data[i][j][k].x;
     y=(double)j+(double)data[i][j][k].y;
     z=(double)k+(double)data[i][j][k].z;
    } else {
      x=(double)i;
      y=(double)j;
      z=(double)k;
    }
    
    i2=(int)floor(x);j2=(int)floor(y);k2=(int)floor(z);
    i3=i2+1;j3=j2+1;k3=k2+1;

    if (Imref->mri[i][j][k]<=0) pixelRef=0;
    else pixelRef=MINI((nbcl*Imref->mri[i][j][k])/max+1, nbcl);
    //else pixelRef=((nbcl2-1)*Imref->mri[i][j][k])/(max2+1)+1;
    
    xrel=x-(double)i2; yrel=y-(double)j2; zrel=z-(double)k2;
    
    if (i2>=0 && i3<wdthReca && j2>=0 && j3<hghtReca && k2>=0 && k3<dpthReca)
    {// tout le voisinage est dans l'image
 
     //on trouve les valeurs du voisinage

     v[0]=imrecaQuantif->mri[i2][j2][k2];
     v[1]=imrecaQuantif->mri[i2][j2][k3];
     v[2]=imrecaQuantif->mri[i2][j3][k2];
     v[3]=imrecaQuantif->mri[i2][j3][k3];
     v[4]=imrecaQuantif->mri[i3][j2][k2];
     v[5]=imrecaQuantif->mri[i3][j2][k3];
     v[6]=imrecaQuantif->mri[i3][j3][k2];
     v[7]=imrecaQuantif->mri[i3][j3][k3];  

     //mise a jour de l'histogramme
     update_histogrammes_VP(pixelRef, donnees_VP, v, w, xrel, yrel, zrel);

    }
    else if (i3>=0 && i2<wdthReca && j3>=0 && j2<hghtReca && k3>=0 && k2<dpthReca)
    {//une partie du voisinage est dans l'image

     if (i2==-1 || j2==-1 || k2==-1) v[0]=0; else v[0]=imrecaQuantif->mri[i2][j2][k2];
     if (i2==-1 || j2==-1 || k3==dpthReca) v[1]=0; else v[1]=imrecaQuantif->mri[i2][j2][k3];
     if (i2==-1 || j3==hghtReca || k2==-1) v[2]=0; else v[2]=imrecaQuantif->mri[i2][j3][k2];
     if (i2==-1 || j3==hghtReca || k3==dpthReca) v[3]=0; else v[3]=imrecaQuantif->mri[i2][j3][k3];
     if (i3==wdthReca || j2==-1 || k2==-1) v[4]=0; else v[4]=imrecaQuantif->mri[i3][j2][k2];
     if (i3==wdthReca || j2==-1 || k3==dpthReca) v[5]=0; else v[5]=imrecaQuantif->mri[i3][j2][k3];
     if (i3==wdthReca || j3==hghtReca || k2==-1) v[6]=0; else v[6]=imrecaQuantif->mri[i3][j3][k2];
     if (i3==wdthReca || j3==hghtReca || k3==dpthReca) v[7]=0; else v[7]=imrecaQuantif->mri[i3][j3][k3];

     //mise a jour de l'histogramme
     update_histogrammes_VP(pixelRef, donnees_VP, v, w, xrel, yrel, zrel);
   
    }
    else
    {//tout le voisinage est hors de l'image

     //on calcule la classe du pixel de l'image non transformee.
     //Il y a une classe 0 pour les voxels nuls qui ne sont pas porteurs d'info
     histo_conjoint[0][pixelRef]+=1.0;

    }

   }

  if (err!=NULL) *err=NO_ERR;

end_func:

 //liberation memoire  
 if (imrecaQuantif) free_grphic3d(imrecaQuantif);
 
 return(1);
}

/*******************************************************************************
**    update_histogrammes_VP(pixelRef, donnees_VP, v, w, cl, xrel, yrel, zrel)
*/
/*!
**    mise a jour de l'histogramme conjoint en fonction des donnees relatives au pixels
**    courants dans imref et imreca transformee
**
**  \param pixelRef : valeur du pixel courant dans imref
**  \param donnees_VP : structure ou mettre les donnees calculees
**  \param v : tableau des valeurs des voisins du pixel courant dans imreca transformee
**  \param w : tableau des poids affectes aux valeurs des voisins du pixel courant dans imreca transformee
**  \param cl : tableau des classes de l'histogramme correspondant
**              aux voisins du pixel courant dans imreca transformee
**  \param xrel, yrel, zrel : distances des coordonnees du pixel courant dans imreca transformee
**                            a celles de son voisin du coin inferieur
**  \retval int : 1 en cas de succes, 0 dans le cas contraire
*******************************************************************************/
void update_histogrammes_VP(int pixelRef, VP_histos * donnees_VP, TYPEMRI3D *v, double *w, double xrel, double yrel, double zrel)
{
 int i;

 //on calcule les poids pour l'histogramme conjoint
 w[0]=(1.0-xrel)*(1.0-yrel)*(1.0-zrel);
 w[1]=(1.0-xrel)*(1.0-yrel)*zrel;
 w[2]=(1.0-xrel)*yrel*(1.0-zrel);
 w[3]=(1.0-xrel)*yrel*zrel;
 w[4]=xrel*(1.0-yrel)*(1.0-zrel);
 w[5]=xrel*(1.0-yrel)*zrel;
 w[6]=xrel*yrel*(1.0-zrel);
 w[7]=xrel*yrel*zrel;

 //update de l'histogramme conjoint
 for (i=0;i<8;i++) donnees_VP->histo_conjoint[v[i]][pixelRef]+=w[i];

}

/*******************************************************************************
**    calcul_histogrammes_VP(Imreca,Imref,champ,donnees_VP)
*/
/*!
**    fonction de calcul des histogrammes utilises dans le calcul du gradient
**    pour le volume partiel
**  \param imreca : image a recaler
**  \param champ : le champ
**  \param donnees_VP : structure ou mettre les donnees calculees
**  \retval int : 1 en cas de succes, 0 dans le cas contraire
*******************************************************************************/
int calcul_histogrammes_gradient_VP(grphic3d *Imreca, grphic3d *Imref, field3d *champ, VP_histos * donnees_VP, ERREUR_RECALAGE *err)
//int calcul_histogrammes_gradient_VP(ptr_distance_param dist_param)
{
 int wdth,hght,dpth;
 int wdthReca,hghtReca,dpthReca;
 int i,j,k,l,i2,j2,k2,i3,j3,k3;
 double x,y,z,xrel,yrel,zrel;
 TYPEMRI3D v[8];
 double dw[8][3];
 int max1,max2,nbcl1,nbcl2;
 vector3d ***data=NULL;
 int p[4] = {0};
 double **** histos_gradient=NULL;
 int pixelRef=0;
 double dxReca, dyReca, dzReca, volVoxelReca;
// grphic3d *Imreca=dist_param->imreca;
// grphic3d *Imref=dist_param->imref;
// field3d *champ=dist_param->champ;
// VP_histos * donnees_VP=dist_param->donnees_VP;
// ERREUR_RECALAGE *err=dist_param->err;

 grphic3d *imrecaQuantif=NULL;

 //////// Initialisations ////////
 wdthReca=Imreca->width;
 hghtReca=Imreca->height;
 dpthReca=Imreca->depth;

 dxReca=(double)Imreca->dx;
 dyReca=(double)Imreca->dy;
 dzReca=(double)Imreca->dz;
 volVoxelReca=dxReca*dyReca*dzReca;
 
 if (donnees_VP==NULL) {
   RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "la structure devant contenant l'histogramme conjoint pour le volume partiel n'est pas allouee dans calcul_histogrammes_gradient_VP\n");
  }
  
  histos_gradient = donnees_VP->histos_gradient;

  wdth=Imref->width;hght=Imref->height;dpth=Imref->depth;

 if (champ==NULL) {
  RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "champ nul dans calcul_histogrammes_gradient_VP\n");
 }

 if (champ->width!=wdth || champ->height!=hght || champ->depth!=dpth)
 {
  RAISE_ERR_TEXT(err, PB_TAILLE_CHAMP, "le champ n'est pas de la bonne taille dans calcul_histogrammes_gradient_VP\n");
 }

 data = champ->raw;

 //determination du nombre de bins de chaque histogramme
 nbcl1=donnees_VP->nb_cl1;
 nbcl2=donnees_VP->nb_cl2;
 max1=donnees_VP->max1;
 max2=donnees_VP->max2;
 
 //quantification de image reca transformee
 imrecaQuantif=cr_grphic3d(Imreca);
 if (!imrecaQuantif) {
  RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation de imrecaQuantif dans calcul_histogrammes_gradient_VP\n");
 }
 imx_copie_param_3d_p(Imreca,imrecaQuantif);

 for (i=0;i<wdthReca;i++)
  for (j=0;j<hghtReca;j++)
   for (k=0;k<dpthReca;k++)
   {
    if (Imreca->mri[i][j][k]<=0) imrecaQuantif->mri[i][j][k]=0;
    else imrecaQuantif->mri[i][j][k]=((nbcl1-1)*Imreca->mri[i][j][k])/(max1+1)+1;
   }

  //initialisation des histos pour le gradient
  for (i=0;i<3;i++)
   for (j=0;j<4;j++)   
    for (k=0;k<nbcl1;k++)
     for (l=0;l<nbcl2;l++)
      histos_gradient[i][j][k][l]=0.0;

//////// calcul des histogrammes ////////

 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    //calcul du voxel transforme
    if (data!=NULL) {
     x=(double)i+data[i][j][k].x;
     y=(double)j+data[i][j][k].y;
     z=(double)k+data[i][j][k].z;
    } else {
      x=(double)i;
      y=(double)j;
      z=(double)k;
    }

    i2=(int)floor(x);j2=(int)floor(y);k2=(int)floor(z);
    i3=i2+1;j3=j2+1;k3=k2+1;

    //calcul du vecteur representant le point 4D transforme par la matrice de transformation 3x4
    p[0]=i;
    p[1]=j;
    p[2]=k;
    p[3]=1;

    if (Imref->mri[i][j][k]<=0) pixelRef=0;
    else pixelRef=((nbcl2-1)*Imref->mri[i][j][k])/(max2+1)+1;

    xrel=(x-(double)i2)*dxReca/volVoxelReca;yrel=(y-(double)j2)*dyReca/volVoxelReca;zrel=(z-(double)k2)*dzReca/volVoxelReca;
    
    if (i2>=0 && i3<wdthReca && j2>=0 && j3<hghtReca && k2>=0 && k3<dpthReca)
    {// tout le voisinage est dans l'image

     //on trouve les valeurs du voisinage
     v[0]=imrecaQuantif->mri[i2][j2][k2];
     v[1]=imrecaQuantif->mri[i2][j2][k3];
     v[2]=imrecaQuantif->mri[i2][j3][k2];
     v[3]=imrecaQuantif->mri[i2][j3][k3];
     v[4]=imrecaQuantif->mri[i3][j2][k2];
     v[5]=imrecaQuantif->mri[i3][j2][k3];
     v[6]=imrecaQuantif->mri[i3][j3][k2];
     v[7]=imrecaQuantif->mri[i3][j3][k3];   

    }
    else if (i3>=0 && i2<wdthReca && j3>=0 && j2<hghtReca && k3>=0 && k2<dpthReca)
    {//une partie du voisinage est dans l'image

     if (i2==-1 || j2==-1 || k2==-1) v[0]=0; else v[0]=imrecaQuantif->mri[i2][j2][k2];
     if (i2==-1 || j2==-1 || k3==dpthReca) v[1]=0; else v[1]=imrecaQuantif->mri[i2][j2][k3];
     if (i2==-1 || j3==hghtReca || k2==-1) v[2]=0; else v[2]=imrecaQuantif->mri[i2][j3][k2];
     if (i2==-1 || j3==hghtReca || k3==dpthReca) v[3]=0; else v[3]=imrecaQuantif->mri[i2][j3][k3];
     if (i3==wdthReca || j2==-1 || k2==-1) v[4]=0; else v[4]=imrecaQuantif->mri[i3][j2][k2];
     if (i3==wdthReca || j2==-1 || k3==dpthReca) v[5]=0; else v[5]=imrecaQuantif->mri[i3][j2][k3];
     if (i3==wdthReca || j3==hghtReca || k2==-1) v[6]=0; else v[6]=imrecaQuantif->mri[i3][j3][k2];
     if (i3==wdthReca || j3==hghtReca || k3==dpthReca) v[7]=0; else v[7]=imrecaQuantif->mri[i3][j3][k3];    
       
    }
    else { v[0]=v[1]=v[2]=v[3]=v[4]=v[5]=v[6]=v[7]=0; }

    //mise a jour des histogrammes de gradient
    update_histogrammes_gradient_VP(pixelRef, donnees_VP, v, dw, p, xrel, yrel, zrel);
  
 }

  if (err!=NULL) *err=NO_ERR;

end_func:

 //liberation memoire
 if (imrecaQuantif) free_grphic3d(imrecaQuantif); 
 
 return(1);
 
}

/*******************************************************************************
**    update_histogrammes_gradient_VP(pixelRef, donnees_VP, v, dw, cl, p, xrel, yrel, zrel)
*/
/*!
**    mise a jour de l'histogramme conjoint en fonction des donnees relatives au pixels
**    courants dans imref et imreca transformee
**
**  \param pixelRef : valeur du pixel courant dans imref
**  \param donnees_VP : structure ou mettre les donnees calculees
**  \param v : tableau des valeurs des voisins du pixel courant dans imreca transformee
**  \param dw : tableau des poids affectes aux valeurs des voisins du pixel courant dans imreca transformee
**  \param cl : tableau des classes de l'histogramme correspondant
**              aux voisins du pixel courant dans imreca transformee
**  \param p : coordonnees du point courant dans imreca transformee
**  \param xrel, yrel, zrel : distances des coordonnees du pixel courant dans imreca transformee
**                            a celles de son voisin du coin inferieur
**  \retval int : 1 en cas de succes, 0 dans le cas contraire
*******************************************************************************/
void update_histogrammes_gradient_VP(int pixelRef, VP_histos * donnees_VP, TYPEMRI3D *v, double dw[8][3], int *p, double xrel, double yrel, double zrel)
{
 int i,l,m;
  
 //on calcule les poids pour l'histogramme conjoint
 dw[0][0]=-(1.0-yrel)*(1.0-zrel); dw[4][0]=-dw[0][0];
 dw[1][0]=-(1.0-yrel)*zrel;       dw[5][0]=-dw[1][0];
 dw[2][0]=-yrel*(1.0-zrel);       dw[6][0]=-dw[2][0];
 dw[3][0]=-yrel*zrel;             dw[7][0]=-dw[3][0];
 dw[0][1]=-(1.0-xrel)*(1.0-zrel); dw[2][1]=-dw[0][1];
 dw[1][1]=-(1.0-xrel)*zrel;       dw[3][1]=-dw[1][1];
 dw[4][1]=-xrel*(1.0-zrel);       dw[6][1]=-dw[4][1];
 dw[5][1]=-xrel*zrel;             dw[7][1]=-dw[5][1];
 dw[0][2]=-(1.0-xrel)*(1.0-yrel); dw[1][2]=-dw[0][2];
 dw[2][2]=-(1.0-xrel)*yrel;       dw[3][2]=-dw[2][2];
 dw[4][2]=-xrel*(1.0-zrel);       dw[5][2]=-dw[4][2];
 dw[6][2]=-xrel*yrel;             dw[7][2]=-dw[6][2];

 //mise a jour des histogrammes de gradient  
 for (l=0;l<3;l++)
  for (m=0;m<4;m++)
   for (i=0;i<8;i++)
    donnees_VP->histos_gradient[l][m][v[i]][pixelRef]+=dw[i][l]*(double)p[m];

}

/*******************************************************************************
**     calcul_info_mutuelle_VP_stochastique_3d(im_type, Imreca, Imref, champ, donnees_VP)
*/
/*!    calcul de la distance se basant sur l'information mutuelle calculee avec le volume partiel et
**     un choix de points de calculs aleatoires
    \param im_type: le type de distance utilisee
        \param Imreca, Imref : les 2 images
        \param champ : le champ de deformation
        \param donnees_VP : structure contenant les donnees pour le volume partiel
        \retval double : l'erreur d'information mut.
*******************************************************************************/
double calcul_info_mutuelle_VP_stochastique_3d(IM3DTYPE im_type, ptr_distance_param dist_param)
{
 int nbcl1,nbcl2;
 int i,j;
 double **histo_conjoint=NULL,*histo1=NULL,*histo2=NULL;
 double info_mutuelle=DBL_MAX, entropie_conjointe=0.0, entropie1=0.0, entropie2=0.0;
 double nbVoxelsOverlap=0.0, nbVoxelsHorsOverlap=0.0;

 VP_histos * donnees_VP=dist_param->donnees_VP;
 field3d *points_calcules=dist_param->points_calcules;
 ERREUR_RECALAGE *err=dist_param->err;

 if (donnees_VP==NULL)
 {
   RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "la structure devant contenant l'histogramme conjoint pour le volume partiel n'est pas allouee dans calcul_info_mutuelle_VP_stochastique_3d\n");
  }

 if (points_calcules==NULL)
 {
   RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "le tableau de points devant etre utilises pour le calcul de l'information mutuelle n'a pas ete cree dans calcul_info_mutuelle_VP_stochastique_3d\n");
 }

 histo_conjoint=donnees_VP->histo_conjoint;
 histo1=donnees_VP->histo_marginal1;
 histo2=donnees_VP->histo_marginal2;

 //appel a la fonction de calcul des histogrammes appropriee
// calcul_histogrammes_VP_stochastique(Imreca, Imref, champ, donnees_VP, points_calcules, err);
 calcul_histogrammes_VP_stochastique(dist_param);
 if ((err)&&(*err!=NO_ERR)) goto end_func;

// nbcl1 = donnees_VP->nb_cl1;
// nbcl2 = donnees_VP->nb_cl2;
 nbcl1=dist_param->nbclReca; nbcl2=dist_param->nbclRef;
// max1=dist_param->maxReca; max2=dist_param->maxRef;

 //on calcule le nb de voxel dans la zone de recouvrement
 for (i=1;i<nbcl1;i++)
  for (j=1;j<nbcl2;j++)
   nbVoxelsOverlap += histo_conjoint[i][j];


 //si aucun recouvrement on renvoie le max de l'energie
 if (nbVoxelsOverlap==0.0)
 {
  info_mutuelle=max_erreur_info_mutuelle(IM, nbcl1, nbcl2);
  RAISE_ERR_TEXT(err, OVERLAP_NUL, "zone d'overlap nulle dans calcul_info_mutuelle_VP_stochastique_3d\n");
 }

 nbVoxelsHorsOverlap = (double)(points_calcules->width*points_calcules->height*points_calcules->depth)-nbVoxelsOverlap;

 info_mutuelle=calc_info_mutuelle_histos(im_type, histo_conjoint, histo1, histo2, nbcl1, nbcl2, nbVoxelsOverlap, nbVoxelsHorsOverlap, &entropie_conjointe, &entropie1, &entropie2);

 //update des donnees
 donnees_VP->entropie_conjointe = entropie_conjointe;
 donnees_VP->entropie1 = entropie1;
 donnees_VP->entropie2 = entropie2;

 if (err!=NULL) *err=NO_ERR;

end_func:

 return info_mutuelle;
}

/*******************************************************************************
**        choix_points_aleatoires(width, height, depth)
**
**        choix de points uniformement repartis dans une image pour les calcul
**        stochastique de la distance inter image
**
**        width, height, depth:dimensions des images a comparer
********************************************************************************/

field3d *choix_points_aleatoires(int width, int height, int depth, ERREUR_RECALAGE *err)
{
 int i,j,k;
 field3d *points_choisis=NULL;
 float positionx,positiony,positionz;

 if ((width<2)||(height<2)||(depth<2))
 {
  RAISE_ERR_TEXT(err, PB_TAILLES_IMAGES, "image de dimensions trop petite dans choix_points_aleatoires\n");
 }

 //les points sont choisis uniformement dans tous les cubes formes par 8 voxels adjacents
 //d'ou la taille du field3d pour les stocker
 points_choisis=cr_field3d(width-1,height-1,depth-1);
 if (!points_choisis)
 {
  RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans choix_points_aleatoires\n");
 }

 for (i=0;i<(width-1);i++)
 {
  for (j=0;j<(height-1);j++)
  {
   for (k=0;k<(depth-1);k++)
   {
    //ici on utilise une astuce pour etre sur que les points ne soient pas sur les bords superieurs des images
    //configuration qui pourrait faire planter certains algos
    positionx = (float)(rand()/RAND_MAX*0.998);
    (points_choisis->raw[i][j][k]).x=(float)(i+positionx+0.001);
    positiony = (float)(rand()/RAND_MAX*0.998);
    (points_choisis->raw[i][j][k]).y=(float)(j+positiony+0.001);
    positionz = (float)(rand()/RAND_MAX*0.998);
    (points_choisis->raw[i][j][k]).z=(float)(k+positionz+0.001);
   }
  }
 }

 if (err) *err=NO_ERR;

end_func:

 return points_choisis;
}

/*******************************************************************************
**        calcul_histogrammes_VP2(Imreca, Imref, points_choisis, champ)
**
**        calcul des histos de l'information mutuelle, methode stochastique
**
**        Imreca, Imref: images a comparer
**        points_choisis: points choisis aleatoirement
**        champ: champs de la transformation appliquee
********************************************************************************/
int calcul_histogrammes_VP_stochastique(ptr_distance_param dist_param)
{
 int wdth,hght,dpth;
 int wdthReca,hghtReca,dpthReca;
 int wdthRef,hghtRef,dpthRef;
 int i,j,k,i2,j2,k2,i3,j3,k3;
 int i22,j22,k22,i32,j32,k32;
 double x,y,z,xrel,yrel,zrel;
 double x2,y2,z2,xrel2,yrel2,zrel2;
 double w[8];
 double w2[8];
 TYPEMRI3D v[8];
 TYPEMRI3D v2[8];
 int max1,max2,nbcl1,nbcl2;
 vector3d ***data=NULL, ***pts_coords=NULL;
 double ** histo_conjoint=NULL;
 grphic3d *imrecaQuantif=NULL, *imrefQuantif=NULL;

 grphic3d *Imreca=dist_param->imreca;
 grphic3d *Imref=dist_param->imref;
 field3d *champ=dist_param->champ;
 VP_histos * donnees_VP=dist_param->donnees_VP;
 field3d *points_choisis=dist_param->points_calcules;
 ERREUR_RECALAGE *err=dist_param->err;

 ///////// Initialisations /////////
 wdthReca=Imreca->width; hghtReca=Imreca->height; dpthReca=Imreca->depth;
 wdthRef=Imref->width; hghtRef=Imref->height; dpthRef=Imref->depth;

 if (donnees_VP==NULL)
 {
   RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "la structure devant contenant l'histogramme conjoint pour le volume partiel n'est pas allouee dans calcul_histogrammes_VP_stochastique\n");
 }

  histo_conjoint = donnees_VP->histo_conjoint;

  wdth=points_choisis->width;hght=points_choisis->height;dpth=points_choisis->depth;

 if (champ==NULL) {
  RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "champ nul dans calcul_histogrammes_VP_stochastique\n");
 }

 if (champ->width!=(wdth+1) || champ->height!=(hght+1) || champ->depth!=(dpth+1))
 {
  RAISE_ERR_TEXT(err, PB_TAILLE_CHAMP, "le champ n'est pas de la bonne taille dans calcul_histogrammes_VP_stochastique\n");
 }

 data = champ->raw;
 
 pts_coords=points_choisis->raw;

 //determination du nombre de bins de chaque histogramme
 nbcl1=dist_param->nbclReca; nbcl2=dist_param->nbclRef;
 max1=dist_param->maxReca; max2=dist_param->maxRef;
 
 //allocation des images quantifiee
 imrefQuantif=cr_grphic3d(Imref);
 if (!imrefQuantif) {
  RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation de imrefQuantif dans calcul_histogrammes_VP_stochastique\n");
 }
 imx_copie_param_3d_p(Imref,imrefQuantif);
 imrecaQuantif=cr_grphic3d(Imreca);
 if (!imrecaQuantif) {
  RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation de imrecaQuantif dans calcul_histogrammes_VP_stochastique\n");
 }
 imx_copie_param_3d_p(Imreca,imrecaQuantif);

 //quantification des images
 for (i=0;i<wdthRef;i++)
  for (j=0;j<hghtRef;j++)
   for (k=0;k<dpthRef;k++)
   {
    if (Imref->mri[i][j][k]<=0) imrefQuantif->mri[i][j][k]=0;
    else imrefQuantif->mri[i][j][k]=MINI(((nbcl2-1)*Imref->mri[i][j][k])/(max2+1)+1, (nbcl2-1));
   }

 for (i=0;i<wdthReca;i++)
  for (j=0;j<hghtReca;j++)
   for (k=0;k<dpthReca;k++)
   {
    if (Imreca->mri[i][j][k]<=0) imrecaQuantif->mri[i][j][k]=0;
    else imrecaQuantif->mri[i][j][k]=MINI(((nbcl1-1)*Imreca->mri[i][j][k])/(max1+1)+1, (nbcl1-1));
   }

 //initialisation de l'histogramme conjoint
 for (i=0;i<nbcl1;i++)
  for (j=0;j<nbcl2;j++)
   histo_conjoint[i][j]=0.0;

///////// calcul des histogrammes ///////////

 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {                                   
    //calcul du voxel transforme
    if (data!=NULL) {
     x=(double)pts_coords[i][j][k].x+(double)data[i][j][k].x;
     y=(double)pts_coords[i][j][k].y+(double)data[i][j][k].y;
     z=(double)pts_coords[i][j][k].z+(double)data[i][j][k].z;
    } else {
      x=(double)pts_coords[i][j][k].x;
      y=(double)pts_coords[i][j][k].y;
      z=(double)pts_coords[i][j][k].z;
    }

    x2=(double)pts_coords[i][j][k].x;
    y2=(double)pts_coords[i][j][k].y;
    z2=(double)pts_coords[i][j][k].z;

    i22=(int)floor(x2);j22=(int)floor(y2);k22=(int)floor(z2);
    i32=i22+1;j32=j22+1;k32=k22+1;
    xrel2=x2-(double)i22; yrel2=y2-(double)j22; zrel2=z2-(double)k22;
    
    v2[0]=imrefQuantif->mri[i22][j22][k22];
    v2[1]=imrefQuantif->mri[i22][j22][k32];
    v2[2]=imrefQuantif->mri[i22][j32][k22];
    v2[3]=imrefQuantif->mri[i22][j32][k32];
    v2[4]=imrefQuantif->mri[i32][j22][k22];
    v2[5]=imrefQuantif->mri[i32][j22][k32];
    v2[6]=imrefQuantif->mri[i32][j32][k22];
    v2[7]=imrefQuantif->mri[i32][j32][k32];

    i2=(int)floor(x);j2=(int)floor(y);k2=(int)floor(z);
    i3=i2+1;j3=j2+1;k3=k2+1;
    xrel=x-(double)i2; yrel=y-(double)j2; zrel=z-(double)k2;

    if (i2>=0 && i3<wdthReca && j2>=0 && j3<hghtReca && k2>=0 && k3<dpthReca)
    {// tout le voisinage est dans l'image

     //on trouve les valeurs du voisinage
     v[0]=imrecaQuantif->mri[i2][j2][k2];
     v[1]=imrecaQuantif->mri[i2][j2][k3];
     v[2]=imrecaQuantif->mri[i2][j3][k2];
     v[3]=imrecaQuantif->mri[i2][j3][k3];
     v[4]=imrecaQuantif->mri[i3][j2][k2];
     v[5]=imrecaQuantif->mri[i3][j2][k3];
     v[6]=imrecaQuantif->mri[i3][j3][k2];
     v[7]=imrecaQuantif->mri[i3][j3][k3];

    }
    else if (i3>=0 && i2<wdthReca && j3>=0 && j2<hghtReca && k3>=0 && k2<dpthReca)
    {//une partie du voisinage est dans l'image

     if (i2==-1 || j2==-1 || k2==-1) v[0]=0; else v[0]=imrecaQuantif->mri[i2][j2][k2];
     if (i2==-1 || j2==-1 || k3==dpthReca) v[1]=0; else v[1]=imrecaQuantif->mri[i2][j2][k3];
     if (i2==-1 || j3==hghtReca || k2==-1) v[2]=0; else v[2]=imrecaQuantif->mri[i2][j3][k2];
     if (i2==-1 || j3==hghtReca || k3==dpthReca) v[3]=0; else v[3]=imrecaQuantif->mri[i2][j3][k3];
     if (i3==wdthReca || j2==-1 || k2==-1) v[4]=0; else v[4]=imrecaQuantif->mri[i3][j2][k2];
     if (i3==wdthReca || j2==-1 || k3==dpthReca) v[5]=0; else v[5]=imrecaQuantif->mri[i3][j2][k3];
     if (i3==wdthReca || j3==hghtReca || k2==-1) v[6]=0; else v[6]=imrecaQuantif->mri[i3][j3][k2];
     if (i3==wdthReca || j3==hghtReca || k3==dpthReca) v[7]=0; else v[7]=imrecaQuantif->mri[i3][j3][k3];

    }
    else
    {//tout le voisinage est hors de l'image

     v[0]=v[1]=v[2]=v[3]=v[4]=v[5]=v[6]=v[7]=0;

    }

    //mise a jour de l'histogramme
    update_histogrammes_VP_stochastique(donnees_VP, v2, v, w2, w, xrel2, yrel2, zrel2, xrel, yrel, zrel);

   }

  if (err) *err=NO_ERR;

end_func:
   
 //libreations memoire  
 if (imrecaQuantif) free_grphic3d(imrecaQuantif);
 if (imrefQuantif) free_grphic3d(imrefQuantif);

 return(0);
}

/*******************************************************************************
**    update_histogrammes_VP_stochastique(donnees_VP, vRef, vReca, wRef, wReca, clRef, clReca, xrelRef, yrelRef, zrelRef, xrelReca, yrelReca, zrelReca)
*/
/*!
**    mise a jour de l'histogramme conjoint en fonction des donnees relatives au pixels
**    courants dans imref et imreca transformee
**
**  \param donnees_VP : structure ou mettre les donnees calculees
**  \param vRef,vReca : tableaux des valeurs des voisins du pixel courant dans imref et imreca transformee alloues
**  \param wRef,wReca : tableaux des poids affectes aux valeurs des voisins du pixel courant dans imref et imreca transformee alloues
**  \param xrelRef, yrelRef, zrelRef, xrelReca, yrelReca, zrelReca : distances des coordonnees du pixel courant dans imref et imreca transformee
**                                                                   a celles de son voisin du coin inferieur
*******************************************************************************/
void update_histogrammes_VP_stochastique(VP_histos * donnees_VP, TYPEMRI3D *vRef, TYPEMRI3D *vReca, double *wRef, double *wReca, double xrelRef, double yrelRef, double zrelRef, double xrelReca, double yrelReca, double zrelReca)
{
 int i,j;

 //on calcule les poids pour l'histogramme conjoint
 wRef[0]=(1.0-xrelRef)*(1.0-yrelRef)*(1.0-zrelRef);
 wRef[1]=(1.0-xrelRef)*(1.0-yrelRef)*zrelRef;
 wRef[2]=(1.0-xrelRef)*yrelRef*(1.0-zrelRef);
 wRef[3]=(1.0-xrelRef)*yrelRef*zrelRef;
 wRef[4]=xrelRef*(1.0-yrelRef)*(1.0-zrelRef);
 wRef[5]=xrelRef*(1.0-yrelRef)*zrelRef;
 wRef[6]=xrelRef*yrelRef*(1.0-zrelRef);
 wRef[7]=xrelRef*yrelRef*zrelRef;

 wReca[0]=(1.0-xrelReca)*(1.0-yrelReca)*(1.0-zrelReca);
 wReca[1]=(1.0-xrelReca)*(1.0-yrelReca)*zrelReca;
 wReca[2]=(1.0-xrelReca)*yrelReca*(1.0-zrelReca);
 wReca[3]=(1.0-xrelReca)*yrelReca*zrelReca;
 wReca[4]=xrelReca*(1.0-yrelReca)*(1.0-zrelReca);
 wReca[5]=xrelReca*(1.0-yrelReca)*zrelReca;
 wReca[6]=xrelReca*yrelReca*(1.0-zrelReca);
 wReca[7]=xrelReca*yrelReca*zrelReca;
 
 for (i=0;i<8;i++)
  for (j=0;j<8;j++)
  {
   donnees_VP->histo_conjoint[vReca[i]][vRef[j]]+=wReca[i]*wRef[j];
  }

}

/*******************************************************************************
**     calcul_info_mutuelle_VP_3d(im_type, Imreca, Imref, champ, donnees_VP)
*/
/*!    calcul de la distance se basant sur l'information mutuelle
**     calculee avec la methode des BSplines
    \param im_type: le type de distance utilisee
        \param Imreca, Imref : les 2 images
        \param champ : le champ de deformation
        \param donnees_VP : structure contenant les donnees pour le volume partiel
        \retval double : l'erreur d'information mut.
*******************************************************************************/
double calcul_info_mutuelle_BSpline_3d(IM3DTYPE im_type, ptr_distance_param dist_param)
{
 int nbclRef,nbclReca;
 int i,j;//,k;
 double **histo_conjoint=NULL,*histo1=NULL,*histo2=NULL;
 double info_mutuelle=DBL_MAX, entropie_conjointe=0.0, entropie1=0.0, entropie2=0.0;
 double nbVoxelsOverlap=0.0, nbVoxelsHorsOverlap=0.0;
 int maxRef, maxReca;
 int wdth,hght,dpth;

 grphic3d *Imreca=dist_param->imreca;
 grphic3d *Imref=dist_param->imref;
 field3d *champ=dist_param->champ;
 ERREUR_RECALAGE *err=dist_param->err;

 wdth=Imref->width; hght=Imref->height; dpth=Imref->depth;

 //calcul des maxpixel dans l'image transformee et l'image ref
 //et des nbs de classes pour l'estimation des densites de proba
// nbclRef=imx_calculer_nombre_de_classes_3d_p(Imref, &maxRef);
// nbclReca=imx_calculer_nombre_de_classes_3d_p(Imreca, &maxReca);
 nbclReca=dist_param->nbclReca; nbclRef=dist_param->nbclRef;
 maxReca=dist_param->maxReca; maxRef=dist_param->maxRef;

 if ((nbclRef<2)||(nbclReca<2))
 {
  RAISE_ERR_TEXT(err, NB_BINS_NUL, "nb de classes nulles dans calcul_info_mutuelle_BSpline_3d\n");
 }

 //allocation memoire et initialisation des histos
 histo1=CALLOC(nbclRef,double);histo2=CALLOC(nbclReca,double);
 histo_conjoint=alloc_dmatrix(nbclRef,nbclReca);
 if ((!histo1)||(!histo2)||(!histo_conjoint))
 {
  RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans calcul_info_mutuelle_BSpline_3d\n");
 }

 //appel a la fonction de calcul des histogrammes appropriee
 calcul_histogramme_BSpline(Imreca, Imref, champ, histo_conjoint, nbclReca, nbclRef, maxReca, maxRef);

 //on calcule le nb de voxel dans la zone de recouvrement
 for (i=1;i<nbclRef;i++)
  for (j=1;j<nbclReca;j++)
   nbVoxelsOverlap += histo_conjoint[i][j];

 //si aucun recouvrement on renvoie le max de l'energie
 if (nbVoxelsOverlap==0.0)
 {
  info_mutuelle=max_erreur_info_mutuelle(im_type, nbclReca, nbclRef);
  RAISE_ERR_TEXT(err, OVERLAP_NUL, "zone d'overlap nulle dans calcul_info_mutuelle_BSpline_3d\n");
 }

 nbVoxelsHorsOverlap = (double)(Imreca->width*Imreca->height*Imreca->depth)-nbVoxelsOverlap;
/*
 //calcul des histogrammes marginaux
 calcul_histos_marginaux(histo_conjoint, histo1, histo2, nbclRef, nbclReca, nbVoxelsHorsOverlap);

 //calcul des entropies
 entropie_conjointe = calcul_entropie_conjointe(histo_conjoint, nbclRef, nbclReca, nbVoxelsOverlap);
 entropie1 = calcul_entropie(histo1, nbclRef, nbVoxelsOverlap);
 entropie2 = calcul_entropie(histo2, nbclReca, nbVoxelsOverlap);

 //calcul de la distance
 switch(im_type)
 {
  case ENTROPIE_CONJOINTE: info_mutuelle=100.0*entropie_conjointe; break;
  case IM: info_mutuelle=calcul_IM(entropie_conjointe, entropie1, entropie2); info_mutuelle=-1000.0*info_mutuelle; break;
  case IMNS: info_mutuelle=calcul_IMNS(entropie_conjointe, entropie1, entropie2); info_mutuelle=-100.0*info_mutuelle; break;
  case IMNM1: info_mutuelle=calcul_IMNM1(entropie_conjointe, entropie1, entropie2); info_mutuelle=-1000.0*info_mutuelle; break;
  case IMNM2: info_mutuelle=calcul_IMNM2(entropie_conjointe, entropie1, entropie2); info_mutuelle=100.0*info_mutuelle; break;
  default: info_mutuelle=100.0*entropie_conjointe; break;
 }
*/
  info_mutuelle=calc_info_mutuelle_histos(im_type, histo_conjoint, histo1, histo2, nbclRef, nbclReca, nbVoxelsOverlap, nbVoxelsHorsOverlap, &entropie_conjointe, &entropie1, &entropie2);

  if (err) *err=NO_ERR;

end_func:

 //liberation memoire
 if (histo_conjoint) free_dmatrix(histo_conjoint,nbclReca,nbclReca);
 if (histo1) FREE(histo1);
 if (histo2) FREE(histo2);

 return info_mutuelle;
}

/*******************************************************************************
**    calcul_histogrammes_VP(Imreca,Imref,champ,donnees_VP)
*/
/*!
**    fonction de calcul de l'histogramme conjoint pour la methode de construction par des BSplines
**    ici on utilise les BSplines cubiques
**
**  \param imreca : image a recaler
**  \param champ : le champ
**  \param donnees_VP : structure ou mettre les donnees calculees
**  \retval int : 1 en cas de succes, 0 dans le cas contraire
*******************************************************************************/
int calcul_histogramme_BSpline(grphic3d *Imreca, grphic3d *Imref, field3d *champ, double **histo_conjoint, int nbclReca, int nbclRef, int maxReca, int maxRef)
{
 int wdth,hght,dpth;
 int i,j,k,l,m;
 double x;
 double val1, val2;
 double valRef[4], valReca[4];
 int indexRef[4], indexReca[4];
 int nbRef, nbReca;
 double tot=0.0;
 TYPEMRI3D mRef, mReca;

///////// Initialisations /////////
 if ((nbclReca<=4)||(nbclRef<=4))
 {
   fprintf (stderr, "pas assez de classes pour le calcul\n");
   return 1;
 }
 wdth=Imref->width;hght=Imref->height;dpth=Imref->depth;

 for (i=0;i<nbclRef;i++)
  for (j=0;j<nbclReca;j++)
   histo_conjoint[i][j]=0.0;

 //calcul de l'histogramme
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
    {
     mRef=Imref->mri[i][j][k]; mReca=Imreca->mri[i][j][k];

     //on ne fait les calculs que sur la zone de recouvrement des images
     if ((mRef>0)&&(mReca>0))
     {

     //calcul des valeurs reelles des images requantifiees dans ]2,(nbcl-2)[
     val1=(double)((mRef+1)*(nbclRef-4))/(double)(maxRef+2)+2.0;
     val2=(double)((mReca+1)*(nbclReca-4))/(double)(maxReca+2)+2.0;

     //calcul des bins affectes par imref et des valeurs correspondantes

     //si val1 est entier seulement 3 bins vont etre affectes sinon c'est 4
     if (floor(val1)==ceil(val1)) nbRef=3; else nbRef=4;

     //les valeurs d'update de l'histogrames sont calulees en appliquant la fonction BSpline cubique
     //a l'ecart entre la valeur reelle de val1 et la valeur du bin update (nombre entier)
     if (nbRef==3)
     {
      indexRef[0]=MAXI((int)floor(val1)-1,0);
      x=val1-floor(val1)+1.0;
      valRef[0]=1.0/6.0;

      indexRef[1]=MAXI((int)floor(val1),0);
      x=val1-floor(val1);
      valRef[1]=2.0/3.0;

      indexRef[2]=MAXI((int)ceil(val1)+1,0);
      x=1.0;
      valRef[2]=1.0/6.0;
     }
     else
     {
      indexRef[0]=MAXI((int)floor(val1)-1,0);
      x=val1-floor(val1)+1.0;
      valRef[0]=(2.0-x)*(2.0-x)*(2.0-x)/6.0;

      indexRef[1]=MAXI((int)floor(val1),0);
      x=val1-floor(val1);
      valRef[1]=2.0/3.0-x*x*(2.0-x)/2.0;

      indexRef[2]=MAXI((int)ceil(val1),0);
      x=ceil(val1)-val1;
      valRef[2]=2.0/3.0-x*x*(2.0-x)/2.0;

      indexRef[3]=MAXI((int)ceil(val1)+1,0);
      x=ceil(val1)+1.0-val1;
      valRef[3]=(2.0-x)*(2.0-x)*(2.0-x)/6.0;
     }

     //calcul des bins affectes par imreca et des valeurs correspondantes

     //si val2 est entier seulement 3 bins vont etre affectes sinon c'est 4
     if (floor(val2)==ceil(val2)) nbReca=3; else nbReca=4;

     //les valeurs d'update de l'histogrames sont calulees en appliquant la fonction BSpline cubique
     //a l'ecart entre la valeur reelle de val1 et la valeur du bin update (nombre entier)
     if (nbReca==3)
     {
      indexReca[0]=MAXI((int)floor(val2)-1,0);
      x=1.0;
      valReca[0]=1.0/6.0;

      indexReca[1]=MAXI((int)floor(val2),0);
      x=val2-floor(val2);
      valReca[1]=2.0/3.0;

      indexReca[2]=MAXI((int)ceil(val2)+1,0);
      x=1.0;
      valReca[2]=1.0/6.0;
     }
     else
     {
      indexReca[0]=MAXI((int)floor(val2)-1,0);
      x=val2-floor(val2)+1.0;
      valReca[0]=(2.0-x)*(2.0-x)*(2.0-x)/6.0;

      indexReca[1]=MAXI((int)floor(val2),0);
      x=val2-floor(val2);
      valReca[1]=2.0/3.0-x*x*(2.0-x)/2.0;

      indexReca[2]=MAXI((int)ceil(val2),0);
      x=ceil(val2)-val2;
      valReca[2]=2.0/3.0-x*x*(2.0-x)/2.0;

      indexReca[3]=MAXI((int)ceil(val2)+1,0);
      x=ceil(val2)+1.0-val2;
      valReca[3]=(2.0-x)*(2.0-x)*(2.0-x)/6.0;
     }

     //update des bin d'histogramme (indiques par indexRef et indexReca)
     //par les valeur calculees (indiquees par valRef et valReca)
     tot=0.0;
     for (l=0;l<nbRef;l++)
      for (m=0;m<nbReca;m++)
      {
       histo_conjoint[indexRef[l]][indexReca[m]]+=valRef[l]*valReca[m];
       tot += valRef[l]*valReca[m];
      }

     }

    }

 return 0;
}

/*******************************************************************************
**    max_erreur_info_mutuelle(im_type, nbcl1, nbcl2)
*/
/*!
**    renvoie le max de l'energie basee sur l'info mutuelle selon les nbs de classes utilisees
**
**  \param im_type: type d'energie utilisee
**  \param nbcl1, nbcl2 : les nombres de classes pour le calcul de l'entropie conjointe
*******************************************************************************/
double max_erreur_info_mutuelle(IM3DTYPE im_type, int nbcl1, int nbcl2)
{
 double max;

 switch (im_type)
 {
   case ENTROPIE_CONJOINTE: max=-log(1.0/(double)(nbcl1*nbcl2)); max=100.0*max; break;
   case IM: max=0.0; max=-1000.0*max; break;
   case IMNS: max=0.0; break;
   case IMNM1: max=0.0; max=-1000.0*max; break;
   case IMNM2: max=-log(1.0/(double)(nbcl1*nbcl2)); max=100.0*max; break;
   default: max=-log(1.0/(double)(nbcl1*nbcl2)); max=100.0*max; break;
 }

 return max;
}

/*******************************************************************************
**    calcul_nb_classes(nb_classes_min, wdth, hght, dpth, max1, max2, nbcl1, nbcl2)
*/
/*!
**    calcul le nombre de classes de valeurs de voxels utilises pour les calcul d'energie selon la resolution de l'image
**
**  \param nb_classes_max: nombre de classes a la resolution min
**  \param wdth, hght, dpth : dimensions des images
**  \param max1, max2 : valeurs max des images
**  \param nbcl1, nbcl2 : nbs de classes calcules
*******************************************************************************/
int calcul_nb_classes(int wdth, int hght, int dpth, int max1, int max2, int *nbcl1, int *nbcl2)
{
  
 *nbcl1=*nbcl2=(int)floor(pow(wdth*hght*dpth,1.0/3.0));
 if (*nbcl1>max1) *nbcl1=(max1+1);
 if (*nbcl2>max2) *nbcl2=(max2+1);

 return 0;
}

/*******************************************************************************
**    imx_calculer_nombre_de_classes_3d_p(image, max_image)
*/
/*!
**    calcul le nombre de classes utilis�es pour le calcul de l'entropie d'une image
**    c'est le nb de classes permettant une estimation non biaisee de la densite de proba
**    d'apres Jenkinson et Smith dans:
**    "A global optimisation method for robust affine registration of brain image"
**
**  \param image: l'image consideree
**  \param max_image : la valeur max de l'image
**  \retval : le nb de classes calcule
*******************************************************************************/
int imx_calculer_nombre_de_classes_3d_p(grphic3d *image, int *max_image)
{
 int wdth, hght, dpth;
 int i,j,k;
 double nbvoxelsOverlap=0.0, sigma=0.0/*, nbcl=0.0*/, w_bin=0.0, val=0.0;
 double moy=0.0,moy2=0.0;
 double eps=1e-4;
 TYPEMRI3D ***imageMRI=NULL;
 TYPEMRI3D val1=0, max1=0;
 int i_nbcl=0;

 wdth=(int)image->width; hght=(int)image->height; dpth=(int)image->depth;
 imageMRI=image->mri;
 for(i=0;i<wdth;i++)
  for(j=0;j<hght;j++)
   for(k=0;k<dpth;k++)
   {
    val1=imageMRI[i][j][k];
    if (val1>0)
    {
     nbvoxelsOverlap+=1.0;
     if (val1>max1) max1=val1;
    }
   }
   
 //cas ou toute l'image est a nul
 if (nbvoxelsOverlap<1.0) { if (max_image!=NULL) (*max_image)=0; return 0; }
//---
// if (max_image!=NULL) (*max_image)=(int)max1;
// i_nbcl=(int)floor(pow(wdth*hght*dpth,1.0/3.0));
// if (i_nbcl>max1) i_nbcl=(max1+1);  return i_nbcl;
//--- 
 for(i=0;i<wdth;i++)
  for(j=0;j<hght;j++)
   for(k=0;k<dpth;k++)
   {
    val=(double)imageMRI[i][j][k];
    if (val>0.0)
    {
     moy+=val/nbvoxelsOverlap;
     moy2+=val*val/nbvoxelsOverlap;
    }
   }

 sigma=sqrt(moy2-moy*moy);
 w_bin=3.49*sigma*pow(nbvoxelsOverlap,-1.0/3.0);
 if (w_bin>eps) i_nbcl=(int)floor(max1/w_bin+0.5); else i_nbcl=0;

 if (i_nbcl>(int)max1) i_nbcl=(int)(max1+1);
 if (max_image!=NULL) (*max_image)=(int)max1;

 return i_nbcl;
}

/*******************************************************************************
**    cr_VP_histos(Imreca, Imref, minimisation_par_gradient)
*/
/*!
**    allocation memoire de la structure devant contenir les donnees necessaires
**    aux calculs lies au volume partiel
**
**  \param imreca : image a recaler
**  \param Imref : image de reference
**  \param minimisation_par_gradient : indique si on utilise une methode de minimisation
**         uilisant le gradient
**  \retval : la structure allouee
*******************************************************************************/
VP_histos * cr_VP_histos(grphic3d *Imreca, grphic3d *Imref,bool minimisation_par_gradient)
{
 int i,j,k,l;
 int nbcl1, nbcl2, max1, max2;
 VP_histos *donnees_VP;

 //calcul des maxpixel dans les deux images
 //et du nombre de classes de valeurs de voxels qu'on va utiliser
 nbcl1=imx_calculer_nombre_de_classes_3d_p(Imreca, &max1);
 nbcl2=imx_calculer_nombre_de_classes_3d_p(Imref, &max2);

 //allocations memoire
 donnees_VP = malloc(sizeof(VP_histos));
 donnees_VP->histo_marginal1 = CALLOC(nbcl1, double);
 donnees_VP->histo_marginal2 = CALLOC(nbcl2, double);
 donnees_VP->histo_conjoint = alloc_dmatrix(nbcl1,nbcl2);
 for (i=0;i<nbcl1;i++)
  for (j=0;j<nbcl2;j++)
   donnees_VP->histo_conjoint[i][j]=0.0;

 if (minimisation_par_gradient)
 {
  donnees_VP->histos_gradient = CALLOC (3, double ***);
  for (i=0;i<3;i++)
   donnees_VP->histos_gradient[i] = CALLOC (4, double **);
  for (i=0;i<3;i++)
   for (j=0;j<4;j++)
    donnees_VP->histos_gradient[i][j] = alloc_dmatrix(nbcl1,nbcl2);
   for (i=0;i<3;i++)
    for (j=0;j<4;j++)
     for (k=0;k<nbcl1;k++)
      for (l=0;l<nbcl2;l++)
       donnees_VP->histos_gradient[i][j][k][l]=0.0;
 }
 else
 { donnees_VP->histos_gradient = NULL; }

 //mises a jour
 donnees_VP->max1 = max1;
 donnees_VP->max2 = max2;
 donnees_VP->nb_cl1 = nbcl1;
 donnees_VP->nb_cl2 = nbcl2;
 donnees_VP->minimisation_par_gradient = minimisation_par_gradient;
 donnees_VP->entropie1 = 0.0;
 donnees_VP->entropie2 = 0.0;
 donnees_VP->entropie_conjointe = 0.0;
 donnees_VP->nb_voxels_recouvrement = 0;

 return donnees_VP;
 
}

/*******************************************************************************
**    free_VP_histos(donnees_VP_ptr)
*/
/*!
**    liberation memoire de la structure devant contenir les donnees necessaires
**    aux calculs lies au volume partiel
**
**  \param donnees_VP_ptr : pointeur sur la structure a liberer
*******************************************************************************/
void free_VP_histos(VP_histos **donnees_VP_ptr)
{
 int i, j, nbcl1, nbcl2;

 nbcl1 = (*donnees_VP_ptr)->nb_cl1;
 nbcl2 = (*donnees_VP_ptr)->nb_cl2;

 //liberations memoire
 FREE((*donnees_VP_ptr)->histo_marginal1);
 FREE((*donnees_VP_ptr)->histo_marginal2);
 free_dmatrix((*donnees_VP_ptr)->histo_conjoint,nbcl1,nbcl2);

 if ((*donnees_VP_ptr)->minimisation_par_gradient)
 {
  for (i=0;i<3;i++)
   for (j=0;j<4;j++)
    free_dmatrix((*donnees_VP_ptr)->histos_gradient[i][j],nbcl1,nbcl2);
  for (i=0;i<3;i++)
   FREE((*donnees_VP_ptr)->histos_gradient[i]);
  FREE((*donnees_VP_ptr)->histos_gradient);
 }

 FREE(*donnees_VP_ptr);

 return;
}

/*******************************************************************************
**        imx_query_erreur_type_3D                                            
*/
/*!                                                                        
**       \brief choix de la fonction de cout                                    
**
**      ouvre une boite de dialogue pour le choix de la fonction de cout
**      \retval int : le numero de la fonction de cout
**          - Quadratict -> 0
**          - Quadratic Robust -> 1
**          - woods -> 2
**          - woods robust -> 3
**          - etc ...
**          
*******************************************************************************/
int imx_query_erreur_type_3D()
{
    char *query[100];
    query[0]="Quadratic";
    query[1]="Quadratic Robust";
    query[2]="woods";
    query[3]="woods robust";
    query[4]="IM Old";
    query[5]="IM";
    query[6]="IM Normalisee Studholme";
    query[7]="IM Normalisee Maes1";
    query[8]="IM Normalisee Maes2";
    query[9]="Entropie conjointe";
    query[10]="Correlation ratio";
    query[11]="IM (VP)";
    query[12]="IM Normalisee Studholme (VP)";
    query[13]="IM Normalisee Maes1 (VP)";
    query[14]="IM Normalisee Maes2 (VP)";
    query[15]="Entropie conjointe (VP)";
    query[16]="vrai critere de Woods";
    query[17]="critere de Woods robuste 2";
    query[18]="Quadratic Robust 2";
    query[19]="IM VP stochastique";
    query[20]="correlation ratio robuste";
    query[21]="IM BSpline";
    query[22]="IM Regionelle";
    query[23]="coefficient de correlation";
    query[24]="Plan de symetrie";
    query[25]="\0";
    query[26]=NULL;
    return GETV_QCM("Distance",(char **)query);
}

/*------------------------------------------------*/
/*!
**  \brief determine la fonction d'erreur
**
**  \param  dist_type : le type de distance choisie
**  \retval dist_func_t : un pointeur sur la fonction d'erreur
*/
/*------------------------------------------------*/

dist_func_t
imx_choose_distance_fct(int dist_type)
{
  dist_func_t distance;
  switch (dist_type)
  {
    case 0: distance=erreur_quad_3d;aff_log("QUAD ");break;
    case 1: distance=erreur_quad_robust_geman_3d;aff_log("QUAD ROBUST ");break;
    case 2: distance=erreur_woods_3d;aff_log("WOODS ");break;
    case 3: distance=erreur_woods_robust_3d;aff_log("WOODS ROBUST ");break;
    case 4: distance=erreur_IM_3d;aff_log("IM ");break;
    case 5: distance=erreur_IM_simple_3d;aff_log("IM SIMPLE ");break;
    case 6: distance=erreur_IM_Normalisee_Studholme_simple_3d;aff_log("IM NORM STUDHOLME ");break;
    case 7: distance=erreur_IM_Normalisee_Maes1_simple_3d;aff_log("IM NORM MAES1 ");break;
    case 8: distance=erreur_IM_Normalisee_Maes2_simple_3d;aff_log("IM NORM MAES2 ");break;
    case 9: distance=erreur_entropie_conjointe_simple_3d;aff_log("Entropie conjointe simple ");break;
    case 10: distance=erreur_CR_3d;aff_log("Correlation ratio ");break;
    case 11: distance=erreur_IM_VP_3d;aff_log("IM (VP) ");break;

    case 12: distance=erreur_IMNS_VP_3d;aff_log("IM NORM STUDHOLME (VP) ");break;
    case 13: distance=erreur_IMNM1_VP_3d;aff_log("IM NORM MAES1 (VP) ");break;
    case 14: distance=erreur_IMNM2_VP_3d;aff_log("IM NORM MAES2 (VP) ");break;
    case 15: distance=erreur_entropie_conjointe_VP_3d;aff_log("Entropie conjointe (VP) ");break;
    case 16: distance=erreur_Woods_vraie_3d;aff_log("vrai critere de Woods ");break;
    case 17: distance=erreur_woods_robust2_3d;aff_log("critere de Woods robuste 2 ");break;
    case 18: distance=erreur_quad_robust2_3d;aff_log("Quadratic Robust 2 ");break;
    case 19: distance=erreur_IM_VP_stochastique_3d;aff_log("IM VP stochastique ");break;
    case 20: distance=erreur_CR_robust_3d;aff_log("correlation ratio robuste ");break;
    case 21: distance=erreur_IM_BSpline_3d;aff_log("IM BSpline ");break;
    case 22: distance=erreur_IM_Region_3d;aff_log("IM Regionelle ");break;
    case 23: distance=erreur_coeff_correlation_d;aff_log("coefficient de correlation ");break;
    case 24: distance=erreur_quad_plandesymetrie_3d;aff_log("Plan de symetrie ");break;
    default:distance=erreur_quad_3d;aff_log("QUAD ");break;
  }
  return distance;
}

double imx_Erreur_3d(int im_src, int im_dst, int err_type)
{
    grphic3d *imsrc,*imdst;
    double res;
    dist_func_t  func_erreur;
    ptr_distance_param dist_param=CALLOC(1, distance_param);

    imsrc = ptr_img_3d(im_src);
    imdst = ptr_img_3d(im_dst);
    func_erreur = imx_choose_distance_fct(err_type);

    dist_param->imreca=imsrc;
    dist_param->imref=imdst;
    res = func_erreur(dist_param);
    FREE(dist_param);
    return res;
}

void Erreur_3d()
{
    int err_type,im_src,im_dst;
    double res;
    im_src = GET_WND3D("image depart");
    im_dst = GET_WND3D("image resultat");
    err_type = imx_query_erreur_type_3D();
    
    res = imx_Erreur_3d(im_src,im_dst,err_type);
    printf("ERREUR = %f\n",res);
}

bool dist_utilisant_VP(dist_func_t distance)
{
  if (!distance) return FALSE;

  if ((distance==erreur_IM_VP_3d)
    ||(distance==erreur_IMNS_VP_3d)
    ||(distance==erreur_IMNM1_VP_3d)
    ||(distance==erreur_IMNM2_VP_3d)
    ||(distance==erreur_entropie_conjointe_VP_3d)
    ||(distance==erreur_IM_VP_stochastique_3d)) return TRUE;

 return FALSE;
}

double calc_info_mutuelle_histos(IM3DTYPE im_type, double **histo_conjoint, double *histo1, double *histo2, int nbcl1, int nbcl2, double nbVoxelsOverlap, double nbVoxelsHorsOverlap, double *EntropieConjointe, double *Entropie1, double *Entropie2)
{
 double info_mutuelle=DBL_MAX,entropie_conjointe=DBL_MAX, entropie1=DBL_MAX, entropie2=DBL_MAX;

 //calcul des histogrammes marginaux
 calcul_histos_marginaux(histo_conjoint, histo1, histo2, nbcl1, nbcl2,nbVoxelsHorsOverlap);

 //calcul des entropies
 entropie_conjointe = calcul_entropie_conjointe(histo_conjoint, nbcl1, nbcl2, nbVoxelsOverlap);
 entropie1 = calcul_entropie(histo1, nbcl1, nbVoxelsOverlap);
 entropie2 = calcul_entropie(histo2, nbcl2, nbVoxelsOverlap);

 //calcul de la distance
 switch(im_type)
 {
  case ENTROPIE_CONJOINTE: info_mutuelle=100.0*entropie_conjointe; break;
  case IM: info_mutuelle=calcul_IM(entropie_conjointe, entropie1, entropie2); info_mutuelle=-1000.0*info_mutuelle; break;
  case IMNS: info_mutuelle=calcul_IMNS(entropie_conjointe, entropie1, entropie2); info_mutuelle=-100.0*info_mutuelle; break;
  case IMNM1: info_mutuelle=calcul_IMNM1(entropie_conjointe, entropie1, entropie2); info_mutuelle=-1000.0*info_mutuelle; break;
  case IMNM2: info_mutuelle=calcul_IMNM2(entropie_conjointe, entropie1, entropie2); info_mutuelle=100.0*info_mutuelle; break;
  default: info_mutuelle=100.0*entropie_conjointe; break;
 }

 *EntropieConjointe=entropie_conjointe;
 *Entropie1=entropie1;
 *Entropie2=entropie2;

 return info_mutuelle;
}


#define N_R_TO_TAILLE_VEC_IM_REGION(N_R) (2*N_R+1)*(2*N_R+1)*(2*N_R+1)

double calcul_info_mutuelle_region_3d(IM3DTYPE im_type, ptr_distance_param dist_param)
{
 grphic3d *imreca=dist_param->imreca;
 grphic3d *imref= dist_param->imref;
 TYPEMRI3D ***imrecaMRI=NULL, ***imrefMRI=NULL, ***imrecaQuantifMRI=NULL, ***imrefQuantifMRI=NULL;
 const int N_R=1;
 double **covar_mat=NULL;
 int taille_vec=N_R_TO_TAILLE_VEC_IM_REGION(N_R);
 int i, j, k;
 int res=0;
 double entropie_conjointe=0.0, entropie1=0.0, entropie2=0.0, info_mutuelle=0.0;
 int nbcl1,nbcl2,max1,max2;
 grphic3d *imrecaQuantif=NULL, *imrefQuantif=NULL;
 int width, height, depth;

 imrecaMRI=imreca->mri; imrefMRI=imref->mri;
 width=imreca->width; height=imreca->height; depth=imreca->depth;
 
 //calcul des maxpixel dans les deux images
 //et du nombre de classes de valeurs de voxels qu'on va utiliser
// nbcl1=imx_calculer_nombre_de_classes_3d_p(imreca, &max1);
// nbcl2=imx_calculer_nombre_de_classes_3d_p(imref, &max2);
 nbcl1=dist_param->nbclReca; nbcl2=dist_param->nbclRef;
 max1=dist_param->maxReca; max2=dist_param->maxRef;

 imrecaQuantif=cr_grphic3d(imreca); imrefQuantif=cr_grphic3d(imref);
 imrecaQuantifMRI=imrecaQuantif->mri; imrefQuantifMRI=imrefQuantif->mri;

 nbcl1=nbcl1-1; max1=max1+1; nbcl2=nbcl2-1; max2=max2+1;
 for (i=0;i<width;i++)
  for (j=0;j<height;j++)
   for (k=0;k<depth;k++)
   {
    if (imrecaMRI[i][j][k]<=0) imrecaQuantifMRI[i][j][k]=0;
    else imrecaQuantifMRI[i][j][k]=(nbcl1*imrecaMRI[i][j][k])/max1+1;
    if (imrefMRI[i][j][k]<=0) imrecaQuantif->mri[i][j][k]=0;
    else imrefQuantifMRI[i][j][k]=(nbcl2*imrefMRI[i][j][k])/max2+1;
   }

 dist_param->imreca=imrecaQuantif;  dist_param->imref=imrefQuantif;
 //allocations memoire
 covar_mat=alloc_dmatrix(2*taille_vec, 2*taille_vec);
 for (i=0;i<(2*taille_vec);i++)
   { for (j=0;j<(2*taille_vec);j++) covar_mat[i][j]=0.0; }

 //calcul de la matrice de covariance
 res=calcul_covar_region_3d(dist_param, N_R, covar_mat);

 //calcul des entropies  
 res=calcul_entropie_covar_region_3d(N_R, covar_mat, &entropie_conjointe, &entropie1, &entropie2);

 //calcul de la distance
 switch(im_type)
 {
  case ENTROPIE_CONJOINTE: info_mutuelle=100.0*entropie_conjointe; break;
  case IM: info_mutuelle=calcul_IM(entropie_conjointe, entropie1, entropie2); info_mutuelle=-1000.0*info_mutuelle; break;
  case IMNS: info_mutuelle=calcul_IMNS(entropie_conjointe, entropie1, entropie2); info_mutuelle=-100.0*info_mutuelle; break;
  case IMNM1: info_mutuelle=calcul_IMNM1(entropie_conjointe, entropie1, entropie2); info_mutuelle=-1000.0*info_mutuelle; break;
  case IMNM2: info_mutuelle=calcul_IMNM2(entropie_conjointe, entropie1, entropie2); info_mutuelle=100.0*info_mutuelle; break;
  default: info_mutuelle=100.0*entropie_conjointe; break;
 }

 dist_param->imreca=imreca;  dist_param->imref=imref;
  
 //liberation memoire
 free_dmatrix(covar_mat, 2*taille_vec, 2*taille_vec); free_grphic3d(imrecaQuantif); free_grphic3d(imrefQuantif);

 return info_mutuelle;
}

int calcul_covar_region_3d(ptr_distance_param dist_param, const int N_R, double **covar_mat)
{
 grphic3d *imreca=dist_param->imreca;
 grphic3d *imref= dist_param->imref;
 TYPEMRI3D ***imrecaMRI=NULL, ***imrefMRI=NULL, ***maskrecaMRI=NULL;
 int i, j, k;
 int l, m, n;
 int index_vec=0;
 double *vec=NULL, *moy_vec=NULL;
 int taille_vec=N_R_TO_TAILLE_VEC_IM_REGION(N_R);
 double nb_count=0.0;
 int width, height, depth;
 double valReca=0.0, valRef=0.0;
 bool negVal=FALSE;

 imrecaMRI=imreca->mri; imrefMRI=imref->mri;
 width=imreca->width; height=imreca->height; depth=imreca->depth;

 maskrecaMRI=(ptr_mask_3d_p(imreca))->mri; 
 
 moy_vec=CALLOC(2*taille_vec, double);
 vec=CALLOC(2*taille_vec, double);

 nb_count=0.0;

 //calcul de la zone de recouvrement
 for (i=N_R; i<(width-N_R); i++)
 {
  for (j=N_R; j<(height-N_R); j++)
  {
   for (k=N_R; k<(depth-N_R); k++)
   {

    negVal=FALSE;
    
    for (l=-N_R;l<=N_R;l++)
    {
     for (m=-N_R;m<=N_R;m++)
     {
      for (n=-N_R;n<=N_R;n++)
      {
        index_vec=(l+N_R)*(2*N_R+1)*(2*N_R+1)+(m+N_R)*(2*N_R+1)+(n+N_R);
        if (imrecaMRI[i+l][j+m][k+n]<=0) { negVal=TRUE; break; }
        if (imrefMRI[i+l][j+m][k+n]<=0) { negVal=TRUE; break; }
      }
     }
    }

    if (!negVal) { nb_count+=1.0; maskrecaMRI[i][j][k]=1; } else { maskrecaMRI[i][j][k]=0; }
    
   }
  }
 }
 
 //calcul sur la zone de recouvrement
 for (i=N_R; i<(width-N_R); i++)
 {
  for (j=N_R; j<(height-N_R); j++)
  {
   for (k=N_R; k<(depth-N_R); k++)
   {

    if (maskrecaMRI[i][j][k]==1)
    {
      
     //remplissage du vecteur de voisinage centre sur la moyenne
     for (l=-N_R;l<=N_R;l++)
     {
      for (m=-N_R;m<=N_R;m++)
      {
       for (n=-N_R;n<=N_R;n++)
       {
        index_vec=(l+N_R)*(2*N_R+1)*(2*N_R+1)+(m+N_R)*(2*N_R+1)+(n+N_R);
        valReca=(double)imrecaMRI[i+l][j+m][k+n];
        vec[index_vec]=valReca; moy_vec[index_vec]+=valReca/nb_count;
        index_vec+=taille_vec;
        valRef=(double)imrefMRI[i+l][j+m][k+n];
        vec[index_vec]=valRef; moy_vec[index_vec]+=valRef/nb_count;
       }
      }
     }

     //mise a jour de la matrice de covariance
     for (l=0;l<(2*taille_vec);l++)
     {
      for (m=0;m<(2*taille_vec);m++)
      {
       covar_mat[l][m]+=vec[l]*vec[m]/nb_count;
      }
     }
    
    }
     
   }
  }
 }

 //normalisation
 for (l=0;l<(2*taille_vec);l++)
 {
  for (m=0;m<(2*taille_vec);m++)
  {
   covar_mat[l][m]=covar_mat[l][m]-moy_vec[l]*moy_vec[m];
  }
 }

 FREE(vec); FREE(moy_vec); 

 return 0; 
}

int calcul_entropie_covar_region_3d(const int N_R, double **covar_mat, double *entropie_conjointe, double *entropie1, double *entropie2)
{
 int i,j;
 double **sous_matrice=NULL;
 int taille_vec=N_R_TO_TAILLE_VEC_IM_REGION(N_R);
 double **covar_mat2=NULL;
 double *diag=NULL;
 double det_conjoint, det1, det2;
 int *index=NULL, *index2=NULL;
 int err=0;

 sous_matrice=alloc_dmatrix(taille_vec, taille_vec);
 covar_mat2=alloc_dmatrix(2*taille_vec, 2*taille_vec);
 index=CALLOC(2*taille_vec, int);
 index2=CALLOC(taille_vec, int);
 diag=CALLOC(2*taille_vec, double);
 
 //calcul du determinant de la matrice de covariance
 for (i=0;i<(2*taille_vec);i++) { for (j=0;j<(2*taille_vec);j++) covar_mat2[i][j]=covar_mat[i][j]; }
 err=choldc(covar_mat2, diag, 2*taille_vec);
 if (err)
 {
   fprintf (stderr, "ERROR\n");
 }
 det_conjoint=1.0;
 for (i=0;i<(2*taille_vec);i++) det_conjoint*=diag[i]*diag[i];
 
 //calcul du determinant de la sous matrice 1
 for (i=0;i<taille_vec;i++) { for (j=0;j<taille_vec;j++) sous_matrice[i][j]=covar_mat[i][j]; }
 err=choldc(sous_matrice, diag, taille_vec);
 if (err)
 {
   fprintf (stderr, "ERROR\n");
 }
 det1=1.0;
 for (i=0;i<taille_vec;i++) det1*=diag[i]*diag[i];
 
 //calcul du determinant de la sous matrice 2
 for (i=0;i<taille_vec;i++) { for (j=0;j<taille_vec;j++) sous_matrice[i][j]=covar_mat[i+taille_vec][j+taille_vec]; }
 err=choldc(sous_matrice, diag, taille_vec);
 if (err)
 {
   fprintf (stderr, "ERROR\n");
 }
 det2=1.0;
 for (i=0;i<taille_vec;i++) det2*=diag[i]*diag[i];
  
 //calcul des entropies
 *entropie_conjointe=log(pow(2.0*PI*exp(1.0), taille_vec)*sqrt(det_conjoint));
 *entropie1=log(pow(2.0*PI*exp(1.0), taille_vec/2.0)*sqrt(det1));
 *entropie2=log(pow(2.0*PI*exp(1.0), taille_vec/2.0)*sqrt(det2));

 //liberation memoire
 free_dmatrix(sous_matrice, taille_vec, taille_vec);
 free_dmatrix(covar_mat2, 2*taille_vec, 2*taille_vec);
 FREE(index); FREE(index2); FREE(diag);
 
 return 0;
}


//---calcul de l'image edgeness---//
void calcul_edgeness_3d()
{
 int im_deb,im_res;
 int taille_voisinage;

 //demande des images
 im_deb=GET_PLACE3D("image a traiter");
 im_res=GET_PLACE3D("image resultat");

 taille_voisinage=GET_INT("taille du voisinage a considerer", 1, &err);

 im_calcul_edgeness_3d(im_deb, im_res, taille_voisinage);
 
}

int im_calcul_edgeness_3d(int im_deb, int im_res, int taille_voisinage)
{
  grphic3d *imdeb,*imres;
  int res;
  
  imdeb=ptr_img_3d(im_deb);
  imres=ptr_img_3d(im_res);

  res=im_calcul_edgeness_3d_p(imdeb, imres, taille_voisinage);

  return res;
}

int im_calcul_edgeness_3d_p(grphic3d *imdeb, grphic3d *imres, int taille_voisinage)
{
 int i,j,k;
 int l,m,n;
 int min_w, max_w, min_h, max_h, min_d, max_d;
 int width, height, depth;
 int val_cour;
 double sum;
 TYPEMRI3D ***imdebMRI=NULL, ***imresMRI=NULL;
 double ***tempIm=NULL;
 
 width=(int)imdeb->width; height=(int)imdeb->height; depth=(int)imdeb->depth;
 imx_copie_3d_p(imdeb, imres); 

 tempIm=alloc_dmatrix_3d(width, height, depth);
 
 imdebMRI=imdeb->mri; imresMRI=imres->mri;
 
 for (i=0;i<width;i++)
 {
  for (j=0;j<height;j++)
  {
   for (k=0;k<depth;k++)
   {
    min_w=MAXI(0, i-taille_voisinage);
    max_w=MINI(width-1, i+taille_voisinage);
    min_h=MAXI(0, j-taille_voisinage);
    max_h=MINI(height-1, j+taille_voisinage);
    min_d=MAXI(0, k-taille_voisinage);
    max_d=MINI(depth-1, k+taille_voisinage);

    val_cour=(int)imdebMRI[i][j][k];
    sum=0.0;
    
    for (l=min_w;l<=max_w;l++)
    {
     for (m=min_h;m<=max_h;m++)
     {
      for (n=min_d;n<=max_d;n++)
      {
       sum+=abs(val_cour-imdebMRI[l][m][n]);
      }
     }
    }

    tempIm[i][j][k]=sum;
           
   } 
  }
 }

 normaliser_3d(imdeb, imres, tempIm, width, height, depth);
 
 free_dmatrix_3d(tempIm);

 return 0;
}

/*! \ingroup     DistanceFonction  @{ */
/*******************************************************************************
**     erreur_ICP_3d(im1,im2)                                        
*/                                                                    
/*!    distance ICP
        \param im1 : carte de distance 
                     im2 : image de r�f�rence binaire 
        \retval ICP
*******************************************************************************/
double erreur_ICP_3d(ptr_distance_param dist_param)
{
 unsigned int i,j,k,Ntot;
 TDimension wdth,hght,dpth;
 double r=DBL_MAX,diff;
 grphic3d *im1=NULL, *im2=NULL;
 TYPEMRI3D ***im1MRI=NULL, ***im2MRI=NULL;
 ERREUR_RECALAGE *err=NULL;

 im1=dist_param->imreca; im2=dist_param->imref;
 im1MRI=im1->mri; im2MRI=im2->mri;
 err=dist_param->err;

 wdth=im1->width;hght=im1->height;dpth=im1->depth;
 if ((wdth!=im2->width)||(hght!=im2->height)||(dpth!=im2->depth))
 {
  RAISE_ERR_TEXT(err, PB_TAILLES_IMAGES, "tailles d'images differentes dans erreur_quad_3d\n");
 }

 r=0.0;
 Ntot=0;
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
    if (im2MRI[i][j][k]>0)
        {
     diff=(double)(im1MRI[i][j][k]);
                 
     r=r+diff;
         Ntot++;
    }
 r=1.0*r/Ntot;

  if (err!=NULL) *err=NO_ERR;
end_func:

 return(r);
}
/*! @} */

/*! \ingroup     DistanceFonction  @{ */
/*******************************************************************************
**     erreur_ICP_3d(im1,im2)                                        
*/                                                                    
/*!    distance ICP
        \param im1 : carte de distance 
                     im2 : image de r�f�rence binaire 
        \retval ICP
*******************************************************************************/
double erreur_ICP_sym_3d(ptr_distance_param dist_param)
{
 unsigned int i,j,k,Ntot;
 TDimension wdth,hght,dpth;
 double r=DBL_MAX,r1=DBL_MAX,r2=DBL_MAX,diff,vmin=0,maskvmin=0,J;
 grphic3d *im1=NULL, *im2=NULL;
 TYPEMRI3D ***im1MRI=NULL, ***im2MRI=NULL;
 ERREUR_RECALAGE *err=NULL;
 double *param;
 double seuil;
 
 param=(double*)malloc(15*sizeof(double));

 im1=dist_param->imreca; im2=dist_param->imref;
 im1MRI=im1->mri; im2MRI=im2->mri;
 err=dist_param->err;

 maskvmin=1.0/im1->mask->rcoeff;
 vmin=1.0/im1->rcoeff;

seuil = 150*im1->rcoeff;

 wdth=im1->width;hght=im1->height;dpth=im1->depth;
 if ((wdth!=im2->width)||(hght!=im2->height)||(dpth!=im2->depth))
 {
  RAISE_ERR_TEXT(err, PB_TAILLES_IMAGES, "tailles d'images differentes dans erreur_quad_3d\n");
 }


/*for (i=0;i<dist_param->transfres->nb_param;i++)
    param[i]=dist_param->transfres->param[i];

 	J=1;

if (dist_param->transfres->typetrans==RIGIDZOOM3D)
    rigidz_to_affine_3d(param);
        
if (dist_param->transfres->typetrans==AFFINEDECOUPLE3D)
    affine_decouple_to_affine_3d(param);
        
if ((dist_param->transfres->typetrans==AFFINE3D)||(dist_param->transfres->typetrans==RIGIDZOOM3D)||(dist_param->transfres->typetrans==AFFINEDECOUPLE3D))
  {
     double a11,a12,a13,a21,a22,a23,a31,a32,a33;
    a11=param[0];a12=param[1];a13=param[2];
  a21=param[3];a22=param[4];a23=param[5];
  a31=param[6];a32=param[7];a33=param[8];
    
    J=fabs(a11*(a22*a33-a32*a23)-a21*(a12*a33-a13*a32)+a31*(a12*a23-a22*a13));
    }
    
free(param);    
    
if ((J>8)||(J<0.125))
    return(HUGE_VAL);

*/
 r1=0.0;
 Ntot=0;
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
    if (im2MRI[i][j][k]>0)
        {
     if (im1MRI[i][j][k]<vmin) 
          diff=(double)(im1->max_pixel);//return(HUGE_VAL);
         else
          diff=(double)(im1MRI[i][j][k]-vmin);
                 
     r1=r1+im2MRI[i][j][k]*diff*diff;
         
                Ntot++;
    }
        
        /*if (Ntot>0)
        r1=1.0*r1/Ntot;
        else 
        r1=0.0;*/
        
        r2=0.0;

        Ntot=0;
        for (i=0;i<wdth;i++)
    for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
            if ((im2->mask->mri[i][j][k]>0)&&(im1->mask->mri[i][j][k]>=maskvmin)&&(im1MRI[i][j][k]>=vmin))
            {   
            
            diff=(double)((im1->mask->mri[i][j][k]-maskvmin));
            
            r2=r2+diff*diff; // /(im1MRI[i][j][k]-vmin+1.0);
            
            Ntot++;
        }

        /*if (Ntot>0)
        r2=1.0*r2/Ntot;
        else
        r2=0.0;*/
        
 r=sqrt(1.0*r1*im1->rcoeff*im1->rcoeff+1.0*r2*im1->mask->rcoeff*im1->mask->rcoeff);

//printf("Eblanc %f  \t Enoir %f \n",r1*im1->rcoeff*im1->rcoeff,r2*im1->mask->rcoeff*im1->mask->rcoeff);

  if (err!=NULL) *err=NO_ERR;
end_func:

 return(r);
}
/*! @} */

/*! \ingroup     DistanceFonction  @{ */
/*******************************************************************************
**     erreur_quad_plandesymetrie_3d(im1,im2)                                        
*/                                                                    
/*!    Fonction de cout pour determiner le plan de symetrie d'une image
*******************************************************************************/
double erreur_quad_plandesymetrie_3d(ptr_distance_param dist_param)
{
 unsigned int i,j,k, tot;
 TDimension wdth,hght,dpth;
 double r=DBL_MAX,diff;//,m1,m2;
 grphic3d *im1=NULL, *im2=NULL, *im1flipped;
 TYPEMRI3D ***im1MRI=NULL, ***im2MRI=NULL;
 ERREUR_RECALAGE *err=NULL;



 im1=dist_param->imreca; im2=dist_param->imref;
 im1MRI=im1->mri; im2MRI=im2->mri;
 err=dist_param->err;

 im1flipped=cr_grphic3d(im1);



 wdth=im1->width;hght=im1->height;dpth=im1->depth;
 if ((wdth!=im2->width)||(hght!=im2->height)||(dpth!=im2->depth))
 {
  RAISE_ERR_TEXT(err, PB_TAILLES_IMAGES, "tailles d'images differentes dans erreur_quad_3d\n");
 }



Miroir_ttaxe_3d_p(im1,im1flipped,1);

 r=0.0;
 tot=0;
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
    {
    if ((im1flipped->mri[i][j][k]==0)&( im1->mri[i][j][k] == 0))
        tot++;
    
     diff=(double)(im1flipped->mri[i][j][k]-im1->mri[i][j][k]);
     r=r+diff*diff;
     r=r-2*im1->mri[i][j][k]*im1->mri[i][j][k];
    }
 //r=sqrt(r);


  if (err!=NULL) *err=NO_ERR;
end_func:

//liberation memoire  
 if (im1flipped) free_grphic3d(im1flipped);

return(r);
}
/*! @} */
