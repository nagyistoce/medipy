/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*******************************************************************************/
/*!	   \file:		detection_anomalies_3d.c
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
#include "noyau/imagix.h"
#include "noyau/io/imx_file.h"
#include "noyau/imx_lang.h"
#include "noyau/io/imx_log.h"
#include "noyau/gui/imx_picture3d.h"
//#include "recalage/chps_3d.h"
#include "recalage/distance_3d.h"
#include "traitement/norm_3d.h"
#include "detection_anomalies_3d.h"

#define SEUIL_Z_MAP 1.0

#ifndef NB_ITER_MAX_LINEMIN_GENERAL_3D
#define NB_ITER_MAX_LINEMIN_GENERAL_3D 30 /*nombre d'iterations max pour la minimisation lineaire*/
#endif
#ifndef CGOLD
#define CGOLD 0.3819660
#endif
#ifndef ZEPS
#define ZEPS 1.0e-10 /*definition de epsilon*/
#endif
#ifndef SGN
#define SGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#endif

//////definitions de types locaux//////

typedef struct s_param_min_TLS
{
 double covariance;
 double varianceRef;
 double varianceSrc;
} param_min_TLS;

double fonction_energie_TLS(double x, param_min_TLS *params_image);
double linemin_general_3d(lin_energie_func_t fonction_energie, void *params_additionnels, double borne_min, double borne_max, double *abscisse_trouve);

//-----------------------------------------------------------------------------//
//                             METHODE Z MAP                                   //
//-----------------------------------------------------------------------------//

/*******************************************************************************
**    selection_Z_Map_3d()
*/
/*!
**    a partir de 2 images spect recalees et segmentees construit des masques englobant
**    la region de reference choisie a partir d'une methode de construction par z Map:
**    "Automated Detection of Local Normalization Areas for Ictal-Interictal Subtraction Brain SPECT"
**
*******************************************************************************/
void selection_Z_Map_3d()
{
 int im_1, im_2;

 im_1=GET_PLACE3D("image 1");
 im_2=GET_PLACE3D("image 2");

 imx_selection_Z_Map_3d(im_1, im_2);

 show_picture_3d(im_1); show_picture_3d(im_2);
}

/*******************************************************************************
**    imx_selection_Z_Map_3d(im_1, im_2)
*/
/*!
**    a partir de 2 images spect recalees et segmentees construit des masques englobant
**    la region de reference choisie a partir d'une methode de construction par z Map:
**    "Automated Detection of Local Normalization Areas for Ictal-Interictal Subtraction Brain SPECT"
**
**  \param im_1, im_2: les spect recales
**	\retval 
*******************************************************************************/
void imx_selection_Z_Map_3d(int im_1, int im_2)
{
 grphic3d *im1, *im2;

 im1=ptr_img_3d(im_1);
 im2=ptr_img_3d(im_2);

 imx_selection_Z_Map_3d_p(im1, im2);
}

/*******************************************************************************
**    imx_selection_Z_Map_3d_p(im1, im2)
*/
/*!
**    a partir de 2 images spect recalees et segmentees construit des masques englobant
**    la region de reference choisie a partir d'une methode de construction par z Map:
**    "Automated Detection of Local Normalization Areas for Ictal-Interictal Subtraction Brain SPECT"
**
**  \param im1, im2: les spect recales
**	\retval : 0 en cas de succes
*******************************************************************************/
int imx_selection_Z_Map_3d_p(grphic3d *im1, grphic3d *im2)
{
 unsigned int i,j,k;
 TDimension wdth,hght,dpth;
 grphic3d  *mask1=NULL, *mask2=NULL;
 double nbvoxelsOverlap1=0.0, nbvoxelsOverlap2=0.0, sigma1=0.0, sigma2=0.0;
 double moy1=0.0, moy2=0.0, var1=0.0, var2=0.0, val1=0.0, val2=0.0, z1=0.0, z2=0.0;

 wdth=im1->width; hght=im1->height; dpth=im1->depth;
 if ((wdth!=im2->width)||(hght!=im2->height)||(dpth!=im2->depth))
 {
  fprintf(stderr,"les dimensions des images ne concordent pas dans imx_selection_Z_Map_3d_p\n");
  return 1;
 }

 //creation des masques
 mask1=cr_mask_3d_ptr(im1);
 mask2=cr_mask_3d_ptr(im2);

 if ((!mask1)||(!mask2))
 {
  fprintf(stderr,"l'allocation des masques a echoue dans imx_selection_Z_Map_3d_p\n");
  return 2;
 }

 //calcul du volume des zones segementees
 for (i=0;i<wdth;i++)
 {
  for (j=0;j<hght;j++)
  {
   for (k=0;k<dpth;k++)
   {

    val1=(double)im1->mri[i][j][k];
    if (val1>0.0)
    {
     nbvoxelsOverlap1+=1.0;
    }

    val2=(double)im2->mri[i][j][k];
    if (val2>0.0)
    {
     nbvoxelsOverlap2+=1.0;
    }

   }
  }
 }

 //calcul de criteres statistiques sur les images
 for (i=0;i<wdth;i++)
 {
  for (j=0;j<hght;j++)
  {
   for (k=0;k<dpth;k++)
   {

    val1=(double)im1->mri[i][j][k];
    if (val1>0.0)
    {
     moy1+=val1/nbvoxelsOverlap1;
     var1+=val1*val1/nbvoxelsOverlap1;
    }

    val2=(double)im2->mri[i][j][k];
    if (val2>0.0)
    {
     moy2+=val2/nbvoxelsOverlap2;
     var2+=val2*val2/nbvoxelsOverlap2;
    }

   }
  }
 }
 sigma1=sqrt(var1-moy1*moy1);
 sigma2=sqrt(var2-moy2*moy2);

 //remplissage des masques par la methode z Map
 for (i=0;i<wdth;i++)
 {
  for (j=0;j<hght;j++)
  {
   for (k=0;k<dpth;k++)
   {

    val1=(double)im1->mri[i][j][k];
    val2=(double)im2->mri[i][j][k];
    if ((val1>0.0)&&(val2>0.0))
    {
     z1=fabs(val1-moy1)/sigma1;
     z2=fabs(val2-moy2)/sigma2;
     if ((z1<SEUIL_Z_MAP)&&(z2<SEUIL_Z_MAP)) mask1->mri[i][j][k]=mask2->mri[i][j][k]=1;
     else mask1->mri[i][j][k]=mask2->mri[i][j][k]=0;
    }
    else mask1->mri[i][j][k]=mask2->mri[i][j][k]=0;

   }
  }
 }

 return 0;
}

//-----------------------------------------------------------------------------//
//                             METHODE DE KOO                                  //
//-----------------------------------------------------------------------------//

/*******************************************************************************
**        soustraction_spect_Koo_3d()
*/
/*!       calcul d'une image soustraction d'un spect ictal et d'un spect interictal recale
**        methode de Koo (http://www.blackwell-synergy.com/links/doi/10.1046/j.1528-1157.2003.29402.x/full/)
**
*******************************************************************************/
void soustraction_spect_Koo_3d()
{
 int im_ictal, im_interictal, im_res ;

 /*choix des images*/
 im_ictal=GET_PLACE3D("choix de l'image ictale");
 im_interictal=GET_PLACE3D("choix de l'image inter ictale");
 im_res=GET_PLACE3D("ou placer l'image soustraction");

 imx_soustraction_spect_Koo_3d(im_ictal, im_interictal, im_res);

 show_picture_3d(im_res);

 return;
}

/*******************************************************************************
**        imx_soustraction_spect_Koo_3d(int im_ictal, int im_interictal, int im_res)
*/
/*!       calcul d'une image soustraction d'un spect ictal et d'un spect interictal recale
**        methode de Koo (http://www.blackwell-synergy.com/links/doi/10.1046/j.1528-1157.2003.29402.x/full/)
**
**        \param im_ictal, im_interictal, im_res: les images
**
*******************************************************************************/
void imx_soustraction_spect_Koo_3d(int im_ictal, int im_interictal, int im_res)
{
 grphic3d *imictal, *iminterictal, *imres;

 imictal=ptr_img_3d(im_ictal);
 iminterictal=ptr_img_3d(im_interictal);
 imres=ptr_img_3d(im_res);

 imx_soustraction_spect_Koo_3d_p(imictal, iminterictal, imres);

 return;
}

/*******************************************************************************
**        soustraction_spect_Koo_3d()
*/
/*!       calcul d'une image soustraction d'un spect ictal et d'un spect interictal recale
**        methode de Koo (http://www.blackwell-synergy.com/links/doi/10.1046/j.1528-1157.2003.29402.x/full/)
**
**        \param imictal, iminterictal, imres: les images
**
*******************************************************************************/
void imx_soustraction_spect_Koo_3d_p(grphic3d *imictal, grphic3d *iminterictal, grphic3d *imres)
{
 int i,j,k;
 int wdth, hght, dpth;
 int diff;
 double moyenne=0.0, variance=0.0, stddev=0.0, minval=0.0;
 double nboverlap=0.0;

 wdth=imictal->width; hght=imictal->height; dpth=imictal->depth;
 imx_copie_param_3d_p(imictal,imres);

 for (i=0;i<wdth;i++)
 {
  for (j=0;j<hght;j++)
  {
   for (k=0;k<dpth;k++)
   {
    if ((imictal->mri[i][j][k]>0)&&(iminterictal->mri[i][j][k]>0))
    {
     diff=imictal->mri[i][j][k]-iminterictal->mri[i][j][k];
     if (diff>0) imres->mri[i][j][k]=diff;
     else imres->mri[i][j][k]=0;
     nboverlap+=1.0;
     moyenne+=(double)diff;
     variance+=(double)(diff*diff);
    }
    else imres->mri[i][j][k]=0;
   }
  }
 }

 moyenne=moyenne/nboverlap;
 variance=variance/nboverlap-moyenne*moyenne;
 stddev=sqrt(variance);
 minval=moyenne+2.0*stddev;

 for (i=0;i<wdth;i++)
 {
  for (j=0;j<hght;j++)
  {
   for (k=0;k<dpth;k++)
   {
    if ((double)imres->mri[i][j][k]<=minval) imres->mri[i][j][k]=0;
   }
  }
 }

 return;
}

///*
//  Cette fonction normalise l'histogramme joint entre deux images (reprise de atlas_3d.c)
//
//  Le principe est le suivant :
//  1 on calcule l'histogramme joint des deux images
//  2 pour chaque valeur de l'image source, on cherche la valeur de l'image
//    destination correspondante : pour ce faire, on calcule la valeur mediane
//    dans la colonne correspondante dans l'histogramme joint
//  3 pour eviter d'avoir des valeurs statistiquement non significatives,
//    on compte sur combien de valeurs a ete calculee chaque colonne de
//    l'histogramme. Si ce n'est pas superieur a un certain seuil, on moyenne
//    la valeur mediane avec les valeurs trouvees par interpolation lineaire
//    des valeurs medianes des colonnes voisines, jusqu'a ce qu'on atteigne
//    le seuil voulu.
//  4 On interpole lineairement les valeurs trouvees pour l'histogramme, ce
//    qui nous donne la fonction de transfert a appliquer a l'image.
//*/
//
//void normalisation_histo_joint_3d(void)
//{
//  int im_ref, im_src, im_dest;
//
//  im_ref = GET_PLACE3D("image de reference");
//  im_src = GET_PLACE3D("image a normaliser");
//  im_dest = GET_PLACE3D("image resultat");
//
//  imx_normalisation_histo_joint_3d(im_ref, im_src, im_dest);
//}
//
//
//int imx_normalisation_histo_joint_3d(int im_ref, int im_src, int im_res)
//{
//  grphic3d *imref;
//  grphic3d *imsrc;
//  grphic3d *imres;
//  int err;
//
//  imref = ptr_img_3d(im_ref);
//  imsrc = ptr_img_3d(im_src);
//  imres = ptr_img_3d(im_res);
//
//  if ((! imref) || (! imsrc) || (! imres))
//    return 1;
//
//  err = imx_normalisation_histo_joint_3d_p(imref, imsrc, imres);
//
//  show_picture_3d(im_res);
//
//  return err;
//}
//
//int imx_normalisation_histo_joint_3d_p(grphic3d *imref, grphic3d *imsrc, grphic3d *imres)
//{
// int nbSrcBins=0, nbRefBins=0;
// long **hist_joint=NULL;
// int *median=NULL;
// int *confiance=NULL;
// double **T_base=NULL;
// int seuil_confiance = 0;
// unsigned int x, y, z;
// TDimension wdth, hght, dpth;
// int i, j;
// TYPEMRI3D T[MAXMRI3D+1]={0};	// tableau pour la table de traduction
// int nb_pts_base=0;
// double val_courante_src, val_courante_ref;
// int min_src=0, max_src=0, min_ref=0, max_ref=0;
// int nbVoxelOverlap=0;
// double coef_courant_src=0.0, coef_courant_ref=0.0;
// int borne_min_src, borne_max_src;
// double *T_base_tmp=NULL;
// //variables pour le filtrage passe-bas de la table de traduction
// int taille_filtre_moyenneur=2;
// double somme_voisins=0.0;
// int nb_somme_voisins=0;
//
// //calcul des bins d'histogramme a affecter au calcul de l'histogramme conjoint
// nbSrcBins=imx_calculer_nombre_de_classes_3d_p(imsrc, &max_src);
// nbRefBins=imx_calculer_nombre_de_classes_3d_p(imref, &max_ref);
//
// //allocations memoire
// hist_joint=CALLOC(nbSrcBins, long*);
// for (i=0;i<nbSrcBins;i++) hist_joint[i]=CALLOC(nbRefBins, long);
// median=CALLOC(nbSrcBins, int);
// confiance=CALLOC(nbSrcBins, int);
// T_base=CALLOC(nbSrcBins, double*);
// for (i=0;i<nbSrcBins;i++) T_base[i]=CALLOC(2, double);
//
// wdth = imref->width; hght = imref->height; dpth = imref->depth;
//
//  if ((imsrc->width != wdth) || (imsrc->height != hght) || (imsrc->depth != dpth))
//  {
//    fprintf (stderr, "tailles d'images incompatibles dans imx_normalisation_histo_joint_3d_p\n");
//    return 1;
//  }
//
//  //calcul des valeurs min positives non nulles des images
//  for (x=0; x<wdth; x++)
//   for (y=0; y<hght; y++)
//    for (z=0; z<dpth; z++)
//    {
//     val_courante_src=imsrc->mri[x][y][z];
//     if ((val_courante_src>0)&&(val_courante_src<min_src)) min_src=(int)val_courante_src;
//     val_courante_ref=imref->mri[x][y][z];
//     if ((val_courante_ref>0)&&(val_courante_ref<min_ref)) min_ref=(int)val_courante_ref;
//    }
//
//  //construction de l'histogramme conjoint
//  for (x=0; x<wdth; x++)
//   for (y=0; y<hght; y++)
//    for (z=0; z<dpth; z++)
//    {
//      if (imsrc->mri[x][y][z]<=0) i=0;
//      else i=((nbSrcBins-1)*imsrc->mri[x][y][z])/(max_src+1-min_src)+1;
//      if (imref->mri[x][y][z]<=0) j=0;
//      else j=((nbRefBins-1)*imref->mri[x][y][z])/(max_ref+1-min_ref)+1;
//      hist_joint[i][j]++;
//    }
//
// //calcul des valeurs medianes de l'histogramme correspondant a chaque valeur requantifiee de imsrc
// for (i=1;i<nbSrcBins;i++)
// {
//  long total=0, demi_total=0;
//
//  for (j=1; j<nbRefBins; j++) total += hist_joint[i][j];
//  confiance[i] = total;
//
//  if (total)
//  {
//   long somme=0;
//   demi_total=total/2;
//
//   for (j=1;j<nbRefBins;j++)
//   {
//    somme+=hist_joint[i][j];
//    if (somme>demi_total) { median[i]=j; break; }
//   }
//  }
//  nbVoxelOverlap += total;
// }
// seuil_confiance=(int)floor(0.01*nbVoxelOverlap);
//
// //calcul des points de bases pour le remplissage de la table de traduction
// coef_courant_src=(double)(max_src-min_src+1)/(double)(nbSrcBins-1);
// coef_courant_ref=(double)(max_ref-min_ref+1)/(double)(nbRefBins-1);
// for (i=1;i<nbSrcBins;i++)
// {
//  //on s'assure que le calcul de la mediane s'est effectuee sur un nombre de voxels statistiquement significatif
//  if (confiance[i]>seuil_confiance)
//  {
//   val_courante_src=(double)i*coef_courant_src+min_src;
//   val_courante_ref=(double)median[i]*coef_courant_ref+min_ref;
//   T_base[nb_pts_base][0]=val_courante_src;
//   T_base[nb_pts_base][1]=val_courante_ref;
//   nb_pts_base++;
//  }
// }
// //on regularise la courbe obtenue plus haut par un filtre moyenneur
// T_base_tmp=CALLOC(nb_pts_base, double);
// for (i=0;i<nb_pts_base;i++)
// {
//  somme_voisins=0.0;
//  nb_somme_voisins=0;
//  for (j=(i-taille_filtre_moyenneur);j<=(i+taille_filtre_moyenneur);j++)
//  {
//   if ((j>0)&&(j<(nb_pts_base-1))) { somme_voisins+=T_base[j][1]; nb_somme_voisins++; }
//  }
//  T_base_tmp[i]=somme_voisins/(double)nb_somme_voisins;
// }
// for (i=0;i<nb_pts_base;i++)
// {
//  T_base[i][1]=T_base_tmp[i];
// }
// FREE(T_base_tmp);
//
// //remplissage de la table de traduction par interpolation lineaire a partir des points de base
// if (nb_pts_base<2)
// {
//  fprintf (stderr, "pas assez de points pour calculer la table de traduction dans imx_normalisation_histo_joint_3d_p\n");
//  return 2;
// }
// coef_courant_src=(double)(T_base[1][1]-T_base[0][1])/(double)(T_base[1][0]-T_base[0][0]);
// borne_min_src=min_src; borne_max_src=(int)floor(T_base[1][0]+0.5);
// for (i=borne_min_src;i<borne_max_src;i++)
// {
//  T[i]=(int)floor(coef_courant_src*(double)(i-T_base[0][0])+T_base[0][1]+0.5);
// }
// for (j=1;j<(nb_pts_base-2);j++)
// {
//  coef_courant_src=(double)(T_base[j+1][1]-T_base[j][1])/(double)(T_base[j+1][0]-T_base[j][0]);
//  borne_min_src=(int)floor(T_base[j][0]+0.5); borne_max_src=(int)floor(T_base[j+1][0]+0.5);
//  for (i=borne_min_src;i<borne_max_src;i++)
//  {
//   T[i]=(int)floor(coef_courant_src*(double)(i-T_base[j][0])+T_base[j][1]+0.5);
//  }
// }
// coef_courant_src=(double)(T_base[nb_pts_base-1][1]-T_base[nb_pts_base-2][1])/(double)(T_base[nb_pts_base-1][0]-T_base[nb_pts_base-2][0]);
// borne_min_src=(int)floor(T_base[nb_pts_base-2][0]+0.5); borne_max_src=max_src;
// for (i=borne_min_src;i<=borne_max_src;i++)
// {
//  T[i]=(int)floor(coef_courant_src*(double)(i-T_base[nb_pts_base-2][0])+T_base[nb_pts_base-2][1]+0.5);
// }
///*
// for (i=0;i<nb_pts_base;i++)
// {
//  printf ("%d,%d\n",(int)floor(T_base[i][0]+0.5),(int)floor(T_base[i][1]+0.5));
// }*/
///*
// for (i=min_src;i<=max_src;i++)
// {
//  printf ("%d,%d\n",i,T[i]);
// }
//*/
//  // et on traduit
//  imx_copie_param_3d_p(imsrc, imres);
//
//  for (x=0; x<wdth; x++)
//   for (y=0; y<hght; y++)
//    for (z=0; z<dpth; z++)
//    {
//      if (imsrc->mri[x][y][z] <= 0) imres->mri[x][y][z] = 0;
//      else imres->mri[x][y][z] = MAXI(T[imsrc->mri[x][y][z]],0);
//    }
//
//  imx_inimaxminpixel_3d_p(imres);
//
// //liberations memoire
// for (i=0;i<nbSrcBins;i++) FREE(T_base[i]);
// FREE(T_base);
// FREE(confiance);
// FREE(median);
// for (i=0;i<nbSrcBins;i++) FREE(hist_joint[i]);
// FREE(hist_joint);
//
//
//  return 0;
//}
//
//void normalisation_lineaire_masque_3d()
//{
// int im_ref, im_src, im_dest;
// int norm_type;
//
// im_ref = GET_PLACE3D("image de reference");
// im_src = GET_PLACE3D("image a normaliser");
// im_dest = GET_PLACE3D("image resultat");
// norm_type = imx_query_normalisation_lineaire_type_3d();
//
// imx_normalisation_lineaire_masque_3d(im_ref, im_src, im_dest, norm_type);
//}
//
//int imx_normalisation_lineaire_masque_3d(int im_ref, int im_src, int im_dest, int norm_type)
//{
// grphic3d *imref, *imsrc, *imres;
// enum NORMLINTYPE normLinType;
// int err;
//
// imref = ptr_img_3d(im_ref);
// imsrc = ptr_img_3d(im_src);
// imres = ptr_img_3d(im_dest);
//
// if ((! imref) || (! imsrc) || (! imres)) return 1;
// normLinType=imx_choose_normalisation_lineaire_type_3d(norm_type);
//
// err = imx_normalisation_lineaire_masque_3d_p(imref, imsrc, imres, normLinType);
//
// show_picture_3d(im_dest);
//
// return err;
//}
//
//int imx_normalisation_lineaire_masque_3d_p(grphic3d *imref, grphic3d *imsrc, grphic3d *imres, enum NORMLINTYPE normLinType)
//{
// double alpha=0.0, beta=0.0;
// TDimension wdth, hght, dpth;
// unsigned int i,j,k;
//
// wdth=imref->width; hght=imref->height; dpth=imref->depth;
//
// if ((wdth!=imsrc->width)||(hght!=imsrc->height)||(dpth!=imsrc->depth))
// {
//  fprintf (stderr,"tailles d'images incompatibles dans imx_normalisation_lineaire_3d_p\n");
//  return 1;
// }
//
// if (!(imref->mask))
// {
//  fprintf (stderr,"le masque de imref manque dans imx_normalisation_lineaire_3d_p\n");
//  return 2;
// }
//
// switch (normLinType)
// {
//  case REGLIN1: imx_regression_lineaire_masque_1_param_3d_p(imref, imsrc, &alpha, &beta); break;
//  case REGLIN2: imx_regression_lineaire_masque_2_param_3d_p(imref, imsrc, &alpha, &beta); break;
//  case REGTLS: imx_regression_TLS_masque_3d_p(imref, imsrc, &alpha, &beta); break;
//  default: imx_regression_lineaire_masque_1_param_3d_p(imref, imsrc, &alpha, &beta); break;
// }
//
// imx_copie_param_3d_p(imsrc, imres);
//
// for (i=0;i<wdth;i++)
// {
//  for (j=0;j<hght;j++)
//  {
//   for (k=0;k<dpth;k++)
//   {
//    if (imsrc->mri[i][j][k]>0)
//    {
//     imres->mri[i][j][k]=(TYPEMRI3D)floor(alpha*(double)imsrc->mri[i][j][k]+beta+0.5);
//    }
//    else imres->mri[i][j][k]=0;
//   }
//  }
// }
//
// imx_inimaxminpixel_3d_p(imres);
//
// return 0;
//}
//
////---------TOOLBOX POUR LA NORMALISATION LINEAIRE-----//
//
///*******************************************************************************
//**    imx_regression_lineaire_masque_1_param_3d_p(grphic3d *imref, grphic3d *imsrc, double *alpha, double *beta)
//*/
///******************************************************************************/
//int imx_regression_lineaire_masque_1_param_3d_p(grphic3d *imref, grphic3d *imsrc, double *alpha, double *beta)
//{
// TDimension wdth, hght, dpth;
// unsigned int i,j,k;
// double covariance=0.0, variance=0.0;
// double valsrc=0.0, valref=0.0, surfMasque=0.0;
//
// wdth=imref->width; hght=imref->height; dpth=imref->depth;
//
// if ((wdth!=imsrc->width)||(hght!=imsrc->height)||(dpth!=imsrc->depth))
// {
//  fprintf (stderr,"tailles d'images incompatibles dans imx_regression_lineaire_masque_1_param_3d_p\n");
//  return 1;
// }
//
// if (!(imref->mask))
// {
//  fprintf (stderr,"le masque de imref manque dans imx_regression_lineaire_masque_1_param_3d_p\n");
//  return 2;
// }
//
// //calcul de la surface du masque
// for (i=0;i<wdth;i++)
// {
//  for (j=0;j<hght;j++)
//  {
//   for (k=0;k<dpth;k++)
//   {
//     if (imref->mask->mri[i][j][k]) surfMasque+=1.0;
//   }
//  }
// }
//
// //calcul des parametres de correlation des images
// for (i=0;i<wdth;i++)
// {
//  for (j=0;j<hght;j++)
//  {
//   for (k=0;k<dpth;k++)
//   {
//    if (imref->mask->mri[i][j][k])
//    {
//     valref=(double)imref->mri[i][j][k];
//     valsrc=(double)imsrc->mri[i][j][k];
//     covariance+=valref*valsrc/surfMasque;
//     variance+=valsrc*valsrc/surfMasque;
//    }
//   }
//  }
// }
//
// //on met les valeurs a jour
// *alpha=covariance/variance;
// *beta=0.0;
//
// printf ("alpha: %f, beta: %f\n", *alpha, *beta);
//
// return 0;
//}
//
//void regression_lineaire_masque_1_param_3d()
//{
//  int im_ref, im_src;
//
//  im_ref = GET_PLACE3D("image de reference");
//  im_src = GET_PLACE3D("image a normaliser");
//
//  imx_regression_lineaire_masque_1_param_3d(im_ref, im_src);
//}
//
//int imx_regression_lineaire_masque_1_param_3d(int im_ref, int im_src)
//{
// double alpha=0.0, beta=0.0;
// grphic3d *imref;
// grphic3d *imsrc;
//
// imref = ptr_img_3d(im_ref);
// imsrc = ptr_img_3d(im_src);
//
// if ((! imref) || (! imsrc)) return 1;
//
// imx_regression_lineaire_masque_1_param_3d_p(imref, imsrc, &alpha, &beta);
//
// return 0;
//}
//
//int imx_regression_lineaire_masque_2_param_3d_p(grphic3d *imref, grphic3d *imsrc, double *alpha, double *beta)
//{
// TDimension wdth, hght, dpth;
// unsigned int i,j,k;
// double covariance=0.0, variance=0.0, moyRef=0.0, moySrc=0.0;
// double valsrc=0.0, valref=0.0, surfMasque=0.0;
//
// wdth=imref->width; hght=imref->height; dpth=imref->depth;
//
// if ((wdth!=imsrc->width)||(hght!=imsrc->height)||(dpth!=imsrc->depth))
// {
//  fprintf (stderr,"tailles d'images incompatibles dans imx_regression_lineaire_masque_2_param_3d_p\n");
//  return 1;
// }
//
// if (!(imref->mask))
// {
//  fprintf (stderr,"le masque de imref manque dans imx_regression_lineaire_masque_2_param_3d_p\n");
//  return 2;
// }
//
// //calcul de la surface du masque
// for (i=0;i<wdth;i++)
// {
//  for (j=0;j<hght;j++)
//  {
//   for (k=0;k<dpth;k++)
//   {
//     if (imref->mask->mri[i][j][k])
//     {
//      surfMasque+=1.0;
//      moyRef+=imref->mri[i][j][k];
//      moySrc+=imsrc->mri[i][j][k];
//     }
//   }
//  }
// }
//
// moyRef/=surfMasque;
// moySrc/=surfMasque;
//
// //calcul des parametres de correlation des images
// for (i=0;i<wdth;i++)
// {
//  for (j=0;j<hght;j++)
//  {
//   for (k=0;k<dpth;k++)
//   {
//    if (imref->mask->mri[i][j][k])
//    {
//     valref=(double)imref->mri[i][j][k];
//     valsrc=(double)imsrc->mri[i][j][k];
//     covariance+=(valref-moyRef)*(valsrc-moySrc)/surfMasque;
//     variance+=(valsrc-moySrc)*(valsrc-moySrc)/surfMasque;
//    }
//   }
//  }
// }
//
// //on met les valeurs a jour
// *alpha=covariance/variance;
// *beta=moyRef-(*alpha)*moySrc;
//
// printf ("alpha: %f, beta: %f\n", *alpha, *beta);
//
// return 0;
//}
//
//void regression_lineaire_masque_2_param_3d()
//{
//  int im_ref, im_src;
//
//  im_ref = GET_PLACE3D("image de reference");
//  im_src = GET_PLACE3D("image a normaliser");
//
//  imx_regression_lineaire_masque_2_param_3d(im_ref, im_src);
//}
//
//int imx_regression_lineaire_masque_2_param_3d(int im_ref, int im_src)
//{
// double alpha=0.0, beta=0.0;
// grphic3d *imref;
// grphic3d *imsrc;
//
// imref = ptr_img_3d(im_ref);
// imsrc = ptr_img_3d(im_src);
//
// if ((! imref) || (! imsrc)) return 1;
//
// imx_regression_lineaire_masque_2_param_3d_p(imref, imsrc, &alpha, &beta);
//
// return 0;
//}
//
//
//void regression_TLS_masque_3d()
//{
//  int im_ref, im_src;
//
//  im_ref = GET_PLACE3D("image de reference");
//  im_src = GET_PLACE3D("image a normaliser");
//
//  imx_regression_TLS_masque_3d(im_ref, im_src);
//}
//
//int imx_regression_TLS_masque_3d(int im_ref, int im_src)
//{
// double alpha=0.0, beta=0.0;
// grphic3d *imref;
// grphic3d *imsrc;
//
// imref = ptr_img_3d(im_ref);
// imsrc = ptr_img_3d(im_src);
//
// if ((! imref) || (! imsrc)) return 1;
//
// imx_regression_TLS_masque_3d_p(imref, imsrc, &alpha, &beta);
//
// return 0;
//}
//
//int imx_regression_TLS_masque_3d_p(grphic3d *imref, grphic3d *imsrc, double *alpha, double *beta)
//{
// TDimension wdth, hght, dpth;
// double covariance=0.0, varianceRef=0.0, varianceSrc=0.0, moyRef=0.0, moySrc=0.0, sigmaRef=0.0, sigmaSrc=0.0;
// int surfMasque=0, taille_covar=0;
// double /*x, y, */r1=0.0, r2=0.0;
// param_min_TLS param_TLS_calc;
//
// wdth=imref->width; hght=imref->height; dpth=imref->depth;
//
// if ((wdth!=imsrc->width)||(hght!=imsrc->height)||(dpth!=imsrc->depth))
// {
//  fprintf (stderr,"tailles d'images incompatibles dans imx_regression_TLS_masque_3d_p\n");
//  return 1;
// }
//
// if (!(imref->mask))
// {
//  fprintf (stderr,"le masque de imref manque dans imx_regression_TLS_masque_3d_p\n");
//  return 2;
// }
//
// //calcul des tailles des masques et des moyennes des images sur les masques
// imx_inimaxminpixel_3d_p(imref);
// imx_inimaxminpixel_3d_p(imsrc);
// surfMasque=imx_calc_moy_im_mask_3d_p(imref,&moyRef);
// surfMasque=imx_calc_moy_im_mask_3d_p(imsrc,&moySrc);
//
// if (surfMasque<=0)
// { fprintf (stderr, "impossible de normaliser pour cause de taille de masque nulle dans imx_regression_TLS_masque_3d_p\n"); return 3; }
//
// //calcul des ecarts types et covariance des images
// imx_calc_sigma_im_mask_3d_p(imref, surfMasque, moyRef, &sigmaRef);
// imx_calc_sigma_im_mask_3d_p(imsrc, surfMasque, moySrc, &sigmaSrc);
// taille_covar=imx_calc_covar_im_mask_3d_p(imref, imsrc, moyRef, moySrc, &covariance);
//
// if (taille_covar<=0)
// { fprintf (stderr, "impossible de normaliser pour cause de taille d'overlap des masques nulle dans imx_regression_TLS_masque_3d_p\n"); return 3; }
//
// varianceSrc=sigmaSrc*sigmaSrc; varianceRef=sigmaRef*sigmaRef;
// param_TLS_calc.covariance=covariance;  param_TLS_calc.varianceSrc=varianceSrc;  param_TLS_calc.varianceRef=varianceRef;
//
///* for (x=-1.0;x<=1.0;x+=0.001)
// {
//  y=fonction_energie_TLS(x, &param_TLS_calc);
//  printf ("%f, %f\n",x, y);
// }*/
//
// //on effecture une minimisation lineaire pour retrouver les parametres de la droite
// linemin_general_3d((lin_energie_func_t)fonction_energie_TLS, &param_TLS_calc, -1.0, 1.0, &r1);
//
// r2=sqrt(1.0-r1*r1);
//
// printf ("XMIN: %f\n", r1);
// //on met les valeurs a jour
// *alpha=-r1/r2;
// *beta=r1/r2*moySrc+moyRef;
//
// printf ("alpha: %f, beta: %f\n", *alpha, *beta);
//
// return 0;
//}

/***    imx_regression_TLS_roi_3d_p()   ------------------------------------------------------
***	NORMALISATION D'UNE IMAGE 3D PAR RAPPORT A UNE AUTRE A L'AIDE DES MASQUES
***
***    normalisation en utilisant la regression lineaire avec la methode total least squares
***
****************************************************************/
void imx_regression_TLS_roi_3d_p(grphic3d *imsrc, grphic3d *imref, grphic3d *imres)
{
 TDimension wdth, hght, dpth;
 double covariance=0.0, varianceRef=0.0, varianceSrc=0.0, moyRef=0.0, moySrc=0.0, sigmaRef=0.0, sigmaSrc=0.0;
 int surfMasque=0, taille_covar=0;
 double r1=0.0, r2=0.0;
 param_min_TLS param_TLS_calc;
 double alpha=0.0, beta=0.0;

 wdth=imref->width; hght=imref->height; dpth=imref->depth;

 if ((wdth!=imsrc->width)||(hght!=imsrc->height)||(dpth!=imsrc->depth))
 {
  fprintf (stderr,"tailles d'images incompatibles dans imx_regression_TLS_masque_3d_p\n");
  return ;
 }

 if (!(imref->mask))
 {
  fprintf (stderr,"le masque de imref manque dans imx_regression_TLS_masque_3d_p\n");
  return ;
 }

 //calcul des tailles des masques et des moyennes des images sur les masques
 imx_inimaxminpixel_3d_p(imref);
 imx_inimaxminpixel_3d_p(imsrc);
 surfMasque=imx_calc_moy_im_mask_3d_p(imref,&moyRef);
 surfMasque=imx_calc_moy_im_mask_3d_p(imsrc,&moySrc);

 if (surfMasque<=0)
 { fprintf (stderr, "impossible de normaliser pour cause de taille de masque nulle dans imx_regression_TLS_masque_3d_p\n"); return ; }

 //calcul des ecarts types et covariance des images
 imx_calc_sigma_im_mask_3d_p(imref, surfMasque, moyRef, &sigmaRef);
 imx_calc_sigma_im_mask_3d_p(imsrc, surfMasque, moySrc, &sigmaSrc);
 taille_covar=imx_calc_covar_im_mask_3d_p(imref, imsrc, moyRef, moySrc, &covariance);

 if (taille_covar<=0)
 { fprintf (stderr, "impossible de normaliser pour cause de taille d'overlap des masques nulle dans imx_regression_TLS_masque_3d_p\n"); return ; }

 varianceSrc=sigmaSrc*sigmaSrc; varianceRef=sigmaRef*sigmaRef;

 param_TLS_calc.covariance=covariance;  param_TLS_calc.varianceSrc=varianceSrc;  param_TLS_calc.varianceRef=varianceRef;

 //on effecture une minimisation lineaire pour retrouver les parametres de la droite
 linemin_general_3d((lin_energie_func_t)fonction_energie_TLS, &param_TLS_calc, -1.0, 1.0, &r1);

 r2=sqrt(1.0-r1*r1);

 printf ("XMIN: %f\n", r1);
 //on met les valeurs a jour
 alpha=-r1/r2;
 beta=r1/r2*moySrc+moyRef;
 printf ("alpha: %f, beta: %f\n", alpha, beta);

 //normalisation lineaire de l'image1
 imx_norm_lin_im_3d_p(imsrc, imres, alpha, beta);

 // Mettre le meme rcoeff pour les deux images
 imx_norm_rcoeff_3d_p(imref,imres);

 return ;
}

///*******************************************************************************
//**    imx_query_normalisation_lineaire_type_3d()
//*/
///*!
//**    interface permettant de choisir une methode pour la normalisation lineaire
//**    retval : l'int correspondant a la methode choisie
//**
//*******************************************************************************/
//int imx_query_normalisation_lineaire_type_3d()
//{
// char *query[100];
// query[0]="regression lineaire a 1 param";
// query[1]="regression lineaire a 2 param";
// query[2]="regression total least squares";
// query[3]=NULL;
// return GETV_QCM("type de normalisation lineaire",(char **)query);
//}
//
///*******************************************************************************
//**    imx_choose_normalisation_lineaire_type_3d(norm_type)
//*/
///*!
//**    renvoie le NORMLINTYPE correspondant a la methode de normalisation lineaire choisie
//**    retval : le bon NORMLINTYPE
//**
//*******************************************************************************/
//enum NORMLINTYPE imx_choose_normalisation_lineaire_type_3d(int norm_type)
//{
// enum NORMLINTYPE normLinType;
// switch (norm_type)
// {
//  case 0: normLinType=REGLIN1;aff_log("REG LIN 1 PARAM ");break;
//  case 1: normLinType=REGLIN2;aff_log("REG LIN 2 PARAM ");break;
//  case 2: normLinType=REGLIN2;aff_log("REG TLS ");break;
//  default:normLinType=REGLIN1;aff_log("REG LIN 1 PARAM ");break;
// }
//
// return normLinType;
//}
//
/*******************************************************************************
**    fonction_energie_TLS(x, params_image)
*/
/*!
**    fonction d'enrgie a minimiser dans la cadre de la normalisation total least squares
**    param x : parametre de minimisation
**    param params_image : les parametres additionnels calcules a partir des images
**    retval : l'energie calculee
**
*******************************************************************************/
double fonction_energie_TLS(double x, param_min_TLS *params_image)
{
 double covariance=0.0, varianceRef=0.0, varianceSrc=0.0, val_cal=0.0;

 covariance=params_image->covariance;
 varianceRef=params_image->varianceRef;
 varianceSrc=params_image->varianceSrc;

 val_cal=x*x*(varianceSrc-varianceRef)+varianceRef+2.0*x*sqrt(1-x*x)*covariance;

 return val_cal;
}

/*******************************************************************************
**    linemin_general_3d(fonction_energie, params_additionnels, borne_min, borne_max, abscisse_trouve)
*/
/*!
**    minimise une fonction d'energie selon un parametre
**    param fonction_energie : la fonction d'energie a minimiser
**    param params_additionnels : les parametres additionneles de la fonction d'energie
**    param borne_min, borne_max, abscisse_trouve:  les bornes min, max et calculee pour la minimisation
**    retval : l'energie min
**
*******************************************************************************/
double linemin_general_3d(lin_energie_func_t fonction_energie, void *params_additionnels, double borne_min, double borne_max, double *abscisse_trouve)
{
 //les cf representent les coefs utilise pour se deplacer sur la droite passant par 'param' et ayant pour direction 'direc'
 //param2[i] = param[i] + coef*direc[i];
 double cfmin=0.0,cfmax=0.0, cfx, cfw, cfv, cfu, cfp, cfq, cfr, cfprec, cfm, cfetemp, cfd=0.0, cfopt, cfe=0.0;
 //parametres de tolerance tenant compte de la precision qu'on desire avoir sur le calcul des parametres optimaux 'prec_param'
 //et des erreurs de calculs dus a la precision machine
 double tol1, tol2;
 //energies calculees pour les parametres correspondant aux coefs: Ew -> cfw, EBorneMin -> cfmin ...
 double Ew, Ev, Ex, Eu, EParam, Eopt;
 //parametres correspondant aux coefs: paramx -> cfx, param_min -> cfmin ...
 int iter;
 int nb_minimisations=0;

 //initisalition de cfmin et cfmax
 cfmin = borne_min; cfmax=borne_max; cfprec=1e-4;

 //calcul de l'energie correspondant aux parametres de depart
 EParam=fonction_energie(*abscisse_trouve, params_additionnels);

 ///////initialisation///////

 cfx=cfw=cfv=0.0;
 Ew=Ev=Ex=EParam;

 ///////minimisation///////

 for (iter=0;iter<NB_ITER_MAX_LINEMIN_GENERAL_3D;iter++)
 {
  cfm=0.5*(cfmin+cfmax);
  tol1=cfprec*fabs(cfx)+ZEPS;
  tol2=2.0*tol1;
  //si on a atteint la precision voulue, on s'arrete
  if (fabs(cfx-cfm)<=(tol2-0.5*(cfmax-cfmin)))
  {
   break;
  }
  //l'intervalle dans lequel on va rechercher le nouveau point et plus grand que tol1:
  //on va tenter une interpolation parabolique avec les points correspondant a cfx, cfv et cfw
  if (fabs(cfe)>tol1)
  {
   cfr=(cfx-cfw)*(Ex-Ev);
   cfq=(cfx-cfv)*(Ex-Ew);
   cfp=(cfx-cfv)*cfq-(cfx-cfw)*cfr;
   cfq=2.0*(cfq-cfr);
   if (cfq>0.0) cfp=-cfp;
   cfq=fabs(cfq);
   cfetemp=cfe;
   cfe=cfd;
   //on teste si l'interpolation parabolique est pertinente
   //ie les 3 points ne doivent pas etre alignes
   if ((fabs(cfp)>=fabs(0.5*cfq*cfetemp))||(cfp<=cfq*(cfmin-cfx))||(cfp>=cfq*(cfmax-cfx)))
    //ca n'a pas marche: on prend la section du nombre d'or dans le plus grand segment
    cfd=CGOLD*(cfe=(cfx >= cfm ? cfmin-cfx : cfmax-cfx));
   else
   {
    //ca a marche: on prend le minimum de la parabole
    cfd=cfp/cfq;
    cfu=cfx+cfd;
    if (((cfu-cfmin)<tol2)||((cfmax-cfu)<tol2))
     cfd=SGN(tol1,cfm-cfx);
   }
  }
  //l'intervalle dans lequel on va rechercher le nouveau point est trop petit:
  //on choisit un autre intervalle
  else
  {
   cfe=(cfx>=cfm ? cfmin-cfx : cfmax-cfx);
   cfd=CGOLD*cfe;
  }

  //calcul du coef correspondant a l'optimum determine, des parametres et de l'energie associes
  cfu=(fabs(cfd) >= tol1 ? cfx+cfd : cfx + SGN(tol1,cfd));
  Eu=fonction_energie(cfu, params_additionnels);

  //si le point determine est meilleur que l'optimum precedent:
  if (Eu<=Ex)
  {
   //on remet les donnees a jour avec le nouvel optimum calcule
   if (cfu>=cfx) cfmin=cfx; else cfmax=cfx;
   cfv=cfw; cfw=cfx; cfx=cfu;
   Ev=Ew; Ew=Ex; Ex=Eu;
   nb_minimisations++;
   if (nb_minimisations>=NB_ITER_MAX_LINEMIN_GENERAL_3D) break;
  }
  //s'il n'est pas meilleur:
  else
  {
   //on se readapte en consequence
   if (cfu<cfx) cfmin=cfu; else cfmax=cfu;
   if ((Eu<=Ew)||(cfw==cfx))
   {
    cfv=cfw;
    cfw=cfu;
    Ev=Ew;
    Ew=Eu;
   }
   else if ((Eu<=Ev)||(cfv=cfx)||(cfv==cfw))
   {
    cfv=cfu;
    Ev=Eu;
   }
  }
 }

 //calcul et renvoie du point optimal
 cfopt=cfx;
 Eopt=fonction_energie(cfopt, params_additionnels);
 *abscisse_trouve=cfopt;

 return Eopt;
}


