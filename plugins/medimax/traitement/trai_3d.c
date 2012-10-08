/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*!**********************************************************************
***
***	\file		trai_3d.c
***
***	project:	Imagix 2.01
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
*************************************************************************/
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_3d.h"
#include "noyau/imx_lang.h"
#include "traitement/trai_3d.h"
#include "math/oper_3d.h"
#include "recalage/mtch_3d.h"
#include "noyau/io/imx_log.h"
#include "outils/function_timer.h"
#include "math/imx_matrix.h"
#include "outils/imx_sort.h"
#include "noyau/mani_3d.h"
#include "math/oper_3d.h"
#include "noyau/io/imx_export_file.h"
//#include "noyau/surf/diag_3d.h"
#include "segmentation/otsu_3d.h"
#include "morpho/morpho_3d.h"
#include "noyau/gui/imx_picture3d.h"
#include "math/ana_3d.h"
#include "segmentation/segm_3d.h"
//#include "serie/seri_3d.h"


/*! \defgroup ToolBoxImage3d
  \brief Boite a outils
*/

#ifndef __DOXYGEN_SHOULD_SKIP_THIS__
extern void cdft (int *which, double *p, double *q, double *t, double *df, int *status, double *bound);
extern void cdfnor (int *which, double *p, double *q, double *x, double *mean, double *sd, int *status, double *bound);
extern void cdff (int *which, double *p, double *q, double *f, double *dfn, double *dfd, int *status, double *bound);
#endif /*__DOXYGEN_SHOULD_SKIP_THIS__*/

#define NORME_EUCLIDIENNE(x,y,z) sqrt((x)*(x)+(y)*(y)+(z)*(z))


/**************************************************
 ** -- d_euclidean_mapping_3d() ------------------------
 **
 **  Transformee de distance d'une image de contour
 **************************************************/


void d_euclidean_mapping_3d(void)
{

  int im_deb,im_res_x,im_res_y,im_res_z,im_res_m;

  im_deb=GET_PLACE3D("Original Image ?");
  im_res_x=GET_PLACE3D("Distance X ?");
  im_res_y=GET_PLACE3D("Distance Y ?");
  im_res_z=GET_PLACE3D("Distance Z ?");
  im_res_m=GET_PLACE3D("Module ?");

  imx_d_euclidean_mapping_3d(im_deb,im_res_x,im_res_y,im_res_z,im_res_m);


  imx_copie_param_3d(im_deb,im_res_x);
  imx_iniparaimg_3d(im_res_x);
  show_picture_3d(im_res_x);

  imx_copie_param_3d(im_deb,im_res_y);
  imx_iniparaimg_3d(im_res_y);
  show_picture_3d(im_res_y);

  imx_copie_param_3d(im_deb,im_res_z);
  imx_iniparaimg_3d(im_res_z);
  show_picture_3d(im_res_z);

  imx_copie_param_3d(im_deb,im_res_m);
  imx_iniparaimg_3d(im_res_m);
  show_picture_3d(im_res_m);

  return;

}




/********************************************
 **    labelcroix_3d()
 */
/*!    Fonction de labelisation 3D
**    A partir de la croissance de region
********************************************/
void	labelcroix_3d(void)
{

  int im_1,im_res;


  /*     Image   				*/
  im_1=GET_PLACE3D(TEXT0018);
  im_res=GET_PLACE3D(TEXT0006);


  imx_labelcroix_3d(im_1,im_res);

  show_picture_3d(im_res);
  return;
}


/*****************************************
 ** -- ordonner_labels_3d()------------------------------
 **    Fonction de rangement (tri) des labels 3D
 **    Si par exemple les labels ont des valeurs non
 **    continues dans N (i.e. 1 2 3 4 ... mais 1 3 8 12..)
 **    cette fonction range tout à partir de 1 et avec un
 **    incrément de 1
 **
 **    Rq : voir imx_ordonner_labels_3d et / ou imx_ordonner_labels_3d_p
 ******************************************/
extern void ordonner_labels_3d(void)
{
  int im_1,im_res;


  /*     Image   				*/
  im_1=GET_PLACE3D(TEXT0018);
  im_res=GET_PLACE3D(TEXT0006);


  printf("Lancement de l'ordonnancement des labels\n");
  imx_ordonner_labels_3d(im_1,im_res);
  printf("Fin de l'ordonnancement des labels\n");


  #ifdef __GTK_GUI
  imx_SetLabelPalette(im_res);
  #endif /*__GTK_GUI*/

  show_picture_3d(im_res);

  return;
}

/*****************************************
 ** -- labelcube_3d()------------------------------
 **    Fonction de labelisation 3D
 **    A partir de la croissance de region
 ******************************************/

int 	labelcube_3d(void)
{
  int im_1,im_res;


  /*     Image   				*/
  im_1=GET_PLACE3D(TEXT0018);
  im_res=GET_PLACE3D(TEXT0006);


  imx_labelcube_3d(im_1,im_res);

#ifdef __GTK_GUI
  imx_SetLabelPalette(im_res);
#endif /*__GTK_GUI*/

  show_picture_3d(im_res);

  return(1);
}


/*****************************************
 ** -- ctr_binaire_3d()------------------------------
 **    Fonction d'extraction des contours en 3D
 **    Utilise un masque en croix
 ******************************************/

void 	ctr_binaire_3d(void)
{
  int im_1,im_res;


  /*     Image   				*/
  im_1=GET_PLACE3D(TEXT0164);
  im_res=GET_PLACE3D(TEXT0006);


  imx_ctr_binaire_3d(im_1,im_res);

  show_picture_3d(im_res);

  return;
}


/*****************************************
 ** -- ctr_labele_3d()------------------------------
 **    Fonction d'extraction des contours en 3D
 **    Utilise voisinage en croix
 ******************************************/

void 	ctr_labele_3d(void)
{
  int im_1,im_res;


  /*     Image   				*/
  im_1=GET_PLACE3D(TEXT0164);
  im_res=GET_PLACE3D(TEXT0006);


  imx_ctr_labele_3d(im_1,im_res);

  show_picture_3d(im_res);

  return;
}


/*****************************************
 ** -- ctr_labele_rempl_3d()------------------------------
 **    Fonction d'extraction des contours en 3D
 **    Utilise voisinage en croix
 ******************************************/

void 	ctr_labele_rempl_3d(void)
{
  int im_1,im_res;


  /*     Image   				*/
  im_1=GET_PLACE3D(TEXT0164);
  im_res=GET_PLACE3D(TEXT0006);


  imx_ctr_labele_rempl_3d(im_1,im_res);

  show_picture_3d(im_res);

  return;
}

/*****************************************
 ** -- sync_labels_3d()------------------------------
 **    Fonction de synchronisation entre deux
 **    images de labels
 ******************************************/

void 	sync_labels_3d(void)
{
  int im_src,im_ref,im_res;


  /*     Image   				*/
  im_src=GET_PLACE3D(TEXT0164);
  im_ref=GET_PLACE3D(TEXT0185);
  im_res=GET_PLACE3D(TEXT0006);


  imx_sync_labels_3d(im_src,im_ref,im_res);

  show_picture_3d(im_res);

  return;
}

/*****************************************
 ** -- sync_labels_sep_3d()------------------------------
 **    Fonction de synchronisation entre deux
 **    images de labels pour la sclérose en plaques
 **    Auteur: F. Rousseau, 2007
 ******************************************/

void 	sync_labels_sep_3d(void)
{
  int im_src,im_ref,im_res;


  /*     Image   				*/
  //im_src=GET_PLACE3D(TEXT0164);
  //im_ref=GET_PLACE3D(TEXT0185);

  im_ref = GET_PLACE3D("Examen 1 (reference)");
  im_src = GET_PLACE3D("Examen 2 (source)");

  imx_sync_labels_sep_3d(im_src,im_ref);

  show_picture_3d(im_res);

  return;
}

/****************************************************
 *** canny_deriche_3d()
 ****************************************************/

void  canny_deriche_3d(void)
{

  int im_deb,im_res,e=0;
  double alpha;
  int window=0;
  int cc_length;

  im_deb=GET_PLACE3D("0riginal Image ?");
  im_res=GET_PLACE3D("Gradient Module Image ?");
  alpha=(double) GET_FLOAT("Alpha parameter (1) ?", 1.0, &e);
  cc_length=(int) GET_INT("Connected component length ?", 5, &e);
  if(e) return;

  window=imx_query_filter_dimension();

  imx_canny_deriche_3d(im_deb,im_res,alpha,window,cc_length);

  imx_copie_param_3d(im_deb,im_res);
  imx_iniparaimg_3d(im_res);
  show_picture_3d(im_res);

}


/*****************************************
 ** -- growing_cerv_3d()------------------------------
 **    Fonction d'extraction du cerveau en 3d
 **    A partir de la croissance de region
 ******************************************/

void 	growing_cerv_3d(void)
{
  int im_1,im_res;
  int *mousepos3d,err=1;
  int x0,y0,z0;        /*Coordonnes du point de depart de la croissance*/

  /*     Image   				*/
  im_res=GET_PLACE3D(TEXT0006);

  PUT_WNDMSG(TEXT0162);
  mousepos3d=GET_MOUSEPOS_3D(&err);
  CLEAR_WNDMSG();
  if (err)
    {
      /*ERREUR dans le GET_MOUSEPOS_3D*/
      PUT_WARN(TEXT0161);
    }
  else
    {
      x0=mousepos3d[3];y0=mousepos3d[4];z0=mousepos3d[5];
      im_1=mousepos3d[2];
      imx_growing_cerv_3d(im_1,im_res,x0,y0,z0);

      show_picture_3d(im_res);
    }

  free(mousepos3d);
  return;
}



/*****************************************
 ** -- extract_cerv_3d()------------------------------
 **    Fonction d'extraction du cerveau en 3d
 **    A partir de ostu avec la labelisation
 ******************************************/

void 	extract_cerv_3d(void)
{
  int im_1,im_res;

  /*     Image   				*/
  im_1=GET_PLACE3D(TEXT0164);
  im_res=GET_PLACE3D(TEXT0006);


  imx_extract_cerv_3d(im_1,im_res);

  show_picture_3d(im_res);

  return;
}


/*******************************************
 ** --  hole_fill_3d() ----------------------------
 **
 **    Fonction remplissage de trous en 3d (d'une image binaire)
 **    A partir de la croissance de region
 ********************************************/
void    hole_fill_3d(void)
{
  int im_1,im_res;

  /*     Image   				*/
  im_1=GET_PLACE3D(TEXT0164);
  im_res=GET_PLACE3D(TEXT0006);


  imx_hole_fill_3d(im_1,im_res);

  show_picture_3d(im_res);

  return;
}

/*******************************************
 ** --  hole_fill_3d_dir() ----------------------------
 **
 **    Fonction remplissage de trous en 3d (d'une image binaire)
 **    A partir de la croissance de region
 ********************************************/
void    hole_fill_3d_dir(void)
{
  int im_1,im_res;

  /*     Image   				*/
  im_1=GET_PLACE3D(TEXT0164);
  im_res=GET_PLACE3D(TEXT0006);


  imx_hole_fill_3d_dir(im_1,im_res);

  show_picture_3d(im_res);

  return;
}


/********************************************
 **   hole_fill_3d_cpc()
 */
/*!
**   \brief Fonction remplissage de trous en 3d (d'une image binaire)
**   coupe par coupe
**   \param  type :
- 1 par rapport X
- 2 par rapport Y
- 3 par rapport Z
- 4 par rapport a tout
********************************************/
int    hole_fill_3d_cpc(int type)
{
  int im_1,im_res;

  /*     Image   				*/
  im_1=GET_PLACE3D(TEXT0164);
  im_res=GET_PLACE3D(TEXT0006);


  imx_hole_fill_3d_cpc(im_1,im_res,type);

  show_picture_3d(im_res);

  return(1);
}



/**************************************************
 ** -- imx_d_euclidean_mapping_3d() ------------------------
 **
 **  Transformee de distance d'une image de contour
 **************************************************/

void imx_d_euclidean_mapping_3d(int im_deb, int im_res_x, int im_res_y, int im_res_z, int im_res_m)
{
  grphic3d *imdeb,*imresx,*imresy,*imresz,*imresm;

  imdeb=ptr_img_3d(im_deb);
  imresx=ptr_img_3d(im_res_x);
  imresy=ptr_img_3d(im_res_y);
  imresz=ptr_img_3d(im_res_z);
  imresm=ptr_img_3d(im_res_m);

  imx_d_euclidean_mapping_3d_p(imdeb,imresx,imresy,imresz,imresm);

  find_sign_3d(imdeb,imresx,imresy,imresz);
  border_processing_3d(imresx,imresy,imresz,imresm);

}


/********************************************
 **    imx_labelcroix_3d()
 */
/*!    Fonction de labelisation 3D
**    A partir de la croissance de region
**	\param im_1 : numero de l'image source
**	\param im_res : numero de l'image resultat
**	\retval 1
********************************************/
int imx_labelcroix_3d(int im_1, int im_res)
{
	grphic3d *im1,*imres;


	im1=ptr_img_3d(im_1);
	imres=ptr_img_3d(im_res);


	imx_labelcroix_3d_p(im1,imres);
#ifdef __GTK_GUI
	imx_SetLabelPalette(im_res);
#endif /*__GTK_GUI*/
	show_picture_3d(im_res);
	return(1);

}


/*****************************************
 ** -- imx_ordonner_labels_3d()------------------------------
 **    Fonction de rangement (tri) des labels 3D
 **    Si par exemple les labels ont des valeurs non
 **    continues dans N (i.e. 1 2 3 4 ... mais 1 3 8 12..)
 **    cette fonction range tout à partir de 1 et avec un
 **    incrément de 1
 **
 **    Rq : voir imx_ordonner_labels_3d_p
 ******************************************/
int imx_ordonner_labels_3d(int im_1, int im_res)
{
  grphic3d *im1,*imres;


  im1=ptr_img_3d(im_1);
  imres=ptr_img_3d(im_res);

  imx_ordonner_labels_3d_p(im1,imres);
  
  return 1;
}

/*******************************************
 ** -- imx_labelcube_3d()------------------------------
 **    Fonction de labelisation 3D
 **    A partir de la croissance de region
 ********************************************/
int imx_labelcube_3d(int im_1, int im_res)
{
  grphic3d *im1,*imres;


  im1=ptr_img_3d(im_1);
  imres=ptr_img_3d(im_res);


  imx_labelcube_3d_p(im1,imres);
  return(1);

}


/*****************************************
 ** -- imx_ctr_binaire_3d()------------------------------
 **    Fonction d'extraction des contours en 3D
 **    Utilise un masque en croix
 ******************************************/
int imx_ctr_binaire_3d(int im_1, int im_res)
{
  grphic3d *im1,*imres;


  im1=ptr_img_3d(im_1);
  imres=ptr_img_3d(im_res);


  imx_ctr_binaire_3d_p(im1,imres);
  return(1);

}

/*****************************************
 ** -- imx_ctr_labele_3d()------------------------------
 **    Fonction d'extraction des contours en 3D
 **    Utilise voisinage en croix
 ******************************************/
int 	imx_ctr_labele_3d(int im_1, int im_res)
{
  grphic3d *im1,*imres;


  im1=ptr_img_3d(im_1);
  imres=ptr_img_3d(im_res);


  imx_ctr_labele_3d_p(im1,imres);
  return(1);

}

/*****************************************
 ** -- imx_ctr_labele_rempl_3d()------------------------------
 **    Fonction d'extraction des contours en 3D
 **    Utilise voisinage en croix
 ******************************************/
int 	imx_ctr_labele_rempl_3d(int im_1, int im_res)
{
  grphic3d *im1,*imres;


  im1=ptr_img_3d(im_1);
  imres=ptr_img_3d(im_res);


  imx_ctr_labele_rempl_3d_p(im1,imres);
  return(1);

}

/*****************************************
 ** -- imx_sync_labels_3d()------------------------------
 **    Fonction de synchronisation entre deux
 **    images de labels
 ******************************************/
int 	imx_sync_labels_3d(int im_src, int im_ref, int im_res)
{
  grphic3d *imsrc,*imref,*imres;


  imsrc=ptr_img_3d(im_src);
  imref=ptr_img_3d(im_ref);
  imres=ptr_img_3d(im_res);


  imx_sync_labels_3d_p(imsrc,imref,imres);
  return(1);

}
/*****************************************
 ** -- imx_sync_labels_sep_3d()------------------------------
 **    Fonction de synchronisation entre deux
 **    images de labels pour la SEP
 **    Auteur : F. Rousseau, 2007
 ******************************************/
int 	imx_sync_labels_sep_3d(int im_src, int im_ref)
{
  grphic3d *imsrc,*imref;


  imsrc=ptr_img_3d(im_src);
  imref=ptr_img_3d(im_ref);

  imx_sync_labels_sep_3d_p(imsrc,imref);
  return(1);

}

/*******************************************
 ** -- imx_growing_cerv_3d()------------------------------
 **    Fonction d'extraction du cerveau en 3d
 **    A partir de la croissance de region
 ********************************************/
int 	imx_growing_cerv_3d(int im_1, int im_res, int x0, int y0, int z0)
{
  grphic3d *im1,*imres;


  im1=ptr_img_3d(im_1);
  imres=ptr_img_3d(im_res);


  imx_growing_cerv_3d_p(im1,imres,x0,y0,z0);
  return(1);

}

/*******************************************
 ** -- imx_extract_cerv_3d()------------------------------
 **    Fonction d'extraction du cerveau en 3d
 **    A partir de la croissance de region
 ********************************************/
int 	imx_extract_cerv_3d(int im_1, int im_res)
{
  grphic3d *im1,*imres;


  im1=ptr_img_3d(im_1);
  imres=ptr_img_3d(im_res);


  imx_extract_cerv_3d_p(im1,imres);
  return(1);

}


/*******************************************
 ** --  imx_hole_fill_3d() ----------------------------
 **
 **    Fonction remplissage de trous en 3d (d'une image binaire)
 **    A partir de la croissance de region
 ********************************************/
int    imx_hole_fill_3d(int im_1, int im_res)
{
  grphic3d  *im1,*imres;

  im1=ptr_img_3d(im_1);
  imres=ptr_img_3d(im_res);

  imx_hole_fill_3d_p(im1,imres);

  return(1);
}


/*******************************************
 ** --  imx_hole_fill_3d_dir() ----------------------------
 **
 **    Fonction remplissage de trous en 3d (d'une image binaire)
 **    A partir de la croissance de region
 ********************************************/
int    imx_hole_fill_3d_dir(int im_1, int im_res)
{
  grphic3d  *im1,*imres;

  im1=ptr_img_3d(im_1);
  imres=ptr_img_3d(im_res);

  imx_hole_fill_3d_dir_p(im1,imres);

  return(1);
}

/*******************************************
 ** --  imx_hole_fill_3d_cpc() ----------------------------
 */
/*!  \ingroup Segmentation
**
**	 \brief   Fonction remplissage de trous en 3d (d'une image binaire)
**    coupe par coupe
**   \param im_1 : numero de l'image source
**   \param im_res : numero de l'image dest
**   \param  type :
- 1 par rapport X
- 2 par rapport Y
- 3 par rapport Z
- 4 par rapport a tout
\retval 1
********************************************/
int    imx_hole_fill_3d_cpc(int im_1, int im_res, int type)
{
  grphic3d  *im1,*imres;

  im1=ptr_img_3d(im_1);
  imres=ptr_img_3d(im_res);

  imx_hole_fill_3d_cpc_p(im1,imres,type);

  return(1);
}


/**************************************************
 ** -- imx_d_euclidean_mapping_3d_p() ------------------------
 **
 **  Transformee de distance d'une image de contour
 **************************************************/
void imx_d_euclidean_mapping_3d_p(grphic3d *imdeb, grphic3d *imresx, grphic3d *imresy, grphic3d *imresz, grphic3d *imresm)
{
  int i,j,k,m;
  int height,width,depth;
  int x[14];
  int y[14];
  int z[14];
  double mod[14];
  double min;


  height=imdeb->height;
  width=imdeb->width;
  depth=imdeb->depth;

  imresx->width=width;
  imresx->height=height;
  imresx->depth=depth;
  imresy->width=width;
  imresy->height=height;
  imresy->depth=depth;
  imresz->width=width;
  imresz->height=height;
  imresz->depth=depth;
  imresm->width=width;
  imresm->height=height;
  imresm->depth=depth;





  /* Mettre les points de contour a zero */

  for(i=0;i<width;i++)
    for(j=0;j<height;j++)
      for(k=0;k<depth;k++)
	{
	  if(imdeb->mri[i][j][k]!=0)
	    {
	      imresx->mri[i][j][k]=0;
	      imresy->mri[i][j][k]=0;
	      imresz->mri[i][j][k]=0;
	    }
	  else if(imdeb->mri[i][j][k]==0)
	    {
	      imresx->mri[i][j][k]=10000;
	      imresy->mri[i][j][k]=10000;
	      imresz->mri[i][j][k]=10000;
	    }
	}


  /********************************** Forward pass ***********************************/

  if(RECA_DEBUG)
    printf("\n\n\n\n Forward pass\n\n");

  for(k=1;k<depth-1;k++)
    {
      for(j=1;j<height-1;j++)
	{
	  for(i=1;i<width-1;i++)
	    {
	      x[0]=imresx->mri[i][j][k];
	      x[1]=imresx->mri[i-1][j][k]+1;
	      x[2]=imresx->mri[i-1][j-1][k]+1;
	      x[3]=imresx->mri[i][j-1][k];
	      x[4]=imresx->mri[i+1][j-1][k]+1;
	      x[5]=imresx->mri[i-1][j-1][k-1]+1;
	      x[6]=imresx->mri[i][j-1][k-1];
	      x[7]=imresx->mri[i+1][j-1][k-1]+1;
	      x[8]=imresx->mri[i-1][j][k-1]+1;
	      x[9]=imresx->mri[i][j][k-1];
	      x[10]=imresx->mri[i+1][j][k-1]+1;
	      x[11]=imresx->mri[i-1][j+1][k-1]+1;
	      x[12]=imresx->mri[i][j+1][k-1];
	      x[13]=imresx->mri[i+1][j+1][k-1]+1;

	      y[0]=imresy->mri[i][j][k];
	      y[1]=imresy->mri[i-1][j][k];
	      y[2]=imresy->mri[i-1][j-1][k]+1;
	      y[3]=imresy->mri[i][j-1][k]+1;
	      y[4]=imresy->mri[i+1][j-1][k]+1;
	      y[5]=imresy->mri[i-1][j-1][k-1]+1;
	      y[6]=imresy->mri[i][j-1][k-1]+1;
	      y[7]=imresy->mri[i+1][j-1][k-1]+1;
	      y[8]=imresy->mri[i-1][j][k-1];
	      y[9]=imresy->mri[i][j][k-1];
	      y[10]=imresy->mri[i+1][j][k-1];
	      y[11]=imresy->mri[i-1][j+1][k-1]+1;
	      y[12]=imresy->mri[i][j+1][k-1]+1;
	      y[13]=imresy->mri[i+1][j+1][k-1]+1;

	      z[0]=imresz->mri[i][j][k];
	      z[1]=imresz->mri[i-1][j][k];
	      z[2]=imresz->mri[i-1][j-1][k];
	      z[3]=imresz->mri[i][j-1][k];
	      z[4]=imresz->mri[i+1][j-1][k];
	      z[5]=imresz->mri[i-1][j-1][k-1]+1;
	      z[6]=imresz->mri[i][j-1][k-1]+1;
	      z[7]=imresz->mri[i+1][j-1][k-1]+1;
	      z[8]=imresz->mri[i-1][j][k-1]+1;
	      z[9]=imresz->mri[i][j][k-1]+1;
	      z[10]=imresz->mri[i+1][j][k-1]+1;
	      z[11]=imresz->mri[i-1][j+1][k-1]+1;
	      z[12]=imresz->mri[i][j+1][k-1]+1;
	      z[13]=imresz->mri[i+1][j+1][k-1]+1;

	      for(m=0;m<14;m++)
		mod[m]=sqrt((double)(x[m]*x[m])+(double)(y[m]*y[m])+(double)(z[m]*z[m]));

	      min=mod[0];
	      imresx->mri[i][j][k]=x[0];
	      imresy->mri[i][j][k]=y[0];
	      imresz->mri[i][j][k]=z[0];

	      for(m=1;m<14;m++)
		if(mod[m]<min)
		  {
		    min=mod[m];
		    imresx->mri[i][j][k]=x[m];
		    imresy->mri[i][j][k]=y[m];
		    imresz->mri[i][j][k]=z[m];
		  }

	    }

	  for(i=width-2;i>=1;i--)
	    {
	      x[0]=imresx->mri[i][j][k];
	      x[1]=imresx->mri[i+1][j][k]+1;

	      y[0]=imresy->mri[i][j][k];
	      y[1]=imresy->mri[i+1][j][k];

	      z[0]=imresz->mri[i][j][k];
	      z[1]=imresz->mri[i+1][j][k];

	      for(m=0;m<2;m++)
		mod[m]=sqrt((double)(x[m]*x[m])+(double)(y[m]*y[m])+(double)(z[m]*z[m]));

	      min=mod[0];
	      imresx->mri[i][j][k]=x[0];
	      imresy->mri[i][j][k]=y[0];
	      imresz->mri[i][j][k]=z[0];

	      for(m=1;m<2;m++)
		if(mod[m]<min)
		  {
		    min=mod[m];
		    imresx->mri[i][j][k]=x[m];
		    imresy->mri[i][j][k]=y[m];
		    imresz->mri[i][j][k]=z[m];
		  }

	    }

	}

      for(j=height-2;j>=1;j--)
	{
	  for(i=width-2;i>=1;i--)
	    {
	      x[0]=imresx->mri[i][j][k];
	      x[1]=imresx->mri[i+1][j][k]+1;
	      x[2]=imresx->mri[i+1][j+1][k]+1;
	      x[3]=imresx->mri[i][j+1][k];
	      x[4]=imresx->mri[i-1][j+1][k]+1;

	      y[0]=imresy->mri[i][j][k];
	      y[1]=imresy->mri[i+1][j][k];
	      y[2]=imresy->mri[i+1][j+1][k]+1;
	      y[3]=imresy->mri[i][j+1][k]+1;
	      y[4]=imresy->mri[i-1][j+1][k]+1;

	      z[0]=imresz->mri[i][j][k];
	      z[1]=imresz->mri[i+1][j][k];
	      z[2]=imresz->mri[i+1][j+1][k];
	      z[3]=imresz->mri[i][j+1][k];
	      z[4]=imresz->mri[i-1][j+1][k];

	      for(m=0;m<5;m++)
		mod[m]=sqrt((double)(x[m]*x[m])+(double)(y[m]*y[m])+(double)(z[m]*z[m]));

	      min=mod[0];
	      imresx->mri[i][j][k]=x[0];
	      imresy->mri[i][j][k]=y[0];
	      imresz->mri[i][j][k]=z[0];

	      for(m=1;m<5;m++)
		if(mod[m]<min)
		  {
		    min=mod[m];
		    imresx->mri[i][j][k]=x[m];
		    imresy->mri[i][j][k]=y[m];
		    imresz->mri[i][j][k]=z[m];
		  }



	    }


	  for(i=1;i<width-1;i++)
	    {
	      x[0]=imresx->mri[i][j][k];
	      x[1]=imresx->mri[i-1][j][k]+1;

	      y[0]=imresy->mri[i][j][k];
	      y[1]=imresy->mri[i-1][j][k];

	      z[0]=imresz->mri[i][j][k];
	      z[1]=imresz->mri[i-1][j][k];

	      for(m=0;m<2;m++)
		mod[m]=sqrt((double)(x[m]*x[m])+(double)(y[m]*y[m])+(double)(z[m]*z[m]));

	      min=mod[0];
	      imresx->mri[i][j][k]=x[0];
	      imresy->mri[i][j][k]=y[0];
	      imresz->mri[i][j][k]=z[0];

	      for(m=1;m<2;m++)
		if(mod[m]<min)
		  {
		    min=mod[m];
		    imresx->mri[i][j][k]=x[m];
		    imresy->mri[i][j][k]=y[m];
		    imresz->mri[i][j][k]=z[m];
		  }

	    }


	}

      if(RECA_DEBUG)
	printf("level %d/%d processed\n",k,depth-2);



    }



  /********************************** Backward pass ***********************************/

  if(RECA_DEBUG)
    printf("\n\n\n\n Backward pass\n\n");

  for(k=depth-2;k>=1;k--)
    {
      for(j=height-2;j>=1;j--)
	{
	  for(i=width-2;i>=1;i--)
	    {
	      x[0]=imresx->mri[i][j][k];
	      x[1]=imresx->mri[i+1][j][k]+1;
	      x[2]=imresx->mri[i+1][j+1][k]+1;
	      x[3]=imresx->mri[i][j+1][k];
	      x[4]=imresx->mri[i-1][j+1][k]+1;
	      x[5]=imresx->mri[i-1][j-1][k+1]+1;
	      x[6]=imresx->mri[i][j-1][k+1];
	      x[7]=imresx->mri[i+1][j-1][k+1]+1;
	      x[8]=imresx->mri[i-1][j][k+1]+1;
	      x[9]=imresx->mri[i][j][k+1];
	      x[10]=imresx->mri[i+1][j][k+1]+1;
	      x[11]=imresx->mri[i-1][j+1][k+1]+1;
	      x[12]=imresx->mri[i][j+1][k+1];
	      x[13]=imresx->mri[i+1][j+1][k+1]+1;

	      y[0]=imresy->mri[i][j][k];
	      y[1]=imresy->mri[i+1][j][k];
	      y[2]=imresy->mri[i+1][j+1][k]+1;
	      y[3]=imresy->mri[i][j+1][k]+1;
	      y[4]=imresy->mri[i-1][j+1][k]+1;
	      y[5]=imresy->mri[i-1][j-1][k+1]+1;
	      y[6]=imresy->mri[i][j-1][k+1]+1;
	      y[7]=imresy->mri[i+1][j-1][k+1]+1;
	      y[8]=imresy->mri[i-1][j][k+1];
	      y[9]=imresy->mri[i][j][k+1];
	      y[10]=imresy->mri[i+1][j][k+1];
	      y[11]=imresy->mri[i-1][j+1][k+1]+1;
	      y[12]=imresy->mri[i][j+1][k+1]+1;
	      y[13]=imresy->mri[i+1][j+1][k+1]+1;

	      z[0]=imresz->mri[i][j][k];
	      z[1]=imresz->mri[i+1][j][k];
	      z[2]=imresz->mri[i+1][j+1][k];
	      z[3]=imresz->mri[i][j+1][k];
	      z[4]=imresz->mri[i-1][j+1][k];
	      z[5]=imresz->mri[i-1][j-1][k+1]+1;
	      z[6]=imresz->mri[i][j-1][k+1]+1;
	      z[7]=imresz->mri[i+1][j-1][k+1]+1;
	      z[8]=imresz->mri[i-1][j][k+1]+1;
	      z[9]=imresz->mri[i][j][k+1]+1;
	      z[10]=imresz->mri[i+1][j][k+1]+1;
	      z[11]=imresz->mri[i-1][j+1][k+1]+1;
	      z[12]=imresz->mri[i][j+1][k+1]+1;
	      z[13]=imresz->mri[i+1][j+1][k+1]+1;


	      for(m=0;m<14;m++)
		mod[m]=sqrt((double)(x[m]*x[m])+(double)(y[m]*y[m])+(double)(z[m]*z[m]));

	      min=mod[0];
	      imresx->mri[i][j][k]=x[0];
	      imresy->mri[i][j][k]=y[0];
	      imresz->mri[i][j][k]=z[0];

	      for(m=1;m<14;m++)
		if(mod[m]<min)
		  {
		    min=mod[m];
		    imresx->mri[i][j][k]=x[m];
		    imresy->mri[i][j][k]=y[m];
		    imresz->mri[i][j][k]=z[m];
		  }

	    }

	  for(i=1;i<width-1;i++)
	    {
	      x[0]=imresx->mri[i][j][k];
	      x[1]=imresx->mri[i-1][j][k]+1;

	      y[0]=imresy->mri[i][j][k];
	      y[1]=imresy->mri[i-1][j][k];

	      z[0]=imresz->mri[i][j][k];
	      z[1]=imresz->mri[i-1][j][k];

	      for(m=0;m<2;m++)
		mod[m]=sqrt((double)(x[m]*x[m])+(double)(y[m]*y[m])+(double)(z[m]*z[m]));

	      min=mod[0];
	      imresx->mri[i][j][k]=x[0];
	      imresy->mri[i][j][k]=y[0];
	      imresz->mri[i][j][k]=z[0];

	      for(m=1;m<2;m++)
		if(mod[m]<min)
		  {
		    min=mod[m];
		    imresx->mri[i][j][k]=x[m];
		    imresy->mri[i][j][k]=y[m];
		    imresz->mri[i][j][k]=z[m];
		  }

	    }

	}

      for(j=1;j<height-1;j++)
	{
	  for(i=1;i<width-1;i++)
	    {
	      x[0]=imresx->mri[i][j][k];
	      x[1]=imresx->mri[i-1][j][k]+1;
	      x[2]=imresx->mri[i-1][j-1][k]+1;
	      x[3]=imresx->mri[i][j-1][k];
	      x[4]=imresx->mri[i+1][j-1][k]+1;

	      y[0]=imresy->mri[i][j][k];
	      y[1]=imresy->mri[i-1][j][k];
	      y[2]=imresy->mri[i-1][j-1][k]+1;
	      y[3]=imresy->mri[i][j-1][k]+1;
	      y[4]=imresy->mri[i+1][j-1][k]+1;

	      z[0]=imresz->mri[i][j][k];
	      z[1]=imresz->mri[i-1][j][k];
	      z[2]=imresz->mri[i-1][j-1][k];
	      z[3]=imresz->mri[i][j-1][k];
	      z[4]=imresz->mri[i+1][j-1][k];

	      for(m=0;m<5;m++)
		mod[m]=sqrt((double)(x[m]*x[m])+(double)(y[m]*y[m])+(double)(z[m]*z[m]));

	      min=mod[0];
	      imresx->mri[i][j][k]=x[0];
	      imresy->mri[i][j][k]=y[0];
	      imresz->mri[i][j][k]=z[0];

	      for(m=1;m<5;m++)
		if(mod[m]<min)
		  {
		    min=mod[m];
		    imresx->mri[i][j][k]=x[m];
		    imresy->mri[i][j][k]=y[m];
		    imresz->mri[i][j][k]=z[m];
		  }



	    }


	  for(i=width-2;i>=1;i--)
	    {
	      x[0]=imresx->mri[i][j][k];
	      x[1]=imresx->mri[i+1][j][k]+1;

	      y[0]=imresy->mri[i][j][k];
	      y[1]=imresy->mri[i+1][j][k];

	      z[0]=imresz->mri[i][j][k];
	      z[1]=imresz->mri[i+1][j][k];

	      for(m=0;m<2;m++)
		mod[m]=sqrt((double)(x[m]*x[m])+(double)(y[m]*y[m])+(double)(z[m]*z[m]));

	      min=mod[0];
	      imresx->mri[i][j][k]=x[0];
	      imresy->mri[i][j][k]=y[0];
	      imresz->mri[i][j][k]=z[0];

	      for(m=1;m<2;m++)
		if(mod[m]<min)
		  {
		    min=mod[m];
		    imresx->mri[i][j][k]=x[m];
		    imresy->mri[i][j][k]=y[m];
		    imresz->mri[i][j][k]=z[m];
		  }

	      imresm->mri[i][j][k]=(int)(0.5+min);

	    }


	}

      if(RECA_DEBUG)
	printf("level %d/%d processed\n",k,depth-2);



    }


} /* end of function */

/************************************************
 **
 *********************************************/
void find_sign_3d(grphic3d *imdeb, grphic3d *imresx, grphic3d *imresy, grphic3d *imresz)
{
  int i,j,k;
  int x,y,z;
  int height,width,depth;
  int flag=0;

  width=imdeb->width;
  height=imdeb->height;
  depth=imdeb->depth;

  for(k=1;k<depth-1;k++)
    for(i=1;i<width-1;i++)
      for(j=1;j<height-1;j++)
	{
	  x=imresx->mri[i][j][k];
	  y=imresy->mri[i][j][k];
	  z=imresz->mri[i][j][k];
	  flag=0;


	  if(i-x>0 && j-y>0 && k-z>0 && flag==0)
	    if(imdeb->mri[i-x][j-y][k-z]!=0)
	      {
		imresx->mri[i][j][k]=-x;
		imresy->mri[i][j][k]=-y;
		imresz->mri[i][j][k]=-z;
		flag=1;
	      }

	  if(i-x>0 && j-y>0 && k+z<depth-1 && flag==0)
	    if(imdeb->mri[i-x][j-y][k+z]!=0)
	      {
		imresx->mri[i][j][k]=-x;
		imresy->mri[i][j][k]=-y;
		imresz->mri[i][j][k]=z;
		flag=1;
	      }

	  if(i-x>0 && j+y<height-1 && k-z>0 && flag==0)
	    if(imdeb->mri[i-x][j+y][k-z]!=0)
	      {
		imresx->mri[i][j][k]=-x;
		imresy->mri[i][j][k]=y;
		imresz->mri[i][j][k]=-z;
		flag=1;
	      }

	  if(i-x>0 && j+y<height-1 && k+z<depth-1 && flag==0)
	    if(imdeb->mri[i-x][j+y][k+z]!=0)
	      {
		imresx->mri[i][j][k]=-x;
		imresy->mri[i][j][k]=y;
		imresz->mri[i][j][k]=z;
		flag=1;
	      }

	  if(i+x<width-1 && j-y>0 && k-z>0 && flag==0)
	    if(imdeb->mri[i+x][j-y][k-z]!=0)
	      {
		imresx->mri[i][j][k]=x;
		imresy->mri[i][j][k]=-y;
		imresz->mri[i][j][k]=-z;
		flag=1;
	      }

	  if(i+x<width-1 && j-y>0 && k+z<depth-1 && flag==0)
	    if(imdeb->mri[i+x][j-y][k+z]!=0)
	      {
		imresx->mri[i][j][k]=x;
		imresy->mri[i][j][k]=-y;
		imresz->mri[i][j][k]=z;
		flag=1;
	      }

	  if(i+x<width-1 && j+y<height-1 && k-z>0 && flag==0)
	    if(imdeb->mri[i+x][j+y][k-z]!=0)
	      {
		imresx->mri[i][j][k]=x;
		imresy->mri[i][j][k]=y;
		imresz->mri[i][j][k]=-z;
		flag=1;
	      }

	  if(i+x<width-1 && j+y<height-1 && k+z<depth-1 && flag==0)
	    if(imdeb->mri[i+x][j+y][k+z]!=0)
	      {
		imresx->mri[i][j][k]=x;
		imresy->mri[i][j][k]=y;
		imresz->mri[i][j][k]=z;
		flag=1;
	      }

	}

}

/******************************************************************************/

void border_processing_3d(grphic3d *imresx, grphic3d *imresy, grphic3d *imresz, grphic3d *imresm)
{
  int i,j,k;
  int width,height,depth;

  width=imresx->width;
  height=imresx->height;
  depth=imresx->depth;


  for(i=0;i<width;i++)
    for(j=0;j<height;j++)
      {
	imresx->mri[i][j][0]=imresx->mri[i][j][1];
	imresy->mri[i][j][0]=imresy->mri[i][j][1];
	imresz->mri[i][j][0]=imresz->mri[i][j][1];
	imresm->mri[i][j][0]=imresm->mri[i][j][1];
	imresx->mri[i][j][depth-1]=imresx->mri[i][j][depth-2];
	imresy->mri[i][j][depth-1]=imresy->mri[i][j][depth-2];
	imresz->mri[i][j][depth-1]=imresz->mri[i][j][depth-2];
	imresm->mri[i][j][depth-1]=imresm->mri[i][j][depth-2];
      }

  for(i=0;i<width;i++)
    for(k=0;k<depth;k++)
      {
	imresx->mri[i][0][k]=imresx->mri[i][1][k];
	imresy->mri[i][0][k]=imresy->mri[i][1][k];
	imresz->mri[i][0][k]=imresz->mri[i][1][k];
	imresm->mri[i][0][k]=imresm->mri[i][1][k];
	imresx->mri[i][height-1][k]=imresx->mri[i][height-2][k];
	imresy->mri[i][height-1][k]=imresy->mri[i][height-2][k];
	imresz->mri[i][height-1][k]=imresz->mri[i][height-2][k];
	imresm->mri[i][height-1][k]=imresm->mri[i][height-2][k];
      }

  for(j=0;j<height;j++)
    for(k=0;k<depth;k++)
      {
	imresx->mri[0][j][k]=imresx->mri[1][j][k];
	imresy->mri[0][j][k]=imresy->mri[1][j][k];
	imresz->mri[0][j][k]=imresz->mri[1][j][k];
	imresm->mri[0][j][k]=imresm->mri[1][j][k];
	imresx->mri[width-1][j][k]=imresx->mri[width-2][j][k];
	imresy->mri[width-1][j][k]=imresy->mri[width-2][j][k];
	imresz->mri[width-1][j][k]=imresz->mri[width-2][j][k];
	imresm->mri[width-1][j][k]=imresm->mri[width-2][j][k];
      }


}

/**************************************************
 ** -- chamfer_distance_3d() ------------------------
 **
 **  Transformee de distance d'une image de contour
 **************************************************/

void chamfer_distance_3d(void)
{
  int im_deb,im_res;

  im_deb=GET_PLACE3D("Original Image ?");
  im_res=GET_PLACE3D("Module ?");

  imx_chamfer_distance_3d(im_deb,im_res);

  imx_copie_param_3d(im_deb,im_res);
  imx_iniparaimg_3d(im_res);
  show_picture_3d(im_res);

  return;

}


/**************************************************
 ** -- imx_chamfer_distance_3d() ------------------------
 **
 **  Transformee de distance d'une image de contour
 **************************************************/

void imx_chamfer_distance_3d(int im_deb, int im_res)
{
  grphic3d *imdeb,*imres;


  imdeb=ptr_img_3d(im_deb);
  imres=ptr_img_3d(im_res);

  imx_chamfer_distance_3d_p(imdeb,imres);

}

/*********************************************************************************/

/**************************************************
 ** -- imx_chamfer_distance_3d_p() ------------------------
 **
 **  Transformee de distance d'une image de contour
 **************************************************/

void imx_chamfer_distance_3d_p(grphic3d *imdeb, grphic3d *imres)
{
  int i,j,k,m;
  int dist[14],min;
  int width,height,depth;
  int d1=3,d2=4,d3=5;

  width=imdeb->width;
  height=imdeb->height;
  depth=imdeb->depth;

imres->rcoeff=1.0;
imres->icomp=0;

  /*************************** TEST ********************************

  for(i=0;i<width;i++)
    for(j=0;j<height;j++)
      for(k=0;k<depth;k++)
           imdeb->mri[i][j][k]=0;
  imdeb->mri[width/2][height/2][depth/2]=1;


  *************************** TEST ********************************/

  /* Mettre les points de contour a zero */

  for(i=0;i<width;i++)
    for(j=0;j<height;j++)
      for(k=0;k<depth;k++)
	{
	  if(imdeb->mri[i][j][k]!=0)
	    imres->mri[i][j][k]=0;
	  else
	    imres->mri[i][j][k]=1000;
	}


  /************* Forward pass **************/

  printf("------ Forward pass ------\n");


  for(k=1;k<depth-1;k++)
    {
      for(i=1;i<width-1;i++)
        for(j=1;j<height-1;j++)
	  {
            dist[0]=imres->mri[i][j][k];
            dist[1]=imres->mri[i][j-1][k]+d1;
            dist[2]=imres->mri[i-1][j-1][k]+d2;
            dist[3]=imres->mri[i-1][j][k]+d1;
            dist[4]=imres->mri[i-1][j+1][k]+d2;
            dist[5]=imres->mri[i-1][j-1][k-1]+d3;
            dist[6]=imres->mri[i-1][j][k-1]+d2;
            dist[7]=imres->mri[i-1][j+1][k-1]+d3;
            dist[8]=imres->mri[i][j-1][k-1]+d2;
            dist[9]=imres->mri[i][j][k-1]+d1;
            dist[10]=imres->mri[i][j+1][k-1]+d2;
            dist[11]=imres->mri[i+1][j-1][k-1]+d3;
            dist[12]=imres->mri[i+1][j][k-1]+d2;
            dist[13]=imres->mri[i+1][j+1][k-1]+d3;


	    min=dist[0];
	    imres->mri[i][j][k]=dist[0];

	    for(m=1;m<14;m++)
	      if(dist[m]<min)
		{
		  min=dist[m];
		  imres->mri[i][j][k]=dist[m];
		}

	  }


      printf("level %d/%d processed \r",k,depth);


    }



  /************* Backward pass **************/

  printf("------ Backward pass ------\n");

  for(k=depth-2;k>0;k--)
    {
      for(i=width-2;i>0;i--)
        for(j=height-2;j>0;j--)
	  {
            dist[0]=imres->mri[i][j][k];
            dist[1]=imres->mri[i][j+1][k]+d1;
            dist[2]=imres->mri[i+1][j+1][k]+d2;
            dist[3]=imres->mri[i+1][j][k]+d1;
            dist[4]=imres->mri[i+1][j-1][k]+d2;
            dist[5]=imres->mri[i-1][j-1][k+1]+d3;
            dist[6]=imres->mri[i-1][j][k+1]+d2;
            dist[7]=imres->mri[i-1][j+1][k+1]+d3;
            dist[8]=imres->mri[i][j-1][k+1]+d2;
            dist[9]=imres->mri[i][j][k+1]+d1;
            dist[10]=imres->mri[i][j+1][k+1]+d2;
            dist[11]=imres->mri[i+1][j-1][k+1]+d3;
            dist[12]=imres->mri[i+1][j][k+1]+d2;
            dist[13]=imres->mri[i+1][j+1][k+1]+d3;


	    min=dist[0];
	    imres->mri[i][j][k]=dist[0];

	    for(m=1;m<14;m++)
	      if(dist[m]<min)
		{
		  min=dist[m];
		  imres->mri[i][j][k]=dist[m];
		}

	  }


      printf("level %d/%d processed \r",k,depth);


    }

  for(i=0;i<width;i++)
    for(j=0;j<height;j++)
      {
	imres->mri[i][j][0]=imres->mri[i][j][1];
	imres->mri[i][j][depth-1]=imres->mri[i][j][depth-2];
      }

  for(i=0;i<width;i++)
    for(k=0;k<depth;k++)
      {
	imres->mri[i][0][k]=imres->mri[i][1][k];
	imres->mri[i][height-1][k]=imres->mri[i][height-2][k];
      }

  for(j=0;j<height;j++)
    for(k=0;k<depth;k++)
      {
	imres->mri[0][j][k]=imres->mri[1][j][k];
	imres->mri[width-1][j][k]=imres->mri[width-2][j][k];
      }

  for(i=0;i<width;i++)
    for(j=0;j<height;j++)
      for(k=0;k<depth;k++)
	imres->mri[i][j][k]/=3;

}

/*********************************************************************************/



/********************************************
 **  fct_grow_label_3d(i0,j0,k0,i,j,k,im1,tab,n)
 */
/*!	Fonction pour tester la connexite
**
**  \param   i0 : point de depart ? a verifier ?
**  \param   j0 : point de depart ? a verifier ?
**  \param   k0 : point de depart ? a verifier ?
**  \param   i : pixel a controler
**  \param   j : pixel a controler
**  \param   k : pixel a controler
**  \param   im1 : pointeur de l'image
**  \param   tab : Tableau contenant le parametres a passer
**  \param   n : nombre de valeurs dans tab
**  \return 1 si reussite, 0 sinon
*********************************************/
int	fct_grow_label_3d(int i0, int j0, int k0, int i, int j, int k, grphic3d *im1, float *tab, int n)
{
  long y,t;

  t=(long)(tab[0]);
  y=im1->mri[i][j][k];
  if(y==t) return(1);
  return(0);

}

/*! \ingroup ToolBoxImage3d @{ */
/********************************************
 **   imx_labelcroix_3d_p()
 */
/*!
**
**    Fonction de labelisation 3D
**    a partir de la croissance de region
** \param im1 : image source
** \param imres : image destination (E/S)
** \retval 1
********************************************/
int     imx_labelcroix_3d_p(grphic3d *im1, grphic3d *imres)
{
  grphic3d *imtemp;
  int x0,y0,z0;      /*point de depart de la croissance*/
  float tab[1];
  long n;
  int w,h,d;

  //ATTENTION - RISQUE DE PLANTAGE
  //le nombre de label ne peut dépasser 65536 !!!

  /* creation de l'image temporaire imtemp */
  imtemp=cr_grphic3d(im1);

  n=1;
  w = im1->width;
  h = im1->height;
  d = im1->depth;

  for (x0=1;x0<w-1;x0++)
    for (y0=1;y0<h-1;y0++)
      for (z0=1;z0<d-1;z0++)
	{
	  if ((im1->mri[x0][y0][z0]>0) && (imtemp->mri[x0][y0][z0]==0))
	    {
	      tab[0]=(float)(im1->mri[x0][y0][z0]);
	      growing_3d_p(x0,y0,z0,im1,imtemp,fct_grow_label_3d,tab,1,n,0);
	      n=n+1;
	    }
	}

  imx_copie_3d_p(imtemp,imres);

  /* liberation de la memoire utilisee par les images temporaires*/
  free_grphic3d (imtemp);

  imres->icomp=0;
  imres->rcoeff=1;
  imres->max_pixel=n-1;
  imres->min_pixel=0;
  imres->cutoff_max=n-1;
  imres->cutoff_min=0;


  return(1);
}
/*! @} */


/*****************************************
 ** -- imx_ordonner_labels_3d_p()------------------------------
 **    Fonction de rangement (tri) des labels 3D
 **    Si par exemple les labels ont des valeurs non
 **    continues dans N (i.e. 1 2 3 4 ... mais 1 3 8 12..)
 **    cette fonction range tout à partir de 1 et avec un
 **    incrément de 1
 **    
 **    Remarque : les labels doivent être de valeur positive
 **   
 **    Principe : 
 **             1 - enregistre dans un tableau toutes les valeurs des labels (via un histogramme)
 **		2 - les valeurs non nulles sont rangées dans un tableau par ordre croissant
 **		3 - il ne reste plus qu'à affecter les valeurs rangées aux voxels de l'image résultante
 ******************************************/
int imx_ordonner_labels_3d_p(grphic3d *im1, grphic3d *imres)
{
 grphic3d *imtemp;
 int valtmp; 
 int nbval = 0;
 int i, j, k; 
 long int *tabOrig;
 int TailleTabOrig = im1->max_pixel + 1;
 int *tabSortedLabels;
  int counter = 0;

 // Calcule le max et le min de im1
 imx_inimaxminpixel_3d_p(im1);
 
 if(im1->min_pixel < 0)
   {
     printf("Le minimum en intensite de l'image d'entree est negatif. Operation d'ordonnacement annulee\n");
     exit(-1);
   }
 

 /************/
 // ETAPE 1 //
 /***********/
 // Création d'un tableau (en fait histogramme des labels)
 // où toutes les valeurs non nulles correspondent à un label


 
 tabOrig = (long int *)malloc(sizeof(long int)*TailleTabOrig);
 memset(tabOrig, 0, im1->max_pixel*sizeof(long int));
 
 for ( i = 0 ; i < im1->width ; i++ )
   { 
     for (j = 0; j < im1->height ; j++)
       {
	 for (k = 0 ; k< im1->depth ; k++) 
	   {
	     tabOrig[im1->mri[i][j][k]]++;
	   }
       }
   }
 
   
 /************/
 // ETAPE 2 //
 /***********/


 //int TailleTabSortedLabels = im1->max_pixel + 1;
 tabSortedLabels = (int *)malloc(sizeof(int)*TailleTabOrig);
 memset(tabSortedLabels, 0, TailleTabOrig*sizeof(int));
 

 for(i=0; i!= TailleTabOrig; i++)
   {
     if(tabOrig[i] != 0)
       {
	 tabSortedLabels[i] = counter;
	 counter++;
       }
   }
  
 // On libère tout de suite ce tableau à présent inutile
 free(tabOrig);
 
 /************/
 // ETAPE 3 //
 /***********/ 
 // On affecte les valeurs rangées dans la nouvelle image
 	for ( i = 0 ; i < im1->width ; i++ )
	{ 
		for (j = 0; j < im1->height ; j++)
		{
			for (k = 0 ; k< im1->depth ; k++) 
			{
			 imres->mri[i][j][k] = tabSortedLabels[im1->mri[i][j][k]];
			}
		}
	} 

 free(tabSortedLabels); 
 
 
 imres->icomp = 0;
 imres->rcoeff = 1;
 imres->max_pixel = counter - 1;
 imres->min_pixel = 0;
 imres->cutoff_max = counter - 1;
 imres->cutoff_min = 0;
  
 return 1;
}



/*******************************************
 ** --  imx_labelcube_3d_p() ----------------------------
 **
 **    Fonction de labelisation 3D
 **    A partir de la croissance de region
 ********************************************/
int     imx_labelcube_3d_p(grphic3d *im1, grphic3d *imres)
{
  grphic3d *imtemp;
  int x0,y0,z0;      /*point de depart de la croissance*/
  float tab[1];
  long n;
  int w,h,d;

  /* creation de l'image temporaire imtemp */
  imtemp=cr_grphic3d(im1);

  n=1;
  w = im1->width;
  h = im1->height;
  d = im1->depth;

  for (x0=1;x0<w-1;x0++)
    for (y0=1;y0<h-1;y0++)
      for (z0=1;z0<d-1;z0++)
	{
	  if ((im1->mri[x0][y0][z0]>0) && (imtemp->mri[x0][y0][z0]==0))
	    {
	      tab[0]=(float)(im1->mri[x0][y0][z0]);
	      growing_3d_p(x0,y0,z0,im1,imtemp,fct_grow_label_3d,tab,1,n,1);
	      n=n+1;
	    }
	}
  imx_copie_3d_p(imtemp,imres);

  /* liberation de la memoire utilisee par les images temporaires*/
  free_grphic3d (imtemp);

  imres->icomp=0;
  imres->rcoeff=1;
  imres->max_pixel=n-1;
  imres->min_pixel=0;
  imres->cutoff_max=n-1;
  imres->cutoff_min=0;


  return(1);
}


/*****************************************
 ** -- imx_ctr_binaire_3d_p()------------------------------
 **    Fonction d'extraction des contours en 3D
 **    Utilise un masque en croix
 ******************************************/
int     imx_ctr_binaire_3d_p(grphic3d *im1, grphic3d *imres)
{
  grphic3d *imtemp;
  int i,j,k,l;
  int w,h,d;
  XPoints_3d tabvois[6];   /*Tableau definissant la position des voisins dans le cas de la croix 3d*/
  long resul;

  /*Definition du tableau tabvois dans le cas de la croix 3d*/
  tabvois[0].x=-1;tabvois[0].y=0;tabvois[0].z=0;
  tabvois[1].x=1;tabvois[1].y=0;tabvois[1].z=0;
  tabvois[2].x=0;tabvois[2].y=-1;tabvois[2].z=0;
  tabvois[3].x=0;tabvois[3].y=1;tabvois[3].z=0;
  tabvois[4].x=0;tabvois[4].y=0;tabvois[4].z=-1;
  tabvois[5].x=0;tabvois[5].y=0;tabvois[5].z=1;

  /* creation de l'image temporaire imtemp */
  imtemp=cr_grphic3d(im1);
  w = im1->width;
  h = im1->height;
  d = im1->depth;

  for (i=1;i<w-1;i++)
    for (j=1;j<h-1;j++)
      for (k=1;k<d-1;k++)
        if ((im1->mri[i][j][k])>0)
	  {
	    resul=im1->mri[i][j][k];
	    for (l=0;l<6;l++) resul=resul && im1->mri[i+tabvois[l].x][j+tabvois[l].y][k+tabvois[l].z];
	    if (resul) imtemp->mri[i][j][k]=0;
            else imtemp->mri[i][j][k]=1;
	  }



  imx_copie_3d_p(imtemp,imres);

  /* liberation de la memoire utilisee par les images temporaires*/
  free_grphic3d (imtemp);

  imres->icomp=0;
  imres->rcoeff=1;
  imx_inimaxminpixel_3d_p(imres);


  return(1);
}

/*****************************************
 ** -- imx_ctr_labele_3d_p()------------------------------
 **    Fonction d'extraction des contours en 3D sur une image labele
 **    Utilise voisinage en croix
 ******************************************/
int     imx_ctr_labele_3d_p(grphic3d *im1, grphic3d *imres)
{
  int i,j,k,s,r,t,result;
  int w,h,d;
  w = im1->width;
  h = im1->height;
  d = im1->depth;

  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      for (k=0;k<d;k++)
	imres->mri[i][j][k]=0;

  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      for (k=0;k<d;k++)
        if ((im1->mri[i][j][k])>0){
	  result=FALSE;
	  for (s= -1; s<=1; s++)
	    for (t= -1; t<=1; t++)
	      for (r= -1; r<=1; r++)
		if ((image_pixel_3d(i+s,j+t,k+r,im1)==TRUE) &&
		    ((s==0 && t==0 && r!=0) || (s==0 && r==0 && t!=0) || (t==0 && r==0 && s!=0))){
		  result=result || ((im1->mri[i+s][j+t][k+r] != im1->mri[i][j][k]) && (imres->mri[i+s][j+t][k+r] != 1));
		  if (result) break;
		}
	  if (result) imres->mri[i][j][k]=1;
	}


  imres->icomp=0;
  imres->rcoeff=1;
  imx_inimaxminpixel_3d_p(imres);
  return(1);
}


/*****************************************
 ** -- imx_ctr_labele_rempl_3d_p()------------------------------
 **    Fonction d'extraction des contours en 3D sur une image labele
 **    Utilise voisinage en croix
 ******************************************/
int     imx_ctr_labele_rempl_3d_p(grphic3d *im1, grphic3d *imres)
{
  int i,j,k,result;
  int w,h,d;
  int o;
  grphic3d *imtemp;

  w = im1->width;
  h = im1->height;
  d = im1->depth;

  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      for (k=0;k<d;k++)
	imres->mri[i][j][k]=0;

  if (im1->max_pixel >255 ){
    PUT_ERROR("Trop d'objets dans l'image");
    return 0;
  }

  for (o=1; o<=im1->max_pixel; o++) {
    printf("Objet = %d /%ld \r",o,im1->max_pixel);
    fflush(stdout);
    imtemp = cr_grphic3d(im1);
    imx_copie_3d_p(im1, imtemp);
    remplit_objet(imtemp, o);


    for (i=1;i<w-1;i++)
      for (j=1;j<h-1;j++)
        for (k=1;k<d-1;k++)
          if ((imtemp->mri[i][j][k])>0)
	    {
	      result=imtemp->mri[i][j][k];
	      result=result && imtemp->mri[i-1][j+0][k+0] && imtemp->mri[i+1][j+0][k+0] && imtemp->mri[i+0][j-1][k+0]
		&& imtemp->mri[i+0][j+1][k+0] && imtemp->mri[i+0][j+0][k-1] && imtemp->mri[i+0][j+0][k+1];
	      if (!result) imres->mri[i][j][k]=o;
	    }

  }

  imx_copie_param_3d_p(im1,imres);
  imres->icomp=0;
  imres->rcoeff=1;
  imx_inimaxminpixel_3d_p(imres);
  return(1);
}

/*****************************************
 ** -- imx_sync_labels_3d_p()------------------------------
 **    Fonction de synchronisation entre deux
 **    images de labels
 **
 **    Principe : 
 **             1 - enregistre dans un tableau de deux dimensions toutes les valeurs des labels
 **                 de l'image source et les valeurs correspondantes de l'image reference (via un histogramme)
 **		2 - enregistre dans un tableau les valeurs des labels reference qui correspondent le mieux aux
 **                 labels source
 **		3 - renumeroter les labels source
 ******************************************/
int     imx_sync_labels_3d_p(grphic3d *imsrc, grphic3d *imref, grphic3d *imres)
{
  int i,j,k;
  int w,h,d;
  int srcLabel, refLabel;
  float overlap, oldOverlap;
  int maxSrc = (int)imsrc->max_pixel+1;
  int maxRef = (int)imref->max_pixel+1;

  int **tabLabels;
  long int *tabLabelSize;
  int *tabNewLabel;


  // arrete si les valeurs ne sont pas entre 0 et 255
  if (imsrc->min_pixel < 0){
    PUT_ERROR("Labels doivent etre positives");
    return 0;
  }
  if (imref->min_pixel < 0){
    PUT_ERROR("Labels doivent etre positives");
    return 0;
  }
  if (imsrc->max_pixel > 255){
    PUT_ERROR("Trop d'objets dans l'image source");
    return 0;
  }
  if (imref->max_pixel > 255){
    PUT_ERROR("Trop d'objets dans l'image reference");
    return 0;
  }


  // allocation de memoire et initialisation des tableaux
  tabLabelSize = (long int *)malloc(sizeof(long int)*maxSrc);
  memset(tabLabelSize, 0, sizeof(long int)*maxSrc);

  tabNewLabel = (int *)malloc(sizeof(int)*maxSrc);
  memset(tabNewLabel, 0, sizeof(int)*maxSrc);

  tabLabels = (int **)malloc(sizeof(int *)*maxSrc);
  for (srcLabel=0; srcLabel<maxSrc; srcLabel++)
    {
      tabLabels[srcLabel] = (int *)malloc(sizeof(int)*maxRef);
      memset(tabLabels[srcLabel], 0, sizeof(int)*maxRef);
    }

  w = (imsrc->width < imref->width ? imsrc->width : imref->width);
  h = (imsrc->height < imref->height ? imsrc->height : imref->height);
  d = (imsrc->depth < imref->depth ? imsrc->depth : imref->depth);

  // enregistrement des histograms des labels dans le tableau 2D
  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      for (k=0;k<d;k++)
	{
	  tabLabels[imsrc->mri[i][j][k]][imref->mri[i][j][k]]++;
	  tabLabelSize[imsrc->mri[i][j][k]]++;
	}
 
  // tester les correspondances entre les labels source et reference
  for (srcLabel=0; srcLabel<maxSrc; srcLabel++)
    {
      oldOverlap = 0;
      for (refLabel=0; refLabel<maxRef; refLabel++)
	{
	  overlap = (float)tabLabels[srcLabel][refLabel]/(float)tabLabelSize[srcLabel];
	  if (overlap > oldOverlap)
	    {
	      tabNewLabel[srcLabel] = refLabel;
	      oldOverlap = overlap;
	    }
	}
      printf("Label %d -> %d (correspondance %3.1f%%)\n", srcLabel, tabNewLabel[srcLabel], oldOverlap*100);
    }

  // renumeroter les labels source
  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      for (k=0;k<d;k++)
	imres->mri[i][j][k] = (int)tabNewLabel[imsrc->mri[i][j][k]];

  imx_copie_param_3d_p(imsrc,imres);
  imres->icomp=0;
  imres->rcoeff=1;
  imx_inimaxminpixel_3d_p(imres);
  return(1);
}

/*****************************************
 ** -- imx_sync_labels_sep_3d_p()------------------------------
 **    Fonction de comptage d'apparition et disparition de composantes connexes
 **    Nécessite des images labélisées en entrée
 **    Principe : 
 **             1 - Compare les labels de l'image source à l'aide des labels de l'image de reference
 **		2 - Permet de détecter des nouveaux labels (ie de nouvelles lésions) et des disparitions de labels
 ******************************************/
int     imx_sync_labels_sep_3d_p(grphic3d *imsrc, grphic3d *imref)
{
  int i,j,k;
  int w,h,d;
  int srcLabel, refLabel;
  float overlap, oldOverlap;
  float dx,dy,dz;
  float d3;
  int maxSrc = (int)imsrc->max_pixel+1;
  int maxRef = (int)imref->max_pixel+1;
  int CurrentNewLabel = maxRef;
    int NbApparitionLesions = 0;
    int NbDisparitionLesions = 0;

    int **tabLabels;
    long int *tabLabelSize;
    long int *tabLabelSizeRef;
    int *tabNewLabel;
    long int VolumeLesionDiff = 0;
	  short dect = 0;
	  	  float NbVoxelRef = 0.0;
	  	    float NbVoxelSrc = 0.0;
	  	    long int DiffVoxel = 0;

  /* On suppose que les labels vont de 0 à max label */
  if (imsrc->min_pixel < 0){
    PUT_ERROR("Labels doivent etre positifs");
    return 0;
  }
  if (imref->min_pixel < 0){
    PUT_ERROR("Labels doivent etre positifs");
    return 0;
  }


  printf("Nombre de labels dans l'image source (image 2) %d\n", maxSrc);
  printf("Nombre de labels dans l'image reference (image 1) %d\n", maxRef);
  printf("Le calcul du volume des labels est effectué sur la résolution de l'image de référence\n");
  printf("Definition d'une apparition ou disparition de label (lesion) : \n");
  printf("-> label dont l'overlap avec le fond (label 0) vaut 1\n");




  // allocation de memoire et initialisation des tableaux
  tabLabelSize = (long int *)malloc(sizeof(long int)*maxSrc);
  memset(tabLabelSize, 0, sizeof(long int)*maxSrc);
  tabLabelSizeRef = (long int *)malloc(sizeof(long int)*maxRef);
  memset(tabLabelSizeRef, 0, sizeof(long int)*maxRef);

  tabNewLabel = (int *)malloc(sizeof(int)*maxSrc);
  memset(tabNewLabel, 0, sizeof(int)*maxSrc);

  tabLabels = (int **)malloc(sizeof(int *)*maxSrc);
  for (srcLabel=0; srcLabel<maxSrc; srcLabel++)
    {
      tabLabels[srcLabel] = (int *)malloc(sizeof(int)*maxRef);
      memset(tabLabels[srcLabel], 0, sizeof(int)*maxRef);
    }

  w = (imsrc->width < imref->width ? imsrc->width : imref->width);
  h = (imsrc->height < imref->height ? imsrc->height : imref->height);
  d = (imsrc->depth < imref->depth ? imsrc->depth : imref->depth);

  dx = imref->dx;
  dy = imref->dy;
  dz = imref->dz;  
  d3 = dx*dy*dz;
  printf("Resolution de l'image de reference : %f %f %f \n",dx,dy,dz);


  if(dx != imsrc->dx)
	printf("Problème : l'image source n'a pas la même résolution en x que l'image de référence ! \n");
  if(dy != imsrc->dy)
	printf("Problème : l'image source n'a pas la même résolution en y que l'image de référence ! \n");
  if(dz != imsrc->dz)
	printf("Problème : l'image source n'a pas la même résolution en z que l'image de référence ! \n");


  // enregistrement des histograms des labels dans le tableau 2D
  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      for (k=0;k<d;k++)
	{
	  tabLabels[imsrc->mri[i][j][k]][imref->mri[i][j][k]]++;
	  tabLabelSize[imsrc->mri[i][j][k]]++;
	  tabLabelSizeRef[imref->mri[i][j][k]]++;
	}
 
  // tester les correspondances entre les labels source et reference
  // on omet le label 0 (fond)
  for (srcLabel=1; srcLabel<maxSrc; srcLabel++)
    {
      oldOverlap = 0;
      for (refLabel=1; refLabel<maxRef; refLabel++)
	{
	  overlap = (float)tabLabels[srcLabel][refLabel]/(float)tabLabelSize[srcLabel];
	  if (overlap > oldOverlap)
	    {
	      tabNewLabel[srcLabel] = refLabel;
	      oldOverlap = overlap;
	    }
	}
      if(oldOverlap==0)
	{
	  //tabNewLabel[srcLabel] = CurrentNewLabel;
	  //CurrentNewLabel++;
	  printf("Nouvelle Lésion : Label %d (image source) de %d voxels", srcLabel, tabLabelSize[srcLabel]);
	  VolumeLesionDiff += tabLabelSize[srcLabel];
	  NbApparitionLesions++;

	  //Connaitre la position d'un point de ce label 

	  i=0;
	  j=0;
	  k=0;
	  while(i<w && dect!=1)
	    {
	      j = 0;
	      while(j<h && dect!=1)
		{
		  k=0;
		  while(k<d && dect!=1)
		    {
		      if(imsrc->mri[i][j][k] == srcLabel)
			{
			  printf("  - Position : %d %d %d \n", i, j, k);
			  dect=1;
			}
		      k++;
		    }
		  j++;
		}
	      i++;
	    }

	}
      //else
      //printf("Label %d -> %d (correspondance %3.1f%%)\n", srcLabel, tabNewLabel[srcLabel], oldOverlap*100);
    }
  
  // tester les correspondances entre les labels reference et source
  // on omet le label 0 (fond)
  for (refLabel=1; refLabel<maxRef; refLabel++)
    {
      oldOverlap = 0;
      for (srcLabel=1; srcLabel<maxSrc; srcLabel++)
	{
	  overlap = (float)tabLabels[srcLabel][refLabel]/(float)tabLabelSizeRef[refLabel];
	  if (overlap > oldOverlap)
	    {
	      //tabNewLabel[srcLabel] = refLabel;
	      oldOverlap = overlap;
	    }
	}
      if(oldOverlap==0)
	{
	  //tabNewLabel[srcLabel] = CurrentNewLabel;
	  //CurrentNewLabel++;
	  printf("Disparition de Lésions : Label %d (image reference) de %d voxels", refLabel, tabLabelSizeRef[refLabel]);
	  VolumeLesionDiff += tabLabelSizeRef[refLabel];
	  NbDisparitionLesions++;

	  //Connaitre la position d'un point de ce label 

	  i=0;
	  j=0;
	  k=0;
	  while(i<w && dect!=1)
	    {
	      j = 0;
	      while(j<h && dect!=1)
		{
		  k=0;
		  while(k<d && dect!=1)
		    {
		      if(imref->mri[i][j][k] == refLabel)
			{
			  printf("  - Position : %d %d %d \n", i, j, k);
			  dect=1;
			}
		      k++;
		    }
		  j++;
		}
	      i++;
	    }

	}
      //else
      //printf("Label %d -> %d (correspondance %3.1f%%)\n", srcLabel, tabNewLabel[srcLabel], oldOverlap*100);
    }

  
  printf("Nombre de lésions apparues %d \n",NbApparitionLesions);
  printf("Nombre de lésions disparues %d \n",NbDisparitionLesions);
  printf("Volume des différences: %d (voxels), %f mm3\n", VolumeLesionDiff, VolumeLesionDiff*d3);


  for (refLabel=1; refLabel<maxRef; refLabel++)
	NbVoxelRef += (float)tabLabelSizeRef[refLabel];
  for (srcLabel=1; srcLabel<maxSrc; srcLabel++)
	NbVoxelSrc += (float)tabLabelSize[srcLabel];

  printf("Volume de lésions dans Image Reference (image 1) : %f voxels, %f mm3 \n", NbVoxelRef,NbVoxelRef*d3);
  printf("Volume de lésions dans Image Source (image 2) : %f voxels, %f mm3 \n", NbVoxelSrc,NbVoxelSrc*d3);
  printf("Différence des volumes (image 1 - image 2): %f voxels, %f mm3\n", (NbVoxelRef - NbVoxelSrc), (NbVoxelRef - NbVoxelSrc)*d3);


  for(i=0;i!=w;i++)
    for(j=0;j!=h;j++)
      for(k=0;k!=d;k++)
	if( ((imref->mri[i][j][k]==0) && (imsrc->mri[i][j][k]!=0)) || ((imref->mri[i][j][k]!=0) && (imsrc->mri[i][j][k]==0)) )
	  DiffVoxel += 1;
  printf("Différence à voxels: %d voxels, %f mm3\n", DiffVoxel, DiffVoxel*d3);

  return(1);

  }
/****************************************************************************************
 **  growing_3d_p(i,j,k,im1,imres,fct_grow_3d,tableau,ntab,valeur,typevois)
 */
/*!	Recherche des points connexe verifiant la fonction fct_grow.
**      Cette fonction recherche en utilisant un tableau les points
**      pour lesquel la fonction fct_grow_3d est vraie (1) .
**      La recherche de connexivite est realise sur une croix (6 points)
**
**      Les points connexes sont mis a la valeur passee en dernier parametre
**	les autres sont mis a zero.
**
**      \param x : coord x du point de depart de la croissance
**		\param y : coord y du point de depart de la croissance
**		\param z : coord z du point de depart de la croissance
**      \param im1 : image a traiter
**      \param imres : image resultat (E/S)
**      \param fct_grow_3d : fonction de test
**      \param tableau : tableau de valeur
**		\param ntab : nbre de valeur du tableau
**      \param valeur: long = valeur a donner a la zone connexe dans imres
**		\param typevois: int =defini le type de voisinage utilise (0:croix 1:cube)
**
**     \remark
**                 Il faut que imres et im1 soit different
**  ATTENTION: La gestion des parametres de imres est a la
**             charge de la procedure appelante
**             maxpixel,icomp,rcoeff
*****************************************************************************************/
int	growing_3d_p(int x, int y, int z, grphic3d *im1, grphic3d *imres, int (*fct_grow_3d) (), float *tableau, int ntab, long int valeur, int typevois)
{
  int i,ii,jj,kk,l,n;
  int hght,wdth,dpth;
  int nbpts,nbpts_sor,nbvois=0;
  XPoints_3d *tab,*tab_sor;
  XPoints_3d tabvois[26];   /*Tableau definissant la position des voisins dans le cas de la croix 3d*/
  double taille_actuelle = 2048;


  /*Definition du tableau tabvois dans le cas de la croix 3d*/
  if (typevois==0)
    {
      nbvois=6;
      tabvois[0].x=-1;tabvois[0].y=0;tabvois[0].z=0;
      tabvois[1].x=1;tabvois[1].y=0;tabvois[1].z=0;
      tabvois[2].x=0;tabvois[2].y=-1;tabvois[2].z=0;
      tabvois[3].x=0;tabvois[3].y=1;tabvois[3].z=0;
      tabvois[4].x=0;tabvois[4].y=0;tabvois[4].z=-1;
      tabvois[5].x=0;tabvois[5].y=0;tabvois[5].z=1;
    }

  /*Definition du tableau tabvois dans le cas du cube 3d*/
  if (typevois==1)
    {
      nbvois=26;
      tabvois[0].x=-1;tabvois[0].y=-1;tabvois[0].z=1;
      tabvois[1].x=0;tabvois[1].y=-1;tabvois[1].z=1;
      tabvois[2].x=1;tabvois[2].y=-1;tabvois[2].z=1;
      tabvois[3].x=-1;tabvois[3].y=0;tabvois[3].z=1;
      tabvois[4].x=0;tabvois[4].y=0;tabvois[4].z=1;
      tabvois[5].x=1;tabvois[5].y=0;tabvois[5].z=1;
      tabvois[6].x=-1;tabvois[6].y=1;tabvois[6].z=1;
      tabvois[7].x=0;tabvois[7].y=1;tabvois[7].z=1;
      tabvois[8].x=1;tabvois[8].y=1;tabvois[8].z=1;
      tabvois[9].x=-1;tabvois[9].y=-1;tabvois[9].z=0;
      tabvois[10].x=0;tabvois[10].y=-1;tabvois[10].z=0;
      tabvois[11].x=1;tabvois[11].y=-1;tabvois[11].z=0;
      tabvois[12].x=-1;tabvois[12].y=0;tabvois[12].z=0;
      tabvois[13].x=1;tabvois[13].y=0;tabvois[13].z=0;
      tabvois[14].x=-1;tabvois[14].y=1;tabvois[14].z=0;
      tabvois[15].x=0;tabvois[15].y=1;tabvois[15].z=0;
      tabvois[16].x=1;tabvois[16].y=1;tabvois[16].z=0;
      tabvois[17].x=-1;tabvois[17].y=-1;tabvois[17].z=-1;
      tabvois[18].x=0;tabvois[18].y=-1;tabvois[18].z=-1;
      tabvois[19].x=1;tabvois[19].y=-1;tabvois[19].z=-1;
      tabvois[20].x=-1;tabvois[20].y=0;tabvois[20].z=-1;
      tabvois[21].x=0;tabvois[21].y=0;tabvois[21].z=-1;
      tabvois[22].x=1;tabvois[22].y=0;tabvois[22].z=-1;
      tabvois[23].x=-1;tabvois[23].y=1;tabvois[23].z=-1;
      tabvois[24].x=0;tabvois[24].y=1;tabvois[24].z=-1;
      tabvois[25].x=1;tabvois[25].y=1;tabvois[25].z=-1;
    }

  tab=(XPoints_3d *)malloc((size_t) 2048*sizeof(XPoints_3d));
  tab_sor=(XPoints_3d *)malloc((size_t) 2048*sizeof(XPoints_3d));

  wdth=im1->width;
  hght=im1->height;
  dpth=im1->depth;
  nbpts=1;
  tab[0].x=x;
  tab[0].y=y;
  tab[0].z=z;
  imres->mri[x][y][z]=(TYPEMRI3D)valeur;
  n=1;



  while (nbpts!=0)
    {
      //Anciennement : 
      /*
      if(nbpts > n*1024)
	{
	  n=n+1;  
	  tab=(XPoints_3d *)
	    realloc((XPoints_3d *) tab,(size_t) n*2048*sizeof(XPoints_3d));
	  tab_sor=(XPoints_3d *)
	    realloc((XPoints_3d *) tab_sor,(size_t) n*2048*sizeof(XPoints_3d)); 

	}
      */

      nbpts_sor=0;
      //Comptage du nb de points
      for (i=0;i<nbpts;i++)
	for (l=0;l<nbvois;l++)
	  {
	    ii=tab[i].x+tabvois[l].x;
	    jj=tab[i].y+tabvois[l].y;
	    kk=tab[i].z+tabvois[l].z;
	    if ((ii>=0) && (ii<wdth) && (jj>=0) && (jj<hght) && (kk>=0) && (kk<=dpth))
	      {
		if (imres->mri[ii][jj][kk]==0 && fct_grow_3d(tab[i].x,tab[i].y,tab[i].z
							     ,ii,jj,kk,im1,tableau,ntab))
		  {
		    nbpts_sor=nbpts_sor+1;
		  }
	      }
	  }
       //Reallocation
      if(nbpts_sor > taille_actuelle)
	{
	  tab_sor=(XPoints_3d *)
	    realloc((XPoints_3d *) tab_sor,(size_t) nbpts_sor*sizeof(XPoints_3d));

	  if(tab_sor==NULL)
	    printf("PB d'allocation de tab_sor ! %d %d %f \n", nbpts_sor,sizeof(XPoints_3d),taille_actuelle);	  
	}


      //----------------

      nbpts_sor=0;
      for (i=0;i<nbpts;i++)
	for (l=0;l<nbvois;l++)
	  {
	    ii=tab[i].x+tabvois[l].x;
	    jj=tab[i].y+tabvois[l].y;
	    kk=tab[i].z+tabvois[l].z;
	    if ((ii>=0) && (ii<wdth) && (jj>=0) && (jj<hght) && (kk>=0) && (kk<=dpth))
	      {
		if (imres->mri[ii][jj][kk]==0 && fct_grow_3d(tab[i].x,tab[i].y,tab[i].z
							     ,ii,jj,kk,im1,tableau,ntab))
		  {
		    imres->mri[ii][jj][kk]=(TYPEMRI3D)valeur;
		    tab_sor[nbpts_sor].x=ii;
		    tab_sor[nbpts_sor].y=jj;
		    tab_sor[nbpts_sor].z=kk;
		    nbpts_sor=nbpts_sor+1;
		  }
	      }
	  }

      //Reallocation
      if(nbpts_sor > taille_actuelle)
	{
	  tab=(XPoints_3d *)
	    realloc((XPoints_3d *) tab,(size_t) nbpts_sor*sizeof(XPoints_3d));
	 
	  if(tab==NULL)
	    printf("PB d'allocation de tab ! %d %d %f \n", nbpts_sor,sizeof(XPoints_3d),taille_actuelle);	  

	  taille_actuelle = nbpts_sor;
	}


      for (i=0;i<nbpts_sor;i++)
	{
	  tab[i].x=tab_sor[i].x;
	  tab[i].y=tab_sor[i].y;
	  tab[i].z=tab_sor[i].z;
	}

      nbpts=nbpts_sor;
      
      if(nbpts_sor > taille_actuelle)
	printf("PB dans growing_3d_p | fichier trai_3d.c: tableau de taille insuffisante ! %d %f\n", nbpts_sor, taille_actuelle);

    }

  free(tab);
  free(tab_sor);
  return(1);
}


/***************************************************************************************
 ** -- growing_cerv_fonction_3d_p(x,y,z,im1,imres,pourcent) ----------------------------
 **	Recherche de la composante connexe associe au cerveau a partir du point
 **      de depart (x,y,z)
 **
 **      im1 : image a traiter
 **      imres : image resultat
 **      pourcent : critere (de 0 a 1) de difference entre deux voisins en pourcentage de ndg.
 **      Remarque :
 **                 Il faut que imres et im1 soit different
 **  ATTENTION: La gestion des parametres de imres est a la charge de la procedure appelante
 **             maxpixel,icomp,rcoeff
 *****************************************************************************************/
int	growing_cerv_fonction_3d_p(int x, int y, int z, grphic3d *im1, grphic3d *imres, float pourcent)
{
  int i,ii,jj,kk,l,n;
  int hght,wdth,dpth;
  int nbpts,nbpts_sor;
  XPoints_3d *tab,*tab_sor;
  XPoints_3d tabvois[26];   /*Tableau definissant la position des voisins dans le cas de la croix 3d*/
  float  pix1,pix2;
  long   volmaxpix,volpix;

  volmaxpix=(long)((float)2000000/im1->dx/im1->dy/im1->dz);
  volpix=0;


  /*Definition du tableau tabvois dans le cas de la croix 3d*/
  tabvois[0].x=-1;tabvois[0].y=-1;tabvois[0].z=1;
  tabvois[1].x=0;tabvois[1].y=-1;tabvois[1].z=1;
  tabvois[2].x=1;tabvois[2].y=-1;tabvois[2].z=1;
  tabvois[3].x=-1;tabvois[3].y=0;tabvois[3].z=1;
  tabvois[4].x=0;tabvois[4].y=0;tabvois[4].z=1;
  tabvois[5].x=1;tabvois[5].y=0;tabvois[5].z=1;
  tabvois[6].x=-1;tabvois[6].y=1;tabvois[6].z=1;
  tabvois[7].x=0;tabvois[7].y=1;tabvois[7].z=1;
  tabvois[8].x=1;tabvois[8].y=1;tabvois[8].z=1;
  tabvois[9].x=-1;tabvois[9].y=-1;tabvois[9].z=0;
  tabvois[10].x=0;tabvois[10].y=-1;tabvois[10].z=0;
  tabvois[11].x=1;tabvois[11].y=-1;tabvois[11].z=0;
  tabvois[12].x=-1;tabvois[12].y=0;tabvois[12].z=0;
  tabvois[13].x=1;tabvois[13].y=0;tabvois[13].z=0;
  tabvois[14].x=-1;tabvois[14].y=1;tabvois[14].z=0;
  tabvois[15].x=0;tabvois[15].y=1;tabvois[15].z=0;
  tabvois[16].x=1;tabvois[16].y=1;tabvois[16].z=0;
  tabvois[17].x=-1;tabvois[17].y=-1;tabvois[17].z=-1;
  tabvois[18].x=0;tabvois[18].y=-1;tabvois[18].z=-1;
  tabvois[19].x=1;tabvois[19].y=-1;tabvois[19].z=-1;
  tabvois[20].x=-1;tabvois[20].y=0;tabvois[20].z=-1;
  tabvois[21].x=0;tabvois[21].y=0;tabvois[21].z=-1;
  tabvois[22].x=1;tabvois[22].y=0;tabvois[22].z=-1;
  tabvois[23].x=-1;tabvois[23].y=1;tabvois[23].z=-1;
  tabvois[24].x=0;tabvois[24].y=1;tabvois[24].z=-1;
  tabvois[25].x=1;tabvois[25].y=1;tabvois[25].z=-1;


  tab=(XPoints_3d *)malloc((size_t) 2048*sizeof(XPoints_3d));
  tab_sor=(XPoints_3d *)malloc((size_t) 2048*sizeof(XPoints_3d));

  wdth=im1->width;
  hght=im1->height;
  dpth=im1->depth;
  nbpts=1;
  tab[0].x=x;
  tab[0].y=y;
  tab[0].z=z;
  imres->mri[x][y][z]=1;
  n=1;
  while (nbpts!=0 && volpix<volmaxpix)
    {
      if(nbpts > n*1024)
	{
	  n=n+1;
	  tab=(XPoints_3d *)
	    realloc((XPoints_3d *) tab,(size_t) n*2048*sizeof(XPoints_3d));
	  tab_sor=(XPoints_3d *)
	    realloc((XPoints_3d *) tab_sor,(size_t) n*2048*sizeof(XPoints_3d));
	}
      nbpts_sor=0;
      for (i=0;i<nbpts;i++)
	for (l=0;l<26;l++)
	  {


	    ii=tab[i].x+tabvois[l].x;
	    jj=tab[i].y+tabvois[l].y;
	    kk=tab[i].z+tabvois[l].z;
	    if ((ii>=0) && (ii<wdth) && (jj>=0) && (jj<hght) && (kk>=0) && (kk<=dpth))
	      {
		pix1=(im1->rcoeff)*(im1->mri[tab[i].x][tab[i].y][tab[i].z]);
		pix2=(im1->rcoeff)*(im1->mri[ii][jj][kk]);
		if (imres->mri[ii][jj][kk]==0 && (fabs((double)(pix1-pix2))<(double)(pourcent*pix1)))
		  {
		    imres->mri[ii][jj][kk]=1;volpix++;
		    tab_sor[nbpts_sor].x=ii;
		    tab_sor[nbpts_sor].y=jj;
		    tab_sor[nbpts_sor].z=kk;
		    nbpts_sor=nbpts_sor+1;
		  }
	      }


	  }

      for (i=0;i<nbpts_sor;i++)
	{
	  tab[i].x=tab_sor[i].x;
	  tab[i].y=tab_sor[i].y;
	  tab[i].z=tab_sor[i].z;
	}
      nbpts=nbpts_sor;
    }

  free(tab);
  free(tab_sor);
  return(1);
}


/*******************************************
 ** --  imx_extract_cerv_3d_p() ----------------------------
 **
 **    Fonction d'extraction du cerveau en 3d
 **    A partir de otsu
 ********************************************/
int     imx_extract_cerv_3d_p(grphic3d *im1, grphic3d *imres)
{
  grphic3d *imtemp;
  long *vol,volmax,ivolmax,i;


  /* creation de l'image temporaire imtemp */
  imtemp=cr_grphic3d(im1);
  imx_copie_param_3d_p(im1,imtemp);

  imx_otsu_brain_3d_p(im1,imtemp);
  imx_labelcroix_3d_p(imtemp,imtemp);

  /*Calcul des volumes*/
  vol=fct_volume_cell_3d_p(imtemp);

  /*Calcul du volume max et de la cellule correspondante*/
  volmax=0;ivolmax=0;
  for (i=1;i<=vol[0];i++)
    {
      if (vol[i]>volmax)
	{volmax=vol[i];ivolmax=i;}
    }

  /*Selection de la cellule ayant le plus grand volume*/
  imx_select_cell_value_3d_p(imtemp,imtemp,ivolmax);

  /*remplissage*/
  imx_hole_fill_3d_p(imtemp,imtemp);

  /*et logique*/
  imx_and_3d_p(im1,imtemp,imtemp);

  imx_copie_3d_p(imtemp,imres);

  /* liberation de la memoire utilisee par les images temporaires*/
  free_grphic3d (imtemp);
  free (vol);

  imx_inimaxminpixel_3d_p(imres);
  return(1);
}

/*******************************************
 ** --  imx_growing_cerv_3d_p() ----------------------------
 **
 **    Fonction d'extraction du cerveau en 3d
 **    A partir de la croissance de region
 ********************************************/
int     imx_growing_cerv_3d_p(grphic3d *im1, grphic3d *imres, int x0, int y0, int z0)
{
  grphic3d *imtemp;


  /* creation de l'image temporaire imtemp */
  imtemp=cr_grphic3d(im1);
  imx_copie_param_3d_p(im1,imtemp);


  /*growing_cerv_fonction_3d_p(x0,y0,z0,im1,imtemp,0.07);*/
  imx_seg_cerv_grow_dis_3d_p(im1,imtemp,x0,y0,z0);

  /*et logique*/
  imx_and_3d_p(im1,imtemp,imtemp);

  imx_copie_3d_p(imtemp,imres);

  /* liberation de la memoire utilisee par les images temporaires*/
  free_grphic3d (imtemp);

  imx_inimaxminpixel_3d_p(imres);

  return(1);
}



/*******************************************
 ** -- fct_grow_hole_fill_3d(i,j,k,im1,tab,n) ----------------------------------
 **	Fonction pour tester la connexite pour la procedure hole_fill_3d
 **      i0,j0,k0 : point dont on cherche les voisins (pas pt de depart)
 **      i,j,k : pixel a controler
 **      im1 : pointeur de l'image
 **      tab : Tableau contenant le parametres a passer
 **      n : nombre de valeurs dans tab
 *********************************************/
int	fct_grow_hole_fill_3d(int i0, int j0, int k0, int i, int j, int k, grphic3d *im1, float *tab, int n)
{


  if((im1->mri[i][j][k])==0) return(1);
  return(0);

}

/*******************************************
 ** -- fct_grow_con_pos_3d(i0,j0,k0,i,j,k,im1,tab,n) -------------------------
 **	Fonction pour tester la connexite pour la procedure sel_con_3d
 **      i0,j0,k0 : point dont on cherche les voisins (pas pt de depart)
 **      i,j,k : pixel a controler
 **      im1 : pointeur de l'image
 **      tab : Tableau contenant le parametres a passer
 **      n : nombre de valeurs dans tab
 *********************************************/
int	fct_grow_con_pos_3d(int i0, int j0, int k0, int i, int j, int k, grphic3d *im1, float *tab, int n)
{


  if((im1->mri[i][j][k])>tab[0]) return(1);
  return(0);

}

/*******************************************
 ** --  imx_hole_fill_3d_p() ----------------------------
 */
/*!
**    Fonction remplissage de trous en 3d (d'une image binaire)
**    A partir de la croissance de region
**	 \param im1 : image source
**	 \param imres : image resultat (E/S)
********************************************/
void    imx_hole_fill_3d_p(grphic3d *im1, grphic3d *imres)
{
  grphic3d *imtemp;


  /* creation de l'image temporaire imtemp */
  imtemp=cr_grphic3d(im1);
  imx_copie_param_3d_p(im1,imtemp);



  growing_3d_p(0,0,0,im1,imtemp,fct_grow_hole_fill_3d,NULL,0,1,0);


  imx_inv_3d_p(imtemp,imtemp);


  imx_copie_3d_p(imtemp,imres);

  /* liberation de la memoire utilisee par les images temporaires*/
  free_grphic3d (imtemp);

  imres->icomp=0;
  imres->rcoeff=1;
  imres->max_pixel=1;
  imres->min_pixel=0;
  imres->cutoff_max=1;
  imres->cutoff_min=0;
}


/*******************************************
 ** --  imx_hole_fill_3d_dir_p() ----------------------------
 **
 **    Fonction remplissage de trous en 3d (d'une image binaire)
 **    A partir de la croissance de region
 ********************************************/
void    imx_hole_fill_3d_dir_p(grphic3d *imdeb, grphic3d *imres)
{
  int i,j,k,l;
  int trou;
  int w,h,d;

  imx_copie_param_3d_p(imdeb,imres);
  imx_copie_3d_p(imdeb,imres);

  w = imdeb->width;
  h = imdeb->height;
  d = imdeb->depth;

  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      for (k=0;k<d;k++)
        if (imdeb->mri[i][j][k]==0)
          {
	    trou=1;
	    l=i;
	    while (trou==1 && l<w && imdeb->mri[l][j][k]==0) {l++;}
	    if (l==w) trou=0;
	    l=i;
	    while (trou==1 && l>=0 && imdeb->mri[l][j][k]==0) {l--;}
	    if (l==-1) trou=0;

	    l=j;
	    while (trou==1 && l<h && imdeb->mri[i][l][k]==0) {l++;}
	    if (l==h) trou=0;
	    l=j;
	    while (trou==1 && l>=0 && imdeb->mri[i][l][k]==0) {l--;}
	    if (l==-1) trou=0;

	    l=k;
	    while (trou==1 && l<d && imdeb->mri[i][j][l]==0) {l++;}
	    if (l==d) trou=0;
	    l=k;
	    while (trou==1 && l>=0 && imdeb->mri[i][j][l]==0) {l--;}
	    if (l==-1) trou=0;

	    if (trou==1) imres->mri[i][j][k]=(TYPEMRI3D)imres->max_pixel;

          }
}


/*! \ingroup ToolBoxImage3d @{*/
/*******************************************
 ** --  imx_hole_fill_3d_cpc_p() ----------------------------
 */
/*!
**
**	\brief   Fonction remplissage de trous en 3d (d'une image binaire)
**    coupe par coupe
**   \param imdeb : image source
**   \param imres : image dest (E/S)
**   \param  type :
- 1 par rapport X
- 2 par rapport Y
- 3 par rapport Z
- 4 par rapport a tout
********************************************/
void    imx_hole_fill_3d_cpc_p(grphic3d *imdeb, grphic3d *imres, int type)
{
  grphic3d *imtemp,*imtdeb;
  int wdth,hght,dpth;
  int i,j,k,nb,xp,yp,xv,yv;
  IMXPoint *voisin;
  int maxDim1, maxDim2;

  wdth=imdeb->width;hght=imdeb->height;dpth=imdeb->depth;

  imtemp=cr_grphic3d(imdeb);
  imx_copie_param_3d_p(imdeb,imtemp);
  imx_copie_3d_p(imdeb,imtemp);
  imtdeb=cr_grphic3d(imdeb);
  imx_copie_param_3d_p(imdeb,imtdeb);
  imx_copie_3d_p(imdeb,imtdeb);

  maxDim1=MAXI(wdth,hght); maxDim2=MAXI(dpth,MINI(wdth,hght));
  voisin=CALLOC(maxDim1*maxDim2,IMXPoint);
  
  /*remplissage suivant x*/
  if (type==1 || type ==4)
    {
      for (i=0;i<wdth;i++)
	{
	  /*on empile les 4 coins */
	  nb=0;
	  if (imtemp->mri[i][0][0]==0)
	    {
	      voisin[nb].x=0;voisin[nb].y=0;
	      imtemp->mri[i][0][0]=1;
	      nb++;
	    }
	  if (imtemp->mri[i][hght-1][0]==0)
	    {
	      voisin[nb].x=hght-1;voisin[nb].y=0;
	      imtemp->mri[i][hght-1][0]=1;
	      nb++;
	    }
	  if (imtemp->mri[i][0][dpth-1]==0)
	    {
	      voisin[nb].x=0;voisin[nb].y=dpth-1;
	      imtemp->mri[i][0][dpth-1]=1;
	      nb++;
	    }
	  if (imtemp->mri[i][hght-1][dpth-1]==0)
	    {
	      voisin[nb].x=hght-1;voisin[nb].y=dpth-1;
	      imtemp->mri[i][hght-1][dpth-1]=1;
	      nb++;
	    }
	  if (nb==0) {
	    PUT_ERROR("Error in imx_hole_fill : No starting point (i)");
	    free_grphic3d(imtemp);free_grphic3d(imtdeb);
	    free(voisin);
	    return ;;
	  }

	  while (nb>=0)
	    {
	      /*on depile un point*/
	      xp=voisin[nb].x;yp=voisin[nb].y;nb--;

	      /*les voisins a zero sont mis a 1 et empiles*/
	      xv=xp-1;yv=yp;
	      if (xv>=0 && xv<hght && yv>=0 && yv<dpth)
		if (imtemp->mri[i][xv][yv]==0)
		  {nb++;imtemp->mri[i][xv][yv]=1;voisin[nb].x=xv;voisin[nb].y=yv;}
	      xv=xp;yv=yp-1;
	      if (xv>=0 && xv<hght && yv>=0 && yv<dpth)
		if (imtemp->mri[i][xv][yv]==0)
		  {nb++;imtemp->mri[i][xv][yv]=1;voisin[nb].x=xv;voisin[nb].y=yv;}
	      xv=xp+1;yv=yp;
	      if (xv>=0 && xv<hght && yv>=0 && yv<dpth)
		if (imtemp->mri[i][xv][yv]==0)
		  {nb++;imtemp->mri[i][xv][yv]=1;voisin[nb].x=xv;voisin[nb].y=yv;}
	      xv=xp;yv=yp+1;
	      if (xv>=0 && xv<hght && yv>=0 && yv<dpth)
		if (imtemp->mri[i][xv][yv]==0)
		  {nb++;imtemp->mri[i][xv][yv]=1;voisin[nb].x=xv;voisin[nb].y=yv;}
	    }
	}

      /*on inverse imtemp et on l'ajoute a imtdeb*/
      for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	  for (k=0;k<dpth;k++)
	    if (imtdeb->mri[i][j][k]>0 || imtemp->mri[i][j][k]==0) {imtemp->mri[i][j][k]=1;imtdeb->mri[i][j][k]=1;}
	    else {imtemp->mri[i][j][k]=0;imtdeb->mri[i][j][k]=0;}
    }

  /*remplissage suivant y*/
  if (type==2 || type==4)
    {
      for (j=0;j<hght;j++)
	{
	  /*on empile les 4 coins */
	  nb=0;
	  if (imtemp->mri[0][j][0]==0)
	    {
	      voisin[nb].x=0;voisin[nb].y=0;
	      imtemp->mri[0][j][0]=1;
	      nb++;
	    }
	  if (imtemp->mri[wdth-1][j][0]==0)
	    {
	      voisin[nb].x=wdth-1;voisin[nb].y=0;
	      imtemp->mri[wdth-1][j][0]=1;
	      nb++;
	    }
	  if (imtemp->mri[0][j][dpth-1]==0)
	    {
	      voisin[nb].x=0;voisin[nb].y=dpth-1;
	      imtemp->mri[0][j][dpth-1]=1;
	      nb++;
	    }
	  if (imtemp->mri[wdth-1][j][dpth-1]==0)
	    {
	      voisin[nb].x=wdth-1;voisin[nb].y=dpth-1;
	      imtemp->mri[wdth-1][j][dpth-1]=1;
	      nb++;
	    }
	  if (nb==0) {
	    PUT_ERROR("Error in imx_hole_fill : No starting point (j) ");
	    free_grphic3d(imtemp);free_grphic3d(imtdeb);
	    free(voisin);
	    return ;
	  }

	  while (nb>=0)
	    {
	      /*on depile un point*/
	      xp=voisin[nb].x;yp=voisin[nb].y;nb--;

	      /*les voisins a zero sont mis a 1 et empiles*/
	      xv=xp-1;yv=yp;
	      if (xv>=0 && xv<wdth && yv>=0 && yv<dpth)
		if (imtemp->mri[xv][j][yv]==0)
		  {nb++;imtemp->mri[xv][j][yv]=1;voisin[nb].x=xv;voisin[nb].y=yv;}
	      xv=xp;yv=yp-1;
	      if (xv>=0 && xv<wdth && yv>=0 && yv<dpth)
		if (imtemp->mri[xv][j][yv]==0)
		  {nb++;imtemp->mri[xv][j][yv]=1;voisin[nb].x=xv;voisin[nb].y=yv;}
	      xv=xp+1;yv=yp;
	      if (xv>=0 && xv<wdth && yv>=0 && yv<dpth)
		if (imtemp->mri[xv][j][yv]==0)
		  {nb++;imtemp->mri[xv][j][yv]=1;voisin[nb].x=xv;voisin[nb].y=yv;}
	      xv=xp;yv=yp+1;
	      if (xv>=0 && xv<wdth && yv>=0 && yv<dpth)
		if (imtemp->mri[xv][j][yv]==0)
		  {nb++;imtemp->mri[xv][j][yv]=1;voisin[nb].x=xv;voisin[nb].y=yv;}
	    }
	}

      /*on inverse imtemp et on l'ajoute a imtdeb*/
      for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	  for (k=0;k<dpth;k++)
	    if (imtdeb->mri[i][j][k]>0 || imtemp->mri[i][j][k]==0) {imtemp->mri[i][j][k]=1;imtdeb->mri[i][j][k]=1;}
	    else {imtemp->mri[i][j][k]=0;imtdeb->mri[i][j][k]=0;}
    }

  /*remplissage suivant z*/
  if (type==3 || type==4)
    {
      for (k=0;k<dpth;k++)
	{
	  /*on empile les 4 coins */
	  nb=0;
	  if (imtemp->mri[0][0][k]==0)
	    {
	      voisin[nb].x=0;voisin[nb].y=0;
	      imtemp->mri[0][0][k]=1;
	      nb++;
	    }
	  if (imtemp->mri[wdth-1][0][k]==0)
	    {
	      voisin[nb].x=wdth-1;voisin[nb].y=0;
	      imtemp->mri[wdth-1][0][k]=1;
	      nb++;
	    }
	  if (imtemp->mri[0][hght-1][k]==0)
	    {
	      voisin[nb].x=0;voisin[nb].y=hght-1;
	      imtemp->mri[0][hght-1][k]=1;
	      nb++;
	    }
	  if (imtemp->mri[wdth-1][hght-1][k]==0)
	    {
	      voisin[nb].x=wdth-1;voisin[nb].y=hght-1;
	      imtemp->mri[wdth-1][hght-1][k]=1;
	      nb++;
	    }
	  if (nb==0) {
	    PUT_ERROR("Error in imx_hole_fill : No starting point (k)");
	    free_grphic3d(imtemp);free_grphic3d(imtdeb);
	    free(voisin);
	    return ;
	  }

	  while (nb>=0)
	    {
	      /*on depile un point*/
	      xp=voisin[nb].x;yp=voisin[nb].y;nb--;

	      /*les voisins a zero sont mis a 1 et empiles*/
	      xv=xp-1;yv=yp;
	      if (xv>=0 && xv<wdth && yv>=0 && yv<hght)
		if (imtemp->mri[xv][yv][k]==0)
		  {nb++;imtemp->mri[xv][yv][k]=1;voisin[nb].x=xv;voisin[nb].y=yv;}
	      xv=xp;yv=yp-1;
	      if (xv>=0 && xv<wdth && yv>=0 && yv<hght)
		if (imtemp->mri[xv][yv][k]==0)
		  {nb++;imtemp->mri[xv][yv][k]=1;voisin[nb].x=xv;voisin[nb].y=yv;}
	      xv=xp+1;yv=yp;
	      if (xv>=0 && xv<wdth && yv>=0 && yv<hght)
		if (imtemp->mri[xv][yv][k]==0)
		  {nb++;imtemp->mri[xv][yv][k]=1;voisin[nb].x=xv;voisin[nb].y=yv;}
	      xv=xp;yv=yp+1;
	      if (xv>=0 && xv<wdth && yv>=0 && yv<hght)
		if (imtemp->mri[xv][yv][k]==0)
		  {nb++;imtemp->mri[xv][yv][k]=1;voisin[nb].x=xv;voisin[nb].y=yv;}
	    }
	}

      /*on inverse imtemp et on l'ajoute a imtdeb*/
      for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	  for (k=0;k<dpth;k++)
	    if (imtdeb->mri[i][j][k]>0 || imtemp->mri[i][j][k]==0) {imtemp->mri[i][j][k]=1;imtdeb->mri[i][j][k]=1;}
	    else {imtemp->mri[i][j][k]=0;imtdeb->mri[i][j][k]=0;}
    }

  imtemp->min_pixel=0;
  imtemp->max_pixel=1;
  imtemp->cutoff_min=0;
  imtemp->cutoff_max=1;
  imtemp->rcoeff=1.0;
  imtemp->icomp=0.0;
  imx_copie_param_3d_p(imtemp,imres);
  imx_copie_3d_p(imtemp,imres);
  free_grphic3d(imtemp);free_grphic3d(imtdeb);
  free(voisin);

}
/*! @} */


/*********************************************************
 ** --  imx_canny_deriche_3d(im_deb,im_res,alpha,window,cc_length) --
 **
 **
 **
 **********************************************************/
void imx_canny_deriche_3d(int im_deb, int im_res, double alpha, int window, int cc_length)
{
  grphic3d *imdeb,*imres;

  imdeb=ptr_img_3d(im_deb);
  imres=ptr_img_3d(im_res);

  imx_canny_deriche_3d_p(imdeb,imres,alpha,window,cc_length);

}

/*********************************************************
 ** --  imx_canny_deriche_3d_p(imdeb,imres,alpha,window) --
 **
 **
 **
 **********************************************************/
void imx_canny_deriche_3d_p(grphic3d *imdeb, grphic3d *imres, double alpha, int window, int cc_length)
{
  grphic3d *gx,*gy,*gz;

  gx=cr_grphic3d(imdeb);
  gy=cr_grphic3d(imdeb);
  gz=cr_grphic3d(imdeb);


  initialize_gradient_3d(gx,gy,gz);

  compute_gx_3d(imdeb,gx,alpha,window);   
  compute_gy_3d(imdeb,gy,alpha,window);   
  compute_gz_3d(imdeb,gz,alpha,window);   

  compute_module_3d(gx,gy,gz,imres);   

  gradient_maxima_3d(gx,gy,gz,imres);  


  imx_hysteresis_thresholding_3d_p(imres,imres,cc_length); 

  FREE(gx);
  FREE(gy);
  FREE(gz);

}

/*********************************************************
 ** ------------ hysteresis_thresholding_3d() ---------------------
 **
 **
 **
 **********************************************************/
void hysteresis_thresholding_3d(void)
{
  int im_deb,im_res;
  int cc_length;
  int e=0;

  im_deb=GET_PLACE3D(TEXT0157);
  im_res=GET_PLACE3D(TEXT0157);
  cc_length=(int) GET_FLOAT("Connected component length ?", 0, &e);


  imx_hysteresis_thresholding_3d(im_deb,im_res,cc_length);

  imx_copie_param_3d(im_deb,im_res);
  imx_iniparaimg_3d(im_res);
  show_picture_3d(im_res);

}

/*********************************************************
 ** ------------ imx_hysteresis_thresholding_3d(im_deb,im_res,cc_length) ---------------------
 **
 **
 **
 **********************************************************/
void imx_hysteresis_thresholding_3d(int im_deb, int im_res, int cc_length)
{

  grphic3d *imdeb,*imres;

  imdeb=ptr_img_3d(im_deb);
  imres=ptr_img_3d(im_res);


  imx_hysteresis_thresholding_3d_p(imdeb,imres,cc_length);

}

/*********************************************************
 ** ------------ imx_hysteresis_thresholding_3d_p(imdeb,imres,cc_length) ---------------------
 **
 **
 **
 **********************************************************/
void imx_hysteresis_thresholding_3d_p(grphic3d *imdeb, grphic3d *imres, int cc_length)
{
  grphic3d *temp1=NULL,*temp2=NULL,*temp3=NULL;
  int i,j,k;
  int width,height,depth;
  int *hl=NULL;
  int *corres=NULL;
  TYPEMRI3D ***temp1MRI=NULL, ***temp2MRI=NULL, ***temp3MRI=NULL, ***imresMRI=NULL;
  TYPEMRI3D val;

  width=imdeb->width;
  height=imdeb->height;
  depth=imdeb->depth;

  //allocation des images temporaires
  temp1=cr_grphic3d(imdeb); temp2=cr_grphic3d(imdeb); temp3=cr_grphic3d(imdeb);
  if ((!temp1)||(!temp2)||(!temp3))
  { fprintf (stderr, "erreur alocation memoire dans imx_hysteresis_thresholding_3d_p\n"); goto end_func; }
  temp1MRI=temp1->mri; temp2MRI=temp2->mri; temp3MRI=temp3->mri;

  imresMRI=imres->mri;
  
  imx_inimaxminpixel_3d_p(imdeb);  
  imx_otsu_thresholding_3d_p(imdeb,temp1,2);

  //--- Deux seuils (bas-haut) ---//

  for(i=0;i<width;i++)
    for(j=0;j<height;j++)
      for(k=0;k<depth;k++)
    	{
    	  temp2MRI[i][j][k]=0;
    	  if(temp1MRI[i][j][k]==2)
   	    {
   	      temp1MRI[i][j][k]=1;
   	      temp2MRI[i][j][k]=1;
   	    }
    	}

  imx_labelcroix_3d_p(temp1,temp3);
  imx_inimaxminpixel_3d_p(temp3);

  hl=alloc_ivector(temp3->max_pixel+1);
  corres=alloc_ivector(temp3->max_pixel+1);

  for(i=0;i<temp3->max_pixel;i++)
    {
      hl[i]=0;
      corres[i]=0;
    }

  for(i=0;i<width;i++)
    for(j=0;j<height;j++)
      for(k=0;k<depth;k++)
    	{
       val=temp3MRI[i][j][k];
    	  if(val!=0)
   	    {
   	      hl[val]++;
   	      if(temp2MRI[i][j][k]==1) corres[val]=1;
   	    }
    	}

  for(i=0;i<width;i++)
    for(j=0;j<height;j++)
      for(k=0;k<depth;k++)
      {
       val=temp3MRI[i][j][k]; 
	     if(hl[val]>cc_length && corres[val]==1) imresMRI[i][j][k]=1;
	     else imresMRI[i][j][k]=0;
      }

end_func:

  if (hl) free_ivector(hl, temp3->max_pixel+1);
  if (corres) free_ivector(corres, temp3->max_pixel+1);
  if (temp1) free_grphic3d(temp1);
  if (temp2) free_grphic3d(temp2);
  if (temp3) free_grphic3d(temp3);
}

/*********************************************************
 ** --  filter_d(x,alpha,c) --
 **
 **
 **
 **********************************************************/
double filter_d(int x, double alpha, double c)
{
  double dx,res;

  dx=(double)x;

  res=-c*dx*exp(-1.0*alpha*fabs(dx));

  return(res);

}

/*********************************************************
 ** --  filter_l(x,alpha,s) --
 **
 **
 **
 **********************************************************/
double filter_l(int x, double alpha, double s)
{
  double dx,res;

  dx=(double)x;

  res=s*(alpha*fabs(dx)+1.0)*exp(-1.0*alpha*fabs(dx));

  return(res);

}

/*********************************************************
 ** --  compute_module_3d(gx,gy,gz,imres) --
 **
 **
 **
 **********************************************************/
void compute_module_3d(grphic3d *gx, grphic3d *gy, grphic3d *gz, grphic3d *imres)
{
  int i,j,k;
  double x,y,z;
  int w,h,d;
  w = gx->width;
  h = gx->height;
  d = gx->depth;

  for(i=0;i<w;i++)
    for(j=0;j<h;j++)
      for(k=0;k<d;k++)
	{
	  x=(double)gx->mri[i][j][k];
	  y=(double)gy->mri[i][j][k];
	  z=(double)gz->mri[i][j][k];

	  imres->mri[i][j][k]=(int)(0.5+sqrt(x*x+y*y+z*z));
	}

}

/*********************************************************
 ** -- initialize_gradient_3d(gx,gy,gz)  --
 **
 **
 **
 **********************************************************/
void initialize_gradient_3d(grphic3d *gx, grphic3d *gy, grphic3d *gz)
{
  int i,j,k;
  int w,h,d;
  w = gx->width;
  h = gx->height;
  d = gx->depth;

  for(i=0;i<w;i++)
    for(j=0;j<h;j++)
      for(k=0;k<d;k++)
	{
	  gx->mri[i][j][k]=0;
	  gy->mri[i][j][k]=0;
	  gz->mri[i][j][k]=0;
	}

}

/*********************************************************
 ** -- compute_gx_3d(imdeb,gx,alpha,wnd)  --
 **
 **
 **
 **********************************************************/
void compute_gx_3d(grphic3d *imdeb, grphic3d *gx, double alpha, int wnd)
{
  int i,j,k;
  int width,height,depth;
  double c,s;
  double expon;
  double norm_d=0.0,norm_l=0.0;
  double *filtre_d=NULL, *filtre_l=NULL;
  double ***tmp_res_x=NULL;
  double *row1DIn=NULL, *row1DOut=NULL;
  TYPEMRI3D ***imdebMRI=NULL, ***gxMRI=NULL;
   
  width=imdeb->width;
  height=imdeb->height;
  depth=imdeb->depth;

  imdebMRI=imdeb->mri; gxMRI=gx->mri;
  
  filtre_d=CALLOC(wnd, double);
  if (!filtre_d) { fprintf(stderr, "memory allocation failed in calcul_norme_gradient_canny_deriche_3d_p\n"); goto end_func; }
  filtre_l=CALLOC(wnd, double);
  if (!filtre_l) { fprintf(stderr, "memory allocation failed in calcul_norme_gradient_canny_deriche_3d_p\n"); goto end_func; }
  tmp_res_x=alloc_dmatrix_3d(width, height, depth);
  if (!tmp_res_x) { fprintf(stderr, "memory allocation failed in calcul_norme_gradient_canny_deriche_3d_p\n"); goto end_func; }

  expon=exp(-1.0*alpha);
  c=(1.0-expon)*(1.0-expon)/expon;
  s=(1.0-expon)*(1.0-expon)/(1.0+2.0*alpha*expon-exp(-2.0*alpha));
  for(i=-wnd/2;i<=wnd/2;i++)
  {
   filtre_d[i+wnd/2]=filter_d(i,alpha,c);
   filtre_l[i+wnd/2]=filter_l(i,alpha,s);
   norm_d+=fabs(filtre_d[i+wnd/2]);
   norm_l+=fabs(filtre_l[i+wnd/2]);
  }
  for(i=-wnd/2;i<=wnd/2;i++)
  {
   filtre_d[i+wnd/2]/=norm_d;
   filtre_l[i+wnd/2]/=norm_l;
  }

 //filtrage en x avec d(x)
 row1DIn=CALLOC(width, double); row1DOut=CALLOC(width, double);
 for (j=0;j<height;j++)
 {
  for (k=0;k<depth;k++)
  {
   for (i=0;i<width;i++) row1DIn[i]=(double)imdebMRI[i][j][k];
   convoluer_filtre_1D_sans_bord(row1DIn, row1DOut, width, filtre_d, wnd);
   for (i=0;i<width;i++) tmp_res_x[i][j][k]=row1DOut[i];
  }
 }
 FREE(row1DIn); row1DIn=NULL; FREE(row1DOut); row1DOut=NULL;

 //filtrage en y avec l(x)
 row1DIn=CALLOC(height, double); row1DOut=CALLOC(height, double);
 for (i=wnd/2;i<(width-wnd/2);i++)
 {
  for (k=0;k<depth;k++)
  {
   for (j=0;j<height;j++) row1DIn[j]=tmp_res_x[i][j][k];
   convoluer_filtre_1D_sans_bord(row1DIn, row1DOut, height, filtre_l, wnd);
   for (j=0;j<height;j++) tmp_res_x[i][j][k]=row1DOut[j];
  }
 }
 FREE(row1DIn); row1DIn=NULL; FREE(row1DOut); row1DOut=NULL;

 //filtrage en z avec l(x)
 row1DIn=CALLOC(depth, double); row1DOut=CALLOC(depth, double);
 for (i=wnd/2;i<(width-wnd/2);i++)
 {
  for (j=wnd/2;j<(height-wnd/2);j++)
  {
   for (k=0;k<depth;k++) row1DIn[k]=tmp_res_x[i][j][k];
   convoluer_filtre_1D_sans_bord(row1DIn, row1DOut, depth, filtre_l, wnd);
   for (k=0;k<depth;k++) gxMRI[i][j][k]=(TYPEMRI3D)floor(0.5+row1DOut[k]);
  }
 }
 FREE(row1DIn); row1DIn=NULL; FREE(row1DOut); row1DOut=NULL;

end_func:

  if (filtre_d) FREE(filtre_d);
  if (filtre_l) FREE(filtre_l);
  if (tmp_res_x) free_dmatrix_3d(tmp_res_x);
}

/*********************************************************
 ** -- compute_gy_3d(imdeb,gy,alpha,wnd)  --
 **
 **
 **
 **********************************************************/

void compute_gy_3d(grphic3d *imdeb, grphic3d *gy, double alpha, int wnd)
{
  int i,j,k;
  int width,height,depth;
  double c,s;
  double expon;
  double norm_d=0.0,norm_l=0.0;
  double *filtre_d=NULL, *filtre_l=NULL;
  double ***tmp_res_y=NULL;
  double *row1DIn=NULL, *row1DOut=NULL;
  TYPEMRI3D ***imdebMRI=NULL, ***gyMRI=NULL;
  
  width=imdeb->width;
  height=imdeb->height;
  depth=imdeb->depth;

  imdebMRI=imdeb->mri; gyMRI=gy->mri;
  
  filtre_d=CALLOC(wnd, double);
  if (!filtre_d) { fprintf(stderr, "memory allocation failed in calcul_norme_gradient_canny_deriche_3d_p\n"); goto end_func; }
  filtre_l=CALLOC(wnd, double);
  if (!filtre_l) { fprintf(stderr, "memory allocation failed in calcul_norme_gradient_canny_deriche_3d_p\n"); goto end_func; }
  tmp_res_y=alloc_dmatrix_3d(width, height, depth);
  if (!tmp_res_y) { fprintf(stderr, "memory allocation failed in calcul_norme_gradient_canny_deriche_3d_p\n"); goto end_func; }

  expon=exp(-alpha);
  c=(1.0-expon)*(1.0-expon)/expon;
  s=(1.0-expon)*(1.0-expon)/(1.0+2.0*alpha*expon-exp(-2.0*alpha));
  for(i=-wnd/2;i<=wnd/2;i++)
  {
   filtre_d[i+wnd/2]=filter_d(i,alpha,c);
   filtre_l[i+wnd/2]=filter_l(i,alpha,s);
   norm_d+=fabs(filtre_d[i+wnd/2]);
   norm_l+=fabs(filtre_l[i+wnd/2]);
  }
  for(i=-wnd/2;i<=wnd/2;i++)
  {
   filtre_d[i+wnd/2]/=norm_d;
   filtre_l[i+wnd/2]/=norm_l;
  }
  
 //filtrage en y avec d(x)
 row1DIn=CALLOC(height, double); row1DOut=CALLOC(height, double);
 for (i=0;i<width;i++)
 {
  for (k=0;k<depth;k++)
  {
   for (j=0;j<height;j++) row1DIn[j]=(double)imdebMRI[i][j][k];
   convoluer_filtre_1D_sans_bord(row1DIn, row1DOut, height, filtre_d, wnd);
   for (j=0;j<height;j++) tmp_res_y[i][j][k]=row1DOut[j];
  }
 }
 FREE(row1DIn); row1DIn=NULL; FREE(row1DOut); row1DOut=NULL;

 //filtrage en x avec l(x)
 row1DIn=CALLOC(width, double); row1DOut=CALLOC(width, double);
 for (j=wnd/2;j<(height-wnd/2);j++)
 {
  for (k=0;k<depth;k++)
  {
   for (i=0;i<width;i++) row1DIn[i]=tmp_res_y[i][j][k];
   convoluer_filtre_1D_sans_bord(row1DIn, row1DOut, width, filtre_l, wnd);
   for (i=0;i<width;i++) tmp_res_y[i][j][k]=row1DOut[i];
  }
 }
 FREE(row1DIn); row1DIn=NULL; FREE(row1DOut); row1DOut=NULL;

 //filtrage en z avec l(x)
 row1DIn=CALLOC(depth, double); row1DOut=CALLOC(depth, double);
 for (i=wnd/2;i<(width-wnd/2);i++)
 {
  for (j=wnd/2;j<(height-wnd/2);j++)
  {
   for (k=0;k<depth;k++) row1DIn[k]=tmp_res_y[i][j][k];
   convoluer_filtre_1D_sans_bord(row1DIn, row1DOut, depth, filtre_l, wnd);
   for (k=0;k<depth;k++) gyMRI[i][j][k]=(TYPEMRI3D)floor(0.5+row1DOut[k]);
  }
 }
 FREE(row1DIn); row1DIn=NULL; FREE(row1DOut); row1DOut=NULL;

end_func:

  if (filtre_d) FREE(filtre_d);
  if (filtre_l) FREE(filtre_l);
  if (tmp_res_y) free_dmatrix_3d(tmp_res_y);
  
}

/*********************************************************
 ** -- compute_gy_3d(imdeb,gy,alpha,wnd)  --
 **
 **
 **
 **********************************************************/
void compute_gz_3d(grphic3d *imdeb, grphic3d *gz, double alpha, int wnd)
{
  int i,j,k;
  int width,height,depth;
  double c,s;
  double expon;
  double norm_d=0.0,norm_l=0.0;
  double *filtre_d=NULL, *filtre_l=NULL;
  double ***tmp_res_z=NULL;
  double *row1DIn=NULL, *row1DOut=NULL;
  TYPEMRI3D ***imdebMRI=NULL, ***gzMRI=NULL;

  width=imdeb->width;
  height=imdeb->height;
  depth=imdeb->depth;

  imdebMRI=imdeb->mri; gzMRI=gz->mri;

  filtre_d=CALLOC(wnd, double);
  if (!filtre_d) { fprintf(stderr, "memory allocation failed in calcul_norme_gradient_canny_deriche_3d_p\n"); goto end_func; }
  filtre_l=CALLOC(wnd, double);
  if (!filtre_l) { fprintf(stderr, "memory allocation failed in calcul_norme_gradient_canny_deriche_3d_p\n"); goto end_func; }
  tmp_res_z=alloc_dmatrix_3d(width, height, depth);
  if (!tmp_res_z) { fprintf(stderr, "memory allocation failed in calcul_norme_gradient_canny_deriche_3d_p\n"); goto end_func; }

  expon=exp(-alpha);
  c=(1.0-expon)*(1.0-expon)/expon;
  s=(1.0-expon)*(1.0-expon)/(1.0+2.0*alpha*expon-exp(-2.0*alpha));
  for(i=-wnd/2;i<=wnd/2;i++)
  {
   filtre_d[i+wnd/2]=filter_d(i,alpha,c);
   filtre_l[i+wnd/2]=filter_l(i,alpha,s);
   norm_d+=fabs(filtre_d[i+wnd/2]);
   norm_l+=fabs(filtre_l[i+wnd/2]);
  }
  for(i=-wnd/2;i<=wnd/2;i++)
  {
   filtre_d[i+wnd/2]/=norm_d;
   filtre_l[i+wnd/2]/=norm_l;
  }

 //filtrage en z avec d(x)
 row1DIn=CALLOC(depth, double); row1DOut=CALLOC(depth, double);
 for (j=0;j<height;j++)
 {
  for (i=0;i<width;i++)
  {
   for (k=0;k<depth;k++) row1DIn[k]=(double)imdebMRI[i][j][k];
   convoluer_filtre_1D_sans_bord(row1DIn, row1DOut, depth, filtre_d, wnd);
   for (k=0;k<depth;k++) tmp_res_z[i][j][k]=row1DOut[k];
  }
 }
 FREE(row1DIn); row1DIn=NULL; FREE(row1DOut); row1DOut=NULL;

 //filtrage en y avec l(x)
 row1DIn=CALLOC(height, double); row1DOut=CALLOC(height, double);
 for (k=wnd/2;k<(depth-wnd/2);k++)
 {
  for (i=0;i<width;i++)
  {
   for (j=0;j<height;j++) row1DIn[j]=tmp_res_z[i][j][k];
   convoluer_filtre_1D_sans_bord(row1DIn, row1DOut, height, filtre_l, wnd);
   for (j=0;j<height;j++) tmp_res_z[i][j][k]=row1DOut[j];
  }
 }
 FREE(row1DIn); row1DIn=NULL; FREE(row1DOut); row1DOut=NULL;

 //filtrage en x avec l(x)
 row1DIn=CALLOC(width, double); row1DOut=CALLOC(width, double);
 for (k=wnd/2;k<(depth-wnd/2);k++)
 {
  for (j=wnd/2;j<(height-wnd/2);j++)
  {
   for (i=0;i<width;i++) row1DIn[i]=tmp_res_z[i][j][k];
   convoluer_filtre_1D_sans_bord(row1DIn, row1DOut, width, filtre_l, wnd);
   for (i=0;i<width;i++) gzMRI[i][j][k]=(TYPEMRI3D)floor(0.5+row1DOut[i]);
  }
 }
 FREE(row1DIn); row1DIn=NULL; FREE(row1DOut); row1DOut=NULL;

end_func:

  if (filtre_d) FREE(filtre_d);
  if (filtre_l) FREE(filtre_l);
  if (tmp_res_z) free_dmatrix_3d(tmp_res_z);

}

/*********************************************************
 ** -- gradient_maxima_3d(gx,gy,gz,imres)  --
 **
 **
 **
 **********************************************************/
void gradient_maxima_3d(grphic3d *gx, grphic3d *gy, grphic3d *gz, grphic3d *imres)
{
  int i,j,k;
  int width,height,depth;
  double cosa,cosb,cosc;
  double x_up=0,y_up=0,z_up=0,x_down=0,y_down=0,z_down=0;
  double t;
  int val_up,val_down;
  double dx,dy,dz;
  double mod;

  width=gx->width;
  height=gx->height;
  depth=gx->depth;


  for(k=3;k<depth-3;k++)
    for(j=3;j<height-3;j++)
      for(i=3;i<width-3;i++)
	if(imres->mri[i][j][k]!=0)
	  {
	    dx=(double)gx->mri[i][j][k];
	    dy=(double)gy->mri[i][j][k];
	    dz=(double)gz->mri[i][j][k];

           // mod=sqrt(dx*dx+dy*dy+dz*dz);
            mod=NORME_EUCLIDIENNE(dx,dy,dz);

            cosa=dx/mod;
            cosb=dy/mod;
            cosc=dz/mod;

            if(fabs(cosa)>=fabs(cosb) && fabs(cosa)>=fabs(cosc))
              {
		t=1.0/cosa;
		x_up=(double)i+1.0;
		y_up=(double)j+t*cosb;
		z_up=(double)k+t*cosc;
		x_down=(double)i-1.0;
		y_down=(double)j-t*cosb;
		z_down=(double)(k)-t*cosc;
              }

	    else if(fabs(cosb)>=fabs(cosa) && fabs(cosb)>=fabs(cosc))
              {
		t=1.0/cosb;
		x_up=(double)i+t*cosa;
		y_up=(double)j+1.0;
		z_up=(double)k+t*cosc;
		x_down=(double)i-t*cosa;
		y_down=(double)j-1.0;
		z_down=(double)(k)-t*cosc;
              }

	    else if(fabs(cosc)>=fabs(cosa) && fabs(cosc)>=fabs(cosb))
              {
		t=1.0/cosc;
		x_up=(double)i+t*cosa;
		y_up=(double)j+t*cosb;
		z_up=(double)k+1.0;
		x_down=(double)i-t*cosa;
		y_down=(double)j-t*cosb;
		z_down=(double)(k)-1.0;
              }


	    val_up=interpol_gradient_3d(gx,gy,gz,x_up,y_up,z_up);
	    val_down=interpol_gradient_3d(gx,gy,gz,x_down,y_down,z_down);

	    if(imres->mri[i][j][k]<val_up || imres->mri[i][j][k]<=val_down)
	      imres->mri[i][j][k]=0;

	  }



}


/*********************************************************
 ** -- interpol_gradient_3d(gx,gy,gz,x,y,z) --
 **
 **
 **
 **********************************************************/
int interpol_gradient_3d(grphic3d *gx, grphic3d *gy, grphic3d *gz, double x, double y, double z)
{
  double val_gx=0,val_gy=0,val_gz=0;
  int x1,x2,y1,y2,z1,z2;
  double gx1,gx2,gx3,gx4,gy1,gy2,gy3,gy4,gz1,gz2,gz3,gz4;
  double d1,d2,d3,d4;
  double w1,w2,w3,w4,w;
  double d_x1,d_x2,d_y1,d_y2,d_z1,d_z2;
  int result=0;
  TYPEMRI3D ***gxMRI=NULL, ***gyMRI=NULL, ***gzMRI=NULL;

  x1=(int)x;
  y1=(int)y;
  z1=(int)z;
  x2=(int)x+1;
  y2=(int)y+1;
  z2=(int)z+1;

  d_x1=(double)x1;
  d_y1=(double)y1;
  d_z1=(double)z1;
  d_x2=(double)x2;
  d_y2=(double)y2;
  d_z2=(double)z2;

  gxMRI=gx->mri; gyMRI=gy->mri; gzMRI=gz->mri;

  if(z==d_z1)
    {
      d1=NORME_EUCLIDIENNE(x-d_x1,y-d_y2,0.0);
      d2=NORME_EUCLIDIENNE(x-d_x2,y-d_y2,0.0);
      d3=NORME_EUCLIDIENNE(x-d_x1,y-d_y1,0.0);
      d4=NORME_EUCLIDIENNE(x-d_x2,y-d_y1,0.0);

      if(d1!=0.0 && d2!=0.0 && d3!=0.0 && d4!=0.0)
    	{
    	  w1=1.0/d1;
    	  w2=1.0/d2;
    	  w3=1.0/d3;
    	  w4=1.0/d4;
    	  w=w1+w2+w3+w4;

    	  gx1=(double)gxMRI[x1][y2][z1];
    	  gx2=(double)gxMRI[x2][y2][z1];
    	  gx3=(double)gxMRI[x1][y1][z1];
    	  gx4=(double)gxMRI[x2][y1][z1];
    	  gy1=(double)gyMRI[x1][y2][z1];
    	  gy2=(double)gyMRI[x2][y2][z1];
    	  gy3=(double)gyMRI[x1][y1][z1];
    	  gy4=(double)gyMRI[x2][y1][z1];
    	  gz1=(double)gzMRI[x1][y2][z1];
    	  gz2=(double)gzMRI[x2][y2][z1];
    	  gz3=(double)gzMRI[x1][y1][z1];
    	  gz4=(double)gzMRI[x2][y1][z1];

    	  val_gx=(w1*gx1+w3*gx3+w2*gx2+w4*gx4)/w;
    	  val_gy=(w1*gy1+w3*gy3+w2*gy2+w4*gy4)/w;
    	  val_gz=(w1*gz1+w3*gz3+w2*gz2+w4*gz4)/w;

    	}

      else if(d1==0.0)
    	{
    	  val_gx=(double)gxMRI[x1][y2][z1];
    	  val_gy=(double)gyMRI[x1][y2][z1];
    	  val_gz=(double)gzMRI[x1][y2][z1];
    	}

      else if(d2==0.0)
    	{
    	  val_gx=(double)gxMRI[x2][y2][z1];
    	  val_gy=(double)gyMRI[x2][y2][z1];
    	  val_gz=(double)gzMRI[x2][y2][z1];
    	}
    
      else if(d3==0.0)
    	{
    	  val_gx=(double)gxMRI[x1][y1][z1];
    	  val_gy=(double)gyMRI[x1][y1][z1];
    	  val_gz=(double)gzMRI[x1][y1][z1];
    	}

      else if(d4==0.0)
    	{
    	  val_gx=(double)gxMRI[x2][y1][z1];
    	  val_gy=(double)gyMRI[x2][y1][z1];
    	  val_gz=(double)gzMRI[x2][y1][z1];
    	}


      result=(int)floor(0.5+NORME_EUCLIDIENNE(val_gx,val_gy,val_gz));
    }


  else if(x==d_x1)
    {
      d1=NORME_EUCLIDIENNE(0.0,y-d_y2,z-d_z1);
      d2=NORME_EUCLIDIENNE(0.0,y-d_y2,z-d_z2);
      d3=NORME_EUCLIDIENNE(0.0,y-d_y1,z-d_z1);
      d4=NORME_EUCLIDIENNE(0.0,y-d_y1,z-d_z2);

      if(d1!=0.0 && d2!=0.0 && d3!=0.0 && d4!=0.0)
    	{
    	  w1=1.0/d1;
    	  w2=1.0/d2;
    	  w3=1.0/d3;
    	  w4=1.0/d4;
    	  w=w1+w2+w3+w4;

    	  gx1=(double)gxMRI[x1][y2][z1];
    	  gx2=(double)gxMRI[x1][y2][z2];
    	  gx3=(double)gxMRI[x1][y1][z1];
    	  gx4=(double)gxMRI[x1][y1][z2];
    	  gy1=(double)gyMRI[x1][y2][z1];
    	  gy2=(double)gyMRI[x1][y2][z2];
    	  gy3=(double)gyMRI[x1][y1][z1];
    	  gy4=(double)gyMRI[x1][y1][z2];
    	  gz1=(double)gzMRI[x1][y2][z1];
    	  gz2=(double)gzMRI[x1][y2][z2];
    	  gz3=(double)gzMRI[x1][y1][z1];
    	  gz4=(double)gzMRI[x1][y1][z2];

    	  val_gx=(w1*gx1+w3*gx3+w2*gx2+w4*gx4)/w;
    	  val_gy=(w1*gy1+w3*gy3+w2*gy2+w4*gy4)/w;
    	  val_gz=(w1*gz1+w3*gz3+w2*gz2+w4*gz4)/w;

    	}

      else if(d1==0.0)
    	{
    	  val_gx=(double)gxMRI[x1][y2][z1];
    	  val_gy=(double)gyMRI[x1][y2][z1];
    	  val_gz=(double)gzMRI[x1][y2][z1];
    	}

      else if(d2==0.0)
    	{
    	  val_gx=(double)gxMRI[x1][y2][z2];
    	  val_gy=(double)gyMRI[x1][y2][z2];
    	  val_gz=(double)gzMRI[x1][y2][z2];
    	}
    
      else if(d3==0.0)
    	{
    	  val_gx=(double)gxMRI[x1][y1][z1];
    	  val_gy=(double)gyMRI[x1][y1][z1];
    	  val_gz=(double)gzMRI[x1][y1][z1];
    	}

      else if(d4==0.0)
    	{
    	  val_gx=(double)gxMRI[x1][y1][z2];
    	  val_gy=(double)gyMRI[x1][y1][z2];
    	  val_gz=(double)gzMRI[x1][y1][z2];
    	}


      result=(int)floor(0.5+NORME_EUCLIDIENNE(val_gx,val_gy,val_gz));
    }


  else if(y==d_y1)
    {
      d1=NORME_EUCLIDIENNE(x-d_x2,0.0,z-d_z1);
      d2=NORME_EUCLIDIENNE(x-d_x2,0.0,z-d_z2);
      d3=NORME_EUCLIDIENNE(x-d_x1,0.0,z-d_z1);
      d4=NORME_EUCLIDIENNE(x-d_x1,0.0,z-d_z2);


      if(d1!=0.0 && d2!=0.0 && d3!=0.0 && d4!=0.0)
    	{
    	  w1=1.0/d1;
    	  w2=1.0/d2;
    	  w3=1.0/d3;
    	  w4=1.0/d4;
    	  w=w1+w2+w3+w4;

    	  gx1=(double)gxMRI[x2][y1][z1];
    	  gx2=(double)gxMRI[x2][y1][z2];
    	  gx3=(double)gxMRI[x1][y1][z1];
    	  gx4=(double)gxMRI[x1][y1][z2];
    	  gy1=(double)gyMRI[x2][y1][z1];
    	  gy2=(double)gyMRI[x2][y1][z2];
    	  gy3=(double)gyMRI[x1][y1][z1];
    	  gy4=(double)gyMRI[x1][y1][z2];
    	  gz1=(double)gzMRI[x2][y1][z1];
    	  gz2=(double)gzMRI[x2][y1][z2];
    	  gz3=(double)gzMRI[x1][y1][z1];
    	  gz4=(double)gzMRI[x1][y1][z2];

    	  val_gx=(w1*gx1+w3*gx3+w2*gx2+w4*gx4)/w;
    	  val_gy=(w1*gy1+w3*gy3+w2*gy2+w4*gy4)/w;
    	  val_gz=(w1*gz1+w3*gz3+w2*gz2+w4*gz4)/w;

    	}
    
      else if(d1==0.0)
    	{
    	  val_gx=(double)gxMRI[x2][y1][z1];
    	  val_gy=(double)gyMRI[x2][y1][z1];
    	  val_gz=(double)gzMRI[x2][y1][z1];
    	}

      else if(d2==0.0)
    	{
    	  val_gx=(double)gxMRI[x2][y1][z2];
    	  val_gy=(double)gyMRI[x2][y1][z2];
    	  val_gz=(double)gzMRI[x2][y1][z2];
    	}

      else if(d3==0.0)
    	{
    	  val_gx=(double)gxMRI[x1][y1][z1];
    	  val_gy=(double)gyMRI[x1][y1][z1];
    	  val_gz=(double)gzMRI[x1][y1][z1];
    	}

      else if(d4==0.0)
    	{
    	  val_gx=(double)gxMRI[x1][y1][z2];
    	  val_gy=(double)gyMRI[x1][y1][z2];
    	  val_gz=(double)gzMRI[x1][y1][z2];
    	}


      result=(int)floor(0.5+NORME_EUCLIDIENNE(val_gx,val_gy,val_gz));
    }



  return(result);

}

/******************************************************************************
 *********
 *********  ----------------- gaussian_filter_3d()
 *********
 ********************************************************************************/
void gaussian_filter_3d(void)
{
  int im_deb,im_res;
  int i,e=0;
  int val,wndX,wndY,wndZ;
  float sigmaX,sigmaY,sigmaZ,FWHM;
  char *quest[16];
  char msg[100];
  grphic3d *imdeb;
  
  im_deb=GET_PLACE3D(TEXT0011);
  im_res=GET_PLACE3D(TEXT0006);

  /*  Initialisation du filtre */
  imdeb=ptr_img_3d(im_deb);
  FWHM=(double) GET_FLOAT("FWHM (mm)?", 0, &e);
  FWHM=FWHM/(float)2.35482005;      /* s=FWHM/sqrt(8*ln2) en mm */
  sigmaX=FWHM/(float)imdeb->dx; /* dimension de l'image 0 */
  sigmaY=FWHM/(float)imdeb->dy; /* dimension de l'image 0 */
  sigmaZ=FWHM/(float)imdeb->dz; /* dimension de l'image 0 */

  /* Taille filtre   */
  for(i=0;i<9;i++)
    quest[i]=CALLOC(80,char);
  strcpy(quest[0],"3");
  strcpy(quest[1],"5");
  strcpy(quest[2],"7");
  strcpy(quest[3],"9");
  strcpy(quest[4],"11");
  strcpy(quest[5],"13");
  strcpy(quest[6],"15");
  strcpy(quest[7],"\0");

  wndX = (int)(4.*sigmaX*2.+1+0.5);
  if (wndX%2==0)
    wndX+=1;
  wndY = (int)(4.*sigmaY*2.+1+0.5);
  if (wndY%2==0)
    wndY+=1;
  wndZ = (int)(4.*sigmaZ*2.+1+0.5);
  if (wndZ%2==0)
    wndZ+=1;
  sprintf(msg,"The optimal value for windows in (X,Y,Z) is (%d,%d,%d)",wndX,wndY,wndZ);
  PUT_WNDMSG(msg);
  val=GETV_QCM("Mask ?",(char **)quest);
  CLEAR_WNDMSG();
  for(i=0;i<9;i++)
    free(quest[i]);

  if(val==0)
    wndX=wndY=wndZ=3;
  else if(val==1)
    wndX=wndY=wndZ=5;
  else if(val==2)
    wndX=wndY=wndZ=7;
  else if(val==3)
    wndX=wndY=wndZ=9;
  else if(val==4)
    wndX=wndY=wndZ=11;
  else if(val==5)
    wndX=wndY=wndZ=13;
  else if(val==6)
    wndX=wndY=wndZ=15;

  /* Filtrage des images par une gaussienne  */
  imx_copie_3d(im_deb,im_res);
  imx_gaussian_filterModif_3d(im_res,wndX,wndY,wndZ,sigmaX,sigmaY,sigmaZ);

  //imx_iniparaimg_3d(im_res);
  show_picture_3d(im_res);

}

/******************************************************************************
 *********
 *********  ----------------- imx_gaussian_filter
 *********
 ********************************************************************************/
void imx_gaussian_filter_3d(int im_deb, int im_res, int wnd, double sigma)
{
  grphic3d *imdeb,*imres;


  imdeb=ptr_img_3d(im_deb);
  imres=ptr_img_3d(im_res);

  imx_gaussian_filter_3d_p(imdeb,imres,wnd,sigma);

}

/******************************************************************************
 *********
 *********  ----------------- imx_gaussian_filter
 *********
 ********************************************************************************/
void imx_gaussian_filter_3d_p(grphic3d *imdeb, grphic3d *imres, int wnd, double sigma)
{
  grphic3d *inter;
  int i,j,k,p;
  double *filter,sum=0.0,coeff;
  double norm=0.0;
  int width,height,depth;

  width=imdeb->width;
  height=imdeb->height;
  depth=imdeb->depth;

  filter=alloc_dvector(wnd);
  inter=cr_grphic3d(imdeb);

  coeff=pow((double) 2.0*PI,(double)0.5);
  for(i=0;i<wnd;i++)
    /* JP 3/99   filter[i]=(1/(2.0*PI*sigma))*exp(-1.0*((double)(i-wnd/2))*((double)(i-wnd/2))/(2.0*sigma*sigma));*/
    /* 2.50662828 correspond a sqrt(2PI)   */
    filter[i]=(1/(coeff*sigma))*exp(-1.0*((double)(i-wnd/2))*((double)(i-wnd/2))/(2.0*sigma*sigma));
  for(i=0;i<wnd;i++)
    norm+=filter[i];

  for(i=0;i<wnd;i++)
    filter[i]/=norm;

  for(i=0;i<width;i++)
    for(j=0;j<height;j++)
      for(k=0;k<depth;k++)
	inter->mri[i][j][k]=0;

  for(k=wnd/2;k<depth-wnd/2;k++)
    for(j=wnd/2;j<height-wnd/2;j++)
      for(i=wnd/2;i<width-wnd/2;i++)
	{
	  sum=0.0;

	  for(p=-wnd/2;p<=wnd/2;p++)
	    sum+=((double)(imdeb->mri[i-p][j][k]))*filter[p+wnd/2];

	  inter->mri[i][j][k]=((int)floor(0.5+sum));
	}

  for(i=0;i<width;i++)
    for(j=0;j<height;j++)
      for(k=0;k<depth;k++)
	imres->mri[i][j][k]=0;


  for(k=wnd/2;k<depth-wnd/2;k++)
    for(i=wnd/2;i<width-wnd/2;i++)
      for(j=wnd/2;j<height-wnd/2;j++)
	{
	  sum=0.0;

	  for(p=-wnd/2;p<=wnd/2;p++)
	    sum+=((double)(inter->mri[i][j-p][k]))*filter[p+wnd/2];

	  imres->mri[i][j][k]=((int)floor(0.5+sum));
	}

  for(i=wnd/2;i<width-wnd/2;i++)
    for(j=wnd/2;j<height-wnd/2;j++)
      for(k=wnd/2;k<depth-wnd/2;k++)
	{
	  sum=0.0;

	  for(p=-wnd/2;p<=wnd/2;p++)
	    sum+=((double)(imres->mri[i][j][k-p]))*filter[p+wnd/2];

	  inter->mri[i][j][k]=((int)floor(0.5+sum));
	}

  for(i=0;i<width;i++)
    for(j=0;j<height;j++)
      for(k=0;k<depth;k++)
	imres->mri[i][j][k]=inter->mri[i][j][k];



  free_dvector(filter,wnd);
  free_grphic3d(inter);

}

/******************************************************************************
 *********
 *********  ----------------- imx_gaussian_filterModif
 *********
 ********************************************************************************/
void imx_gaussian_filterModif_3d(int im_fil, int wndX, int wndY, int wndZ,
                                 double sigmaX, double sigmaY, double sigmaZ)
{
  grphic3d *imfil;


  imfil=ptr_img_3d(im_fil);
 
  imx_brukermax_update_3d(imfil);
  imx_gaussian_filterModif_3d_p(imfil,wndX,wndY,wndZ,sigmaX,sigmaY,sigmaZ);
  imx_inimaxminpixel_3d_p(imfil);

}

/******************************************************************************
 *********
 *********  ----------------- imx_gaussian_filter modifier
 *********
 ********************************************************************************/
void imx_gaussian_filterModif_3d_p(grphic3d *imfmri,
                                   int wndX       ,
                                   int wndY       ,
                                   int wndZ       ,
                                   double sigmaX ,
                                   double sigmaY ,
                                   double sigmaZ)
{
  int i,j,k;
  int wdth,hght,dpth;
  float	*timecrs;
  double *filter,coeff;
  double norm=0.0;
  int swap;
  /* Debut Initialisation variables */

  wdth = imfmri->width;
  hght = imfmri->height;
  dpth = imfmri->depth;

  /*Construction du filtre pour les X*/
  filter=(double*)malloc(wndX*sizeof(double));
  if (filter == NULL)
    PUT_ERRFAT( ERR_MEM_ALLOC);
  swap = (wndX-1)/2; /*wnd est impair*/
  coeff=pow((double) 2.0*PI,(double)0.5);
  for(i=0;i<wndX;i++)
    filter[i]=(1/(coeff*sigmaX))

              *exp(-1.0*((double)(i-swap))*((double)(i-swap))/((double)2.0*sigmaX*sigmaX));

  for(i=0;i<wndX;i++)
    norm+=filter[i];
  for(i=0;i<wndX;i++)
    filter[i]/=norm;

  /* Dans le sens des X */
  if( (timecrs = (float*)malloc(sizeof(float)*wdth)) == NULL)
  {
    PUT_WNDMSG( ERR_MEM_ALLOC);
    return;
  }

  for(j=0 ; j<hght ; j++)
    for(k=0 ; k<dpth ; k++)
    {
      for (i=0;i<wdth;i++)
        timecrs[i]=(float)imfmri->mri[i][j][k];

      sc_gaussian_filterModif(timecrs,wdth,wndX,sigmaX,filter);

      for (i=0;i<wdth;i++)
        imfmri->mri[i][j][k]=(int)(timecrs[i]+0.5);
    }

  free(timecrs);



  /*Construction du filtre pour les Y*/
  free(filter);
  filter=(double*)malloc(wndY*sizeof(double));
  if (filter == NULL)
    PUT_ERRFAT( ERR_MEM_ALLOC);
  swap = (wndY-1)/2; /*wnd est impair*/
  coeff=pow((double) 2.0*PI,(double)0.5);
  for(i=0;i<wndY;i++)
    filter[i]=(1/(coeff*sigmaY))
              *exp(-1.0*((double)(i-swap))*((double)(i-swap))/((double)2.0*sigmaY*sigmaY));
  for(i=0;i<wndY;i++)
    norm+=filter[i];
  for(i=0;i<wndY;i++)
    filter[i]/=norm;


  /* Dans le sens des Y */
  if( (timecrs = (float*)malloc(sizeof(float)*hght)) == NULL)
  {
    PUT_WNDMSG( ERR_MEM_ALLOC);
    return;
  }
  for(i=0 ; i<wdth ; i++)
    for(k=0 ; k<dpth ; k++)
    {
      for (j=0;j<hght;j++)
        timecrs[j]=imfmri->mri[i][j][k];

      sc_gaussian_filterModif(timecrs,hght,wndY,sigmaY,filter);

      for (j=0;j<hght;j++)
        imfmri->mri[i][j][k]=(int)(timecrs[j]+0.5);
    }

  free(timecrs);


  /*Construction du filtre pour les Z*/
  free(filter);
  filter=(double*)malloc(wndZ*sizeof(double));
  if (filter == NULL)
    PUT_ERRFAT( ERR_MEM_ALLOC);
  swap = (wndZ-1)/2;
  coeff=pow((double) 2.0*PI,(double)0.5);
  for(i=0;i<wndZ;i++)
    filter[i]=(1/(coeff*sigmaZ))
              *exp(-1.0*((double)(i-swap))*((double)(i-swap))/((double)2.0*sigmaZ*sigmaZ));
  for(i=0;i<wndZ;i++)
    norm+=filter[i];
  for(i=0;i<wndZ;i++)
    filter[i]/=norm;


  /* Dans le sens des Z */
  if( (timecrs = (float*)malloc(sizeof(float)*dpth)) == NULL)
  {
    PUT_WNDMSG( ERR_MEM_ALLOC);
    return;
  }
  for(i=0 ; i<wdth ; i++)
    for(j=0 ; j<hght ; j++)
    {
      for (k=0;k<dpth;k++)
        timecrs[k]=imfmri->mri[i][j][k];

      sc_gaussian_filterModif(timecrs,dpth,wndZ,sigmaZ,filter);

      for (k=0;k<dpth;k++)
        imfmri->mri[i][j][k]=(int)(timecrs[k]+0.5);
    }

  free(timecrs);
  free(filter);
}

void sc_gaussian_filterModif(float *y, long int nbpts, int wnd, double sigma,double * filter)
{
  int /*unused variable i,*/l,p,swap;
  double norm=0.0;
  float *yres;

  swap = (wnd-1)/2; /*wnd est impair*/

  /* Allocation du tableau de  yres */
  if( (yres = (float*)malloc(sizeof(float)*nbpts)) == NULL)
  {
    PUT_WNDMSG( ERR_MEM_ALLOC);
    return;
  }


  for (l=0;l<nbpts;l++)
    yres[l]=0.;

  for (l=0;l<swap;l++)    /*Probleme de bord a gauche*/
  {
    if (y[l]!=0)
    {
      norm = 0.;
      for(p=-swap;p<=swap;p++)
      {
        if ((l+p)>=0)
        {
          if (y[l+p]!=0)
          {
            yres[l]+=(double)y[l+p]*filter[p+swap];
            norm += filter[p+swap];
          }
        }
      }
      if (norm>0)
        yres[l]/=norm;
    }
  }

  for (l=nbpts-swap;l<nbpts;l++)   /*Probleme de bord a droite*/
  {
    if (y[l]!=0)
    {
      norm = 0.;
      for(p=-swap;p<=swap;p++)
      {
        if ((l+p)<nbpts)
        {
          if (y[l+p]!=0)
          {
            yres[l]+=(double)y[l+p]*filter[p+swap];
            norm += filter[p+swap];
          }
        }
      }
      if (norm>0)
        yres[l]/=norm;
    }
  }

  for (l=swap;l<nbpts-swap;l++)
  {
    if (y[l]!=0)
    {
      norm = 0;
      for(p=-swap;p<=swap;p++)
      {
        if (y[l+p]!=0)
        {
          yres[l]+=(double)y[l+p]*filter[p+swap];
          norm+=filter[p+swap];
        }
      }
      if (norm>0)
        yres[l]/=norm;
    }
  }



  for (l=0;l<nbpts;l++)
    y[l]=yres[l];


  free(yres);
}

/******************************************************************************
 **    moyen_filter_3d() :
 **
 */
/*! filtrage moyenneur
************************************************************************/
void moyen_filter_3d(void)
{
  int im_deb,im_res,wnd=0;

  im_deb=GET_PLACE3D(TEXT0011);
  im_res=GET_PLACE3D(TEXT0006);

  wnd=imx_query_filter_dimension();

  imx_moyen_filter_3d(im_deb,im_res,wnd);

  //imx_iniparaimg_3d(im_res);
  show_picture_3d(im_res);
}

/******************************************************************************
 **
 **   imx_moyen_filter_3d()
 */
/*!   filtrage moyenneur
**	 \param im_deb : numero de l'image source
**   \param im_res : numero de l'image resultat
**   \param wnd : taille de la fenetre
**
*************************************************************************/


int imx_moyen_filter_3d(int im_deb, int im_res, int wnd)
{
  grphic3d *imdeb,*imres;
  int err=0;

  imdeb=ptr_img_3d(im_deb);
  imres=ptr_img_3d(im_res);

  err=imx_moyen_filter_3d_p(imdeb,imres,wnd);

  return err;
}


/******************************************************************************
 **
 **  imx_moyen_filter_3d_p()
 */
/*!  \ingroup ToolBoxImage3d
**
**  filtrage moyenneur
**  \param imdeb :image source
**  \param imres : image resultat (E/S)
**  \param wnd : taille de la fenetre
**
****************************************************************************/
int imx_moyen_filter_3d_p(grphic3d *imdeb, grphic3d *imres, int wnd)
{
  int i,j,k,l,m,n,p;
  int width,height,depth;
  double val;
  int err=0;
  int i1, i2, j1, j2, k1, k2;

  TYPEMRI3D ***imdebMRI=NULL, ***imtempMRI=NULL;
  grphic3d *imtemp;

 if ((!imdeb)||(!imres))
  { fprintf(stderr,"les deux images doivent etre non nulles dans imx_median_3d_p\n"); return 1; }

 if (wnd%2!=1)
 { fprintf(stderr,"la taille du filtre doit etre impaire dans imx_median_3d_p\n"); return 2; }


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
     p=0;
	val=0;
     k1=MAXI(k-wnd, 0); k2=MINI(depth-1, k+wnd);

     for(l=i1;l<=i2;l++)
      for(m=j1;m<=j2;m++)
       for(n=k1;n<=k2;n++)
       {
       val+=imdebMRI[l][m][n];
        p++;
       }
       imtempMRI[i][j][k]=(TYPEMRI3D)(val/(float)p);
    }
   }
  }

 imx_copie_3d_p(imtemp,imres);
 imx_inimaxminpixel_3d_p(imres);
 free_grphic3d(imtemp);

end_func:


 return err;
}






/******************************************************************************
 **    median_3d() :
 **
 */
/*! filtrage median
************************************************************************/
void median_3d(void)
{
  int im_deb,im_res,wnd=0;

  im_deb=GET_PLACE3D(TEXT0011);
  im_res=GET_PLACE3D(TEXT0006);

  wnd=imx_query_filter_dimension();

  imx_median_3d(im_deb,im_res,wnd);

  //imx_iniparaimg_3d(im_res);
  show_picture_3d(im_res);
}

/******************************************************************************
 **
 **   imx_median_3d()
 */
/*!   filtrage median
**	 \param im_deb : numero de l'image source
**   \param im_res : numero de l'image resultat
**   \param wnd : taille de la fenetre
**
*************************************************************************/


int imx_median_3d(int im_deb, int im_res, int wnd)
{
  grphic3d *imdeb,*imres;
  int err=0;

  imdeb=ptr_img_3d(im_deb);
  imres=ptr_img_3d(im_res);

  err=imx_median_3d_p(imdeb,imres,wnd);

  return err;
}


/******************************************************************************
 **
 **  imx_median_3d_p()
 */
/*!  \ingroup ToolBoxImage3d
**
**  filtrage median
**  \param imdeb :image source
**  \param imres : image resultat (E/S)
**  \param wnd : taille de la fenetre
**
****************************************************************************/
int imx_median_3d_p(grphic3d *imdeb, grphic3d *imres, int wnd)
{
  int i,j,k,l,m,n,p;
  int *tab;
  int width,height,depth;
  double val;
  int err=0;
  int i1, i2, j1, j2, k1, k2;
  TYPEMRI3D ***imdebMRI=NULL, ***imtempMRI=NULL;
  grphic3d *imtemp;

 if ((!imdeb)||(!imres))
  { fprintf(stderr,"les deux images doivent etre non nulles dans imx_median_3d_p\n"); return 1; }

 if (wnd%2!=1)
 { fprintf(stderr,"la taille du filtre doit etre impaire dans imx_median_3d_p\n"); return 2; }

  tab=alloc_ivector(wnd*wnd*wnd);
  if (!tab) { fprintf (stderr, "erreur allocation memoire dans imx_median_3d_p\n"); goto end_func; }

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
     p=0;

     k1=MAXI(k-wnd, 0); k2=MINI(depth-1, k+wnd);

     for(l=i1;l<=i2;l++)
      for(m=j1;m<=j2;m++)
       for(n=k1;n<=k2;n++)
       {
        tab[p]=(int)(imdebMRI[l][m][n]);
        p++;
       }
       val=quick_select_integer(((i2-i1+1)*(j2-j1+1)*(k2-k1+1)/2)+1,tab,(i2-i1+1)*(j2-j1+1)*(k2-k1+1));
       imtempMRI[i][j][k]=(TYPEMRI3D)val;
    }
   }
  }

 imx_copie_3d_p(imtemp,imres);
 imx_inimaxminpixel_3d_p(imres);
 free_grphic3d(imtemp);

end_func:

 if (tab) free_ivector(tab,wnd*wnd*wnd);

 return err;
}


/* -- imx_diffusion_filter_3d_p() --------------------------------------------
 *    Filtre approximant une equation de diffusion
 *
 *    Usage   imx_diffusion_filter_dif_3d_p()
 */
void
imx_diffusion_filter_3d_p(grphic3d *imdeb, grphic3d *imres, unsigned int nb_it)
{
  unsigned int width, height, depth;
  unsigned int index, i, j, k;
  grphic3d *imtemp;
  unsigned char binaire;

  imtemp = cr_grphic3d(imdeb);
  imx_copie_param_3d_p(imdeb, imtemp);

  width = imdeb->width;
  height = imdeb->height;
  depth = imdeb->depth;

  binaire = (imdeb->max_pixel == 1);

  if (binaire){
    for (i=0; i<width; i++)
      for (j=0; j<height; j++)
	for (k=0; k<depth; k++)
	  imtemp->mri[i][j][k] = (imdeb->mri[i][j][k] == 0) ? 0 : 32767;
  }else{
    for (i=0; i<width; i++)
      for (j=0; j<height; j++)
	for (k=0; k<depth; k++)
	  imtemp->mri[i][j][k] = imdeb->mri[i][j][k];
  }

  /* On copie tout dans l'image resultat */
  imx_copie_3d_p(imtemp, imres);

  /* On applique le filtre de diffusion directement sur l'image */
  for (index=0; index<nb_it; index++){
    for (i=1; i<(width-1); i++)
      for (j=1; j<(height-1); j++)
	for (k=1; k<(depth-1); k++)
	  imres->mri[i][j][k] =
	    imres->mri[i][j][k]/2 + (imres->mri[i+1][j][k] + imres->mri[i-1][j][k]
				     + imres->mri[i][j+1][k] + imres->mri[i][j-1][k]
				     + imres->mri[i][j][k+1] + imres->mri[i][j][k-1])/12;
  }

  if (binaire){
    for (i=0; i<width; i++)
      for (j=0; j<height; j++)
	for (k=0; k<depth; k++)
	  imres->mri[i][j][k] = (imres->mri[i][j][k] > 16383) ? 1 : 0 ;
  }

  free_grphic3d(imtemp);
}



/* -- imx_diffusion_filter_3d() --------------------------------------------
**    Filtre par diffusion
**
**    Usage   imx_filtre_dif_3d(im_deb, im_res, nb_it)
*/
void
imx_diffusion_filter_3d(int im_deb, int im_res, unsigned int nb_it)
{
  grphic3d *imdeb, *imres;

  imdeb = ptr_img_3d(im_deb);
  imres = ptr_img_3d(im_res);

  imx_diffusion_filter_3d_p(imdeb, imres, nb_it);
}



/* -- diffusion_filter_3d() --------------------------------------------
**    Fonction d'entree pour le filtre par diffusion
**
**    Usage   diffusion_filter_3d()
*/
void diffusion_filter_3d(void)
{
  int im_deb, im_res;
  int nb_it;
  int err=0;

  im_deb = GET_PLACE3D( TEXT0001);
  im_res = GET_PLACE3D( TEXT0006);

  /* On demande le nombre d'iterations */
  do{
    nb_it = GET_INT(TEXT0007, 0, &err);
    if (nb_it == 0) break;
  }while (nb_it < 1);


  imx_diffusion_filter_3d(im_deb, im_res, nb_it);

  imx_iniparaimg_3d(im_res);
  show_picture_3d(im_res);
}




/*************************************************************
 ** sub_3d()
 **
 ** difference pixels a pixels entre deux images
 **
 *****************************************************************/

void  sub_3d(void)
{
  int im_1,im_2,im_res;

  im_1=GET_PLACE3D(TEXT0183);
  im_2=GET_PLACE3D(TEXT0184);
  im_res=GET_PLACE3D(TEXT0006);

  imx_sub_3d(im_1,im_2,im_res);
  show_picture_3d(im_res);

  return ;
}

/****subabs_3d()****************************************************
 **
 **
 ** difference absolue entre deux images
 **
 *****************************************************************/

void  subabs_3d(void)
{
  int im_1,im_2,im_res;

  im_1=GET_PLACE3D(TEXT0183);
  im_2=GET_PLACE3D(TEXT0184);
  im_res=GET_PLACE3D(TEXT0006);

  imx_subabs_3d(im_1,im_2,im_res);
  show_picture_3d(im_res);

  return ;
}

/****sub_normamax2_3d()****************************************************
 **
 **
 ** difference normalise par rapport a une image.
 ** Par definition le resultat sera une image comprise entre -255 et +255.
 **
 *****************************************************************/

void  sub_normamax2_3d(void)
{
  int im_1,im_2,im_res;

  im_1=GET_PLACE3D(TEXT0183);
  im_2=GET_PLACE3D(TEXT0184);
  im_res=GET_PLACE3D(TEXT0006);

  imx_sub_normamax2_3d(im_1,im_2,im_res);
  show_picture_3d(im_res);

  return ;
}

/****sub_normpixel2_3d()****************************************************
 **
 **
 ** difference normalise par rapport a une image.
 ** Par definition le resultat sera une image comprise entre -255 et +255.
 **
 *****************************************************************/

void  sub_normpixel2_3d(void)
{
  int im_1,im_2,im_res;

  im_1=GET_PLACE3D(TEXT0183);
  im_2=GET_PLACE3D(TEXT0184);
  im_res=GET_PLACE3D(TEXT0006);

  imx_sub_normpixel2_3d(im_1,im_2,im_res);
  show_picture_3d(im_res);

  return ;
}

/****sub_normpixelmean2_3d()****************************************************
 **
 */
/*! difference normalise par rapport a une image.
** Par definition le resultat sera une image comprise entre -255 et +255.
**
*****************************************************************/

void  sub_normpixelmean2_3d(void)
{
  int im_1,im_2,im_res;

  im_1=GET_PLACE3D(TEXT0183);
  im_2=GET_PLACE3D(TEXT0184);
  im_res=GET_PLACE3D(TEXT0006);

  imx_sub_normpixelmean2_3d(im_1,im_2,im_res);
  show_picture_3d(im_res);

  return;
}

/*******************************************
 ** --  imx_sub_normamax2_3d() ----------------------------
 **     Soustractive normalise par rapport a l'image soustraite
 **     im_1,im_2 : im_1 operation with im_2 /maxpixel(im_2)
 **     im_res : Result image
 ********************************************/
int	imx_sub_normamax2_3d(int im_1, int im_2, int im_res)
{
  grphic3d *im1, *im2, *imres;
  im1=ptr_img_3d(im_1);
  im2=ptr_img_3d(im_2);
  imres=ptr_img_3d(im_res);
  return (imx_sub_normamax2_3d_p(im1, im2, imres));
}

/*******************************************
 ** --  imx_sub_normamax2_3d() ----------------------------
 **     Soustractive normalise par rapport a l'image soustraite
 **     im_1,im_2 : im_1 operation with im_2 /maxpixel(im_2)
 **     im_res : Result image
 ********************************************/
int	imx_sub_normamax2_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{
  int i,j,k,wdth,hght,dpth;
  long icomp;
  float amp1,amp2,ampres;
  float amax2,amaxres,aminres;
  float rcoeff,rcoeff1,rcoeff2;


  hght=MINI(im1->height,im2->height);
  wdth=MINI(im1->width,im2->width);
  dpth=MINI(im1->depth,im2->depth);
  rcoeff1=im1->rcoeff;
  rcoeff2=im2->rcoeff;

  amax2=(float) im2->max_pixel*rcoeff2;

  /*   Recherche amaxres et aminres  */
  amaxres=0.; aminres=(float) LONGMAX;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++) {
        amp1=rcoeff1*im1->mri[i][j][k];
        amp2=rcoeff2*im2->mri[i][j][k];
        ampres=amp1-amp2;
        amaxres=MAXI(amaxres,ampres);
        aminres=MINI(aminres,ampres);
      }
  amaxres=(float)(100.*amaxres/amax2);
  aminres=(float)(100.*aminres/amax2);

  imx_brukermax_3d(amaxres,aminres,imres);
  rcoeff=imres->rcoeff;
  icomp=(long)floor(imres->icomp);

  /*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=rcoeff1*im1->mri[i][j][k];
	  amp2=rcoeff1*im2->mri[i][j][k];
	  ampres= (float)(100.*(amp1-amp2)/(amax2*rcoeff));
	  imres->mri[i][j][k]=(TYPEMRI3D)floor(ampres);
	}

  imx_copie_param_3d_p(im1,imres);
  imres->max_pixel=(long)floor(amaxres/rcoeff);
  imres->min_pixel=(long)floor(aminres/rcoeff);
  imres->cutoff_max=(int)floor(amaxres/rcoeff);
  imres->cutoff_min=(int)floor(aminres/rcoeff);
  imres->icomp=(float)icomp;
  imres->rcoeff=rcoeff;
  if (MESG_DEBUG) printf(" fin sub\n");
  return(1);

}

/*******************************************
 ** --  imx_sub_normpixel2_3d() ----------------------------
 **     Soustractive normalise pixel a pixel par rapport a l'image soustraite
 **     im_1,im_2 : im_1 operation with im_2 / im_2
 **     im_res : Result image
 ********************************************/
int	imx_sub_normpixel2_3d(int im_1, int im_2, int im_res)
{
  grphic3d *im1,*im2,*imres;
  int i,j,k,wdth,hght,dpth;
  int err;
  long icomp;
  float amp1,amp2,ampres,amaxres,aminres;
  float rcoeff,rcoeff1,rcoeff2;

  im1=ptr_img_3d(im_1);
  im2=ptr_img_3d(im_2);
  imres=ptr_img_3d(im_res);

  hght=MINI(im1->height,im2->height);
  wdth=MINI(im1->width,im2->width);
  dpth=MINI(im1->depth,im2->depth);
  rcoeff1=im1->rcoeff;
  rcoeff2=im2->rcoeff;

  /*   Recherche amaxres et aminres  */
  amaxres=0.; aminres=(float) LONGMAX;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++) {
        amp1=rcoeff1*im1->mri[i][j][k];
        amp2=rcoeff2*im2->mri[i][j][k];
        ampres=(amp2==0) ? (float) 1. : (float)((amp1-amp2)/amp2);
        amaxres=MAXI(amaxres,ampres);
        aminres=MINI(aminres,ampres);
      }
  amaxres=(float)(100.*amaxres);
  aminres=(float)(100.*aminres);

  err=imx_brukermax_3d(amaxres,aminres,imres);
  rcoeff=imres->rcoeff;
  icomp=(long)imres->icomp;

  /*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++) {
	amp1=rcoeff1*im1->mri[i][j][k];
	amp2=rcoeff2*im2->mri[i][j][k];
	if (im2->mri[i][j][k]==0)
	  ampres=(float)(100./rcoeff);
	else
	  ampres=(float)( 100.*(amp1-amp2)/(amp2*rcoeff));
	imres->mri[i][j][k]=(TYPEMRI3D)floor(ampres);
      }

  imx_copie_param_3d_p(im1,imres);
  imres->max_pixel=(long)floor(amaxres/rcoeff);
  imres->min_pixel=(long)floor(aminres/rcoeff);
  imres->cutoff_max=(int)floor(amaxres/rcoeff);
  imres->cutoff_min=(int)floor(aminres/rcoeff);
  imres->icomp=(float)icomp;
  imres->rcoeff=rcoeff;
  if (MESG_DEBUG) printf(" fin sub\n");
  return(1);

}

/*******************************************
 ** --  imx_sub_normpixelmean2_3d() ----------------------------
 */
/*!    Soustractive normalise pixel a pixel par rapport a la valeur du pixel moyen
**     \param im_1,im_2 : im_1 operation with im_2 /maxpixel(im_2)
**     \param im_res : Result image
**	   \retval 1
**	   \remark : modifie l'image correspondant a im_res
********************************************/
int	imx_sub_normpixelmean2_3d_p(grphic3d * im1, grphic3d * im2, grphic3d * imres)
{
  int i,j,k,wdth,hght,dpth;
  int err,nb;
  long icomp;
  float amp1,amp2,ampres;
  float amax2,amaxres=0,aminres=0;
  float rcoeff,rcoeff1,rcoeff2;


  hght=MINI(im1->height,im2->height);
  wdth=MINI(im1->width,im2->width);
  dpth=MINI(im1->depth,im2->depth);
  rcoeff1=im1->rcoeff;
  rcoeff2=im2->rcoeff;

  /* Recherche du pixel moyen et de amaxres et aminres   */
  amax2=0.;   nb=0;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++) {
        amp1=rcoeff1*im1->mri[i][j][k];
        amp2=rcoeff2*im2->mri[i][j][k];
        ampres=amp1-amp2;
        amaxres=MAXI(amaxres,ampres);
        aminres=MINI(aminres,ampres);
        if (im2->mri[i][j][k]!=0) {
	  amax2=amax2+(float) im2->mri[i][j][k];
	  nb=nb+1;
	}
      }
  amax2=(float) rcoeff2*amax2/nb;
  amaxres=(float)(100.*amaxres/amax2);
  aminres=(float)(100.*aminres/amax2);

  err=imx_brukermax_3d(amaxres,aminres,imres);
  rcoeff=imres->rcoeff;
  icomp=(long)imres->icomp;

  /*   Calcul de l'image   */
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	{
	  amp1=rcoeff1*im1->mri[i][j][k];
	  amp2=rcoeff2*im2->mri[i][j][k];
	  ampres= (float)(100.*(amp1-amp2)/(amax2*rcoeff));
	  imres->mri[i][j][k]=(TYPEMRI3D)floor(ampres);
	}

  imx_copie_param_3d_p(im1,imres);
  imres->max_pixel=(int)floor(amaxres/rcoeff);
  imres->min_pixel=(int)floor(aminres/rcoeff);
  imres->cutoff_max=(int)floor(amaxres/rcoeff);
  imres->cutoff_min=(int)floor(aminres/rcoeff);
  imres->icomp=(float)icomp;
  imres->rcoeff=rcoeff;
  if (MESG_DEBUG) printf(" fin sub\n");
  return(1);

}

/*******************************************
 ** --  imx_sub_normpixelmean2_3d() ----------------------------
 */
/*!    Soustractive normalise pixel a pixel par rapport a la valeur du pixel moyen
**     \param im_1,im_2 : im_1 operation with im_2 /maxpixel(im_2)
**     \param im_res : Result image
**	   \retval 1
**	   \remark : modifie l'image correspondant a im_res
********************************************/
int	imx_sub_normpixelmean2_3d(int im_1, int im_2, int im_res)
{
  grphic3d *im1,*im2,*imres;

  im1=ptr_img_3d(im_1);
  im2=ptr_img_3d(im_2);
  imres=ptr_img_3d(im_res);

  imx_sub_normpixelmean2_3d_p(im1, im2, imres);

  return(1);
}



/*****compmanuel_3d()************************************************
 **    Methode de vraisemblance
 **    methode statistique pour eliminer le bruit
 **
 **
 **   modelisation constante par zone I=cte+bruit
 **   Choix du seuil
 **
 *****************************************************************/

int compmanuel_3d(void)
{
  int im_1,im_2,im_res;
  int N;
  double seuil;
  int err=0;


  im_1=GET_PLACE3D(TEXT0185);
  im_2=GET_PLACE3D(TEXT0186 );
  im_res=GET_PLACE3D(TEXT0006);
  N= GET_INT(TEXT0187, 0, &err);

  if (N>0) {
    seuil= GET_DOUBLE(TEXT0188, 0, &err);
    imx_compmanuel_3d(im_1,im_2,im_res,N,seuil);
    show_picture_3d(im_res);
  }
  else PUT_WNDMSG(TEXT0189);


  return 0;
}

/****imx_compmanuel_3d()*********************************************
 **  im_1,im_2: comparaison des images im_1 et im_2
 **  im_res: image resultat
 **  N: taille de la fenetre d'operation
 **  seuil: valeur du seuil de difference entre les 2 images
 *****************************************************************/
int imx_compmanuel_3d(int im_1, int im_2, int im_res, int N, double seuil)
{
  grphic3d *im1,*im2,*imres;

  im1=ptr_img_3d(im_1);
  im2=ptr_img_3d(im_2);
  imres=ptr_img_3d(im_res);

  imx_compmanuel_3d_p(im1,im2,imres,N,seuil);

  return 0;
}


/****imx_compmanuel_3d_p()**********************************************
 ** Sous programme
 **  Compare la difference entre deux images en utilisant les
 **  pointeurs d'image
 *****************************************************************/
int imx_compmanuel_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres, int N, double seuil)
{
  grphic3d *imrestemp;
  int width,height,depth;
  int w;
  int i,j,k,f,g,h;
  double moy1,moy2;
  double diffmoy;
  double sigma1,sigma2;          /* variance */
  double var1,var2;
  int width2,height2,depth2;
  double rcoeff1,rcoeff2;
  double rcoeff;
  int N3;
  double M;
  double aux1,aux2,aux3;
  double max=0,min=LONGMAX;
  double diff;
  double I1,I2;
  int inf1,inf2,inf3,sup1,sup2,sup3;
  int err;
  double R;

  imrestemp=cr_grphic3d(im1);
  width=im1->width;
  height=im1->height;
  depth=im1->depth;
  w=N/2;
  width2=width-w;
  height2=height-w;
  depth2=depth-w;
  rcoeff1=(float)im1->rcoeff;
  rcoeff2=(float)im2->rcoeff;
  N3=N*N*N;
  M= pow((double) N,1.5)/2;


  for (i=0;i<width;i++)
    for (j=0;j<height;j++)
      for (k=0;k<depth;k++)
	{
	  imrestemp->mri[i][j][k]=0;

	  diff=(double) rcoeff1*(im1->mri[i][j][k])-rcoeff2*(im2->mri[i][j][k]);
	  if (diff>max) max=diff;
	  min=MINI(min,diff);
	}


  /*   Calcul de imres->maxpixel et imres->icomp et imres->rcoeff a partir de max*/
  err=imx_brukermax_3d((float)max,(float)min,imrestemp);
  rcoeff=imrestemp->rcoeff;


  for(i=w;i<width2;i++)
    for(j=w;j<height2;j++)
      for (k=w;k<depth2;k++)
	{
	  moy1=0;
	  moy2=0;
	  sigma2=0;
	  sigma1=0;
	  var1=0;
	  var2=0;
	  diff=0;
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
		  moy1 += rcoeff1*(im1->mri[f][g][h]);
		  moy2 += rcoeff2*(im2->mri[f][g][h]);
		}
	  moy1=moy1/(N3);
	  moy2=moy2/(N3);
	  diffmoy=fabs(moy1-moy2);
	  for(f=inf1;f<=sup1;f++)
	    for(g=inf2;g<=sup2;g++)
	      for(h=inf3;h<=sup3;h++)
		{
		  I1=rcoeff1*(im1->mri[f][g][h])-moy1;
		  I2=rcoeff2*(im2->mri[f][g][h])-moy2;
		  var1 += I1*I1;
		  var2 += I2*I2;
		}
	  var1=var1/N3;
	  var2=var2/N3;
	  sigma1=sqrt(var1);
	  sigma2=sqrt(var2);

	  diff=rcoeff1*(im1->mri[i][j][k])-rcoeff2*(im2->mri[i][j][k]);
	  diff= diff/rcoeff;

	  if ((sigma1==0 || sigma2==0) && diffmoy!=0)  imrestemp->mri[i][j][k]=(TYPEMRI3D)floor(diff);
	  else
	    if ((sigma1==0 || sigma2==0) && diffmoy==0) imrestemp->mri[i][j][k]=0;
	    else {
	      aux1=sigma1/(2*sigma2);
	      aux2=sigma2/(2*sigma1);
	      aux3=(diffmoy*diffmoy)/(2*sigma1*sigma2);
	      R=aux1+aux2+aux3;
	      if (R>seuil) imrestemp->mri[i][j][k]=(TYPEMRI3D)floor(diff);
	    }

	}


  imx_copie_param_3d_p(imrestemp,imres);

  for (i=0;i<width;i++)
    for (j=0;j<height;j++)
      for (k=0;k<depth;k++) imres->mri[i][j][k]=imrestemp->mri[i][j][k];



  imres->min_pixel=0;
  imres->width=width;
  imres->height=height;
  imres->depth=depth;



  free_grphic3d(imrestemp);

  return 0;
}


/******compauto_3d()************************************************
 **    Methode de vraisemblance 2
 **    methode statistique pour eliminer le bruit
 **
 **
 **   modelisation constante par zone I=cte+bruit
 **   Seuil AUTOMATIQUE
 **
 *****************************************************************/

int compauto_3d(void)
{
  int im_1,im_2,im_res;
  int N;
  float K=0.0;
  char *quest[3];
  int answer_nr;
  int err=0,i;

  for(i=0;i<3;i++)
    quest[i]=CALLOC(80,char);
  strcpy(quest[0],"ecarts types differents dans les 2 fenetres" );
  strcpy(quest[1],"meme ecart type dans les 2 fenetres");
  strcpy(quest[2],"\0");

  im_1=GET_PLACE3D(TEXT0185);
  im_2=GET_PLACE3D(TEXT0186);
  im_res=GET_PLACE3D(TEXT0006);
  N= GET_INT(TEXT0187, 0, &err);
  /* Put a question box at position 150,200. */
  answer_nr=GETV_QCM(TEXT0043,(char **)quest);

  if (answer_nr==0)
    K= GET_FLOAT("ckoix de K tel que (mu1-mu2)>K(sigma1+sigma2)/2", 0, &err);
  if (answer_nr==1)
    K= GET_FLOAT("ckoix de K tel que (mu1-mu2)>K*sigma", 0, &err);

  if (N>0) {
    imx_compauto_3d(im_1,im_2,im_res,N,answer_nr,K);
    show_picture_3d(im_res);
  }
  else PUT_WNDMSG(TEXT0189);

  for(i=0;i<3;i++) FREE(quest[i]);
  return 0;
}

/*****imx_compauto_3d()********************************************
 **  im_1,im_2: comparaison des images im_1 et im_2
 **  im_res: image resultat
 **  N: taille de la fenetre d'operation
 *****************************************************************/
int imx_compauto_3d(int im_1, int im_2, int im_res, int N, int answer_nr, float K)
{
  grphic3d *im1,*im2,*imres;

  im1=ptr_img_3d(im_1);
  im2=ptr_img_3d(im_2);
  imres=ptr_img_3d(im_res);


  imx_compauto_3d_p(im1,im2,imres,N,answer_nr,K);

  return 0;
}

/*****imx_compauto_3d_p()*******************************************
 ** Sous programme
 **  Compare la difference entre deux images en utilisant les
 **  pointeurs d'image
 *****************************************************************/
int imx_compauto_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres, int N, int answer_nr, float K)
{
  grphic3d *imrestemp;
  int width,height,depth;
  int w;
  int i,j,k,f,g,h;
  double moy1,moy2;
  double diffmoy;
  double sigma1,sigma2;          /* variance */
  double var1,var2;
  double cond;
  int width2,height2,depth2;
  double rcoeff1,rcoeff2;
  double rcoeff;
  int N3;
  double M;
  double aux1,aux2,aux3;
  double max=0,min=LONGMAX;
  double diff;
  double I1,I2;
  int inf1,inf2,inf3,sup1,sup2,sup3;
  int err;
  double R;


  imrestemp=cr_grphic3d(im1);
  width=im1->width;
  height=im1->height;
  depth=im1->depth;
  w=N/2;
  width2=width-w;
  height2=height-w;
  depth2=depth-w;
  rcoeff1=(float)im1->rcoeff;
  rcoeff2=(float)im2->rcoeff;
  N3=N*N*N;
  M= pow((double) N,1.5)/2;


  for (i=0;i<width;i++)
    for (j=0;j<height;j++)
      for (k=0;k<depth;k++)
	{
	  imrestemp->mri[i][j][k]=0;

	  diff=(double) rcoeff1*(im1->mri[i][j][k])-rcoeff2*(im2->mri[i][j][k]);
	  if (diff>max) max=diff;
	  min=MINI(min,diff);
	}

  /*   Calcul de imres->maxpixel et imres->icomp et imres->rcoeff a partir de max*/
  err=imx_brukermax_3d((float)max,(float)min,imrestemp);
  rcoeff=imrestemp->rcoeff;



  for(i=w;i<width2;i++)
    for(j=w;j<height2;j++)
      for (k=w;k<depth2;k++)
	{
	  moy1=0;
	  moy2=0;
	  sigma2=0;
	  sigma1=0;
	  var1=0;
	  var2=0;
	  diff=0;
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
		  moy1 += rcoeff1*(im1->mri[f][g][h]);
		  moy2 += rcoeff2*(im2->mri[f][g][h]);
		}
	  moy1=moy1/(N3);
	  moy2=moy2/(N3);
	  diffmoy=fabs(moy1-moy2);
	  for(f=inf1;f<=sup1;f++)
	    for(g=inf2;g<=sup2;g++)
	      for(h=inf3;h<=sup3;h++)
		{
		  I1=rcoeff1*(im1->mri[f][g][h])-moy1;
		  I2=rcoeff2*(im2->mri[f][g][h])-moy2;
		  var1 += I1*I1;
		  var2 += I2*I2;
		}
	  var1=var1/N3;
	  var2=var2/N3;
	  sigma1=sqrt(var1);
	  sigma2=sqrt(var2);

	  diff=rcoeff1*(im1->mri[i][j][k])-rcoeff2*(im2->mri[i][j][k]);
	  diff= diff/rcoeff;

	  switch (answer_nr) {
	  case 1: if (sigma1!= 0) R=M*diffmoy/sigma1;
	  else R=0;
	    cond=M*K;
	    if (R>cond) imrestemp->mri[i][j][k]=(TYPEMRI3D)floor(diff);
	    break;
	  case 0: if ((sigma1==0 || sigma2==0) && diffmoy!=0)  imrestemp->mri[i][j][k]=(TYPEMRI3D)floor(diff);
	  else
	    if ((sigma1==0 || sigma2==0) && diffmoy==0) imrestemp->mri[i][j][k]=0;
	    else {
	      aux1=sigma1/(2*sigma2);
	      aux2=sigma2/(2*sigma1);
	      aux3=(diffmoy*diffmoy)/(2*sigma1*sigma2);
	      R=aux1+aux2+aux3;
	      cond=aux1+aux2+K*K*(1+aux1+aux2)/4 ;
	      if (R>cond) imrestemp->mri[i][j][k]=(TYPEMRI3D)floor(diff);
	    }
	    break;
	  default: break;
	  }
	}

  imx_copie_param_3d_p(imrestemp,imres);

  for (i=0;i<width;i++)
    for (j=0;j<height;j++)
      for (k=0;k<depth;k++) imres->mri[i][j][k]=imrestemp->mri[i][j][k];


  imres->min_pixel=0;
  imres->width=width;
  imres->height=height;
  imres->depth=depth;

  free_grphic3d(imrestemp);
  return 0;
}



/********pix_compte_3d()********************************************
 ***
 ***	Compte le nombre de pixels de l'image NON NULS
 ***
 *****************************************************************/

void pix_compte_3d(void)
{
  int im_1;
  grphic3d *im1;
  int i,j,k;
  int width,height,depth;
  int compte;

  im_1=GET_PLACE3D("Image a considerer");
  im1=ptr_img_3d(im_1);

  width=im1->width;
  height=im1->height;
  depth=im1->depth;
  compte=0;

  for (i=0;i<width;i++)
    for (j=0;j<height;j++)
      for (k=0;k<depth;k++)
	{
	  if (im1->mri[i][j][k] != 0.0) compte++;
	}

  PUTI("Nombre de pixels non nuls  =",compte," ");

  return ;
}




/*************************** student_3d() *******************************/
/*                                                                      */
/*                          Test de Student                             */
/*                                                                      */
/*      Equivalent au test de vraisemblance g��alis�                 */
/*      avec mod�isation constante de la fonction intensit�           */
/*      et utilisant la 'variance en commun' comme m�e estimation      */
/*      des variances dans les deux fen�res                            */
/*                                                                      */
/************************************************************************/

void student_3d(void)
{
  int im_1,im_2,im_res;
  int N;
  float fausse_alarme;
  char *quest[3];
  int err=0,i,type_remplissage;

  for(i=0;i<3;i++) quest[i]=CALLOC(80,char);
  strcpy(quest[0],TEXT0227);
  strcpy(quest[1],TEXT0228);
  strcpy(quest[2],"\0");

  im_1=GET_PLACE3D(TEXT0185);
  im_2=GET_PLACE3D(TEXT0186);
  im_res=GET_PLACE3D(TEXT0006);

  while((N= GET_INT(TEXT0187, 0, &err))<=0)
    PUT_WNDMSG(TEXT0189);
  CLEAR_WNDMSG();

  while((fausse_alarme= GET_FLOAT(TEXT0221, 0, &err))<0 || fausse_alarme>1)
    PUT_WNDMSG(TEXT0222);
  CLEAR_WNDMSG();

  type_remplissage=GETV_QCM(TEXT0218,(char **)quest);
  imx_student_3d(im_1,im_2,im_res,N,fausse_alarme,type_remplissage);
  show_picture_3d(im_res);

  for(i=0;i<3;i++) FREE(quest[i]);
  return ;
}

/*****imx_student_3d()********************************************
 **  im_1,im_2: comparaison des images im_1 et im_2
 **  im_res: image resultat
 **  N: taille de la fenetre d'operation
 *****************************************************************/
int imx_student_3d(int im_1, int im_2, int im_res, int N, float fausse_alarme, int type_remplissage)
{
  grphic3d *im1,*im2,*imres;

  im1=ptr_img_3d(im_1);
  im2=ptr_img_3d(im_2);
  imres=ptr_img_3d(im_res);

  imx_student_3d_p(im1,im2,imres,N,fausse_alarme,type_remplissage);

  return 0;
}

/*****imx_student_3d_p()*******************************************
 ** Sous programme
 **  Compare la difference entre deux images en utilisant les
 **  pointeurs d'image
 *****************************************************************/
int imx_student_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres, int N, float fausse_alarme, int type_remplissage)
{

  float *mri_temp;
  int width,height,depth;
  int w;
  int i,j,k,f,g,h;
  float moy1,moy2;
  int width2,height2,depth2;
  float rcoeff1,rcoeff2;
  float rcoeff;
  int N3;
  float M;
  float max=0,min=LONGMAX;
  float diff;
  float I1,I2;
  int inf1,inf2,inf3,sup1,sup2,sup3;
  int err;
  float somme;
  float t,t2,seuil2;
  int which,status;
  double p,q,seuil,degree_of_freedom,bound;


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
  M= (float)(sqrt((double) N3*(N3-1)));

  if (( mri_temp = (float *)malloc(sizeof(double)*width*height*depth))==NULL)
    {
      printf("failled to allocate memory\n");
      //return NULL;
      return (-1);   // Il faut un entier !!! RC. Nov 2001
    }

  /* Inverse Student Distribution with (2n-2) degrees of freedom */
  /* Calculates threshold parameter given probability of false alarm */
  which=2;
  q=fausse_alarme/2.0;
  p=1-q;
  degree_of_freedom=N3+N3-2;
  cdft(&which, &p, &q, &seuil, &degree_of_freedom, &status, &bound);
  seuil2=(float)(seuil*seuil);

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
	  moy1=0;
	  moy2=0;
	  diff=0;
	  somme=0;
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
		  moy1 += im1->mri[f][g][h];
		  moy2 += im2->mri[f][g][h];
		}
	  moy1=rcoeff1*moy1/(N3);
	  moy2=rcoeff2*moy2/(N3);

	  for(f=inf1;f<=sup1;f++)
	    for(g=inf2;g<=sup2;g++)
	      for(h=inf3;h<=sup3;h++)
		{
		  I1=rcoeff1*(im1->mri[f][g][h])-moy1;
		  I2=rcoeff2*(im2->mri[f][g][h])-moy2;
		  somme += I1*I1+I2*I2;
		}
	  somme = (float)MAXI(somme,1.0e-9);

	  diff=rcoeff1*(im1->mri[i][j][k])-rcoeff2*(im2->mri[i][j][k]);

	  t=(float)(M*(moy1-moy2)/sqrt(somme));
	  t2=t*t;

	  if (t2>seuil2)
	    {
	      if (type_remplissage == 0)
		{
		  if((mri_temp[i*height2*depth2+j*depth2+k]=diff)>max)
		    max=diff;
		  min=MINI(min,diff);
		}
	      if (type_remplissage == 1)
		{
		  if((mri_temp[i*height2*depth2+j*depth2+k]=(float)log(t2))>max)
		    max=(float)log(t2);
		  min=(float)MINI(min,log(t2));
		}
	    }

       	}

  /* Calcul de imres->maxpixel et imres->icomp et imres->rcoeff a partir de max*/
  err=imx_brukermax_3d(max,min,imres);
  rcoeff=imres->rcoeff;

  for (i=w;i<width2;i++)
    for (j=w;j<height2;j++)
      for (k=w;k<depth2;k++)
	imres->mri[i][j][k]=(TYPEMRI3D)floor(mri_temp[i*height2*depth2+j*depth2+k]/rcoeff);


  imres->min_pixel=0;
  imres->width=width;
  imres->height=height;
  imres->depth=depth;

  free(mri_temp);

  return 0;
}




/************************* matched_filter_3d() *****************************/
/*                 Methode de vraisemblance g��alis�                   */
/*                 'Generalized Likelihood Ratio Test'                    */
/*                                                                        */
/*            Modelisation constante de l'intensit�I=cte+bruit           */
/*                   	  Variance du bruit connue                        */
/**************************************************************************/

void matched_filter_3d(void)
{
  int im_1,im_2,im_res;
  int N;
  double fausse_alarme;
  double sigma,sigma1,sigma2;
  int err=0,i,type_remplissage;
  char *quest[3];

  for(i=0;i<3;i++)
    quest[i]=CALLOC(80,char);

  strcpy(quest[0],TEXT0227);
  strcpy(quest[1],TEXT0228);
  strcpy(quest[2],"\0");

  im_1=GET_PLACE3D(TEXT0185);
  while((sigma1= GET_DOUBLE(TEXT0219, 0, &err))<=0)
    PUT_WNDMSG(TEXT0220);
  CLEAR_WNDMSG();

  im_2=GET_PLACE3D(TEXT0186);
  while((sigma2= GET_DOUBLE(TEXT0219, 0, &err))<=0)
    PUT_WNDMSG(TEXT0220);
  CLEAR_WNDMSG();

  im_res=GET_PLACE3D(TEXT0006);

  while((N= GET_INT(TEXT0187, 0, &err))<=0)
    PUT_WNDMSG(TEXT0189);
  CLEAR_WNDMSG();

  while((fausse_alarme= GET_DOUBLE(TEXT0221, 0, &err))<0 || fausse_alarme>1)
    PUT_WNDMSG(TEXT0222);
  CLEAR_WNDMSG();


  type_remplissage=GETV_QCM(TEXT0218,(char **)quest);

  sigma=sqrt(sigma1*sigma1+sigma2*sigma2);

  imx_matched_filter_3d(im_1,im_2,im_res,N,fausse_alarme,type_remplissage,sigma);
  show_picture_3d(im_res);


  for(i=0;i<3;i++)
    FREE(quest[i]);

  return;
}

/**************************************************************************/
/*                          Comparaison par rang                          */
/**************************************************************************/
void comparaison_par_rang(void) {
  int im_ref, im_cmp, im_res, methode, i;
  char *quest[3];
  for(i=0;i<3;i++)
    quest[i]=CALLOC(40,char);

  strcpy(quest[0],"comparaison par integrale discrete");
  strcpy(quest[1],"comparaison par rang");
  strcpy(quest[2],"\0");

  /* choix des images a traiter */
  printf("comparaison par rang...\n");
  im_ref=GET_PLACE3D(TEXT0185);
  im_cmp=GET_PLACE3D(TEXT0186);
  im_res=GET_PLACE3D(TEXT0006);

  methode=GETV_QCM("Choix de la methode",(char **)quest);

  /* calcul de l'image de comparaison */
  imx_comparaison_par_rang(im_ref, im_cmp, im_res, methode);
}

void imx_comparaison_par_rang(int im_ref, int im_cmp, int im_res, int methode) {
  grphic3d *imref, *imcmp, *imres;

  /* recuperation des pointeurs des images a traiter */
  imref = ptr_img_3d(im_ref);
  imcmp = ptr_img_3d(im_cmp);
  imres = ptr_img_3d(im_res);

  /* calcul de l'image de comparaison */
  if(methode==0)
    imx_comparaison_par_int_dis_p(imref, imcmp, imres);
  else
    imx_comparaison_par_rang_p(imref, imcmp, imres);

  /* raffrechissement de l'affichage de l'image resultat */
  show_picture_3d(im_res);
}

typedef int (*compfn)(const void*, const void*);

typedef struct s_point3d {
  int x, y, z;
  float value;
} point3d;

int compare(point3d *premiere, point3d *deuxieme) {
  float diff = premiere->value - deuxieme->value;
  return(diff>0?1:(diff==0?0:-1));
}

void imx_rechelonage01_p(grphic3d *imref, grphic3d *imres) {
  int i,j,k, cptVal=0;
  int ratio;
  point3d image[256*256*256];

  /* en chaque point de l'image... */
  for(i=0;i<imref->width;i++)
    for(j=0;j<imref->height;j++)
      for(k=0;k<imref->depth;k++) {
        /* si dans le masque */
        if(!imref->mask || imref->mask->mri[i][j][k]!=0) {
          image[cptVal].x = i;
          image[cptVal].y = j;
          image[cptVal].z = k;
          image[cptVal].value = imref->mri[i][j][k] * imref->rcoeff;
          cptVal++;
        }
      }

  /* trie du tableau */
  qsort(
        (void*)&image,
        cptVal,
        sizeof(point3d),
        (compfn)compare);

  /* recopie des parametres de l'image de reference dans l'image resultat */
  imx_copie_param_3d_p(imref, imres);

  /* initialisation du rcoeff */
  ratio = cptVal / MAXMRI3D +1;
  imres->rcoeff = ((float)ratio)/((float)cptVal);

  /* initaialisation de la nouvelle image */
  for(i=0;i<imref->width;i++)
    for(j=0;j<imref->height;j++)
      for(k=0;k<imref->depth;k++)
        imres->mri[i][j][k] = 0;

  /* affectation des nouvelles valeurs par ordre croissant */
  for(i=0; i<cptVal; i++) {
    imres->mri[image[i].x][image[i].y][image[i].z] = (TYPEMRI3D)(i/ratio);
  }

  /* mise a jour du min et du max de l'image */
  //  imx_inimaxminpixel_3d_p(imres);
}

void imx_comparaison_par_rang_p(grphic3d *imref, grphic3d *imcmp, grphic3d *imres) {
  grphic3d *tmpref, *tmpcmp;

  tmpref = cr_grphic3d(imref);
  tmpcmp = cr_grphic3d(imcmp);
  
  printf("traitement de l\'image de reference...");fflush(stdout);
  imx_rechelonage01_p(imref, tmpref);
  printf(" ok\ntraitement de l\'image de comparaison...");fflush(stdout);
  imx_rechelonage01_p(imcmp, tmpcmp);
  printf(" ok\n");

  imx_sub_3d_p(tmpref, tmpcmp, imres);

  free_grphic3d(tmpref);
  free_grphic3d(tmpcmp);
}

void imx_comparaison_par_int_dis_p(grphic3d *imref, grphic3d *imcmp, grphic3d *imres) {
  const int size=256;
  const int scale=100;
  int historef[256], histocmp[256];
  int i,j,k, minref, mincmp, ampref, ampcmp;

  /* calcul des histogramme en total cumule decroissant des images a comparer */
  histogramme_cumule_dec_p(historef, imref, size, &minref, &ampref);
  histogramme_cumule_dec_p(histocmp, imcmp, size, &mincmp, &ampcmp);

  /* recopie des parametres de l'image de reference dans l'image resultat */
  imx_copie_param_3d_p(imref, imres);

  /* initialisation du rcoeff */
  imres->rcoeff = 1./(float)scale;

  /* en chaque point de l'image... */
  for(i=0;i<imref->width;i++)
    for(j=0;j<imref->height;j++)
      for(k=0;k<imref->depth;k++) {
        /* si en dehors du masque : voxel resultat positionne a zero */
        if(imref->mask && imref->mask->mri[i][j][k]==0) {
          imres->mri[i][j][k] = 0;
        }
        /* sinon voxel resultat egale a la difference du rang */
        else {
          imres->mri[i][j][k] =
            historef[((imref->mri[i][j][k] - minref) * (size-1)) / ampref] * scale / historef[0]
            - histocmp[((imcmp->mri[i][j][k] - mincmp) * (size-1)) / ampcmp] * scale / histocmp[0];
        }
      }
  /* mise a jour du min et du max de l'image */
  imx_inimaxminpixel_3d_p(imres);
}

void histogramme_cumule_dec_p(int *histo, grphic3d *image, int size, int *min, int *amplitude) {

  int i,j,k;
  int indice;
  int max;
  int nbval =0;

  /* initialisation de l'histogramme a zero */
  for(i=0;i<size;i++) histo[i]=0;

  /* initialisation du min et du max a zero */
  *min = 0;
  max = 0;

  /* calcul du min et du max de l'image */
  for(i=0;i<image->width;i++)
    for(j=0;j<image->height;j++)
      for(k=0;k<image->depth;k++)
        if(!image->mask || image->mask->mri[i][j][k]!=0) {
          if(image->mri[i][j][k]>max) max = image->mri[i][j][k];
          if(image->mri[i][j][k]<*min) *min = image->mri[i][j][k];
        }

  /* calcul de l'amplitude de l'image */
  *amplitude = max - *min;

  /* calcul de l'histogramme de l'image */
  for(i=0;i<image->width;i++)
    for(j=0;j<image->height;j++)
      for(k=0;k<image->depth;k++)
        if(!image->mask || image->mask->mri[i][j][k]!=0) {
          histo[((image->mri[i][j][k] - *min) * (size-1)) / *amplitude]++;
          nbval++;
        }

  /* calcul du total cumule decroissant */
  for(i=size-2;i>=0;i--) histo[i]+=histo[i+1];
}



/************************ unit_matched_filter_3d() ************************/
/*                 Methode de vraisemblance g��alis�                   */
/*                 'Generalized Likelihood Ratio Test'                    */
/*                                                                        */
/*            Modelisation constante de l'intensit�I=cte+bruit           */
/*                   	  Variance du bruit = 1                           */
/**************************************************************************/

void unit_matched_filter_3d(void)
{
  int im_1,im_2,im_res;
  int N;
  double sigma, fausse_alarme;
  int err=0,i,type_remplissage;
  char *quest[3];

  for(i=0;i<3;i++)
    quest[i]=CALLOC(80,char);

  strcpy(quest[0],TEXT0227);
  strcpy(quest[1],TEXT0228);
  strcpy(quest[2],"\0");

  PUT_WARN(TEXT0235);
  sigma=sqrt(2.0);

  im_1=GET_PLACE3D(TEXT0185);
  im_2=GET_PLACE3D(TEXT0186);
  im_res=GET_PLACE3D(TEXT0006);

  while((N= GET_INT(TEXT0187, 0, &err))<=0)
    PUT_WNDMSG(TEXT0189);
  CLEAR_WNDMSG();

  while((fausse_alarme= GET_DOUBLE(TEXT0221, 0, &err))<0 || fausse_alarme>1)
    PUT_WNDMSG(TEXT0222);
  CLEAR_WNDMSG();

  type_remplissage=GETV_QCM(TEXT0218,(char **)quest);

  imx_matched_filter_3d(im_1,im_2,im_res,N,fausse_alarme,type_remplissage,sigma);
  show_picture_3d(im_res);


  for(i=0;i<3;i++)
    FREE(quest[i]);

  return ;
}



/********************** imx_matched_filter_3d() ******************************/
/*               im_1,im_2: comparaison des images im_1 et im_2              */
/*                          im_res: image resultat                           */
/*                    N: taille de la fenetre d'operation                    */
/*         fausse_alarme: probabilit�tol�� des fausses alrmes             */
/*****************************************************************************/
int imx_matched_filter_3d(int im_1, int im_2, int im_res, int N, double fausse_alarme, int type_remplissage, double sigma)
{
  grphic3d *im1,*im2,*imres;

  im1=ptr_img_3d(im_1);
  im2=ptr_img_3d(im_2);
  imres=ptr_img_3d(im_res);


  imx_matched_filter_3d_p(im1,im2,imres,N,fausse_alarme,type_remplissage,sigma);

  return 0;
}

/************************ imx_matched_filter_3d_p() **************************/
/*                                                                           */
/*                         pointer based program                             */
/*****************************************************************************/

int imx_matched_filter_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres, int N, double fausse_alarme, int type_remplissage, double sigma)
{
  int width,height,depth;
  int w;
  int i,j,k,f,g,h;
  double moy1,moy2;
  double diffmoy;
  int width2,height2,depth2;
  double rcoeff1,rcoeff2;
  double rcoeff;
  int N3;
  double max=0.0,min=LONGMAX;
  double diff;
  int inf1,inf2,inf3,sup1,sup2,sup3;
  int err;
  double R;
  double *mri_temp;
  int which,status;
  double p,q,seuil,mean,sd,bound;


  width=MINI(im1->width,im2->width);
  height=MINI(im1->height,im2->height);
  depth=MINI(im1->depth,im2->depth);
  w=N/2;
  width2=width-w;
  height2=height-w;
  depth2=depth-w;
  rcoeff1=(float)im1->rcoeff;
  rcoeff2=(float)im2->rcoeff;
  N3=N*N*N;

  if (( mri_temp = (double *)malloc(sizeof(double)*width*height*depth))==NULL)
    {
      printf("failled to allocate memory\n");
      return (-1);
    }

  /* Inverse normal distribution (mean=0 , variance=1)*/
  /* Calculates threshold parameter given probability of false alarm */
  which=2;
  q=fausse_alarme/2.0;
  p=1-q;
  mean=0.0;
  sd=1.0;
  cdfnor(&which,&p,&q,&seuil,&mean,&sd,&status,&bound);

  for (i=0;i<width;i++)
    for (j=0;j<height;j++)
      for (k=0;k<depth;k++)
	mri_temp[i*height*depth+j*depth+k]=0.0;

  for(i=w;i<width2;i++)
    for(j=w;j<height2;j++)
      for (k=w;k<depth2;k++)
	{
	  moy1=0;
	  moy2=0;
	  diff=0;
	  inf1=i-w;
	  sup1=i+w;
	  inf2=j-w;
	  sup2=j+w;
	  inf3=k-w;
	  sup3=k+w;

	  diff=rcoeff2*im2->mri[i][j][k]-rcoeff1*im1->mri[i][j][k];
	  /*Modif JP 9/99 	    diff=fabs(diff);  */

	  for(f=inf1;f<=sup1;f++)
	    for(g=inf2;g<=sup2;g++)
	      for(h=inf3;h<=sup3;h++)
		{
		  moy1 += rcoeff1*(im1->mri[f][g][h]);
		  moy2 += rcoeff2*(im2->mri[f][g][h]);
		}
	  moy1=moy1/(N3);
	  moy2=moy2/(N3);
	  diffmoy=fabs(moy1-moy2);

	  R=sqrt(N3/2.0)*diffmoy/sigma;
	  if (R>seuil)
	    {
	      if (type_remplissage == 0)
		{
		  if((mri_temp[i*height2*depth2+j*depth2+k]=diff)>max) max=diff;
	 	  min=MINI(min,diff);
		}
	      if (type_remplissage == 1)
		{
	          if((mri_temp[i*height2*depth2+j*depth2+k]=log(R))>max) max=log(R);
	          min=MINI(min,log(R));
		}
	    }
	}


  /* Calcul de imres->maxpixel et imres->icomp et imres->rcoeff a partir de max*/
  err=imx_brukermax_3d((float)max,(float)min,imres);
  rcoeff=imres->rcoeff;

  for (i=w;i<width2;i++)
    for (j=w;j<height2;j++)
      for (k=w;k<depth2;k++)
        imres->mri[i][j][k]=(TYPEMRI3D)floor(mri_temp[i*height2*depth2+j*depth2+k]/rcoeff);

  imres->min_pixel=(long)floor(min);
  imres->width=width;
  imres->height=height;
  imres->depth=depth;


  free(mri_temp);
  return 0;
}



/************************* vraisemblance_3d() *****************************/
/*                 Methode de vraisemblance g��alis�                   */
/*                 'Generalized Likelihood Ratio Test'                    */
/*                                                                        */
/*            Modelisation constante de l'intensit�I=cte+bruit           */
/**************************************************************************/

void vraisemblance_3d(void)
{
  int im_1,im_2,im_res;
  int N;
  double fausse_alarme;
  int err=0,i,type_remplissage;
  char *quest[3];

  for(i=0;i<3;i++)
    quest[i]=CALLOC(80,char);

  strcpy(quest[0],TEXT0227);
  strcpy(quest[1],TEXT0228);
  strcpy(quest[2],"\0");

  im_1=GET_PLACE3D(TEXT0185);
  im_2=GET_PLACE3D(TEXT0186);
  im_res=GET_PLACE3D(TEXT0006);

  while((N= GET_INT(TEXT0187, 0, &err))<=0)
    PUT_WNDMSG(TEXT0189);
  CLEAR_WNDMSG();

  while((fausse_alarme= GET_DOUBLE(TEXT0221, 0, &err))<0 || fausse_alarme>1)
    PUT_WNDMSG(TEXT0222);
  CLEAR_WNDMSG();

  type_remplissage=GETV_QCM(TEXT0218,(char **)quest);
  imx_vraisemblance_3d(im_1,im_2,im_res,N,fausse_alarme,type_remplissage);
  show_picture_3d(im_res);


  for(i=0;i<3;i++)
    FREE(quest[i]);

  return ;
}


/********************** imx_vraisemblance_3d() *******************************/
/*               im_1,im_2: comparaison des images im_1 et im_2              */
/*                          im_res: image resultat                           */
/*                    N: taille de la fenetre d'operation                    */
/*         fausse_alarme: probabilit�tol�� des fausses alrmes             */
/*****************************************************************************/
int imx_vraisemblance_3d(int im_1, int im_2, int im_res, int N, double fausse_alarme, int type_remplissage)
{
  grphic3d *im1,*im2,*imres;

  im1=ptr_img_3d(im_1);
  im2=ptr_img_3d(im_2);
  imres=ptr_img_3d(im_res);


  imx_vraisemblance_3d_p(im1,im2,imres,N,fausse_alarme,type_remplissage);

  return 0;
}


/************************ imx_vraisemblance_3d_p() ***************************/
/*                                                                           */
/*                         pointer based program                             */
/*****************************************************************************/

int imx_vraisemblance_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres, int N, double fausse_alarme, int type_remplissage)
{
  int width,height,depth;
  int w;
  int i,j,k,f,g,h;
  double moy1,moy2;
  double diffmoy;
  double sigma1,sigma2;
  double var1,var2,covar;
  int width2,height2,depth2;
  double rcoeff1,rcoeff2;
  double rcoeff;
  int N3;
  double max=0.0,min=LONGMAX;
  double diff;
  double I1,I2;
  int inf1,inf2,inf3,sup1,sup2,sup3;
  int err;
  double R;
  double *mri_temp;
  int which,status;
  double p,q,seuil,degree_of_freedom_n,degree_of_freedom_d,bound;


  width=MINI(im1->width,im2->width);
  height=MINI(im1->height,im2->height);
  depth=MINI(im1->depth,im2->depth);
  w=N/2;
  width2=width-w;
  height2=height-w;
  depth2=depth-w;
  rcoeff1=(float)im1->rcoeff;
  rcoeff2=(float)im2->rcoeff;
  N3=N*N*N;

  if (( mri_temp = (double *)malloc(sizeof(double)*width*height*depth))==NULL)
    {
      printf("failled to allocate memory\n");
      return (-1);
    }

  /* Inverse Fisher Distribution (1,2n-2) degrees of freedom */
  /* Calculates threshold parameter given probability of false alarm */
  which=2;
  q=fausse_alarme;
  p=1-q;
  degree_of_freedom_n=1;
  degree_of_freedom_d=N3+N3-2;
  cdff(&which,&p,&q,&seuil,&degree_of_freedom_n,&degree_of_freedom_d,&status,&bound);

  for (i=0;i<width;i++)
    for (j=0;j<height;j++)
      for (k=0;k<depth;k++)
	{
	  imres->mri[i][j][k]=0;
	  mri_temp[i*height*depth+j*depth+k]=0.0;
	}

  for(i=w;i<width2;i++)
    for(j=w;j<height2;j++)
      for (k=w;k<depth2;k++)
	{
	  moy1=0;
	  moy2=0;
	  sigma2=0;
	  sigma1=0;
	  var1=0;
	  var2=0;
	  diff=0;
	  inf1=i-w;
	  sup1=i+w;
	  inf2=j-w;
	  sup2=j+w;
	  inf3=k-w;
	  sup3=k+w;

	  diff=rcoeff2*im2->mri[i][j][k]-rcoeff1*im1->mri[i][j][k];
	  /*Modif JP 9/99 diff=fabs(diff);   */

	  for(f=inf1;f<=sup1;f++)
	    for(g=inf2;g<=sup2;g++)
	      for(h=inf3;h<=sup3;h++)
		{
		  moy1 += rcoeff1*(im1->mri[f][g][h]);
		  moy2 += rcoeff2*(im2->mri[f][g][h]);
		}
	  moy1=moy1/(N3);
	  moy2=moy2/(N3);
	  diffmoy=fabs(moy1-moy2);

	  for(f=inf1;f<=sup1;f++)
	    for(g=inf2;g<=sup2;g++)
	      for(h=inf3;h<=sup3;h++)
		{
		  I1=rcoeff1*(im1->mri[f][g][h])-moy1;
		  I2=rcoeff2*(im2->mri[f][g][h])-moy2;
		  var1 += I1*I1;
		  var2 += I2*I2;
		}

	  var1=var1/N3;
	  var2=var2/N3;

	  var1=MAXI(var1,1.00e-9);
	  var2=MAXI(var2,1.00e-9);

	  sigma1=sqrt(var1);
	  sigma2=sqrt(var2);
	  covar=sigma1*sigma2;


	  R = (var1+var2)/2.0;
	  R = R + (diffmoy*diffmoy/4.0);
	  R = R - (covar);
	  R = R/(covar);
	  R = (N3-1)*R;

	  if (R>seuil)
	    {
	      if (type_remplissage == 0)
		{
		  if((mri_temp[i*height2*depth2+j*depth2+k]=diff)>max) max=diff;
	 	  min=MINI(min,diff);
		}
	      if (type_remplissage == 1)
		{
	          if((mri_temp[i*height2*depth2+j*depth2+k]=log(R))>max) max=log(R);
	          min=MINI(min,log(R));
		}
	    }
	}


  /* Calcul de imres->maxpixel et imres->icomp et imres->rcoeff a partir de max*/
  err=imx_brukermax_3d((float)max,(float)min,imres);
  rcoeff=imres->rcoeff;

  for (i=w;i<width2;i++)
    for (j=w;j<height2;j++)
      for (k=w;k<depth2;k++)
        imres->mri[i][j][k]=(TYPEMRI3D)floor(mri_temp[i*height2*depth2+j*depth2+k]/rcoeff);

  imres->min_pixel=0;
  imres->width=width;
  imres->height=height;
  imres->depth=depth;



  free(mri_temp);
  return 0;
}



/**************************** kolmogorov_3d() ********************************/
/*                      Kolmogorov test for change detection                 */
/*                                                                           */
/*      Given two images for comparison, and probability of false alarm,     */
/*      this routine calculates the Kolmogorov-Smirnov statistic and then    */
/*      decides for change detection.                                        */
/*                                                                           */
/*****************************************************************************/

void kolmogorov_3d(void)
{
  int im_1,im_2,im_res;
  int N;
  int err=0;
  int type_remplissage;
  float fausse_alarme;
  char *quest[] = {TEXT0227, TEXT0228, "" };

  im_1=GET_PLACE3D(TEXT0185);
  im_2=GET_PLACE3D(TEXT0186);
  im_res=GET_PLACE3D(TEXT0006);

  while((N= GET_INT(TEXT0187, 0, &err))<=0)
    PUT_WNDMSG(TEXT0189);
  CLEAR_WNDMSG();

  while((fausse_alarme= GET_FLOAT(TEXT0221, 0, &err))<0 || fausse_alarme>1)
    PUT_WNDMSG(TEXT0222);
  CLEAR_WNDMSG();

  type_remplissage=GETV_QCM(TEXT0218,(char **)quest);
  imx_kolmogorov_3d(im_1,im_2,im_res,N,fausse_alarme,type_remplissage);
  show_picture_3d(im_res);

  return ;
}


/************************* imx_kolmogorov_3d() *******************************/
/*                          imx style program                                */
/*****************************************************************************/

int imx_kolmogorov_3d(int im_1, int im_2, int im_res, int N,
		      float fausse_alarme, int type_remplissage)
{
  grphic3d *im1,*im2,*imres;

  im1=ptr_img_3d(im_1);
  im2=ptr_img_3d(im_2);
  imres=ptr_img_3d(im_res);

  if (fausse_alarme !=0.001 && fausse_alarme !=0.01 && fausse_alarme != 0.05) fausse_alarme = (float)0.001;
  imx_kolmogorov_3d_p(im1,im2,imres,N,fausse_alarme,type_remplissage);

  return 0;
}


/************************* imx_kolmogorov_3d_p() *****************************/
/*                         pointer based program                             */
/*****************************************************************************/

int imx_kolmogorov_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres, int N, float fausse_alarme, int type_remplissage)
{
  float kolmogorov_proba(float lambda);
  float *mri_temp;
  int width,height,depth;
  int w;
  int i,j,k,f,g,h;
  int j1,j2;
  int N2,N3;
  int width2,height2,depth2;
  int inf1,inf2,inf3,sup1,sup2,sup3;
  int err;
  float rcoeff,rcoeff1,rcoeff2;
  float diff, max=0.0,min=LONGMAX;
  float *data1, *data2, *data1_head, *data2_head;
  float d1,d2,d_temp,d,cumul1,cumul2;
  float seuil, confidence;
  TYPEMRI3D ***image1, ***image2;

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
  N2=N*N;

  if (( mri_temp = (float *)malloc(sizeof(double)*width*height*depth))==NULL)
    {
      PUT_ERRFAT("failled to allocate memory\n");
      return (-1);
    }

  /* Distribution de Kolmogorov */

  if (fausse_alarme==1.0e-8) seuil=(float)3.0914;
  if (fausse_alarme==0.001) seuil=(float)1.950;
  if (fausse_alarme==0.01) seuil=(float)1.628;
  if (fausse_alarme==0.05) seuil=(float)1.359;


  data1 = (float *)calloc((size_t) N3,(size_t) sizeof(float));
  data2 = (float *)calloc((size_t) N3,(size_t) sizeof(float));
  data1_head = data1;
  data2_head = data2;

  bzero(mri_temp, sizeof(float) * height * width * depth);
/*   bzero(imres->mri[0][0], sizeof(TYPEMRI3D) * imres->mri_width */
/* 	* imres->mri_depth * imres->mri_height); */

  image1 = im1->mri;
  image2 = im2->mri;
  {
    START_TIMER;
    for(i=w;i<width2;i++)
      for(j=w;j<height2;j++)
	for (k=w;k<depth2;k++)
	  {
	    inf1=i-w;
	    sup1=i+w;
	    inf2=j-w;
	    sup2=j+w;
	    inf3=k-w;
	    sup3=k+w;

	    diff = rcoeff1 * image1[i][j][k] - rcoeff1 * image1[i][j][k];
	    if((diff=(float)fabs(diff))>max) max = diff;

	    for(f=inf1;f<=sup1;f++)
	      for(g=inf2;g<=sup2;g++)
		for(h=inf3;h<=sup3;h++)
		  {
		    *(data1++) = rcoeff1 * image1[f][g][h];
		    *(data2++) = rcoeff2 * image2[f][g][h];
		  }

	    data1=data1_head;
	    data2=data2_head;

	    shell_sort(N3,data1);
	    shell_sort(N3,data2);

	    j1=1;
	    j2=1;
	    cumul1=0.0;
	    cumul2=0.0;
	    d=0.0;
	    while (j1<=N3 && j2<=N3)
	      {
		if ((d1=data1[j1-1]) <= (d2=data2[j2-1])) cumul1=j1++/(float)N3;
		if (d2 <= d1) cumul2=j2++/(float)N3;
		if ((d_temp=(float)fabs(cumul1-cumul2)) > d) d=d_temp;
	      }

	    d=(float)(d*sqrt(N3/2.0));
	    confidence=kolmogorov_proba(d);

	    if(confidence <= fausse_alarme)
	      {
		if (type_remplissage == 0)
		  {
		    if((mri_temp[i*height2*depth2+j*depth2+k]=diff)>max)
		      max=diff;
		    min=MINI(min,diff);
		  }
		if (type_remplissage == 1)
		  {
		    if((mri_temp[i*height2*depth2+j*depth2+k]=(float)log(d))>max)
		      max=(float)log(d);
		    min=(float)MINI(min,log(d));
		  }
	      }
	    else
	      mri_temp[i*height2*depth2+j*depth2+k] = 0.0;
	  }

    /* Calcul de imres->maxpixel et imres->icomp et imres->rcoeff a partir de max*/
    err=imx_brukermax_3d(max,min,imres);
    rcoeff=imres->rcoeff;

    for (i=w;i<width2;i++)
      for (j=w;j<height2;j++)
	for (k=w;k<depth2;k++)
	  imres->mri[i][j][k]=(TYPEMRI3D)floor(mri_temp[i*height2*depth2+j*depth2+k]/rcoeff);
    END_TIMER;
  }

  free(mri_temp);
  free(data1);
  free(data2);

  return 0;

}



/************************ shell_sort() *******************************
 *       Sorts an array data[1..n] into ascending numerical order
 *                       by Shell's method
 *       data is replaced on output by its sorted rearrangement
 **********************************************************************/

void shell_sort(int n, float *data)
{
  int i,j,inc;
  float v;
  float *data_temp;

  data_temp=(float *)calloc((size_t) (n+1),(size_t) sizeof(float));

  for(i=0;i<n;i++)
    data_temp[i+1] = data[i];

  inc=1;

  do {
    inc *= 3;
    inc++;
  } while (inc <= n);

  do {
    inc /= 3;
    for (i=inc+1;i<=n;i++)
      {
       	v=data_temp[i];
	j=i;
	while (data_temp[j-inc] > v)
	  {
	    data_temp[j] = data_temp[j-inc];
	    j -= inc;
	    if (j<=inc) break;
	  }
	data_temp[j]=v;
      }
  } while (inc > 1);

  for(i=0;i<n;i++)
    data[i] = data_temp[i+1];

  free(data_temp);
}

/*************** kolmogorov_proba(float lambda) *****************************/
/*               Kolmogorov probability function                            */
/*                                                                          */
/*      Given  the Kolmogorov statistic lambda, this routine                */
/*      returns corresponding Kolmogorov probability number                 */
/*                                                                          */
/****************************************************************************/

float kolmogorov_proba(float lambda)
{
  int j;
  float term1,term2,sum=0.0,termbf=0.0;
  float sign=1.0,eps1=(float)1.0e-3, eps2=(float)1.0e-8;

  term1 =(float)( -2.0*lambda*lambda);
  for (j=1;j<=100;j++)
    {
      term2 = (float)(2.0*sign*exp(term1*j*j));
      sum += term2;
      if ( fabs(term2) <= eps1*termbf || fabs(term2) <= eps2*sum) return sum;
      sign = -sign;
      termbf = (float)fabs(term2);
    }
  return 1.0;
}


/************************    sigma_median_3d()   *****************************/
/*								 	     */
/*****************************************************************************/
void sigma_median_3d(void)
{
  int im_ref;
  float median;

  im_ref=GET_PLACE3D(TEXT0185);

  imx_sigma_median_3d(im_ref,&median);

  PUTF("Noise Std. Div. (median value) =",median," ");

}

/*****************************************************************************/
void imx_sigma_median_3d(int im_ref, float *median)
{
  grphic3d *imref;

  imref=ptr_img_3d(im_ref);
  imx_sigma_median_3d_p(imref,median);
}


/*****************************************************************************/
void imx_sigma_median_3d_p(grphic3d *imref, float *median)
{
  int i,j,k;
  int start_x,start_y,start_z,size,edge,edge_2;
  int width,height,depth;
  float *sigma,*sigma_head;

  if (( sigma = (float *)malloc(sizeof(float)*64))==NULL)
    {
      printf("failled to allocate memory\n");
    }

  edge=0;
  size=5;
  edge_2=2*size+edge;
  width=imref->width;
  height=imref->height;
  depth=imref->depth;
  sigma_head=sigma;

  /*** corner 000 ***/
  start_x=edge;
  for (i=0;i<2;i++)
    {
      start_y=edge;
      for (j=0;j<2;j++)
	{
	  start_z=edge;
	  for (k=0;k<2;k++)
	    {
	      *sigma++=imx_sigma_from_VOI_p(imref,start_x,start_y,start_z,size);
	      start_z+=size;
	    }
	  start_y+=size;
	}
      start_x+=size;
    }


  /*** corner 001 ***/
  start_x=edge;
  for (i=0;i<2;i++)
    {
      start_y=edge;
      for (j=0;j<2;j++)
	{
	  start_z=depth-edge_2;
	  for (k=0;k<2;k++)
	    {
	      *sigma++=imx_sigma_from_VOI_p(imref,start_x,start_y,start_z,size);
	      start_z+=size;
	    }
	  start_y+=size;
	}
      start_x+=size;
    }


  /*** corner 010 ***/
  start_x=edge;
  for (i=0;i<2;i++)
    {
      start_y=height-edge_2;
      for (j=0;j<2;j++)
	{
	  start_z=edge;
	  for (k=0;k<2;k++)
	    {
	      *sigma++=imx_sigma_from_VOI_p(imref,start_x,start_y,start_z,size);
	      start_z+=size;
	    }
	  start_y+=size;
	}
      start_x+=size;
    }


  /*** corner 011 ***/
  start_x=edge;
  for (i=0;i<2;i++)
    {
      start_y=height-edge_2;
      for (j=0;j<2;j++)
	{
	  start_z=depth-edge_2;
	  for (k=0;k<2;k++)
	    {
	      *sigma++=imx_sigma_from_VOI_p(imref,start_x,start_y,start_z,size);
	      start_z+=size;
	    }
	  start_y+=size;
	}
      start_x+=size;
    }


  /*** corner 100 ***/
  start_x=width-edge_2;
  for (i=0;i<2;i++)
    {
      start_y=edge;
      for (j=0;j<2;j++)
	{
	  start_z=edge;
	  for (k=0;k<2;k++)
	    {
	      *sigma++=imx_sigma_from_VOI_p(imref,start_x,start_y,start_z,size);
	      start_z+=size;
	    }
	  start_y+=size;
	}
      start_x+=size;
    }


  /*** corner 101 ***/
  start_x=width-edge_2;
  for (i=0;i<2;i++)
    {
      start_y=edge;
      for (j=0;j<2;j++)
	{
	  start_z=depth-edge_2;
	  for (k=0;k<2;k++)
	    {
	      *sigma++=imx_sigma_from_VOI_p(imref,start_x,start_y,start_z,size);
	      start_z+=size;
	    }
	  start_y+=size;
	}
      start_x+=size;
    }


  /*** corner 110 ***/
  start_x=width-edge_2;
  for (i=0;i<2;i++)
    {
      start_y=height-edge_2;
      for (j=0;j<2;j++)
	{
	  start_z=edge;
	  for (k=0;k<2;k++)
	    {
	      *sigma++=imx_sigma_from_VOI_p(imref,start_x,start_y,start_z,size);
	      start_z+=size;
	    }
	  start_y+=size;
	}
      start_x+=size;
    }


  /*** corner 111 ***/
  start_x=width-edge_2;
  for (i=0;i<2;i++)
    {
      start_y=height-edge_2;
      for (j=0;j<2;j++)
	{
	  start_z=depth-edge_2;
	  for (k=0;k<2;k++)
	    {
	      *sigma++=imx_sigma_from_VOI_p(imref,start_x,start_y,start_z,size);
	      start_z+=size;
	    }
	  start_y+=size;
	}
      start_x+=size;
    }


  /*** sorting data ***/
  sigma=sigma_head;
  shell_sort(64,sigma);


  /*** median ***/
  *median=(float)((sigma[31]+sigma[32])/2.0);

  free(sigma);
}

/*****************************************************************/
float imx_sigma_from_VOI_p(grphic3d *imref, int start_x, int start_y, int start_z, int size)
{
  int i,j,k,m;
  float moy,sigma;

  m=0;
  moy = 0.0;
  for (i=start_x;i<start_x+size;i++)
    for (j=start_y;j<start_y+size;j++)
      for (k=start_z;k<start_z+size;k++)
	{
	  m++;
	  moy += ((imref->rcoeff)*(imref->mri[i][j][k]));
	}

  moy=moy/(m);
  sigma=(float)(moy*sqrt(2.0/PI));

  return (sigma);
}


/************************  standardize_noise_3d() ****************************/
/*								 	     */
/*****************************************************************************/
void standardize_noise_3d(void)
{
  int im_ref,im_res;
  float median;

  im_ref=GET_PLACE3D(TEXT0185);
  im_res=GET_PLACE3D(TEXT0006);

  imx_sigma_median_3d(im_ref,&median);
  PUTF("Noise Std. Div. (median value) =",median," ");

  if (median!=0)
    imx_standardize_noise_3d(im_ref,im_res,median);
  else
    PUT_ERROR("function aborted: zero noise variance!!");

  show_picture_3d(im_res);

}

/*****************************************************************************/
int imx_standardize_noise_3d(int im_ref, int im_res, float sigma)
{
  grphic3d *imref,*imres;

  imref=ptr_img_3d(im_ref);
  imres=ptr_img_3d(im_res);
  imx_standardize_noise_3d_p(imref,imres,sigma);

  return 0;
}


/*****************************************************************************/
void imx_standardize_noise_3d_p(grphic3d *imref, grphic3d *imres, float sigma)
{
  int width,height,depth,i,j,k,err;
  float max,min;
  float rcoeff1,rcoeff2;

  width=imref->width;
  height=imref->height;
  depth=imref->depth;

  imx_copie_param_3d_p(imref,imres);

  /* Calcul de imres->maxpixel, imres->icomp et imres->rcoeff */
  rcoeff1=imref->rcoeff;
  max=(rcoeff1*imref->max_pixel)/sigma;
  min=(rcoeff1*imref->min_pixel)/sigma;
  err=imx_brukermax_3d(max,min,imres);
  rcoeff2=imres->rcoeff;

  for (i=0;i<width;i++)
    for (j=0;j<height;j++)
      for (k=0;k<depth;k++)
	imres->mri[i][j][k]=(TYPEMRI3D)floor(rcoeff1*imref->mri[i][j][k]/(sigma*rcoeff2));
}


/********************** roc_normal() *****************************/
/*								 */
/*****************************************************************/
void roc_normal(void)
{
  int which,status,err=0;
  double p1,q1,seuil,mean1,sd1,bound;
  double p2,q2,mean2,sd2;

  while((mean2= GET_FLOAT("SNR", 0, &err))<=0)
    PUT_WNDMSG("use a positive SNR");
  CLEAR_WNDMSG();

  while((q1= GET_DOUBLE(TEXT0221, 0, &err))<0 || q1>1)
    PUT_WNDMSG(TEXT0222);
  CLEAR_WNDMSG();

  PUTF("\nSNR=",(float)mean2," ");

  which=2;
  q1/=2.0;
  p1=1-q1;
  mean1=0.0;
  sd1=1.0;
  cdfnor(&which,&p1,&q1,&seuil,&mean1,&sd1,&status,&bound);

  which=1;
  sd2=1.0;
  cdfnor(&which,&p2,&q2,&seuil,&mean2,&sd2,&status,&bound);

  PUTF("PFA=",(float)(2*q1),"_");
  PUTF("\tPD=",(float)q2,"_");
  PUTF("\tSeuil=",(float)seuil,"");
}



/******************************************************************************
 ** -- Hist_Eq() ---------------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July 22nd, 1999
 **
 ******************************************************************************/
void Hist_Eq_3d(void)
{
  grphic3d *ia,*im1,*imres;
  float *x,*hist,xmin,xmax;
  int i,j,k,wdth,hght,dpth,im_1,im_res;

  im_1=GET_PLACE3D( "Histogram equalization from");
  im_res=GET_PLACE3D( "to :");
  im1=ptr_img_3d(im_1);
  imres=ptr_img_3d(im_res);
  if ((ia=(grphic3d *)cr_grphic3d(NULL))==NULL) {
    PUT_ERRFAT(
	       "Can't allocate temporary image in imx_user()");
    _imagix_error=IMX_NOT_ENOUGH_CORE;
    return;
  }
  x= CALLOC(257,float);
  hist= CALLOC(257,float);
  wdth=im1->width;
  hght=im1->height;
  dpth=im1->depth;
  xmin=0.;
  xmax=im1->rcoeff*im1->max_pixel;
  for (i=0;i<wdth;i++)
    for (j=0;j<hght;j++)
      for (k=0;k<dpth;k++)
	ia->mri[i][j][k] = (im1->mri[i][j][k] > 0) ? 1 : 0 ;

  for (i=0;i<256;i++) hist[i]=0.0;
  imx_copie_param_3d_p(im1,imres);
  imx_histogram_3d_p(im1,ia,xmin,xmax,(long int)256,x,hist);
  imx_histogram_equalization_3d_p(im1,imres,ia,hist,(long int)256);
  for (i=0;i<wdth;i++)
    for (j=0;j<hght;j++)
      for (k=0;k<dpth;k++)
	if ( imres->mri[i][j][k] < 0 ) imres->mri[i][j][k]=0;
  imx_inimaxminpixel_3d_p(imres); /* Positionne max,min,   ...*/
  show_picture_3d(im_res);
  free_grphic3d(ia);
  FREE(x);
  FREE(hist);
}
/*-----------------  image_pixel_3d	--------------------------*/
/*
  return TRUE if point (x,y,z) is within the bounds of image im1resc
*/
int image_pixel_3d(int x, int y, int z, grphic3d *im1resc)
{
  if ( (x>=0) && (x<(int)im1resc->width) && (y>=0) && (y<(int)im1resc->height) && (z>=0)  && (z<(int)im1resc->depth) )
    return(TRUE);
  else
    return(FALSE);
}

/*******************************************************************************
 **        convoluer_filtre_1D(ligneIn, ligneOut, tailleLigne, filtre, tailleFiltre)
 */
/*!       filtrage lineaire 1D, les effets de bord sont geres en repliant le veteur en entree sur lui meme
**
**        \param ligneIn, ligneOut, tailleLigne: vecteurs en entree, en sortie et taille du vecteur
**        \param filtre, tailleFiltre: filtre de convolution 1D et taille du filtre
**        \retval : 0 en cas de succes
**
*******************************************************************************/
int convoluer_filtre_1D(double *ligneIn, double *ligneOut, int tailleLigne, double *filtre, int tailleFiltre)
{
  int i, j, p;
  double sum;

  for (i=tailleFiltre/2; i<(tailleLigne-tailleFiltre/2); i++)
    {
      sum=0.0;
      for (p=-tailleFiltre/2; p<=tailleFiltre/2; p++) sum+=ligneIn[i-p]*filtre[p+tailleFiltre/2];
      ligneOut[i]=sum;
    }

  //on remplit les bords en repliant le ligne sur elle meme
  for (i=0;i<tailleFiltre/2;i++)
    {
      sum=0.0;
      for (p=-tailleFiltre/2; p<=tailleFiltre/2; p++)
	{
	  j=i-p; if (j<0) j+=tailleLigne;
	  sum+=ligneIn[j]*filtre[p+tailleFiltre/2];
	}
      ligneOut[i]=sum;
    }

  for (i=(tailleLigne-tailleFiltre/2);i<tailleLigne;i++)
    {
      sum=0.0;
      for (p=-tailleFiltre/2; p<tailleFiltre/2; p++)
	{
	  j=i-p; if (j>=tailleLigne) j%=tailleLigne;
	  sum+=ligneIn[j]*filtre[p+tailleFiltre/2];
	}
      ligneOut[i]=sum;
    }

  return 0;
}

/*******************************************************************************
 **        convoluer_filtre_1D_sans_bord(ligneIn, ligneOut, tailleLigne, filtre, tailleFiltre)
 */
/*!       filtrage lineaire 1D, les bords sont ignores
**
**        \param ligneIn, ligneOut, tailleLigne: vecteurs en entree, en sortie et taille du vecteur
**        \param filtre, tailleFiltre: filtre de convolution 1D et taille du filtre
**        \retval : 0 en cas de succes
**
*******************************************************************************/
int convoluer_filtre_1D_sans_bord(double *ligneIn, double *ligneOut, int tailleLigne, double *filtre, int tailleFiltre)
{
  int i,p;
  double sum;

  for (i=tailleFiltre/2; i<(tailleLigne-tailleFiltre/2); i++)
    {
      sum=0.0;
      for (p=-tailleFiltre/2; p<=tailleFiltre/2; p++) sum+=ligneIn[i-p]*filtre[p+tailleFiltre/2];
      ligneOut[i]=sum;
    }

  return 0;
}

/*******************************************************************************
 **        imx_filtre_gaussien_anisotrope_3d_p(imdeb, imres, wnd, sigmax, sigmay, sigmaz)
 */
/*!       filtrage gaussien d'une image avec valaur d'ecart type differentes selon les directions
**
**        \param imdeb, imres: les images a filtrer et resultat
**        \param wnd: taille de fenetre du filtrage
**        \param sigmax, sigmay, sigmaz: ecarts types selon les directions
**        \retval : 0 en cas de succes
**
*******************************************************************************/
int imx_filtre_gaussien_anisotrope_3d_p(grphic3d *imdeb, grphic3d *imres, unsigned int wnd, double sigmax, double sigmay, double sigmaz)
{
  unsigned int i,j,k;
  TDimension width,height,depth;
  double coeff;
  double normx=0.0, normy=0.0, normz=0.0;
  double *filterx, *filtery, *filterz;
  double *ligneIn=NULL, *ligneOut=NULL;
  float ***tmp_res;

  width=imdeb->width; height=imdeb->height; depth=imdeb->depth;

  if ((imres->width!=width)||(imres->height!=height)||(imres->depth!=depth))
    {
      fprintf (stderr,"l'image res n'a pas la meme taille que l'image a filtrer dans imx_filtre_gaussien_anisotrope_3d_p\n");
      return 1;
    }

  if (2*(wnd/2)==wnd) { fprintf(stderr,"la taille de fenetre doit etre impaire dans imx_filtre_gaussien_anisotrope_3d_p\n"); return 2; }

  //allocations memoire
  filterx=CALLOC(wnd, double); filtery=CALLOC(wnd, double); filterz=CALLOC(wnd, double);
  // image temporaire en double
  tmp_res=alloc_matrix_3d(width, height, depth);
  if (tmp_res == (float ***)NULL) { aff_log("Unable to allocate memory"); return 3; }

  //remplissage des filtres
  coeff=sqrt(2.0*PI);
  if (sigmax!=0)
    {
      for(i=0;i<wnd;i++) { filterx[i]=(1/(coeff*sigmax))*exp(-((double)((i-wnd/2)*(i-wnd/2)))/(2.0*sigmax*sigmax)); normx+=filterx[i]; }
      for(i=0;i<wnd;i++) { filterx[i]/=normx; }
    }
  else { filterx[wnd/2]=1.0; }
  if (sigmay!=0)
    {
      for(i=0;i<wnd;i++) { filtery[i]=(1/(coeff*sigmay))*exp(-((double)((i-wnd/2)*(i-wnd/2)))/(2.0*sigmay*sigmay)); normy+=filtery[i]; }
      for(i=0;i<wnd;i++) { filtery[i]/=normy; }
    }
  else { filtery[wnd/2]=1.0; }
  if (sigmaz!=0)
    {
      for(i=0;i<wnd;i++) { filterz[i]=(1/(coeff*sigmaz))*exp(-((double)((i-wnd/2)*(i-wnd/2)))/(2.0*sigmaz*sigmaz)); normz+=filterz[i]; }
      for(i=0;i<wnd;i++) { filterz[i]/=normz; }
    }
  else { filterz[wnd/2]=1.0; }

  //init
  for(i=0;i<width;i++)
    for(j=0;j<height;j++)
      for(k=0;k<depth;k++)
	tmp_res[i][j][k]=(float)imdeb->mri[i][j][k];

  //filtrage selon x
  ligneIn=CALLOC(width, double);
  ligneOut=CALLOC(width, double);
  for(k=0;k<depth;k++)
    for(j=0;j<height;j++)
      {
	for (i=0;i<width;i++) ligneIn[i]=tmp_res[i][j][k];
	convoluer_filtre_1D(ligneIn, ligneOut, width, filterx, wnd);
	for (i=0;i<width;i++) tmp_res[i][j][k]=(float)ligneOut[i];
      }
  FREE(ligneIn); FREE(ligneOut);

  //filtrage selon y
  ligneIn=CALLOC(height, double);
  ligneOut=CALLOC(height, double);
  for(k=wnd/2;k<depth-wnd/2;k++)
    for(i=wnd/2;i<width-wnd/2;i++)
      {
	for (j=0;j<height;j++) ligneIn[j]=tmp_res[i][j][k];
	convoluer_filtre_1D(ligneIn, ligneOut, height, filtery, wnd);
	for (j=0;j<height;j++) tmp_res[i][j][k]=(float)ligneOut[j];
      }
  FREE(ligneIn); FREE(ligneOut);

  //filtrage selon z
  ligneIn=CALLOC(depth, double);
  ligneOut=CALLOC(depth, double);
  for(i=wnd/2;i<width-wnd/2;i++)
    for(j=wnd/2;j<height-wnd/2;j++)
      {
	for (k=0;k<depth;k++) ligneIn[k]=tmp_res[i][j][k];
	convoluer_filtre_1D(ligneIn, ligneOut, depth, filterz, wnd);
	for (k=0;k<depth;k++) tmp_res[i][j][k]=(float)ligneOut[k];
      }
  FREE(ligneIn); FREE(ligneOut);

  //remplissage de l'image finale
  imx_copie_param_3d_p(imdeb, imres);

  for(i=0;i<width;i++)
    for(j=0;j<height;j++)
      for(k=0;k<depth;k++)
	imres->mri[i][j][k]=(int)floor(tmp_res[i][j][k]+0.5);

  imx_inimaxminpixel_3d_p(imres);

  //liberations memoire
  FREE(filterx); FREE(filtery); FREE(filterz);
  free_matrix_3d(tmp_res);

  return 0;
}


/******************************************************************************
 **    filtre_moyenneur_3d() :
 **
 */
/*! filtrage moyenneur
************************************************************************/
void filtre_moyenneur_3d()
{
  int im_src, im_dst, taille_filtre=0;

  //choix des images
  im_src=GET_PLACE3D("Image Source");
  im_dst=GET_PLACE3D("Image Destination");

  //choix de la taille du filtre
  taille_filtre=imx_query_filter_dimension();

  imx_filtre_moyenneur_3d(im_src, im_dst, taille_filtre);

  show_picture_3d(im_dst);
}

/******************************************************************************
 **
 **   imx_filtre_moyenneur_3d()
 */
/*!   filtrage moyenneur
**	 \param im_src : numero de l'image source
**   \param im_dst : numero de l'image resultat
**   \param taille_filtre : taille de la fenetre
**
*************************************************************************/
int imx_filtre_moyenneur_3d(int im_src, int im_dst, int taille_filtre)
{
  grphic3d *imsrc=NULL, *imdst=NULL;
  int res=0;

  imsrc=ptr_img_3d(im_src);
  imdst=ptr_img_3d(im_dst);

  res=imx_filtre_moyenneur_3d_p(imsrc, imdst, taille_filtre);

  return res;
}

/******************************************************************************
 **
 **  imx_median_3d_p()
 */
/*!  \ingroup ToolBoxImage3d
**
**  filtrage moyenneur
**  \param imsrc :image source
**  \param imdst : image resultat (E/S)
**  \param taille_filtre : taille de la fenetre
**
****************************************************************************/
int imx_filtre_moyenneur_3d_p(grphic3d *imsrc, grphic3d *imdst, int taille_filtre)
{
 int i,j,k, l,m,n;
 int min_i, max_i, min_j, max_j, min_k, max_k;
 int err=0;
 long sum=0;
 int aire_sum=0;
 int width, height, depth;
 TYPEMRI3D ***imsrcMRI=NULL, ***imdstMRI=NULL;
 grphic3d *imtmp=NULL;

 if ((!imsrc)||(!imdst))
  { fprintf(stderr,"les deux images doivent etre non nulles dans imx_filtre_moyenneur_3d_p\n"); return 1; }

 if (taille_filtre%2!=1)
 { fprintf(stderr,"la taille du filtre doit etre impaire dans imx_filtre_moyenneur_3d_p\n"); return 2; }

 //pour le cas ou on filtrerait l'image sur elle meme
  if (imdst==imsrc)
  {
    imtmp=cr_grphic3d(imdst);
    if (!imtmp) { fprintf(stderr, "erreur allocation memoire dans imx_median_3d_p\n"); err=1; goto end_func; }
    imdstMRI=imtmp->mri;
  }
  //la copie de l'image permet aussi de s'assurer que le buffer de imres a la bonne taille
  else
  {
   err=imx_copie_3d_p(imsrc, imdst); if (err) goto end_func;
   imdstMRI=imdst->mri;
  }

// //recopie des parametres dans l'image destination
// imx_copie_param_3d_p(imsrc, imdst);
 //la copie de l'image permet aussi de s'assurer que le buffer de imres a la bonne taille
 err=imx_copie_3d_p(imsrc, imdst); if (err) return err;

 width=(int)imsrc->width; height=(int)imsrc->height; depth=(int)imsrc->depth;

 //filtrage
 taille_filtre=taille_filtre/2; imsrcMRI=imsrc->mri;
 for (i=0;i<width;i++)
 {
  min_i=MAXI(0, i-taille_filtre); max_i=MINI(width-1, i+taille_filtre);
  for (j=0;j<height;j++)
  {
   min_j=MAXI(0, j-taille_filtre); max_j=MINI(height-1, j+taille_filtre);
   for (k=0;k<depth;k++)
   {
    min_k=MAXI(0, k-taille_filtre); max_k=MINI(depth-1, k+taille_filtre);

    sum=0; aire_sum=0;
    for (l=min_i;l<=max_i;l++)
     for (m=min_j;m<=max_j;m++)
      for (n=min_k;n<=max_k;n++)
       { sum +=imsrcMRI[l][m][n]; aire_sum++; }

    imdstMRI[i][j][k]=(TYPEMRI3D)floor((double)sum/(double)aire_sum + 0.5);
   }
  }
 }

 if (imtmp) imx_copie_3d_p(imtmp, imdst);

 imx_inimaxminpixel_3d_p(imdst);

end_func:

 if (imtmp) free_grphic3d(imtmp);

 return err;
}

/******************************************************************************
 **    norme_sqr_quad_gradient() :
 **
 */
/*! calcule et affiche la norme quadratique carree du gradient d'une image
************************************************************************/
void norme_sqr_quad_gradient()
{
  int im_src, im_dst;

  //choix des images
  im_src=GET_PLACE3D("Image Source");
  im_dst=GET_PLACE3D("Image Destination");

  imx_norme_sqr_quad_gradient_3d(im_src, im_dst);

  show_picture_3d(im_dst);
}

/******************************************************************************
 **    imx_norme_sqr_quad_gradient_3d() :
 **
 */
/*! calcule et affiche la norme quadratique carree du gradient d'une image
**
**  \param im_src :image source
**  \param im_dst : image resultat (E/S)
**
************************************************************************/
int imx_norme_sqr_quad_gradient_3d(int im_src, int im_dst)
{
 grphic3d *imsrc=NULL, *imdst=NULL;
 double *d_image=NULL, *d_gradient=NULL;
 int res=0;
 TDimension width, height, depth;
 unsigned int i,j,k;
 double max=0.0;
 double rcoeff;
 double iVal;

 imsrc=ptr_img_3d(im_src);
 imdst=ptr_img_3d(im_dst);

 width=imsrc->width; height=imsrc->height; depth=imsrc->depth;

 //allocations memoire
 d_image=CALLOC(width*height*depth, double);
 d_gradient=CALLOC(width*height*depth, double);

 //on remplit l'image double
 rcoeff=imsrc->rcoeff;
 for (i=0;i<width;i++)
  for (j=0;j<height;j++)
   for (k=0;k<depth;k++)
    d_image[i*height*depth+j*depth+k]=(double)rcoeff*imsrc->mri[i][j][k];

 //calcul de la norme carree du gradient
 res=imx_norme_sqr_quad_gradient_3d_p(d_image, d_gradient, width, height, depth);

 //si erreur on libere la memoire et on quitte
 if (res) { FREE(d_image); FREE(d_gradient); return res; }

 //creation du grphic3d resultat
 imx_copie_param_3d_p(imsrc, imdst);

 max=0.0;
 for (i=0;i<width;i++)
  for (j=0;j<height;j++)
   for (k=0;k<depth;k++)
   {
    iVal = d_gradient[i*height*depth+j*depth+k]+0.5;
    if (max<iVal) max=iVal;
   }

 rcoeff=max/MAXMRI3D;
 imdst->rcoeff=(float)rcoeff;
 for (i=0;i<width;i++)
  for (j=0;j<height;j++)
   for (k=0;k<depth;k++)
    imdst->mri[i][j][k]=(TYPEMRI3D)floor(d_gradient[i*height*depth+j*depth+k]/rcoeff+0.5);

 imx_inimaxminpixel_3d_p(imdst);

 //liberations memoire
 FREE(d_image); FREE(d_gradient);

 return res;
}

/******************************************************************************
 **    filtre_moyenneur_3d() :
 **
 */
/*! calcule et affiche la norme quadratique carree du gradient d'une image
**
**  \param image :image source
**  \param gradient : image gradient (E/S)
**  \param width, height, depth : dimensions de l'image
**
************************************************************************/
int imx_norme_sqr_quad_gradient_3d_p(double *image, double *gradient, int width, int height, int depth)
{
 int i,j,k;
 double gx, gy, gz;
 int n_Deriv=0;

 for (i=0;i<width;i++)
 {
  for (j=0;j<height;j++)
  {
   for (k=0;k<depth;k++)
   {
    //derivee selon x
    n_Deriv=0; gx=0.0;
    if (i<(width-1)) {gx+=(image[(i+1)*height*depth+j*depth+k]-image[i*height*depth+j*depth+k])/2.0; n_Deriv++; }
    if (i>0) {gx+=(image[i*height*depth+j*depth+k]-image[(i-1)*height*depth+j*depth+k])/2.0; n_Deriv++; }
    gx/=n_Deriv;
    //derivee selon y
    n_Deriv=0; gy=0.0;
    if (j<(height-1)) {gy+=(image[i*height*depth+(j+1)*depth+k]-image[i*height*depth+j*depth+k])/2.0; n_Deriv++; }
    if (j>0) {gy+=(image[i*height*depth+j*depth+k]-image[i*height*depth+(j-1)*depth+k])/2.0; n_Deriv++; }
    gy/=n_Deriv;
    //derivee selon z
    n_Deriv=0; gz=0.0;
    if (k<(depth-1)) {gz+=(image[i*height*depth+j*depth+k+1]-image[i*height*depth+j*depth+k])/2.0; n_Deriv++; }
    if (k>0) {gz+=(image[i*height*depth+j*depth+k]-image[i*height*depth+j*depth+k-1])/2.0; n_Deriv++; }
    gz/=n_Deriv;
    //norme du gradient
    gradient[i*height*depth+j*depth+k]=gx*gx+gy*gy+gz*gz;
   }
  }
 }

 return 0;
}

/******************************************************************************
 **    filtre_diffusion_anisotrope_3d() :
 **
 */
/*! filtrage de diffusion anisotrope tel que preconise dans
**  "Magnetic Resonance Image Tissue Claissification Using a Partial Volume Model"
**  Shattuck et al.
**  avec kd=5.0 et t0=1/8
**
************************************************************************/
void filtre_diffusion_anisotrope_3d()
{
  int im_src, im_dst, nb_iter=0;
  int err=0;

  //choix des images
  im_src=GET_PLACE3D("Image Source");
  im_dst=GET_PLACE3D("Image Destination");

  //choix de la taille du filtre
  nb_iter=GET_INT("Nombre d'iterations", 0, &err);

  imx_filtre_diffusion_anisotrope_3d(im_src, im_dst, nb_iter);

  show_picture_3d(im_dst);
}

/******************************************************************************
 **    imx_filtre_diffusion_anisotrope_3d() :
 **
 */
/*! filtrage de diffusion anisotrope tel que preconise dans
**  "Magnetic Resonance Image Tissue Claissification Using a Partial Volume Model"
**  Shattuck et al.
**  avec kd=5.0 et t0=1/8
**
**  \param im_src :image source
**  \param im_dst : image resultat (E/S)
**  \param taille_filtre : nb_iter : nombre d'iterations
**
****************************************************************************/
int imx_filtre_diffusion_anisotrope_3d(int im_src, int im_dst, int nb_iter)
{
  grphic3d *imsrc=NULL, *imdst=NULL;
  int res=0;

  imsrc=ptr_img_3d(im_src);
  imdst=ptr_img_3d(im_dst);

  res=imx_filtre_diffusion_anisotrope_3d_p(imsrc, imdst, nb_iter);

  return res;
}

/******************************************************************************
 **    imx_filtre_diffusion_anisotrope_3d_p() :
 **
 */
/*! filtrage de diffusion anisotrope tel que preconise dans
**  "Magnetic Resonance Image Tissue Claissification Using a Partial Volume Model"
**  Shattuck et al.
**  avec kd=5.0 et t0=1/8
**
**  \param imsrc :image source
**  \param imdst : image resultat (E/S)
**  \param nb_iter : nombre d'iterations
**
****************************************************************************/
int imx_filtre_diffusion_anisotrope_3d_p(grphic3d *imsrc, grphic3d *imdst, int nb_iter)
{
 unsigned int i,j,k;
 TDimension width, height, depth;
 double *tmp_res_n0=NULL, *tmp_res_n1=NULL, *norm_grad=NULL;
 int iter, res=0;
 double kd=5.0, t0=1.0/8.0, val, term_f=0.0, max, iVal, rcoeff;

 if ((!imsrc)||(!imdst))
  { fprintf(stderr,"les deux images doivent etre non nulles dans imx_filtre_diffusion_anisotrope_3d_p\n"); return 1; }

 width=imsrc->width; height=imsrc->height; depth=imsrc->depth;

 //allocations memoire
 tmp_res_n0=CALLOC(width*height*depth, double);
 tmp_res_n1=CALLOC(width*height*depth, double);
 norm_grad=CALLOC(width*height*depth, double);

 //on remplit l'image double
 rcoeff=imsrc->rcoeff;
 for (i=0;i<width;i++)
  for (j=0;j<height;j++)
   for (k=0;k<depth;k++)
    tmp_res_n0[i*height*depth+j*depth+k]=(double)rcoeff*imsrc->mri[i][j][k];

 //filtrage de l'image
 for (iter=0;iter<nb_iter;iter++)
 {
  //calcul de la norme du gradient de l'image
  res=imx_norme_sqr_quad_gradient_3d_p(tmp_res_n0, norm_grad, width, height, depth);

  //calcul de l'image des coefficients de diffusion
  for (i=0;i<width;i++)
   for (j=0;j<height;j++)
    for (k=0;k<depth;k++)
    {
     val=exp(-norm_grad[i*height*depth+j*depth+k]/(kd*kd));
     norm_grad[i*height*depth+j*depth+k]=val;
    }

  //calcul de l'image filtree
  for (i=0;i<width;i++)
   for (j=0;j<height;j++)
    for (k=0;k<depth;k++)
    {
     //calcul du terme de diffusion
     term_f=0.0;
     val=tmp_res_n0[i*height*depth+j*depth+k];
     if ((i>0)&&(i<(width-1)))
     {
      term_f+=norm_grad[(i-1)*height*depth+j*depth+k]*(tmp_res_n0[(i+1)*height*depth+j*depth+k]-val);
      term_f+=norm_grad[(i+1)*height*depth+j*depth+k]*(tmp_res_n0[(i-1)*height*depth+j*depth+k]-val);
     }
     if ((j>0)&&(j<(height-1)))
     {
      term_f+=norm_grad[i*height*depth+(j-1)*depth+k]*(tmp_res_n0[i*height*depth+(j+1)*depth+k]-val);
      term_f+=norm_grad[i*height*depth+(j+1)*depth+k]*(tmp_res_n0[i*height*depth+(j-1)*depth+k]-val);
     }
     if ((k>0)&&(k<(depth-1)))
     {
      term_f+=norm_grad[i*height*depth+j*depth+k-1]*(tmp_res_n0[i*height*depth+j*depth+k+1]-val);
      term_f+=norm_grad[i*height*depth+j*depth+k+1]*(tmp_res_n0[i*height*depth+j*depth+k-1]-val);
     }
     term_f*=t0;
     //update de l'image res
     tmp_res_n1[i*height*depth+j*depth+k]=val+term_f;
    }

  //on reboucle si on n'a pas atteint le nombre d'iterations voulues
  if (iter<(nb_iter-1))
  {
   for (i=0;i<width;i++)
    for (j=0;j<height;j++)
     for (k=0;k<depth;k++)
      {
       tmp_res_n0[i*height*depth+j*depth+k]=tmp_res_n1[i*height*depth+j*depth+k];
      }
  }
 }

 //creation du grphic3d resultat
 imx_copie_param_3d_p(imsrc, imdst);

 max=0.0;
 for (i=0;i<width;i++)
  for (j=0;j<height;j++)
   for (k=0;k<depth;k++)
   {
    iVal = tmp_res_n1[i*height*depth+j*depth+k]+0.5;
    if (max<iVal) max=iVal;
   }

 rcoeff=max/MAXMRI3D;
 imdst->rcoeff=(float)rcoeff;
 for (i=0;i<width;i++)
  for (j=0;j<height;j++)
   for (k=0;k<depth;k++)
    imdst->mri[i][j][k]=(TYPEMRI3D)floor(tmp_res_n1[i*height*depth+j*depth+k]/rcoeff+0.5);

 imx_inimaxminpixel_3d_p(imdst);


 FREE(tmp_res_n0); FREE(tmp_res_n1); FREE(norm_grad);

 return res;
}

/******************************************************************************
 **    imx_query_filter_dimension() :
 **
 */
/*! demande de la taille du filtre pour les filtrages convolutif
**
**  \retval : taille du filtre (impaire)
****************************************************************************/
int imx_query_filter_dimension()
{
 char *quest[7];
 int i;
 int val;

 for(i=0;i<7;i++) quest[i]=CALLOC(80,char);

 strcpy(quest[0],"3 x 3 x 3");
 strcpy(quest[1],"5 x 5 x 5");
 strcpy(quest[2],"7 x 7 x 7");
 strcpy(quest[3],"9 x 9 x 9");
 strcpy(quest[4],"11 x 11 x 11");
 strcpy(quest[5],"13 x 13 x 13");
 strcpy(quest[6],"\0");

 val=GETV_QCM("Taille du filtre ?",(char **)quest);

 val=val*2+3;

 for(i=0;i<7;i++) FREE(quest[i]);

 return val;
}

/* -- remplit_objet() -------------------------------------------
**      Fonction qui remplit les objets labelises pour eviter
**      d'avoir deux surfaces confondues quand on a des objets
**      contenus l'un dans l'autre.
**
**      Usage   remplit_objet(imobj, label)
**
**    imobj  : grphic3d* : l'image contenant l'objet
**    label  : int       : label de l'objet dans l'image
**
*/
void remplit_objet(grphic3d *imobj, int label_obj)
{

  grphic3d *imtemp, *imtemp2;
  int i, j, k, l;

  /* On isole le label de l'objet dans l'image */
  for(i = 0; i < imobj->width; i++)
    for(j = 0; j < imobj->height; j++)
      for(k = 0; k < imobj->depth; k++)
        imobj->mri[i][j][k] = (TYPEMRI)(imobj->mri[i][j][k] == label_obj);

  imobj->icomp = 0;
  imobj->rcoeff = 1;
  imobj->max_pixel = 1;
  imobj->min_pixel = 0;
  imobj->cutoff_max = 1;
  imobj->cutoff_min = 0;

  imtemp = cr_grphic3d(imobj);
  imtemp2 = cr_grphic3d(imobj);

  /* On separe l'objet en composantes connexes par labelisation */
  imx_labelcroix_3d_p(imobj, imtemp);

  /* RAZ_mri(imobj); */
  for(i = 0; i < imobj->width; i++)
    for(j = 0; j < imobj->height; j++)
      for(k = 0; k < imobj->depth; k++)
        imobj->mri[i][j][k] = 0;

  for (l = 0; l < imtemp->max_pixel; l++)
  {
    for(i = 0; i < imtemp->width; i++)
      for(j = 0; j < imtemp->height; j++)
        for(k = 0; k < imtemp->depth; k++)
          imtemp2->mri[i][j][k] = (TYPEMRI)(imtemp->mri[i][j][k] == l+1);

    /* On comble les trous de la composante connexe */
    imx_hole_fill_3d_p(imtemp2, imtemp2);

    for(i = 0; i < imobj->width; i++)
      for(j = 0; j < imobj->height; j++)
        for(k = 0; k < imobj->depth; k++)
        {
          if (imobj->mri[i][j][k] == 0)
            if (imtemp2->mri[i][j][k]) imobj->mri[i][j][k] = label_obj;
        }
  }

  free_grphic3d(imtemp);
  free_grphic3d(imtemp2);

}

