/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*************************************************************************
**/	
/*!	\file:		ana_3d.c
***
***	project:	Imagix 1.01
***			
***
***	\brief description:	Contient les appel menu de l'analyse
***	
***	
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***     Last user action by: Mr. ARMSPACH on Sept 15th 1993
***
**************************************************************************/
#include <config.h>
#include <stdio.h> 
#include <stdlib.h> 	/* Needed when you use GET_FLOAT	*/
#include <math.h> 	/* Needed when you use math function	*/
#include <string.h> 	

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_3d.h"
#include "noyau/gui/imx_picture3d.h"
#include "noyau/imx_lang.h"
#include "math/imx_matrix.h"
#include "noyau/io/imx_export_file.h"
#include "noyau/io/imx_file.h"
#include "math/oper_3d.h"
#include "traitement/trai_3d.h"
#include "math/ana_3d.h" 
#include "outils/imx_sort.h"
#include "math/imx_math.h"




/********************************************
**   select_cell_3d() 
*/
/*! Selection de cellule en fonction de critere
**
********************************************/
void	select_cell_3d(void)
{ 
  int i;
  int choix;
  char *quest[5];
  int im_1,im_res;
  im_1=GET_PLACE3D(TEXT0059);
  im_res=GET_PLACE3D(TEXT0060);

  for(i=0;i<5;i++)
    quest[i]=CALLOC(80,char);
  strcpy(quest[0],TEXT0062);
  strcpy(quest[1],TEXT0063);
  strcpy(quest[2],TEXT0116);
  strcpy(quest[3],TEXT0117);
  strcpy(quest[4],"\0");
  choix=GETV_QCM("Votre choix",(const char **)quest);

  imx_select_cell_3d(im_1,im_res,choix);

  for(i=0;i<5;i++)
    FREE(quest[i]);
  
}

/********************************************
**   imx_select_cell_3d() 
*/
/*! Selection de cellule en fonction de critere
**  \param im_1 : numero de l'image source
**  \param im_2 : numero de l'image resultat
**  \param choix : critere de selection
					\li  0 : Surface
					\li  1 : Position centre de gravite
					\li  2 : La plus grande surface
					\li  3 : Par la valeur
**  
********************************************/
int 	imx_select_cell_3d(int im_1, int im_2,int choix)
{
	grphic3d *im1, *imres;
  	im1=ptr_img_3d(im_1);
  	imres=ptr_img_3d(im_2);
	
	return (imx_select_cell_3d_p(im1, imres,choix));
}


/********************************************
**   select_cell_3d() 
*/
/*! \ingroup Segmentation 
**  Selection de cellule sur une image labelisee en fonction de critere
**	\param im1 : image de depart
**  \param imres : image resultat (E/S)
**  \param choix : critere de selection
					\li  0 : Surface
					\li  1 : Position centre de gravite
					\li  2 : La plus grande surface
					\li  3 : Par la valeur
	\retval : 1
********************************************/
int	imx_select_cell_3d_p(grphic3d *im1, grphic3d *imres, int choix)
{ 
  int i,j,k,l,err=0;
  int wdth,hght,dpth;
  long surfinf,surfsup,limindexinf,limindexsup;
  int *index;
  double *tabsurf;
  long numerocell/*,isurf,icengra,jcengra*/;
  long max,nbr;

  tabsurf=alloc_dvector(LIM_LABEL);
  index  =alloc_ivector(LIM_LABEL);
  for (i=0;i<LIM_LABEL;i++) tabsurf[i]=0;

  wdth=im1->width;
  hght=im1->height;
  dpth=im1->depth;


  switch ((int)choix) {
  
	case 0: 	/*  Surface */
		surfinf=GET_INT(TEXT0064,2,&err);
		surfsup=GET_INT( TEXT0065,2,&err);
		if(err!=0) {
		  free_dvector(tabsurf,LIM_LABEL);
		  free_ivector(index,LIM_LABEL);
		  return 0; }
                nbr=0;
                for (i=0;i<wdth;i++)
                  for (j=0;j<hght;j++)
                    for (k=0;k<dpth;k++) {
		      numerocell=im1->mri[i][j][k];
                      if ( numerocell > 0  && numerocell < LIM_LABEL) {
                         tabsurf[numerocell]++;
			 nbr=MAXI(nbr,numerocell);
			 }
                      }

                tabsurf[0]=0.0; nbr=nbr+1;
				indexed_quicksort(tabsurf,nbr,index);

                limindexinf=0;
		        limindexsup=0;
                for (i=0;i<nbr;i++) {
		  if( tabsurf[index[i]] <= surfinf) limindexinf=i;
		  if( tabsurf[index[i]] <= surfsup) limindexsup=i;
		  }

                for (i=0;i<wdth;i++)
                  for (j=0;j<hght;j++)
                    for (k=0;k<dpth;k++) {
		      numerocell=im1->mri[i][j][k];
		      if (numerocell!=0) {
		        imres->mri[i][j][k]=0;
		        for(l=limindexinf;l<=limindexsup;l++) {
                          if ( index[l] == numerocell )
			    imres->mri[i][j][k]=(TYPEMRI3D)floor(numerocell);
                          }
                        }
                      else
			imres->mri[i][j][k]=0;
                      }

                imx_copie_param_3d_p(im1,imres);
                show_picture_3d(imres->pos);
		break;

	case 1:		/*  Centre de gravite */

		break;
	
	case 2:	  /* La plus grande surface   */
                nbr=0;
                for (i=0;i<wdth;i++)
                  for (j=0;j<hght;j++)
                    for (k=0;k<dpth;k++) {
		      numerocell=im1->mri[i][j][k];
                      if ( numerocell > 0  && numerocell < LIM_LABEL) {
                         tabsurf[numerocell]++;
			 nbr=MAXI(nbr,numerocell);
			 }
                      }

                tabsurf[0]=0.; nbr=nbr+1;
                indexed_quicksort(tabsurf,nbr,index);
 
                for (i=0;i<wdth;i++)
                  for (j=0;j<hght;j++)
                    for (k=0;k<dpth;k++) {
		      numerocell=im1->mri[i][j][k];
		      if (numerocell!=0) {
		        imres->mri[i][j][k]=0;
                          if ( index[nbr-1] == numerocell )
			    imres->mri[i][j][k]=(TYPEMRI3D)floor(numerocell);
                        }
                      else
			imres->mri[i][j][k]=0;
                      }

                imx_copie_param_3d_p(im1,imres);
		show_picture_3d(imres->pos);
		break;
	
	case 3:		/* La valeur d'une zone */
		max=GET_INT(TEXT0117,0,&err);
                imx_select_cell_value_3d_p(im1,imres,max);
                imx_copie_param_3d_p(im1,imres);
		show_picture_3d(imres->pos);
		break;
	
	default:
		fprintf(stderr,"Warning: in Cell selection\n");
		break;

  }

  free_dvector(tabsurf,LIM_LABEL);
  free_ivector(index,LIM_LABEL);
  
  return 1;

}

/********************************************
**  imx_select_cell_value_3d()
*/
/*! \brief Fonction qui ne garde que la cellule donc la valeur est
**    egale valeur
**  \param im_1 : numero de l'image source
**  \param im_res : numero de l'image dest
**  \param valeur : la valeur
**  \retval 1  
********************************************/
int    imx_select_cell_value_3d(int im_1, int im_res, long int valeur)
{
  grphic3d  *im1,*imres;
  
  im1=ptr_img_3d(im_1);
  imres=ptr_img_3d(im_res);
  
  imx_select_cell_value_3d_p(im1,imres,valeur);  

  return(1);
}


/*********************************************
**  imx_select_cell_value_3d_p()
*/
/*!  \brief Fonction qui ne garde que la cellule donc la valeur est
**   egale a valeur
**  \param im1 : image source
**  \param imres : image dest (E/S)
**  \param valeur : la valeur   
**  \retval 1
********************************************/
int    imx_select_cell_value_3d_p(grphic3d *im1, grphic3d *imres, long int valeur)
{
  int i,j,k;
  int wdth,hght,dpth;
  
  wdth=im1->width;
  hght=im1->height;
  dpth=im1->depth;

  imx_copie_param_3d_p(im1,imres); 
   
  for (i=0;i<wdth;i++)
    for (j=0;j<hght;j++)
      for (k=0;k<dpth;k++)
      {
       if ((im1->mri[i][j][k])==valeur)
          {imres->mri[i][j][k]=1;}
        else
          {imres->mri[i][j][k]=0;}
      }
  
  imres->max_pixel=1;
  imres->min_pixel=1;
  imres->cutoff_max=1;
  imres->cutoff_min=0;
  imres->icomp=0;
  imres->rcoeff=1.0;
  return(1);
}


/*******************************************
** --  fct_volume_cell_3d_p() ----------------------------
**/
/*!  Fonction qui calcul le volume des celulles.
**   Renvoie un tableau (pointeur sur un tableau) dont le ieme element est egal
**   au nombre de voxels de l'image ayant la valeur (valeur numerique, pas 
**   valeur reelle) i.
**   Le premier element du tableau indique le nombre d'element dans le tableau.
**   ATTENTION: le nombre d'element dans le tableau est limite a LIM_LABEL.
**   La place memoire du tableau est allouee dans cette fonction mais la liberation
**   de cette place est a la charge de la fonction appelante
**    
********************************************/
long    *fct_volume_cell_3d_p(grphic3d *im1)
{
  int i,j,k;
  long *tab;
  long max;
  int wdth,hght,dpth;

  /*allocation de la memoire du tableau a sa taille max possible et mise a zero*/
  tab=CALLOC(_lim_label+1,long);
 
  wdth=im1->width;
  hght=im1->height;
  dpth=im1->depth;
  
  max=0;
  for (i=0;i<wdth;i++)
    for (j=0;j<hght;j++)
      for (k=0;k<dpth;k++)
      {
       if (im1->mri[i][j][k]<=_lim_label)
       {
         (tab[im1->mri[i][j][k]])++;
         if (im1->mri[i][j][k]>max) {max=im1->mri[i][j][k];}
       }
      }   
  tab[0]=max;
  
  /*reallocation de la memoire du tableau a la bonne taille*/
  tab=REALLOC(tab,max+1,long);
  
  return(tab);
}

/*******************************************
** --  imx_statistique_3d() ----------------------------
*/
/*! Calcul statistique sur une image a partir d'une zone d'interet
********************************************/
int   imx_statistique_3d(int im_1, int mask3d_, float aminh, float amaxh, int *inb, float *moy, float *ecty, float *skew, float *kurt, float *mediane, float *per5quantile, float *per95quantile)
{
  grphic3d *im1, *mask;
  im1 = ptr_img_3d(im_1);
  mask = ptr_img_3d(mask3d_);
  
  return imx_statistique_3d_p(im1, mask, aminh, amaxh, inb, moy, ecty, skew, kurt, mediane, per5quantile, per95quantile);
}

int   imx_statistique_3d_p(grphic3d* im1, grphic3d* mask3d, float aminh, float amaxh, int *inb, float *moy, float *ecty, float *skew, float *kurt, float *mediane, float *per5quantile, float *per95quantile)

{ 
  int i,j,k;
  int wdth,hght,dpth;
  long *ibuf;
  float rcoeff;
  double a1;
  float tmp;
  double std_tmp, moy_tmp;

  wdth=im1->width;
  hght=im1->height;
  dpth=im1->depth;
  rcoeff=im1->rcoeff;

  ibuf=CALLOC(_lim_label,long int);

/*   Calcul du nombre de points, de la moyenne et de l'ecart type  */
  *moy=0.0;
  *ecty=0.0;
  std_tmp=0.0;
  moy_tmp=0.0;
  *inb=0;
  for (i=0;i<wdth;i++)
    for (j=0;j<hght;j++)
     for (k=0;k<dpth;k++)
      {
      if (mask3d->mri[i][j][k]!=0)
        {
        a1=im1->mri[i][j][k]*rcoeff;
        if (a1>=aminh&&a1<=amaxh)
	  {
	  // *moy=*moy+a1;
	  // *ecty=*ecty+a1*a1;
	  moy_tmp+=a1;
	  std_tmp+=a1*a1;
	  *inb=*inb+1;
	  if (*inb<=_lim_label) ibuf[*inb]=im1->mri[i][j][k];
	  }
        }
      }
   moy_tmp/=(double)(*inb);
   std_tmp/=(double)(*inb);
   std_tmp = std_tmp - moy_tmp*moy_tmp;
   std_tmp=(float)pow((double)std_tmp,(double) 0.5);

   *moy = (float) moy_tmp;
   *ecty=(float) std_tmp;
  // *moy=(*moy)/(*inb);
  // *ecty=(*ecty)/(*inb);
  // *ecty=*ecty-(*moy)*(*moy);
  // *ecty=(float)pow((double)*ecty,(double) 0.5);

/*   Calcul de la mediane, et quantile 5%   */
  if ((*inb)<=_lim_label)
    {
    tri_rapide(ibuf,0,*inb);
    *mediane=ibuf[*inb/2]*rcoeff;
    i=(int)floor(5.0*(*inb)/100.0);
    *per5quantile=ibuf[i]*rcoeff;
    *per95quantile=ibuf[*inb-i]*rcoeff;
    }
  else
    {
    if(MESG_DEBUG)PUT_WARN(" Impossible to calculate median, too much points\n");
    *mediane=0.;
    *per5quantile=0.;
    *per95quantile=0.;
    }
  
/*   Calcul Skewness et Kurtosis  */
  *skew=0.0;
  *kurt=0.0;
  for (i=0;i<wdth;i++)
    for (j=0;j<hght;j++)
     for (k=0;k<dpth;k++)
      {
      if (mask3d->mri[i][j][k]!=0)
        {
        a1=im1->mri[i][j][k]*rcoeff;
        if (a1>=aminh&&a1<=amaxh)
	  {
	  tmp=(a1-*moy);
	  *skew=(float)(*skew+pow(tmp,(double) 3.0));
	  *kurt=(float)(*kurt+pow(tmp,(double) 4.0));
	  }
        }
      }
  if (*ecty==0)
    {
    *skew=0.0;
    *kurt=0.0;
    }
  else
    {
    *skew=*skew/((*inb)*(*ecty)*(*ecty)*(*ecty));
    *kurt=(float)(*kurt/((*inb)*(*ecty)*(*ecty)*(*ecty)*(*ecty))-3.);
    }

  FREE(ibuf);

  return(1);

}

/*******************************************
** --  imx_histogram_3d_p() ----------------------------
*/
/*! Recherche histogramme de l'image 3D
**
**  \param im_1      :	image
**  \param roi       : 	masque de traitement (mettre a 1 pour totalite image)
**  \param xmin,xmax :	minimun et maximum pour l'histogramme
**  \param Nbpas     :	nombre de pas de l'histogramme
**  \param x,y       :	Tableu resultat de l'histogramme
**              Cumul l'histogramme dans y Il faut donc remmetre y
**              a zero si l'on souhaite un histogramme non cumule
********************************************/
int	imx_histogram_3d_p(grphic3d *im1, grphic3d *roi, float xmin, float xmax, long int Nbpas, float *x, float *y)
{ 
  int 	i,j,k,wdth,hght,dpth;
  float reso,a1;
  float rcoeff;

  wdth=im1->width;
  hght=im1->height;
  dpth=im1->depth;
  rcoeff=im1->rcoeff;

/*   Creation du tableu Y   */
  reso= ynist(xmin,xmax,Nbpas);
  if (reso==0) { return(0); } 
  for (i=0;i<wdth;i++)
    for (j=0;j<hght;j++)
      for (k=0;k<dpth;k++)
      {
      if (roi->mri[i][j][k]!=0)
      {
        a1=im1->mri[i][j][k]*rcoeff;
	    if (a1>=xmin&&a1<=xmax)
        {
	    ystog(y,Nbpas,xmin,reso,a1,1.);
      }
    }
  }

/*  Creation du tableau X   */
  for (i=0;i<Nbpas;i++)
    {
    x[i]=xmin+i*reso;
    }

  return(1);
}

/******************************************
**  programme cardiac_mednuc_volume_3d()
*/
/*!
**  Permet de calculer le volume d'une image de Medecine
**  nucleaire cardiac.
**  ----------
**  |       / |
**  |      /  |
**  | S1  /   |
**  |    /\ S2|
**  |   /  \  |
**  |  /    \ |
**  | /  S3  \|
**  ----------
******************************************/
void cardiac_mednuc_volume_3d()
{
 int i,j,k;
 int im_1;
 grphic3d *im1;
 long S1,S2,S3;
 double a;
 int w,h,d;

 im_1=GET_PLACE3D(TEXT0001);
 
 im1=ptr_img_3d(im_1);
 
 a=(double) im1->height/ (double) im1->width;
 S1=S2=S3=0;

 w = im1->width;
 h = im1->height;
 d = im1->depth;

 for (i=0;i<w; i++) 
   for (j=0;j<h; j++)
     {
	 if (j<(h-a*i) ) { 
	   for (k=0;k <d; k++)
         if(im1->mri[i][j][k] != 0 ) S1=S1+1; 
       }
     if (a*i>j && j>(h-a*i) ) { 
        for (k=0;k <d; k++)
	       if(im1->mri[i][j][k] != 0 ) S2=S2+1; 
       }
     if (a*i<j && j>(h-a*i) ) { 
         for (k=0;k <d; k++)
	       if(im1->mri[i][j][k] != 0 ) S3=S3+1; 
       }
	 }  
  PUTI(" Surface (pixels) S1 =_",S1," _");
  PUTI("S2 =_",S2," _");
  PUTI("S3 =_",S3," ");
	    
}

/********************************************
**   select_multi_cell_3d()
*/
/*! Selection de plusieurs cellules
    en fonction de critere
**
********************************************/
void  select_multi_cell_3d(void)
{
  int i;
  int choix;
  char *quest[3];
  int im_1,im_res;
  im_1=GET_PLACE3D(TEXT0059);
  im_res=GET_PLACE3D(TEXT0060);

  for(i=0;i<3;i++)
    quest[i]=CALLOC(80,char);
  strcpy(quest[0],TEXT0062);
  strcpy(quest[1],TEXT0117);
  strcpy(quest[2],"\0");
  choix=GETV_QCM("Votre choix",(const char **)quest);
  imx_select_multi_cell_3d(im_1,im_res,choix);

  for(i=0;i<3;i++)
    FREE(quest[i]);

}
/********************************************
**   imx_select_multi_cell_3d()
*/
/*! Selection de plusieurs cellules
**  im_1 :numero de l'image source
**  im_2 :numero de l'image resultat
**  choix : parametre
**      0:surface
**      1:valeur
********************************************/
int   imx_select_multi_cell_3d(int im_1, int im_2,int choix)
{
  grphic3d *im1, *imres;
  im1=ptr_img_3d(im_1);
  imres=ptr_img_3d(im_2);

  return (imx_select_multi_cell_3d_p(im1, imres,choix));
}

int  imx_select_multi_cell_3d_p(grphic3d *im1, grphic3d *imres, int choix)
{
  int val,err=0;
  int i,j,k;
  int choix2;
  BOOL loop=TRUE;
//char *quest[3];
char quest[3][256];
  int wdth,hght,dpth;
    
//  for(i=0;i<3;i++)
	memset((char*)quest,0,256*3);
//quest[i]=CALLOC(256,char);

  strcpy(quest[0],TEXT1740);//oui
  strcpy(quest[1],TEXT1750);//non
  strcpy(quest[2],"\0");
  
  wdth=im1->width;
  hght=im1->height;
  dpth=im1->depth;

  imx_copie_param_3d_p(im1,imres);
  for (i=0;i<wdth;i++)
    for (j=0;j<hght;j++)
      for (k=0;k<dpth;k++)
        imres->mri[i][j][k]=0;
     
   
  switch(choix){
    case 0: /*  surface      */
    
    break;

    case 1: /*  valeur    */
       while(loop){
        val=GET_INT(TEXT0117,0,&err);//valeur
        for (i=0;i<wdth;i++)
          for (j=0;j<hght;j++)
            for (k=0;k<dpth;k++)
              if ((im1->mri[i][j][k])==val)
                imres->mri[i][j][k]=val;
        imx_inimaxminpixel_3d_p(imres);
        show_picture_3d(imres->pos);      
        choix2=GETV_QCM("Nouvelle valeur",(const char **)quest);
        if(choix2==1)
          loop=FALSE;
       }
    break;
  }

//  for(i=0;i<3;i++)
//    FREE(quest[i]);

  return 1;
}

