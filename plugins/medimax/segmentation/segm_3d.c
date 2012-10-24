/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*!************************************************************************
***	
***	\file:		segm_3d.c
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
#include <stdio.h> 
#include <stdlib.h> 	/* Needed when you use GET_FLOAT	*/
#include <math.h> 	/* Needed when you use math function	*/
#include <time.h>
#include <string.h>

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_3d.h"
#include "noyau/gui/imx_picture3d.h"
#include "noyau/imx_lang.h" 
#include "traitement/trai_3d.h"
#include "math/ana_3d.h"
#include "math/oper_3d.h"
#include "noyau/mani_3d.h"
#include "segmentation/segm_3d.h"
#include "segmentation/segm_seuillage_3d.h"
#include "segmentation/otsu_3d.h"

/* DOXYGEN SECTION */
/*! \defgroup Segmentation
	\brief module de segmentation d'image medecine nucleaire
*/

/* END DOXYGEN SECTION */



/********************************************
**   imx_threh_img_3d() 
*/
/*!
** Seuillage d'une image a l'aide des valeurs vrais
** Garde la valeur de tout ce qui est > bas et < haut
** \param bas : seuil bas (valeur vraie de l'image)
** \param haut: seuil haut (valeur vraie de l'image)
** \param im_1 : image seuille
** \param im_res : resultat du seuillage (E/S)
** \retval 1
********************************************/
int	imx_threh_img_3d(int im_1, int im_res, float bas, float haut)
{
  grphic3d *im1,*imres;

  im1=ptr_img_3d(im_1);
  imres=ptr_img_3d(im_res);

  imx_threh_img_3d_p(im1,imres,bas,haut);

  return(1);
} 

/********************************************
**   imx_threh_img_3d_p() 
*/
/*! \ingroup Segmentation
**
** \brief Segmentation de l'image avec 2 seuils
**
** Seuillage d'une image a l'aide de 2 seuils (bas, haut).
** Garde la valeur de tout ce qui est > bas et < haut
** \param bas : seuil bas (valeur vraie de l'image)
** \param haut: seuil haut (valeur vraie de l'image)
** \param im_1 : image seuille
** \param im_res : resultat du seuillage (E/S)
** \retval 1
********************************************/
int	imx_threh_img_3d_p(grphic3d *im1, grphic3d *imres, float bas, float haut)
{
  int i,j,k;
  int w,h,d;
  TYPEMRI3D bas_mri, haut_mri;
  TYPEMRI3D ***im1MRI=NULL, ***imresMRI=NULL;

  //calcul des seuils sur le tableau image
  if ((bas/im1->rcoeff)<(double)MINMRI3D) bas_mri=MINMRI3D;
  else if ((bas/im1->rcoeff)>(double)MAXMRI3D) bas_mri=MAXMRI3D;
  else bas_mri  = (TYPEMRI3D) floor(bas/im1->rcoeff);

  if ((haut/im1->rcoeff)>(double)MAXMRI3D) haut_mri=MAXMRI3D;
  else if ((haut/im1->rcoeff)<(double)MINMRI3D) haut_mri=MINMRI3D;
  else  haut_mri = (TYPEMRI3D) ceil(haut/im1->rcoeff);

  //copie des parametres
  imx_copie_param_3d_p(im1,imres);

  //recuperation des donnees image
  w = (int)im1->width; h = (int)im1->height; d = (int)im1->depth;
  im1MRI=im1->mri; imresMRI=imres->mri;

  //seuillage
  for (i=0;i<w;i++) 
    for (j=0;j<h;j++) 
      for (k=0;k<d;k++)
        {
      	  if ((im1MRI[i][j][k]>=bas_mri)&&(im1MRI[i][j][k]<=haut_mri))
			      imresMRI[i][j][k]=im1MRI[i][j][k];
        	else
			      imresMRI[i][j][k]=0;
	      }

    return(1);

} 

/*************************************************
** -- opt_threh_3d() -------------------------------
**
**  Optimisation automatique du seuillage
**
**
*************************************************/
void   opt_threh_3d(void)
{

    int im_1,im_res;

    im_1=GET_PLACE3D(TEXT0015);
    im_res=GET_PLACE3D(TEXT0006);

    imx_opt_threh_3d(im_1,im_res);
/*    imx_iniparaimg_3d(im_res);*/     
    show_picture_3d(im_res);

    return;
}

/*************************************************
** -- imx_opt_threh_3d() -----------------------------
**
**  Optimisation automatique du seuillage
**  im_1      image a filtrer
**  im_res    image resultat
**
*************************************************/
int      imx_opt_threh_3d(int im_1, int im_res)
{
     grphic3d *im1,*imres;
     int x,y,z,i,j,flag,thres,GO,GB,sum;
     int w,h,d;
     float ihist[256],hist[256],sum1;

     imres=ptr_img_3d(im_res);
     im1=ptr_img_3d(im_1);
     imx_copie_param_3d_p(im1,imres);

     /* Histogramme   */
     w = im1->width;
     h = im1->height;
     d = im1->depth;

	 memset(ihist,0,256);
     sum = 0;
     for (x=0;x<w;x++)
       for (y=0;y<h;y++)
         for (z=0;z<d;z++)
	     {
	     j=im1->mri[x][y][z]*255/im1->max_pixel;
         j=MAXI(j,0);
		 ihist[j]=ihist[j]+1;
	     sum++;
         }
     
     for (i=0;i<=255;i++)
      {  hist[i]=ihist[i]/sum; }

     GO=255;
     GB=0;

	 for (y=0;y<=255;y++)
	 {
	 	sum1=0.0;
	 	j=0;
	 	for (x=-15;x<=15;x++)
	    {
	 	    j++;
	 	    if ((y-x)>=0)
			sum1=sum1+hist[y-x];    
       	}
     	hist[y]=sum1/j;
	 }
     y=2;
     flag=0;
     thres=0;
     while (flag==0 && y<254)
	 {
	 	if ((hist[y-1]>=hist[y]) && (hist[y]<hist[y+1]))
        {
	    	 flag=1;
	    	 thres=y;
        }
        y++;
     }
	 printf("thres=%d\n",thres);

     for (x=0;x<w;x++)
      for(y=0;y<h;y++)
		for (z=0;z<d;z++)
	    	{
	 		if ((im1->mri[x][y][z]*255/im1->max_pixel)>thres)
			imres->mri[x][y][z]=GO;
            else imres->mri[x][y][z]=GB;  /*norm  GB*/
            }
  /*   copie presque conforme rem quequiestimp????.*/

  imx_inimaxminpixel_3d_p(imres);
  imres->x3dc=im1->x3dc;
  imres->y3dc=im1->y3dc;
  imres->z3dc=im1->z3dc;
  
     return(1);
}

/*******************************************
** --  raz_img_partiel_3d_p() ----------------------------
**
**  Remise a zero partiel de l image pointez par im1
********************************************/
/*int	raz_img_partiel_3d_p(grphic3d *im1, int i_old, int j_old, int k_old)
{
  int i,j,k,wdth,hght,dpth;

  wdth=MINI(im1->width,MAX_WIDTH_3D);
  hght=MINI(im1->height,MAX_HEIGHT_3D);
  dpth=MINI(im1->depth,MAX_DEPTH_3D);

  for (i=0;i<wdth;i++)
    if (i<i_old) {
      for (j=0;j<hght;j++)
	if (j<j_old)
          for (k=k_old;k<dpth;k++)
            im1->mri[i][j][k]=0;
        else
	  for (k=0;k<dpth;k++)
	    im1->mri[i][j][k]=0;
       }
     else
       for (j=0;j<hght;j++)
	 for (k=0;k<dpth;k++)
	   im1->mri[i][j][k]=0;

   return(0);
}*/

/*********************************************************************/
/****** SEGMENTATION DU CERVEAU PAR RECHERCHE DE DISCONTINUITE *******/
/*********************************************************************/

/**********************************
**-- seg_cerv_grow_dis_3d()
** segmentation du cerveau par croissance
** et recherche de discontinuite
************************************/
int seg_cerv_grow_dis_3d(void)
{
 int im_1,im_res;
 int *mousepos3d,err=1,x0=0,y0=0,z0=0;

  im_res=GET_PLACE3D(TEXT0006);
 
 /*     Image   				*/
  PUT_WNDMSG(TEXT0162);
  mousepos3d=XGET_MOUSEPOS_3D(&err);
  CLEAR_WNDMSG();
  if (err)
    {
    /*ERREUR dans le GET_MOUSEPOS_3D*/
    PUT_WARN(TEXT0161);
	return (0);
    }
  
  x0=mousepos3d[3];y0=mousepos3d[4];z0=mousepos3d[5];
  im_1=mousepos3d[2];
  
  imx_seg_cerv_grow_dis_3d(im_1,im_res,x0,y0,z0);
  
  show_picture_3d(im_res);
  
  FREE(mousepos3d);
  return(1);
  
}

/*************************************
**-- imx_seg_cerv_grow_dis_3d()
** segmentation du cerveau par croissance
** et recherche de discontinuite
**************************************/
int imx_seg_cerv_grow_dis_3d(int im_deb, int im_res, int x0, int y0, int z0)
{
 grphic3d *imdeb,*imres;
 
 imdeb=ptr_img_3d(im_deb);
 imres=ptr_img_3d(im_res);
 

 imx_seg_cerv_grow_dis_3d_p(imdeb,imres,x0,y0,z0); 
 
 imx_inimaxminpixel_3d_p(imres);
 
 return(1);
 
}

/*************************************
**-- imx_moyenne_3d_p()
** moyenne d'une image
**************************************/
void imx_moyenne_3d_p(grphic3d *imdeb, grphic3d *imres)
{
 int i,j,k,t,ii,jj,kk,x,y,z,nb;
 int w,h,d;
 long moy;
 grphic3d *imtemp;
 
  /* creation de l'image temporaire imtemp */
  imtemp=cr_grphic3d(imdeb);
  imx_copie_param_3d_p(imdeb,imtemp);

 w = imdeb->width;
 h = imdeb->height;
 d = imdeb->depth;

 t=2;
  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      for (k=0;k<d;k++)
       {
        nb=0;moy=0;
        for (ii=-t;ii<=t;ii++)
          {
           x=i+ii;
           if (x>=0 && x<w)
           for (jj=-t;jj<=t;jj++)
             {
              y=j+jj;
              if (y>=0 && y<h)
              for (kk=-t;kk<=t;kk++)
              {
               z=k+kk;
               if (k>=0 && k<d)
               {moy=moy+imdeb->mri[x][y][z];nb++;}
              }
             }   
          }
        imtemp->mri[i][j][k]=(TYPEMRI3D)floor(moy/(long)nb);   
       }          
 
 
 imx_copie_param_3d_p(imtemp,imres);
 imx_copie_3d_p(imtemp,imres);
 free_grphic3d(imtemp);
 
}


/*************************************
**-- imx_seg_cerv_grow_dis_3d_p()
** segmentation du cerveau par croissance
** et recherche de discontinuite
*************************************/
int imx_seg_cerv_grow_dis_3d_p(grphic3d *imdeb, grphic3d *imres, int x0, int y0, int z0)
{
 grphic3d* imtemp,*imtemp2;
 int i,j,k;
 int w,h,d;
 float pasparam,param1,param2,vol1,vol2;
 float paramdeb,paramfin,voldeb,volfin;
 float difvolmax;
 float dx,dy,dz;
 float tab[4];
 long  nb;
 int  time1,time2,timedeb,timefin;
 
 dx=imdeb->dx;dy=imdeb->dy;dz=imdeb->dz;
  
 time1=time(NULL);timedeb=time1;
  
 /* creation de l'image temporaire imtemp */
  imtemp=cr_grphic3d(imdeb);
  imx_copie_param_3d_p(imdeb,imtemp);
 
 /* creation de l'image temporaire imtemp2 */
  imtemp2=cr_grphic3d(imdeb);
  imx_copie_param_3d_p(imdeb,imtemp2);


 /*moyenne de l'image de depart*/
 /*imx_median_3d_p(imdeb,imtemp);*/
 /*imx_moyenne_3d_p(imdeb,imtemp);*/
 imx_copie_3d_p(imdeb,imtemp);

 time2=time(NULL);
 printf("temps initialisation de l'image : %d \n",time2-time1);
 

  paramdeb=0.0;voldeb=0.0,volfin=0.0; 
  paramfin=5;
   
  pasparam=(float)0.002;
 do
  {
   param1=paramdeb;vol1=voldeb;
   difvolmax=0;   
   do
    {
      param2=param1+pasparam;
      if (param2!=paramfin)
       {
        /*mise a zero de imtemp2*/
          w = imtemp2->width;
          h = imtemp2->height;
          d = imtemp2->depth;

          for (i=0;i<w;i++)
           for (j=0;j<h;j++)
             for (k=0;k<d;k++)
                imtemp2->mri[i][j][k]=0;
  
         tab[0]=param2;
         growing_cerv_fonction_3d_p(x0,y0,z0,imtemp,imtemp2,param2);
    
        /*calcul du volume*/
         nb=0;
          w = imtemp2->width;
          h = imtemp2->height;
          d = imtemp2->depth;
          for (i=0;i<w;i++)
            for (j=0;j<h;j++)
              for (k=0;k<d;k++)
                 if (imtemp2->mri[i][j][k]>0) 
                   {nb++;}  
          vol2=dx*dy*dz*(float)nb;
        }
        else {vol2=volfin;}  
       
       if (vol2-vol1>difvolmax && vol1>500000 && vol1<1500000) {paramdeb=param1;voldeb=vol1;difvolmax=vol2-vol1;}
   
       /*printf("param1:%f vol1:%f param2:%f vol2:%f \n",param1,vol1,param2,vol2);*/
       param1=param2;vol1=vol2;
       
    }while (param2<paramfin && vol2<1500000);
    paramfin=paramdeb+pasparam;volfin=voldeb+difvolmax;
   /*printf("RESULTAT: paramdeb:%f voldeb:%f paramfin:%f volfin:%f \n",paramdeb,voldeb,paramfin,volfin);*/ 
   pasparam=(float)((paramfin-paramdeb)/4.0);  
  }while (pasparam>0.0001);
              
  
 /*calcul final avec param=paramdeb*/ 
   /*printf("Calcul final avec %f \n",paramdeb);*/

   /*mise a zero de imtemp2*/
   w = imtemp2->width;
   h = imtemp2->height;
   d = imtemp2->depth;
     for (i=0;i<w;i++)
       for (j=0;j<h;j++)
         for (k=0;k<d;k++)
            imtemp2->mri[i][j][k]=0;
  
    tab[0]=paramdeb;
    growing_cerv_fonction_3d_p(x0,y0,z0,imtemp,imtemp2,paramdeb);

 time1=time(NULL);
 printf("temps pour recherche discontinuite: %d \n",time1-time2);

 
 
  imx_copie_param_3d_p(imtemp2,imres);
  imx_copie_3d_p(imtemp2,imres); 
  
  free_grphic3d(imtemp);free_grphic3d(imtemp2);
  
  timefin=time(NULL);
  
  printf("temps total: %d \n",timefin-timedeb);
 
  return(1);
}


/**************************************************************************/
/********* SEGMENTATION DU CERVEAU PAR CONTOURS SPHERIQUES ****************/
/**************************************************************************/


/**********************************
**-- grow_spher_cerv_3d()
** extraction de la ROI correspondante au cerveau
**   option: 0 normal
**           1 suppression de la moelle
************************************/
int grow_spher_cerv_3d(int option)
{
 int im_1,im_res;
 int *mousepos3d,err=1,x0,y0,z0;

  im_res=GET_PLACE3D(TEXT0006);
 
 /*     Image   				*/
  PUT_WNDMSG(TEXT0162);
  mousepos3d=XGET_MOUSEPOS_3D(&err);
  CLEAR_WNDMSG();
  if (err)
    {
    /*ERREUR dans le GET_MOUSEPOS_3D*/
    PUT_ERROR(TEXT0161);
    return(0);
	}
  
  x0=mousepos3d[3];y0=mousepos3d[4];z0=mousepos3d[5];
  im_1=mousepos3d[2];

  imx_grow_spher_cerv_3d(im_1,im_res,x0,y0,z0,option);
  
  show_picture_3d(im_res);
  
  FREE(mousepos3d);
  return(1);  
}

/*************************************
**-- imx_grow_spher_cerv_3d()
**   extraction de la ROI correspondante au cerveau
**   option: 0 normal
**           1 suppression de la moelle
**************************************/
int imx_grow_spher_cerv_3d(int im_deb, int im_res, int x0, int y0, int z0, int option)
{
 grphic3d *imdeb,*imres,*imtemp;
 
 imdeb=ptr_img_3d(im_deb);
 imres=ptr_img_3d(im_res);
 
 /* creation de l'image temporaire imtemp */
  imtemp=cr_grphic3d(imdeb);
  imx_copie_param_3d_p(imdeb,imtemp);

 imx_cerveau_to_ROI_3d_p(imdeb,imtemp,x0,y0,z0,option);
 
 imx_mul_3d_p(imdeb,imtemp,imres);
 
 imx_inimaxminpixel_3d_p(imres);
 
 free_grphic3d(imtemp);
 return(1);
 
}


/**********************************
**-- cerveau_to_ROI_3d()
** extraction de la ROI correspondante au cerveau
**   option: 0 normal
**           1 suppression de la moelle
************************************/
int cerveau_to_ROI_3d(int option)
{
 int im_1;
 int *mousepos3d,err=1,x0,y0,z0;

 
 /*     Image   				*/
  PUT_WNDMSG(TEXT0162);
  mousepos3d=XGET_MOUSEPOS_3D(&err);
  CLEAR_WNDMSG();
  if (err)
  {
   /*ERREUR dans le GET_MOUSEPOS_3D*/
   PUT_ERROR(TEXT0161);
   return(0);
  }

  x0=mousepos3d[3];y0=mousepos3d[4];z0=mousepos3d[5];
  im_1=mousepos3d[2];
  
  imx_cerveau_to_ROI_3d(im_1,0,x0,y0,z0,option);
  
  FREE(mousepos3d);
  return(1);  
}


/*************************************
**-- imx_cerveau_to_ROI_3d()
**   extraction de la ROI correspondante au cerveau
**   x0,y0,z0: point de depart appartenant au cerveau
**   option: 0 normal
**           1 suppression de la moelle
**************************************/
int imx_cerveau_to_ROI_3d(int im_deb, int im_res, int x0, int y0, int z0, int option)
{
 grphic3d *imdeb,*imres;
 
 imdeb=ptr_img_3d(im_deb);
 imres=ptr_img_3d(im_res);

 imx_cerveau_to_ROI_3d_p(imdeb,imres,x0,y0,z0,option);
 
 return 0;
}

/*************************************
**-- imx_cerveau_to_ROI_3d_p()
**   extraction de la ROI correspondante au cerveau
**   x0,y0,z0: point de depart appartenant au cerveau
**   option: 0 normal
**           1 suppression de la moelle
**************************************/
int imx_cerveau_to_ROI_3d_p(grphic3d *imdeb, grphic3d *imres, int x0, int y0, int z0, int option)
{
 grphic3d *imtemp,*imtemp2;
 int i,j,k,xg,yg,zg;
 int w,h,d;
 long moyk,moyj,moyi;
 long nb;
 double teta,phi,r,pi,rmoy,rpas;
 double coefi,coefj,coefk;
 int iteta,iphi,nteta,nphi,tfilt,rmax;
 int iiteta,iiphi;
 double *rayon,*rayon2,*tabcteta,*tabsteta,*tabcphi,*tabsphi;
 double rayonmax,rayonmin;
 int l,m;
 int iphi1,iphi2,iteta1,iteta2,dphi;
 
 double pasteta,pasphi;
 
 pi=(double)PI;
 
 
 /* creation de l'image temporaire imtemp */
  imtemp=cr_grphic3d(imdeb);
  imx_copie_param_3d_p(imdeb,imtemp);
 
  /*EXTRACTION DU CERVEAU*/ 
  /*growing_cerv_fonction_3d_p(x0,y0,z0,imdeb,imtemp,0.07);*/
  imx_seg_cerv_grow_dis_3d_p(imdeb,imtemp,x0,y0,z0); 
 
  
 /* creation de l'image temporaire imtemp2 */
  imtemp2=cr_grphic3d(imtemp); 
  imx_copie_param_3d_p(imdeb,imtemp2); 
  
 /*EXTRACTION DE LA SURFACE DU CERVEAU*/
 imx_ctr_binaire_3d_p(imtemp,imtemp2); 
     
 w = imtemp2->width;
 h = imtemp2->height;
 d = imtemp2->depth;

     
 /*CALCUL DU CENTRE DE GRAVITE DE LA SURFACE DU CERVEAU*/
 moyi=0;moyj=0;moyk=0;nb=0;
 for (i=0;i<w;i++)
   for (j=0;j<h;j++)
     for (k=0;k<d;k++)
      if (imtemp2->mri[i][j][k]>0)
        {
          nb++;
          moyi=moyi+i;
          moyj=moyj+j;
          moyk=moyk+k;
        }
        
 xg=(int)(moyi/nb);
 yg=(int)(moyj/nb);
 zg=(int)(moyk/nb);
 
 
 /*liberation de la memoire pour l'image temporaire imtemp2*/
 free_grphic3d(imtemp2);
 w = imdeb->width;
 h = imdeb->height;
 d = imdeb->depth;
   /*MISE EN SPHERIQUE DE LA SURFACE*/
 
   /* premiere mise en spherique pour le calcul du vrai rmax */
      /*definition de rmax par rapport a la distance aux bords de l'image*/
         rmax=xg;
         if (yg>rmax) {rmax=yg;}
         if (zg>rmax) {rmax=zg;}
         if (w-xg>rmax) {rmax=imdeb->width-xg;}
         if (w-yg>rmax) {rmax=imdeb->height-yg;}
         if (w-zg>rmax) {rmax=imdeb->depth-zg;}
         printf("rmax:%d \n",rmax);
   
      /*definition de nteta et nphi*/ 
         nteta=rmax/2;nphi=2*nteta;
         pasteta=pi/(double)nteta;pasphi=2*pi/(double)nphi; 
         printf("nteta: %d  nphi: %d  \n",nteta,nphi);
  
      /*allocation et calcul des tableau tabteta et tabphi*/
         tabcteta=CALLOC((nteta+1),double);
         tabsteta=CALLOC((nteta+1),double);
         tabcphi=CALLOC((nphi+1),double);
         tabsphi=CALLOC((nphi+1),double);
         for (iteta=0;iteta<=nteta;iteta++)
           {
             teta=pasteta*(double)iteta;
             tabcteta[iteta]=cos(teta);
             tabsteta[iteta]=sin(teta);
           }

         for (iphi=0;iphi<=nphi;iphi++)
           {
             phi=pasphi*(double)iphi;
             tabcphi[iphi]=cos(phi);
             tabsphi[iphi]=sin(phi);
           }
  
      /*calcul de rayonmax*/
         rayonmax=0;
   for (iteta=0;iteta<=nteta;iteta++)
     for (iphi=0;iphi<=nphi;iphi++)
       {
          coefi=tabsteta[iteta]*tabcphi[iphi];
          coefj=tabsteta[iteta]*tabsphi[iphi];
          coefk=tabcteta[iteta];
    
          r=(double)rmax; 
          w = imtemp->width;
          h = imtemp->height;
          d = imtemp->depth;
  
          do
             {
               i=xg+(int)(r*coefi);
               j=yg+(int)(r*coefj);
               k=zg+(int)(r*coefk);
               r=r-1;
             } while (i<0 || i>=w || j<0 || j>=h || k<0 || k>=d ); 
    
          rpas=1;
          do
             {
               do
                  {
                    r=r-rpas;
                    i=xg+(int)(r*coefi);
                    j=yg+(int)(r*coefj);
                    k=zg+(int)(r*coefk);
                  } while(imtemp->mri[i][j][k]==0 && r>0);
               r=r+rpas;
              rpas=rpas/2;
             } while(rpas>0.5 && r>0);
           if (r>rayonmax) rayonmax=r;
       }
   
   rayonmax=1.2*rayonmax;
   rmax=(int)rayonmax;
 
   FREE(tabcteta);FREE(tabsteta);FREE(tabcphi);FREE(tabsphi);
 
 /*deuxieme mise en spherique avec la bonne valeur de rmax*/
   /*definition de nteta et nphi*/ 
      nteta=4*rmax;nphi=2*nteta;
      pasteta=pi/(double)nteta;pasphi=2*pi/(double)nphi; 
      tfilt=(int)(nteta/40);
      printf("nteta: %d  nphi: %d  taille filtre: %d \n",nteta,nphi,tfilt);
 
   /*allocation dynamique de la place pour les 2 tableaux de rayons en spherique*/
      rayon=CALLOC((nteta+1)*(nphi+1),double);
      rayon2=CALLOC((nteta+1)*(nphi+1),double);
 
   /*allocation et calcul des tableau tabteta et tabphi*/
      tabcteta=CALLOC((nteta+1),double);
      tabsteta=CALLOC((nteta+1),double);
      tabcphi=CALLOC((nphi+1),double);
      tabsphi=CALLOC((nphi+1),double);
      for (iteta=0;iteta<=nteta;iteta++)
        {
         teta=pasteta*(double)iteta;
         tabcteta[iteta]=cos(teta);
         tabsteta[iteta]=sin(teta);
        }
  
      for (iphi=0;iphi<=nphi;iphi++)
        {
         phi=pasphi*(double)iphi;
         tabcphi[iphi]=cos(phi);
         tabsphi[iphi]=sin(phi);
        }
  
 
    for (iteta=0;iteta<=nteta;iteta++)
      for (iphi=0;iphi<=nphi;iphi++)
        {
          coefi=tabsteta[iteta]*tabcphi[iphi];
          coefj=tabsteta[iteta]*tabsphi[iphi];
          coefk=tabcteta[iteta];
    
          r=(double)rmax; 
          w = imtemp->width;
          h = imtemp->height;
          d = imtemp->depth;
          do
            {
               i=xg+(int)(r*coefi);
               j=yg+(int)(r*coefj);
               k=zg+(int)(r*coefk);
               r=r-1;
            } while (i<0 || i>=w || j<0 || j>=h || k<0 || k>=d ); 
    
          rpas=2;
          do
            {
              do
                {
                   r=r-rpas;
                   i=xg+(int)(r*coefi);
                   j=yg+(int)(r*coefj);
                   k=zg+(int)(r*coefk);
                } while(imtemp->mri[i][j][k]==0 && r>0);
              r=r+rpas;
              rpas=rpas/2;
            } while(rpas>0.5 && r>0);
          rayon[iteta+iphi*(nteta+1)]=r;
        }


 /*ELIMINATION DE LA MOELLE EPINIERE*/
 iiteta=0;iiphi=0;
 if (option==1)
   {
     printf("%s \n","Elimination de la moelle epiniere");
     /*recherche de l'angle correspondant a la moelle*/
        rayonmax=0;
        for (iteta=nteta/4;iteta<=3*nteta/4;iteta++)
          for (iphi=nphi/8;iphi<=3*nphi/8;iphi++)
            {
              if (rayon[iteta+iphi*(nteta+1)]>rayonmax) {rayonmax=rayon[iteta+iphi*(nteta+1)];iiteta=iteta;iiphi=iphi;}
            }
       printf("iiteta:%d  iiphi:%d \n",iiteta,iiphi);
     /*elimination de la moelle*/
       for (m=-1;m<=1;m=m+2)
       {
        iteta=iiteta;
          do
          { 
           /*Calcul du rayon mini, moyen et maxi pour ce iteta*/
             rmoy=0;rayonmax=0;rayonmin=1000;nb=0;
             for (iphi=0;iphi<=nphi;iphi++)
               {
                 r=rayon[iteta+iphi*(nteta+1)];
                 rmoy=rmoy+r;
                 if (r>rayonmax) rayonmax=r;
                 if (r<rayonmin) rayonmin=r;
                 nb++;
               }
             rmoy=rmoy/(double)nb;
            
            /*choix du rayon seuil en fonction de rmin,rmoy et rmax*/
            /*if (rayonmax>1.2*rmoy)
            rmoy=rayonmin+(rmoy-rayonmin)/2; 
            else 
            rmoy=rayonmax;*/
            rmoy=rayonmin+(rmoy-rayonmin)/2;
            
           
           /*recherche de iphi1 et iphi2*/
             iphi=iiphi;
             do
               {
                iphi--;
                r=rayon[iteta+iphi*(nteta+1)];
               } while (r>=rmoy && iphi>0);
             iphi1=iphi;
             iphi=iiphi;
             do
               {
                iphi++;
                r=rayon[iteta+iphi*(nteta+1)];
               } while (r>=rmoy && iphi<nphi/2);
             iphi2=iphi;
             
           /*Exploitation de la symetrie de la moelle*/
              dphi=MINI(iiphi-iphi1,iphi2-iiphi);
              iphi1=iiphi-dphi;iphi2=iiphi+dphi;
           
          /*suppression de la partie relative a la moelle dans la surface*/
          
             for (iphi=iphi1;iphi<=iphi2;iphi++)
               {
                 rayon[iteta+iphi*(nteta+1)]=rmoy;
               }
            
               
           iteta=iteta+m;
          } while (dphi>1 && iteta>=0 && iteta<=nteta);
        if (m==-1) {iteta1=iteta+1;} else {iteta2=iteta-1;}   
       }
   }
  
 
 /*FILTRAGE MAX*/
   for (iteta=0;iteta<=nteta;iteta++)
     for (iphi=0;iphi<=nphi;iphi++)
       {
         rmoy=0;
         for (i=-tfilt;i<=tfilt;i++)
           for (j=-tfilt;j<=tfilt;j++)
             {
               iiteta=iteta+i;
               iiphi=iphi+j;
               if (iiteta<0 || iiteta>nteta || iiphi<0 || iiphi>nphi)
                   {
                     if (iiteta<0) {iiteta=-iiteta;iiphi=iiphi-nphi/2;}
                     if (iiteta>nteta-1) {iiteta=iiteta-nteta;iiphi=iiphi-nphi/2;}
                     if (iiphi<0) {iiphi=nphi-iiphi;}
                     if (iiphi>nphi-1) {iiphi=iiphi-nphi;}
                   }
               if (rayon[iiteta+iiphi*(nteta+1)]>rmoy) {rmoy=rayon[iiteta+iiphi*(nteta+1)];}
             }
         rayon2[iteta+iphi*(nteta+1)]=rmoy; 
       }

    
 /*FILTRAGE MOYEN*/
   for (iteta=0;iteta<=nteta;iteta++)
     for (iphi=0;iphi<=nphi;iphi++)
       {
         rmoy=0;nb=0;
         for (i=-tfilt;i<=tfilt;i++)
           for (j=-tfilt;j<=tfilt;j++)
             {
               iiteta=iteta+i;
               iiphi=iphi+j;
               if (iiteta<0 || iiteta>nteta || iiphi<0 || iiphi>nphi)
                   {
                     if (iiteta<0) {iiteta=-iiteta;iiphi=iiphi-nphi/2;}
                     if (iiteta>nteta-1) {iiteta=iiteta-nteta;iiphi=iiphi-nphi/2;}
                     if (iiphi<0) {iiphi=nphi-iiphi;}
                     if (iiphi>nphi-1) {iiphi=iiphi-nphi;}
                   }
               nb++;
               rmoy=rmoy+rayon2[iiteta+iiphi*(nteta+1)];
             }
         rayon[iteta+iphi*(nteta+1)]=rmoy/(double)nb;
       }


/*MISE A ZERO DE LA ROI*/
  w = imtemp->width;
  h = imtemp->height;
  d = imtemp->depth;

  for (i=0;i<w;i++)
   for (j=0;j<h;j++)
     for (k=0;k<d;k++)
     imtemp->mri[i][j][k]=0; 
 
 
 /*MISE EN PLACE DE LA SURFACE  3D DU CERVEAU DANS LA ROI*/
  for (iteta=0;iteta<nteta;iteta++)
   for (iphi=0;iphi<nphi;iphi++)
   {
    r=rayon[iteta+iphi*(nteta+1)];
    for (l=(int)r-4;l<=(int)r;l++)
    {
     i=xg+(int)((double)l*tabsteta[iteta]*tabcphi[iphi]);
     j=yg+(int)((double)l*tabsteta[iteta]*tabsphi[iphi]);
     k=zg+(int)((double)l*tabcteta[iteta]);
     imtemp->mri[i][j][k]=1;
    } 
   }
 
 imx_hole_fill_3d_p(imtemp,imtemp);
 imx_copie_3d_p(imtemp,imres);
 imx_copie_param_3d_p(imtemp,imres);
 
 free_grphic3d(imtemp);
 FREE(rayon2);
 FREE(rayon);
 FREE(tabcteta);FREE(tabsteta);FREE(tabcphi);FREE(tabsphi);

 return(1); 
 
}

/**************************************************************
**
**  Segmentation du crane par otsu + labelisation + zone la plus grande
**  + remplir trou 
*************************************************/

/*************************************
**   otsu_head_3d()
**   
**************************************/
void otsu_head_3d(void)
{
    int im_deb,im_res;

    im_deb=GET_PLACE3D("Original Image ? ");
    im_res=GET_PLACE3D("Result Image ? ");

    imx_otsu_head_3d(im_deb,im_res);

    show_picture_3d(im_res);

    return;

}

/****************************************************
**  imx_otsu_head_3d(im_deb,im_res)
**
**
**
****************/
void imx_otsu_head_3d(int im_deb, int im_res)
{
    grphic3d *imdeb,*imres;

    imdeb=ptr_img_3d(im_deb);
    imres=ptr_img_3d(im_res);

    imx_otsu_head_3d_p(imdeb,imres);
}

/****************************************************
**  imx_otsu_head_3d_p(imdeb,imres)
**
**
**
**************************************************/
void imx_otsu_head_3d_p(grphic3d *imdeb, grphic3d *imres)
{
  int threshold[MAX_THRESHOLD_NB];
  int i,j,k;
  int wdth,hght,dpth;
  int wdthm1,hghtm1,dpthm1;
  int thr_nb,n;
  long *vol,volmax,ivolmax;
  grphic3d *imtemp1,*imtemp2;
  float tab[1];


  wdth=imdeb->width  ;  wdthm1=wdth-1;
  hght=imdeb->height ;  hghtm1=hght-1;
  dpth=imdeb->depth  ;  dpthm1=dpth-1;

  imtemp1=cr_grphic3d(imdeb);
  imtemp2=cr_grphic3d(imdeb);


/*  Segmentation otsu 2 classe puis image binaire */
  thr_nb=1;
  find_threshold_3d(imdeb,threshold,thr_nb,TRUE);
  for(i=0;i<wdth;i++)
     for(j=0;j<hght;j++)
        for(k=0;k<dpth;k++)
           {
           if(imdeb->mri[i][j][k]>threshold[0])
                imtemp1->mri[i][j][k]=1;
           else 
                imtemp1->mri[i][j][k]=0;
            }
	    
/*  labelisation   */
  n=1; 
  for (i=1;i<wdthm1;i++) 
    for (j=1;j<hghtm1;j++) 
      for (k=1;k<dpthm1;k++)
        {
      	if ((imtemp1->mri[i][j][k]>0) && (imtemp2->mri[i][j][k]==0)) {
          tab[0]=(float)(imtemp1->mri[i][j][k]);
	  growing_3d_p(i,j,k,imtemp1,imtemp2,fct_grow_label_3d,tab,1,n,0);
	  n=n+1;
	  }
        }

/*   Recherche de la plus grande zone connexe */
  /*Calcul des volumes*/
  vol=fct_volume_cell_3d_p(imtemp2);
  
  /*Calcul du volume max et de la cellule correspondante*/
  volmax=0;ivolmax=0;
  for (i=1;i<=vol[0];i++)
  {
    if (vol[i]>volmax)
     {volmax=vol[i];ivolmax=i;}
  }
  
  /*Selection de la cellule ayant le plus grand volume*/
  imx_select_cell_value_3d_p(imtemp2,imtemp2,ivolmax);
  
  /*remplissage*/
  imx_hole_fill_3d_cpc_p(imtemp2,imtemp2,4);
  
  /*et logique*/
  for (i=0;i<wdth;i++) 
    for (j=0;j<hght;j++) 
      for (k=0;k<dpth;k++)
        {
		imres->mri[i][j][k]=imdeb->mri[i][j][k]*imtemp2->mri[i][j][k];
        }

  imx_copie_param_3d_p(imdeb,imres);
  imx_inimaxminpixel_3d_p(imres);
          
/* liberation de la memoire utilisee par les images temporaires*/
  free_grphic3d (imtemp1);
  free_grphic3d (imtemp2);
  free (vol);
}

/********************************************
**   imx_threh_img_outliers_3d_p()
*/
/*! \ingroup Segmentation
**
** \brief Segmentation de l'image avec 2 seuils
**
** Seuillage d'une image a l'aide de 2 seuils (bas, haut).
** Garde la valeur de tout ce qui est < bas et ce qui est > haut
** \param bas : seuil bas (valeur vraie de l'image)
** \param haut: seuil haut (valeur vraie de l'image)
** \param im_1 : image seuille
** \param im_res : resultat du seuillage (E/S)
** \retval 1
********************************************/
int imx_threh_img_outliers_3d_p (grphic3d *im1, grphic3d *imres, float bas, float haut)
{
  int i,j,k;
  int w,h,d;
  TYPEMRI3D bas_mri, haut_mri;
  TYPEMRI3D ***im1MRI=NULL, ***imresMRI=NULL;

  //calcul des seuils sur le tableau image
  bas_mri  = (TYPEMRI3D) floor(bas/im1->rcoeff);
  haut_mri = (TYPEMRI3D) ceil(haut/im1->rcoeff);

  //copie des parametres
  imx_copie_param_3d_p(im1,imres);

  //recuperation des donnees image
  w = (int)im1->width; h = (int)im1->height; d = (int)im1->depth;
  im1MRI=im1->mri; imresMRI=imres->mri;

  //seuillage
  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      for (k=0;k<d;k++)
        {
      	  if ((im1MRI[i][j][k]<=bas_mri)||(im1MRI[i][j][k]>=haut_mri))
			      imresMRI[i][j][k]=im1MRI[i][j][k];
        	else
			      imresMRI[i][j][k]=0;
	      }

    return 0;
}

