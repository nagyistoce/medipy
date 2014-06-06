/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*!--------------------------------------------------------------------------
***
***	\file:		morpho_3d.c
***
***	project:	Imagix 2.01
***
***
***	\brief description:    Fichier pour les operateurs morphologiques 3D
***
***
***	Copyright (c) 1993-2000, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***
***--------------------------------------------------------------------------*/

#include <config.h>
#include <stdio.h>
#include <stdlib.h> 	/* Needed when you use GET_FLOAT	*/
#include <math.h> 	/* Needed when you use math function	*/
#include <string.h>



#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/gui/imx_picture3d.h"
#include "noyau/imx_3d.h"
#include "noyau/mani_3d.h"
#include "noyau/imx_lang.h"
#include "outils/queues_3d.h"
#include "math/oper_3d.h"
#include "morpho/morpho_3d.h"
#include "noyau/io/imx_file.h"







/* FIFOs */

extern void create_p_fifo_3d(FIFO *fifo);

extern void insert_p_fifo_3d(FIFO *fifo, TPINFO info);

extern int fifo_empty_3d(FIFO fifo);

extern void free_p_fifo_3d(FIFO *fifo);

extern TPINFO fifo_first_3d(FIFO *fifo);



/*******************************************

** --  eros_3d() ----------------------------

**

**  procedure erosion 3D

**  maj par Cyril Feruglio

**  derniere maj le 16-04-1998

********************************************/



int imx_eros_3d (int im_1, int im_res, int numas, int niter);

int imx_eros_3d_p (grphic3d *im1, grphic3d *imres, int numas, int niter);

extern XPoints_3d *masque_3d (int numas);

extern int imx_inimaxminpixel_3d_p (grphic3d *im1);

int imx_dilat_3d (int im_1, int im_res, int numas, int niter);

int imx_dilat_3d_p (grphic3d *im1, grphic3d *imres, int numas, int niter);

int imx_GradM_3d (int im_1, int im_res, int numas, int numgrad, int niter);

int imx_GradM_3d_p (grphic3d *im1, grphic3d *imres, int numas, int numgrad, int niter);

int imx_open_3d (int im_1, int im_res, int numas, int niter);

int imx_open_3d_p (grphic3d *im1, grphic3d *imres, int numas, int niter);

int imx_close_3d (int im_1, int im_res, int numas, int niter);

int imx_close_3d_p (grphic3d *im1, grphic3d *imres, int numas, int niter);

int imx_tophat_3d (int im_1, int im_res, int numas, int black_white_tophat, int niter);

int imx_tophat_3d_p (grphic3d *im1, grphic3d *imres, int numas, int black_white_tophat, int niter);

int imx_reconstruction_geod_dilat_3d (int im_1, int im_2, int im_res, int answer_nr);

int imx_reconstruction_geod_dilat_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres, int numas);

int imx_reconstruction_geod_eros_3d (int im_1, int im_2, int im_res, int answer_nr);

int imx_reconstruction_geod_eros_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres, int numas);

int imx_sharpen_morpho_3d (int im_1, int im_res, int numas, int niter);

int imx_sharpen_morpho_3d_p (grphic3d *im1, grphic3d *imres, int numas, int niter);

int imx_contrast_morpho_3d (int im_1, int im_res, int numas, int niter);

int imx_contrast_morpho_3d_p (grphic3d *im1, grphic3d *imres, int numas, int niter);



void     eros_3d(void)

{

  int im_1,im_res,niter,i,err=0;

  char *quest[12];

  int answer_nr;



  /*     Image   				*/

  im_1=GET_PLACE3D(TEXT0011);

  im_res=GET_PLACE3D(TEXT0006);

  niter= GET_INT(TEXT0007, 2, &err);



  /* Question concerning type of masque  */

  if ((niter>0) && (!err))

    {

      for (i=0;i<12;i++)

	quest[i]=CALLOC(80,char);

      strcpy(quest[0],TEXT0159);

      strcpy(quest[1],TEXT0160);

      strcpy(quest[2],TEXT0191);

      strcpy(quest[3],TEXT0238);

      strcpy(quest[4],TEXT0474);

      strcpy(quest[5],TEXT0475);

      strcpy(quest[6],TEXT0476);

      strcpy(quest[7],TEXT0477);

      strcpy(quest[8],TEXT0478);

      strcpy(quest[9],TEXT0479);

      strcpy(quest[10],TEXT0481);

      strcpy(quest[11],"\0");



      /* Put a question box at position 150,200. */

      answer_nr=GETV_QCM(TEXT0043,(char **)quest);

      for (i=0;i<12;i++)

	FREE(quest[i]);



      if (answer_nr != 10) {

	imx_eros_3d(im_1,im_res,answer_nr,niter);

	show_picture_3d(im_res);

      }

      else

	PUT_MESG("Annulation du traitement");

    }

  else

    PUT_MESG("Erreur (niter)");



  return;

}



/*******************************************

** --  imx_eros_3d() ----------------------------

**

** Erosion d'une image

**

********************************************/

int     imx_eros_3d(int im_1, int im_res, int numas, int niter)

{

  grphic3d *im1,*imres;



  im1=ptr_img_3d(im_1);

  imres=ptr_img_3d(im_res);



  imx_eros_3d_p(im1,imres,numas,niter);

  return(1);

}



/*********************************************************************************/





/*******************************************

** --  imx_eros_3d_p() ----------------------------

**

** Erosion d'une image

**

** reecrit par Feruglio Cyril

** derniere maj le :  16 04 1998

********************************************/

int     imx_eros_3d_p(grphic3d *im1, grphic3d *imres, int numas, int niter)

{



  int i;                 /* indice de deplacement dans le masque */

  int imas,jmas,kmas; /* indices de deplacement du masque dans l'image */

  int taille_mas=0;        /* taille du masque (cube:26, plans:18, croix:6) */

  int ii, bord=1;

  XPoints_3d *mas;

  grphic3d *imtemp;

  long min;

  int w,h,d;



  /* creation de l'image temporaire imtemp */



  /* allocation memoire pour l'image temporaire */

  imtemp=cr_grphic3d(im1);



  /* Copie du contenu de l'image im1 dans imtemp*/

  imx_copie_3d_p(im1,imtemp);



  /* Copie des parametres de im1 (nom,tailles ...) dans l'image resultat */

  imx_copie_param_3d_p(im1,imres);



  w = imtemp->width;

  h = imtemp->height;

  d = imtemp->depth;





  /* initialisation de l'image resultat a 0 */

  for(imas=0; imas<w; imas++)

    for(jmas=0; jmas<h; jmas++)

      for(kmas=0; kmas<d; kmas++)

	{imres->mri[imas][jmas][kmas]=0;}



  /* allocation de l'espace memoire pour le masque (mas) d'apres numas(type de masque) */

  mas= masque_3d(numas);



  /* initialisation de la taille du masque */

  switch(numas){

  case 0 : taille_mas= 26; 	/* cas du voisinage cube (26 elements)  */

    break;



  case 1 : taille_mas= 6; 	/* cas du voisinage croix (6 elements)  */

    break;



  case 2 : taille_mas= 18; 	/* cas du voisinage plans (18 elements) */

    break;



  case 3 : taille_mas= 8; 	/* cas du voisinage croix (8 elements)  */

    break;



  case 4 : taille_mas=80; bord=2;

    break;



  case 5 : taille_mas=13;

    break;



  case 6 : taille_mas=3;

    break;



  case 7 : taille_mas=9;

    break;



  case 8 : taille_mas=4;

    break;



  case 9 : taille_mas=40; bord=2;

    break;

  }



  /* Calcul de l'erosion */

  for (i=0;i<niter;i++)

    {

      for(imas=bord;imas<(w-bord);imas++)

        for(jmas=bord;jmas<(h-bord);jmas++)

          for(kmas=bord;kmas<(d-bord);kmas++)

            {

	      min=imtemp->mri[imas][jmas][kmas];

	      for(ii=0;ii<=taille_mas;ii++)

	 	{

		  if (imtemp->mri[imas+mas[ii].x][jmas+mas[ii].y][kmas+mas[ii].z]<min)

		    {

		      min=imtemp->mri[imas+mas[ii].x][jmas+mas[ii].y][kmas+mas[ii].z];

		    }

		}

	      imres->mri[imas][jmas][kmas]=(TYPEMRI3D)min;

	    }

      imx_copie_3d_p(imres,imtemp);

    }





  /* liberation de la memoire utilisee par les images temporaires*/

  free_grphic3d (imtemp);



  /* sauvegarde de l'image erodee */

  imx_inimaxminpixel_3d_p(imres);



  /* liberation de la memoire utilisee par le masque */

  free(mas);



  return(1);

}



/*******************************************

** --  dilat_3d() ----------------------------

**

**  procedure dilatation 3D

**  maj par Cyril Feruglio

**  derniere maj le 16-04-1998

********************************************/

void    dilat_3d(void)

{

  int im_1,im_res,niter,i,err=0;

  char *quest[12];

  int answer_nr;



  /*     Image   				*/

  im_1=GET_PLACE3D(TEXT0011);

  im_res=GET_PLACE3D(TEXT0006);

  niter= GET_INT(TEXT0007, 2, &err);



  /* Question concerning type of masque  */

  if ((niter>0) && (!err))

    {

      for (i=0;i<12;i++)

	quest[i]=CALLOC(80,char);

      strcpy(quest[0],TEXT0159);

      strcpy(quest[1],TEXT0160);

      strcpy(quest[2],TEXT0191);

      strcpy(quest[3],TEXT0238);

      strcpy(quest[4],TEXT0474);

      strcpy(quest[5],TEXT0475);

      strcpy(quest[6],TEXT0476);

      strcpy(quest[7],TEXT0477);

      strcpy(quest[8],TEXT0478);

      strcpy(quest[9],TEXT0479);

      strcpy(quest[10],TEXT0481);

      strcpy(quest[11],"\0");



      /* Put a question box at position 150,200. */

      answer_nr=GETV_QCM(TEXT0043,(char **)quest);

      for (i=0;i<12;i++)

	FREE(quest[i]);



      if ( answer_nr != 10) {

	imx_dilat_3d(im_1,im_res,answer_nr,niter);

	show_picture_3d(im_res);

      }

      else

	PUT_MESG("Annulation du traitement");

    }

  else

    PUT_MESG("Erreur (niter)");



  return;

}



/*******************************************

** --  imx_dilat_3d() ----------------------------

**

** Dilatation d'une image

**

********************************************/

int     imx_dilat_3d(int im_1, int im_res, int numas, int niter)

{

  grphic3d *im1,*imres;



  im1=ptr_img_3d(im_1);

  imres=ptr_img_3d(im_res);



  imx_dilat_3d_p(im1,imres,numas,niter);

  return(1);

}



/*******************************************

** --  imx_dilat_3d_p() ----------------------------

**

** Dilatation d'une image

**

** reecrit par Feruglio Cyril

** derniere maj le :  16 04 1998

********************************************/

int     imx_dilat_3d_p(grphic3d *im1, grphic3d *imres, int numas, int niter)

{



  int i;                 /* indice de deplacement dans le masque */

  int imas,jmas,kmas; /* indices de deplacement du masque dans l'image */

  int taille_mas=0;        /* taille du masque (cube:26, plans:18, croix:6) */

  int ii, bord=1;

  XPoints_3d *mas;

  grphic3d *imtemp;

  long max;

  int w,h,d;



  /* creation de l'image temporaire imtemp */

  /* allocation memoire pour l'image temporaire */

  imtemp=cr_grphic3d(im1);



  /* Copie du contenu de l'image im1 dans imtemp*/

  imx_copie_3d_p(im1,imtemp);



  /* Copie des parametres de im1 (nom,tailles ...) dans l'image resultat */

  imx_copie_param_3d_p(im1,imres);



  w=imtemp->width;

  h=imtemp->height;

  d=imtemp->depth;



  /* initialisation de l'image resultat a 0 */

  for(imas=0; imas<(w); imas++)

    for(jmas=0; jmas<(h); jmas++)

      for(kmas=0; kmas<(d); kmas++)

	{imres->mri[imas][jmas][kmas]=0;}



  /* allocation de l'espace memoire pour le masque (mas) d'apres numas(type de masque) */

  mas= masque_3d(numas);



  /* initialisation de la taille du masque */

  switch(numas){

  case 0 : taille_mas= 26; 	/* cas du voisinage cube (26 elements)  */

    break;



  case 1 : taille_mas= 6; 	/* cas du voisinage croix (6 elements)  */

    break;



  case 2 : taille_mas= 18; 	/* cas du voisinage plans (18 elements) */

    break;



  case 3 : taille_mas= 8; 	/* cas du voisinage croix (8 elements)  */

    break;



  case 4 : taille_mas=80; bord=2;

    break;



  case 5 : taille_mas=13;

    break;



  case 6 : taille_mas=3;

    break;



  case 7 : taille_mas=9;

    break;



  case 8 : taille_mas=4;

    break;



  case 9 : taille_mas=40; bord=2;

    break;

  }





  /* dilatation proprement dite */

  for (i=0;i<niter;i++)

    {

      for(imas=bord;imas<(w-bord);imas++)

        for(jmas=bord;jmas<(h-bord);jmas++)

          for(kmas=bord;kmas<(d-bord);kmas++)

            {

	      max=imtemp->mri[imas][jmas][kmas];

	      for(ii=0;ii<=taille_mas;ii++)

	 	{

		  if



		    (imtemp->mri[imas+mas[ii].x][jmas+mas[ii].y][kmas+mas[ii].z]>max)

		    {

		      max=imtemp->mri[imas+mas[ii].x][jmas+mas[ii].y][kmas+mas[ii].z];

		    }

		}

	      imres->mri[imas][jmas][kmas]=(TYPEMRI3D)max;

	    }

      imx_copie_3d_p(imres,imtemp);

    }







  /* liberation de la memoire utilisee par les images temporaires*/

  free_grphic3d (imtemp);



  /* sauvegarde de l'image dilatee */

  imx_inimaxminpixel_3d_p(imres);



  /* liberation de la memoire utilisee par le masque */

  free(mas);



  return(1);

}





/*******************************************

** --  GradM_3d() ----------------------------

**

**  procedure Gradient Morphologique 3D

**   par Fahima DJAOUI 99

**

********************************************/

void    GradM_3d(void)

{

  int im_1,im_res,niter,i,err=0;

  char *quest1[12];

  char *quest2[5];

  int answer_nr1;

  int answer_nr2;



  /*     Image   				*/

  im_1=GET_PLACE3D(TEXT0462);

  im_res=GET_PLACE3D(TEXT0006);





  niter= GET_INT(TEXT0007, 2, &err);



  /* Question concerning type of masque  */

  if ((niter>0) && (!err))

    {

      for (i=0;i<12;i++)

	quest1[i]=CALLOC(80,char );

      strcpy(quest1[0],TEXT0159);

      strcpy(quest1[1],TEXT0160);

      strcpy(quest1[2],TEXT0191);

      strcpy(quest1[3],TEXT0238);

      strcpy(quest1[4],TEXT0474);

      strcpy(quest1[5],TEXT0475);

      strcpy(quest1[6],TEXT0476);

      strcpy(quest1[7],TEXT0477);

      strcpy(quest1[8],TEXT0478);

      strcpy(quest1[9],TEXT0479);

      strcpy(quest1[10],TEXT0481);

      strcpy(quest1[11],"\0");

      /* Put a question box at position 150,200. */

      answer_nr1=GETV_QCM(TEXT0043,(char **)quest1);

      for (i=0;i<12;i++)

	FREE(quest1[i]);



      if (answer_nr1 != 10) {

	for (i=0;i<5;i++)

	  quest2[i]=CALLOC(80,char );

	strcpy(quest2[0],TEXT0459);

	strcpy(quest2[1],TEXT0460);

	strcpy(quest2[2],TEXT0461);

	strcpy(quest2[3],TEXT0481);

	strcpy(quest2[4],"\0");

        /* Put a question box at position 150,200. */

	answer_nr2=GETV_QCM(TEXT0043,(char **)quest2);

	for (i=0;i<5;i++)

	  FREE(quest2[i]);



	if (answer_nr2 != 3) {

	  imx_GradM_3d(im_1,im_res,answer_nr1, answer_nr2,niter);

	  show_picture_3d(im_res);

	}

	else

	  PUT_MESG("Annulation du traitement");

      }

      else

	PUT_MESG("Annulation du traitement");

    }

  else

    PUT_MESG("Erreur (niter)");



  return;

}



/*******************************************

** --  imx_GradM_3d() ----------------------------

**

** Gradient Morphologique d'une image 3D

**

********************************************/

int     imx_GradM_3d(int im_1, int im_res, int numas, int numgrad, int niter)

{

  grphic3d *im1,*imres;



  im1=ptr_img_3d(im_1);

  imres=ptr_img_3d(im_res);



  imx_GradM_3d_p(im1,imres,numas,numgrad,niter);

  return(1);

}



/*******************************************

** --  imx_GradM_3d_p() ----------------------------

**

** Gradient morphologique d'une image 3D

** ecrit le:   03 06 1999 by F.Djaoui

********************************************/

int     imx_GradM_3d_p(grphic3d *im1, grphic3d *imres, int numas, int numgrad, int niter)

{



  int i;                 /* indice de deplacement dans le masque */

  int imas,jmas,kmas; /* indices de deplacement du masque dans l'image */

  int taille_mas=0;        /* taille du masque (cube:26, plans:18, croix:6) */

  int ii, bord=1;

  XPoints_3d *mas;

  grphic3d *imtemp;

  long max, min;

  int w,h,d;



  /* creation de l'image temporaire imtemp */



  /* allocation memoire pour l'image temporaire */

  imtemp=cr_grphic3d(im1);



  /* Copie du contenu de l'image im1 dans imtemp*/

  imx_copie_3d_p(im1,imtemp);



  /* Copie des parametres de im1 (nom,tailles ...) dans l'image resultat */

  imx_copie_param_3d_p(im1,imres);



  w = imtemp->width;

  h = imtemp->height;

  d = imtemp->depth;





  /* initialisation de l'image resultat a 0 */

  for(imas=0; imas<(w); imas++)

    for(jmas=0; jmas<(h); jmas++)

      for(kmas=0; kmas<(d); kmas++)

	{imres->mri[imas][jmas][kmas]=0;

	}



  /* allocation de l'espace memoire pour le masque (mas) d'apres numas(type de masque) */

  mas= masque_3d(numas);



  /* initialisation de la taille du masque */

  switch(numas){

  case 0 : taille_mas= 26; 	/* cas du voisinage cube (26 elements)  */

    break;



  case 1 : taille_mas= 6; 	/* cas du voisinage croix (6 elements)  */

    break;



  case 2 : taille_mas= 18; 	/* cas du voisinage plans (18 elements) */

    break;



  case 3 : taille_mas= 8; 	/* cas du voisinage croix (8 elements)  */

    break;



  case 4 : taille_mas=80; bord=2;

    break;



  case 5 : taille_mas=13;

    break;



  case 6 : taille_mas=3;

    break;



  case 7 : taille_mas=9;

    break;



  case 8 : taille_mas=4;

    break;



  case 9 : taille_mas=40; bord=2;

    break;

  }





  /* Calcul du gradient */

  for (i=0;i<niter;i++)

    {

      for(imas=bord;imas<(w-bord);imas++)

        for(jmas=bord;jmas<(h-bord);jmas++)

          for(kmas=bord;kmas<(d-bord);kmas++)

            {max=imtemp->mri[imas][jmas][kmas];

	    min=imtemp->mri[imas][jmas][kmas];

	    for(ii=0;ii<=taille_mas;ii++)

	      {

	        if(imtemp->mri[imas+mas[ii].x][jmas+mas[ii].y][kmas+mas[ii].z]>max)

		  {

		    max=imtemp->mri[imas+mas[ii].x][jmas+mas[ii].y][kmas+mas[ii].z];

		  }



		if(imtemp->mri[imas+mas[ii].x][jmas+mas[ii].y][kmas+mas[ii].z]<min)

		  {

		    min=imtemp->mri[imas+mas[ii].x][jmas+mas[ii].y][kmas+mas[ii].z];

		  }



	      }

	    if (numgrad==0) /*dilatation - erosion*/

	      imres->mri[imas][jmas][kmas]=(TYPEMRI3D)floor(max- min);

	    if (numgrad==1) /*image initial - erosion*/

	      imres->mri[imas][jmas][kmas]=(TYPEMRI3D)floor(imtemp->mri[imas][jmas][kmas]- min);

	    if (numgrad==2) /*dilatation - image initial*/

	      imres->mri[imas][jmas][kmas]=(TYPEMRI3D)floor(max-imtemp->mri[imas][jmas][kmas]);

	    }

      /* sauvegarde de l'image  */

      imx_copie_3d_p(imres,imtemp);

    }





  /* liberation de la memoire utilisee par les images temporaires*/

  free_grphic3d (imtemp);



  /* sauvegarde de l'image resultat */

  imx_inimaxminpixel_3d_p(imres);



  /* liberation de la memoire utilisee par le masque */

  free(mas);



  return(1);

}



/*******************************************

** --  open_3d() ----------------------------

**

**  procedure ouverture 3D

**   par Nicodeme PAUL

**

********************************************/

void    open_3d(void)

{

  int im_1,im_res,niter,i,err=0;

  char *quest[12];

  int answer_nr;



  /*     Image   				*/

  im_1=GET_PLACE3D(TEXT0465);

  im_res=GET_PLACE3D(TEXT0006);

  niter= GET_INT(TEXT0007, 2, &err);



  /* Question concerning type of masque  */

  if ((niter>0) && (!err))

    {

      for (i=0;i<12;i++)

	quest[i]=CALLOC(80,char);

      strcpy(quest[0],TEXT0159);

      strcpy(quest[1],TEXT0160);

      strcpy(quest[2],TEXT0191);

      strcpy(quest[3],TEXT0238);

      strcpy(quest[4],TEXT0474);

      strcpy(quest[5],TEXT0475);

      strcpy(quest[6],TEXT0476);

      strcpy(quest[7],TEXT0477);

      strcpy(quest[8],TEXT0478);

      strcpy(quest[9],TEXT0479);

      strcpy(quest[10],TEXT0481);

      strcpy(quest[11],"\0");

      /* Put a question box at position 150,200. */

      answer_nr=GETV_QCM(TEXT0043,(char **)quest);

      for (i=0;i<12;i++)

	FREE(quest[i]);



      if (answer_nr != 10) {

	imx_open_3d(im_1,im_res,answer_nr,niter);

	show_picture_3d(im_res);

      }

      else

	PUT_MESG("Annulation du traitement");

    }

  else

    PUT_MESG("Erreur (niter)");



  return;

}



/*******************************************

** --  imx_open_3d() ----------------------------

**

** Ouverture d'une image

**  reecrite par Nicodeme PAUL Mars 2000

********************************************/

int     imx_open_3d(int im_1, int im_res, int numas, int niter)

{

  grphic3d *im1,*imres;



  im1=ptr_img_3d(im_1);

  imres=ptr_img_3d(im_res);



  imx_open_3d_p(im1,imres,numas,niter);



  return(1);

}



/*******************************************

** --  imx_open_3d_p() ----------------------------

**

** Ouverture d'une image

** ecrite par Nicodeme PAUL Mars 2000

********************************************/

int     imx_open_3d_p(grphic3d *im1, grphic3d *imres, int numas, int niter)

{

  grphic3d *imtemp;



  /* allocation memoire pour l'image temporaire */

  imtemp=cr_grphic3d(im1);



  /*eroder de im1 et mettre le resultat dans imtemp*/

  imx_eros_3d_p(im1,imtemp,numas,niter);



  /*dilater imtemp et mettre le resultat dans imres*/

  imx_dilat_3d_p(imtemp,imres,numas,niter);



  /* liberation de la memoire utilisee par les images temporaires*/

  free_grphic3d (imtemp);



  return(1);

}





/*******************************************

** --  close_3d() ----------------------------

**

**  procedure fermeture 3D

**   par Fahima DJAOUI 99

**

********************************************/

void    close_3d(void)

{

  int im_1,im_res,niter,i,err=0;

  /*

  int im_2,im_res2;

  */

  char *quest[12];

  int answer_nr;



  /*     Image   				*/

  im_1=GET_PLACE3D(TEXT0466);

  im_res=GET_PLACE3D(TEXT0006);

  niter= GET_INT(TEXT0007, 2, &err);



  /* Question concerning type of masque  */

  if ((niter>0) && (!err))

    {

      for (i=0;i<12;i++)

	quest[i]=CALLOC(80,char);

      strcpy(quest[0],TEXT0159);

      strcpy(quest[1],TEXT0160);

      strcpy(quest[2],TEXT0191);

      strcpy(quest[3],TEXT0238);

      strcpy(quest[4],TEXT0474);

      strcpy(quest[5],TEXT0475);

      strcpy(quest[6],TEXT0476);

      strcpy(quest[7],TEXT0477);

      strcpy(quest[8],TEXT0478);

      strcpy(quest[9],TEXT0479);

      strcpy(quest[10],TEXT0481);

      strcpy(quest[11],"\0");

      /* Put a question box at position 150,200. */

      answer_nr=GETV_QCM(TEXT0043,(char **)quest);

      for (i=0;i<12;i++)

	FREE(quest[i]);



      if (answer_nr != 10) {

	imx_close_3d(im_1,im_res,answer_nr,niter);

	show_picture_3d(im_res);

      }

      else

	PUT_MESG("Annulation du traitement");

    }

  else

    PUT_MESG("Erreur (niter)");



  return;

}





/*******************************************

** --  imx_close_3d() ----------------------------

**

** Fermeture d'une image

** By DJAOUI 1999

********************************************/

int     imx_close_3d(int im_1, int im_res, int numas, int niter)

{

  grphic3d *im1,*imres;



  im1=ptr_img_3d(im_1);

  imres=ptr_img_3d(im_res);



  imx_close_3d_p(im1,imres,numas,niter);



  return(1);

}





/*******************************************

** --  imx_close_3d_p() ----------------------------

**

** Ouverture d'une image

** ecrite par Nicodeme PAUL Mars 2000

********************************************/

int     imx_close_3d_p(grphic3d *im1, grphic3d *imres, int numas, int niter)

{

  grphic3d *imtemp;



  /* allocation memoire pour l'image temporaire */

  imtemp=cr_grphic3d(im1);



  /*eroder de im1 et mettre le resultat dans imtemp*/

  imx_dilat_3d_p(im1,imtemp,numas,niter);



  /*dilater imtemp et mettre le resultat dans imres*/

  imx_eros_3d_p(imtemp,imres,numas,niter);



  /* liberation de la memoire utilisee par les images temporaires*/

  free_grphic3d (imtemp);



  return(1);

}



/*******************************************

** --  tophat_3d() ----------------------------

**    difference entre une image et son ouverture ou

**    difference entre une image et sa fermeture

**   par Nicodeme PAUL Fevrier 2000

**

********************************************/

void    tophat_3d(void)

{

  int im_1,im_res,i,niter,err=0;

  char *quest1[12];

  char *quest2[4];

  int answer_nr1;

  int answer_nr2;



  /*     Image   				*/

  im_1=GET_PLACE3D(TEXT0467);

  im_res=GET_PLACE3D(TEXT0006);

  niter= GET_INT(TEXT0007, 2, &err);



  /* Question concerning type of masque  */

  if ((niter>0) && (!err))

    {

      for (i=0;i<12;i++)

	quest1[i]=CALLOC(80,char );

      strcpy(quest1[0],TEXT0159);

      strcpy(quest1[1],TEXT0160);

      strcpy(quest1[2],TEXT0191);

      strcpy(quest1[3],TEXT0238);

      strcpy(quest1[4],TEXT0474);

      strcpy(quest1[5],TEXT0475);

      strcpy(quest1[6],TEXT0476);

      strcpy(quest1[7],TEXT0477);

      strcpy(quest1[8],TEXT0478);

      strcpy(quest1[9],TEXT0479);

      strcpy(quest1[10],TEXT0481);

      strcpy(quest1[11],"\0");

      /* Put a question box at position 150,200. */

      answer_nr1=GETV_QCM(TEXT0043,(char **)quest1);

      for (i=0;i<12;i++)

	FREE(quest1[i]);



      if (answer_nr1 != 10) {

	for (i=0;i<4;i++)

          quest2[i]=CALLOC(80,char );

	strcpy(quest2[0],TEXT0463);

	strcpy(quest2[1],TEXT0464);

	strcpy(quest2[2],TEXT0481);

	strcpy(quest2[3],"\0");

	/* Put a question box at position 150,200. */

	answer_nr2=GETV_QCM(TEXT0043,(char **)quest2);

	for (i=0;i<4;i++)

	  FREE(quest2[i]);



	if (answer_nr2 != 2) {

	  imx_tophat_3d(im_1,im_res,answer_nr1,answer_nr2,niter);

	  show_picture_3d(im_res);

	}

	else

	  PUT_MESG("Annulation du traitement");

      }

      else

	PUT_MESG("Annulation du traitement");

    }

  else

    PUT_MESG("Erreur (niter)");



  return;

}



/*******************************************

** --  imx_tophat_3d() ----------------------------

**     difference entre une image et son ouverture ou

**     difference entre une image et sa fermeture

**     par Nicodeme PAUL Fevrier 2000

********************************************/

int     imx_tophat_3d(int im_1, int im_res, int numas, int black_white_tophat, int niter)

{

  grphic3d *im1,*imres;



  im1=ptr_img_3d(im_1);

  imres=ptr_img_3d(im_res);



  imx_tophat_3d_p(im1,imres,numas,black_white_tophat,niter);



  return(1);

}





/*******************************************

** --  imx_tophat_3d_p() ----------------------------

**     difference entre une image et son ouverture ou

**     difference entre une image et sa fermeture

**     par Nicodeme PAUL Fevrier 2000

********************************************/

int     imx_tophat_3d_p(grphic3d *im1, grphic3d *imres, int numas, int black_white_tophat, int niter)

{

  grphic3d *imtemp;

  int i,j,k;

  int signe=1;

  int w,h,d;



  /* allocation memoire pour l'image temporaire */

  imtemp=cr_grphic3d(im1);



  if (black_white_tophat==0) {

    imx_eros_3d_p(im1,imtemp,numas,niter);

    imx_dilat_3d_p(imtemp,imres,numas,niter);

  }

  else

    {

      imx_dilat_3d_p(im1,imtemp,numas,niter);

      imx_eros_3d_p(imtemp,imres,numas,niter);

      signe=-1;

    }



  w = imres->width;

  h = imres->height;

  d = imres->depth;



  /* liberation de la memoire utilisee par les images temporaires*/

  free_grphic3d (imtemp);



  /* difference d'images */

  for(i=0;i<(w);i++)

    for(j=0;j<(h);j++)

      for(k=0;k<(d);k++)

	{

	  imres->mri[i][j][k]=((im1->mri[i][j][k])-(imres->mri[i][j][k]))*signe;

	}



  /* sauvegarde de l'image resultat */

  imx_inimaxminpixel_3d_p(imres);



  return(1);

}





/*******************************************

** --  reconstruction_geod_dilat_3d() ----------------------------

**   reconstruction morphologique d'une image

**   a partir d'une autre image par dilatation

**   par Nicodeme PAUL Fevrier 2000

**

********************************************/



void reconstruction_geod_dilat_3d(void)

{

  int im_1,im_2,im_res,i;

  char *quest1[12];

  int answer_nr;



  /*     Image   				*/

  im_1=GET_PLACE3D(TEXT0472);

  im_2=GET_PLACE3D(TEXT0473);

  im_res=GET_PLACE3D(TEXT0006);



  /* Question concerning type of masque  */

  for (i=0;i<12;i++)

    quest1[i]=CALLOC(80,char );

  strcpy(quest1[0],TEXT0159);

  strcpy(quest1[1],TEXT0160);

  strcpy(quest1[2],TEXT0191);

  strcpy(quest1[3],TEXT0238);

  strcpy(quest1[4],TEXT0474);

  strcpy(quest1[5],TEXT0475);

  strcpy(quest1[6],TEXT0476);

  strcpy(quest1[7],TEXT0477);

  strcpy(quest1[8],TEXT0478);

  strcpy(quest1[9],TEXT0479);

  strcpy(quest1[10],TEXT0481);

  strcpy(quest1[11],"\0");

  /* Put a question box at position 150,200. */

  answer_nr=GETV_QCM(TEXT0043,(char **)quest1);

  for (i=0;i<12;i++)

    FREE(quest1[i]);



  if (answer_nr != 10) {

    imx_reconstruction_geod_dilat_3d(im_1,im_2,im_res,answer_nr);

    show_picture_3d(im_res);

  }

  else

    PUT_MESG("Annulation du traitement");



  return;

}



/*******************************************

** --  imx_reconstruction_geod_dilat_3d() ----------------------------

**   reconstruction morphologique d'une image

**   a partir d'une autre image par dilatation

**   par Nicodeme PAUL Fevrier 2000

**

********************************************/

int     imx_reconstruction_geod_dilat_3d(int im_1, int im_2, int im_res, int answer_nr)

{

  grphic3d *im1, *im2, *imres;



  im1=ptr_img_3d(im_1);

  im2=ptr_img_3d(im_2);

  imres=ptr_img_3d(im_res);



  imx_reconstruction_geod_dilat_3d_p(im1,im2,imres,answer_nr);



  return(1);

}





/*******************************************

** --  imx_reconstruction_geod_dilat_3d_p() ----------------------------

**   reconstruction morphologique d'une image

**   a partir d'une autre image par dilatation

**   par Nicodeme PAUL Fevrier 2000

**

********************************************/

int     imx_reconstruction_geod_dilat_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres, int numas)

{

  int i, j, k;

  int taille_mas=0, taille=0;        /* taille du masque (cube:26, plans:18, croix:6) */

  int ii, bord=1;

  int width,height,depth;

  char reponse='n';

  XPoints_3d *mas;

  grphic3d *imtemp;

  long max;



  width=im1->width;

  height=im1->height;

  depth=im1->depth;



  /* allocation memoire pour l'image temporaire */

  imtemp=cr_grphic3d(im2);



  /* Copie du contenu de l'image im2 dans imres*/

  imx_copie_3d_p(im2,imres);



  /* Copie du contenu de l'image im1 dans imtemp*/

  imx_copie_3d_p(im2,imtemp);



  /* allocation de l'espace memoire pour le masque (mas) d'apres numas(type de masque) */

  mas= masque_3d(numas);



  /* initialisation de la taille du masque */

  switch(numas){

  case 0 : taille_mas= 26; 	/* cas du voisinage cube (26 elements)  */

    break;



  case 1 : taille_mas= 6; 	/* cas du voisinage croix (6 elements)  */

    break;



  case 2 : taille_mas= 18; 	/* cas du voisinage plans (18 elements) */

    break;



  case 3 : taille_mas= 8; 	/* cas du voisinage croix (8 elements)  */

    break;



  case 4 : taille_mas=80; bord=2;

    break;



  case 5 : taille_mas=13;

    break;



  case 6 : taille_mas=3;

    break;



  case 7 : taille_mas=9;

    break;



  case 8 : taille_mas=4;

    break;



  case 9 : taille_mas=40; bord=2;

    break;

  }



  taille=taille_mas/2;



  /* Repeter jusqu'a stabilite : imtemp=imres*/

  do {



    /************* Forward pass **************/



    printf("------ Forward pass ------\n");





    for(k=bord;k<depth-bord;k++)

      for(j=bord;j<height-bord;j++)

        for(i=bord;i<width-bord;i++)

	  {

	    max=imres->mri[i][j][k];

            for(ii=0;ii<taille;ii++)

	      max=MAXI(imres->mri[i+mas[ii].x][j+mas[ii].y][k+mas[ii].z],max);

	    imres->mri[i][j][k]=(TYPEMRI3D)MINI(max,im1->mri[i][j][k]);

	  }



    /************* Backward pass **************/



    printf("------ Backward pass ------\n");



    for(k=depth-bord-1;k>bord-1;k--)

      for(j=height-bord-1;j>bord-1;j--)

        for(i=width-bord-1;i>bord-1;i--)

	  {

	    max=imres->mri[i][j][k];

            for(ii=taille;ii<taille_mas;ii++)

	      max=MAXI(imres->mri[i+mas[ii].x][j+mas[ii].y][k+mas[ii].z],max);

	    imres->mri[i][j][k]=(TYPEMRI3D)MINI(max,im1->mri[i][j][k]);

	  }

    /* verification de imres-imtemp */

    for (k=0; k< depth; k++)

      for (j=0; j<height; j++)

	for (i=0; i<width; i++) {

	  if (imres->mri[i][j][k] != imtemp->mri[i][j][k]) {

	    reponse='y';

	    /* Copie du contenu de l'image imres dans imtemp*/

	    imx_copie_3d_p(imres,imtemp);

	    break;

	  }

	  else

	    {

	      reponse='n';

	    }

	}

  }

  while (reponse == 'y');



  /* liberation de la memoire utilisee par les images temporaires*/

  free_grphic3d (imtemp);



  /* sauvegarde de l'image resultat */

  imx_inimaxminpixel_3d_p(imres);



  /* liberation de la memoire utilisee par le masque */

  free(mas);



  return(1);

}





/*******************************************

** --  imx_reconstruction_geod_dilat_bin_3d_p() ----------------------------

**   reconstruction morphologique d'une composante connexe

**   d'une image binaire a partir d'un marqueur

**   par Nicodeme PAUL Fevrier 2000

**

********************************************/

void     imx_reconstruction_geod_dilat_bin_3d_p(grphic3d *imtemp2, grphic3d *imres, int x3d, int y3d, int z3d, int im_res, int numas)

{

  int taille_mas=0;        /* taille du masque cube:26 */

  int ii, bord=1;

  int width,height,depth;

  int x,y,z,grey_level;

  XPoints_3d *mas;

  TPINFO info;

  FIFO  fifo;



  /* initialisation */

  /*  fifo.head=NULL;

  fifo.tail=NULL;

  */



  width=imtemp2->width;

  height=imtemp2->height;

  depth=imtemp2->depth;



  /* allocation de l'espace memoire pour le masque (mas) d'apres numas(type de masque) */

  mas= masque_3d(numas);



  /* initialisation de la taille du masque */

  switch(numas){

  case 0 : taille_mas= 26; 	/* cas du voisinage cube (26 elements)  */

    break;



  case 1 : taille_mas= 6; 	/* cas du voisinage croix (6 elements)  */

    break;



  case 2 : taille_mas= 18; 	/* cas du voisinage plans (18 elements) */

    break;



  case 3 : taille_mas= 8; 	/* cas du voisinage croix (8 elements)  */

    break;



  case 4 : taille_mas=80; bord=2;

    break;



  case 5 : taille_mas=13;

    break;



  case 6 : taille_mas=3;

    break;



  case 7 : taille_mas=9;

    break;



  case 8 : taille_mas=4;

    break;



  case 9 : taille_mas=40; bord=2;

    break;

  }





  create_p_fifo_3d(&fifo);

  info.coord_x=x3d;

  info.coord_y=y3d;

  info.coord_z=z3d;

  info.grey_level=imtemp2->mri[x3d][y3d][z3d];

  imres->mri[x3d][y3d][z3d]=(TYPEMRI3D)info.grey_level;

  imtemp2->mri[x3d][y3d][z3d]=0;

  insert_p_fifo_3d(&fifo, info);





  /* Repeter */

  while (!((fifo.head==NULL)&&(fifo.tail==NULL))) {

    info = fifo_first_3d(&fifo);

    x=info.coord_x;

    y=info.coord_y;

    z=info.coord_z;

    grey_level=info.grey_level;

    if (x<bord || y<bord || z<bord || x>=width-bord || y>=height-bord || z>=depth-bord)

      printf("Je suis arrive au bord\n");

    else {

      for(ii=0;ii<taille_mas;ii++) {

	if (imtemp2->mri[x+mas[ii].x][y+mas[ii].y][z+mas[ii].z]>0) {

	  info.coord_x=x+mas[ii].x;

	  info.coord_y=y+mas[ii].y;

	  info.coord_z=z+mas[ii].z;

	  info.grey_level=imtemp2->mri[x+mas[ii].x][y+mas[ii].y][z+mas[ii].z];

	  insert_p_fifo_3d(&fifo, info);

	  imres->mri[x+mas[ii].x][y+mas[ii].y][z+mas[ii].z]=(TYPEMRI3D)info.grey_level;

	  imtemp2->mri[x+mas[ii].x][y+mas[ii].y][z+mas[ii].z]=0;

	}

      }

    }

  }





  /* liberation de la memoire utilisee par le masque */

  free(mas);

  free_p_fifo_3d(&fifo);



  show_picture_3d(im_res);

}



/*******************************************

** --  reconstruction_geod_eros_3d() ----------------------------

**   reconstruction morphologique d'une image

**   a partir d'une autre image par erosion

**   par Nicodeme PAUL Fevrier 2000

**

********************************************/



void reconstruction_geod_eros_3d(void)

{

  int im_1,im_2,im_res,i;

  char *quest1[12];

  int answer_nr;



  /*     Image   				*/

  im_1=GET_PLACE3D(TEXT0472);

  im_2=GET_PLACE3D(TEXT0473);

  im_res=GET_PLACE3D(TEXT0006);



  /* Question concerning type of masque  */

  for (i=0;i<12;i++)

    quest1[i]=CALLOC(80,char );

  strcpy(quest1[0],TEXT0159);

  strcpy(quest1[1],TEXT0160);

  strcpy(quest1[2],TEXT0191);

  strcpy(quest1[3],TEXT0238);

  strcpy(quest1[4],TEXT0474);

  strcpy(quest1[5],TEXT0475);

  strcpy(quest1[6],TEXT0476);

  strcpy(quest1[7],TEXT0477);

  strcpy(quest1[8],TEXT0478);

  strcpy(quest1[9],TEXT0479);

  strcpy(quest1[10],TEXT0481);

  strcpy(quest1[11],"\0");

  /* Put a question box at position 150,200. */

  answer_nr=GETV_QCM(TEXT0043,(char **)quest1);

  for (i=0;i<12;i++)

    FREE(quest1[i]);



  if (answer_nr != 10) {

    imx_reconstruction_geod_eros_3d(im_1,im_2,im_res,answer_nr);

    show_picture_3d(im_res);

  }

  else

    PUT_MESG("Annulation du traitement");



  return;

}



/*******************************************

** --  imx_reconstruction_geod_eros_3d() ----------------------------

**   reconstruction morphologique d'une image

**   a partir d'une autre image par erosion

**   par Nicodeme PAUL Fevrier 2000

**

********************************************/

int     imx_reconstruction_geod_eros_3d(int im_1, int im_2, int im_res, int answer_nr)

{

  grphic3d *im1, *im2, *imres;



  im1=ptr_img_3d(im_1);

  im2=ptr_img_3d(im_2);



  imres=ptr_img_3d(im_res);



  imx_reconstruction_geod_eros_3d_p(im1,im2,imres,answer_nr);



  return(1);

}





/*******************************************

** --  imx_reconstruction_geod_eros_3d_p() ----------------------------

**   reconstruction morphologique d'une image

**   a partir d'une autre image par erosion

**   par Nicodeme PAUL Fevrier 2000

**

********************************************/

int     imx_reconstruction_geod_eros_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres, int numas)

{

  int i, j, k;

  int taille_mas=0, taille=0;        /* taille du masque (cube:26, plans:18, croix:6) */

  int ii, bord=1;

  int width,height,depth;

  char reponse='n';

  XPoints_3d *mas;

  grphic3d *imtemp;

  long min;



  width=im1->width;

  height=im1->height;

  depth=im1->depth;



  /* allocation memoire pour l'image temporaire */

  imtemp=cr_grphic3d(im2);



  /* Copie du contenu de l'image im2 dans imres*/

  imx_copie_3d_p(im2,imres);



  /* Copie du contenu de l'image im1 dans imtemp*/

  imx_copie_3d_p(im2,imtemp);



  /* allocation de l'espace memoire pour le masque (mas) d'apres numas(type de masque) */

  mas= masque_3d(numas);



  /* initialisation de la taille du masque */

  switch(numas){

  case 0 : taille_mas= 26; 	/* cas du voisinage cube (26 elements)  */

    break;



  case 1 : taille_mas= 6; 	/* cas du voisinage croix (6 elements)  */

    break;



  case 2 : taille_mas= 18; 	/* cas du voisinage plans (18 elements) */

    break;



  case 3 : taille_mas= 8; 	/* cas du voisinage croix (8 elements)  */

    break;



  case 4 : taille_mas=80; bord=2;

    break;



  case 5 : taille_mas=13;

    break;



  case 6 : taille_mas=3;

    break;



  case 7 : taille_mas=9;

    break;



  case 8 : taille_mas=4;

    break;



  case 9 : taille_mas=40; bord=2;

    break;

  }



  taille=taille_mas/2;



  /* Repeter jusqu'a stabilite : imtemp=imres*/

  do {



    /************* Forward pass **************/



    printf("------ Forward pass ------\n");





    for(k=bord;k<depth-bord;k++)

      for(j=bord;j<height-bord;j++)

        for(i=bord;i<width-bord;i++)

	  {

	    min=imres->mri[i][j][k];

            for(ii=0;ii<taille;ii++)

	      min=MINI(imres->mri[i+mas[ii].x][j+mas[ii].y][k+mas[ii].z],min);

	    imres->mri[i][j][k]=MAXI(min,im1->mri[i][j][k]);

	  }



    /************* Backward pass **************/



    printf("------ Backward pass ------\n");



    for(k=depth-bord-1;k>bord-1;k--)

      for(j=height-bord-1;j>bord-1;j--)

        for(i=width-bord-1;i>bord-1;i--)

	  {

	    min=imres->mri[i][j][k];

            for(ii=taille;ii<taille_mas;ii++)

	      min=MINI(imres->mri[i+mas[ii].x][j+mas[ii].y][k+mas[ii].z],min);

	    imres->mri[i][j][k]=MAXI(min,im1->mri[i][j][k]);

	  }

    /* verification de imres-imtemp */

    for (k=0; k< depth; k++)

      for (j=0; j<height; j++)

	for (i=0; i<width; i++) {

	  if (imres->mri[i][j][k] != imtemp->mri[i][j][k]) {

	    reponse='y';

	    /* Copie du contenu de l'image imres dans imtemp*/

	    imx_copie_3d_p(imres,imtemp);

	    break;

	  }

	  else

	    {

	      reponse='n';

	    }

	}

  }

  while (reponse == 'y');



  /* liberation de la memoire utilisee par les images temporaires*/

  free_grphic3d (imtemp);



  /* sauvegarde de l'image resultat */

  imx_inimaxminpixel_3d_p(imres);



  /* liberation de la memoire utilisee par le masque */

  free(mas);



  return(1);

}





/*******************************************

** --  sharpen_morpho_3d() ----------------------------

**

**  procedure ameliorant les contours d'une image 3D

**        PAUL Nicodeme

**        05-05-2000

********************************************/

void    sharpen_morpho_3d(void)

{

  int im_1,im_res,niter,i,err=0;

  char *quest[12];

  int answer_nr;



  /*     Image   				*/

  im_1=GET_PLACE3D(TEXT0164);

  im_res=GET_PLACE3D(TEXT0006);

  niter= GET_INT(TEXT0007, 2, &err);



  /* Question concerning type of masque  */

  if ((niter>0) && (!err))

    {

      for (i=0;i<12;i++)

	quest[i]=CALLOC(80,char);

      strcpy(quest[0],TEXT0159);

      strcpy(quest[1],TEXT0160);

      strcpy(quest[2],TEXT0191);

      strcpy(quest[3],TEXT0238);

      strcpy(quest[4],TEXT0474);

      strcpy(quest[5],TEXT0475);

      strcpy(quest[6],TEXT0476);

      strcpy(quest[7],TEXT0477);

      strcpy(quest[8],TEXT0478);

      strcpy(quest[9],TEXT0479);

      strcpy(quest[10],TEXT0481);

      strcpy(quest[11],"\0");

      /* Put a question box at position 150,200. */

      answer_nr=GETV_QCM(TEXT0043,(char **)quest);

      for (i=0;i<12;i++)

	FREE(quest[i]);



      if (answer_nr != 10) {

	imx_sharpen_morpho_3d(im_1,im_res,answer_nr,niter);

	show_picture_3d(im_res);

      }

      else

	PUT_MESG("Annulation du traitement");

    }

  else

    PUT_MESG("Erreur (niter)");



  return;

}



/*******************************************

** --  imx_sharpen_morpho_3d() ----------------------------

**

**     Amelioration les contours de l'image

**     par la  morphologie mathematique

********************************************/

int     imx_sharpen_morpho_3d(int im_1, int im_res, int numas, int niter)

{

  grphic3d *im1,*imres;



  im1=ptr_img_3d(im_1);

  imres=ptr_img_3d(im_res);



  imx_sharpen_morpho_3d_p(im1,imres,numas,niter);

  return(1);

}



/*******************************************

** --  imx_sharpen_morpho_3d_p() ----------------------------

**

**    Amelioration les contours de l'image

**          PAUL Nicodeme

**          05-05-2000

********************************************/

int     imx_sharpen_morpho_3d_p(grphic3d *im1, grphic3d *imres, int numas, int niter)

{



  int i;                 /* indice de deplacement dans le masque */

  int imas,jmas,kmas; /* indices de deplacement du masque dans l'image */

  int taille_mas=0;        /* taille du masque (cube:26, plans:18, croix:6) */

  int ii, bord=1;

  XPoints_3d *mas;

  grphic3d *imtemp;

  long max, min;

  int w,h,d;



  /* creation de l'image temporaire imtemp */

  /* allocation memoire pour l'image temporaire */

  imtemp=cr_grphic3d(im1);



  /* Copie du contenu de l'image im1 dans imtemp*/

  imx_copie_3d_p(im1,imtemp);



  /* Copie des parametres de im1 (nom,tailles ...) dans l'image resultat */

  imx_copie_param_3d_p(im1,imres);



  w = imtemp->width;

  h = imtemp->height;

  d = imtemp->depth;



  /* initialisation de l'image resultat a 0 */

  for(imas=0; imas<(w); imas++)

    for(jmas=0; jmas<(h); jmas++)

      for(kmas=0; kmas<(d); kmas++)

	{imres->mri[imas][jmas][kmas]=0;}



  /* allocation de l'espace memoire pour le masque (mas) d'apres numas(type de masque) */

  mas= masque_3d(numas);



  /* initialisation de la taille du masque */

  switch(numas){

  case 0 : taille_mas= 26; 	/* cas du voisinage cube (26 elements)  */

    break;



  case 1 : taille_mas= 6; 	/* cas du voisinage croix (6 elements)  */

    break;



  case 2 : taille_mas= 18; 	/* cas du voisinage plans (18 elements) */

    break;



  case 3 : taille_mas= 8; 	/* cas du voisinage croix (8 elements)  */

    break;



  case 4 : taille_mas=80; bord=2;

    break;



  case 5 : taille_mas=13;

    break;



  case 6 : taille_mas=3;

    break;



  case 7 : taille_mas=9;

    break;



  case 8 : taille_mas=4;

    break;



  case 9 : taille_mas=40; bord=2;

    break;

  }





  /* Traitement */

  for (i=0;i<niter;i++)

    {

      for(imas=bord;imas<(w-bord);imas++)

        for(jmas=bord;jmas<(h-bord);jmas++)

          for(kmas=bord;kmas<(d-bord);kmas++)

            {

	      max=min=imtemp->mri[imas][jmas][kmas];

	      for(ii=0;ii<=taille_mas;ii++)

		{

		  max=MAXI(max,imtemp->mri[imas+mas[ii].x][jmas+mas[ii].y][kmas+mas[ii].z]);

		  min=MINI(min,imtemp->mri[imas+mas[ii].x][jmas+mas[ii].y][kmas+mas[ii].z]);

		}

	      if ((max-imtemp->mri[imas][jmas][kmas])<(imtemp->mri[imas][jmas][kmas]-min))

		imres->mri[imas][jmas][kmas]=(TYPEMRI3D)max;

	      else

		if ((max-imtemp->mri[imas][jmas][kmas])>(imtemp->mri[imas][jmas][kmas]-min))

		  imres->mri[imas][jmas][kmas]=(TYPEMRI3D)min;

		else

		  imres->mri[imas][jmas][kmas]=imtemp->mri[imas][jmas][kmas];

	    }

      imx_copie_3d_p(imres,imtemp);

    }





  /* liberation de la memoire utilisee par les images temporaires*/

  free_grphic3d (imtemp);



  /* sauvegarde de l'image dilatee */

  imx_inimaxminpixel_3d_p(imres);



  /* liberation de la memoire utilisee par le masque */

  free(mas);



  return(1);

}





/*******************************************

** --  contrast_morpho_3d() ----------------------------

**

**  procedure ameliorant le contraste d'une image 3D

**        PAUL Nicodeme

**        05-05-2000

********************************************/

void    contrast_morpho_3d(void)

{

  int im_1,im_res,niter,i,err=0;

  char *quest[12];

  int answer_nr;



  /*     Image   				*/

  im_1=GET_PLACE3D(TEXT0164);

  im_res=GET_PLACE3D(TEXT0006);

  niter= GET_INT(TEXT0007, 2, &err);



  /* Question concerning type of masque  */

  if ((niter>0) && (!err))

    {

      for (i=0;i<12;i++)

	quest[i]=CALLOC(80,char);

      strcpy(quest[0],TEXT0159);

      strcpy(quest[1],TEXT0160);

      strcpy(quest[2],TEXT0191);

      strcpy(quest[3],TEXT0238);

      strcpy(quest[4],TEXT0474);

      strcpy(quest[5],TEXT0475);

      strcpy(quest[6],TEXT0476);

      strcpy(quest[7],TEXT0477);

      strcpy(quest[8],TEXT0478);

      strcpy(quest[9],TEXT0479);

      strcpy(quest[10],TEXT0481);

      strcpy(quest[11],"\0");

      /* Put a question box at position 150,200. */

      answer_nr=GETV_QCM(TEXT0043,(char **)quest);

      for (i=0;i<12;i++)

	FREE(quest[i]);



      if (answer_nr != 10) {

	imx_contrast_morpho_3d(im_1,im_res,answer_nr,niter);

	show_picture_3d(im_res);

      }

      else

   	PUT_MESG("Annulation du traitement");

    }

  else

    PUT_MESG("Erreur (niter)");



  return;

}



/*******************************************

** --  imx_contrast_morpho_3d() ----------------------------

**

**     Amelioration du contraste de l'image

**     par la  morphologie mathematique

********************************************/

int     imx_contrast_morpho_3d(int im_1, int im_res, int numas, int niter)

{

  grphic3d *im1,*imres;



  im1=ptr_img_3d(im_1);

  imres=ptr_img_3d(im_res);



  imx_contrast_morpho_3d_p(im1,imres,numas,niter);

  return(1);

}



/*******************************************

** --  imx_sharpen_morpho_3d_p() ----------------------------

**

**    Amelioration du contraste de l'image

**          PAUL Nicodeme

**          05-05-2000

********************************************/

int     imx_contrast_morpho_3d_p(grphic3d *im1, grphic3d *imres, int numas, int niter)

{



  int i,j,k; /* indices de deplacement du masque dans l'image */

  int taille_mas=0;        /* taille du masque (cube:26, plans:18, croix:6) */

  int ii, bord=1;

  XPoints_3d *mas;

  grphic3d *imtemp, *imgrad;

  double max, min;

  double sum, sum_pond, c1, c2, quotient;

  double *data;

  int w,h,d;



  /* creation de l'image temporaire imtemp */

  /* allocation memoire pour l'image temporaire */

  imtemp=cr_grphic3d(im1);



  /* Copie du contenu de l'image im1 dans imtemp*/

  imx_copie_3d_p(im1,imtemp);



  /* Copie des parametres de im1 (nom,tailles ...) dans l'image resultat */

  imx_copie_param_3d_p(im1,imres);



  /* allocation memoire pour l'image temporaire */

  imgrad=cr_grphic3d(im1);



  w = imtemp->width;

  h = imtemp->height;

  d = imtemp->depth;



  /* initialisation de l'image resultat a 0 */

  for(i=0; i<w; i++)

    for(j=0; j<h; j++)

      for(k=0; k<d; k++)

	{imres->mri[i][j][k]=0;}



  /* allocation de l'espace memoire pour le masque (mas) d'apres numas(type de masque) */

  mas= masque_3d(numas);



  /* initialisation de la taille du masque */

  switch(numas){

  case 0 : taille_mas= 26; 	/* cas du voisinage cube (26 elements)  */

    break;



  case 1 : taille_mas= 6; 	/* cas du voisinage croix (6 elements)  */

    break;



  case 2 : taille_mas= 18; 	/* cas du voisinage plans (18 elements) */

    break;



  case 3 : taille_mas= 8; 	/* cas du voisinage croix (8 elements)  */

    break;



  case 4 : taille_mas=80; bord=2;

    break;



  case 5 : taille_mas=13;

    break;



  case 6 : taille_mas=3;

    break;



  case 7 : taille_mas=9;

    break;



  case 8 : taille_mas=4;

    break;



  case 9 : taille_mas=40; bord=2;

    break;

  }



  /* Calcul du gradient de l'image */

  imx_GradM_3d_p(im1,imgrad,numas,0,1);



  data = ( double *) calloc

    ((im1->width-2*bord+1)*(im1->height-2*bord+1)*(im1->depth-2*bord+1),sizeof(double));



  w = im1->width;

  h = im1->height;

  d = im1->depth;

  /* Traitement */

  for (i=0;i<niter;i++)

    {

      for(i=bord;i<(w-bord);i++)

        for(j=bord;j<(h-bord);j++)

          for(k=bord;k<(d-bord);k++)

            {

	      sum_pond=0; sum=0;

	      for(ii=0;ii<=taille_mas;ii++)

		{

		  sum += (imgrad->rcoeff*

			  imgrad->mri[i+mas[ii].x][j+mas[ii].y][k+mas[ii].z]);

		  sum_pond += (imtemp->rcoeff * imgrad->rcoeff *

			       imtemp->mri[i+mas[ii].x][j+mas[ii].y][k+mas[ii].z]*

			       imgrad->mri[i+mas[ii].x][j+mas[ii].y][k+mas[ii].z]);

		}

	      if (sum != 0.0)

		quotient=(double) (sum_pond/sum);

	      else

		quotient=0.0;

	      if (imtemp->rcoeff*imtemp->mri[i][j][k] > quotient) {

		c1=(double) ((imtemp->rcoeff * imtemp->mri[i][j][k])-quotient)/

		  ((imtemp->rcoeff * imtemp->mri[i][j][k])+quotient);

		c2=sqrt(c1);

		data[i*(im1->height-2*bord)*(im1->depth-2*bord)+j*(im1->depth-2*bord)+k-bord]=(quotient*(1+c2))/(1-c2);

	      }

	      else {

		if ((imtemp->rcoeff * imtemp->mri[i][j][k])+quotient == 0.0)

		  c1=0.0;

		else {

		  c1=(double)	(quotient-(imtemp->rcoeff * imtemp->mri[i][j][k]))/

		    ((imtemp->rcoeff * imtemp->mri[i][j][k])+quotient);

		}

		c2=sqrt(c1);

		data[i*(im1->height-2*bord)*(im1->depth-2*bord)+j*(im1->depth-2*bord)+k-bord]=(quotient*(1-c2))/(1+c2);

	      }

	    }

      max=0.0; min=10000.0;

      for(i=0;i<(w-2*bord);i++)

        for(j=0;j<(h-2*bord);j++)

          for(k=0;k<(d-2*bord);k++)	{

	    max=MAXI(max,data[i*(im1->height-2*bord)*(im1->depth-2*bord)+j*(im1->depth-2*bord)+k]);

	    min=MINI(min,data[i*(im1->height-2*bord)*(im1->depth-2*bord)+j*(im1->depth-2*bord)+k]);

	  }



      /* calcul de rcoeff */

      imx_brukermax_3d((float)max,(float)min,imres);



      w = imres->width;

      h = imres->height;

      d = imres->depth;

      /* calcul de l'image resultat */

      for (i = bord; i < w-bord; i++)

	for (j = bord; j < h-bord; j++)

	  for (k = bord; k < d-bord; k++)  {

	    imres->mri[i][j][k]=(TYPEMRI3D)

	      floor(data[i*(im1->height-2*bord)*(im1->depth-2*bord)+j*(im1->depth-2*bord)+k-bord]/(imres->rcoeff));

	  }

      imx_copie_3d_p(imres,imtemp);

    }





  /* liberation de la memoire utilisee par les images temporaires*/

  free_grphic3d (imtemp);



  /* liberation de la memoire utilisee par les images temporaires*/

  free_grphic3d (imgrad);



  free (data);



  /* sauvegarde de l'image dilatee */

  imx_inimaxminpixel_3d_p(imres);



  /* liberation de la memoire utilisee par le masque */

  free(mas);



  return(1);

}







/*******************************************

** --  eros_bin_3d() ----------------------------

**

**  procedure erosion binaire 3D

********************************************/

int     eros_bin_3d(void)

{

  int im_1,im_res,niter,i,err=0;

  char *quest[3];

  int answer_nr;



  /*     Image   				*/

  im_1=GET_PLACE3D(TEXT0011);

  im_res=GET_PLACE3D(TEXT0006);

  niter= GET_INT(TEXT0007, 2, &err);



  /* Question concerning type of masque  */

  if (!err)

    {

      for (i=0;i<3;i++)

	quest[i]=CALLOC(80,char);

      strcpy(quest[0],TEXT0159);

      strcpy(quest[1],TEXT0160);

      strcpy(quest[2],"\0");

      /* Put a question box at position 150,200. */

      answer_nr=GETV_QCM(TEXT0043,(char **)quest);

      for (i=0;i<3;i++)

	FREE(quest[i]);



      /*protoize:???*/

      imx_eros_bin_3d(im_1,im_res,answer_nr,0/*???*/);

      show_picture_3d(im_res);

    }

  return(1);

}



/*******************************************

** --  imx_eros_bin_3d() ----------------------------

**

** Erosion binaire d'une image

**

********************************************/

int     imx_eros_bin_3d(int im_1, int im_res, int numas, int niter)

{

  grphic3d *im1,*imres;



  im1=ptr_img_3d(im_1);

  imres=ptr_img_3d(im_res);



  imx_eros_bin_3d_p(im1,imres,numas);

  return(1);

}





/*******************************************

** --  imx_eros_bin_3d_p() ----------------------------

**

** Erosion binaire d'une image (equivaut a enlever a une image

** le contour

**

********************************************/

int     imx_eros_bin_3d_p(grphic3d *im1, grphic3d *imres, int numas)

{

  int i,j,k,l;

  long resul;

  grphic3d *imtemp;

  int wdth,hght,dpth,wdthm1,hghtm1,dpthm1;

  XPoints_3d tabvois[8];





  /* creation de l'image temporaire imtemp */

  imtemp=cr_grphic3d(im1);



  imx_copie_param_3d_p(im1,imres);



  wdth=im1->width; wdthm1=wdth-1;

  hght=im1->height;hghtm1=hght-1;

  dpth=im1->depth; dpthm1=dpth-1;



  /*  masque_3d(numas,mas); */



  /*Definition du tableau tabvois dans le cas de la croix 3d*/

  tabvois[0].x=-1;tabvois[0].y=0;tabvois[0].z=0;

  tabvois[1].x=1;tabvois[1].y=0;tabvois[1].z=0;

  tabvois[2].x=0;tabvois[2].y=-1;tabvois[2].z=0;

  tabvois[3].x=0;tabvois[3].y=1;tabvois[3].z=0;

  tabvois[4].x=0;tabvois[4].y=0;tabvois[4].z=-1;

  tabvois[5].x=0;tabvois[5].y=0;tabvois[5].z=1;







  for (i=1;i<wdthm1;i++)

    for (j=1;j<hghtm1;j++)

      for (k=1;k<dpthm1;k++)

        {

	  if (im1->mri[i][j][k]!=0) {

	    resul=im1->mri[i][j][k];

	    for (l=0;l<6;l++)

	      resul=resul && im1->mri[i+tabvois[l].x][j+tabvois[l].y][k+tabvois[l].z];

	    if (resul) imtemp->mri[i][j][k]=im1->mri[i][j][k];

	    else imtemp->mri[i][j][k]=0;

          }

        }

  imx_copie_3d_p(imtemp,imres);



  /* liberation de la memoire utilisee par les images temporaires*/

  free_grphic3d (imtemp);



  imx_inimaxminpixel_3d_p(imres);



  return(1);

}





/*******************************************

** --  FAS_3d() ----------------------------

**

**  procedure filtre alterne 3D

**   par Fahima DJAOUI 99

**

********************************************/

void    FAS_3d(void)

{

  int im_1,im_res,niter,i,err=0;

  /*

  int im_2,im_res2;

  */

  char *quest[5];

  int answer_nr;



  /*     Image   				*/

  im_1=GET_PLACE3D(TEXT0005);

  im_res=GET_PLACE3D(TEXT0006);





  niter= GET_INT(TEXT0007, 2, &err);



  /* Question concerning type of masque  */

  if (!err)

    {

      for (i=0;i<5;i++)

	quest[i]=CALLOC(80,char);

      strcpy(quest[0],TEXT0159);

      strcpy(quest[1],TEXT0160);

      strcpy(quest[2],TEXT0191);

      strcpy(quest[3],TEXT0238);

      strcpy(quest[4],"\0");

      /* Put a question box at position 150,200. */

      answer_nr=GETV_QCM(TEXT0043,(char **)quest);

      for (i=0;i<5;i++)

	FREE(quest[i]);



      imx_open_3d(im_1,im_res,answer_nr,niter);

      imx_close_3d(im_res,im_1,answer_nr,niter);



      imx_copie_param_3d(im_1,im_res);



      show_picture_3d(im_res);

    }

  return;

}



/*******************************************

** --  FASS_3d() ----------------------------

**

**  procedure filtre alterne sequentiel 3D

**   par Fahima DJAOUI 99

**

********************************************/

void    FASS_3d(void)

{

  int im_1,im_res,niter,i,err=0;



  char *quest[5];

  int answer_nr;



  /*     Image   				*/

  im_1=GET_PLACE3D(TEXT0005);

  im_res=GET_PLACE3D(TEXT0006);





  niter= GET_INT(TEXT0007, 2, &err);



  /* Question concerning type of masque  */

  if (!err)

    {

      for (i=0;i<5;i++)

	quest[i]=CALLOC(80,char);

      strcpy(quest[0],TEXT0159);

      strcpy(quest[1],TEXT0160);

      strcpy(quest[2],TEXT0191);

      strcpy(quest[3],TEXT0238);

      strcpy(quest[4],"\0");

      /* Put a question box at position 150,200. */

      answer_nr=GETV_QCM(TEXT0043,(char **)quest);

      for (i=0;i<5;i++)

	FREE(quest[i]);





      for (i=1;i<= niter;i++)

	{

	  imx_open_3d(im_1,im_res,answer_nr,i);



	  imx_close_3d(im_res,im_1,answer_nr,i);

	}



      imx_copie_param_3d(im_1,im_res);

      show_picture_3d(im_res);



    }

  return;

}





/*******************************************

** --  masque_3d() ----------------------------

**

**   Creation d'un masque pour l'erosion, la dilatation

**   numas : numero du masque

**   mas  masque proprement dir

**

**   Generation du  tableau "mas"

**   reecrite par Cyril Feruglio

**   derniere maj le :  16-04-1998

********************************************/

XPoints_3d  *masque_3d(int numas)

{

  XPoints_3d *mas=NULL;



  switch (numas)

    {



    case 0 : 					/* cas du voisinage cube */

      mas=CALLOC(26,XPoints_3d);	/* 26 elements dans le voisinage */

      mas[0].x= -1; mas[0].y= -1; mas[0].z= -1;

      mas[1].x=-1; mas[1].y= -1; mas[1].z= 0;

      mas[2].x= -1; mas[2].y=-1; mas[2].z= 1;

      mas[3].x= 0; mas[3].y= -1; mas[3].z=-1;

      mas[4].x=0; mas[4].y=-1; mas[4].z= 0;

      mas[5].x= 0; mas[5].y=-1; mas[5].z=1;

      mas[6].x=1; mas[6].y= -1; mas[6].z=-1;

      mas[7].x=1; mas[7].y=-1; mas[7].z=0;

      mas[8].x= 1; mas[8].y=-1; mas[8].z= 1;

      mas[9].x= -1; mas[9].y= 0; mas[9].z= -1;

      mas[10].x= -1; mas[10].y= 0; mas[10].z= 0;

      mas[11].x= -1; mas[11].y=0; mas[11].z= 1;

      mas[12].x= 0; mas[12].y= 0; mas[12].z=-1;

      mas[13].x= 0; mas[13].y= 0; mas[13].z=1;

      mas[14].x=1; mas[14].y= 0; mas[14].z= -1;

      mas[15].x= 1; mas[15].y=0; mas[15].z= 0;

      mas[16].x=1; mas[16].y=0; mas[16].z= 1;

      mas[17].x= -1; mas[17].y= 1; mas[17].z= -1;

      mas[18].x= -1; mas[18].y= 1; mas[18].z=0;

      mas[19].x= -1; mas[19].y= 1; mas[19].z= 1;

      mas[20].x= 0; mas[20].y=1; mas[20].z= -1;

      mas[21].x= 0; mas[21].y= 1; mas[21].z= 0;

      mas[22].x=0; mas[22].y= 1; mas[22].z= 1;

      mas[23].x= 1; mas[23].y=1; mas[23].z=-1;

      mas[24].x=1; mas[24].y= 1; mas[24].z=0;

      mas[25].x=1; mas[25].y=1; mas[25].z= 1;

      break;



    case 1 :				/* cas du voisinage croix */

      mas=CALLOC(6, XPoints_3d); 	/* 6 elements dans le voisinage */

      mas[0].x= 0; mas[0].y= -1; mas[0].z= 0;

      mas[1].x= -1; mas[1].y= 0; mas[1].z=0;

      mas[2].x= 0; mas[2].y= 0; mas[2].z= -1;

      mas[3].x= 0; mas[3].y=0; mas[3].z= 1;

      mas[4].x= 1; mas[4].y= 0; mas[4].z= 0;

      mas[5].x=0; mas[5].y= 1; mas[5].z= 0;

      break;



    case 2 :					/* cas du voisinage plans */

      mas=CALLOC(18, XPoints_3d);  /* 18 elements dans le voisinage */

      mas[0].x=-1; mas[0].y=-1; mas[0].z= 0;

      mas[1].x= 0; mas[1].y= -1; mas[1].z= -1;

      mas[2].x= 0; mas[2].y= -1; mas[2].z= 0;

      mas[3].x= 0; mas[3].y= -1; mas[3].z= 1;

      mas[4].x= 1; mas[4].y=-1; mas[4].z= 0;

      mas[5].x= -1; mas[5].y= 0; mas[5].z=-1;

      mas[6].x= -1; mas[6].y= 0; mas[6].z=0;

      mas[7].x=-1; mas[7].y= 0; mas[7].z= 1;

      mas[8].x= 0; mas[8].y=0; mas[8].z= -1;

      mas[9].x= 0; mas[9].y= 0; mas[9].z= 1;

      mas[10].x= 1; mas[10].y= 0; mas[10].z=-1;

      mas[11].x= 1; mas[11].y= 0; mas[11].z= 0;

      mas[12].x= 1; mas[12].y=0; mas[12].z= 1;

      mas[13].x= -1; mas[13].y= 1; mas[13].z= 0;

      mas[14].x=0; mas[14].y= 1; mas[14].z= -1;

      mas[15].x= 0; mas[15].y=1; mas[15].z=0;

      mas[16].x=0; mas[16].y= 1; mas[16].z=1;

      mas[17].x=1; mas[18].y= 1; mas[18].z= 0;



      break;



    case 3 :				/* cas du voisinage coin */

      mas=CALLOC(8, XPoints_3d);  /* 8 elements dans le voisinage */

      mas[0].x=-1; mas[0].y=-1; mas[0].z= -1;

      mas[1].x=-1; mas[1].y=-1; mas[1].z= 1;

      mas[2].x=1; mas[2].y=-1; mas[2].z=-1;

      mas[3].x=1; mas[3].y=-1; mas[3].z=1;

      mas[4].x=-1; mas[4].y=1; mas[4].z=-1;

      mas[5].x=-1; mas[5].y=1; mas[5].z=1;

      mas[6].x=1; mas[6].y=1; mas[6].z=-1;

      mas[7].x=1; mas[7].y=1; mas[7].z=1;



      break;



    case 4 :

      mas=CALLOC(80, XPoints_3d);

      mas[0].x=-1; mas[0].y=-2; mas[0].z= -1;

      mas[1].x=-1; mas[1].y=-2; mas[1].z= 0;

      mas[2].x=-1; mas[2].y=-2; mas[2].z= 1;

      mas[3].x=0; mas[3].y=-2; mas[3].z= -1;

      mas[4].x=0; mas[4].y=-2; mas[4].z= 0;

      mas[5].x=0; mas[5].y=-2; mas[5].z= 1;

      mas[6].x=1; mas[6].y=-2; mas[6].z= -1;

      mas[7].x=1; mas[7].y=-2; mas[7].z= 0;

      mas[8].x=1; mas[8].y=-2; mas[8].z= 1;

      mas[9].x=-2; mas[9].y=-1; mas[9].z= -1;

      mas[10].x=-2; mas[10].y=-1; mas[10].z= 0;

      mas[11].x=-2; mas[11].y=-1; mas[11].z= 1;

      mas[12].x=-1; mas[12].y=-1; mas[12].z= -2;

      mas[13].x=-1; mas[13].y=-1; mas[13].z= -1;

      mas[14].x=-1; mas[14].y=-1; mas[14].z= 0;

      mas[15].x=-1; mas[15].y=-1; mas[15].z= 1;

      mas[16].x=-1; mas[16].y=-1; mas[16].z= 2;

      mas[17].x=0; mas[17].y=-1; mas[17].z= -2;

      mas[18].x=0; mas[18].y=-1; mas[18].z= -1;

      mas[19].x=0; mas[19].y=-1; mas[19].z= 0;

      mas[20].x=0; mas[20].y=-1; mas[20].z= 1;

      mas[21].x=0; mas[21].y=-1; mas[21].z= 2;

      mas[22].x=1; mas[22].y=-1; mas[22].z= -2;

      mas[23].x=1; mas[23].y=-1; mas[23].z= -1;

      mas[24].x=1; mas[24].y=-1; mas[24].z= 0;

      mas[25].x=1; mas[25].y=-1; mas[25].z= 1;

      mas[26].x=1; mas[26].y=-1; mas[26].z= 2;

      mas[27].x=2; mas[27].y=-1; mas[27].z= -1;

      mas[28].x=2; mas[28].y=-1; mas[28].z= 0;

      mas[29].x=2; mas[29].y=-1; mas[29].z= 1;

      mas[30].x=-2; mas[30].y=0; mas[30].z= -1;

      mas[31].x=-2; mas[31].y=0; mas[31].z= 0;

      mas[32].x=-2; mas[32].y=0; mas[32].z= 1;

      mas[33].x=-1; mas[33].y=0; mas[33].z= -2;

      mas[34].x=-1; mas[34].y=0; mas[34].z= -1;

      mas[35].x=-1; mas[35].y=0; mas[35].z= 0;

      mas[36].x=-1; mas[36].y=0; mas[36].z= 1;

      mas[37].x=-1; mas[37].y=0; mas[37].z= 2;

      mas[38].x=0; mas[38].y=0; mas[38].z= -2;

      mas[39].x=0; mas[39].y=0; mas[39].z= -1;

      mas[40].x=0; mas[40].y=0; mas[40].z= 1;

      mas[41].x=0; mas[41].y=0; mas[41].z= 2;

      mas[42].x=1; mas[42].y=0; mas[42].z= -2;

      mas[43].x=1; mas[43].y=0; mas[43].z= -1;

      mas[44].x=1; mas[44].y=0; mas[44].z= 0;

      mas[45].x=1; mas[45].y=0; mas[45].z= 1;

      mas[46].x=1; mas[46].y=0; mas[46].z= 2;

      mas[47].x=2; mas[47].y=0; mas[47].z= -1;

      mas[48].x=2; mas[48].y=0; mas[48].z= 0;

      mas[49].x=2; mas[49].y=0; mas[49].z= 1;

      mas[50].x=-2; mas[50].y=1; mas[50].z= -1;

      mas[51].x=-2; mas[51].y=1; mas[51].z= 0;

      mas[52].x=-2; mas[52].y=1; mas[52].z= 1;

      mas[53].x=-1; mas[53].y=1; mas[53].z= -2;

      mas[54].x=-1; mas[54].y=1; mas[54].z= -1;

      mas[55].x=-1; mas[55].y=1; mas[55].z= 0;

      mas[56].x=-1; mas[56].y=1; mas[56].z= 1;

      mas[57].x=-1; mas[57].y=1; mas[57].z= 2;

      mas[58].x=0; mas[58].y=1; mas[58].z= -2;

      mas[59].x=0; mas[59].y=1; mas[59].z= -1;

      mas[60].x=0; mas[60].y=1; mas[60].z= 0;

      mas[61].x=0; mas[61].y=1; mas[61].z= 1;

      mas[62].x=0; mas[62].y=1; mas[62].z= 2;

      mas[63].x=1; mas[63].y=1; mas[63].z= -2;

      mas[64].x=1; mas[64].y=1; mas[64].z= -1;

      mas[65].x=1; mas[65].y=1; mas[65].z= 0;

      mas[66].x=1; mas[66].y=1; mas[66].z= 1;

      mas[67].x=1; mas[67].y=1; mas[67].z= 2;

      mas[68].x=2; mas[68].y=1; mas[68].z= -1;

      mas[69].x=2; mas[69].y=1; mas[69].z= 0;

      mas[70].x=2; mas[70].y=1; mas[70].z= 1;

      mas[71].x=-1; mas[71].y=2; mas[71].z= -1;

      mas[72].x=-1; mas[72].y=2; mas[72].z= 0;

      mas[73].x=-1; mas[73].y=2; mas[73].z= 1;

      mas[74].x=0; mas[74].y=2; mas[74].z= -1;

      mas[75].x=0; mas[75].y=2; mas[75].z= 0;

      mas[76].x=0; mas[76].y=2; mas[76].z= 1;

      mas[77].x=1; mas[77].y=2; mas[77].z= -1;

      mas[78].x=1; mas[78].y=2; mas[78].z= 0;

      mas[79].x=1; mas[79].y=2; mas[79].z= 1;

      break;



    case 5 : mas=CALLOC(13,XPoints_3d);

      mas[0].x= -1; mas[0].y= -1; mas[0].z= -1;

      mas[1].x=-1; mas[1].y= -1; mas[1].z= 0;

      mas[2].x= -1; mas[2].y=-1; mas[2].z= 1;

      mas[3].x= 0; mas[3].y= -1; mas[3].z=-1;

      mas[4].x=0; mas[4].y=-1; mas[4].z= 0;

      mas[5].x= 0; mas[5].y=-1; mas[5].z=1;

      mas[6].x=1; mas[6].y= -1; mas[6].z=-1;

      mas[7].x=1; mas[7].y=-1; mas[7].z=0;

      mas[8].x= 1; mas[8].y=-1; mas[8].z= 1;

      mas[9].x= -1; mas[9].y= 0; mas[9].z= -1;

      mas[10].x= -1; mas[10].y= 0; mas[10].z= 0;

      mas[11].x= -1; mas[11].y=0; mas[11].z= 1;

      mas[12].x= 0; mas[12].y= 0; mas[12].z=-1;

      break;



    case 6 :

      mas=CALLOC(3, XPoints_3d);

      mas[0].x= 0; mas[0].y= -1; mas[0].z= 0;

      mas[1].x= -1; mas[1].y= 0; mas[1].z=0;

      mas[2].x= 0; mas[2].y= 0;  mas[2].z= -1;

      break;



    case 7 :

      mas=CALLOC(9, XPoints_3d);

      mas[0].x=-1; mas[0].y=-1; mas[0].z= 0;

      mas[1].x= 0; mas[1].y= -1; mas[1].z= -1;

      mas[2].x= 0; mas[2].y= -1; mas[2].z= 0;

      mas[3].x= 0; mas[3].y= -1; mas[3].z= 1;

      mas[4].x= 1; mas[4].y=-1; mas[4].z= 0;

      mas[5].x= -1; mas[5].y= 0; mas[5].z=-1;

      mas[6].x= -1; mas[6].y= 0; mas[6].z=0;

      mas[7].x=-1; mas[7].y= 0; mas[7].z= 1;

      mas[8].x= 0; mas[8].y=0; mas[8].z= -1;

      break;



    case 8 :

      mas=CALLOC(4, XPoints_3d);

      mas[0].x=-1; mas[0].y=-1; mas[0].z= -1;

      mas[1].x=-1; mas[1].y=-1; mas[1].z= 1;

      mas[2].x=1; mas[2].y=-1; mas[2].z=-1;

      mas[3].x=1; mas[3].y=-1; mas[3].z=1;

      break;



    case 9 :

      mas=CALLOC(40, XPoints_3d);

      mas[0].x=-1; mas[0].y=-2; mas[0].z= -1;

      mas[1].x=-1; mas[1].y=-2; mas[1].z= 0;

      mas[2].x=-1; mas[2].y=-2; mas[2].z= 1;

      mas[3].x=0; mas[3].y=-2; mas[3].z= -1;

      mas[4].x=0; mas[4].y=-2; mas[4].z= 0;

      mas[5].x=0; mas[5].y=-2; mas[5].z= 1;

      mas[6].x=1; mas[6].y=-2; mas[6].z= -1;

      mas[7].x=1; mas[7].y=-2; mas[7].z= 0;

      mas[8].x=1; mas[8].y=-2; mas[8].z= 1;

      mas[9].x=-2; mas[9].y=-1; mas[9].z= -1;

      mas[10].x=-2; mas[10].y=-1; mas[10].z= 0;

      mas[11].x=-2; mas[11].y=-1; mas[11].z= 1;

      mas[12].x=-1; mas[12].y=-1; mas[12].z= -2;

      mas[13].x=-1; mas[13].y=-1; mas[13].z= -1;

      mas[14].x=-1; mas[14].y=-1; mas[14].z= 0;

      mas[15].x=-1; mas[15].y=-1; mas[15].z= 1;

      mas[16].x=-1; mas[16].y=-1; mas[16].z= 2;

      mas[17].x=0; mas[17].y=-1; mas[17].z= -2;

      mas[18].x=0; mas[18].y=-1; mas[18].z= -1;

      mas[19].x=0; mas[19].y=-1; mas[19].z= 0;

      mas[20].x=0; mas[20].y=-1; mas[20].z= 1;

      mas[21].x=0; mas[21].y=-1; mas[21].z= 2;

      mas[22].x=1; mas[22].y=-1; mas[22].z= -2;

      mas[23].x=1; mas[23].y=-1; mas[23].z= -1;

      mas[24].x=1; mas[24].y=-1; mas[24].z= 0;

      mas[25].x=1; mas[25].y=-1; mas[25].z= 1;

      mas[26].x=1; mas[26].y=-1; mas[26].z= 2;

      mas[27].x=2; mas[27].y=-1; mas[27].z= -1;

      mas[28].x=2; mas[28].y=-1; mas[28].z= 0;

      mas[29].x=2; mas[29].y=-1; mas[29].z= 1;

      mas[30].x=-2; mas[30].y=0; mas[30].z= -1;

      mas[31].x=-2; mas[31].y=0; mas[31].z= 0;

      mas[32].x=-2; mas[32].y=0; mas[32].z= 1;

      mas[33].x=-1; mas[33].y=0; mas[33].z= -2;

      mas[34].x=-1; mas[34].y=0; mas[34].z= -1;

      mas[35].x=-1; mas[35].y=0; mas[35].z= 0;

      mas[36].x=-1; mas[36].y=0; mas[36].z= 1;

      mas[37].x=-1; mas[37].y=0; mas[37].z= 2;

      mas[38].x=0; mas[38].y=0; mas[38].z= -2;

      mas[39].x=0; mas[39].y=0; mas[39].z= -1;

      break;



    default:

      PUT_MESG("Annulation du traitement");

      break;



    }



  return(mas);

}

