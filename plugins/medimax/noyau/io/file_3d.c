/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*!***********************************************************************
***    
***    \file:        imx_3d.c
***
***    project:    Imagix 1.01
***            
***
***    \brief description:    Menu de visualisation et traitement 3d
***
***     Remarque : L'utilisation de la 3D necessite au minimum
***                de 128 Moctets de memoire. En deca le programme tourne
***                mais ne peut plus executer des programmes
***                tels que xgraph, impression, ...
***    
***    Copyright (c) 1993, ULP-IPB Strasbourg.
***    All rights are reserved.
***
***     Last user action by: Mr. ARMSPACH on Jun 23th 1995
***     Last user action by: Mr. ARNOLD on Mar 1996
***     Last user action by: Mr. G.ENGER on July 1996
***     Last user action by: Mr. GOUVERNEUR on Nov 24th 1998
***     Last user action by: Mr. FREYMANN on Dec 2007
***
**************************************************************************/

#include <config.h>
#include <stdio.h> 
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <limits.h>
#ifdef WIN32
# include <io.h>  // pour _sleep
#endif
#include <string.h>

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_3d.h"
#include "noyau/gui/imx_picture3d.h"
#include "noyau/imx_lang.h"
#include "noyau/io/file_3d.h"
#include "noyau/io/imx_file.h"
#include "math/imx_matrix.h"
#include "math/imx_math.h"
#include "outils/imx_misc.h"
#include "noyau/mani_3d.h"
#include "recalage/reca_3d.h"
#include "talairach/tala_3d.h"
#ifdef __GTK_GUI
#include "noyau/io/imx__atm.h"
#include "noyau/zone_3d.h"
#include "noyau/gui/curve_3d.h"
#include "math/oper_3d.h"
#include "math/operation.h"
#include "noyau/io/imx_export_file.h"
#include "noyau/io/lect_analyze_3d.h"

#ifndef COMPILE_FOR_MEDIPY
#include "noyau/annotations/annotations_export.h"
#endif

#ifdef HAS_CORRECT_BOOST
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "noyau/io/io_xml.h"

#include "noyau/propertymanager.h"
#endif


#endif /* __GTK_GUI */


#ifdef HAS_NIFTI
#include "noyau/io/file_nifti.h"
#include "nifti1_io.h"
#endif

/*   Declaration des fonctions externes */
#ifndef __DOXYGEN_SHOULD_SKIP_THIS__
extern int IMXColor_BindColormap_3d_ptr(ptr_grphic3d ptrPic, int nNewCMap);
extern int access(const char *pathname, int mode);


/** Fonction declared in this file:
*** ------------------------------ */
grphic3d    *load_mri_bruker_3d (const char *file, grphic3d *graphic, int dpth, int numfirst);
float       pix_value_3d_p      (grphic3d *im1, float *X, int TYPE);
int 		find_whd_from_orient_p(grphic3d *img, int wdth, int hght, int dpth, int read_format_3d);
int 		fill_mri_32b		(grphic3d *img,long *buff32b,int i, int j, int k, int pos,  float coeff, int read_format, char is_filebigEndian);
int 		fill_mri_16b		(grphic3d *img,short *buff16b,int i, int j, int k, int pos, float coeff, int read_format,  char is_filebigEndian);
int 		fill_mri_8b			(grphic3d *img,unsigned char *buff8b,int i, int j, int k, int pos, float coeff, int read_format);
/* ATTENTION : Lors de la recuperation des modifications apportees
   il faut proceder autrement pour eviter la redeclaration de fonctions */

int         save_mri_ipb_3d     (char *file, int start, int num_images, int step);
int 	    save_pictures_3d	(char *file, int start, int num_images, int step);
int         load_mri_ipb_3d     (const char *file, grphic3d *image, int number);

long getPosFile(const char *file, int numfirst, int filewdth, int filehght, int filedpth, int bytesppixel, int format);
#endif /*__DOXYGEN_SHOULD_SKIP_THIS__*/



#ifdef __GTK_GUI
/* --load_picture_3d() ------------------------------------------------
**
**    Chargement d'une image 3d
**    pos : numero image 3d
**    file : nom du fichier
**    max : nombre de coupes
**    step : pas pour la recherche dans fichiers multiples
**    numfirst : numero de la premiere image
**
**    max et step ne sont actuellement utilise que pour les format 
**    BRUKER_FORMAT et SIEMENS_FORMAT
-----------------------------------------------------------------*/
int    load_picture_3d_cst(int pos, const char *file, int max, int step, int numfirst)
{
  grphic3d *graphic;
  int ismri;
  int temp1,temp2,zoom,old_zoom;
       
  ismri = is_mri(file,NULL);

  if (!(ismri == INTERFILE_FORMAT || ismri == BRUKER_FORMAT || ismri == SIEMENS_FORMAT	    ||ismri == PARAVISION_FORMAT
        || ismri == IPB_FORMAT || ismri == VANDERBILT_FORMAT 
	|| ismri == IPB_3D_FORMAT || ismri == IPB_2D_FORMAT 
	|| ismri == ACRNEMA2L_FORMAT || ismri == ACRNEMA2B_FORMAT || ismri == ACRNEMA1L_FORMAT || ismri == ACRNEMA1B_FORMAT
	|| ismri == SMIS_FORMAT || ismri == ECAT7_FORMAT 
	|| ismri == DICMBIMP_FORMAT || ismri == DICMLIMP_FORMAT
	|| ismri == DICMBEXP_FORMAT || ismri == DICMLEXP_FORMAT
	|| ismri == DICM3DBIMP_FORMAT || ismri == DICM3DLIMP_FORMAT
	|| ismri == DICM3DBEXP_FORMAT || ismri == DICM3DLEXP_FORMAT
	|| ismri == ANALYZE_FORMAT_LITTLE_ENDIAN || ismri == ANALYZE_FORMAT_BIG_ENDIAN
	|| ismri == XML_FILEFORMAT
	)) return(-1);  
  
  /*BERST gestion images temporaires des automates
    if(pos>MAX_PICTURES_3D || pos<-MAX_PICTURES_3D) {
    PUT_WARN("Error Loading Image |position|>MAX_PICTURES_3D\n"); 
    return(-1);
    }
  */    
  graphic=ptr_img_3d(pos);
  
  /******* Liberation de la memoire reservee pour la labellisation *****/    
  if(graphic)
  {
    if(graphic->tala.labelbox != NULL)
    {
      FREE(graphic->tala.labelbox); // labelbox de l'ancienne image lib�� 
    }
  }
  /*********************************************************************/
  
  graphic->pos=pos;
  SetImg3d(pos,graphic);
  	
  if (graphic->mask)
    ptr_mask_activate_3d_ptr(graphic, FALSE);

  load_mri_3d(file,graphic,numfirst,max,step);
  
  /*  Initialisation de la table de couleur pour affichage */
  IMXColor_BindColormap_3d(pos, graphic->cmapinfo.noColormap); 
  IMXColor_UpdateColormap();  
	
	
  /*  Recherche du zoom de visu pour toujours utiliser l'espace maximum */
  old_zoom = graphic->zoom;
  temp1= (int) MAX_WIDTH_3D/graphic->width;
  temp2= (int) MAX_HEIGHT_3D/graphic->height;
  zoom=MINI(temp1,temp2);
  temp1= (int) MAX_DEPTH_3D/graphic->depth;
  zoom=MINI(temp1,zoom);
  if (old_zoom < zoom)
    graphic->zoom = zoom; 
                
   show_picture_3d(pos);     
   
  return(1);                                 
}

int    load_picture_3d(int pos, char *file, int max, int step, int numfirst)
{
  grphic3d *graphic;
  int ismri;
  int temp1,temp2,zoom,old_zoom;
  char *buff;
  int  n=0, n2=0;
  char *file2;
       
  
  if(_bIsParalActivated)
    { // TODO:PARAL
      n    = strlen(file);
      buff = CALLOC(n + 3, char);
      strcpy(buff, file);
        
      if(buff[n - 4] == '.')
        {
	  buff[n - 4] = '\0';
	  n2 = 4;
        }
       
      strcat(buff, ".*");

      file2 = IMXAtm_RCP(buff);

      n = strlen(file2) + n2;
      file = CALLOC(n + 1, char);
      sprintf(file, "%s.%s", file2, &file[strlen(file) - n]);
 
      FREE(buff);
      FREE(file2);
    } 

  ismri = is_mri(file,NULL);

  if (!(ismri == INTERFILE_FORMAT || ismri == BRUKER_FORMAT || ismri == SIEMENS_FORMAT	    ||ismri == PARAVISION_FORMAT
        || ismri == IPB_FORMAT || ismri == VANDERBILT_FORMAT 
	|| ismri == IPB_3D_FORMAT || ismri == IPB_2D_FORMAT 
	|| ismri == ACRNEMA2L_FORMAT || ismri == ACRNEMA2B_FORMAT || ismri == ACRNEMA1L_FORMAT || ismri == ACRNEMA1B_FORMAT
	|| ismri == SMIS_FORMAT || ismri == ECAT7_FORMAT 
	|| ismri == DICMBIMP_FORMAT || ismri == DICMLIMP_FORMAT
	|| ismri == DICMBEXP_FORMAT || ismri == DICMLEXP_FORMAT
	|| ismri == DICM3DBIMP_FORMAT || ismri == DICM3DLIMP_FORMAT
	|| ismri == DICM3DBEXP_FORMAT || ismri == DICM3DLEXP_FORMAT
	|| ismri == ANALYZE_FORMAT_LITTLE_ENDIAN || ismri == ANALYZE_FORMAT_BIG_ENDIAN
	|| ismri == XML_FILEFORMAT || ismri == NIFTI_FORMAT || ismri == NIFTI_GZ_FORMAT
	)) return(-1);  
  
  /*BERST gestion images temporaires des automates
    if(pos>MAX_PICTURES_3D || pos<-MAX_PICTURES_3D) {
    PUT_WARN("Error Loading Image |position|>MAX_PICTURES_3D\n"); 
    return(-1);
    }
  */    
  graphic=ptr_img_3d(pos);
  
  /******* Liberation de la memoire reservee pour la labellisation *****/    
  if(graphic)
  {
    if(graphic->tala.labelbox != NULL)
    {
      FREE(graphic->tala.labelbox); // labelbox de l'ancienne image lib�� 
    }
  }
  /*********************************************************************/
  
  graphic->pos=pos;
  SetImg3d(pos,graphic);
  	
  if (graphic->mask)
    ptr_mask_activate_3d_ptr(graphic, FALSE);

  load_mri_3d(file,graphic,numfirst,max,step);
  
  /*  Initialisation de la table de couleur pour affichage */
  IMXColor_BindColormap_3d(pos, graphic->cmapinfo.noColormap); 
  IMXColor_UpdateColormap();  
	
	
  /*  Recherche du zoom de visu pour toujours utiliser l'espace maximum */
  old_zoom = graphic->zoom;
  temp1= (int) MAX_WIDTH_3D/graphic->width;
  temp2= (int) MAX_HEIGHT_3D/graphic->height;
  zoom=MINI(temp1,temp2);
  temp1= (int) MAX_DEPTH_3D/graphic->depth;
  zoom=MINI(temp1,zoom);
  if (old_zoom < zoom)
    graphic->zoom = zoom; 
                
   show_picture_3d(pos);     
   
  return(1);                                 
}

/* --imx_patient_identity() ------------------------------------------------
**	retrieve patient identity & fill struct grphic3d
**	
**
**  file : nom du fichier 
**  graphic : structure graphic3d a remplir
** ----------------------------------------------------------------*/
int imx_patient_identity(const char* file, grphic3d *graphic)
{
  char *tmp;
  tmp=CALLOC(128,char);

  if (getmri_name(file))	strcpy(tmp,getmri_name(file));
  else strcpy(tmp,"");
  if (tmp != (char*)NULL) strcpy(graphic->patient.name,tmp);
  else bzero(graphic->patient.name,sizeof(graphic->patient.name));        

  if (getmri_dbirth(file)) strcpy(tmp,getmri_dbirth(file));
  else strcpy(tmp,"");
  if (tmp != (char*)NULL) strcpy(graphic->patient.d_birth,tmp); 
  else bzero(graphic->patient.d_birth,sizeof(graphic->patient.d_birth));       

  if (getmri_dexamen(file)) strcpy(tmp,getmri_dexamen(file));
  else strcpy(tmp,"");
  if (tmp != (char*)NULL) strcpy(graphic->patient.d_examen,tmp);
  else bzero(graphic->patient.d_examen,sizeof(graphic->patient.d_examen));       

  if (getmri_medecin(file)) strcpy(tmp,getmri_medecin(file));
  else strcpy(tmp,"");
  if (tmp != (char*)NULL) strcpy(graphic->patient.medecin,tmp);
  else bzero(graphic->patient.medecin,sizeof(graphic->patient.medecin));


  if(tmp) FREE(tmp);
  return 1;
}

/************** reorient_If_Negative_Voxel_Dim()***********
**
**   Detecte, signale et reoriente les axes de l'image
**   entres en argument si ceux-ci ont une dimension
**   de voxel negative.
**
**   retourne : nombre de dimensions negatives trouvees
**********************************************************/
int reorient_if_negative_voxel_dim(grphic3d *graphic)
{
  int NegativeVoxDim = 0; // += 1 si dx dy ou dy est n�atif
  
  /* Dans le cas d'une inversion de l'image sur X */
  if (graphic->dx < 0)
  {
    Miroir_ttaxe_3d_p(graphic,graphic,1); // l'argument "1" indique que le mirroir doit etre fait sur x
    NegativeVoxDim++;
  }

  /* Dans le cas d'une inversion de l'image sur Y */
  if (graphic->dy < 0)
  {
    Miroir_ttaxe_3d_p(graphic,graphic,2); // l'argument "2" indique que le mirroir doit etre fait sur y
    NegativeVoxDim++;
  }

  /* Dans le cas d'une inversion de l'image sur Z */
  if (graphic->dz < 0)
  {
    Miroir_ttaxe_3d_p(graphic,graphic,3); // l'argument "3" indique que le mirroir doit etre fait sur z
    NegativeVoxDim++;
  }

  /* Signalement du fait que l'image a retourne */
  if (NegativeVoxDim != 0)
    PUT_ERROR("dx<0 ou dy<0 ou dz<0 !");

  /* mise a jour des dimensions */
  graphic->dx = fabs(graphic->dx); 
  graphic->dy = fabs(graphic->dy); 
  graphic->dz = fabs(graphic->dz); 

  return NegativeVoxDim;
}

/* --load_mri_3d() ------------------------------------------------
**
**    Chargement d'une image 3d
**    pos : numero image 3d
**    file : nom du fichier
**    max : nombre de coupes 
**    step : pas pour la recherche dans fichiers multiples
**    numfirst : numero de la premiere image
**
**    max et step ne sont actuellement utilise 
**    que pour BRUKER_FORMAT et SIEMENS_FORMAT
**
**    ATTENTION : Bien verifier si width, height, et depth de graphic
**                correspondent a ceux de file surtout pour les formats
**                PARAVISION, VANDERBILT et INTERFILE car il y a des inversions
**                entre width, height et depth.
**                Les fonctions get_width, get_height, et get_depth 
**                donnent les bonnes valeurs de width height et depth
**                qui seront dans mri.
-----------------------------------------------------------------*/
int    load_mri_3d(const char *file, grphic3d *graphic, int numfirst, int max, int step)
{
  char file1[FILE_LEN];
  char s[64];
  char *file_mri;
  char *tmp=NULL;
  char v[80];				/* pour le format PARAVISION */
  char answer[128];
  int wdth, hght, dpth;
  int i, j, k, taille;
  int retour, err, format;
  int numfirst_sav;
  int minpixel,maxpixel;
  //	int read_format;
  int nbImg=0;
  int temp,zoom;
  float amax=0.,amin=0.;
  int readformat;
  float value;

  //	grphic3d *gtmp;

    
  strcpy(file1, file);
  numfirst_sav = numfirst;
  format = is_mri(file, NULL);

  retour=1;
  taille=2; /* Valeur par defaut */


  wdth = getmri_width(file, 1, MODE_3D);
  hght = getmri_height(file, 1, MODE_3D);
  dpth = getmri_depth(file, 1, MODE_3D);
  /*	if (!wdth || !hght || !dpth)
	{
	char tmp[128];
	sprintf(tmp,"probleme de dimension, wdth=%d, hght=%d, dpth=%d\n",wdth,hght,dpth);
	PUT_WARN( tmp);
	return (-1);
	}
  */
  graphic->min_pixel = getmri_minpixel(file, 1, MODE_3D);
  graphic->max_pixel = getmri_maxpixel(file, 1, MODE_3D);
  graphic->cutoff_max = getmri_cutoff_max(file, 1, MODE_3D);
  graphic->cutoff_min = getmri_cutoff_min(file, 1, MODE_3D);
  graphic->icomp = getmri_icomp(file, 1, MODE_3D);
  graphic->rcoeff = (double)pow((double)2, (double)graphic->icomp);
  graphic->dx = getmri_dx(file, 1, MODE_3D);
  graphic->dy = getmri_dy(file, 1, MODE_3D);
  graphic->dz = getmri_dz(file, 1, MODE_3D);
  graphic->mask_type = getmri_mask_type(file, 1, MODE_3D);
  graphic->zoom   =1;
  graphic->zoom_x =0;
  graphic->zoom_y =0;
  graphic->zoom_z =0;

  graphic->dti.is_dti = getmri_dti_is_dti(file, 1, MODE_3D);

	if (graphic->dti.is_dti) /* collecte des infos dti */
		{
	  graphic->dti.nb_directions = getmri_dti_nb_directions(file, 1, MODE_3D);
	  graphic->dti.gx = getmri_dti_gx(file, 1, MODE_3D);
	  graphic->dti.gy = getmri_dti_gy(file, 1, MODE_3D);
	  graphic->dti.gz = getmri_dti_gz(file, 1, MODE_3D);
	  graphic->dti.b = getmri_dti_b(file, 1, MODE_3D);
		}


  /*  Lecture fichiers format BRUKER 2D, Siemens 2D, SMIS 2D et DICOM 2D*/
  switch (format)
    {
    case (BRUKER_FORMAT)	: 
    case (SIEMENS_FORMAT)	:
    case (ACRNEMA2L_FORMAT)	:
    case (ACRNEMA2B_FORMAT)	:
    case (ACRNEMA1L_FORMAT)	:
    case (ACRNEMA1B_FORMAT)	:
    case (SMIS_FORMAT)     	:
    case (DICMBIMP_FORMAT) 	:
    case (DICMLIMP_FORMAT)	:
    case (DICMBEXP_FORMAT) 	:
    case (DICMLEXP_FORMAT) 	:

      graphic->width = wdth ;
      graphic->height = hght ;
      if (format != SMIS_FORMAT) graphic->depth = MINI(max, MAX_DEPTH_3D);
      else graphic->depth = dpth ;

      /* Recherche idendite du patient */

      imx_patient_identity(file, graphic);

      /* recherche du nombre maxi de fichiers pouvant etre lu */
      nbImg=graphic->depth;
      for (dpth=0;dpth<graphic->depth;dpth++)
	{
#ifndef WIN32
	  if ((access(file1, 0)!=0)) nbImg--;
#else // WIN32
	  if ((_access(file1, 0)!=0)) nbImg--;
#endif
	  incfic_img(file1,step,&numfirst);
	}
      strcpy(file1, file);
      numfirst = numfirst_sav;

      /*  Recherche du minimum et maximum de l'image 3D */
      /*        amax = 0.; 
		amin=LONGMAX;
		for (dpth=0; dpth<nbImg; dpth++) {
		amax = MAXI(getmri_maxpixel(file1, 1, MODE_3D) *((double)pow((double)2, (double)getmri_icomp(file1,
                1, MODE_3D))), amax);
		amin = MINI(getmri_minpixel(file1, 1, MODE_3D) *((double)pow((double)2, (double)getmri_icomp(file1,
                1, MODE_3D))), amin);
		incfic_img(file1, step, &numfirst);
		}
		numfirst = numfirst_sav;
		if (MESG3D_DEBUG) printf(" amax=%f \n", amax);
      */

      /*  Chargement de l'image 3D    
      **  apres calcul de icomp et de rcoeff max_pixel de l'image 3D  
      **  afin de tenir compte de la taille de mri (SHORT, LONG, ..)   */
      /*		if ( MAXI(amax,fabs((double) amin))  > MAXMRI3D )
			{
			err = imx_brukermax_3d(amax,amin,graphic);
			}
      */
      /* Nouvelle gestion du maxpixel JPA Mai 2002 */
      if ( sizeof(TYPEMRI3D) < getmri_dattype(file, 0, MODE_3D) ) 
	{
	  amax=graphic->rcoeff*graphic->max_pixel;
	  amin=graphic->rcoeff*graphic->min_pixel;
	  err=imx_brukermax_3d(amax,amin,graphic);
	}
        
      /* Lecture des fichiers */
      temp = MINI(nbImg, MAX_DEPTH_3D);
			
      graphic->depth = temp;

      strcpy(file1, file);
      for(dpth=0; dpth<temp; dpth++)  
	{
	  if((load_mri_bruker_3d(file1, (grphic3d *) graphic, dpth, numfirst)) == NULL) 
	    {
	      _imagix_error = IMX_BAD_POSITION;
	      _imagix_err_ext = errno;
	      retour=-1;				
	    }
	  incfic_img(file1, step, &numfirst);		
	}

      temp=MAXI(wdth,hght);
      zoom=1;
			
      if (temp>MAX_WIDTH_3D || temp>MAX_HEIGHT_3D)
	zoom=ceil((float)temp/((float)MAXI(MAX_WIDTH_3D,MAX_HEIGHT_3D)));
      graphic->dx = graphic->dx*zoom;
      graphic->dy = graphic->dy*zoom;
      graphic->dz = graphic->dz*zoom;
      graphic->width=graphic->width/zoom;
      graphic->height=graphic->height/zoom;
      graphic->depth = MINI(dpth, MAX_DEPTH_3D);

      /* Recherche du max et min pixel de l'image 3D */
      /*MISE EN PLACE DES BONNES VALEURS DANS WIDHT HEIGHT et DEPTH, dx, dy, dz */
      wdth = ceil(((double)wdth)/((double)zoom));
      hght = ceil(((double)hght)/((double)zoom));
      find_whd_from_orient_p(graphic, wdth, hght, dpth, getmri_readformat(file));

      // Recherche du max et min pixel de l'image 3D */
      minpixel=amax;
      maxpixel=0;
      for (i=0; i<graphic->width; i++)
	for (j=0; j<graphic->height; j++)
	  for (k=0; k<graphic->depth; k++)
	    {
	      if (graphic->mri[i][j][k]>maxpixel) maxpixel=graphic->mri[i][j][k];
	      if (graphic->mri[i][j][k]<minpixel) minpixel=graphic->mri[i][j][k];
	    }
      graphic->min_pixel = minpixel;
      graphic->max_pixel = maxpixel;
      graphic->cutoff_max = maxpixel;
      graphic->cutoff_min = minpixel;
      break;
      /* Fin chargement image BRUKER SIEMENS DICOM 2D (multi-fichier) */ 

      /* chargement DICOM 3D */
    case (DICM3DBIMP_FORMAT)	:
    case (DICM3DLIMP_FORMAT)	:
    case (DICM3DBEXP_FORMAT)	:
    _is_filebigEndian = getmri_endian(file);
    case (DICM3DLEXP_FORMAT)	:
		
      /* Recherche idendite du patient */
      imx_patient_identity(file, graphic);

      /* fin de la recherche de l'identite du patient */		
      graphic->width = MINI(wdth, MAX_WIDTH_3D);
      graphic->height = MINI(hght, MAX_HEIGHT_3D);
      graphic->depth = MINI(dpth, MAX_DEPTH_3D);
      fprintf(stderr, "file:%s\n", file);
      taille = getmri_dattype(file, 1,0);
      _read_format_3d = getmri_readformat(file);    /*format Siemens Vanderbilt*/
      _is_filebigEndian = getmri_endian(file);
      read_raw_3d(file, numfirst, wdth, hght, dpth, graphic, taille,format);
      /* resolution de maxpixel et minpixel lint    load_mri_3d(char *file, grphic3d *graphic, orsque les tags manquent
	 a l'appel  */
      if (graphic->max_pixel == 0) 
	imx_inimaxminpixel_3d_p(graphic);
      fprintf(stderr, "numfirst: %d, width: %d, height: %d, depth: %d, taille: %d\n", numfirst, wdth, hght, dpth,
	      taille);
      break;
      /* Fin chargement image DICOM 3D */ 


      /* chargement image type Siemens Vanderbilt */
    case (VANDERBILT_FORMAT )	:
    
      graphic->width = MINI(wdth, MAX_WIDTH_3D);
      graphic->height = MINI(hght, MAX_HEIGHT_3D);
      graphic->depth = MINI(dpth, MAX_DEPTH_3D);
      fprintf(stderr, "file:%s\n", file);
      taille = getmri_dattype(file, 1,0);
      readformat = _read_format_3d;
      fprintf(stderr, "readformat:'%d'", readformat);
      _read_format_3d = getmri_readformat(file);    /*format Siemens Vanderbilt*/
      file_mri = getfilename_img(file);
      fprintf(stderr, "file_img a lire:%s\n", file_mri);
      _is_filebigEndian = getmri_endian(file);
      read_raw_3d(file_mri, numfirst, wdth, hght, dpth, graphic, taille,format);
      fprintf(stderr, "numfirst: %d, width: %d, height: %d, depth: %d, taille: %d\n", numfirst, wdth, hght, dpth,
	      taille);
      free(file_mri);
      break;
      /* fin siemens Vanderbilt */

      /* chargement image type ECAT7 */
    case (ECAT7_FORMAT)	:
      wdth = getmri_width(file, numfirst, MODE_3D);
      hght = getmri_height(file, numfirst, MODE_3D);
      dpth = getmri_depth(file, numfirst, MODE_3D);
      /* dans ce cas, max ne represente pas le nombre de coupes
	 mais le numero de la matrice dans le fichier ECAT */
      max=dpth;
      graphic->min_pixel = getmri_minpixel(file, 1, MODE_3D);
      graphic->max_pixel = getmri_maxpixel(file, 1, MODE_3D);
      graphic->cutoff_max = getmri_cutoff_max(file, 1, MODE_3D);
      graphic->cutoff_min = getmri_cutoff_min(file, 1, MODE_3D);
      graphic->icomp = getmri_icomp(file, numfirst, MODE_3D);
      graphic->rcoeff = (double)pow((double)2, (double)graphic->icomp);
      if ((graphic->dx = getmri_dx(file, numfirst, MODE_3D)) < 0)
      {
        graphic->dx = - (graphic->dx);
      }
      graphic->dy = getmri_dy(file, numfirst, MODE_3D);
      graphic->dz = getmri_dz(file, numfirst, MODE_3D);
      graphic->width = MINI(wdth, MAX_WIDTH_3D);
      graphic->height = MINI(hght, MAX_HEIGHT_3D);
      graphic->depth = MINI(dpth, MAX_DEPTH_3D);

      /* lecture des informations du patient ... */
      imx_patient_identity(file, graphic);
			
      graphic->icomp=getmri_icomp(file,1,MODE_3D);
      graphic->rcoeff=(double)pow(2,(graphic->icomp));

      /* chargement du buffer */
      _read_format_3d = getmri_readformat(file);
      _is_filebigEndian = getmri_endian(file);
      read_raw_3d(file,numfirst,wdth,hght,dpth,graphic,taille,format);

      free(tmp);
      break;
      /* fin ECAT7 */


      /*  Chargement image type interfile */  
    case (INTERFILE_FORMAT):
      // Recherche idendite du patient 
      imx_patient_identity(file, graphic);

      // Recherche nom fichier data 
      file_mri = getfilename_img(file);
      _is_filebigEndian = getmri_endian(file);
      taille = getmri_dattype(file, 1,0);
      _read_format_3d = getmri_readformat(file);
      read_raw_3d(file_mri, numfirst, wdth, hght, dpth, graphic, taille,format);
      free(file_mri);
      break;  
      /* fin interfile */
	
      /*  Chargement image 3D type PARAVISION   */
    case (PARAVISION_FORMAT) :
      /* les champs dy et dz sont inverses pour l'angiographie */
      file_mri = getfilename_imnd(file);
      strcpy(v, getheader_interfile(file_mri,
				    "##$IMND_slice_scheme=", ASCII8B, 0));
      if (!strcmp(v, "Angio"))
	{
	  value = graphic->dy;
	  graphic->dy = graphic->dz;
	  graphic->dz = value;
	}
      free(file_mri);
      _read_format_3d = getmri_readformat(file);
      file_mri = getfilename_2dseq(file);
      taille = getmri_dattype(file, 1,0);
      /* lecture des informations du patient ... */
      _is_filebigEndian = getmri_endian(file);
      read_raw_3d(file_mri, numfirst, wdth, hght, dpth, graphic, taille,format);
      free(file_mri);
      break;
      /* fin image 3D type PARAVISION */
      imx_patient_identity(file, graphic);

      /*  Chargement image 3D type IPB  (vieux format)*/
    case (IPB_FORMAT) 	:

      file_mri = getfilename_img(file);
      taille = getmri_dattype(file, 1,0);
      readformat = _read_format_3d;
      _read_format_3d = getmri_readformat(file); /*format IPB perso*/
      _is_filebigEndian = getmri_endian(file);
      /*ATTENTION la fonction read_raw_3d ne permet pas de lire une image
        a une certaine position dans
        le fichier. Tout le reste (header avec indexage des parametres)
        est prevu pour pouvoir stocker
        et lire plusieurs images 3d. Pour cela il faudra modifier read_raw_3d
        en ajoutant un parametre
        relatif a la position dans le fichier de l'image que l'on souhaite lire */
      read_raw_3d(file_mri, numfirst, wdth, hght, dpth, graphic, taille,format);
      free(file_mri);
      break;
      /* fin image 3D type IPB  (vieux format) */

      /*  Chargement image 3D a partir d'une suite d'image 2D type IPB  */
    case (IPB_2D_FORMAT)	:
      file_mri = getfilename_img(file);
      taille = getmri_dattype(file, 1,0);
      readformat = _read_format_3d;
      _read_format_3d = getmri_readformat(file); /*format IPB perso*/
      /*ATTENTION la fonction read_raw_3d ne permet pas de lire une image
        a une certaine position dans
        le fichier. Tout le reste (header avec indexage des parametres)
        est prevu pour pouvoir stocker
        et lire plusieurs images 3d. Pour cela il faudra modifier read_raw_3d
        en ajoutant un parametre
        relatif a la position dans le fichier de l'image que l'on souhaite lire */
      read_raw_3d(file_mri, numfirst, wdth, hght, dpth, graphic, taille,format);
      free(file_mri);
      break;
      /* fin image 3D a partir d'une suite d'image 2D type IPB */
    /***************** lecture du xml ***********************
     
     Pour fonctionner, les fonctions relatives au traitement des fichiers atlas
     en xml font appel a deux classes : IO_xml et PropertyManager.
     Ces classes utilisent tres fortement les deux conteneurs "any" et "variant"
     de Boost (www.boost.org).
     Les classes IO_xml et PropertyManager sont �mentionner dans les Makefile
     respectifs noyau/io/Makefile et noyau/Makefile, en ajoutant completant
     ainsi : SRCXX = io_xml.cpp           pour noyau/io/Makefile
             SRCXX = propertymanager.cpp  pour noyau/Makefile
     
     Concernant les conteneurs "any" et "variant", il faut telecharger la
     librairie boost_1_32_0, la decompresser et l'installer en suivant les
     instructions presentes sur le site a l'url:
     http://www.boost.org/more/getting_started.html
     
     Ceci fait il ne reste plus qu'a modifier le fichier "configure" de sorte
     que les conteneurs soient trouves en lancant la compilation avec l'option :
     --with-boost= chemin d'acces vers le repertoire qui contient 
     le dossier "boost"
     ex: --with-boost=/home/henn/lib/boost-1_32/  
    
     */ 
   
#ifdef HAS_CORRECT_BOOST
   case (XML_FILEFORMAT)	:
      {
/*Olexa long taille;
        void* ioxml = newIOxml("","");
        int fileIsBigEndian;        

        if(xmlFileGetDataInTree(ioxml, file, graphic, numfirst, &fileIsBigEndian) != 0)
          return (-1);

        taille = (long)(graphic->bitppixel/8);
   
        read_raw_xml_3d(file, numfirst,
                        (long)graphic->width, (long)graphic->height, (long)graphic->depth,
                        graphic, taille, fileIsBigEndian);

        FREE(ioxml);
        return(1);
*/
	return(load_mri_atlas_3d((char*)file, graphic, numfirst, (int) 0));
    
	}
#endif
#ifdef HAS_NIFTI
    case (NIFTI_GZ_FORMAT) :
    case (NIFTI_FORMAT) :
      {
	if(!read_data_in_file_nifti(file, graphic)){
	  sprintf(s,"Unable to load 3d mri from file \"%s\"\n", file);
	  return -1;
	}

        return 1;
      }
#endif

      /********************************************************/
      /*   Lecture fichier format IPB_3D_FORMAT   */
    case (IPB_3D_FORMAT)	:
      {
	int num_tmp=1;
	//Dans le cas d'une image avec un mask la N-eme image n'est pas forcement a l'offset N
	//Il faut donc recuperer sa position dans le fichier a l'aide du tag  ORIGINAL_OFFSET
	strcpy(answer,getheader_interfile(file,"ORIGINAL_OFFSET=",ENTIER_IPB_2D,numfirst));
	if (strlen(answer)>0)
	  num_tmp = atoi(answer);
	else
	  num_tmp = numfirst;

         if (!graphic->mri)
	  graphic->mri = cr_mri_3d(wdth, hght, dpth); 
  
        if (load_mri_ipb_3d(file, graphic, num_tmp) != 1) 
	  {
	  sprintf(s,"Unable to load 3d mri from file \"%s\"\n", file);
	  PUT_ERROR(s);
	  return(-1);
	  }	  
 
	// lecture du mask de l'image s'il existe
	strcpy(answer,getheader_interfile(file,"MASK_OFFSET=",ENTIER_IPB_2D,numfirst));
	if(strlen(answer)>0)
	  {			
	    num_tmp = atoi(answer);
	    cr_mask_3d_ptr(graphic);			
	    graphic->mask_active=TRUE;

	    if (load_mri_ipb_3d(file, graphic->mask, num_tmp) != 1)
	      {
		sprintf(s,"Unable to load 3d mri from file \"%s\"\n", file);
		PUT_ERROR(s);
		return(-1);
	      }
	    //  Initialisation du mask 
	    if (graphic->mask)
	      {	
		init_mask_3d_ptr(graphic);
		if (graphic->pos!=99) // BERST si pos = 99 ne pas afficher l'image ????
		  ptr_mask_activate_3d_ptr(graphic, TRUE);
	      }	
	  }
	// Pour placer le repère correctement
	graphic->x3dc = getmri_x3dc(file, 1, MODE_3D);
	graphic->y3dc = getmri_y3dc(file, 1, MODE_3D);
	graphic->z3dc = getmri_z3dc(file, 1, MODE_3D);
	IMX3D_SetXYZ(graphic->pos,graphic->x3dc,graphic->y3dc,graphic->z3dc);
	IMX3D_RefreshObjects3D(graphic->pos);
	return(1);
      }
      break;
      /* fin format IPB_3D_FORMAT  */

      /*   Lecture fichier format ANALYZE_FORMAT   */
    case (ANALYZE_FORMAT_LITTLE_ENDIAN)	:
    case (ANALYZE_FORMAT_BIG_ENDIAN)	:
      if ( load_mri_analyze_3d(file, graphic, numfirst) != 1 ) 
	    {
	      sprintf(s,"Unable to load 3d mri from file \"%s\"\n", file);
	      PUT_ERROR(s);
	      return(-1);
	    }
        return(1);
      break;
		
    default : 
     break;
    }
  /* fin format format ANALYZE_FORMAT   */

  /* Suite et partie commune a la lecture de quelques formats */
  strcpy(graphic->filename, file);

  graphic->x3dc = graphic->width/2;
  graphic->y3dc = graphic->height/2;
  graphic->z3dc = graphic->depth/2;

  graphic->zoom = 1;          /* Not determine yet !  */
  graphic->zoom_x = 0;        /* Not determine yet !  */
  graphic->zoom_y = 0;        /* Not determine yet !  */
  graphic->zoom_z = 0;        /* Not determine yet !  */
  graphic->img_type = IMAGE_NORMAL;
  graphic->type = 'b';
  graphic->bitppixel = 32;
  
  printf("For file width=%d height=%d depth=%d bytes par pixel=%d format de lecture=%d\n \
	        retenu : width=%d height=%d depth=%d\n", wdth, hght, dpth, taille,_read_format_3d
	 , graphic->width, graphic->height, graphic->depth);

  /*  RAZ du reste de l'image a visualiser sur l'ecran */
  if (graphic->pos != 99) 
    {
      for (i=graphic->width+1; i<MAX_WIDTH_3D; i++)
        for (j=graphic->height+1; j<MAX_HEIGHT_3D; j++)
          for (k=graphic->depth+1; k<MAX_DEPTH_3D; k++)
	    graphic->mri[i][j][k] = 0;
    }		

  /*  Si le maxpixel est nul bien verifier s'il n'y pas eu d'erreur de lecture du maxpixel 
      en recherchant le max et le min dans l'image . */
  if (graphic->max_pixel == 0) 
    imx_inimaxminpixel_3d_p(graphic);   

  reorient_if_negative_voxel_dim(graphic);
  
  return(retour);
}


/* --load_picture_raw_3d() ------------------------------------------------
**
**    Chargement d'une image raw 3d
**    pos : numero image 3d
**    file : nom du fichier
**    max : nombre de coupes
**    step : pas pour la recherche dans fichiers multiples
**    numfirst : numero de la premiere image
**
-----------------------------------------------------------------*/
int    load_picture_raw_3d(const char *file, int numfirst, int wdth, int hght, int dpth, int pos, int taille)
{
  grphic3d *graphic;
  int i,j,k,size;
  //    int err;  		unused variable
  long ia,maxpixel=0;
  //    float amax;   	unused variable
  short *buff16b;
  long  *buff;
  unsigned char *buff8b;
  FILE *fp;

  graphic=ptr_img_3d(pos);
  SetImg3d(pos,graphic);

  graphic->width =MINI(wdth,MAX_WIDTH_3D);
  graphic->height=MINI(hght,MAX_HEIGHT_3D);
  graphic->depth=MINI(dpth,MAX_DEPTH_3D);
  printf( " width file=%d height file =%d depth file=%d \n retenu : width=%d height=%d depth=%d\n", wdth,hght,dpth,
	  graphic->width,graphic->height,graphic->depth);
  graphic->min_pixel=0;
  graphic->icomp=0;
  graphic->rcoeff=(double)pow((double)2,(double)graphic->icomp);

  graphic->dx =1;
  graphic->dy =1;
  graphic->dz =1;
  graphic->zoom = 1;
  graphic->zoom_x = 0;
  graphic->zoom_y = 0;
  graphic->zoom_z = 0;

  /*  Recherche du maximum de l'image 3D */
  fp=fopen(file,"rb");
  if (fp==NULL) {
    PUT_WARN("Unknown file in load_picture_raw\n");
    return (0);
  }
	  
  switch(taille)
    {
    case 1: /*On lit un fichier 8 bits*/    
      if((buff8b=CALLOC(wdth*hght,unsigned char))==NULL) { /* Cannot alloc... */
        fclose(fp);
        return(0);
      }
      for (k=0;k<dpth;k++) {
        size=fread((unsigned char*)buff8b,sizeof(unsigned char),wdth*hght,fp);
        for (i=0;i<wdth*hght;i++) {
	  ia=buff8b[i];
	  maxpixel=MAXI(maxpixel,ia);
        }
      }
      break;
    case 2: /*On lit un fichier 16 bits*/    
      if((buff16b=CALLOC(wdth*hght,short))==NULL) { /* Cannot alloc... */
        fclose(fp);
        return(0);
      }
      for (k=0;k<dpth;k++) {
        size=FREAD_SHORT((short*)buff16b,wdth*hght,fp);
        for (i=0;i<wdth*hght;i++) {
	  ia=buff16b[i];
	  maxpixel=MAXI(maxpixel,ia);
        }
      }
      break;
    case 4: /*On lit un fichier 32 bits*/
      if((buff=CALLOC(wdth*hght,long))==NULL) { /* CannotIMXColor_BindColormap_3d alloc... */
        fclose(fp);
        return(0);
      }
      for (k=0;k<dpth;k++) {
        size=FREAD_LONG((long*)buff,wdth*hght,fp);
        for (i=0;i<wdth*hght;i++) {
	  ia=buff[i];
	  maxpixel=MAXI(maxpixel,ia);
        }
      }
      break;
    }  
  fclose(fp);
  graphic->max_pixel=maxpixel;
  graphic->cutoff_max=maxpixel;

  /*  Chargement de l'image 3D*/
  read_raw_3d ( file,numfirst,wdth,hght,dpth,graphic,taille,0);

  // graphic->depth=MINI(dpth,MAX_DEPTH_3D);
  graphic->x3dc=graphic->width/2;
  graphic->y3dc=graphic->height/2;
  graphic->z3dc=graphic->depth/2;

  graphic->zoom=1;        /* Not determine yet ! */
  graphic->zoom_x=0;        /* Not determine yet ! */
  graphic->zoom_y=0;        /* Not determine yet ! */
  graphic->zoom_z=0;        /* Not determine yet ! */
  graphic->pos=-1;        /* Not determine yet ! */
  graphic->img_type =IMAGE_NORMAL;
  graphic->type='b';
  graphic->bitppixel=32;
  strcpy (graphic->filename,file);

  /*  RAZ du reste de l'image avisualiser sur l'ecran */
  for (i=graphic->width+1; i<MAX_WIDTH_3D; i++)
    for (j=graphic->height+1; j<MAX_HEIGHT_3D; j++)
      for (k=graphic->depth+1; k<MAX_DEPTH_3D; k++)
	graphic->mri[i][j][k] = 0;

  /* Mise en place table de couleur */
  IMXColor_BindColormap_3d_ptr(graphic, 1);
  IMXColor_UpdateColormap();

  show_picture_3d(pos);

  return(1);
}

/* --load_mri_bruker_3d() ----------------------------------------------------    
**
**    This function loads a mri 3D.
**
**    Usage:    load_mri_bruker_3d(file,graphic,dpth,numfirst)
**        file: char *: Designs file to load.
**        graphic: grphic:Where to stock file.
**
*/
grphic3d *load_mri_bruker_3d(const char *file, grphic3d *graphic, int dpth, int numfirst)
{
  long *buff;
  //	long max,i;    unused variable
  int number_img;
  float anc_rcoef,coeff;
  long maxpix,minpix;
  long format=0,read_format;
  int i,j,l,tmp,zoom;
  int tmp_wdth,tmp_hght;

  //    maxpix=graphic->max_pixel;
  //    amax=graphic->rcoeff*maxpix;
  maxpix = 0;
  minpix = 0;
  zoom = 1;
  read_format = getmri_readformat(file);

  format=is_mri(file,NULL);
  if(format==INTERFILE_FORMAT)
    {
      if ((number_img= atol(getheader_interfile(file,"!total number of images:=",ENTIER,0)) )<numfirst)
        {
	  char s[256];
	  sprintf(s,"Sorry, cannot load image number %d only %d images \n" ,numfirst,number_img);
	  PUT_MESG(s);
	  _imagix_error=IMX_PICTURE_NOT_FOUND;
	  _imagix_err_ext=errno;
	  return(NULL);
        }
    }

  if(format==BRUKER_FORMAT 
     || format==SIEMENS_FORMAT 
     || format==INTERFILE_FORMAT 
     || format==VANDERBILT_FORMAT
     || format==ACRNEMA2L_FORMAT
     || format==ACRNEMA2B_FORMAT
     || format==ACRNEMA1L_FORMAT
     || format==ACRNEMA1B_FORMAT
     || format==SMIS_FORMAT
     || format==DICMBIMP_FORMAT  || format==DICMLIMP_FORMAT 
     || format==DICMBEXP_FORMAT  || format==DICMLEXP_FORMAT)
    { 

      errno=0;
      strncpy(graphic->filename,file,256);
      if((buff=(long *)getmri(file,numfirst))==NULL) 
	{
	  char s[256];
	  sprintf(s,"Sorry, cannot load %s...(*errno:%d) \n",file,errno);
	  PUT_MESG(s);
	  _imagix_error=IMX_PICTURE_NOT_FOUND;
	  _imagix_err_ext=errno;
	  return(NULL);
        }
      anc_rcoef=pow((double)2,(double)getmri_icomp(file, 1, MODE_3D));
      coeff = anc_rcoef/((float) graphic->rcoeff);
		
      tmp = MAXI(graphic->width, graphic->height);
      tmp_wdth = graphic->width;
      tmp_hght = graphic->height;
		
      l=0;
      if (tmp>MAX_WIDTH_3D || tmp>MAX_HEIGHT_3D)
	{
	  zoom=ceil((float)tmp/((float)MAXI(MAX_WIDTH_3D,MAX_HEIGHT_3D)));
	  tmp_wdth = tmp_wdth/zoom;
	  tmp_hght = tmp_hght/zoom;
	}
		
      for(i=0;i<tmp_wdth;i++)
	{
	  for(j=0;j<tmp_hght;j++)
	    {	
	      fill_mri_32b(graphic, buff, i, j, dpth, l, coeff, read_format, 1);
	      l+=zoom;
	      // 	a voir plus tard  ... faire une moyenne des valeurs entre l et l+zoom
	      // parceque pour l'instant on omet les valeurs ]l;l+zoom[
				
	    }	
	  if (zoom>1)
	    l=l+(tmp_hght*zoom);
	}		
		
      FREE(buff);    /* Free width*height*4 + 18432 */
    }
  return (graphic);
}

#endif /* __GTK_GUI */
/* --load_mri_ipb_3d() ---------------------------------------------------------
**
**    load_mri_ipb_3d() loads a picture in IPB 3D file format.
**
**    Usage ...
*/
int load_mri_ipb_3d(const char *file, grphic3d *image, int number)
{
  int data_size;
  int len;
  char headerfile[256];
  char datafile[256];
  char datafilegz[256];
  char syscall[1000];
  int iscompressed=0;
  char crvfile[PATH_LEN];
  struct stat statbuff;
  char *s;
  char ss[100];
  char v[16]; 
  unsigned long offset;
  char *buffchar;
  short *buffshort;
  long *bufflong;
  void *buffer;
  int i, j, k, l, num_items;
  FILE *fp;
  int w_step, h_step, d_step, step, origWidth;
  int idCMap1;
  int is_filebigEndian = FALSE;
  //	int idCMap2;  unused variable
  int n1;

  data_size = sizeof(TYPEMRI3D);

  if(number <= 0)
    {
      sprintf(ss,"[load_mri_ipb_3d] no such image in file (number=%d)!\n", number);
      PUT_ERROR(ss);
      number = 1;
    }

  if(getmri_numberimg(file) < number) {
    sprintf(ss,"[load_mri_ipb_3d] no such image in file (number=%d)!\n", number);
    PUT_ERROR(ss);
    return(0);
  }

  if((len = strlen(file)) < 4) {
    sprintf(ss,"[load_mri_ipb_3d] Invalid file name \"%s\" !\n", file);
    PUT_ERROR(ss);
    return(0);
  }

  strcpy(headerfile, file);
  strcpy(datafile, file);
  strcpy(datafilegz, file);
  strcpy(crvfile, file);
  strcpy(datafile + len - 4, ".img");
  strcpy(datafilegz + len - 4, ".img.gz");
  strcpy(crvfile + len - 4, "");


  if(stat(headerfile, &statbuff) == -1) { /* File doesn't exist */
    sprintf(ss,"[load_mri_ipb_3d] file \"%s\" doesn't exists !\n", headerfile);
    PUT_ERROR(ss);
    return(0);
  }

  if(stat(datafile, &statbuff) == -1) 
    { /* File doesn't exist */
      if(stat(datafilegz, &statbuff) == -1) 
	{
	  sprintf(ss,"[load_mri_ipb_3d] file \"%s\" doesn't exists !\n", datafile);
	  PUT_ERROR(ss);
	  return(0);
	}
      iscompressed=1;
      sprintf(datafile,"/tmp/imx_tmp_decompress.img");
      sprintf(syscall,"gunzip -c %s > %s",datafilegz,datafile);
      system(syscall);
      if(stat(datafile, &statbuff) == -1) 
	{ 
	  PUT_ERROR("file decompress failed");
	}
    }

  if((s = getmri_name(file)) == NULL) image->patient.name[0] = '\0';
  else strncpy(image->patient.name, s,100);
  if (s) {free(s); s=NULL;}
  if((s = getmri_dbirth(file)) == NULL) image->patient.d_birth[0] = '\0';
  else strncpy(image->patient.d_birth, s,100);
  if (s) {free(s); s=NULL;}
  if((s = getmri_medecin(file)) == NULL) image->patient.medecin[0] = '\0';
  else strncpy(image->patient.medecin, s,100);
  if (s) {free(s); s=NULL;}
  if((s = getmri_dexamen(file)) == NULL) image->patient.d_examen[0] = '\0';
  else strncpy(image->patient.d_examen, s,100);
  if (s) {free(s); s=NULL;}
  strcpy(image->filename, headerfile);

  image->width    = getmri_width(file, number, MODE_3D);
  image->height   = getmri_height(file, number, MODE_3D);
  image->depth    = getmri_depth(file, number, MODE_3D);

  if(getmri_numberimg(headerfile) == 1) {
    int file_size;

    file_size = image->width * image->height * image->depth * getmri_dattype(headerfile, 1, MODE_3D);

    if(statbuff.st_size == 2 * file_size) {
      //            char s[256];		unused variable

      sprintf(ss,"[load_mri_ipb_3d] file size doesn't match : %s!\nPlease use ipb3dconv ...\n",file);
      PUT_ERROR(ss);
      return(0);
    }
    else {
      if(statbuff.st_size != file_size) {
	PUT_ERROR("[load_mri_ipb_3d] file size doesn't match !\n");
	return(0);
      }
    }
  }

  image->icomp = getmri_icomp(file, number, MODE_3D);
  /* Utiliser le rcoeff du header ? */
  image->rcoeff = (double)pow((double)2, (double)image->icomp);

  image->max_pixel = getmri_maxpixel(file, number, MODE_3D);
  image->min_pixel = getmri_minpixel(file, number, MODE_3D);
  image->cutoff_max = getmri_cutoff_max(file, number, MODE_3D);
  if(image->cutoff_max == 0) image->cutoff_max = image->max_pixel;
  image->cutoff_min = getmri_cutoff_min(file, number, MODE_3D);
  if(image->cutoff_min == 0) image->cutoff_min = image->min_pixel;
  image->dx = getmri_dx(file, number, MODE_3D);
  image->dy = getmri_dy(file, number, MODE_3D);
  image->dz = getmri_dz(file, number, MODE_3D);

#ifdef __GTK_GUI
  idCMap1 = getmri_nCMap(file, number, MODE_3D);

  if(idCMap1 == 0)
    n1 = IMXColor_NewColormap(7, 0, NULL, 0);
  else
    n1 = IMXColor_NewColormap(idCMap1, 0, NULL, 0);
 
  list_moveto(&_lstColormaps, n1);
  image->cmapinfo.noColormap  = n1;
  image->cmapinfo.cmap        = (ptrColormap) list_item(&_lstColormaps);
  image->cmapinfo.dspOverBlancking = getmri_dspOverBlancking(file, number, MODE_3D);
  //    IMXColor_BindColormap_3d_ptr(image, n1 );  modif mask
  //    IMXColor_UpdateColormap();   modif mask
  /*#else 
    getmri_nCMap(file, number, MODE_3D);*/
#endif /* __GTK_GUI */

  strcpy(v, getheader_interfile(headerfile, "x3dc=", ENTIER_IPB_2D, number));
  image->x3dc = atol(v);
  strcpy(v, getheader_interfile(headerfile, "y3dc=", ENTIER_IPB_2D, number));
  image->y3dc = atol(v);
  strcpy(v, getheader_interfile(headerfile, "z3dc=", ENTIER_IPB_2D, number));
  image->z3dc = atol(v);
  strcpy(v, getheader_interfile(headerfile, "zoom=", ENTIER_IPB_2D, number));
  image->zoom = atol(v);
  strcpy(v, getheader_interfile(headerfile, "zoom_x=", ENTIER_IPB_2D, number));
  image->zoom_x = atol(v);
  strcpy(v, getheader_interfile(headerfile, "zoom_y=", ENTIER_IPB_2D, number));
  image->zoom_y = atol(v);
  strcpy(v, getheader_interfile(headerfile, "zoom_z=", ENTIER_IPB_2D, number));
  image->zoom_z = atol(v);

  strcpy(v, getheader_interfile(headerfile, "AC_x=", REEL_IPB_2D, number));
  image->tala.AC.x = atof(v);
  strcpy(v, getheader_interfile(headerfile, "AC_y=", REEL_IPB_2D, number));
  image->tala.AC.y = atof(v);
  strcpy(v, getheader_interfile(headerfile, "AC_z=", REEL_IPB_2D, number));
  image->tala.AC.z = atof(v);
  strcpy(v, getheader_interfile(headerfile, "PC_x=", REEL_IPB_2D, number));
  image->tala.PC.x = atof(v);
  strcpy(v, getheader_interfile(headerfile, "PC_y=", REEL_IPB_2D, number));
  image->tala.PC.y = atof(v);
  strcpy(v, getheader_interfile(headerfile, "PC_z=", REEL_IPB_2D, number));
  image->tala.PC.z = atof(v);
  strcpy(v, getheader_interfile(headerfile, "AP_x=", REEL_IPB_2D, number));
  image->tala.AP.x = atof(v);
  strcpy(v, getheader_interfile(headerfile, "AP_y=", REEL_IPB_2D, number));
  image->tala.AP.y = atof(v);
  strcpy(v, getheader_interfile(headerfile, "AP_z=", REEL_IPB_2D, number));
  image->tala.AP.z = atof(v);
  strcpy(v, getheader_interfile(headerfile, "PP_x=", REEL_IPB_2D, number));
  image->tala.PP.x = atof(v);
  strcpy(v, getheader_interfile(headerfile, "PP_y=", REEL_IPB_2D, number));
  image->tala.PP.y = atof(v);
  strcpy(v, getheader_interfile(headerfile, "PP_z=", REEL_IPB_2D, number));
  image->tala.PP.z = atof(v);
  strcpy(v, getheader_interfile(headerfile, "RM_x=", REEL_IPB_2D, number));
  image->tala.RM.x = atof(v);
  strcpy(v, getheader_interfile(headerfile, "RM_y=", REEL_IPB_2D, number));
  image->tala.RM.y = atof(v);
  strcpy(v, getheader_interfile(headerfile, "RM_z=", REEL_IPB_2D, number));
  image->tala.RM.z = atof(v);
  strcpy(v, getheader_interfile(headerfile, "LM_x=", REEL_IPB_2D, number));
  image->tala.LM.x = atof(v);
  strcpy(v, getheader_interfile(headerfile, "LM_y=", REEL_IPB_2D, number));
  image->tala.LM.y = atof(v);
  strcpy(v, getheader_interfile(headerfile, "LM_z=", REEL_IPB_2D, number));
  image->tala.LM.z = atof(v);
  strcpy(v, getheader_interfile(headerfile, "SP_x=", REEL_IPB_2D, number));
  image->tala.SP.x = atof(v);
  strcpy(v, getheader_interfile(headerfile, "SP_y=", REEL_IPB_2D, number));
  image->tala.SP.y = atof(v);
  strcpy(v, getheader_interfile(headerfile, "SP_z=", REEL_IPB_2D, number));
  image->tala.SP.z = atof(v);
  strcpy(v, getheader_interfile(headerfile, "IP_x=", REEL_IPB_2D, number));
  image->tala.IP.x = atof(v);
  strcpy(v, getheader_interfile(headerfile, "IP_y=", REEL_IPB_2D, number));
  image->tala.IP.y = atof(v);
  strcpy(v, getheader_interfile(headerfile, "IP_z=", REEL_IPB_2D, number));
  image->tala.IP.z = atof(v);
  strcpy(v, getheader_interfile(headerfile, "IC_x=", REEL_IPB_2D, number));
  image->tala.IC.x = atof(v);
  strcpy(v, getheader_interfile(headerfile, "IC_y=", REEL_IPB_2D, number));
  image->tala.IC.y = atof(v);
  strcpy(v, getheader_interfile(headerfile, "IC_z=", REEL_IPB_2D, number));
  image->tala.IC.z = atof(v);

  strcpy(v, getheader_interfile(headerfile, "Tu_x=", REEL_IPB_2D, number));
  image->tala.Tu.x = atof(v);
  strcpy(v, getheader_interfile(headerfile, "Tu_y=", REEL_IPB_2D, number));
  image->tala.Tu.y = atof(v);
  strcpy(v, getheader_interfile(headerfile, "Tu_z=", REEL_IPB_2D, number));
  image->tala.Tu.z = atof(v);
  strcpy(v, getheader_interfile(headerfile, "Tv_x=", REEL_IPB_2D, number));
  image->tala.Tv.x = atof(v);
  strcpy(v, getheader_interfile(headerfile, "Tv_y=", REEL_IPB_2D, number));
  image->tala.Tv.y = atof(v);
  strcpy(v, getheader_interfile(headerfile, "Tv_z=", REEL_IPB_2D, number));
  image->tala.Tv.z = atof(v);
  strcpy(v, getheader_interfile(headerfile, "Tw_x=", REEL_IPB_2D, number));
  image->tala.Tw.x = atof(v);
  strcpy(v, getheader_interfile(headerfile, "Tw_y=", REEL_IPB_2D, number));
  image->tala.Tw.y = atof(v);
  strcpy(v, getheader_interfile(headerfile, "Tw_z=", REEL_IPB_2D, number));
  image->tala.Tw.z = atof(v);
  strcpy(v, getheader_interfile(headerfile, "x_AP=", REEL_IPB_2D, number));
  image->tala.x_AP = atof(v);
  strcpy(v, getheader_interfile(headerfile, "x_PC=", REEL_IPB_2D, number));
  image->tala.x_PC = atof(v);
  strcpy(v, getheader_interfile(headerfile, "x_PP=", REEL_IPB_2D, number));
  image->tala.x_PP = atof(v);
  strcpy(v, getheader_interfile(headerfile, "z_IP=", REEL_IPB_2D, number));
  image->tala.z_IP = atof(v);
  strcpy(v, getheader_interfile(headerfile, "z_SP=", REEL_IPB_2D, number));
  image->tala.z_SP = atof(v);
  strcpy(v, getheader_interfile(headerfile, "y_RM=", REEL_IPB_2D, number));
  image->tala.y_RM = atof(v);
  strcpy(v, getheader_interfile(headerfile, "y_LM=", REEL_IPB_2D, number));
  image->tala.y_LM = atof(v);

  if ((image->tala.AC.x != 0.0)&&(image->tala.AC.y != 0.0)&&(image->tala.AC.z != 0.0))
    image->tala.etat = 1;
  else image->tala.etat = 0;
  image->tala.p_sagi_ok = 0;
  image->tala.grid_id = NULL;
  image->tala.grid_on = 0;

  strcpy(v, getheader_interfile(headerfile, "is_dti=", REEL_IPB_2D, number));
  image->dti.is_dti = atoi(v);
	strcpy(v, getheader_interfile(headerfile, "nb_directions=", REEL_IPB_2D, number));
  image->dti.nb_directions = atoi(v);
  strcpy(v, getheader_interfile(headerfile, "gx=", REEL_IPB_2D, number));
  image->dti.gx = atof(v);
 	strcpy(v, getheader_interfile(headerfile, "gy=", REEL_IPB_2D, number));
  image->dti.gy = atof(v);
 	strcpy(v, getheader_interfile(headerfile, "gz=", REEL_IPB_2D, number));
  image->dti.gz = atof(v);
 	strcpy(v, getheader_interfile(headerfile, "b=", REEL_IPB_2D, number));
  image->dti.b = atof(v);

#ifndef COMPILE_FOR_MEDIPY
	read_annotations(image,headerfile,number);
#endif
	
#ifdef __GTK_GUI
  imx_rdcurv(crvfile, number, image);
#endif /* __GTK_GUI */

  image->bitppixel = getmri_dattype(headerfile, number, MODE_3D) * 8;
  is_filebigEndian = getmri_endian(file);

  strcpy(v, getheader_interfile(headerfile, "offset in raw=", ENTIER_IPB_2D, number));
  offset = atol(v);

  num_items = image->width * image->height;
  origWidth = image->width;

  w_step = (image->width + MAX_WIDTH_3D - 1) / MAX_WIDTH_3D;
  h_step = (image->height + MAX_HEIGHT_3D - 1) / MAX_HEIGHT_3D;
  d_step = (image->depth + MAX_DEPTH_3D - 1) / MAX_DEPTH_3D;
  step = (w_step > h_step)? w_step: h_step;
  step = (step > d_step)? step: d_step;
  if(step > 1) {
    sprintf(ss,"Warning : image has been reduced by a factor %d\n",step);
    PUT_WNDMSG(ss);
    image->width    = image->width/step;
    image->height   = image->height/step;
    image->depth    = image->depth/step;
    image->dx = image->dx * step;
    image->dy = image->dy * step;
    image->dz = image->dz * step;
    image->x3dc = image->x3dc/step;
    image->y3dc = image->y3dc/step;
    image->z3dc = image->z3dc/step;
  }

  if((buffer = malloc(num_items * data_size)) == NULL) {
    sprintf(ss,"[load_mri_ipb_3d] memory allocation error !\n");
    PUT_ERROR(ss);
    return(0);
  }

  if((fp = fopen(datafile, "rb")) == NULL) 
    {
      sprintf(ss,"[load_mri_ipb_3d] open error on file \"%s\" !\n", datafile);
      PUT_ERROR(ss);
      free(buffer);
      return(0);
    }

  fseek(fp, offset, SEEK_SET);

  switch(data_size) {
  case sizeof(char):
    for(k=0; k<image->depth; k++) {

      if(fread(buffer, data_size, num_items, fp) != num_items) {
	sprintf(ss,"[load_mri_ipb_3d] read error (1)!\n");
	PUT_ERROR(ss);
	free(buffer);
	fclose(fp);
	return(0);
      }

      buffchar = buffer;
      l=0;
      for(j=0; j<image->height; j++)
	{
	  for(i=0; i<image->width; i++)
	    image->mri[i][j][k] =  buffchar[l];l=l+step;
//	  bug - doesn't work for width not divisible by step
//	  l=l+image->width*step*(step-1);
	  l=l+step*(origWidth-image->width);	
	}	
      for(i=0; i<step-1; i++) {
	if(fread(buffer, data_size, num_items, fp) != num_items) {
	  sprintf(ss,"[load_mri_ipb_3d] read error (2)!\n");
	  PUT_ERROR(ss);
	  free(buffer);
	  fclose(fp);
	  return(0);
	}
      }
    }
    break;
  case sizeof(short):
    for(k=0; k<image->depth; k++) {

      if(fread(buffer, data_size, num_items, fp) != num_items) {
	sprintf(ss,"[load_mri_ipb_3d] read error (2)!\n");
	PUT_ERROR(ss);
	free(buffer);
	fclose(fp);
	return(0);
      }

      buffshort = buffer;
      l=0;
      for(j=0; j<image->height; j++)
	{
	  for(i=0; i<image->width; i++)
	    if (is_filebigEndian == _is_bigEndian) {image->mri[i][j][k] = buffshort[l];l=l+step;}
	    else	{image->mri[i][j][k] = shEndianConversion(buffshort[l]);l=l+step;}
//	  bug - doesn't work for width not divisible by step
//	  l=l+image->width*step*(step-1);
	  l=l+step*(origWidth-image->width);
	}		
      for(i=0; i<step-1; i++) {
	if(fread(buffer, data_size, num_items, fp) != num_items) {
	  sprintf(ss,"[load_mri_ipb_3d] read error (2)!\n");
	  PUT_ERROR(ss);
	  free(buffer);
	  fclose(fp);
	  return(0);
	}
      }
    }
    break;
  case sizeof(long):
    for(k=0; k<image->depth; k++) {

      if(fread(buffer, data_size, num_items, fp) != num_items) {
	sprintf(ss,"[load_mri_ipb_3d] read error (3)!\n");
	PUT_ERROR(ss);
	free(buffer);
	fclose(fp);
	return(0);
      }

      bufflong = buffer;
      l=0;
      for(j=0; j<image->height; j++) 
	{
	  for(i=0; i<image->width; i++)
	    if (is_filebigEndian == _is_bigEndian) {image->mri[i][j][k] = bufflong[l];l=l+step;}
	    else	{image->mri[i][j][k] = longEndianConversion(bufflong[l]);l=l+step;}
//	  bug - doesn't work for width not divisible by step
//	  l=l+image->width*step*(step-1);
	  l=l+step*(origWidth-image->width);
	}		
      for(i=0; i<step-1; i++) {
	if(fread(buffer, data_size, num_items, fp) != num_items) {
	  sprintf(ss,"[load_mri_ipb_3d] read error (2)!\n");
	  PUT_ERROR(ss);
	  free(buffer);
	  fclose(fp);
	  return(0);
	}
      }	
    }
    break;
  }

  free(buffer);
  fclose(fp);

  if(iscompressed)
    {
      sprintf(syscall,"rm %s",datafile);
      system(syscall);
    }
    
#ifndef COMPILE_FOR_MEDIPY
  loadLabels(image);
#endif 
  return(1);
}

#ifdef __GTK_GUI
/* --save_pictures_3d() ---------------------------------------------------------
**
**    save_pictures_3d() saves more pictures in IPB 3D file format.
**
**    Usage ...
**   file :       nom du fichier
**   start :      numero de la premier image a sauver 
**   num_images : nombre d'images a sauver 
**   step :       pas si plusieurs images
**
*/
int save_pictures_3d(char *file, int start, int num_images, int step)
{
  int err;
  int ipb_format;

  if (check_file_3d(file)!=0 ) return (0);

  ipb_format = check_ipb_format_3d(file);
  if (ipb_format && ipb_format!=IPB_WITHOUT_MASK_FORMAT) 
    { 
      PUT_WARN( "[save_pictures_3d] warning: detected ipb with mask format !\n");
      return (0);
    }
  err=save_mri_ipb_3d(file, start, num_images, step);

  return(err);
}

/* --save_pictures_3d_std() ---------------------------------------------------------
**
**    save_pictures_3d_std() saves more pictures in IPB 3D file format.
**
**    Usage ...
**   file :       nom du fichier
**   start :      numero de la premier image a sauver 
**   num_images : nombre d'images a sauver 
**   step :       pas si plusieurs images
**
**    Rajoute pour eviter dans les automates les question d'exostence
*/
int save_pictures_3d_std(char *file, int start, int num_images, int step)
{
  int err;
  
  err=save_mri_ipb_3d(file, start, num_images, step);

  return(err);
}

/* --save_pictures_3d_erase() ---------------------------------------------------------
**
**    save_pictures_3d_erase() saves more pictures in IPB 3D file format.
**
**    Usage ...
**   file :       nom du fichier
**   start :      numero de la premier image a sauver 
**   num_images : nombre d'images a sauver 
**   step :       pas si plusieurs images
**
**    Rajoute pour eviter dans les automates les question d'exostence. Ecrase automatiquement les fichiers existants
*/
int save_pictures_3d_erase(char *file, int start, int num_images, int step)
{
  //Ecrase le fichier si il existe : 
  char file_without_ext[PATH_LEN];
  char headerfile[PATH_LEN];
  char datafile[PATH_LEN];
  char s[64];
  struct stat statbuff;

  remove_file_extension(file,file_without_ext);
  strcpy(headerfile, file_without_ext);
  strcpy(datafile, file_without_ext);
  strcat(headerfile, ".ipb");
  strcat(datafile, ".img");

  if(stat(headerfile, &statbuff) == 0) 
    { 
      /* On supprime les fichiers .ipb, et .img */

      if(remove(headerfile) == -1) 
	{
	  sprintf(s, "[save_header] unable to remove file \"%s\"\n",
		  headerfile);	 
	  PUT_ERROR( s);
	  return(1); /* Erreur */
	}  
         
      if(remove(datafile) == -1) 
	{
	  sprintf(s,"[save_header] unable to remove file \"%s\"\n",
		  datafile);
	  PUT_ERROR( s);
	  return(1); /* Erreur */
	}	 		 
    }
  //Fin de : Ecrase le fichier si il existe  
  
  //sauvegarde du fichier
  int err;
  err=save_mri_ipb_3d(file, start, num_images, step);

  return(err);
}

/* --save_mri_ipb_3d() ---------------------------------------------------------
**
**    save_mri_ipb_3d() saves a picture in IPB 3D file format.
**
**    Usage ...
**   file :       nom du fichier
**   start :      numero de la premier image a sauver 
**   num_images : nombre d'images a sauver 
**   step :       pas si plusieurs images
**
**   return : 1 if Ok  !=1 if error
*/
int save_mri_ipb_3d(char *file, int start, int num_images, int step)
{
  grphic3d *image;
  int end;
  char s[256];
    
  end = (num_images - 1) * step + start;
  if(end > MAX_PICTURES_3D) {
    end = MAX_PICTURES_3D;
    PUT_WARN("[save_mri_ipb_3d] warning: cannot save all images !\n");
  }
  
  while(start <= end) {
    if((image = ptr_img_3d(start)) == NULL) {
      sprintf(s,"[save_mri_ipb_3d] unable to get image %d info!\n", start);
      PUT_ERROR(s);
      return(0);
    }
  
    if ( save_mri_ipb_3d_p(file, image) == 0) 
      return(0) ;
    start += step;
  }

  return(1);
        
}

#endif /* __GTK_GUI */
/* --save_mri_ipb_3d_p() ---------------------------------------------------------
**
**    save_mri_ipb_3d_p() saves a picture in IPB 3D file format.
**
**    Usage ...
**   file :       nom du fichier
**   image : pointeur sur struicture grphic3d 
**
**   return : 0 if Ok  !=0 if error
*/
int save_mri_ipb_3d_p(char *file, grphic3d *image)
{
  char s[256];
  int data_size;
  char headerfile[PATH_LEN];
  char datafile[PATH_LEN];
  struct stat statbuff;
  int index;
  long entier;
  double reel;
  long offset;
  char *buffchar;
  short *buffshort;
  long *bufflong;
  void *buffer;
  int i, j, k, num_items;
  FILE *fp;
  HEADER header;

  data_size = sizeof(TYPEMRI3D);

  remove_file_extension(file,headerfile);
  remove_file_extension(file,datafile);
  put_file_extension(headerfile,".ipb",headerfile);
  put_file_extension(datafile,".img",datafile);

  index = 0;
  offset = 0;

  if(stat(headerfile, &statbuff) == 0) { /* File exists */
    if((index = getmri_numberimg(headerfile)) == 0) {
      fprintf(stderr, "[save_mri_ipb_3d] Error reading header (1)!\n");
      return(0);
    }

    for(i=0; i<index; i++) {
      int width, height, depth, data_size;

      width = getmri_width(headerfile, i+1, MODE_3D);
      height = getmri_height(headerfile, i+1, MODE_3D);
      depth = getmri_depth(headerfile, i+1, MODE_3D);
      data_size = getmri_dattype(headerfile, i+1, MODE_3D);
      offset += width * height * depth * data_size;
    }	
  }

  index ++;

  if((fp = fopen(datafile, "ab")) == NULL) {
    sprintf(s,"Unable to save file %s !\n", file);
    PUT_ERROR(s);
    return(0);
  }

  if(load_header(headerfile, &header) == 0) {
    PUT_ERROR( "[save_mri_ipb_3d] Error reading header (3)!\n");
    fclose(fp);
    return(0);
  }

  put_header(&header, "IPB3", STRING, "", 0);
  put_header(&header, "patient name=",  STRING, image->patient.name,     0);
  put_header(&header, "date of birth=", STRING, image->patient.d_birth,  0);
  put_header(&header, "doctor name=",   STRING, image->patient.medecin,  0);
  put_header(&header, "date of exam=",  STRING, image->patient.d_examen, 0);
  put_header(&header, "number of images=", ENTIER, &index, 0);
#ifdef __GTK_GUI
  entier = (image->cmapinfo.cmap ? image->cmapinfo.cmap->noMode : 0);
  put_header(&header, "colormap=", ENTIER, &entier, index);
  put_header(&header, "display style=", ENTIER, &image->cmapinfo.dspStyle, index);
  put_header(&header, "overblancking=", ENTIER, &image->cmapinfo.dspOverBlancking, index);
#endif /* __GTK_GUI */
  entier = image->width;	
  put_header(&header, "width=", ENTIER, &entier, index);
  entier = image->height;	
  put_header(&header, "height=", ENTIER, &entier, index);
  entier = image->depth;	
  put_header(&header, "depth=", ENTIER, &entier, index);
  reel   = image->rcoeff;	
  put_header(&header, "rcoeff=", REEL, &reel, index);
  reel   = image->icomp;	
  put_header(&header, "icomp=", REEL, &reel, index);
  entier = image->max_pixel;  
  put_header(&header, "max pixel=", ENTIER, &entier, index);
  entier = image->min_pixel;  
  put_header(&header, "min pixel=", ENTIER, &entier, index);
  entier = image->cutoff_max; 
  put_header(&header, "cutoff max=", ENTIER, &entier, index);
  entier = image->cutoff_min; 
  put_header(&header, "cutoff min=", ENTIER, &entier, index);
  entier = data_size * 8;	
  put_header(&header, "bits per pixel=", ENTIER, &entier, index);
  reel   = image->dx; 	
  put_header(&header, "dx=", REEL, &reel, index);
  reel   = image->dy; 	
  put_header(&header, "dy=", REEL, &reel, index);
  reel   = image->dz; 	
  put_header(&header, "dz=", REEL, &reel, index);
  entier = image->x3dc;	
  put_header(&header, "x3dc=", ENTIER, &entier, index);
  entier = image->y3dc;	
  put_header(&header, "y3dc=", ENTIER, &entier, index);
  entier = image->z3dc;	
  put_header(&header, "z3dc=", ENTIER, &entier, index);
  entier = image->zoom;	
  put_header(&header, "zoom=", ENTIER, &entier, index);
  entier = image->zoom_x;	
  put_header(&header, "zoom_x=", ENTIER, &entier, index);
  entier = image->zoom_y;	
  put_header(&header, "zoom_y=", ENTIER, &entier, index);
  entier = image->zoom_z;	
  put_header(&header, "zoom_z=", ENTIER, &entier, index);
  entier = offset;		
  put_header(&header, "offset in raw=", ENTIER, &entier, index);

  put_header(&header, "--- Atlas de Talairach ---", STRING, "", 0);

  reel   = image->tala.AC.x;  
  put_header(&header, "AC_x=", REEL, &reel, index);
  reel   = image->tala.AC.y;  
  put_header(&header, "AC_y=", REEL, &reel, index);
  reel   = image->tala.AC.z;  
  put_header(&header, "AC_z=", REEL, &reel, index);
  reel   = image->tala.PC.x;  
  put_header(&header, "PC_x=", REEL, &reel, index);
  reel   = image->tala.PC.y;  
  put_header(&header, "PC_y=", REEL, &reel, index);
  reel   = image->tala.PC.z;  
  put_header(&header, "PC_z=", REEL, &reel, index);
  reel   = image->tala.AP.x;  
  put_header(&header, "AP_x=", REEL, &reel, index);
  reel   = image->tala.AP.y;  
  put_header(&header, "AP_y=", REEL, &reel, index);
  reel   = image->tala.AP.z;  
  put_header(&header, "AP_z=", REEL, &reel, index);
  reel   = image->tala.PP.x;  
  put_header(&header, "PP_x=", REEL, &reel, index);
  reel   = image->tala.PP.y;  
  put_header(&header, "PP_y=", REEL, &reel, index);
  reel   = image->tala.PP.z;  
  put_header(&header, "PP_z=", REEL, &reel, index);
  reel   = image->tala.RM.x;  
  put_header(&header, "RM_x=", REEL, &reel, index);
  reel   = image->tala.RM.y;  
  put_header(&header, "RM_y=", REEL, &reel, index);
  reel   = image->tala.RM.z;  
  put_header(&header, "RM_z=", REEL, &reel, index);
  reel   = image->tala.LM.x;  
  put_header(&header, "LM_x=", REEL, &reel, index);
  reel   = image->tala.LM.y;  
  put_header(&header, "LM_y=", REEL, &reel, index);
  reel   = image->tala.LM.z;  
  put_header(&header, "LM_z=", REEL, &reel, index);
  reel   = image->tala.SP.x;  
  put_header(&header, "SP_x=", REEL, &reel, index);
  reel   = image->tala.SP.y;  
  put_header(&header, "SP_y=", REEL, &reel, index);
  reel   = image->tala.SP.z;  
  put_header(&header, "SP_z=", REEL, &reel, index);
  reel   = image->tala.IP.x;  
  put_header(&header, "IP_x=", REEL, &reel, index);
  reel   = image->tala.IP.y;  
  put_header(&header, "IP_y=", REEL, &reel, index);
  reel   = image->tala.IP.z;  
  put_header(&header, "IP_z=", REEL, &reel, index);
  reel   = image->tala.IC.x;  
  put_header(&header, "IC_x=", REEL, &reel, index);
  reel   = image->tala.IC.y;  
  put_header(&header, "IC_y=", REEL, &reel, index);
  reel   = image->tala.IC.z;  
  put_header(&header, "IC_z=", REEL, &reel, index);

  reel   = image->tala.Tu.x;  
  put_header(&header, "Tu_x=", REEL, &reel, index);
  reel   = image->tala.Tu.y;  
  put_header(&header, "Tu_y=", REEL, &reel, index);
  reel   = image->tala.Tu.z;  
  put_header(&header, "Tu_z=", REEL, &reel, index);
  reel   = image->tala.Tv.x;  
  put_header(&header, "Tv_x=", REEL, &reel, index);
  reel   = image->tala.Tv.y;  
  put_header(&header, "Tv_y=", REEL, &reel, index);
  reel   = image->tala.Tv.z;  
  put_header(&header, "Tv_z=", REEL, &reel, index);
  reel   = image->tala.Tw.x;  
  put_header(&header, "Tw_x=", REEL, &reel, index);
  reel   = image->tala.Tw.y;  
  put_header(&header, "Tw_y=", REEL, &reel, index);
  reel   = image->tala.Tw.z;  
  put_header(&header, "Tw_z=", REEL, &reel, index);

  reel   = image->tala.x_AP;  
  put_header(&header, "x_AP=", REEL, &reel, index);
  reel   = image->tala.x_PC;  
  put_header(&header, "x_PC=", REEL, &reel, index);
  reel   = image->tala.x_PP;  
  put_header(&header, "x_PP=", REEL, &reel, index);
  reel   = image->tala.z_IP;  
  put_header(&header, "z_IP=", REEL, &reel, index);
  reel   = image->tala.z_SP;  
  put_header(&header, "z_SP=", REEL, &reel, index);
  reel   = image->tala.y_RM;  
  put_header(&header, "y_RM=", REEL, &reel, index);
  reel   = image->tala.y_LM;  
  put_header(&header, "y_LM=", REEL, &reel, index);

  put_header(&header, "--- DTI ---", STRING, "", 0);

	entier = image->dti.is_dti;	
  put_header(&header, "is_dti=", ENTIER, &entier, index);
	entier = image->dti.nb_directions;	
  put_header(&header, "nb_directions=", ENTIER, &entier, index);
	reel   = image->dti.gx;  
  put_header(&header, "gx=", REEL, &reel, index);
	reel   = image->dti.gy;  
  put_header(&header, "gy=", REEL, &reel, index);
	reel   = image->dti.gz;  
  put_header(&header, "gz=", REEL, &reel, index);
	reel   = image->dti.b;  
  put_header(&header, "b=", REEL, &reel, index);

#ifndef COMPILE_FOR_MEDIPY
 if(image->annotations)
    {
      write_annotations(image,&header,index);
    }
#endif

#ifdef __GTK_GUI
  imx_updcurv(&header, index, file, image);
#endif /* __GTK_GUI */

  offset += image->width * image->height * image->depth * data_size;

  num_items = image->width * image->height;

  switch(data_size)
    {
    case sizeof(char):
      buffchar = (char *)malloc(num_items);
      if(buffchar == NULL) {
	PUT_ERROR(
		  "[save_mri_ipb_3d] memory allocation error (2)!\n");
	fclose(fp);
	return(0);
      }
      buffer = buffchar;

      for(k=0; k<image->depth; k++)
	{
	  buffchar = buffer;

	  for(j=0; j<image->height; j++)
	    for(i=0; i<image->width; i++)
	      *buffchar ++ = image->mri[i][j][k];

	  if(fwrite(buffer, data_size, num_items, fp) != num_items) {
	    PUT_ERROR(
		      "[save_mri_ipb_3d] error writing data (2)!\n");
	    fclose(fp);
	    return(0);
	  }
	}
      break;

    case sizeof(short):
      buffshort = (short *)malloc(num_items * sizeof(short));
      if(buffshort == NULL) {
	PUT_ERROR(
		  "[save_mri_ipb_3d] memory allocation error (4)!\n");
	fclose(fp);
	return(0);
      }
      buffer = buffshort;

      for(k=0; k<image->depth; k++)
	{
	  buffshort = buffer;

	  for(j=0; j<image->height; j++)
	    for(i=0; i<image->width; i++)
	      //		    if (_is_bigEndian)
              *buffshort ++ = image->mri[i][j][k];
	  //			else
	  //              *buffshort ++ = shEndianConversion(image->mri[i][j][k]);

	  if(fwrite(buffer, data_size, num_items, fp) != num_items) {
	    PUT_ERROR(
		      "[save_mri_ipb_3d] error writing data (4)!\n");
	    fclose(fp);
	    return(0);
	  }
	}
      break;

    case sizeof(long):
      bufflong = (long *)malloc(num_items * sizeof(long));
      if(bufflong == NULL) {
	PUT_ERROR(
		  "[save_mri_ipb_3d] memory allocation error (6)!\n");
	fclose(fp);
	return(0);
      }

      buffer = bufflong;

      for(k=0; k<image->depth; k++)
	{
	  bufflong = buffer;

	  for(j=0; j<image->height; j++)
	    for(i=0; i<image->width; i++)
	      //		    if (_is_bigEndian)
              *bufflong ++ = image->mri[i][j][k];
	  //			else
	  //              *bufflong ++ = longEndianConversion(image->mri[i][j][k]);

	  if(fwrite(buffer, data_size, num_items, fp) != num_items) {
	    PUT_ERROR(
		      "[save_mri_ipb_3d] error writing data (6)!\n");
	    fclose(fp);
	    return(0);
	  }
	}
      break;

    default:
      PUT_ERROR(
		"[save_mri_ipb_3d] invalid data size\n");
      fclose(fp);
      return(0);
    }

  free(buffer);

  index ++;

  entier = index - 1;
  put_header(&header, "number of images=", ENTIER, &entier, 0);

  if (_is_bigEndian)
    put_header(&header, "endian=", STRING, "bigEndian", 0);
  else
    put_header(&header, "endian=", STRING, "littleEndian", 0);

  if(save_header(&header, headerfile) == 0) {
    sprintf(s,"Unable to save file header %s !\n", file);
    PUT_ERROR(s);
    fclose(fp);
    return(0);
  }

  fclose(fp);
  return(1);
}

/******************************************************************************
 ** -- update_tala() -----------------------------------------------------------
 **
 **	update talairach data in a picture
 **
 **	>>	Last user action by: F. BOULEAU on 01st Dec., 1999
 **		- bug fixed: saving position in file
 **
 ******************************************************************************/

int update_tala(int pos, int index, const char *file)
{
  HEADER header;
  grphic3d  *image;
  double reel;
  char s[64];

  if(load_header(file, &header) == 0) {
    PUT_ERROR( "[update_tala] Error reading header (1)!\n");
    return(0);
  }

  if((image = ptr_img_3d(pos)) == NULL) {
    sprintf(s,"[update_tala] unable to get image %d info!\n", pos);
    PUT_ERROR(s);
    return(0);
  }

  put_header(&header, "--- Atlas de Talairach ---", STRING, "", 0);
  reel   = image->tala.AC.x;  
  put_header(&header, "AC_x=", REEL, &reel, index);
  reel   = image->tala.AC.y;  
  put_header(&header, "AC_y=", REEL, &reel, index);
  reel   = image->tala.AC.z;  
  put_header(&header, "AC_z=", REEL, &reel, index);
  reel   = image->tala.PC.x;  
  put_header(&header, "PC_x=", REEL, &reel, index);
  reel   = image->tala.PC.y;  
  put_header(&header, "PC_y=", REEL, &reel, index);
  reel   = image->tala.PC.z;  
  put_header(&header, "PC_z=", REEL, &reel, index);
  reel   = image->tala.AP.x;  
  put_header(&header, "AP_x=", REEL, &reel, index);
  reel   = image->tala.AP.y;  
  put_header(&header, "AP_y=", REEL, &reel, index);
  reel   = image->tala.AP.z;  
  put_header(&header, "AP_z=", REEL, &reel, index);
  reel   = image->tala.PP.x;  
  put_header(&header, "PP_x=", REEL, &reel, index);
  reel   = image->tala.PP.y;  
  put_header(&header, "PP_y=", REEL, &reel, index);
  reel   = image->tala.PP.z;  
  put_header(&header, "PP_z=", REEL, &reel, index);
  reel   = image->tala.RM.x;  
  put_header(&header, "RM_x=", REEL, &reel, index);
  reel   = image->tala.RM.y;  
  put_header(&header, "RM_y=", REEL, &reel, index);
  reel   = image->tala.RM.z;  
  put_header(&header, "RM_z=", REEL, &reel, index);
  reel   = image->tala.LM.x;  
  put_header(&header, "LM_x=", REEL, &reel, index);
  reel   = image->tala.LM.y;  
  put_header(&header, "LM_y=", REEL, &reel, index);
  reel   = image->tala.LM.z;  
  put_header(&header, "LM_z=", REEL, &reel, index);
  reel   = image->tala.SP.x;  
  put_header(&header, "SP_x=", REEL, &reel, index);
  reel   = image->tala.SP.y;  
  put_header(&header, "SP_y=", REEL, &reel, index);
  reel   = image->tala.SP.z;  
  put_header(&header, "SP_z=", REEL, &reel, index);
  reel   = image->tala.IP.x;  
  put_header(&header, "IP_x=", REEL, &reel, index);
  reel   = image->tala.IP.y;  
  put_header(&header, "IP_y=", REEL, &reel, index);
  reel   = image->tala.IP.z;  
  put_header(&header, "IP_z=", REEL, &reel, index);
  reel   = image->tala.IC.x;  
  put_header(&header, "IC_x=", REEL, &reel, index);
  reel   = image->tala.IC.y;  
  put_header(&header, "IC_y=", REEL, &reel, index);
  reel   = image->tala.IC.z;  
  put_header(&header, "IC_z=", REEL, &reel, index);
  reel   = image->tala.Tu.x;  
  put_header(&header, "Tu_x=", REEL, &reel, index);
  reel   = image->tala.Tu.y;  
  put_header(&header, "Tu_y=", REEL, &reel, index);
  reel   = image->tala.Tu.z;  
  put_header(&header, "Tu_z=", REEL, &reel, index);
  reel   = image->tala.Tv.x;  
  put_header(&header, "Tv_x=", REEL, &reel, index);
  reel   = image->tala.Tv.y;  
  put_header(&header, "Tv_y=", REEL, &reel, index);
  reel   = image->tala.Tv.z;  
  put_header(&header, "Tv_z=", REEL, &reel, index);
  reel   = image->tala.Tw.x;  
  put_header(&header, "Tw_x=", REEL, &reel, index);
  reel   = image->tala.Tw.y;  
  put_header(&header, "Tw_y=", REEL, &reel, index);
  reel   = image->tala.Tw.z;  
  put_header(&header, "Tw_z=", REEL, &reel, index);
  reel   = image->tala.x_AP;  
  put_header(&header, "x_AP=", REEL, &reel, index);
  reel   = image->tala.x_PC;  
  put_header(&header, "x_PC=", REEL, &reel, index);
  reel   = image->tala.x_PP;  
  put_header(&header, "x_PP=", REEL, &reel, index);
  reel   = image->tala.z_IP;  
  put_header(&header, "z_IP=", REEL, &reel, index);
  reel   = image->tala.z_SP;  
  put_header(&header, "z_SP=", REEL, &reel, index);
  reel   = image->tala.y_RM;  
  put_header(&header, "y_RM=", REEL, &reel, index);
  reel   = image->tala.y_LM;  
  put_header(&header, "y_LM=", REEL, &reel, index);

  if(save_header(&header, file) == 0) {
    sprintf(s,"Unable to save file %s !\n", file);
    PUT_ERROR(s);
    return(0);
  }

  return(1);
}


/* --buff3pic() -------------------------------------------------------    
**    Transforms a buffer 
**      to a picture[MAX_WIDTH_3D][MAX_HEIGHT_3D][MAX_DEPTH_3D].
**
*/
long    buff3pic(long int *buff, grphic3d *graphic, int dpth, float anc_rcoef)
{
  long i,j,k;
  int wdth,hght;
  double coef;
  int zoom,tmp;    

  wdth=graphic->width;
  hght=graphic->height;
  k=j=i=(long)0;
  coef=((double) anc_rcoef)/((double) graphic->rcoeff);

  tmp=MAXI(wdth,hght);
  if (tmp>MAX_WIDTH_3D || tmp>MAX_HEIGHT_3D)
    {
      zoom=ceil((float)tmp/((float)MAXI(MAX_WIDTH_3D,MAX_HEIGHT_3D)));
      for(i=0;i<wdth/zoom;i++)
	for(j=0;j<hght/zoom;j++){	
	  graphic->mri[i][j][dpth]=buff[i*wdth*zoom+j*zoom]*coef;
	  k++;
	}			 
    }
  else
    {
      for(i=0;i<wdth;i++)
	for(j=0;j<hght;j++){
	  graphic->mri[i][j][dpth]=(short)buff[k++]*coef;
	}
    }


  return(k);
}

#ifdef __GTK_GUI
/**********************************************************
 *** int fill_mri_8b()
 **	
 **	remplir une structure mri avec des donnee 8bits suivant une orientation
 **  a partir d'un buffer 8b et d'une position dans le buffer
 **
 **********************************************************/
int fill_mri_8b(grphic3d *img, 
		unsigned char *buff8b, 
		int i, int j, int k, int pos,float coeff,
		int read_format)
{
  int wdth,hght,dpth;

  wdth = img->width;
  hght = img->height;
  dpth = img->depth;
	
  switch (read_format) 
    {
    case 0: /* On ne change rien */
      img->mri[i][j][k]=coeff*buff8b[pos];
      break;
    case 1:  /* Coupe transversale par le sommet du crane */
      img->mri[i][k][j]=coeff*buff8b[pos];
      break;
    case 2:   /* IPB et paravision sagittal */
      img->mri[k][j][i]=coeff*buff8b[pos];
      break;
    case 3:     /*  Roufach 1 */
      img->mri[j][i][k]=coeff*buff8b[pos];
      break;
    case 4:  /*  MN dicom */
      img->mri[k][j][wdth-i-1]=coeff*buff8b[pos];
      break;
    case 6: /* Coupe transversale par le cou */
      img->mri[i][dpth-k-1][j]=coeff*buff8b[pos];
      break;
    case 7:  /* pour ECAT7 */
      img->mri[wdth-i-1][k][j]=coeff*buff8b[pos];						
      break;
    case 8:  /* pour  ANALYZE transverse unflipped */
      img->mri[i][k][hght-1-j]=coeff*buff8b[pos];
      break;
    case 9:  /*  */
      img->mri[j][k][i]=coeff*buff8b[pos];
      break;
    case 10:  /* */
      img->mri[k][i][j]=coeff*buff8b[pos];
      break;
    case 11:  /* */
      img->mri[i][j][dpth-k]=coeff*buff8b[pos];
      break;
    case 12:  /* */
      img->mri[i][dpth-k-1][hght-j-1]=coeff*buff8b[pos];
      break;
    default : /* On ne change rien */
      img->mri[i][j][k]=coeff*buff8b[pos];
      break;
    }
  return 1;
}

/**********************************************************
 *** int fill_mri_16b()
 **	
 **	remplir une structure mri avec des donnee 8bits suivant une orientation
 **  a partir d'un buffer 8b et d'une position dans le buffer
 **
 **********************************************************/
int fill_mri_16b(grphic3d *img,
		 short *buff16b, 
		 int i, int j, int k, int pos, float coeff,
		 int read_format,
		 char is_filebigEndian)
{
  int wdth,hght,dpth;
	
  wdth = img->width;
  hght = img->height;
  dpth = img->depth;
	

		
  switch (read_format) 
    {
    case 0: /* On ne change rien */
      if(is_filebigEndian)
	img->mri[i][j][k]=coeff*buff16b[pos];
      else
	img->mri[i][j][k] = (TYPEMRI3D) shEndianConversion(buff16b[pos]);   
      break;
    case 1:  /* Coupe transversale par le sommet du crane */
      if(is_filebigEndian)
	img->mri[i][k][j]=coeff*buff16b[pos];
      else
	img->mri[i][k][j] = coeff* (TYPEMRI3D) shEndianConversion(buff16b[pos]);   
      break;
    case 2:  /*  IPB, paravision sagital, et siemens Vanderbilt */
      /* Attention pour scheiber 3/99 modif k-dpth-1 */
      if(is_filebigEndian)
	img->mri[k][j][i]=coeff*buff16b[pos];
      else
	img->mri[k][j][i] = (TYPEMRI3D) shEndianConversion(buff16b[pos]);   
      break;
    case 3:  /*  Roufach 1 */
      if(is_filebigEndian)
	img->mri[j][i][k]=coeff*buff16b[pos];
      else
	img->mri[j][i][k] = (TYPEMRI3D) shEndianConversion(buff16b[pos]);   
      break;
    case 4:  /*  MN dicom */
      if(is_filebigEndian)
	img->mri[k][j][wdth-i-1]=coeff*buff16b[pos];
      else
	img->mri[k][j][wdth-i-1] = (TYPEMRI3D) shEndianConversion(buff16b[pos]);   
      break;
    case 6:  /* Coupe transversale par le cou */
      if(is_filebigEndian)
	img->mri[i][dpth-k-1][j]=coeff*buff16b[pos];
      else
	img->mri[i][dpth-k-1][j] = (TYPEMRI3D) shEndianConversion(buff16b[pos]);   
      break;
    case 7:  /* pour  ECAT7 */
      if(is_filebigEndian)
	img->mri[wdth-i-1][k][j]=coeff*buff16b[pos];
      else
	img->mri[wdth-i-1][k][j] = (TYPEMRI3D) shEndianConversion(buff16b[pos]);   
      break;
    case 8:  /* pour  ANALYZE transverse unflipped */
      if(is_filebigEndian)
	img->mri[i][k][hght-1-j]=coeff*buff16b[pos];
      else
	img->mri[i][k][hght-1-j] = (TYPEMRI3D) shEndianConversion(buff16b[pos]);   
      break;
    case 9:  /*  */
      if(is_filebigEndian)
	img->mri[j][k][i]=coeff*buff16b[pos];
      else
	img->mri[j][k][i] = (TYPEMRI3D) shEndianConversion(buff16b[pos]);   
      break;
    case 10:  /* */
      if(is_filebigEndian)
	img->mri[k][i][j]=coeff*buff16b[pos];
      else
	img->mri[k][i][j] = (TYPEMRI3D) shEndianConversion(buff16b[pos]);   
      break;
    case 11:  /* */
      if(is_filebigEndian)
	img->mri[i][j][dpth-k]=coeff*buff16b[pos];
      else
	img->mri[i][j][dpth-k] = (TYPEMRI3D) shEndianConversion(buff16b[pos]);   
      break;
    case 12:/*  */
      if(is_filebigEndian)
	img->mri[i][dpth-k-1][hght-j-1]=coeff*buff16b[pos];
      else
	img->mri[i][dpth-k-1][hght-j-1] = (TYPEMRI3D) shEndianConversion(buff16b[pos]);
      break;
    default: /* On ne change rien */	
      if(is_filebigEndian)
	img->mri[i][j][k]=coeff*buff16b[pos];
      else
	img->mri[i][j][k] = (TYPEMRI3D) shEndianConversion(buff16b[pos]);   
      break;
    }
  return 1;
}

int fill_mri_16b2(grphic3d *img,
		 short *buff16b,
		 int i, int j, int k, int pos, float coeff,
		 int read_format,
		 char is_filebigEndian)
{
  int wdth,hght,dpth;
  TYPEMRI3D valPix=0;

  wdth = img->width;
  hght = img->height;
  dpth = img->depth;

  if (is_filebigEndian) valPix=(TYPEMRI3D)floor(coeff*shEndianConversion(buff16b[pos])+0.5);
  else valPix=(TYPEMRI3D)buff16b[pos];
  //valPix=(TYPEMRI3D)floor(coeff*buff16b[pos]+0.5);

  switch (read_format)
    {
    case 0: // On ne change rien
      img->mri[i][j][k]=valPix;
      break;
    case 1:  // Coupe transversale par le sommet du crane 
      img->mri[i][k][j]=valPix;
      break;
    case 2:  //  IPB, paravision sagital, et siemens Vanderbilt 
      // Attention pour scheiber 3/99 modif k-dpth-1 
      img->mri[k][j][i]=valPix;
      break;
    case 3:  //  Roufach 1 
      img->mri[j][i][k]=valPix;
      break;
    case 4:  //  MN dicom 
      img->mri[k][j][wdth-i-1]=valPix;
      break;
    case 6:  // Coupe transversale par le cou 
      img->mri[i][dpth-k-1][j]=valPix;
      break;
    case 7:  // pour  ECAT7 
      img->mri[wdth-i-1][k][j]=valPix;
      break;
    case 8:  // pour  ANALYZE transverse unflipped 
      img->mri[i][k][hght-1-j]=valPix;
      break;
    case 9:  //  
      img->mri[j][k][i]=valPix;
      break;
    case 10:  // 
      img->mri[k][i][j]=valPix;
      break;
    case 11:  // 
      img->mri[i][j][dpth-k]=valPix;
      break;
    case 12://  
      img->mri[i][dpth-k-1][hght-j-1]=valPix;
      break;
    default: // On ne change rien 
      img->mri[i][j][k]=valPix;
      break;
    }
  return 1;
}

/**********************************************************
 *** int fill_mri_32b()
 **	
 **	remplir une structure mri avec des donnee 8bits suivant une orientation
 **  a partir d'un buffer 8b et d'une position dans le buffer
 **
 **********************************************************/
int fill_mri_32b(grphic3d *img, 
		 long *buff32b, 
		 int i, int j, int k, int pos, float coeff,
		 int read_format, 
		 char is_filebigEndian)
{
  int wdth,hght,dpth;
	
  wdth = img->width;
  hght = img->height;
  dpth = img->depth;
	
  //	if (coeff*buff32b[pos]!=0)
  //	 printf ("info %d %d %d %d\n",i,j,k,read_format);
	 
  switch (read_format) 
    {
    case 0: /* On ne change rien */
      img->mri[i][j][k]=coeff*buff32b[pos];
      break;
    case 1: /* Coupe transversale par le sommet du crane */
      img->mri[i][k][j]=coeff*buff32b[pos];
      break;
    case 2:  /*  IPB et paravision sagital */
      img->mri[k][j][i]=coeff*buff32b[pos];
      break;
    case 3: /* Roufach 1 */
      img->mri[j][i][k]=coeff*buff32b[pos];
      break;
    case 4:  /*  MN dicom */
      img->mri[k][j][wdth-i-1]=coeff*buff32b[pos];
      break;
    case 6:/* Coupe transversale par le cou */
      img->mri[i][dpth-k-1][j]=coeff*buff32b[pos];
      break;
    case 7:  /* pour ECAT7 */
      img->mri[wdth-i-1][k][j]=coeff*buff32b[pos];						
      break;
    case 8:  /* pour  ANALYZE transverse unflipped */
      img->mri[i][k][hght-1-j]=coeff*buff32b[pos];
      break;
    case 9:  /*  */
      img->mri[j][k][i]=coeff*buff32b[pos];
      break;
    case 10:  /* */
      img->mri[k][i][j]=coeff*buff32b[pos];
      break;
    case 11:  /* */
      img->mri[i][j][dpth-k]=coeff*buff32b[pos];
      break;
    case 12:/*  */
      img->mri[i][dpth-k-1][hght-j-1]=coeff*buff32b[pos];
      break;
    default: /* On ne change rien */
      img->mri[i][j][k]=coeff*buff32b[pos];
      break;
    }
  return 1;
}

int fill_mri_32b2(grphic3d *img,
		 long *buff32b,
		 int i, int j, int k, int pos, float coeff,
		 int read_format,
		 char is_filebigEndian)
{
  int wdth,hght,dpth;
  TYPEMRI3D valPix=0;

  wdth = img->width;
  hght = img->height;
  dpth = img->depth;

/*  if (is_filebigEndian) valPix=(TYPEMRI3D)floor(coeff*buff32b[pos]+0.5);
  else (TYPEMRI3D) valPix=(TYPEMRI3D)floor(coeff*longEndianConversion(buff32b[pos])+0.5); */
  if (is_filebigEndian) valPix=(TYPEMRI3D)floor(coeff*longEndianConversion(buff32b[pos])+0.5);
  else valPix=(TYPEMRI3D)floor(coeff*buff32b[pos]+0.5);

  //	if (coeff*buff32b[pos]!=0)
  //	 printf ("info %d %d %d %d\n",i,j,k,read_format);

  switch (read_format)
    {
    case 0: /* On ne change rien */
      img->mri[i][j][k]=valPix;
      break;
    case 1: /* Coupe transversale par le sommet du crane */
      img->mri[i][k][j]=valPix;
      break;
    case 2:  /*  IPB et paravision sagital */
      img->mri[k][j][i]=valPix;
      break;
    case 3: /* Roufach 1 */
      img->mri[j][i][k]=valPix;
      break;
    case 4:  /*  MN dicom */
      img->mri[k][j][wdth-i-1]=valPix;
      break;
    case 6:/* Coupe transversale par le cou */
      img->mri[i][dpth-k-1][j]=valPix;
      break;
    case 7:  /* pour ECAT7 */
      img->mri[wdth-i-1][k][j]=valPix;
      break;
    case 8:  /* pour  ANALYZE transverse unflipped */
      img->mri[i][k][hght-1-j]=valPix;
      break;
    case 9:  /*  */
      img->mri[j][k][i]=valPix;
      break;
    case 10:  /* */
      img->mri[k][i][j]=valPix;
      break;
    case 11:  /* */
      img->mri[i][j][dpth-k]=valPix;
      break;
    case 12:/*  */
      img->mri[i][dpth-k-1][hght-j-1]=valPix;
      break;
    default: /* On ne change rien */
      img->mri[i][j][k]=valPix;
      break;
    }
  return 1;
}

/*********************************************************
 **
 ** find_whd_from_orient_p()
 **
 ** met les bonnes valeurs de img->width, img->height etc
 ** en fonction de l'orientation 
 **
 **
 ********************************************************/

int find_whd_from_orient_p(grphic3d *img, int wdth, int hght, int dpth, int read_format_3d)
{
  switch (read_format_3d)
    {
    case 0:/*on echange rien*/
      img->width=wdth;
      img->height=hght;
      img->depth=dpth;
      break;
    case 1: /*on echange k et j*/
      img->width=wdth;
      img->height=dpth;
      img->depth=hght;
      SWAPF(img->dy,img->dz);
      SWAP(img->y3dc,img->z3dc);
      break;
    case 2:/*on echange k et i*/
      img->width=dpth;
      img->height=hght;
      img->depth=wdth;
      SWAPF(img->dx,img->dz);
      SWAP(img->x3dc,img->z3dc);
      break;
    case 3:/*on echange i et j*/
      img->width=hght;
      img->height=wdth;
      img->depth=dpth;
      SWAPF(img->dx,img->dy);
      SWAP(img->x3dc,img->y3dc);
      break;
    case 4:/*on echange k et i + miroir i*/
      img->width=dpth;
      img->height=hght;
      img->depth=wdth;
      SWAPF(img->dx,img->dz);
      SWAP(img->x3dc,img->z3dc);
      break;
    case 6: /*on echange k et j + miroir  j*/
      img->width=wdth;
      img->height=dpth;
      img->depth=hght;
      SWAPF(img->dy,img->dz);
      SWAP(img->y3dc,img->z3dc);
      break;
    case 7:/*on echange k et j + miroir i*/
      img->width=wdth;
      img->height=dpth;
      img->depth=hght;
      SWAPF(img->dy,img->dz);
      SWAP(img->y3dc,img->z3dc);
      break;
    case 8:/*on echange k et j + miroir k*/
      img->width=wdth;
      img->height=dpth;
      img->depth=hght;
      SWAPF(img->dy,img->dz);
      SWAP(img->y3dc,img->z3dc);
      break;
    case 9:/*on echange k et j + miroir k*/
      img->width=hght;
      img->height=dpth;
      img->depth=wdth;
      SWAPF(img->dx,img->dy);
      SWAPF(img->dy,img->dz);
      SWAP(img->x3dc,img->y3dc);
      SWAP(img->y3dc,img->z3dc);
      break;
    case 10:/*on echange k et j + miroir k*/
      img->width=dpth;
      img->height=wdth;
      img->depth=hght;
      SWAPF(img->dx,img->dz);
      SWAPF(img->dy,img->dz);
      SWAP(img->x3dc,img->z3dc);
      SWAP(img->y3dc,img->z3dc);
      break;
    case 11:/*on echange rien*/
      img->width=wdth;
      img->height=hght;
      img->depth=dpth;
      break;
    case 12:/*on echange k et j + miroir k + miroir j*/
      img->width=wdth;
      img->height=dpth;
      img->depth=hght;
      SWAPF(img->dy,img->dz);
      SWAP(img->y3dc,img->z3dc);
      break;
    default:/*on echange rien*/
      img->width=wdth;
      img->height=hght;
      img->depth=dpth;
      break;
    }
  return 1;
}

/*********************************************************
 **
 ** replace_whd_from_orient()
 **
 ** met les bonnes valeurs de img->width, img->height etc
 ** en fonction de l'orientation 
 **
 **
 ********************************************************/

int replace_whd_from_orient(int *wdth, int *hght, int *dpth, int read_format_3d)
{
  int tmp,tmp2;
	
  switch (read_format_3d)
    {
    case 0:/*on echange rien*/
      *wdth=*wdth;
      *hght=*hght;
      *dpth=*dpth;
      break;
    case 1: /*on echange k et j*/
      *wdth=*wdth;
      tmp  =*hght;
      *hght=*dpth;
      *dpth=tmp;
      break;
    case 2:/*on echange k et i*/
      tmp  =*wdth;
      *wdth=*dpth;
      *dpth=tmp;
      *hght=*hght;
      break;
    case 3:/*on echange i et j*/
      tmp  =*wdth;
      *wdth=*hght;
      *hght=tmp;
      *dpth=*dpth;
      break;
    case 4:/*on echange k et i + miroir i*/
      tmp  =*wdth;
      *wdth=*dpth;
      *hght=*hght;
      *dpth=tmp;
      break;
    case 6: /*on echange k et j + miroir  j*/
      *wdth=*wdth;
      tmp  =*hght;
      *hght=*dpth;
      *dpth=tmp;
      break;
    case 7:/*on echange k et j + miroir i*/
      *wdth=*wdth;
      tmp  =*hght;
      *hght=*dpth;
      *dpth=tmp;
      break;
    case 8:/*on echange k et j + miroir k*/
      *wdth=*wdth;
      tmp  =*hght;
      *hght=*dpth;
      *dpth=tmp;
      break;
    case 9:/*on echange i->j j->k k->i miroir k*/
      tmp  =*wdth;
      tmp2 =*dpth;
      *wdth=*hght;
      *dpth=tmp;
      *hght=tmp2;
      break;
    case 10:/*on echange k et j + miroir k*/
      tmp  =*wdth;
      tmp2 =*hght;
      *wdth=*dpth;
      *hght=tmp;
      *dpth=tmp2;
      break;
    case 11:/*on echange rien*/
      *wdth=*wdth;
      *hght=*hght;
      *dpth=*dpth;
      break;
    case 12:/*on echange k et j + miroir k + miroir j*/
      tmp  =*wdth;
      tmp2 =*hght;
      *wdth=*dpth;
      *hght=tmp;
      *dpth=tmp2;
      break;
    default:/*on echange rien*/
      *wdth=*wdth;
      *hght=*hght;
      *dpth=*dpth;
      break;
    }
  return 1;
}

/**********************************************************
 *** int read_raw_3d ()
 ***
 ***    Permet de lire une image 3D rawdata quelque soit 
 ***     taille de celle ci.
 ***
 ***     PARAMETRES
 ***        file : nom de l'image
 ***        numfirst : numero de l'image dans le fichier raw data
 ***        img  : pointeur de l'image        
 ***        dimx,dimy,dimz : dimension de l'image
 ***        bytesppixel : nombre d'octets par point (1byte, 2bytes etc)
 ***             bytesppixel doit etre inferieure a 5 !!!
 ***        format : resultat de la commande is_mri (type de fichier)
 ***    REMARQUE
 ***        La fonction teste la compatibilite entre 
 ***        dimx,dimy,dimz et taille d'une part et la taille
 ***        en octet de l'image d'autre part.
 ***
 ***             Il faut que graphix->max_pixel soit positionne
 ***********************************************************/

int read_raw_3d(const char *file, int numfirst, int filewdth, int filehght, int filedpth, grphic3d *img, int bytesppixel, int format)
{
  int i,j,k,l,size,err;
  long setposfile;
  float amax,coeff,amin;

#ifdef HAS_NIFTI
  znzFile fp;
#else
  FILE * fp;
#endif
  short *buff16b;
  long *buff32b;
  unsigned char *buff8b;
  int wdth=0,hght=0,dpth=0;
  int pas,pasw,pash,pasd;
  //  int tmp,dst;					unused variable
  //  char* m_int=(char*)malloc(2);	unused variable
  //  char is_filebigEndian=FALSE;
	
#ifdef HAS_NIFTI
  if(format == NIFTI_GZ_FORMAT){
	if((fp=znzopen(file, "rb", 1)) == NULL){
		printf("Le fichier %s n'existe pas \n",file);
    		PUT_MESG( "The file does'nt exist (read_raw_3d)\n");
    		return(-1);
	}
  }else{
  	if ( (fp=znzopen(file,"rb", 0))==NULL )  {
  	  printf("Le fichier %s n'existe pas \n",file);
  	  PUT_MESG( "The file does'nt exist (read_raw_3d)\n");
  	  return(-1);
  	}
  }
  znzseek(fp,0l,SEEK_END);
  size=znztell(fp);
#else
  if ( (fp=fopen(file,"rb"))==NULL )  {
  	  printf("Le fichier %s n'existe pas \n",file);
  	  PUT_MESG( "The file does'nt exist (read_raw_3d)\n");
  	  return(-1);
  }
  fseek(fp, 01, SEEK_END);
  size=ftell(fp);
#endif

  if(format != NIFTI_GZ_FORMAT){
	  if ( fp==NULL || size < numfirst*filewdth*filehght*filedpth*bytesppixel) {
	    PUT_MESG( "Cannot open the raw file or size too small\n");
	    return (-1);
	  }
  }

  /*  Recherche de l'offset en cas de header dans le fichier raw data */
  setposfile=getPosFile(file,numfirst,filewdth,filehght,filedpth,bytesppixel,format);

#ifdef HAS_NIFTI
  znzseek(fp,setposfile,SEEK_SET);
#else
  fseek(fp, setposfile, SEEK_SET);
#endif
	
  /*MISE EN PLACE DES MAX DANS LES IMAGES EN FONCTION DES bytesppixelS DE LECTURE*/
  switch (_read_format_3d)
    {
    case 0:  /*on echange rien*/
      wdth=MINI(filewdth,MAX_WIDTH_3D);
      hght=MINI(filehght,MAX_HEIGHT_3D);
      dpth=MINI(filedpth,MAX_DEPTH_3D);
      break;
    case 1:  /*on echange k et j*/
      wdth=MINI(filewdth,MAX_WIDTH_3D);
      hght=MINI(filehght,MAX_DEPTH_3D);
      dpth=MINI(filedpth,MAX_HEIGHT_3D);
      break;
    case 2:  /*on echange k et i*/
      wdth=MINI(filewdth,MAX_DEPTH_3D);
      hght=MINI(filehght,MAX_HEIGHT_3D);
      dpth=MINI(filedpth,MAX_WIDTH_3D);
      break;
    case 3:  /*on echange i et j*/
      wdth=MINI(filewdth,MAX_HEIGHT_3D);
      hght=MINI(filehght,MAX_WIDTH_3D);
      dpth=MINI(filedpth,MAX_DEPTH_3D);
      break;
    case 4:  /*on echange k et i et +miroir k*/
      wdth=MINI(filewdth,MAX_DEPTH_3D);
      hght=MINI(filehght,MAX_HEIGHT_3D);
      dpth=MINI(filedpth,MAX_WIDTH_3D);
      break;
    case 6:  /*on echange k et j + miroir  j*/
      wdth=MINI(filewdth,MAX_WIDTH_3D);
      hght=MINI(filehght,MAX_DEPTH_3D);
      dpth=MINI(filedpth,MAX_HEIGHT_3D);
      break;
    case 7:  /*on echange k et j + miroir i*/
      wdth=MINI(filewdth,MAX_WIDTH_3D);
      hght=MINI(filehght,MAX_DEPTH_3D);
      dpth=MINI(filedpth,MAX_HEIGHT_3D);
      break;
    case 8:  /*on echange k et j + miroir k*/
      wdth=MINI(filewdth,MAX_WIDTH_3D);
      hght=MINI(filehght,MAX_DEPTH_3D);
      dpth=MINI(filedpth,MAX_HEIGHT_3D);
      break;
    case 9:  /*i->k ; j->i ; k->j*/
      wdth=MINI(filewdth,MAX_HEIGHT_3D);
      hght=MINI(filehght,MAX_DEPTH_3D);
      dpth=MINI(filedpth,MAX_WIDTH_3D);
      break;
    case 10:  /*i->j ; j->k ; k->i*/
      wdth=MINI(filewdth,MAX_DEPTH_3D);
      hght=MINI(filehght,MAX_WIDTH_3D);
      dpth=MINI(filedpth,MAX_HEIGHT_3D);
      break;
    case 11:  /*i->i ; j->j ; k->k*/
      wdth=MINI(filewdth,MAX_WIDTH_3D);
      hght=MINI(filehght,MAX_HEIGHT_3D);
      dpth=MINI(filedpth,MAX_DEPTH_3D);
      break;
    case 12:  /*on echange k et j + miroir k*/
      wdth=MINI(filewdth,MAX_WIDTH_3D);
      hght=MINI(filehght,MAX_DEPTH_3D);
      dpth=MINI(filedpth,MAX_HEIGHT_3D);
      break;
    default :  /*on echange rien*/
      wdth=MINI(filewdth,MAX_WIDTH_3D);
      hght=MINI(filehght,MAX_HEIGHT_3D);
      dpth=MINI(filedpth,MAX_DEPTH_3D);
      break;
    }

  /*Calcul du pas de lecture si l'image est plus grande suivant une direction */
  /*cela permet un ajustement automatique de la taille de l'image a la place disponible a l'arrivee*/
  pasw=1;
  pash=1;
  pasd=1;
  if (filewdth>wdth || filehght>hght || filedpth>dpth)
    {
      if (wdth && filewdth%wdth==0) {
	pasw=filewdth/wdth;
      } else {
	pasw=filewdth/wdth+1;
      }
      if (pasw<1) {
	pasw=1;
      }
      if (hght && filehght%hght==0) {
	pash=filehght/hght;
      } else {
	pash=filehght/hght+1;
      }
      if (pash<1) {
	pash=1;
      }
      if (dpth && filedpth%dpth==0) {
	pasd=filedpth/dpth;
      } else {
	pasd=filedpth/dpth+1;
      }
      if (pasd<1) {
	pasd=1;
      }
    }
  pas=MAXI(pasw,pash);
  pas=MAXI(pas,pasd);
  if (pas>1) {
    printf("ATTENTION l'image est trop grande: le pas de lecture est %d \n",pas);
    PUT_WARN("L'image est trop grande.");
    img->dx *= pasw;
    img->dy *= pash;
    img->dz *= pasd;
  }
   
  wdth=filewdth/pasw;
  hght=filehght/pash;
  dpth=filedpth/pasd;
  /* Pour mettre tous les pas a la meme valeur et on reajuste les dimensions
     enlever le commentaire suivant */
  /*   Commentaire
       if (pasw!=pas) {
       pasw=pas;
       wdth=wdth/pas;
       }
       if (pash!=pas) {
       pash=pas;
       hght=hght/pas;
       }
       if (pasd!=pas) {
       pasd=pas;
       dpth=dpth/pas;
       }
       Fin commentaire */
  img->width = wdth;
  img->height= hght;
  img->depth = dpth;
  if (!img->mri)
    img->mri = cr_mri_3d(wdth, hght, dpth);
  switch(bytesppixel)
    {
    case 1: /*On lit un fichier 8 bits*/
      /*allocation du buffer a la bonne taille*/
      if((buff8b=CALLOC(filewdth*filehght,unsigned char))==NULL) { /* Cannot alloc... */
#ifdef HAS_NIFTI
	znzclose(fp);
#else
	fclose(fp);
#endif
	return(0);
      }

      /* Recherche maxpixel */
      coeff=1;
      for (k=0;k<dpth;k++) {
#ifdef HAS_NIFTI
	size=znzread((unsigned char*)buff8b,sizeof(unsigned char),filewdth*filehght,fp);
#else
	size=fread((unsigned char*)buff8b,sizeof(unsigned char),filewdth*filehght,fp);
#endif
	l=0;
	for (j=0;j<hght;j++)
	  {
	    for (i=0;i<wdth;i++) 
	      {
		fill_mri_8b(img, buff8b, i, j, k, l, coeff,_read_format_3d);
		l=l+pasw;
	      }

	    if( pash!=1  && filewdth%2)l++;
	    l=l+filewdth*(pash-1);
	  }
	for (i=0;i<pasd-1;i++){
#ifdef HAS_NIFTI
	size=znzread((unsigned char*)buff8b,sizeof(unsigned char),filewdth*filehght,fp);
#else
	size=fread((unsigned char*)buff8b,sizeof(unsigned char),filewdth*filehght,fp);
#endif
	}
      }
		

      free(buff8b);
      break;

    case 2: /*On lit un fichier 16 bits*/
      /*allocation du buffer a la bonne taille*/
      if((buff16b=CALLOC(filewdth*filehght,short))==NULL) { /* Cannot alloc... */
#ifdef HAS_NIFTI
	znzclose(fp);
#else
	fclose(fp);
#endif
	return(0);
      }

      /* Lecture de l'image */	

       if (sizeof(TYPEMRI3D) < 2 )
       {
        coeff=img->rcoeff;
        amax=img->rcoeff*img->max_pixel;
        amin=img->rcoeff*img->min_pixel;
        err=imx_brukermax_3d(amax,amin,img);
        coeff=coeff/img->rcoeff;
       }
       else
        coeff=1;
 
      for (k=0;k<dpth;k++)
      {
#ifdef HAS_NIFTI
	size=znzread((short*)buff16b,(size_t)sizeof(short),(size_t)filewdth*filehght,fp);
#else
	size=FREAD_SHORT2((short*)buff16b,filewdth*filehght,fp);
#endif
       l=0;

       for (j=0;j<hght;j++)
       {
        for (i=0;i<wdth;i++) 
        {
   //      fill_mri_16b(img, buff16b, i, j, k, l, coeff, _read_format_3d, _is_filebigEndian);
         fill_mri_16b2(img, buff16b, i, j, k, l, coeff, _read_format_3d, _is_filebigEndian);
         l=l+pasw;
        }
        l=l+filewdth*(pash-1);
       }
       /*
       for (i=0;i<pasd-1;i++)
       {
        size=FREAD_SHORT((short*)buff16b,filewdth*filehght,fp);
       }*/
#ifdef HAS_NIFTI
       znzseek(fp,(pasd-1)*filewdth*filehght*sizeof(short),SEEK_CUR);
#else
	fseek(fp,(pasd-1)*filewdth*filehght*sizeof(short),SEEK_CUR);
#endif
      }


      free(buff16b);
      break;

    case 4: /*On lit un fichier 32 bits*/
      /*allocation du buffer a la bonne taille*/
      if((buff32b=CALLOC(filewdth*filehght,long))==NULL) { /* Cannot alloc... */
#ifdef HAS_NIFTI
	znzclose(fp);
#else
	fclose(fp);
#endif
	return(0);
      }

      if (sizeof(TYPEMRI3D) < 4 ) 
	{
	  coeff=img->rcoeff;
	  amax=img->rcoeff*img->max_pixel;
	  amin=img->rcoeff*img->min_pixel;
	  err=imx_brukermax_3d(amax,amin,img);
	  coeff=coeff/img->rcoeff;
	}
      else
	coeff=1;

      /* Lecture de l'image */
      for (k=0;k<dpth;k++)
      {
#ifdef HAS_NIFTI
	size=znzread(buff32b,(size_t)sizeof(long),filewdth*filehght,fp);
#else
	size=FREAD_LONG2(buff32b,filewdth*filehght,fp);
#endif
       l=0;
       for (j=0;j<hght;j++)
       {
        for (i=0;i<wdth;i++) 
        {
         fill_mri_32b2(img, buff32b, i, j, k, l, coeff,_read_format_3d, _is_filebigEndian);
         l=l+pasw;
        }
        l=l+filewdth*(pash-1);
       }
      /* for (i=0;i<pasd-1;i++)
       {
        size=FREAD_LONG2(buff32b,filewdth*filehght,fp);
       }*/
#ifdef HAS_NIFTI
       znzseek(fp,(pasd-1)*filewdth*filehght*sizeof(long),SEEK_CUR);
#else
	fseek(fp,(pasd-1)*filewdth*filehght*sizeof(long),SEEK_CUR);
#endif
      }


      free(buff32b);
      break;

    }

  /*MISE EN PLACE DES BONNES VALEURS DANS WIDHT HEIGHT et DEPTH, dx, dy, dz */
  find_whd_from_orient_p(img, wdth, hght, dpth, _read_format_3d);


  /*  Fermeture fichier image  */
#ifdef HAS_NIFTI
  znzclose(fp);
#else
  fclose(fp);
#endif

  return(0);
}
#endif

/*********************************************
 **  --  pix_value_3d_p()  ---------------------------------
 **
 **
 ***********************************************/
float    pix_value_3d_p(grphic3d *im1, float *X, int TYPE)
{
  //    int i,j,k; 	unused variable
  int i1,j1,k1;
  float val;

  switch (TYPE)
    {
    case 1:   /*Pixel le plus proche   */
      i1=X[0]+.5;
      j1=X[1]+.5;
      k1=X[2]+.5;
      if( i1<0 || i1>=im1->width || j1<0 || j1>=im1->height || k1<0 || k1>=im1->depth )
	val=0;
      else 
	val=im1->mri[i1][j1][k1];
      return((float) val);
      break;

    case 2:   /* Interpolation bilineaire a faire  */
      break;
    }

  return((float) 0.);

}

#ifdef __GTK_GUI


/***********************************************
 **   -- interp_3d()  -----------------------
 **
 **   Interpolation lineaire 3d sur une dimension(depth)
 **
 **
 ***********************************************/
int     interp_3d(int im_deb, int im_res, int TYPINTERP)
{
  //   grphic3d *imdeb,*imres;	unused variable

  imx_interp_3d(im_deb,im_res,TYPINTERP);
  show_picture_3d(im_res);
  return(1);
}

/***********************************************
 **   -- imx_interp_3d()  -----------------------
 **
 **   Interpolation lineaire 3d sur une dimension(depth)
 **
 **
 ***********************************************/
int     imx_interp_3d(int im_deb, int im_res, int TYPINTERP)
{
  grphic3d *imdeb,*imres;
  int i,j,k,z;
  int cptdeb,cptres;
  int hght,wdth,dpth,dpth_res;
  float vary,varx,val;
  float dx,dy,dz;
  int tab[MAX_DEPTH_3D];
  float sum=0.0,sum1;
   
  imdeb=ptr_img_3d(im_deb);
  imres=ptr_img_3d(im_res);

  dx=imdeb->dx;
  dy=imdeb->dy;
  dz=imdeb->dz;

  hght=imdeb->height;
  wdth=imdeb->width;
  dpth=imdeb->depth;

    
  dpth_res=(int)(0.5+dz/dx*((double)dpth));
  dpth_res=MINI(MAX_DEPTH_3D,dpth_res);

  for (i=0;i<MAX_DEPTH_3D;i++) tab[i]=0.;

  switch(TYPINTERP)
    {

    case 1:  /* Interpolation lineaire   */
      for(i=0;i<wdth;i++)
	for(j=0;j<hght;j++)
	  {
	    for(k=0;k<dpth;k++)
	      tab[k]=imdeb->mri[i][j][k];
	    k=0;
	    cptdeb=0;
	    imres->mri[i][j][k]=tab[k];
	    k++;
	    while(k<MAX_DEPTH_3D && (cptdeb+1)<dpth)
	      {
		if( (k*dx) > ((cptdeb+1)*dz) ) cptdeb++;
		varx=k*dx-cptdeb*dz;
		vary=tab[cptdeb+1]-tab[cptdeb];
		imres->mri[i][j][k]=tab[cptdeb]+(varx/dz)*vary;
		k++;
	      }
	  }
      break;

    case 2:   /* Interpolation sinus cardinal  */
		   
      for(i=0;i<wdth;i++)
	for(j=0;j<hght;j++)
	  {
	    for(k=0;k<dpth;k++)
	      tab[k]=imdeb->mri[i][j][k];
	       
	    k=0;
	           
	    while(k*dx<=2*dz) 
	      k++;
	      
	    cptres=0;
	    cptdeb=2;
	   
	    while(cptres<MAX_DEPTH_3D && cptdeb<dpth-2)
	      {
		if(k*dx>(cptdeb+1)*dz) 
		  cptdeb++;
	                 
		val=0.0;
		sum=0.0;

		for(z=cptdeb-2;z<=cptdeb+2;z++)
		  {
		    sum1=sin_car( (float) (dx*k-dz*z));
		    val+=((double)tab[z])*fabs(sum1);
		    sum+=fabs(sum1);
		  }

		imres->mri[i][j][cptres]=(int)(0.5+val/sum);
                        
		k++;
   
		cptres++;
	      }
	  }

      break;

  
  
  
    }


  SetImg3d(im_res,imres);
  imres->width =imdeb->width;
  imres->height=imdeb->height;
  imres->depth=dpth_res;
  imres->min_pixel=imdeb->min_pixel;
  imres->max_pixel=imdeb->max_pixel;
  imres->cutoff_min=imdeb->cutoff_min;
  imres->cutoff_max=imdeb->cutoff_max;
  imres->icomp =imdeb->icomp;
  imres->rcoeff=imdeb->rcoeff;
  imres->pos=-1;
  imres->img_type =IMAGE_NORMAL;
  imres->type='b';
  imres->bitppixel=32;
  imres->x3dc=imdeb->x3dc;
  imres->y3dc=imdeb->y3dc;
  imres->z3dc=imdeb->z3dc;
  imres->dx=imdeb->dx;
  imres->dy=imdeb->dy;
  imres->dz=imdeb->dx;
  return(1);
}

/*****************************************************
 **
 ** --  sin_car() -------------------------------------
 **
 ** Fonction sinus cardinal pour les interpolations
 **
 ******************************************************/
float sin_car(float nb)
{
  if(nb==0) return( (float) 1);
  else return( (float) sin( (double) PI*nb)/(PI*nb) );
}

/***********************************************
 **   -- interp_tri_3d()  -----------------------
 **
 **   Interpolation trilineaire 3d
 **
 **
 ***********************************************/
int     interp_tri_3d(int im_deb, int im_res)
{
  //    grphic3d *imdeb,*imres;	unused variable

  imx_interp_tri_3d(im_deb,im_res);
  show_picture_3d(im_res);
  return(1);
}

/***********************************************
 **   -- imx_interp_tri_3d()  -----------------------
 **
 **   Interpolation trilineaire 3d 
 **
 **
 ***********************************************/
int     imx_interp_tri_3d(int im_deb, int im_res)
{
  grphic3d *imdeb,*imres,*temp;
  double dx,dy,dz;
  int wdth,hght,dpth;


  imdeb=ptr_img_3d(im_deb);
  imres=ptr_img_3d(im_res);

  dx=imdeb->dx;
  dy=imdeb->dy;
  dz=MAXI(dx,dy);

  wdth=imdeb->width;
  hght=imdeb->height;
  dpth=((double)imdeb->depth)*(imdeb->dz)/dz;


  temp=cr_grphic3d(imdeb);
  temp->width=imdeb->width; 
  temp->height=imdeb->height; 
  temp->depth=dpth; 
  temp->dx=imdeb->dx;
  temp->dy=imdeb->dy;
  temp->dz=dz;

  imx_zoom_image_3d_p(imdeb,temp,imres,TRILINEAR_ROTATION_3D);

  free_grphic3d(temp);

  return(1);
}



/******************************************************************************
 ** -- Open_3d() ---------------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July, 23rd 1999
 **
 ******************************************************************************/

void Open_3d(void)
{
  int e;
  int is_targz=0;

  int format;

  /* Load 3d files.             */
  /*   1st, Call file selection box and:    */
  /*   2nd, Make severals files.         */
  char *file;
  int pos,max=0,step=0,numfirst=1;
  file=GET_FILE_IMG(NULL, &numfirst, &e);
  if (e) return;
  if(file==NULL || strlen(file)<=0)
    PUT_WARN("Sorry, cannot read picture !\n");
  else {

	is_targz=untargz_Ipb_file(file);

    format=is_mri(file,NULL);
    printf("format : %d\n",format);
    switch(format)
      {
#ifdef HAS_NIFTI
      case NIFTI_GZ_FORMAT :
      case NIFTI_FORMAT :
	pos=GET_WND3D(TEXT0001);
	if(load_picture_3d(pos,file,max,step,numfirst)==-1)
	  PUT_ERROR("Cannot read picture.\n");
	break;
#endif

      case BRUKER_FORMAT :
	max= GET_INT(TEXT0004, 1, &e);
	if (e) break;
	step= GET_INT(TEXT0024, 1, &e);
	if (e) break;
	pos=GET_WND3D(TEXT0001);
	if(load_picture_3d(pos,file,max,step,numfirst)==-1)
	  PUT_ERROR("Cannot read picture.\n");
	break;

      case SIEMENS_FORMAT :
	max= GET_INT(TEXT0004, 0, &e);
	if (e) break;
	step= GET_INT(TEXT0024, 0, &e);
	if (e) break;
	pos=GET_WND3D(TEXT0001);

	if(load_picture_3d(pos,file,max,step,numfirst)==-1)
	  PUT_ERROR("Cannot read picture.\n");
	imx_inimaxminpixel_3d(pos);
	show_picture_3d(pos);
	break;


      case VANDERBILT_FORMAT :
	max= GET_INT(TEXT0004, 0, &e);
	if (e) break;
	step= GET_INT(TEXT0024, 0, &e);
	if (e) break;
	pos=GET_WND3D(TEXT0001);

	if(load_picture_3d(pos,file,max,step,numfirst)==-1)
	  PUT_ERROR("Cannot read picture.\n");
	imx_inimaxminpixel_3d(pos);
	show_picture_3d(pos);
	break;
                
      case INTERFILE_FORMAT :
	max= GET_INT(TEXT0004, 1, &e);
	if (e) break;
	step= GET_INT(TEXT0024, 1, &e);
	if (e) break;
	pos=GET_WND3D(TEXT0001);
	if(load_picture_3d(pos,file,max,step,numfirst)==-1)
	  PUT_ERROR("Cannot read picture.\n");
	break;

      case PARAVISION_FORMAT :
	pos=GET_WND3D(TEXT0001);
	if(load_picture_3d(pos,file,max,step,numfirst) == -1)
	  PUT_ERROR("Cannot read picture.\n");
	break;

      case IPB_FORMAT :
	pos=GET_WND3D(TEXT0001);
	if(load_picture_3d(pos,file,max,step,numfirst) == -1)
	  PUT_ERROR("Cannot read picture.\n");
	break;
      case IPB_2D_FORMAT :
	max= GET_INT(TEXT0004, 1, &e);
	if (e) break;
	step= GET_INT(TEXT0024, 1, &e);
	if (e) break;
	pos=GET_WND3D(TEXT0001);
	if(load_picture_3d(pos,file,max,step,numfirst) == -1)
	  PUT_ERROR("Cannot read picture.\n");
	break;
      case IPB_3D_FORMAT :
	pos=GET_WND3D(TEXT0001);

	if(load_picture_3d(pos,file,0,0,numfirst) == -1)
	  PUT_ERROR("Cannot read picture.\n");
	break;

      case ACRNEMA2L_FORMAT :
	max= GET_INT(TEXT0004, 0, &e);
	if (e) break;
	step= GET_INT(TEXT0024, 0, &e);
	if (e) break;
	pos=GET_WND3D(TEXT0001);
	if(load_picture_3d(pos,file,max,step,numfirst)==-1)
	  PUT_ERROR("Cannot read picture.\n");
	break;
				
      case ACRNEMA2B_FORMAT :
	max= GET_INT(TEXT0004, 0, &e);
	if (e) break;
	step= GET_INT(TEXT0024, 0, &e);
	if (e) break;
	pos=GET_WND3D(TEXT0001);
	if(load_picture_3d(pos,file,max,step,numfirst)==-1)
	  PUT_ERROR("Cannot read picture.\n");
	break;
				
      case ACRNEMA1L_FORMAT :
	max= GET_INT(TEXT0004, 0, &e);
	if (e) break;
	step= GET_INT(TEXT0024, 0, &e);
	if (e) break;
	pos=GET_WND3D(TEXT0001);
	if(load_picture_3d(pos,file,max,step,numfirst)==-1)
	  PUT_ERROR("Cannot read picture.\n");
	break;
				
      case ACRNEMA1B_FORMAT :
	max= GET_INT(TEXT0004, 0, &e);
	if (e) break;
	step= GET_INT(TEXT0024, 0, &e);
	if (e) break;
	pos=GET_WND3D(TEXT0001);
	if(load_picture_3d(pos,file,max,step,numfirst)==-1)
	  PUT_ERROR("Cannot read picture.\n");
	break;
				
      case DICMBIMP_FORMAT : 
      case DICMLIMP_FORMAT : 
      case DICMBEXP_FORMAT : 
      case DICMLEXP_FORMAT : 
	max= GET_INT(TEXT0004, 0, &e);
	if (e) break;
	step= GET_INT(TEXT0024, 0, &e);
	if (e) break;
	pos=GET_WND3D(TEXT0001);
	if(load_picture_3d(pos,file,max,step,numfirst)==-1)
	  PUT_ERROR("Cannot read picture.\n");
	break;

      case SMIS_FORMAT :
	//                    max= GET_INT(TEXT0004, 0, &e);
	//                    if (e) break;
	//                    step= GET_INT(TEXT0024, 0, &e);
	//                    if (e) break;
	pos=GET_WND3D(TEXT0001);
	if(load_picture_3d(pos,file,0,1,numfirst)==-1)
	  PUT_ERROR("Cannot read picture.\n");
	break;

      case ECAT7_FORMAT :			
	pos=GET_WND3D(TEXT0001);
	if(load_picture_3d(pos,file,max,step,numfirst)==-1)
	  PUT_ERROR("Cannot read picture.\n");
	break;

      case ANALYZE_FORMAT_LITTLE_ENDIAN :
      case ANALYZE_FORMAT_BIG_ENDIAN :				
	pos=GET_WND3D(TEXT0001);
	if(load_picture_3d(pos,file,max,step,numfirst) == -1)
	  PUT_ERROR("Cannot read picture.\n");
	break;
      case XML_FILEFORMAT	:
         pos=GET_WND3D(TEXT0001);

	if(load_picture_3d(pos,file,0,0,numfirst) == -1)
	  PUT_ERROR("Cannot read picture.\n");
	break;
      default : 
	pos=GET_WND3D(TEXT0001);
	if(load_picture_3d(pos,file,max,step,numfirst) == -1)
	  PUT_ERROR("Cannot read picture.\n");
	break;
      }
  
  if (is_targz)
  	{
	delete_ipb_file(file);
	}
	
  } 
}

/******************************************************************************
 ** -- nOpen_3d() ---------------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July, 23rd 1999
 **
 ******************************************************************************/

void nOpen_3d(void)
{
  int e;
int is_targz=0;
	
  char *file;
  // unused variable			char **quest,answer[80];
  int pos,i,max,step,numfirst;
  // unused variable			int k;
			
  file=GET_FILE_IMG(NULL, &numfirst, &e);
  if (e) return;
  if(file==NULL || strlen(file)<=0)
    PUT_WARN("Sorry, cannot read picture !\n");
  else {
    
    is_targz=untargz_Ipb_file(file);

    switch(is_mri(file,NULL))
      {
      case IPB_3D_FORMAT :
	max= GET_INT(TEXT0004, 4, &e);
	if (e) break;
	step= GET_INT(TEXT0024, 1, &e);
	if (e) break;
	pos=GET_WND3D(TEXT0001);
                    
	for(i=0; i<max; i++) {
	  if(load_picture_3d(pos,file,0,0,numfirst) == -1) {
	    PUT_ERROR("Cannot read picture.\n");
	    break;
	  }
	  pos ++;
	  numfirst ++;
	}
	break;

      case PARAVISION_FORMAT :
	max= GET_INT(TEXT0004, 4, &e);
	if (e) break;
	step= GET_INT(TEXT0024, 1, &e);
	if (e) break;
	pos=GET_WND3D(TEXT0001);
                    
	for(i=0; i<max; i++) {
	  if(load_picture_3d(pos,file,0,0,numfirst) == -1) {
	    PUT_ERROR("Cannot read picture.\n");
	    break;
	  }
	  pos ++;
	  numfirst ++;
	}
	break;

      default : 
	break;                 
      }
   if (is_targz)
  	{
	delete_ipb_file(file);
	}
  }
}

/******************************************************************************
 ** -- Save_3d_xml() -----------------------------------------------------------
 ******************************************************************************/
#ifdef HAS_CORRECT_BOOST
void Save_3d_xml(void)
{

  char *questAction[3], **questImgName;
  int choixAction, i, err;
  int tAtlas;
  char *fileName, imgName[32];
  char datafile[PATH_LEN];

  grphic3d* im1 = ptr_img_3d(GET_WND3D("image a ajouter à l'atlas"));
  
  for (i=0; i<3; i++)
    questAction[i] = CALLOC(25,char);

  strcpy(questAction[0],"atlas existant");
  strcpy(questAction[1],"nouvel atlas");
  strcpy(questAction[2],"\0");

  choixAction = GETV_QCM("ecrire dans :",(const char**)questAction);
  for (i=0; i<3; i++) FREE(questAction[i]);

  if(choixAction == 0)
  {
    fileName = GET_FILE("fichier atlas", &err);

    remove_file_extension(fileName,datafile);
    put_file_extension(datafile,".xml",datafile);

    if(err!=0)  return;
    tAtlas = get_type_atlas(datafile);
    if(tAtlas < 0 || tAtlas >= XML_TYPE_ATLAS_SIZE) return;

    questImgName = (char **)XML_NOM_IMAGE[tAtlas];
    char titleImgName[64];
    strcpy(titleImgName,"nom de l'image\ndans atlas ");
    strcat(titleImgName,XML_NOM_ATLAS[tAtlas]);
    strcpy(imgName,questImgName[GETV_QCM(titleImgName,(const char**) questImgName)]);

    err = saveInExistingAtlas(im1, fileName, imgName);
  }
  else
  {
    tAtlas = GETV_QCM("Type de l'atlas",(const char**) XML_NOM_ATLAS);

    questImgName = (char **)XML_NOM_IMAGE[tAtlas];
    strcpy(imgName,questImgName[GETV_QCM("nom de l'image",(const char**) questImgName)]);

    fileName = GET_FILE("nouveau fichier atlas", &err);
    
    if(err!=0) return;

    remove_file_extension(fileName,datafile);
    put_file_extension(datafile,".xml",datafile);
    
    FILE *fp = fopen (datafile, "r");

    if (fp == NULL) // ecrase-t-on un fichier ? si non on continue
    {
      err = saveInNewAtlas(im1, fileName, tAtlas, imgName);
    }
    else // si oui on n'écrit pas
    {
      PUT_ERROR("Ce fichier existe deja");
      fclose(fp);
    }
  }
  if(err!=0) PUT_ERROR("Erreur de sauvegarde de l'image");

}
#endif
/******************************************************************************
 ** -- Save_3d_Ipb() -----------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July, 23rd 1999
 **
 ******************************************************************************/

void Save_3d_Ipb(void)
{
  char	file[PATH_LEN];
  int 	pos, e;
	
  strcpy(file,GET_FILE(TEXT0420, &e));
  put_file_extension(file,".ipb",file);
	
  if(file == NULL || strlen(file) <= 0 || e) {
    PUT_WARN( TEXT0421);
  }else 
    {
      pos = GET_WND3D( NULL);
      /*if (e) 
	return;*/
      if (save_pictures_3d(file, pos, 1, 1) == 0)
	PUT_ERROR(TEXT0422);
      if (e) 
	return;
    }
}


/******************************************************************************
 ** -- Save_3d_Ipb_targz() -----------------------------------------------------------
 **
 ** >> Last user action by: V. Noblet 07/11/07
 **
 ******************************************************************************/

void Save_3d_Ipb_targz(void)
{
  char	file[PATH_LEN];
  int 	pos, e;
	
  strcpy(file,GET_FILE(TEXT0420, &e));
  put_file_extension(file,".ipb",file);
	
  if(file == NULL || strlen(file) <= 0 || e) {
    PUT_WARN( TEXT0421);
  }else 
    {
      pos = GET_WND3D( NULL);
      /*if (e) 
	return;*/
      if (save_pictures_3d(file, pos, 1, 1) == 0)
	PUT_ERROR(TEXT0422);
	else
	{
	 targz_Ipb_file(file);
	}
      if (e) 
	return;
    }
}

/******************************************************************************
 ** -- targz_Ipb_file() -----------------------------------------------------------
 **
 ** >> Last user action by: V. Noblet 07/11/07
 **
 ******************************************************************************/
void targz_Ipb_file(char *file)
{
char commande_unix[PATH_LEN];
char nom_archive[PATH_LEN];
char nom_fichier[PATH_LEN];
char nom_fichier_a_effacer[PATH_LEN];
char* nom_fichier_basename;
char* nom_fichier_dirname;
char *ext;
	

nom_fichier_basename = basename(file);
nom_fichier_dirname = dirname(file);


strcpy(nom_fichier,nom_fichier_basename);

ext=strstr(nom_fichier,".ipb");
if (ext!=NULL) *ext='\0';

strcpy(nom_archive,nom_fichier);

ext=strstr(nom_archive,".tar.gz");
if (ext!=NULL) 
	{*ext='\0';
	
	// On renomme les fichiers images car elles sont enregistrees avec l'extension tar.gz
	strcpy(commande_unix,"cd ");
	strcat(commande_unix,nom_fichier_dirname);
	strcat(commande_unix,"; mv ");
	strcat(commande_unix,nom_fichier);
	strcat(commande_unix,".ipb ");
	strcat(commande_unix,nom_archive);
	strcat(commande_unix,".ipb ");
	strcat(commande_unix,"; mv ");
	strcat(commande_unix,nom_fichier);
	strcat(commande_unix,".img ");
	strcat(commande_unix,nom_archive);
	strcat(commande_unix,".img ");
	//printf("%s \n",commande_unix);

	system(commande_unix);
	}


strcpy(nom_fichier,nom_archive);
strcat(nom_archive,".tar.gz ");
	
strcpy(commande_unix,"cd ");
strcat(commande_unix,nom_fichier_dirname);
strcat(commande_unix,"; tar czf ");
strcat(commande_unix,nom_archive);
strcat(nom_fichier,".i*");
strcat(commande_unix,nom_fichier);

/* creation de l'archive */
//printf("%s \n",commande_unix);

system(commande_unix);
	
strcpy(commande_unix,"rm ");
strcat(commande_unix,nom_fichier_dirname);
strcat(commande_unix,"/");
strcat(commande_unix,nom_fichier);

//printf("%s \n",commande_unix);
	
/* effacement des fichiers ipb et img */
system(commande_unix);
}



/******************************************************************************
 ** -- untargz_Ipb_file() -----------------------------------------------------------
 **
 ** >> Last user action by: V. Noblet 07/11/07
 **
 ******************************************************************************/
int untargz_Ipb_file(char *file)
{
char *ext;
int res=0;
	
ext=strstr(file,".tar.gz");
if (ext!=NULL) 
	{
	char commande_unix[PATH_LEN];
	char* file_basename;
	*ext='\0';

	strcpy(commande_unix,"tar -xzf ");
	strcat(commande_unix,file);
	strcat(commande_unix,".tar.gz ");
	strcat(commande_unix,"-C /tmp/");

	/* dezippage de l'archive */	
	system(commande_unix);
	
	file_basename = basename (file);
	strcpy(file,"/tmp/");
	strcat(file,file_basename);
	strcat(file,".ipb");
		
	res=1;
	}
return(res);
}


/******************************************************************************
 ** -- delete_ipb_file() -----------------------------------------------------------
 **
 ** >> Last user action by: V. Noblet 07/11/07
 **
 ******************************************************************************/
void delete_ipb_file(char* file)
{
char *ext;
	
ext=strstr(file,".ipb");
if (ext!=NULL) 
	{
	char datafile[PATH_LEN];
	char header[PATH_LEN];
	*ext='\0';

	strcpy(datafile,file);
	strcat(datafile,".img");
	
	strcpy(header,file);
	strcat(header,".ipb");
	
	remove(datafile);
	remove(header);
	}
}

/******************************************************************************
 ** -- is_targz(file) -----------------------------------------------------------
 **
 ** >> Last user action by: V. Noblet 07/11/07
 **
 ******************************************************************************/
int is_targz(char* file)
{
char *ext;

if (file == NULL)
 return(0);
 
ext=strstr(file,".tar.gz");
if (ext!=NULL) 
	return(1);
else 
	return(0);

}



/******************************************************************************
 ** -- Save_3d_Ipb_std() -----------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July, 23rd 1999
 ** >> Last user action by: F. BOULEAU on July, 23rd 1999
 **  Fct creer pour avoir toujours les memes questions que le fichier exist ou pas
 **  Fct pas visible du menu
 ******************************************************************************/
void Save_3d_Ipb_std(void)
{
  char	file[PATH_LEN];
  int 	pos, e;
	
  strcpy(file,GET_FILE(TEXT0420, &e));
  put_file_extension(file,".ipb",file);
	
  if(file == NULL || strlen(file) <= 0 || e) {
    PUT_WARN( TEXT0421);
  }else 
    {
      pos = GET_WND3D( NULL);
      if (save_pictures_3d_std(file, pos, 1, 1) == 0)
	PUT_ERROR(TEXT0422);
      if (e) 
	return;
    }
}

/******************************************************************************
 ** -- Save_3d_Ipb_erase() -----------------------------------------------------------
 **
 ** >> Last user action by: F. ROUSSEAU on Sept, 4th 2006
 ** >> Last user action by: F. ROUSSEAU on Sept, 4th 2006
 **  Fct creer qui permet de sauver un fichier et d'ecraser automatiquement le fichier existant
 **  Fct pas visible du menu
 ******************************************************************************/
void Save_3d_Ipb_erase(void)
{
  char	file[PATH_LEN]; 
  int 	pos, e;
	
  strcpy(file,GET_FILE(TEXT0420, &e));
  put_file_extension(file,".ipb",file);
	
  if(file == NULL || strlen(file) <= 0 || e) {
    PUT_WARN( TEXT0421);
  }else 
    {
      pos = GET_WND3D( NULL);
      if (save_pictures_3d_erase(file, pos, 1, 1) == 0)
	PUT_ERROR(TEXT0422);
      if (e) 
	return;
    }
}


/******************************************************************************
 ** -- Save_3d_Ipb_mask_erase() -----------------------------------------------------------
 **
 ** >> Last user action by: JPA on Sep 8th 2008
 **
 ******************************************************************************/

void Save_3d_Ipb_mask_erase(void)
{
  char	file[PATH_LEN];
  int 	pos, e;


  //Ecrase le fichier si il existe : 
  char file_without_ext[PATH_LEN];
  char headerfile[PATH_LEN];
  char datafile[PATH_LEN];
  char s[64];
  struct stat statbuff;


  strcpy(file,GET_FILE(TEXT0420, &e));
  put_file_extension(file,".ipb",file);
  pos = GET_WND3D( NULL);

  if (!ptr_has_mask_3d(pos))
    {
      PUT_WARN("ERROR : no mask in  Save_3d_Ipb_mask!!\n");
      return;   
    }
	
  /* Ecrase s il existe */
  remove_file_extension(file,file_without_ext);
  strcpy(headerfile, file_without_ext);
  strcpy(datafile, file_without_ext);
  strcat(headerfile, ".ipb");
  strcat(datafile, ".img");

  if(stat(headerfile, &statbuff) == 0) 
    { 
      /* On supprime les fichiers .ipb, et .img */

      if(remove(headerfile) == -1) 
	{
	  sprintf(s, "[save_header] unable to remove file \"%s\"\n",
		  headerfile);	 
	  PUT_ERROR( s);
	  return; /* Erreur */
	}  
         
      if(remove(datafile) == -1) 
	{
	  sprintf(s,"[save_header] unable to remove file \"%s\"\n",
		  datafile);
	  PUT_ERROR( s);
	  return; /* Erreur */
	}	 		 
    }
  //Fin de : Ecrase le fichier si il existe  

  // save the *;ipb file
  if(file == NULL || strlen(file) <= 0 || e) 
    {
      PUT_WARN( TEXT0421);
      return;
    }
  else 
    {
      Save_3d_Ipb_mask_p(file, ptr_img_3d(pos));    
    }
	
  //update the header file
  if (!add_mask_Ipb_3d(file,pos))
    {
      PUT_ERROR(TEXT0422);
      return ;
    }	

}


/******************************************************************************
 ** -- check_ipb_format_3d() -----------------------------------------------------------
 **
 ** >> Last user action by: T. BERST on April, 15th 2003
 **
 ******************************************************************************/

int  check_ipb_format_3d(char * file)
{
  struct stat statbuff;
  char file_without_ext[FILE_LEN];
  char headerfile[FILE_LEN];
  char answer[64];
	
  remove_file_extension(file,file_without_ext);
  strcpy(headerfile, file_without_ext);
  strcat(headerfile, ".ipb");
  if(stat(headerfile, &statbuff) == 0) 
    {
      strcpy(answer,getheader_interfile(headerfile,"MASK_OFFSET=",ENTIER_IPB_2D,1));
      if (strlen(answer)==0)
	{
	  return IPB_WITHOUT_MASK_FORMAT;
	}
      else 
	{
	  return IPB_WITH_MASK_FORMAT;
	}
    }
  return 0;

}

/******************************************************************************
 ** -- Save_3d_Ipb_mask_p() -----------------------------------------------------------
 **
 ** >> Last user action by: T. BERST on July, 8th 2002
 **
 ******************************************************************************/

int  Save_3d_Ipb_mask_p(char * file, grphic3d *img)
{
  int retVal=0;
  int format=0;
  //faire un check si le file existe pour verifier que c'est un format ipb+mask
  format = check_ipb_format_3d(file);
  if (format && format != IPB_WITH_MASK_FORMAT)
    {
      PUT_ERROR("Erreur : ne peut pas ajouter une image+masque a une image sans masque");
      return 0;
    }
  if (check_file_3d(file)!=0 ) return (0);

  save_mri_ipb_3d_p(file,img);
  if (img->mask)
    save_mri_ipb_3d_p(file,img->mask);
  else
    retVal = 1;
  return retVal;
}

/******************************************************************************
 ** -- Save_3d_Ipb_mask() -----------------------------------------------------------
 **
 ** >> Last user action by: T. BERST on July, 8th 2002
 **
 ******************************************************************************/

void Save_3d_Ipb_mask(void)
{
  char	file[PATH_LEN];
  int 	pos, e;

  strcpy(file,GET_FILE(TEXT0420, &e));
  put_file_extension(file,".ipb",file);
	
  check_file_3d(file);
  pos = GET_WND3D( NULL);
  if (!ptr_has_mask_3d(pos))
    {
      PUT_WARN("ERROR : no mask in  Save_3d_Ipb_mask!!\n");
      return;   
    }
	
  // save the *;ipb file
  if(file == NULL || strlen(file) <= 0 || e) 
    {
      PUT_WARN( TEXT0421);
      return;
    }
  else 
    {
      Save_3d_Ipb_mask_p(file, ptr_img_3d(pos));    
    }
	
  //update the header file
  if (!add_mask_Ipb_3d(file,pos))
    {
      PUT_ERROR(TEXT0422);
      return ;
    }	

}

/******************************************************************************
 ** -- Save_3d_Ipb_mask_targz() -----------------------------------------------------------
 **
 ** >> Last user action by: V.Noblet 07/11/07
 **
 ******************************************************************************/

void Save_3d_Ipb_mask_targz(void)
{
  char	file[PATH_LEN];
  int 	pos, e;

  strcpy(file,GET_FILE(TEXT0420, &e));
  put_file_extension(file,".ipb",file);
	
  check_file_3d(file);
  pos = GET_WND3D( NULL);
  if (!ptr_has_mask_3d(pos))
    {
      PUT_WARN("ERROR : no mask in  Save_3d_Ipb_mask!!\n");
      return;   
    }
	
  // save the *;ipb file
  if(file == NULL || strlen(file) <= 0 || e) 
    {
      PUT_WARN( TEXT0421);
      return;
    }
  else 
    {
      Save_3d_Ipb_mask_p(file, ptr_img_3d(pos));    
    }
	
  //update the header file
  if (!add_mask_Ipb_3d(file,pos))
    {
      PUT_ERROR(TEXT0422);
      return ;
    }	
    
 targz_Ipb_file(file);

}

/* --save_mri_img_3d() ---------------------------------------------------------
**
**    save_mri_img_3d() saves a picture in img 3D file format.
**
**    Usage ...
**   file :       nom du fichier  !attention : aucun ajout d'extension au nom
**   start :      numero de la premier image a sauver 
**   num_images : nombre d'images a sauver 
**   step :       pas si plusieurs images
**
**   return : 1 if Ok  !=1 if error
*/
int save_mri_img_3d(const char *file, int start, int num_images, int step)
{
  grphic3d *image;
  int end;
  char s[256];
    
  end = (num_images - 1) * step + start;
  if(end > MAX_PICTURES_3D) {
    end = MAX_PICTURES_3D;
    PUT_WARN("[save_mri_ipb_3d] warning: cannot save all images !\n");
  }
  
  while(start <= end) {
    if((image = ptr_img_3d(start)) == NULL) {
      sprintf(s,"[save_mri_img_3d] unable to get image %d info!\n", start);
      PUT_ERROR(s);
      return(0);
    }
  
    if ( save_mri_img_3d_p(file, image) == 0) 
      return(0) ;
    start += step;
  }

  return(1);
        
}

/* --save_mri_img_3d_p() ---------------------------------------------------------
**
**    save_mri_img_3d_p() saves a picture in IPB 3D file format.
**
**    Usage ...
**   file :       nom du fichier
**   image : pointeur sur struicture grphic3d 
**
**   return : 0 if Ok  !=0 if error
*/

int save_mri_img_3d_p(const char *file_img, grphic3d *image)
{
  char s[256];
  int data_size;
  int index;
  long offset;
  char *buffchar;
  short *buffshort;
  long *bufflong;
  void *buffer;
  int i, j, k, num_items;
  FILE *fp;

  data_size = sizeof(TYPEMRI3D);

  index = 0;
  offset = 0;

  index ++;

  if((fp = fopen(file_img, "ab")) == NULL) {
    sprintf(s,"Unable to save file %s !\n", file_img);
    PUT_ERROR(s);
    return(0);
  }


  offset += image->width * image->height * image->depth * data_size;

  num_items = image->width * image->height;

  switch(data_size)
    {
    case sizeof(char):
      buffchar = (char *)malloc(num_items);
      if(buffchar == NULL) {
	PUT_ERROR("[save_mri_ipb_3d] memory allocation error (2)!\n");
	fclose(fp);
	return(0);
      }
      buffer = buffchar;

      for(k=0; k<image->depth; k++)
	{
	  buffchar = buffer;

	  for(j=0; j<image->height; j++)
	    for(i=0; i<image->width; i++)
	      *buffchar ++ = image->mri[i][j][k];

	  if(fwrite(buffer, data_size, num_items, fp) != num_items) {
	    PUT_ERROR("[save_mri_ipb_3d] error writing data (2)!\n");
	    fclose(fp);
	    return(0);
	  }
	}
      break;

    case sizeof(short):
      buffshort = (short *)malloc(num_items * sizeof(short));
      if(buffshort == NULL) {
	PUT_ERROR("[save_mri_ipb_3d] memory allocation error (4)!\n");
	fclose(fp);
	return(0);
      }
      buffer = buffshort;

      for(k=0; k<image->depth; k++)
	{
	  buffshort = buffer;

	  for(j=0; j<image->height; j++)
	    for(i=0; i<image->width; i++)
	      if (_is_bigEndian) *buffshort ++ = image->mri[i][j][k];
	      else	*buffshort ++ = shEndianConversion(image->mri[i][j][k]);

	  if(fwrite(buffer, data_size, num_items, fp) != num_items) {
	    PUT_ERROR("[save_mri_ipb_3d] error writing data (4)!\n");
	    fclose(fp);
	    return(0);
	  }
	}
      break;

    case sizeof(long):
      bufflong = (long *)malloc(num_items * sizeof(long));
      if(bufflong == NULL) {
	PUT_ERROR("[save_mri_ipb_3d] memory allocation error (6)!\n");
	fclose(fp);
	return(0);
      }

      buffer = bufflong;

      for(k=0; k<image->depth; k++)
	{
	  bufflong = buffer;

	  for(j=0; j<image->height; j++)
	    for(i=0; i<image->width; i++)
	      if (_is_bigEndian) *bufflong ++ = image->mri[i][j][k];
	      else	*bufflong ++ = longEndianConversion(image->mri[i][j][k]);

	  if(fwrite(buffer, data_size, num_items, fp) != num_items) {
	    PUT_ERROR("[save_mri_ipb_3d] error writing data (6)!\n");
	    fclose(fp);
	    return(0);
	  }
	}
      break;

    default:
      PUT_ERROR("[save_mri_ipb_3d] invalid data size\n");
      fclose(fp);
      return(0);
    }

  free(buffer);

  fclose(fp);
  return(1);
}


/* --add_mask_Ipb_3d_3d() ---------------------------------------------------------
**
**    add_mask_Ipb_3d() adds mask item in the *.ipb.
**
**    Usage ...
**   file :       nom du fichier
**   pos :        position de l'image dont on veut sauver le masque 
**   
**	return 0 si echec 1 sinon 
*/

int add_mask_Ipb_3d(char * file, int pos)
{
  int index;
  HEADER header;
  int entier;
  char s[PATH_LEN];
  char file_with_extension[PATH_LEN];
	
  put_file_extension(file,".ipb",file_with_extension);
	
  if(load_header(file_with_extension, &header) == 0) {
    PUT_ERROR( "[save_mri_ipb_3d] Error reading header (3)!\n");
    return(0);
  }
	
  if((index = getmri_numberimg(file_with_extension)) == 0) 
    {
      fprintf(stderr, "[save_mri_ipb_3d] Error reading header (1)!\n");
      return(0);
    }

  entier = index-1;
  put_header(&header, "ORIGINAL_OFFSET=", ENTIER, &entier, index/2);
  entier = index;
  put_header(&header, "MASK_OFFSET=", ENTIER, &entier, index/2);

  if(save_header(&header, file_with_extension) == 0) 
    {
      sprintf(s,"Unable to save file header %s !\n", file);
      PUT_ERROR(s);
      return(0);
    }

  /*
    char mask_file[PATH_LEN];
    char file_with_extension[PATH_LEN];
    int i;
	
    put_file_extension(file,".ipb",file_with_extension);
    for (i=strlen(file);i>0 && file[i]!='/';i--)
    ;
    i++;  // on est alle un cran trop loin,file[i]='premiere lettre du nom du fichier'
	
    remove_file_extension(&file[i],mask_file);
    put_file_extension(mask_file,"_mask.ipb",mask_file);


    put_header(&header, "MASK_FILE=", STRING, mask_file, 0);
    put_header(&header, "MASK_TYPE=", ENTIER, &(ptr_img_3d(pos)->mask_type), 0);

    if(save_header(&header, file_with_extension) == 0) {
    sprintf(s,"Unable to save file header %s !\n", file);
    PUT_ERROR(s);
    return(0);
    }
  */
  return(1);
}

/******************************************************************************
*******************************************************************************
** -- Save_Ipbr_3d() ----------------------------------------------------------
**
** >> Last user action by: F. BOULEAU on Aug., 10th 1999
**
******************************************************************************/

void Save_Ipbr_3d(void)
{
  char     *file;
  int       pos, e;
  grphic3d *tmp;
	
  file = GET_FILE("3D File ?", &e);
	
  if (file == NULL || strlen(file) <= 0 || e)
    PUT_WARN( "Sorry, cannot save!\n");
  else 
    {
      pos = GET_WND3D( NULL);

      tmp = cr_grphic3d(ptr_img_3d(pos));
      imx_copie_3d_p(ptr_img_3d(pos), tmp);
      scaleimg_3d_p(tmp, ptr_img_3d(pos));
		
      if (e) 
	return;
      if (save_pictures_3d(file, pos, 1, 1) == 0)
	PUT_ERROR(
		  "erreur enregistrement format IPB 3D\n");

      imx_copie_3d_p(tmp, ptr_img_3d(pos));
      free_grphic3d(tmp);
    }
}

/******************************************************************************
 ** -- nSave_3d_Ipb() ----------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July, 23rd 1999
 ** >> Last user action by: R. CHEVRIER on Dec, 03rd 2001
 **     Bug fix. var file_to_save was truncated.
 ******************************************************************************/

void nSave_3d_Ipb(void)
{
  /*	const */char *file_to_save;
  int pos, max, step, e;
  file_to_save=strdup(GET_FILE("3D File", &e));
  if(file_to_save == NULL || strlen(file_to_save) <= 0 || e)
    PUT_WARN("Sorry, cannot save!\n");
  else {
    pos=GET_WND3D(NULL);
    if (e) return;
    max= GET_INT(TEXT0004, 0, &e);
    if (e) return;
    step= GET_INT(TEXT0024, 0, &e);
    if (e) return;
    if (save_pictures_3d(file_to_save,pos,max,step) == 0)
      PUT_ERROR(
		"erreur enregistrement format IPB 3D\n");
  }
}


/******************************************************************************
 ** -- nSave_3d_Ipb_targz() ----------------------------------------------------------
 **
 ** >> Last user action by: V. Noblet 07/11/07
 ******************************************************************************/

void nSave_3d_Ipb_targz(void)
{
  /*	const */char *file_to_save;
  int pos, max, step, e;
  file_to_save=strdup(GET_FILE("3D File", &e));
  if(file_to_save == NULL || strlen(file_to_save) <= 0 || e)
    PUT_WARN("Sorry, cannot save!\n");
  else {
    pos=GET_WND3D(NULL);
    if (e) return;
    max= GET_INT(TEXT0004, 0, &e);
    if (e) return;
    step= GET_INT(TEXT0024, 0, &e);
    if (e) return;
    if (save_pictures_3d(file_to_save,pos,max,step) == 0)
      PUT_ERROR(
		"erreur enregistrement format IPB 3D\n");
  }
 targz_Ipb_file(file_to_save);
  
}


/******************************************************************************
 ** -- Int_Lin_3d() ------------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July, 23rd 1999
 **
 ******************************************************************************/

void Int_Lin_3d(void)
{ 
  int im_deb,im_res;

  im_deb=GET_PLACE3D(TEXT0001);
  im_res=GET_PLACE3D(TEXT0006);
  interp_3d(im_deb,im_res,1);
}

/******************************************************************************
 ** -- Int_Tri_3d() ------------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July, 23rd 1999
 **
 ******************************************************************************/

void Int_Tri_3d(void)
{ 
  int im_deb,im_res;

  im_deb=GET_PLACE3D(TEXT0001);
  im_res=GET_PLACE3D(TEXT0006);
  interp_tri_3d(im_deb,im_res);
}

/******************************************************************************
 ** -- Int_Sinc_3d() -----------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July, 23rd 1999
 **
 ******************************************************************************/

void Int_Sinc_3d(void)
{ 
  int im_deb,im_res;

  im_deb=GET_PLACE3D(TEXT0001);
  im_res=GET_PLACE3D(TEXT0006);
  interp_3d(im_deb,im_res,2);
}

/*
**  -- Open_raw_3d() ----------------------------------------------------------
**
**  >> Last user action by: JC. PELOUTIER on May, 17th 2000
**
*/

void Open_raw_3d(void)
{

  int pos,taille,dimx,dimy,dimz,e;

  char *file;
  const char *bppxl[] = { "8  bits", "16 bits", "       ",
		    "32 bits", "Cancel", "" };
  
  const char *echanges[] = { "Pas de modification", "echange k et j", "echange k et i",
		       "echange i et j", "echange k et i, miroir k", "echange k et j, miroir j", 
		       "echange k et j, miroir i", "echange k et j, miroir k", 
		       "i->k ; j->i ; k->j", "i->j ; j->k ; k->i", "i->i ; j->j ; k->k", ""};
  const char *endianess[] = { "LittleEndian", "BigEndian", ""};

  file = GET_FILE(NULL, &e);  
  if (e) return; 
  
  /* choix taille du pixel en byte */
  taille = GETV_QCM( TEXT0043,(const char **)bppxl);
  
  if (taille==4) return;
  
  /* Choix du format (read_format) de lecture */
  _read_format_3d=GETV_QCM( "read format",(const char **)echanges);

  _is_filebigEndian=GETV_QCM(
			     "Endianess?",(const char **)endianess);

  dimx = GET_INT(TEXT0125, 0, &e);

  if (e) return;

  dimy = GET_INT(TEXT0126, 0, &e);

  if (e) return;

  dimz = GET_INT("Dimension Z", 0, &e);

  if (e) return;

  pos = GET_WND3D(NULL);

  load_picture_raw_3d(file,1,dimx,dimy,dimz,pos,taille+1);
}



/******************************************************************************
 ** -- nOpen_raw_3d() -----------------------------------------------------------
 **
 ** >> Last user action by: JC. PELOUTIER on May, 17th 2000
 **
 ******************************************************************************/

void nOpen_raw_3d(void)
{ 
  UnAvailable();
}

#endif /* __GTK_GUI */

/******************************************************************************
 ** -- scaleimg_3d_p() ---------------------------------------------------------
 **
 ** >> Derniere revision par F. BOULEAU le  4 Aout 1999
 **
 **  Redimensionne une image 3d selon une methode parmi plusieurs.
 **
 **  Entree : imres : image resultat
 **
 **  Sortie : imres est modifie
 **
 ******************************************************************************/

int scaleimg_3d_p(grphic3d *imdeb, grphic3d *imres)
{
  int meth, /* methode de redimensionnement      */
	factx, facty, factz,		/* facteurs de reduction             */
	rep; /* reponse au qcm */     
    
  int err; /* erreur get_int */
  const char	*qcm1[] = {"1:2", "1:4", "Personnalise", ""};
  const char  *qcm2[] = {"Saut de pixels", "Moyenne des pixels",
		   "Mediane des pixels", "MIP", "" };

  rep = GETV_QCM( "Reduction ?", (const char **) qcm1);


  switch(rep)
    {
    case 0:
      factx = facty = factz = 2;
      break;

    case 1:		/* 1:2 ou 1:4 */
      factx = facty = factz = 4;
      break;

    case 2:		/* Personnalise */
      factx =  GET_INT("X", 0, &err);

      if(!err)
	{
	  PUT_ERROR("Typing error!\n");
	  return ERR_TYPING;
	} 

      facty =  GET_INT("Y", 0, &err);

      if(!err)
	{
	  PUT_ERROR("Typing error!\n");
	  return ERR_TYPING;
	} 

      factz =  GET_INT("Z", 0, &err);

      if(!err)
	{
	  PUT_ERROR("Typing error!\n");
	  return ERR_TYPING;
	} 
      break;

    default:
      PUT_ERROR("Typing error!\n");
      return ERR_TYPING;
    }

 

  meth = GETV_QCM("Type de reduction ?", (const char **) qcm2);


  return imx_scaleimg_3d_p(imdeb, imres, meth, factx, facty, factz);
}

/******************************************************************************
 ** -- imx_scaleimg_3d_p() -----------------------------------------------------
 **
 ** >> Derniere revision par F. BOULEAU le  4 Aout 1999
 **
 **  Redimensionne une image 3d selon une methode parmi plusieurs.
 **
 **  Entree : imdeb : image de reference
 **           imres : image resultat
 **           meth  : methode a utiliser
 **           factx : facteur de reduction en x
 **           facty : facteur de reduction en y
 **           factz : facteur de reduction en z
 **
 **  Sortie : imres est modifie
 **
 ******************************************************************************/

int imx_scaleimg_3d_p(grphic3d *imdeb, grphic3d *imres, int meth, int factx, int facty, int factz)
     /* image de reference                */
     /* image resultat                    */
     /* methode de redimensionnement      */
     /* facteurs de reduction             */
{
  int  i, j, k;						/* indices de parcours de l'image    */
  int  resi, resj, resk;				/* indices des pixels de l'image res */
  int  mini=0, maxi=0, once = 0;
  int  res=0;							/* valeur du pixel                  */
  int  v;								/* calculs intermediaires           */
  int  n;								/* indice pour le tableau "tri"     */
  int  tot = factx * facty * factz;	/* taille de la fenetre             */
  int *tri;							/* tableau trie des valeurs         */
  int  nbpix;							/* nombre de pixels dans "tri"      */
  int  ii, jj, kk;					/* fenetre de parcours              */
  int  tmp;
  TYPEMRI3D *pix;						/* ligne a parcourir                */

  tri = CALLOC(tot, int);

  /* -- calcul des dimensions de la nouvelle image ----------------------- */

  imx_copie_param_3d_p(imdeb, imres);

  imres->width = imdeb->width / factx 
    + (imdeb->width % factx > 0 ? 1 : 0);
  imres->height = imdeb->height / facty 
    + (imdeb->height % facty > 0 ? 1 : 0);
  imres->depth = imdeb->depth / factz 
    + (imdeb->depth % factz > 0 ? 1 : 0);

  imres->x3dc = imres->width / 2;
  imres->y3dc = imres->height / 2;
  imres->z3dc = imres->depth / 2;

  imres->dx = imres->dx * factx;
  imres->dy = imres->dy * facty;
  imres->dz = imres->dy * factz;

  /* -- calcul des pixels ------------------------------------------------ */

  for(i = resi = 0 ; i < imdeb->width ; i += factx, resi++)
    for(j = resj = 0 ; j < imdeb->height ; j += facty, resj++)
      for(k = resk = 0 ; k < imdeb->depth ; k += factz, resk++)
	{
	  switch(meth)
	    {
	    case REDUC_JUMP:	/* un pixel sur n */
	      res = imdeb->mri[i][j][k];
	      break;
		
	    case REDUC_MEAN:	/* moyenne des pixels */
	      v = 0;
	      for(ii = 0 ; ii < factx ; ii++)
		for(jj = 0 ; jj < facty ; jj++)
		  for(kk = 0 ; kk < factz ; kk++)
		    v += imdeb->mri[i + ii][j + jj][k + kk];
	      res = v / tot;
	      break;
			
	    case REDUC_MEDIAN:	/* valeur mediane des pixels */
	      nbpix = 0;

	      for(ii = 0 ; ii < factx ; ii++)
		for(jj = 0 ; jj < facty ; jj++)
		  {
		    pix = imdeb->mri[i + ii][j + jj];
		    for(kk = 0 ; kk < factz ; kk++)
		      {
			v = pix[k + kk];

			/* insertion avec tri de la valeur dans le tableau   */

			for(n = 0 ; n < nbpix ; n++)
			  {
			    /* si la valeur stockee est plus grande que la   */
			    /* valeur a inserer, on echange les deux et on   */
			    /* continue le parcours                          */

			    if (tri[n] > v)
			      {

				tmp = v;
				v = tri[n];
				tri[n] = tmp;
			      }
			  }

			/* stockage de la valeur conservee dans la variable  */
			/* (la plus grande) en derniere position             */

			tri[n] = v;
			nbpix++;
		      }
		  }		/* fin boucle ii et jj */

	      res = tri[tot / 3];			/* recuperation de la valeur mediane */
	      break;
					
	    case REDUC_MIP:		/* valeur du pixels d'intensite maximale */
	      res = 0;
	      for(ii = 0 ; ii < factx ; ii++)
		for(jj = 0 ; jj < facty ; jj++)
		  for(kk = 0 ; kk < factz ; kk++)
		    res = MAXI(res, imdeb->mri[i + ii][j + jj][k + kk]);
	      break;
	    }

	  imres->mri[resi][resj][resk] = res;

	  if (once == 0)
	    {
	      mini = maxi = res;
	      once = 1;
	    }
	  else
	    {
	      mini = MINI(mini, res);
	      maxi = MAXI(maxi, res);
	    }
	
	}	/* fin boucle i, j et k */

  imres->cutoff_min = mini;
  imres->cutoff_max = maxi;

  free(tri);
			
  return 0;
}

#ifdef __GTK_GUI
/********************************************************
 * getPosFile()
 * permet de se placer au bon offset pour la lecture d'un
 * ficher 3D en raw data
 *******************************************************/
long getPosFile(const char *file, int numfirst, int filewdth, int filehght, int filedpth, int bytesppixel, int format)
{
  char* m_int=(char*)malloc(2),tmp;	
  int dst,i;
  long taille,offset;
  FILE *fp=fopen(file,"rb");
   

  switch(format) {

  /* Retourne l'offset car il y a un header dans le fichier NIFTI */
#ifdef HAS_NIFTI
  case NIFTI_GZ_FORMAT :
  case NIFTI_FORMAT :
	return 352 + nifti_extension_size(nifti_image_read(file, 0));
	break;
#endif
    /*  Recherche de l'offset car header DICOM avant raw data */
  case DICM3DBIMP_FORMAT :
  case DICM3DBEXP_FORMAT :
  case DICM3DLIMP_FORMAT :
  case DICM3DLEXP_FORMAT :
    /* Recherche taille et offset des donne du MetaElementGroup */
    taille=Acrnema_getGroupSizeandOffset(file, 0x7fe0, 0x10, &offset, format);
    return(offset);
    break;

    /*  Recherche de l'offset car header ECAT avant raw data */
  case ECAT7_FORMAT :
    if(numfirst<32){
      fseek(fp,0x206+numfirst*0x10,SEEK_SET);
      fread(m_int,2,1,fp);
      if(!_is_bigEndian){
	tmp=m_int[0];
	m_int[0]=m_int[1];
	m_int[1]=tmp;
      }
      dst=*(unsigned short*)m_int;
    }else{
      dst=1;
      for(i=1;i<(numfirst/32);++i){
	fseek(fp,0x1F6+0x200*dst,SEEK_SET);
	fread(m_int,2,1,fp);
	if(!_is_bigEndian){
	  tmp=m_int[0];
	  m_int[0]=m_int[1];
	  m_int[1]=tmp;
	}
	dst=*(unsigned short*)m_int;
      }
      fseek(fp,(numfirst%32)*0x10,SEEK_CUR);
      fread(m_int,2,1,fp);
      if(!_is_bigEndian){
	tmp=m_int[0];
	m_int[0]=m_int[1];
	m_int[1]=tmp;
      }
      dst=*(unsigned short*)m_int;
    }
    free(m_int);
    fclose(fp);
    return(0x200*(dst+1));
    break;
	
    /*  Recherche de l'offset si plusieurs images dans le fichier */
  default :	
    free(m_int); 
    fclose(fp);
    return ((numfirst-1)*filewdth*filehght*filedpth*bytesppixel);
  }

}

/* temporary, this is just for testing */
/* #ifndef HAS_ANNOTATIONS */
/* void draw_annotations (grphic3d *image){;} */
/* void read_annotations (grphic3d *image,const char *filename,int imageNumber){;} */
/* void write_annotations(grphic3d *image,HEADER *header,      int imageNumber){;} */
/* void annotations_test(){;} */
/* void add_annotation(){;} */
/* void delete_annotation(){;} */
/* void edit_annotation(){;} */
/* #endif // HAS_ANNOTATIONS */


#endif /* __GTK_GUI */


/*
**  read_raw_float_3d 
**    Load a raw float image
**   
*/
void read_raw_float_3d(const char* file, grphic3d* img, int endianess)
{
  int i,j,k;
  float min = INT_MAX, max = INT_MIN;
  float *end,*data;
  FILE* input = fopen(file, "rb");
  float* buffer = malloc(sizeof(float) * img->depth * img->width * img->height);
  fread(buffer, sizeof(float), img->depth  * img->width * img->height, input);
  
  
  end = &buffer[img->depth * img->width * img->height];
 
  for (data = buffer; data != end; data++)
    {
      if (*data < min) min = *data;
      if (*data > max) max = *data;
    }
    
  min = 1 / (max - min) * MAXMRI3D;
  for (data = buffer, i = 0; i < img->width; i++)
    for (j = 0; j < img->height; j++)
      for (k = 0; k < img->depth; k++, data++)
	img->mri[i][j][k] = (TYPEMRI3D)rint((*data) * min);  

  imx_brukermax_3d(max, min, img);
  free(buffer);
  fclose(input);
}


