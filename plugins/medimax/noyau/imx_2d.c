/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*
*** imx_2d.c
*** Login : <wiest@localhost.localdomain>
*** Started on  Wed Jan 14 09:09:59 2004 Nicolas Wiest
*** $Id: imx_2d.c,v 1.7 2007/06/11 11:25:04 bonnet Exp $
*** 	
***	Copyright (c) 2004, ULP-IPB Strasbourg.
***	All rights are reserved.
***
*/

#include <config.h>
#include <errno.h> 
#include <string.h>

#include "imagix.h" 
#include "noyau/gtkimagix.h"
#include"noyau/imx_2d.h"
#include "noyau/imx_lang.h"
#ifdef __GTK_GUI
#include "noyau/gui/imx_logical_cmap.h"
#include "noyau/gui/imx_picture_common.h"
#endif /* __GTK_GUI */
#include "noyau/gui/imx_picture2d.h"
#include "noyau/gui/imx_picture3d.h"
#ifdef __GTK_GUI
#include "noyau/gui/imx_text.h"
#include "outils/imx_misc.h"
#include "noyau/gui/interface.h"
#include "noyau/imx_zone.h"
#include "noyau/gui/imx_cmap.h"

//void (*IMXPict2d_ShowBar)(int noCMap, ptrColormap cmap, UCHAR *pDest);
#endif /* __GTK_GUI */

/* -- ptr_img() -------------------------------------------------------    
 */
/*!  Permet d'obtenir le pointeur d'une structure s_grphic Ã partir
**  du numero de l'image (nr)
**	\param nr : numero de l'image
**	\retval : un pointeur sur la structure grphic concernee, NULL en cas d'echec
*/
grphic *ptr_img(int nr)
{
	#ifdef __GTK_GUI
	if(abs(nr) > 0 && abs(nr)<=MAX_PICTURES ) 
	{
		if(nr > 0) 
		return((grphic *)&_imx_img[nr-1]);
		else
		return ptr_mask(abs(nr));
	}
	
	/*
		if(nr==0)  return((grphic *)&_imx_roi);
	*/
	if(abs(nr)>MAX_PICTURES)
	{  
		PUT_ERROR("BAD_POSITION nr >MAX_PICTURES in ptr_img");
		return((grphic *)NULL);
	}
	
	#endif /* __GTK_GUI */

	return NULL;
	
}

/***************************************************
 *** --  cr_grphic()   ---------------------------
 *** 
 ***  Permet de creer un zone memoire grphic
 ***  c'est a dire un image temporaire avec
 ***  un calloc
 ***  si temp est un pointeur d'image l'image cree prendra 
 ***  les caracteristiques de cette image.
 ***
 ***  Diffenrece avec 3d. La taille de l'image est toujours MAX_WIDTH, MAX_HEIGHT
 ***  >> LAST USER ACTION BY R. CHEVRIER on Aug, 21st 2001
 ***  >> Ajout de la creation du masque correspondant a la ROI.
 ****************************************************/

grphic *cr_grphic(grphic *temp)
{
	grphic *ia;
	
	#ifdef __GTK_GUI
	ptrColormap cmap = NULL;
	
	if(list_nbitems(&_lstColormaps))
	{
		list_movefirst(&_lstColormaps);
		cmap = list_item(&_lstColormaps);
	}
	#endif /* __GTK_GUI */
	
	if((ia = CALLOC(1, grphic)) == NULL)
	{
		PUT_ERRFAT("Can't allocate s_grphic in cr_grphic()");
		return((grphic *) NULL);  /* Cannot alloc... */
	} 
	
	/*
		#ifdef __GTK_GUI
	ia->displayable_image_data = CALLOC(_pic_size_width * _pic_size_height * _nBytesPerPixel, UCHAR);
	#endif *//* __GTK_GUI */
	
	if (temp == NULL) 
	{
		ia->bitppixel   = 8;
		ia->ref_count   = 1;
		ia->img_type    = IMAGE_NORMAL;
		ia->mask_type   = MASK_BINARY;
		ia->width       = MAX_WIDTH;
		ia->height      = MAX_HEIGHT;
		ia->icomp       = 0; 
		ia->rcoeff      = 1; 
		ia->cutoff_min  = 0;
		ia->pos         = 99;
		ia->zoom        = 1;
		#ifdef __GTK_GUI
		ia->cmapinfo.noColormap  = 1;
		ia->cmapinfo.dspStyle    = DSP_POSITIVE;
		ia->cmapinfo.dspOverBlancking = OVERBLANCKING_NONE;
	
		ia->cmapinfo.cmap        = cmap;
		ia->bar                  = 0;
		ia->indexes     = TEXT_NONE;
		#endif /* __GTK_GUI */
		ia->mask                 = NULL;
		ia->mask_active          = FALSE;
	
		#ifdef __GTK_GUI
			IMXColor_Colormap_AddRef(1);
		#endif /* __GTK_GUI */
	}
	else 
	{
		#ifdef __GTK_GUI
		ia->bar         = 0;
		ia->indexes     = TEXT_NONE;
		#endif /* __GTK_GUI */
		ia->ref_count   = 1;
		ia->img_type    = temp->img_type;
		ia->mask_type   = temp->mask_type;
		ia->width       = temp->width;
		ia->height      = temp->height;
		ia->icomp       = temp->icomp; 
		ia->rcoeff      = temp->rcoeff; 
		ia->max_pixel   = temp->max_pixel;
		ia->min_pixel   = temp->min_pixel;
		ia->cutoff_max  = temp->cutoff_max;
		ia->cutoff_min  = temp->cutoff_min;
		ia->zoom        = temp->zoom;                                        
		ia->zoom_x      = temp->zoom_x;
		ia->zoom_y      = temp->zoom_y;
		ia->bitppixel   = temp->bitppixel;
		ia->pos         = 99;
		#ifdef __GTK_GUI
		ia->cmapinfo.noColormap  = temp->cmapinfo.noColormap;
		ia->cmapinfo.dspStyle    = temp->cmapinfo.dspStyle;
		ia->cmapinfo.dspOverBlancking = temp->cmapinfo.dspOverBlancking;
		
		ia->cmapinfo.cmap        = temp->cmapinfo.cmap;
		ia->bar                  = temp->bar;
		ia->mask                 = NULL;
		//B.MARAT+JP ARMSPACH 20/06/2002 bug si temp->displayable_image_data==NULL
		//memcpy( ia->displayable_image_data, temp->displayable_image_data, _pic_size_width * _pic_size_height * _nBytesPerPixel);
		
		IMXColor_Colormap_AddRef(temp->cmapinfo.noColormap);
		#endif /* __GTK_GUI */
	}
	
	if ((ia->mri = cr_mri(ia->width, ia->height)) == NULL) 
	{
		PUT_ERRFAT("Can't allocate mri in cr_grphic()");
		return(NULL);
	}
	
	return(ia);
}

/***************************************************
 *** --  cr_grphic_modif()   ---------------------------
 *** 
 ***  Permet de creer un zone memoire grphic
 ***  en saissant ses parametres et en tenant compte 
 ***  de width, height pour la creation de mri 
 ***  
 ***   
 ***  
 ****************************************************/
grphic *cr_grphic_modif(int width, int height, float icomp, float rcoeff, int maxpixel)
{
	grphic *ia;
	#ifdef __GTK_GUI
	ptrColormap cmap = NULL;
	
	if (list_nbitems(&_lstColormaps))
	{
		list_movefirst(&_lstColormaps);
		cmap = list_item(&_lstColormaps);
	}
	#endif /*__GTK_GUI */
	
	if ((ia = CALLOC(1, grphic)) == NULL)
	{
		PUT_ERRFAT("Can't allocate s_grphic in cr_grphic()");
		return((grphic *) NULL);  /* Cannot alloc... */
	} 
		
	ia->width       = width;
	ia->height      = height;
	ia->icomp       = icomp; 
	ia->rcoeff      = rcoeff; 
	ia->max_pixel   = maxpixel;
	ia->min_pixel   = 0;
	ia->cutoff_max  = ia->max_pixel;
	ia->cutoff_min  = ia->min_pixel;
	ia->zoom        = 1;  				      
	ia->zoom_x      = 0;
	ia->zoom_y      = 0;
	ia->pos         = 99;
	ia->ref_count   = 1;
	#ifdef __GTK_GUI
	ia->cmapinfo.noColormap  = 1;
	ia->cmapinfo.cmap        = cmap;
	ia->cmapinfo.dspStyle    = DSP_POSITIVE;
	ia->mask_active          = FALSE;
	ia->bar         = 0;
	#endif /*__GTK_GUI */
	
	if ((ia->mri = cr_mri(ia->width, ia->height)) == NULL) 
	{
		PUT_ERRFAT("Can't allocate mri in cr_grphic()");
		FREE(ia);
		return((grphic *)NULL);
	}
	
	return(ia);
}

/*** free_mri  ***********************************************************
 **
 **  Libere la memoire creer par cr_mri
 **
 *************************************************************************/
void free_mri(long int **mri)
{
  FREE(mri[0]);
  FREE(mri);
}


/*
**  -- free_grphic() ----------------------------------------------------------
**
**  Libere la memoire creer par cr_grphic
**
*/

void free_grphic(grphic *temp)
{
  free_mri(temp->mri);
  free_mask_p(temp);
  FREE(temp); 
}  


#ifdef __GTK_GUI
/***************************************************
 *** --  cr_grphic_array()   ---------------------------
 *** 
 ***  Permet de creer un tableau de grphic
 ***  c'est a dire num  image temporaire avec
 ***  
 ***

****************************************************/
grphic *cr_grphic_array(int num, int width, int height)
{
	grphic  *tmp;
	grphic  *modele;
	int     i;
	
	modele              = cr_grphic(NULL);
	modele->width       = width;
	modele->height      = height;
	modele->ref_count   = 1;
	modele->bar         = 0;
	modele->indexes     = TEXT_NONE;
	
	if ((tmp = (grphic*)calloc(num, sizeof(grphic))) == NULL)
	{
		PUT_ERRFAT("Can't allocate s_grphic in cr_grphic_array()");
			
		return((grphic *)NULL);  // Cannot alloc...
	}
	
	for (i = 0 ; i < num ; i++)
	{
		if ((tmp[i].mri = cr_mri(width, height)) == NULL) 
		{
			FREE(tmp);
			PUT_ERRFAT("Can't allocate s_grphic in cr_grphic_array()");
			return((grphic *)NULL);  
		}
		
		imx_copie_p(modele, &tmp[i]);
	}
	
	free_grphic(modele);
	
	return tmp;
}
#endif /* __GTK_GUI */


/*
**  -- free_grphic_array() ----------------------------------------------------
**
**  Libere la memoire creer par cr_grphic_array
**
*/

void free_grphic_array(grphic *array, int n)
{
	int i;
		
	for(i = 0 ; i < n ; i++)  
	{
		free_mri(array[i].mri);
		free_mask_p(&array[i]);
	}
		
	FREE(array);
}


/***************************************************
 *** --  cr_mri()   ---------------------------
 *** 
 ***  Permet de creer le tableau mri pour la 
 **   structure grphic.
 ***  
 ***
 ****************************************************/
TYPEMRI **cr_mri(int width, int height)
{
	TYPEMRI	**t3;
	int	i;
		
	if( (t3=(TYPEMRI**)calloc((width+1),sizeof(TYPEMRI*))) == NULL)
		return (TYPEMRI**)NULL;
		
	if( (t3[0]=(TYPEMRI*)calloc((width*height),sizeof(TYPEMRI))) == NULL)
	{
		FREE(t3);
		return (TYPEMRI**)NULL;
	}
		
	for(i=1 ; i<width ; i++)
		t3[i]=t3[i-1]+height; 
		
	return t3;
}


//
//  ---- ptr_mask() -----------------------------------------------------------
//
//  Description:
//      ...
//
//  >>  Created on Sep 06th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//

grphic *ptr_mask(int pos)
{
  char		str[80];
  grphic	*img = NULL;

  if(pos < 1 || pos > MAX_PICTURES) 
    {
      sprintf(str, "ptr_mask(): cannot get mask at #%d", pos);
      PUT_ERROR(str);
      return(0);
    }
  img = ptr_img(pos);
  return ptr_mask_p(img);
}

/*
** ---- ptr_mask_p() ---------------------------------------------------------
*/

grphic	*ptr_mask_p(grphic* img)
{
  char		str[80];
  grphic	*mask = NULL;

  if (!img)
    {
      sprintf(str, "ptr_mask_p(): cannot get mask ");
      PUT_ERROR(str);
      return(0);    
    }
  if((mask = img->mask) == NULL)
    {
      img->mask = cr_grphic(NULL);
      mask = img->mask;
      mask->rcoeff = 1;
    }
  
  return mask;
}

/*
** ---- ptr_has_mask() -------------------------------------------------------
*/
int ptr_has_mask(int nr)
{
  grphic *img;
  int retVal=0;
  
  img = ptr_img(abs(nr));
  if (img && img->mask)
    retVal=1;
  return retVal;
}

/*
** ---- ptr_has_mask_p() -----------------------------------------------------
*/

int ptr_has_mask_p(grphic *img)
{
  if (img && img->mask)
    return 1;
  return 0;
}


/* FIXME manque par rapport a la 3d
** ---- cr_mask_ptr() --------------------------------------------------------
*/

/* FIXME manque par rapport a la 3d
** ---- init_mask_ptr() ------------------------------------------------------
*/

/*---- free_mask_p() ---------------------------------------------------*/
/*!
//  Description: libere la zone memoire du mask d'une image 2d
//  
//
//  >>  Created on Sep 06th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
// \param im : l'image dont on veut liberer le mask (E/S)
// \remark : ne libere effectivement la memoire que si mask->ref_count vaut 0
// et que mask->pos vaut -1.
*********************************/
void free_mask_p(grphic *im)
{
  ptr_grphic mask = im->mask;

  if(!mask) return;
    
  if(mask->ref_count > 0) 
    mask->ref_count--;
    
  if(mask->ref_count == 0 && mask->pos == -1)
    {
      free_grphic(mask);
      im->mask = NULL;
      im->mask_active = FALSE;
    }
  else
    {
      im->mask = NULL;
      im->mask_active = FALSE;	
    }
}

/*
**  -- info_mri_2d() ----------------------------------------------------------
**
**  Show information for a 2D picture
**
**  >>  Created by: on
**  >>  Last user action by: F. BOULEAU on Sep 29th, 00
**      - added informations: binded colormaps and display mode
**  >>  Last user action by: on
**
*/

void imx_info_mri_2d(int nr)
{ 
#ifdef __GTK_GUI
  GtkWidget   *window;
  GtkWidget *str;
    
  grphic      *p;
  char        msg[4096];
  ptrColormap cmap = NULL;

  if((p = ptr_img(nr)) == NULL)
    return ;

  window = gtk_window_new( GTK_WINDOW_TOPLEVEL );
  gtk_window_set_transient_for( GTK_WINDOW(window), GTK_WINDOW(_wnd_main));

  if(p->cmapinfo.noColormap)
    {
      list_moveto(&_lstColormaps, p->cmapinfo.noColormap);
      cmap = (t_Colormap *)list_item(&_lstColormaps);
    }

  sprintf(msg, "Fenetre 3D no %d\n", abs(nr));
  sprintf(msg,"%s\n", msg);

  sprintf(msg, "%sType d'image : %s\n", msg, IMXPict_DescImageType(p->img_type));
  sprintf(msg, "%sDimensions (pixels) : %dx%d\n", msg, p->width, p->height);
  sprintf(msg, "%sBits par pixel : %d\n", msg, p->bitppixel);
  sprintf(msg, "%sZoom : %d\n", msg, p->zoom);
  sprintf(msg, "%sOffsets : x=%d, y=%d\n", msg, p->zoom_x, p->zoom_y);
  sprintf(msg, "%sCutoff : min=%ld max=%ld\n", msg, p->cutoff_min, p->cutoff_max);
  sprintf(msg, "%sValeurs de MRI : min=%ld max=%ld\n", msg, p->min_pixel, p->max_pixel);
  sprintf(msg, "%sCoefficient de valeur reelle : rcoeff=%f icomp=%f\n", msg, p->rcoeff, p->icomp);
  sprintf(msg, "%sResolution (mm) : dx=%f dy=%f dz=%f\n", msg, p->dx, p->dy, p->dz);
  sprintf(msg,"%s\n", msg);
    
  sprintf(msg, "%sPalette de couleurs : %s (no %d, id %d)\n", msg, cmap->szName, p->cmapinfo.noColormap, cmap->noMode);
  sprintf(msg, "%sMode d'affichage : %s\n", msg, imx_GetStrDisplayMode(p->cmapinfo.dspStyle));
  sprintf(msg, "%sOverBlancking : %s\n", msg, imx_GetStrOverBlanckingMode(p->cmapinfo.dspOverBlancking));
  sprintf(msg, "%sType de palette : %s\n", msg, IMXColor_DescCMapType(cmap->noType));
    
//USELESS
//  if(cmap->noType == MAP_FIXED)
//    sprintf(msg, "%sType de palier : %s\n", msg, IMXColor_DescStageType(cmap->noStageType));
    
  sprintf(msg,"%s\n", msg);

  sprintf(msg, "%sMasque : %s (%s)\n", msg, IMXPict_DescMaskType(nr), (p->mask_active ? "active" : "desactive"));
  sprintf(msg,"%s\n", msg);

  sprintf(msg, "%sAffichage : %s %s\n", msg, (p->bar ? "palette" : "-"), (p->bar && p->indexes != TEXT_NONE ? "et valeurs" : ""));
  sprintf(msg, "%sReferences : %d\n", msg, p->ref_count);
  sprintf(msg,"%s\n", msg);

  sprintf(msg, "%sNom du fichier : %s\n", msg, p->filename);
    
  if(p->type)
    sprintf(msg, "%sType : %c\n", msg, p->type);
    
  sprintf(msg,"%s\n", msg);

  sprintf(msg, "%sPatient : %s (%s)\n", msg, p->patient.name, p->patient.d_birth);
  sprintf(msg, "%sMedecin : %s\n", msg, p->patient.medecin);
  sprintf(msg, "%sDate d'examen : %s\n", msg, p->patient.d_examen);
  sprintf(msg,"%s\n", msg);

  str = gtk_label_new( msg );
  gtk_label_set_justify(GTK_LABEL(str), GTK_JUSTIFY_LEFT);
  gtk_container_add( GTK_CONTAINER(window), str );

  gtk_widget_show_all( window );

#endif /* __GTK_GUI */
}


//
//  ---- IMXPict2d_CleanPicture() -----------------
//
//  Cleans a 2D picture (developped for Therese =)
//  Created by R. CHEVRIER on Dec 6th, 2001
//

void IMXPict2d_CleanPicture( void ) {
	int pos = GET_PLACE( "Picture to reset" );

	IMXPict2d_ResetPicture( pos );

	show_picture( pos );
}



//
//  ---- IMXPict2d_ResetPicture() -----------------
//
//  Resets a 2D MRI and unactivates mask.
//  Created by R. CHEVRIER on Dec 6th, 2001
//

void IMXPict2d_ResetPicture( int pos ) {
	grphic *img = ptr_img( pos );
	int i, j;

	for( i=0; i<img->width; i++ )
		for( j=0; j<img->height; j++ )
			img->mri[i][j] = 0l;

	img->mask_active = FALSE;
	img->zoom = 1;
	img->zoom_x = img->zoom_y = 0;
}


/*******************************************
** --  imx_brukermax() ---------------------
* Sous programme determine MAXPIXEL et ICOMP a partir de AMAX
* MAXPIXEL sera le plus proche possible de 4194304
********************************************/
int imx_brukermax(float amax, int *maxpixel, float *icomp, float *rcoeff)
{ 
        double aicomp,rcoeff2,two;

/*
*       Sous programme determinent MAXPIXEL et ICOMP a partir de AMAX
*       MAXPIXEL sera le plus proche possible de 4194304
*/
       two=2;
      *maxpixel=4194304;
      
      if (amax==0) aicomp=0.;
      else {
    rcoeff2=amax/(*maxpixel);
           aicomp=log(rcoeff2)/log(two);
    }
      
      if (aicomp<=0) *icomp=(float)((int)aicomp);
      else *icomp=(float)((int)(aicomp+1));
      
      aicomp=*icomp;
      *rcoeff=pow(two,aicomp);
      *maxpixel=amax/(*rcoeff);
      return(1);
}

/* FIXME  manque par rapport a la 3d
** ---- imx_iniimg3d() -------------------------------------------------------
*/

/*******************************************
** --  imx_inimaxpixel() ---------------------
**  Sous programme 
**   -  Cheche maxpixel dans l'image
**
********************************************/
int imx_inimaxpixel(int im_1)
{ 
 grphic *im1;

/*
*       Recherche du maximum vrai
*/
  im1=ptr_img(im_1);
  imx_inimaxpixel_p(im1);

  return(1);
}
/*******************************************
** --  imx_inimaxpixel_p() ---------------------
**  Sous programme 
**   -  Cheche maxpixel dans l'image
**
********************************************/
int imx_inimaxpixel_p(grphic *im1)
{ 
 long maxpixel;
 int i,j;

/*
*       Recherche du maximum vrai
*/

  maxpixel=0;
  for(i=0;i<im1->width;i++)
    for(j=0;j<im1->height;j++)
      {
      maxpixel=MAXI(maxpixel,im1->mri[i][j]);
/*      if (im1->mri[i][j]>maxpixel)
  maxpixel=im1->mri[i][j];*/
      }

  im1->max_pixel=maxpixel;
  im1->cutoff_max=maxpixel;

  return(1);
}

/*******************************************
** --  imx_iniminpixel() ---------------------
**  Sous programme 
**   -  Cheche maxpixel dans l'image
**
********************************************/
int imx_iniminpixel(int im_1)
{ 
 grphic *im1;

/*
*       Recherche du maximum vrai
*/
  im1=ptr_img(im_1);
  imx_iniminpixel_p(im1);

  return(1);
}
/*******************************************
** --  imx_iniminpixel_p() ---------------------
**  Sous programme 
**   -  Cheche maxpixel dans l'image
**
********************************************/
int imx_iniminpixel_p(grphic *im1)
{ 
 long minpixel;
 int i,j;

/*
*       Recherche du maximum vrai
*/

  minpixel=1000000000;
  for(i=0;i<im1->width;i++)
    for(j=0;j<im1->height;j++)
      minpixel=MINI(minpixel,im1->mri[i][j]);

  im1->min_pixel=minpixel;

  return(1);
}

/*******************************************
** --  imx_inimaxminpixel() ---------------------
**  Sous programme 
**   -  Cheche maxpixel et minpixel dans l'image
**
********************************************/
int imx_inimaxminpixel(int im_1)
{ 
 grphic *im1;

/*
*       Recherche du maximum et minimum 
*/
  im1=ptr_img(im_1);
  imx_inimaxminpixel_p(im1);

  return(1);
}

/*******************************************
** --  imx_inimaxminpixel_p() ---------------------
**  Sous programme 
**   -  Cheche maxpixel et minpixel  dans l'image
**
********************************************/
int imx_inimaxminpixel_p(grphic *im1)
{ 
 long maxpixel,minpixel;
 int i,j;

/*
*       Recherche du maximum vrai
*/

  maxpixel=0;
  minpixel=1000000000;
  for(i=0;i<im1->width;i++)
    for(j=0;j<im1->height;j++) {
      maxpixel=MAXI(maxpixel,im1->mri[i][j]);
      minpixel=MINI(minpixel,im1->mri[i][j]);
      }

  im1->max_pixel=maxpixel;
  im1->min_pixel=minpixel;
  im1->cutoff_max=maxpixel;
  im1->cutoff_min=minpixel;

  return(1);
}

/*******************************************
** --  imx_iniparam() ---------------------
**  Sous programme 
**   -  Initialisation des parametres
**      max-pixel, min_pixel, rcoeff, icomp
**
********************************************/
int imx_iniparam(int im_1, long int max_pixel, long int min_pixel, float icomp, float rcoeff)
{ 
 grphic *im1;

  im1=ptr_img(im_1);

  im1->max_pixel=max_pixel;
  im1->min_pixel=min_pixel;
  im1->cutoff_max=max_pixel;
  im1->cutoff_min=min_pixel;
  im1->icomp=icomp;
  im1->rcoeff=rcoeff;

  return(1);
}
/*******************************************
** --  imx_iniparaimg() ---------------------
**  Sous programme 
**   -  Initialisation des parametres
**      max-pixel, min_pixel, rcoeff, icomp
**    avec recherche de max_pixel et min_pixel
**
********************************************/
int imx_iniparaimg(int im_1)
{ 
  grphic *im1;

  im1=ptr_img(im_1);
  imx_iniparaimg_p(im1);

  return(1);
}

/*******************************************
** --  imx_iniparaimg_p() ---------------------
**  Sous programme 
**   -  Initialisation des parametres
**      max-pixel, min_pixel, rcoeff, icomp
**    avec recherche de max_pixel et min_pixel
**
********************************************/
int imx_iniparaimg_p(grphic *im1)
{ 

  imx_inimaxminpixel_p(im1);
  im1->icomp=0;
  im1->rcoeff=1;

  return(1);
}
/*******************************************
** --  imx_copie_param() ---------------------
**  Sous programme 
**   -  Copie les parametres d'une image dans l'autre
**
********************************************/
int imx_copie_param(int im_1, int im_res)
{ 
 grphic *im1,*imres;

  im1=ptr_img(im_1);
  imres=ptr_img(im_res);
  imx_copie_param_p(im1,imres);

  return(1);

}
/*******************************************
** --  imx_copie_param_p() ---------------------
**  Sous programme 
**   -  Copie les parametres d'une image dans l'autre
**      en utilisant les pointeur d'image
********************************************/
int imx_copie_param_p(grphic *im1, grphic *imres)
{ 
  imres->max_pixel    = im1->max_pixel;
  imres->min_pixel    = im1->min_pixel;
  imres->cutoff_max   = im1->cutoff_max;
  imres->cutoff_min   = im1->cutoff_min;
  imres->icomp        = im1->icomp;
  imres->rcoeff       = im1->rcoeff;
  imres->height       = im1->height;
  imres->width        = im1->width;
  imres->zoom         = im1->zoom;
  imres->zoom_x       = im1->zoom_x;
  imres->zoom_y       = im1->zoom_y;
  imres->dx           = im1->dx;
  imres->dy           = im1->dy;
  imres->dz           = im1->dz;
  imres->type         = im1->type;
  imres->ref_count    = 1;
  imres->mask_active  = im1->mask_active;
#ifdef __GTK_GUI
  imres->img_type     = im1->img_type;
  imres->mask_type    = im1->mask_type;
  imres->bitppixel    = im1->bitppixel;  
  imres->cmapinfo.noColormap = im1->cmapinfo.noColormap ;  // Colormap index
  imres->cmapinfo.dspStyle    =   im1->cmapinfo.dspStyle    ;    // colormap curve style to use
  imres->cmapinfo.dspOverBlancking = im1->cmapinfo.dspOverBlancking;
  imres->cmapinfo.cmap     =   im1->cmapinfo.cmap     ;
  imres->bar       =   im1->bar   ;
  imres->indexes      =   im1->indexes  ;
#endif /* __GTK_GUI */
  strcpy(imres->filename, im1->filename);
  strcpy(imres->patient.name, im1->patient.name);
  strcpy(imres->patient.d_birth, im1->patient.d_birth);
  strcpy(imres->patient.medecin, im1->patient.medecin);
  strcpy(imres->patient.d_examen, im1->patient.d_examen);

  return(1);
}
