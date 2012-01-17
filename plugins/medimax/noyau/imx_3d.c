/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/********************************************************************************
**/	
/*!	   \file:        imx_3d.c
***
***    project:    Imagix 2.01
***
***
***    \brief description:    creation/destruction structure grphic3d 
***
***
***---------------------------------------------------------------------*/


#include <config.h>
#include <string.h> 	/* Needed for strcpy	*/




#include "imagix.h"
#include "noyau/gtkimagix.h"
#ifdef __GTK_GUI
#include "noyau/gui/imx_logical_cmap.h"
#include "noyau/gui/imx_picture_common.h"
#include "noyau/gui/interface.h"
#endif // __GTK_GUI
#include "noyau/imx_3d.h"
#include "noyau/gui/imx_picture3d.h"
#include "noyau/zone_3d.h"
#include "noyau/imx_tmp_grphic3d.h"

extern int IMXColor_BindColormap_3d_ptr(ptr_grphic3d ptrPic, int nNewCMap);



/* -- ptr_img_3d() -------------------------------------------------------
 */
/*!
**  Permet d'obtenir le pointeur d'une structure s_grphic3d �partir
**  du numero de l'image (nr)
**  \param nr : le numero de l'image
**  \retval un pointeur vers la structure grphic3d concernee
*/
grphic3d *ptr_img_3d(int nr)
{
#ifdef __GTK_GUI
  if(abs(nr) > 0 && abs(nr) <= MAX_PICTURES_3D)
    {
      if(nr > 0)
	return((grphic3d *)&_imx_img3d[nr-1]);
      else
	return ptr_mask_3d(abs(nr));
    }
  else
    {
      if(nr==0)
	return((grphic3d *)&_imx_roi3d);
      else
	if (nr>0)     // <- acces a une image temporaire
	  return (ptr_img_tmp_3d(nr));
	else
	  return (ptr_mask_3d(abs(nr)));
    }

  if(nr>MAX_PICTURES_3D){  
    PUT_ERROR("BAD_POSITION nr >MAX_PICTURES in ptr_img_3d");
    return((grphic3d *)NULL);
  }
  return((grphic3d *)NULL);


#else  /* __GTK_GUI */
  return NULL;
#endif /* __GTK_GUI */
}

/***************************************************
 *** --  cr_grphic3d()   ---------------------------
 **/
/*! \brief creation d'une objet grphic3d*/
/*!
  Permet de creer un zone memoire grphic3d
  c'est a dire un image temporaire avec
  un calloc
  \param temp : les parametres de temp (par ex width, height etc...) sont copies dans la structure
  nouvellement cree. Si temp vaut NULL alors la structure est cree avec la taille maximal et d'autre
  parametres par defaut.
  \retval  la structure cree si reussite , NULL si echec
  \remark  pour detruire l'objet ainsi cree on fera appel a ::free_grphic3d
*/
/****************************************************/
grphic3d  *cr_grphic3d(grphic3d *temp)
{
	grphic3d *ia;
	
	#ifdef __GTK_GUI
		ptrColormap cmap = NULL;
		
		if(list_nbitems(&_lstColormaps))
		{
			list_movefirst(&_lstColormaps);
			cmap = list_item(&_lstColormaps);
		}
	#endif /*__GTK_GUI */
		
	if((ia = CALLOC(1, grphic3d)) == NULL)
	{
		PUT_ERRFAT("Can't allocate s_grphic3d in cr_grphic3d()");
		return((grphic3d *) NULL);  /* Cannot alloc... */
	}

	if (temp == NULL)
	{
		ia->img_type	= IMAGE_NORMAL;
		ia->mask_type	= MASK_BINARY;
		ia->width		= MAX_WIDTH_3D;
		ia->height		= MAX_HEIGHT_3D;
		ia->depth		= MAX_DEPTH_3D;
		ia->icomp		= 0;
		ia->rcoeff		= 1;
		ia->pos			= 99;
		ia->zoom		= 1;
		ia->zoom_x		= 0;
		ia->zoom_y		= 0;
		ia->zoom_z		= 0;
		ia->mask		= NULL;
		#ifdef __GTK_GUI
		ia->cmapinfo.noColormap	= 1;
		ia->cmapinfo.cmap		= cmap;
		ia->cmapinfo.dspStyle	= DSP_POSITIVE;
		ia->mask_active			= FALSE;
		ia->bar					= 0;
		ia->annotations			= NULL;
		#endif /*__GTK_GUI */
/*JPA*/
printf("temp = NULL %d %d %d \n",ia->width,ia->height,ia->depth);
	}
	
	else
	{
		ia->bar         = 0;
		ia->width       = temp->width;
		ia->height      = temp->height;
		ia->depth       = temp->depth;
		ia->x3dc        = temp->x3dc;
		ia->y3dc        = temp->y3dc;
		ia->z3dc        = temp->z3dc;
		ia->icomp       = temp->icomp;
		ia->rcoeff      = temp->rcoeff;
		ia->max_pixel   = temp->max_pixel;
		ia->min_pixel   = temp->min_pixel;
		ia->cutoff_max  = temp->cutoff_max;
		ia->cutoff_min  = temp->cutoff_min;
		ia->zoom        = temp->zoom;
		ia->zoom_x      = temp->zoom_x;
		ia->zoom_y      = temp->zoom_y;
		ia->zoom_z      = temp->zoom_z;
		ia->dx          = temp->dx;
		ia->dy          = temp->dy;
		ia->dz          = temp->dz;
		ia->pos         = 99;
		ia->mask        = NULL;
		ia->img_type    = temp->img_type;
		ia->mask_type   = temp->mask_type;
		#ifdef __GTK_GUI
		ia->cmapinfo.noColormap  = temp->cmapinfo.noColormap;
		ia->cmapinfo.dspStyle    = temp->cmapinfo.dspStyle;
		ia->cmapinfo.cmap        = temp->cmapinfo.cmap;
		ia->mask_active          = FALSE;
		ia->annotations			= NULL;
		#endif /*__GTK_GUI */
	
		strcpy(ia->filename, temp->filename);
		strcpy(ia->patient.name, temp->patient.name);
		strcpy(ia->patient.d_birth, temp->patient.d_birth);
		strcpy(ia->patient.medecin, temp->patient.medecin);
		strcpy(ia->patient.d_examen, temp->patient.d_examen);
/*JPA*/
printf("%d %d %d \n",ia->width,ia->height,ia->depth);
	}
	
	if ((ia->mri = cr_mri_3d(ia->width, ia->height, ia->depth)) == NULL)
	{
		PUT_ERRFAT("Can't allocate mri in cr_grphic3d()");
		free(ia);
		return(NULL);
	}
	
	ia->mri_width		= ia->width;
	ia->mri_height		= ia->height;
	ia->mri_depth		= ia->depth;
	
	/* initialisation des champs pour la coupe curviligne */
	ia->curv.state		= 0;             /* donnees inexistantes              */
	ia->curv.nb_params	= 0;             /* courbe de degre 0                 */
	
	return(ia);
}

/***************************************************
 *** --  cr_grphic3d_modif()   ---------------------------
 **/
/*! \brief  Permet de creer un zone memoire grphic3d */
/*!  c'est a dire un image temporaire avec
  un calloc

  \param width : ...
  \param height : ...
  \param depth : ...
  \param icomp : ...
  \param rcoeff : ...
  \param maxpixel : ....
  \retval  la structure cree si reussite , NULL si echec
  \remark  pour detruire l'objet ainsi cree on fera appel a ::free_grphic3d
*/
/****************************************************/
grphic3d  *cr_grphic3d_modif(UINT width, UINT height, UINT depth, float icomp, float rcoeff, UINT maxpixel)
{
	grphic3d *ia;
	
	#ifdef __GTK_GUI
	ptrColormap cmap = NULL;
	
	if (list_nbitems(&_lstColormaps))
	{
		list_movefirst(&_lstColormaps);
		cmap = list_item(&_lstColormaps);
	}
	#endif /*__GTK_GUI */
	
	if ((ia = CALLOC(1, grphic3d)) == NULL)
	{
		PUT_ERRFAT("Can't allocate s_grphic3d in cr_grphic3d_modif()");
		return((grphic3d *) NULL);  /* Cannot alloc... */
	}
	
	ia->img_type    = IMAGE_NORMAL;
	ia->mask_type   = MASK_BINARY;
	ia->width       = width;
	ia->height      = height;
	ia->depth       = depth;
	ia->icomp       = icomp;
	ia->rcoeff      = rcoeff;
	ia->max_pixel   = maxpixel;
	ia->min_pixel   = 0;
	ia->cutoff_max  = ia->max_pixel;
	ia->cutoff_min  = ia->min_pixel;
	ia->zoom        = 1;
	ia->zoom_x      = 0;
	ia->zoom_y      = 0;
	ia->zoom_z      = 0;
	ia->pos         = 99;
	ia->mask        = NULL;
	#ifdef __GTK_GUI
	ia->cmapinfo.noColormap  = 1;
	ia->cmapinfo.cmap        = cmap;
	ia->cmapinfo.dspStyle    = DSP_POSITIVE;
	ia->mask_active          = FALSE;
	ia->bar                  = 0;
	#endif /*__GTK_GUI */
	
	if ((ia->mri = cr_mri_3d(ia->width, ia->height, ia->depth)) == NULL)
	{
		PUT_ERRFAT("Can't allocate mri in cr_grphic3d_modif()");
		free(ia);
		return(NULL);
	}
	
	ia->mri_width       = ia->width;
	ia->mri_height      = ia->height;
	ia->mri_depth       = ia->depth;
	
	/* initialisation des champs pour la coupe curviligne */
	ia->curv.state     = 0;             /* donnees inexistantes              */
	ia->curv.nb_params = 0;             /* courbe de degre 0                 */
	
	return(ia);
}


/*** free_mri_3d  **********************************************************/
/**/
/*! \brief Libere la memoire creer par cr_mri_3d*/
/*! \param mri : le tableau a liberer (E/S)
  \retval void
  \remark s'utilise rarement seul, utiliser de preferecence ::free_grphic3d
*/
/************************************************************************/
void free_mri_3d(TYPEMRI3D ***mri)
{
  if (mri)
    {
      FREE(mri[0][0]);
      FREE(mri[0]);
      FREE(mri);
    }
}

/*** free_grphic3d  *******************************************************/
/**/
/*! \brief Libere la memoire creer par ::cr_grphic3d
  \param temp : la structure a liberer (E/S)
  \retval void */
/*************************************************************************/
void free_grphic3d(grphic3d *temp)
{
	if (temp)
	{
		free_mri_3d(temp->mri);
		
		if (temp->curv.state != 0)	/* coupe curviligne ? */
		{
			FREE(temp->curv.r);
			if (temp->curv.nb_params > 0) 
        FREE(temp->curv.vect);
		}
		
		free_mask_3d_p(temp);
		FREE(temp);
	}
}

/****************************************************
 **** --  cr_grphic3d_array()   ---------------------------
 **/
/*!  \brief Permet de creer un tableau de grphic3d
  c'est a dire num  image temporaire avec
  \param num : nombre d'images a creer
  \param width : ...
  \param height : ...
  \param depth : ...
  \retval grphic3d* : un tableau de grphic3d  */
/****************************************************/
grphic3d *cr_grphic3d_array(UINT num, UINT width, UINT height, UINT depth)
{
	grphic3d	*tmp;
	unsigned int		i;
	
	#ifdef __GTK_GUI
	ptrColormap cmap = NULL;
	
	if(list_nbitems(&_lstColormaps))
	{
		list_movefirst(&_lstColormaps);
		cmap = list_item(&_lstColormaps);
	}
	#endif /*__GTK_GUI*/
	
	if( (tmp = (grphic3d*)calloc(num, sizeof(grphic3d))) == NULL)
	{
		PUT_ERRFAT("Can't allocate s_grphic3d in cr_grphic3d_array()");
		return((grphic3d *)NULL);  /* Cannot alloc... */
	}
	
	for(i=0 ; i<num ; i++)
	{
		if( (tmp[i].mri = cr_mri_3d(width, height, depth)) == NULL)
		{
			FREE(tmp);
			PUT_ERRFAT("Can't allocate s_grphic3d in cr_grphic3d_array()");
			return((grphic3d *)NULL);  /* Cannot alloc... */
		}
		
		tmp[i].width		= width;
		tmp[i].height		= height;
		tmp[i].depth		= depth;
		tmp[i].mri_width	= width;
		tmp[i].mri_height	= height;
		tmp[i].mri_depth	= depth;
		tmp[i].pos			= 99;
		tmp[i].bar			= 0;
		tmp[i].icomp		= 0;
		tmp[i].rcoeff		= 1;
		#ifdef __GTK_GUI
		tmp[i].cmapinfo.noColormap	= 1;
		tmp[i].cmapinfo.dspStyle	= DSP_POSITIVE;
		tmp[i].cmapinfo.cmap		= cmap;
		tmp[i].mask					= NULL; //cr_grphic3d(tmp[i].mask);
		#endif /*__GTK_GUI*/
		
		/* initialisation des champs pour la coupe curviligne */
		
		tmp[i].curv.state = 0;			/* donnees inexistantes              */
		tmp[i].curv.nb_params = 0;		/* courbe de degre 0                 */
	}
	
	return tmp;
}

/****************************************************
 **** --  cr_ptr_grphic3d_array()   ---------------------------
 **/
/*!  \brief Permet de creer un tableau de grphic3d* pointant sur des structure grphic3d valides
  c'est a dire num  image temporaire avec
  \param num : nombre d'images a creer
  \param width : ...
  \param height : ...
  \param depth : ...
  \retval grphic3d** : un tableau de grphic3d*  */
/****************************************************/
grphic3d **cr_ptr_grphic3d_array(UINT num, UINT width, UINT height, UINT depth)
{
	grphic3d	**tmp;
	unsigned int		i;
	
	#ifdef __GTK_GUI
	ptrColormap cmap = NULL;
	
	if (list_nbitems(&_lstColormaps))
	{
		list_movefirst(&_lstColormaps);
		cmap = list_item(&_lstColormaps);
	}
	#endif /*__GTK_GUI*/
	
	if ((tmp = (grphic3d**)calloc(num, sizeof(grphic3d*))) == NULL)
	{
		PUT_ERRFAT("Can't allocate s_grphic3d in cr_grphic3d_array()");
		return((grphic3d **)NULL);  /* Cannot alloc... */
	}
	
	for (i=0 ; i<num ; i++)
	{
		tmp[i] = cr_grphic3d_modif(width, height, depth, 0., 1.,0);
		
		if (tmp[i] == NULL) 
		{
			FREE(tmp);
			PUT_ERRFAT("Can't allocate s_grphic3d in cr_grphic3d_array()");
			return((grphic3d **)NULL);  /* Cannot alloc... */
		}
		
		tmp[i]->pos	= 99;
		tmp[i]->bar	= 0;
	
		#ifdef __GTK_GUI
		tmp[i]->cmapinfo.noColormap	= 1;
		tmp[i]->cmapinfo.dspStyle	= DSP_POSITIVE;
		tmp[i]->cmapinfo.cmap		= cmap;
		tmp[i]->mask				= NULL;
		#endif /*__GTK_GUI*/
	
		/* initialisation des champs pour la coupe curviligne */
	
		tmp[i]->curv.state = 0;			/* donnees inexistantes              */
		tmp[i]->curv.nb_params = 0;		/* courbe de degre 0                 */
	}
	
	return tmp;
}

/*** free_grphic3d_array  ************************************************
 **/
/*! \brief Libere la memoire creer par cr_grphic3d_array
  \param array : le tableau de grphic3d a liberer
  \param n : la taille du tableau
  \retval void
*/
/************************************************************************/
void free_grphic3d_array(grphic3d *array, int n)
{
	int i;
	
	for(i=0 ; i<n ; i++)
	{
		free_mri_3d(array[i].mri);
		free_mask_3d_p(&array[i]);
	}
	
	FREE(array);
}

/*** free_ptr_grphic3d_array  ************************************************
 **/
/*! \brief Libere la memoire creer par cr_grphic3d_array
  \param array : le tableau de grphic3d a liberer
  \param n : la taille du tableau
  \retval void
*/
/************************************************************************/
void free_ptr_grphic3d_array(grphic3d **array, int n)
{
	int i;
	
	for(i=0 ; i<n ; i++)
	{
		free_grphic3d(array[i]);
		free_mask_3d_p(array[i]);
	}
	
	FREE(array);
}
/***************************************************
 *** --  cr_mri_3d()   ---------------------------
 */
/*! \brief  Alloue le tableau mri pour la structure grphic3d.

\param width : ...
\param height : ...
\param depth : ...
\retval  TYPEMRI3D *** : un tableau mri, NULL si echec

****************************************************/

TYPEMRI3D ***cr_mri_3d(UINT width, UINT height, UINT depth)
{
	TYPEMRI3D *t1, ***t3;
	unsigned int       i, j;
	
	if ((t3 = CALLOC((width + 1), TYPEMRI3D**)) == NULL)
		return (TYPEMRI3D***)NULL;
	
	if ((t3[0] = CALLOC((width * height), TYPEMRI3D*)) == NULL)
	{
		FREE(t3);
		return (TYPEMRI3D***)NULL;
	}
	
	for(i = 1 ; i < width ; i++)
		t3[i] = t3[i - 1] + height;
	
	if ((t3[0][0] = CALLOC((width * height * depth), TYPEMRI3D)) == NULL)
	{
		FREE(t3[0]);
		FREE(t3);
		return (TYPEMRI3D***)NULL;
	}
	
	t1 = t3[0][0];
	
	for (i=0 ; i < width ; i++)
	{
		for(j = 0 ; j < height ; j++)
		{
			t3[i][j] =  t1;
			t1       += depth;
		}
	}
	
	return t3;
}

/*********************** Gestion du Masque *********************/

//
//  ---- ptr_mask_3d() --------------------------------------------------------
//
//  Description:
//    Retourne le pointeur du masque. S'il n'existe pas il cree le masque.
//
//  >>  Created on Sep 06th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//
grphic3d *ptr_mask_3d(int pos3d)
{
	grphic3d *img, *mask = NULL;
	char str[80];
	int img_tran,img_sagi,img_coro,img_curv;
	
	img = ptr_img_3d(pos3d);
	
	if (!img)
	{
		sprintf(str, "ptr_mask_3d(): cannot get mask at #%d", pos3d);
		PUT_ERROR(str);
		return(0);
	}
	
	if((mask = img->mask) == NULL)
	{
		img->mask = cr_grphic3d_modif(img->mri_width, img->mri_height, img->mri_depth, 0., 1.,0);
		(img->mask)->rcoeff = 1;
		(img->mask)->pos = -pos3d;
		mask = img->mask;
		
		if (pos3d > 0 && pos3d <= MAX_PICTURES_3D)
		{
			imx_iniimg3d(-pos3d,&img_coro,&img_sagi,&img_tran,&img_curv);
			(img->mask)->img_coro=img_coro;
			(img->mask)->img_sagi=img_sagi;
			(img->mask)->img_tran=img_tran;
			(img->mask)->img_curv=img_curv;
			(img->mask)->zoom=1;
		}
	}
	
	return mask;
}

//
//  ---- ptr_mask_3d_p() --------------------------------------------------------
//
//  Description:
//      ...
//
//  >>  Created on 2004 JPA
//
grphic3d *ptr_mask_3d_p(grphic3d *img)
{
	grphic3d  *mask = NULL;
	char str[80];
	int img_tran,img_sagi,img_coro,img_curv;
	
	
	if (!img)
	{
		sprintf(str, "ptr_mask_3d_p(): cannot get mask ");
		PUT_ERROR(str);
		return(0);
	}
	
	if((mask = img->mask) == NULL)
	{
		img->mask = cr_grphic3d_modif(img->mri_width, img->mri_height, img->mri_depth, 0., 1.,0);
		// NWD: Why the width/depth/height are not set??? do it!
		(img->mask)->width = img->width;
		(img->mask)->height = img->height;
		(img->mask)->depth = img->depth;
	
		(img->mask)->rcoeff = 1;
		(img->mask)->pos = img->pos;
		mask = img->mask;
		
		if (img->pos > 0 && img->pos <= MAX_PICTURES_3D)
		{
			imx_iniimg3d(img->pos,&img_coro,&img_sagi,&img_tran,&img_curv);
			(img->mask)->img_coro=-img_coro;
			(img->mask)->img_sagi=-img_sagi;
			(img->mask)->img_tran=-img_tran;
			(img->mask)->img_curv=-img_curv;
			(img->mask)->zoom=1;
		}
	}
	
	return mask;
}



//  ---- ptr_has_mask_3d() --------------------------------------------------------
//  indique si un masque existe
int ptr_has_mask_3d(int pos)
{
  grphic3d *img;
  int retVal = 0;
  img = ptr_img_3d(abs(pos));
  if (img && img->mask)
    retVal = 1;
  return retVal;
}

//  ---- ptr_has_mask_3d_p() --------------------------------------------------------
//  indique si un masque existe
int ptr_has_mask_3d_p(grphic3d *img)
{
  int retVal = 0;
  if (img && img->mask)
    retVal = 1;
  return retVal;
}

//  ---- cr_mask_3d_ptr() --------------------------------------------------------
//
grphic3d *cr_mask_3d_ptr(grphic3d *img)
{

  if( img->mask == NULL)
    {
      img->mask = cr_grphic3d_modif(img->mri_width, img->mri_height, img->mri_depth, 0., 1.,0);
      (img->mask)->rcoeff = 1;
      (img->mask)->zoom=1;
    }
  return img->mask;
}

//  ---- init_mask_3d_ptr() --------------------------------------------------------
//
grphic3d *init_mask_3d_ptr(grphic3d *img)
{
  int img_tran,img_sagi,img_coro,img_curv;

  if( img->mask == NULL)
    img->mask = cr_grphic3d(NULL);
  //(img->mask)->rcoeff = 1;  Hennequin: g��e bug �l'enregistrement ipb+mask ?
  (img->mask)->pos = img->pos;
  if (img->pos!=99) // il ne faut pas afficher qd pos=99 ???? a verifier
    {
      imx_iniimg3d(img->pos,&img_coro,&img_sagi,&img_tran,&img_curv);
      (img->mask)->img_coro=-img_coro;
      (img->mask)->img_sagi=-img_sagi;
      (img->mask)->img_tran=-img_tran;
      (img->mask)->img_curv=-img_curv;
      (img->mask)->zoom=img->zoom;
    }
  return img->mask;
}

//
//  ---- free_mask_3d_p() -----------------------------------------------------
//
/*!  Description: libere la memoire alloue au mask d'une image
//      ...
//
//  >>  Created on Sep 06th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
// \param im : l'image dont le masque va etre libere (E/S)
// \attention Si l'image est affiche il peut etre util de faire appel
//  a free_display_mask_3d_p() pour liberer la memoire utilise par l'affichage
*/
void free_mask_3d_p(grphic3d *im)
{
  ptr_grphic3d mask = im->mask;

  if(!mask) return;

  free_grphic3d(mask);
  im->mask = NULL;
  im->mask_active = FALSE;
}

/*
**  -- info_mri_3d() ----------------------------------------------------------
**/
/*! \brief Show information for a 3D picture
  \param  nr : numero de l'image 3d
  \retval void
*/
/*
**  >>  Created by: on
**  >>  Last user action by: F. BOULEAU on Sep 29th, 00
**      - added informations: binded colormaps and display mode
**  >>  Last user action by: R. CHEVRIER on Aug 29th, 2001
**      - translating to GTK.
**
*/

void imx_info_mri_3d(int nr)
{
#ifdef __GTK_GUI
  GtkWidget   *window;
  GtkLabel *str;

  grphic3d    *p;
  char        msg[4096];
  ptrColormap cmap = NULL;

  if((p = ptr_img_3d(nr)) == NULL)
    return;

  window = gtk_window_new( GTK_WINDOW_TOPLEVEL );
  gtk_window_set_transient_for( GTK_WINDOW(window), GTK_WINDOW(_wnd_main) );

  if(p->cmapinfo.noColormap)
    {
      list_moveto(&_lstColormaps, p->cmapinfo.noColormap);
      cmap = (t_Colormap *)list_item(&_lstColormaps);
    }

  sprintf(msg, "Fenetre 3D no %d\n", abs(nr));
  sprintf(msg,"%s\n", msg);

  sprintf(msg, "%sImages 2D associees : coro=%d sagi=%d tran=%d curv=%d\n", msg, p->img_coro, p->img_sagi, p->img_tran, p->img_curv);
  sprintf(msg, "%sCoordonnees des coupes : %d,%d,%d\n", msg, p->x3dc, p->y3dc, p->z3dc);
  sprintf(msg, "%s\n", msg);

  sprintf(msg, "%sType d'image : %s\n", msg, IMXPict_DescImageType(p->img_type));
  sprintf(msg, "%sDimensions (pixels) : %dx%dx%d\n", msg, p->width, p->height, p->depth);
  sprintf(msg, "%sBits par pixel : %d\n", msg, p->bitppixel);
  sprintf(msg, "%sZoom : %d\n", msg, p->zoom);
  sprintf(msg, "%sOffsets : x=%d, y=%d, z=%d\n", msg, p->zoom_x, p->zoom_y, p->zoom_z);
  sprintf(msg, "%sCutoff : min=%d max=%d\n", msg, p->cutoff_min, p->cutoff_max);
  sprintf(msg, "%sValeurs de MRI : min=%ld max=%ld\n", msg, p->min_pixel, p->max_pixel);
  sprintf(msg, "%sCoefficient de valeur reelle : rcoeff=%f icomp=%f\n", msg, p->rcoeff, p->icomp);
  sprintf(msg, "%sValeurs de l'image : min=%f max=%f\n", msg, p->min_pixel*p->rcoeff, p->max_pixel*p->rcoeff);
  sprintf(msg, "%sResolution (mm) : dx=%f dy=%f dz=%f\n", msg, p->dx, p->dy, p->dz);
  sprintf(msg,"%s\n", msg);

  sprintf(msg, "%sPalette de couleurs : %s (no %d, id %d)\n", msg, cmap->szName, p->cmapinfo.noColormap, cmap->noMode);
  sprintf(msg, "%sMode d'affichage : %s\n", msg, imx_GetStrDisplayMode(p->cmapinfo.dspStyle));
  sprintf(msg, "%sType de palette : %s\n", msg, IMXColor_DescCMapType(cmap->noType));

  sprintf(msg,"%s\n", msg);

  sprintf(msg, "%sMasque : %s (%s)\n", msg, IMXPict_DescMaskType(p->mask_type), (p->mask_active ? "active" : "desactive"));
  sprintf(msg,"%s\n", msg);

  sprintf(msg, "%sAffichage : %s\n", msg, (p->bar ? "palette" : "-"));
  sprintf(msg,"%s\n", msg);

  sprintf(msg, "%sNom du fichier : %s\n", msg, p->filename);

  if(p->type)
    sprintf(msg, "%sType : %c\n", msg, p->type);

  sprintf(msg,"%s\n", msg);

  sprintf(msg, "%sPatient : %s (%s)\n", msg, p->patient.name, p->patient.d_birth);
  sprintf(msg, "%sMedecin : %s\n", msg, p->patient.medecin);
  sprintf(msg, "%sDate d'examen : %s\n", msg, p->patient.d_examen);
  sprintf(msg,"%s\n", msg);

  sprintf(msg, "%sTalairach : %s\n", msg, (p->tala.grid_on ? "active" : "desactive"));
  sprintf(msg, "%sId grille : %p\n", msg, p->tala.grid_id);

  if(p->tala.grid_on)
    {
      sprintf(msg,"%s\n", msg);
      sprintf(msg, "%sAC : (%f, %f, %f)\n", msg, p->tala.AC.x, p->tala.AC.y, p->tala.AC.z);
      sprintf(msg, "%sPC : (%f, %f, %f)\n", msg, p->tala.PC.x, p->tala.PC.y, p->tala.PC.z);
      sprintf(msg, "%sAP : (%f, %f, %f)\n", msg, p->tala.AP.x, p->tala.AP.y, p->tala.AP.z);
      sprintf(msg, "%sPP : (%f, %f, %f)\n", msg, p->tala.PP.x, p->tala.PP.y, p->tala.PP.z);
      sprintf(msg, "%sRM : (%f, %f, %f)\n", msg, p->tala.RM.x, p->tala.RM.y, p->tala.RM.z);
      sprintf(msg, "%sLM : (%f, %f, %f)\n", msg, p->tala.LM.x, p->tala.LM.y, p->tala.LM.z);
      sprintf(msg, "%sSP : (%f, %f, %f)\n", msg, p->tala.SP.x, p->tala.SP.y, p->tala.SP.z);
      sprintf(msg, "%sIP : (%f, %f, %f)\n", msg, p->tala.IP.x, p->tala.IP.y, p->tala.IP.z);
    }

  sprintf(msg,"%s\n", msg);

  sprintf(msg, "%siCoupe curviligne : %s\n", msg, (p->curv.state ? "active" : "desactive"));

  if(p->curv.state)
    {
      sprintf(msg, "%scentre de gravite : %d,%d,%d\n", msg, p->curv.xg, p->curv.yg, p->curv.zg);
      sprintf(msg, "%stheta : %f\n", msg, p->curv.theta);
      sprintf(msg, "%stheta : min=%f max=%f\n", msg, p->curv.theta_min, p->curv.theta_max);
      sprintf(msg, "%sphi : %f\n", msg, p->curv.phi);
      sprintf(msg, "%sphi : min=%f max=%f\n", msg, p->curv.phi_min, p->curv.phi_max);
      sprintf(msg, "%srayon : %f\n", msg, p->curv.rayon);
      sprintf(msg, "%szoom : %f\n", msg, p->curv.zoom);
    }

  sprintf(msg,"%s\n", msg);

  
  if(p->dti.is_dti)
    {
		  sprintf(msg, "%sInfo DTI\n", msg);
    	sprintf(msg, "%snombre de direction : %d\n", msg, p->dti.nb_directions);
      sprintf(msg, "%sgradient de diffusion : [ %f ; %f ; %f ]\n", msg, p->dti.gx,p->dti.gy,p->dti.gz);
      sprintf(msg, "%scoefficient de diffusion b : %f\n", msg, p->dti.b);
    }
  str = (GtkLabel*)gtk_label_new( msg );
  gtk_label_set_justify(str, GTK_JUSTIFY_LEFT);
  gtk_container_add( GTK_CONTAINER(window),GTK_WIDGET( str) );

  gtk_widget_show_all( window );

#endif /*__GTK_GUI*/
}

//
//  ---- IMXPict3d_CleanPicture() ------------------------
//
/*! \brief  Cleans a 3D picture  (developped for Therese =) */
//  Created by R. CHEVRIER on Dec 6th, 2001.
//
void IMXPict3d_CleanPicture( void )
{
  int pos3d = GET_PLACE3D( "Picture to reset" );

  IMXPict3d_ResetPicture( pos3d );

  show_picture_3d( pos3d );
}

//
//  ---- IMXPict3d_ResetPicture() ------------------------
//
/*! \brief  Resets 3D MRI and unactivates mask. */
/*! \param pos3d : numero de l'image*/
//  Created by R. CHEVRIER on Dec 6th, 2001.
//
void IMXPict3d_ResetPicture( int pos3d )
{
	#ifdef __GTK_GUI
	grphic3d *img = ptr_img_3d( pos3d );
	int i, j, k;
	
	for( i=0; i<img->width; i++ )
		for( j=0; j<img->height; j++ )
			for( k=0; k<img->depth; k++ )
				img->mri[i][j][k] = 0l;
	
	ptr_img(img->img_coro)->mask_active = FALSE;
	ptr_img(img->img_tran)->mask_active = FALSE;
	ptr_img(img->img_sagi)->mask_active = FALSE;
	
	img->mask_active = FALSE;
	img->zoom = 1;
	img->zoom_x = img->zoom_y = img->zoom_z = 0;
	#endif /*__GTK_GUI*/
}


/* --imx_brukermax_3d() ----------------------------------------------
 */
/*! \brief Calcul de icomp et rcoeff afin de tenir compte de la
**   modification de la taille du pixel
**
**   (8 bits, 16 bits, 32 bits)
**
**   \param amax : valeur maximum vrai de l'image
**   \param amin : valeur mimimun vrai de l'image
**   \param graphic : pointeur sur la structure image 3D
**
**    ecriture dans graphic de icomp, rcoeff, max_pixel, min_pixel
**	\retval 1
***********************************************************************/
int      imx_brukermax_3d(float amax, float amin, grphic3d *graphic)
{
  double aicomp,temp1,two;
  long maxpix,minpix;
  float max;
  double eps=1e-5;

  max=(float) MAXI(amax,fabs((double) amin));

  maxpix=MAXMRI3D;

  two=2.;
  if(max==0.0)
    aicomp=0.;
  else
    {
      temp1=max/maxpix;
      aicomp=log(temp1)/log(two);
    }

  if(aicomp<=0.0)
    {
      //     graphic->icomp=(float)((int)floor(aicomp));
      //difference significative on prend la valeur au dessus
      //sinon on prend l'arrondi
      if (fabs(aicomp-floor(aicomp+0.5))>eps) graphic->icomp=(float)ceil(aicomp);
      else graphic->icomp=(float)floor(aicomp+0.5);
    }
  else
    {
      //     graphic->icomp=(float)((int)floor((aicomp+1.0)));
      //difference significative on prend la valeur au dessus
      //sinon on prend l'arrondi
      if (fabs(aicomp-floor(aicomp+0.5))>eps) graphic->icomp=(float)ceil(aicomp+1.0);
      else graphic->icomp=(float)floor(aicomp+1.0+0.5);
    }
   
  aicomp=graphic->icomp;
  graphic->rcoeff=(float)pow(two,aicomp);
  maxpix=(long)floor(amax/graphic->rcoeff);
  graphic->max_pixel=maxpix;
  graphic->cutoff_max=maxpix;
  minpix=(long)floor(amin/graphic->rcoeff);
  graphic->min_pixel=minpix;
  graphic->cutoff_min=minpix;
#ifdef __GTK_GUI
  switch (graphic->cmapinfo.dspStyle)
    {
    case DSP_POSITIVE : graphic->cutoff_min=0; break;
    case DSP_NEGATIVE : graphic->cutoff_max=0; break;
    default 			: break;
    }
#endif /*__GTK_GUI*/


  if(MESG3D_DEBUG) printf("maxpix=%ld max=%f icomp=%f rcoef=%f \n"
			  ,graphic->max_pixel,max,graphic->icomp,graphic->rcoeff);

  return(1);

}

/* --imx_brukermax_update_3d_p() ----------------------------------------------
 */
/*! \brief Met a jour l'image pour avoir une dynamique de calcul maximum
**
**   (8 bits, 16 bits, 32 bits)
**
**   \param im1 : pointeur sur la structure image 3D a mettre a jour
**
***********************************************************************/
void  imx_brukermax_update_3d(grphic3d *im)
{
 float rcoeff, rapport;
 int i,j,k;
 
 rcoeff=im->rcoeff;
 imx_inimaxminpixel_3d_p(im);
 imx_brukermax_3d(im->max_pixel*rcoeff,im->min_pixel*rcoeff,im);

rapport=rcoeff/im->rcoeff;

 for(i=0;i<im->width;i++)
 for(j=0;j<im->height;j++)
 for(k=0;k<im->depth;k++)
 	{
	im->mri[i][j][k]=(int)(im->mri[i][j][k]*rapport);
	}


}

/***************************************************************
 **  imx_iniimg3d()
 */
/*!
**   \brief Initialise les variable img_tran, img_sagi, img_coro
**   pour la visualisation multi-planar 3D
**   \param  pos : position de l'image 3d
**   \param  img_coro : pointeur sur la position de l'image coro correspondante (E/S)
**   \param  img_sagi : pointeur sur la position de l'image sagi correspondante (E/S)
**   \param  img_tran : pointeur sur la position de l'image tran correspondante (E/S)
**   \param  img_curv : pointeur sur la position de l'image curv correspondante (E/S)
*****************************************************************/
void imx_iniimg3d(int pos, int *img_coro, int *img_sagi, int *img_tran, int *img_curv)
{
  int nb_3d_ligne,deb_ligne;


  if (pos >100 && pos <128)
    { // images temporaire (dans les automates par ex.)
      return ;
    }
  else
    { // images a afficher dans l'interface
      nb_3d_ligne = MINI(MAX_PIC_IN_LINE/2,MAX_PIC_3D_IN_LINE);
      deb_ligne = ((int) ((abs(pos)-1)/nb_3d_ligne))*2*MAX_PIC_IN_LINE;
      *img_coro = ((abs(pos)-1)%nb_3d_ligne)*2+1+deb_ligne;
      *img_sagi = ((abs(pos)-1)%nb_3d_ligne)*2+2+deb_ligne;
      *img_tran = ((abs(pos)-1)%nb_3d_ligne)*2+1+deb_ligne+MAX_PIC_IN_LINE;
      *img_curv = *img_tran+1;
      *img_coro = *img_coro*(SIGN(pos));
      *img_sagi = *img_sagi*(SIGN(pos));
      *img_tran = *img_tran*(SIGN(pos));
      *img_curv = *img_curv*(SIGN(pos));
      if (*img_coro > MAX_PICTURES || *img_sagi >  MAX_PICTURES || *img_tran > MAX_PICTURES || *img_curv > MAX_PICTURES)
	printf ( "#####  ATTENTION Pas assez d'images 2D pour la 3D \n");
    }
}

/*******************************************
 ** --  imx_inimaxpixel_3d() ---------------------
 **  Sous programme
 **   -  Cheche maxpixel dans l'image
 **
 ********************************************/
int	imx_inimaxpixel_3d(int im_1)
{
  grphic3d *im1;

  /*
   *       Recherche du maximum vrai
   */
  im1=ptr_img_3d(im_1);
  imx_inimaxpixel_3d_p(im1);

  return(1);
}
/*******************************************
 ** --  imx_inimaxpixel_3d_p() ---------------------
 */
/*!  Sous programme
**   -  Cheche maxpixel dans l'image et positionne max_pixel et cutoff_max de l'image
**	\param im1 : image 3d (E/S)
**  \retval 1
********************************************/
int	imx_inimaxpixel_3d_p(grphic3d *im1)
{
  long maxpixel;
  int i,j,k;
  int wdth,hght,dpth;

  /*
   *       Recherche du maximum vrai
   */

  wdth=im1->width;
  hght=im1->height;
  dpth=im1->depth;

  maxpixel=0;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++)
	maxpixel=MAXI(maxpixel,im1->mri[i][j][k]);

  im1->max_pixel=maxpixel;
  im1->cutoff_max=maxpixel;

  return(1);
}

/*******************************************
 ** --  imx_inimaxminpixel_3d() ---------------------*/
/*!
**  Sous programme
**   -  Cherche maxpixel et minpixel dans l'image
**	\param im_1 : numero de l'image a modifier
**	\retval 1
********************************************/
int	imx_inimaxminpixel_3d(int im_1)
{
  grphic3d *im1;

  /*
   *       Recherche du maximum vrai
   */
  im1=ptr_img_3d(im_1);
  imx_inimaxminpixel_3d_p(im1);

  return(1);
}

/*******************************************
 ** --  imx_inimaxminpixel_3d_p() ---------------------
 */
/*!  Sous programme
**   -  Cheche maxpixel et minpixel  dans l'image
**  \param im1 : l'image concernee
**  \retval : 1
**  \remark im1 est modifie (positionnement de minpixel et maxpixel)
********************************************/
int	imx_inimaxminpixel_3d_p(grphic3d *im1)
{
  long maxpixel,minpixel;
  int i,j,k;
  int wdth,hght,dpth;

  /*
   *       Recherche du maximum vrai
   */

  wdth=im1->width;
  hght=im1->height;
  dpth=im1->depth;

  maxpixel=0;
  minpixel=LONGMAX;
  for(i=0;i<wdth;i++)
    for(j=0;j<hght;j++)
      for(k=0;k<dpth;k++) {
	if (maxpixel < im1->mri[i][j][k] )
	  /*			  printf ("i=%d j=%d k=%d \n",i,j,k);*/
	  maxpixel=MAXI(maxpixel,im1->mri[i][j][k]);
	minpixel=MINI(minpixel,im1->mri[i][j][k]);
      }

  im1->max_pixel=maxpixel;
  im1->min_pixel=minpixel;
  im1->cutoff_max=maxpixel;
  im1->cutoff_min=minpixel;

  return(1);
}

/*******************************************
 ** --  imx_iniminpixel_3d() ---------------------
 **  Sous programme
 **   -  Cheche minpixel dans l'image
 **
 ********************************************/
int	imx_iniminpixel_3d(int im_1)
{
  grphic3d *im1;

  /*
   *       Recherche du minimum vrai
   */
  im1=ptr_img_3d(im_1);
  imx_iniminpixel_3d_p(im1);

  return(1);
}

/*******************************************
 ** --  imx_iniminpixel_3d_p() ---------------------
 **  Sous programme
 **   -  Cheche minpixel dans l'image
 **
 ********************************************/
int	imx_iniminpixel_3d_p(grphic3d *im1)
{
  long minpixel;
  int i,j,k;
  int w,h,d;


  /*
   *       Recherche du minimum vrai
   */

  minpixel=0;
  w = im1->width;
  h = im1->height;
  d = im1->depth;
  for(i=0;i<w;i++)
    for(j=0;j<h;j++)
      for(k=0;k<d;k++)
	minpixel=MINI(minpixel,im1->mri[i][j][k]);

  im1->min_pixel=minpixel;

  return(1);
}


/*******************************************
 ** --  imx_iniparam_3d() ---------------------
 **  Sous programme
 **   -  Initialisation des parametres
 **      max-pixel, min_pixel, rcoeff, icomp
 **
 ********************************************/
int	imx_iniparam_3d(int im_1, long int max_pixel, long int min_pixel, float icomp, float rcoeff)
{
  grphic3d *im1;

  im1=ptr_img_3d(im_1);

  im1->max_pixel=max_pixel;
  im1->min_pixel=min_pixel;
  im1->cutoff_max=max_pixel;
  im1->cutoff_min=min_pixel;
  im1->icomp=icomp;
  im1->rcoeff=rcoeff;

  return(1);
}

/*******************************************
 ** --  imx_iniparaimg_3d() ---------------------
 */
/*!  Sous programme
**   -  Initialisation des parametres
**      max-pixel, min_pixel, rcoeff, icomp
**    avec recherche de max_pixel et min_pixel
**	\param im_1 : numero de l'image a traiter
**  \retval 1
********************************************/
int	imx_iniparaimg_3d(int im_1)
{
  grphic3d *im1;

  im1=ptr_img_3d(im_1);
  imx_iniparaimg_3d_p(im1);

  return(1);
}

/*******************************************
 ** --  imx_iniparaimg_3d_p() ---------------------
 **  Sous programme
 **   -  Initialisation des parametres
 **      max-pixel, min_pixel, rcoeff, icomp
 **    avec recherche de max_pixel et min_pixel
 **
 ********************************************/
int	imx_iniparaimg_3d_p(grphic3d *im1)
{

  imx_inimaxpixel_3d_p(im1);
  im1->min_pixel=0;
  im1->icomp=0;
  im1->rcoeff=1;

  return(1);
}

/*******************************************
 **   imx_copie_param_3d()
 */
/*! \brief Copie des parametres de im1 dans im2 (3D)
**  \param im_1 : numero de l'image image source
**	\param im_2 : numero de l'image destination
********************************************/
int        imx_copie_param_3d(int im_1, int im_2)
{
  grphic3d *im1,*im2;
  int img_tran,img_sagi,img_coro,img_curv;

  im1=ptr_img_3d(im_1);
  im2=ptr_img_3d(im_2);

  if (im_2==0) {
    im2->img_tran=0;
    im2->img_sagi=0;
    im2->img_coro=0;
    im2->img_curv=0;
  }
  else {
    imx_iniimg3d(im_2,&img_coro,&img_sagi,&img_tran,&img_curv);
    im2->img_tran=img_tran;
    im2->img_sagi=img_sagi;
    im2->img_coro=img_coro;
    im2->img_curv=img_curv;
  }

  imx_copie_param_3d_p(im1,im2);
#ifdef __GTK_GUI
  IMXColor_BindColormap_3d(im_2, im1->cmapinfo.noColormap);
#endif /*__GTK_GUI*/

  return(1);
}

/*******************************************
 ** --  imx_copie_param_3d_p() ----------------------------
 */
/*! \brief Copie des parametres lies au patient de im1 dans im2 (3D)
**	\param im1 : ptr sur l'image source
**	\param im2 : ptr sur l'image destination
**	\retval 1, im2 est modifiee
**  \remark im2 doit avoir ete alloue precedement
********************************************/
int imx_copie_patient_params_3d_p(grphic3d *im1, grphic3d *im2)
{
  strcpy(im2->filename, im1->filename);
  strcpy(im2->patient.name, im1->patient.name);
  strcpy(im2->patient.d_birth, im1->patient.d_birth);
  strcpy(im2->patient.medecin, im1->patient.medecin);
  strcpy(im2->patient.d_examen, im1->patient.d_examen);

  return(1);
}

/*******************************************
 ** --  imx_copie_param_3d_p() ----------------------------
 */
/*! \brief Copie des parametres concernant la visualisation de im1 dans im2 (3D)
**	\param im1 : ptr sur l'image source
**	\param im2 : ptr sur l'image destination
**	\retval 1, im2 est modifiee
**  \remark im2 doit avoir ete alloue precedement
********************************************/
int imx_copie_visual_params_3d_p(grphic3d *im1, grphic3d *im2)
{

  im2->zoom            = im1->zoom;
  im2->zoom_x          = im1->zoom_x;
  im2->zoom_y          = im1->zoom_y;
  im2->zoom_z          = im1->zoom_z;
  im2->x3dc            = im1->x3dc;
  im2->y3dc            = im1->y3dc;
  im2->z3dc            = im1->z3dc;
#ifdef __GTK_GUI
  im2->cmapinfo.dspStyle  =   im1->cmapinfo.dspStyle    ; 	  // colormap curve style to use
  im2->cmapinfo.cmap  	  =   im1->cmapinfo.cmap  	  ;
  im2->bar	              =   im1->bar	  ;
  IMXColor_BindColormap_3d_ptr(im2, im1->cmapinfo.noColormap);
  //modif BERST : mise a jour des cmap pour chaque vue 2D de l'image 3D
  // ne serait il pas judicieux d'inclure la partie suivante dans IMXColor_BindColormap_3d_ptr????
  if (im2->img_tran)
    IMXColor_BindColormap(im2->img_tran, im1->cmapinfo.noColormap);
  if (im2->img_sagi)
    IMXColor_BindColormap(im2->img_sagi, im1->cmapinfo.noColormap);
  if (im2->img_coro)
    IMXColor_BindColormap(im2->img_coro, im1->cmapinfo.noColormap);
  if (im2->img_curv)
    IMXColor_BindColormap(im2->img_curv, im1->cmapinfo.noColormap);
  // fin modif BERST
#endif /*__GTK_GUI*/
  return(1);
}

/*******************************************
 ** --  imx_copie_apparance_3d_p() ----------------------------
 */
/*! \brief Copie des parametres concernant la visualisation de im1 dans im2 (3D)
**	\param im1 : ptr sur l'image source
**	\param im2 : ptr sur l'image destination
**	\retval 1, im2 est modifiee
**  \remark im2 doit avoir ete alloue precedement
********************************************/
int imx_copie_apparance_params_3d_p(grphic3d *im1, grphic3d *im2)
{
  im2->zoom            = im1->zoom;
  im2->zoom_x          = im1->zoom_x;
  im2->zoom_y          = im1->zoom_y;
  im2->zoom_z          = im1->zoom_z;
  im2->x3dc            = im1->x3dc;
  im2->y3dc            = im1->y3dc;
  im2->z3dc            = im1->z3dc;
  return 1;
}

/*******************************************
 ** --  imx_copie_image_params_3d_p() ----------------------------
 */
/*! \brief Copie des parametres concernant l'image
**	\param im1 : ptr sur l'image source
**	\param im2 : ptr sur l'image destination
**	\retval 1, im2 est modifiee
**  \remark im2 doit avoir ete alloue precedement
********************************************/
int  imx_copie_image_params_3d_p(grphic3d *im1, grphic3d *im2)
{
  im2->min_pixel       = im1->min_pixel;
  im2->max_pixel       = im1->max_pixel;
  im2->cutoff_min      = im1->cutoff_min;
  im2->cutoff_max      = im1->cutoff_max;
  im2->icomp           = im1->icomp;
  im2->rcoeff          = im1->rcoeff;
  im2->img_type        = IMAGE_NORMAL;
  im2->mask_type       = im1->mask_type;
  im2->type            = 'b';
  im2->bitppixel       = im1->bitppixel;

  return(1);
}

/*******************************************
 ** --  imx_copie_dimension_params_3d_p() ----------------------------
 */
/*! \brief Copie des parametres de dimension l'image
**	\param im1 : ptr sur l'image source
**	\param im2 : ptr sur l'image destination
**	\retval 1, im2 est modifiee
**  \remark im2 doit avoir ete alloue precedement
********************************************/
int  imx_copie_dimension_params_3d_p(grphic3d *im1, grphic3d *im2)
{
  im2->dx              = im1->dx;
  im2->dy              = im1->dy;
  im2->dz              = im1->dz;
  im2->width           = im1->width;
  im2->height          = im1->height;
  im2->depth           = im1->depth;
  im2->x3dc            = im1->x3dc;
  im2->y3dc            = im1->y3dc;
  im2->z3dc            = im1->z3dc;


  return(1);
}

/*******************************************
 ** --  imx_copie_dti_params_3d_p() ----------------------------
 */
/*! \brief Copie des parametres relatif au dti de l'image
**	\param im1 : ptr sur l'image source
**	\param im2 : ptr sur l'image destination
**	\retval 1, im2 est modifiee
**  \remark im2 doit avoir ete alloue precedement
********************************************/
int  imx_copie_dti_params_3d_p(grphic3d *im1, grphic3d *im2)
{
  im2->dti.is_dti              	= im1->dti.is_dti;
	im2->dti.nb_directions      	= im1->dti.nb_directions;
	im2->dti.gx              			= im1->dti.gx;
	im2->dti.gy              			= im1->dti.gy;
	im2->dti.gz              			= im1->dti.gz;
	im2->dti.b              			= im1->dti.b;

  return(1);
}

/*******************************************
 ** --  imx_copie_param_3d_p() ----------------------------
 */
/*! \brief Copie des parametres de im1 dans im2 (3D)
**	\param im1 : ptr sur l'image source
**	\param im2 : ptr sur l'image destination
**	\retval 1, im2 est modifiee
**  \remark im2 doit avoir ete alloue precedement
********************************************/
int imx_copie_param_3d_p(grphic3d *im1, grphic3d *im2)
{
  imx_copie_image_params_3d_p(im1, im2);
  imx_copie_dimension_params_3d_p(im1, im2);
  imx_copie_visual_params_3d_p(im1, im2);
  imx_copie_patient_params_3d_p(im1, im2);
  imx_copie_dti_params_3d_p(im1, im2);

  return(1);
}

/* -- SetImg3d() -----------------------------------------------------
**
**   Initialise les variable img_tran, img_sagi, img_coro
**   pour la visualisation multi-planar 3D
**------------------------------------------------------------*/
void SetImg3d(int pos, grphic3d *graphic)
{
  int img_tran,img_sagi,img_coro,img_curv;

  /* Positionnement de img_tran, img_sagi,img_coro */
  if (pos==0) { /*  Initialisation roi 3D */
		//PUT_WARN("Error SetImg3d : pos=0\n");
    return ;
  }
  else {
    imx_iniimg3d(pos,&img_coro,&img_sagi,&img_tran,&img_curv);
    graphic->img_coro=img_coro;
    graphic->img_sagi=img_sagi;
    graphic->img_tran=img_tran;
    graphic->img_curv=img_curv;
    graphic->width=MAX_WIDTH_3D;
    graphic->height=MAX_HEIGHT_3D;
    graphic->depth=MAX_DEPTH_3D;
    graphic->zoom=1;
    graphic->zoom_x=0;
    graphic->zoom_y=0;
    graphic->zoom_z=0;
    graphic->pos=pos;
#ifdef __GTK_GUI
    graphic->cmapinfo.noColormap = IMXColor_NewColormap(7, 0, NULL, 0); /* Palette grise par defaut */
    list_moveto(&_lstColormaps, graphic->cmapinfo.noColormap );
    graphic->cmapinfo.cmap = (ptrColormap) list_item(&_lstColormaps);
#endif /*__GTK_GUI*/
  }

}

/* -- imx_pos_2dto3d() -----------------------------------------------------
**
**   Calcul de la position 3d en fonction de la position 2d
**
**------------------------------------------------------------*/
int  imx_pos_2dto3d(int pos)
{
  int nbli2d,nbco2d,nbli3d,nbco3d,pos3d;
  int sign = SIGN(pos);
  int abs_pos = abs(pos);
  
  nbli2d=(int)((abs_pos-1)/MAX_PIC_IN_LINE)+1;
  nbco2d=(int)((abs_pos-1)%MAX_PIC_IN_LINE)+1;
  nbli3d=(nbli2d-1)/2+1;
  nbco3d=(nbco2d-1)/2+1;
  pos3d=(int)((nbli3d-1)*MAX_PIC_3D_IN_LINE)+nbco3d;
  pos3d=pos3d*sign;
  return(pos3d);
}

/* -- imx_2dto3d() -----------------------------------------------------
**
**   Calcul de la position x,y,z 3d en fonction de la position x,y 2d
**   Return : pos3d
**
**------------------------------------------------------------*/
int  imx_ecran2dto3d(int pos2d, int x2d, int y2d, int *x3d, int *y3d, int *z3d)
{
	int pos3d,im_t,im_s,im_c,im_curv;
	grphic3d *img;
	int w,h,d;
	pos3d=imx_pos_2dto3d(pos2d);
	
	img=ptr_img_3d(pos3d);
	im_c=img->img_coro;
	im_s=img->img_sagi;
	im_t=img->img_tran;
	im_curv=img->img_curv;
	
	if (pos2d==im_c) {
		*x3d=(x2d)/img->zoom+img->zoom_x;
		*y3d=(y2d)/img->zoom+img->zoom_y;
		*z3d=img->z3dc;
	}
	if (pos2d==im_s) {
		*x3d=img->x3dc;
		*y3d=(y2d)/img->zoom+img->zoom_y;
		*z3d=(x2d)/img->zoom+img->zoom_z;
	}
	if (pos2d==im_t) {
		*x3d=(x2d)/img->zoom+img->zoom_x;
		*y3d=img->y3dc;
		*z3d=(y2d)/img->zoom+img->zoom_z;
	}
	if (pos2d==im_curv) {
		*x3d=(x2d)/img->zoom;
		*y3d=(y2d)/img->zoom;
		*z3d=img->z3dc;;
	}
	
	w = img->width;
	h = img->height;
	d = img->depth;
	
	(*x3d)=CLAMP(*x3d,0,w-1);
	(*y3d)=CLAMP(*y3d,0,h-1);
	(*z3d)=CLAMP(*z3d,0,d-1);
	
	return(pos3d);
}

/* -- imx_3dto2d() -----------------------------------------------------
**
**   Calcul de la position x,y en 2d en fonction de la position x,y z et de l'orientation
**   Return : pos3d
**
**------------------------------------------------------------*/
int  imx_ecran3dto2d(int pos2d, int *x2d, int *y2d, int x3d, int y3d, int z3d)
{
	int pos3d,im_t,im_s,im_c,im_curv;
	grphic3d *img;
	
	pos3d=imx_pos_2dto3d(pos2d);
	
	img=ptr_img_3d(pos3d);
	im_c=img->img_coro;
	im_s=img->img_sagi;
	im_t=img->img_tran;
	im_curv = img->img_curv;
	
	if (pos2d==im_c) {
		*x2d=(x3d)/img->zoom+img->zoom_x;
		*y2d=(y3d)/img->zoom+img->zoom_y;
	}
	if (pos2d==im_s) {
		*y2d=(y3d)/img->zoom+img->zoom_y;
		*x2d=(z3d)/img->zoom+img->zoom_z;
	}
	if (pos2d==im_t) {
		*x2d=(x3d)/img->zoom+img->zoom_x;
		*y2d=(z3d)/img->zoom+img->zoom_z;
	}
	if (pos2d==im_curv) {
		*x2d=(x3d)/img->zoom;
		*y2d=(y3d)/img->zoom;
	}
	
	return(pos3d);
}

int resize_grphic3d_buffer(grphic3d *img, TDimension new_width, TDimension new_height, TDimension new_depth)
{
 //le buffer est deja a la bonne taille: presque rien a faire
 if ((img->mri)&&(img->mri_width>=new_width)&&(img->mri_height>=new_height)&&(img->mri_depth>=new_depth))
 {
  img->width=new_width;
  img->height=new_height;
  img->depth=new_depth;
  return 0;
 }

 //libere l'ancien buffer
 if (img->mri) free_mri_3d(img->mri);

 //aloue le nouveau buffer
 img->mri=cr_mri_3d(new_width, new_height, new_depth);
 if (!img->mri) { fprintf (stderr, "memory allocation failed in resize_grphic3d_buffer\n"); return 1; }
 //met a jour les dimensions de l'image
 img->mri_width=new_width;
 img->mri_height=new_height;
 img->mri_depth=new_depth;
 img->width=new_width;
 img->height=new_height;
 img->depth=new_depth;

 return 0;
}

int copie_image_param_adapt_buffer(grphic3d *imdeb, grphic3d *imres)
{
 int err=0;

 imx_copie_param_3d_p(imdeb,imres);
 err=resize_grphic3d_buffer(imres, imdeb->width, imdeb->height, imdeb->depth);

 return err;
}
