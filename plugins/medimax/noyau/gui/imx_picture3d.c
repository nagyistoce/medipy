/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*!
//  \file    : imx_picture3d.c
//  Project : GTK-Imagix
//
//  \brief Description: gui fonctions interface
//      ...
//
//  Copyright (c) 1993, ULP-IPB Strasbourg
//  All rights are reserved
//
//  >>  Created on ... ..th, .... by F. BOULEAU
//  >>  Last user action by Remy CHEVRIER on  Aug 29th, 2001.
//
*/

#include <config.h>
#include <string.h>

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_3d.h"
#include "noyau/gui/imx_picture3d.h"
#include "noyau/zone_3d.h"
#include "noyau/imx_lang.h"
#ifdef __GTK_GUI
#include "noyau/imx_zone.h"
#endif // __GTK_GUI
extern int IMXColor_BindColormap_3d_ptr(ptr_grphic3d ptrPic, int nNewCMap);

/******************************************************************
 **   show_picture_3d()
 */
/*!
** \brief  This function show a 3d picture with 3 images.
**
**
**  \param  pos:  index of picture to show
**
**  \retval -1 if an error occurred, pos otherwise
**
*/
int show_picture_3d(int pos)
{
#ifdef  __GTK_GUI
  grphic3d *img3d;

  if(pos > MAX_PICTURES)
    {
    _imagix_error   = IMX_BAD_POSITION;
    _imagix_err_ext = pos;
    return(-1);
    }

  if(MESG_DEBUG)
    printf("Show> %d\n", pos);

  if (pos<0)
    pos = abs(pos);

  if(pos <= MAX_PICTURES_3D)
    {
      if(MESG_DEBUG)
        printf("Show picture at %d\n", pos);

      img3d = ptr_img_3d(pos);
      picture_3dto2d(img3d);

      if(img3d->mask_active)
        {
        img3d->mask->zoom   = img3d->zoom;
        img3d->mask->zoom_x = img3d->zoom_x;
        img3d->mask->zoom_y = img3d->zoom_y;
        img3d->mask->zoom_z = img3d->zoom_z;
        if(MESG_DEBUG)
          fprintf(stderr, "picture mask: %d\n", pos);
        picture_mask_3dto2d(img3d);
        ptr_img(img3d->mask->img_curv)->bar = img3d->mask->bar;
        ptr_img(img3d->mask->img_tran)->indexes
           = (img3d->mask->bar ? TEXT_ALIGN_RIGHT : TEXT_NONE);
        }

      ptr_img(img3d->img_curv)->bar = img3d->bar;
      ptr_img(img3d->img_tran)->indexes
         = (img3d->bar ? TEXT_ALIGN_RIGHT : TEXT_NONE);

      show_picture(img3d->img_coro);
      show_picture(img3d->img_sagi);
      show_picture(img3d->img_tran);
      show_picture(img3d->img_curv);
    }
#endif /* GTK_GUI*/

  return(pos);
}

/******************************************************************************
 ** -- Refresh_3d() ------------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July, 23rd 1999
 **
 ******************************************************************************/

void Refresh_3d(void)
{
#ifdef __GTK_GUI
  /* Redraw all pictures... */
  int i;
  for(i=1;i<=MAX_PICTURES_3D;i++)
    show_picture_3d(i);
  if(_imx_roi3d.pos>0&& _imx_roi3d.pos<=MAX_PICTURES_3D)
    show_roi_3d(_imx_roi3d.pos);
#endif /*__GTK_GUI*/
}

/*
**  -- picture_3dto2d() -------------------------------------------------------
**
**  This function transform a 3d picture to 3 images 2d.
**
**  Arguments:
**  - img : Which picture ?
**
**               ----------> i (x)         ----------> k (z)
**              |            |
**              | img_coro            | img_sagi
**              |            |
**              v            v
**
**              j            j
**             (y)               (y)
**
**               ----------> i (x)
**              |
**              | img_tran
**              |
**              v
**
**              k
**             (z)
**
*/

int picture_3dto2d(grphic3d *img)
{
#ifdef  __GTK_GUI
  grphic *imt,*ims,*imc,*imcurv;

  imc=ptr_img(img->img_coro);
  ims=ptr_img(img->img_sagi);
  imt=ptr_img(img->img_tran);
  imcurv=ptr_img(img->img_curv);

  picture_3dto2d_param(img,imt,ims,imc,imcurv);
  picture_3dto2d_copy(img,imt,ims,imc,imcurv);
#endif  //__GTK_GUI

  return(1);
}

int picture_3dto2d_param(grphic3d *img, grphic *imt, grphic *ims, grphic *imc, grphic *imcurv)
{
#ifdef  __GTK_GUI
  if (imc)
    {
      imc->cmapinfo = img->cmapinfo;
      imc->alpha_channel = img->alpha_channel;
      imc->img_type = img->img_type;
      imc->mask_type = img->mask_type;

      imc->max_pixel=img->max_pixel;
      imc->min_pixel=img->min_pixel;
      imc->cutoff_max=img->cutoff_max;
      imc->cutoff_min=img->cutoff_min;
      imc->icomp=img->icomp;
      imc->rcoeff=img->rcoeff;
      imc->height=img->height;
      imc->width=img->width;
      imc->zoom=img->zoom;
      imc->zoom_x=img->zoom_x*img->zoom;
      imc->zoom_y=img->zoom_y*img->zoom;
      imc->dx=img->dx;
      imc->dy=img->dy;
      strcpy(imc->filename,img->filename);
      strcpy(imc->patient.name,img->patient.name);
      strcpy(imc->patient.d_birth,img->patient.d_birth);
      strcpy(imc->patient.medecin,img->patient.medecin);
      strcpy(imc->patient.d_examen,img->patient.d_examen);
    }

  if (ims)
    {
      ims->cmapinfo = img->cmapinfo;
      ims->alpha_channel = img->alpha_channel;
      ims->img_type = img->img_type;
      ims->mask_type = img->mask_type;
      ims->max_pixel=img->max_pixel;
      ims->min_pixel=img->min_pixel;
      ims->cutoff_max=img->cutoff_max;
      ims->cutoff_min=img->cutoff_min;
      ims->icomp=img->icomp;
      ims->rcoeff=img->rcoeff;
      ims->height=img->height;
      ims->width=img->depth;
      ims->zoom=img->zoom;
      ims->zoom_x=img->zoom_z*img->zoom;
      ims->zoom_y=img->zoom_y*img->zoom;
      ims->dx=img->dz;
      ims->dy=img->dy;
      strcpy(ims->filename,img->filename);
      strcpy(ims->patient.name,img->patient.name);
      strcpy(ims->patient.d_birth,img->patient.d_birth);
      strcpy(ims->patient.medecin,img->patient.medecin);
      strcpy(ims->patient.d_examen,img->patient.d_examen);
    }

  if (imt)
    {
      imt->cmapinfo = img->cmapinfo;
      imt->alpha_channel = img->alpha_channel;
      imt->img_type = img->img_type;
      imt->mask_type = img->mask_type;
      imt->max_pixel=img->max_pixel;
      imt->min_pixel=img->min_pixel;
      imt->cutoff_max=img->cutoff_max;
      imt->cutoff_min=img->cutoff_min;
      imt->icomp=img->icomp;
      imt->rcoeff=img->rcoeff;
      imt->height=img->depth;
      imt->width=img->width;
      imt->zoom=img->zoom;
      imt->zoom_x=img->zoom_x*img->zoom;
      imt->zoom_y=img->zoom_z*img->zoom;
      imt->dx=img->dx;
      imt->dy=img->dz;
      strcpy(imt->filename,img->filename);
      strcpy(imt->patient.name,img->patient.name);
      strcpy(imt->patient.d_birth,img->patient.d_birth);
      strcpy(imt->patient.medecin,img->patient.medecin);
      strcpy(imt->patient.d_examen,img->patient.d_examen);
    }

  if (imcurv)
    {
      // partial copy ... it's enough
      imcurv->cutoff_min = img->cutoff_min;
      imcurv->cutoff_max = img->cutoff_max;
      imcurv->max_pixel  = img->max_pixel;
      imcurv->min_pixel  = img->min_pixel;
      imcurv->rcoeff     = img->rcoeff;
      imcurv->cmapinfo   = img->cmapinfo;
      imcurv->mask_type  = img->mask_type;
      imcurv->alpha_channel =img->alpha_channel;
    }

#endif  //  __GTK_GUI
  return(1);
}

int picture_3dto2d_copy(grphic3d *img, grphic *imt, grphic *ims, grphic *imc, grphic *imcurv)
{
#ifdef  __GTK_GUI
  int i,j,k;
  int wdth,hght,dpth;
  int debi,debj,debk;
  int x3dc,y3dc,z3dc;

  wdth=img->width;
  hght=img->height;
  dpth=img->depth;

  x3dc=img->x3dc;
  y3dc=img->y3dc;
  z3dc=img->z3dc;

  debi=debj=debk=0;

  if(z3dc<img->depth)
  {
	  for (i=debi;i<wdth;i++)
	  for (j=debj;j<hght;j++)
      imc->mri[i][j]=img->mri[i][j][z3dc];
  }
  else
  {
	  printf("WARNING: picture_3dto2d_copy: cursor coordinates are outside of image!\n");
  }

  if(x3dc<img->width)
  {
	  for (j=debj;j<hght;j++)
	  for (k=debk;k<dpth;k++)
      ims->mri[k][j]=img->mri[x3dc][j][k];
  }
  else
  {
	  printf("WARNING: picture_3dto2d_copy: cursor coordinates are outside of image!\n");
  }

  if(y3dc<img->height)
  {
	  for (i=debi;i<wdth;i++)
	  for (k=debk;k<dpth;k++)
      imt->mri[i][k]=img->mri[i][y3dc][k];
  }
  else
  {
	  printf("WARNING: picture_3dto2d_copy: cursor coordinates are outside of image!\n");
  }

#endif  //  __GTK_GUI
  return(1);
}


int picture_mask_3dto2d(grphic3d *img)
{
#ifdef  __GTK_GUI
  grphic *imt,*ims,*imc,*imcurv;
  grphic3d *mask;

  mask=ptr_mask_3d_p(img);
  imc=ptr_mask(img->img_coro);
  ims=ptr_mask(img->img_sagi);
  imt=ptr_mask(img->img_tran);
  imcurv=ptr_mask(img->img_curv);
  mask->x3dc = img->x3dc;
  mask->y3dc = img->y3dc;
  mask->z3dc = img->z3dc;

  picture_3dto2d_param(mask, imt, ims, imc, imcurv);
  picture_3dto2d_copy(mask, imt, ims, imc, imcurv);
#endif  //  __GTK_GUI

  return(1);
}





/******************************************************************************
 ** -- imx_reinitialise_visual_params_3d_p() ------------------------------------------------------------
 **
 **
 **
 ******************************************************************************/
int imx_reinitialise_visual_params_3d_p(grphic3d *imres)
{
#ifdef  __GTK_GUI
  int max_dimension;
  int max_visualisation;

  if (!imres)
    return 1;

  imres->zoom_x = 0;
  imres->zoom_y = 0;
  imres->zoom_z = 0;

  max_dimension = MAXI(imres->height,imres->width);
  max_dimension = MAXI(max_dimension,imres->depth);
  max_visualisation = MAXI(DEFAULT_IMAGIX_HEIGHT,DEFAULT_IMAGIX_WIDTH);
  printf("%d   %d\n",max_dimension,max_visualisation);
  imres->zoom = max_visualisation/max_dimension;
  imres->zoom = MAXI(1,imres->zoom);
  imres->x3dc = imres->width/2;
  imres->y3dc = imres->height/2;
  imres->z3dc = imres->depth/2;

#endif /*  __GTK_GUI*/
  return 0;
}



#ifdef  __GTK_GUI
/******************************************************************************
 ** -- Zoom_Vi_3d() ------------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July, 23rd 1999
 **
 ******************************************************************************/
void Zoom_Vi_3d(void)
{ /* Zoom for visualisation... */

  char *quest[16];
  int pos, i;
  long nzoom;
  //			long wdth;  unused variable

  pos=GET_WND3D(NULL);
  for(i=0;i<16;i++)
    quest[i]=CALLOC(80,char);

  strcpy(quest[0],"No reduction");
  strcpy(quest[1],"Zoom 2");
  strcpy(quest[2],"Zoom 4");
  strcpy(quest[3],"\0");
  nzoom=GETV_QCM(
		 TEXT0043,(char **)quest);
  nzoom = (nzoom==0)? 1 : (nzoom==1) ? 2 : 4;
  for(i=0;i<16;i++)
    FREE(quest[i]);

  if ( nzoom==1) { /* Pas de deplacement de l'image */
    ptr_img_3d(pos)->zoom_x=0;
    ptr_img_3d(pos)->zoom_y=0;
    ptr_img_3d(pos)->zoom_z=0;
    ptr_img_3d(pos)->zoom=nzoom;
    if ( ptr_img_3d(pos)->mask )
      {
	ptr_img_3d(-pos)->zoom_x=0;
	ptr_img_3d(-pos)->zoom_y=0;
	ptr_img_3d(-pos)->zoom_z=0;
	ptr_img_3d(-pos)->zoom=nzoom;
      }
    show_picture_3d(pos);
  }
  else /* Zoom avec deplacement   */
    { // BERST : a quoi sert ce else???
      ptr_img_3d(pos)->zoom=nzoom;
      ptr_img_3d(pos)->zoom_x=0;
      ptr_img_3d(pos)->zoom_y=0;
      ptr_img_3d(pos)->zoom_z=0;
      if ( ptr_img_3d(pos)->mask )
	{
	  ptr_img_3d(-pos)->zoom_x=0;
	  ptr_img_3d(-pos)->zoom_y=0;
	  ptr_img_3d(-pos)->zoom_z=0;
	  ptr_img_3d(-pos)->zoom=nzoom;
	}
      show_picture_3d(pos);
      // TODO:GTK IMXMouse_MoveImage(pos);
    }

}

/******************************************************************************
 ** -- Zoom_Vin_3d() -----------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July, 23rd 1999
 **
 ******************************************************************************/
void Zoom_Vin_3d(void)
{ /* Zoom for visualisation... */

  int pos,max,i,e=0;
  char *quest[16];
  long nzoom;
  pos=GET_WND3D(NULL);
  max= GET_INT(TEXT0004, 0, &e);
  if (!e) {
    for(i=0;i<16;i++)
      quest[i]=CALLOC(80,char);
    strcpy(quest[0],"No reduction");
    strcpy(quest[1],"Zoom 2");
    strcpy(quest[2],"Zoom 4");
    strcpy(quest[3],"\0");
    nzoom=GETV_QCM(TEXT0043,(char **)quest);
    nzoom = (nzoom==0)? 1 : (nzoom==1) ? 2 : 4;
    for(i=0;i<16;i++)
      FREE(quest[i]);

    if(MESG_DEBUG) printf("  nzoom=%d max=%d %d \n",(int)nzoom
			  ,(int) (MAX_WIDTH/ptr_img_3d(pos)->width)
			  ,(int) (-2>2) );

    for (i=0;i<max;i++) {
      if(-2>2l) {
	if(MESG_DEBUG) printf("  nzoom=%d\n",(int)nzoom);
	ptr_img_3d(pos+i)->zoom=MAX_WIDTH/ptr_img_3d(pos+i)->width;
	ptr_img_3d(pos+i)->zoom_x = 0;
	ptr_img_3d(pos+i)->zoom_y = 0;
	ptr_img_3d(pos+i)->zoom_z = 0;
	if ( ptr_img_3d(pos+i)->mask )
	  {
	    ptr_img_3d(-pos-i)->zoom_x=MAX_WIDTH/ptr_img_3d(pos+i)->width;
	    ptr_img_3d(-pos-i)->zoom_y=0;
	    ptr_img_3d(-pos-i)->zoom_z=0;
	    ptr_img_3d(-pos-i)->zoom=nzoom;
	  }
      }
      else
	{
	  ptr_img_3d(pos+i)->zoom=nzoom;
	  ptr_img_3d(pos+i)->zoom_x = 0;
	  ptr_img_3d(pos+i)->zoom_y = 0;
	  ptr_img_3d(pos+i)->zoom_z = 0;
	  if ( ptr_img_3d(pos+i)->mask )
	    {
	      ptr_img_3d(-pos-i)->zoom_x=nzoom;
	      ptr_img_3d(-pos-i)->zoom_y=0;
	      ptr_img_3d(-pos-i)->zoom_z=0;
	      ptr_img_3d(-pos-i)->zoom=nzoom;
	    }
	}
      show_picture_3d(pos+i);
    }
  }
}

/******************************************************************************
 ** -- Zoom_Vin_3d() -----------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July, 23rd 1999
 **
 ******************************************************************************/
void Zoom_Vim_3d(void)
{ /* Zoom for visualisation max ..image */

  int pos,max,i;
  char *quest[16];
  long nzoom;

  max=MAX_PICTURES_3D;
  pos=1;
  for(i=0;i<5;i++)
    quest[i]=CALLOC(80,char);
  strcpy(quest[0],"No reduction");
  strcpy(quest[1],"Zoom 2");
  strcpy(quest[2],"Zoom 4");
  strcpy(quest[3],"\0");
  nzoom=GETV_QCM(TEXT0043,(char **)quest);
  nzoom = (nzoom==0)? 1 : (nzoom==1) ? 2 : 4;
  for(i=0;i<5;i++)
    FREE(quest[i]);

  if(MESG_DEBUG) printf("  nzoom=%d max=%d %d \n",(int)nzoom
			,(int)(MAX_WIDTH/ptr_img_3d(pos)->width)
			,(int)(-2>2) );

  for (i=0;i<max;i++) {
    if(-2>2l) {
      if(MESG_DEBUG) printf("  nzoom=%d\n",(int)nzoom);
      ptr_img_3d(pos+i)->zoom=MAX_WIDTH/ptr_img_3d(pos+i)->width;
      ptr_img_3d(pos+i)->zoom_x=0;
      ptr_img_3d(pos+i)->zoom_y=0;
      ptr_img_3d(pos+i)->zoom_z=0;
      if ( ptr_img_3d(pos+i)->mask )
	{
	  ptr_img_3d(-pos-i)->zoom_x=MAX_WIDTH/ptr_img_3d(pos+i)->width;
	  ptr_img_3d(-pos-i)->zoom_y=0;
	  ptr_img_3d(-pos-i)->zoom_z=0;
	  ptr_img_3d(-pos-i)->zoom=nzoom;
	}
    }
    else
      {
	ptr_img_3d(pos+i)->zoom=nzoom;
	ptr_img_3d(pos+i)->zoom_x=0;
	ptr_img_3d(pos+i)->zoom_y=0;
	ptr_img_3d(pos+i)->zoom_z=0;
	if ( ptr_img_3d(pos+i)->mask )
	  {
	    ptr_img_3d(-pos-i)->zoom_x=nzoom;
	    ptr_img_3d(-pos-i)->zoom_y=0;
	    ptr_img_3d(-pos-i)->zoom_z=0;
	    ptr_img_3d(-pos-i)->zoom=nzoom;
	  }
      }
    show_picture_3d(pos+i);
  }
}

/******************************************************************************
 ** -- Zoom_Vi_3d_pos_xyz() ------------------------------------------------------------
 **
 ** >> Last user action by: T. BERST on 05/11/02
 **
 ******************************************************************************/
void Zoom_Vi_3d_pos_xyz(int x, int y, int z, int pos ,int nzoom)
{
  int wdth,hght,dpth;
  int zoom_x,zoom_y,zoom_z;
  int dflt_size=DEFAULT_IMAGIX_WIDTH;

  /*	for(i=0;i<16;i++)
    	quest[i]=CALLOC(80,char);

	strcpy(quest[0],"No reduction");
	strcpy(quest[1],"Zoom 2");
	strcpy(quest[2],"Zoom 4");
	strcpy(quest[3],"Zoom 8");
	strcpy(quest[4],"Zoom 16");
	strcpy(quest[5],"Cancel");
	strcpy(quest[6],"\0");
	nzoom=GETV_QCM(TEXT0043,(char **)quest);
	for(i=0;i<16;i++)
	FREE(quest[i]);

	if (nzoom == 5) return ;
  */
  wdth = ptr_img_3d(pos)->width;
  hght = ptr_img_3d(pos)->height;
  dpth = ptr_img_3d(pos)->depth;

  ptr_img_3d(pos)->zoom	=	nzoom;

  // calcul de la position du coin superieur gauche

  zoom_x=x-dflt_size/(nzoom *2);

  if ( (zoom_x+wdth/nzoom)>wdth )
    zoom_x = wdth-wdth/nzoom;
  ptr_img_3d(pos)->zoom_x	=	MAXI(0,zoom_x);

  zoom_y=y-dflt_size/(nzoom *2);

  if ( (zoom_y+hght/nzoom)>hght )
    zoom_y = hght-hght/nzoom;
  ptr_img_3d(pos)->zoom_y	=	MAXI(0,zoom_y);

  zoom_z=z-dflt_size/(nzoom *2);

  if ( (zoom_z+dpth/nzoom)>dpth )
    zoom_z = dpth-dpth/nzoom;
  ptr_img_3d(pos)->zoom_z	=	MAXI(0,zoom_z);

  //---------------------------------//
  //-mise a jour du masque si existe //
  if (ptr_img_3d(pos)->mask)
    {
      ptr_img_3d(pos)->mask->zoom		= nzoom;
      ptr_img_3d(pos)->mask->zoom_x	= MAXI(0,zoom_x);
      ptr_img_3d(pos)->mask->zoom_y	= MAXI(0,zoom_y);
      ptr_img_3d(pos)->mask->zoom_z	= MAXI(0,zoom_z);

    }

  //----------------------------------------

  /*
    printf("%d %d %d %d\n",ptr_img_3d(pos)->zoom_x,ptr_img_3d(pos)->zoom_y,ptr_img_3d(pos)->zoom_z,nzoom);
    printf("%d %d %d\n",x,y,z);
  */
  show_picture_3d(pos);
}


/******************************************************************************
 ** -- Zoom_Vi_3d_pos() ------------------------------------------------------------
 **
 ** >> Last user action by: T. BERST on 05/11/02
 **
 ******************************************************************************/
void Zoom_Vi_3d_pos(void)
{

  int *tab,e=1;

  tab=XGET_MOUSEPOS_3D(&e);
  if (e)
    {	printf("error\n"); return ; }

  Zoom_Vi_3d_pos_xyz(tab[3], tab[4], tab[5], tab[2],ptr_img_3d(tab[2])->zoom * 2);

  FREE(tab);

}

#endif  //  __GTK_GUI
