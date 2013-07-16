/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*******************************************************************************/
/*!	   \file:		transformations_3d.c
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
#include "noyau/gtkimagix.h"
#include "noyau/imx_3d.h"
#include "noyau/imx_lang.h"
#include "noyau/mani_3d.h"
#include "noyau/gui/imx_picture3d.h"
#include "recalage/chps_3d.h"
#include "recalage/mtch_3d.h"

#include "recalage/transformations_3d.h"


/******************************************************************************
** -- imx_zoom_isotrope_3d_p() -----------------------------------------------------------
*/
/*! Zoom l'image pour que dx,dy,dz soient egaux.
**
**  \param imsrc : image a zoomer
**  \param imdst : image destination
**	\param res_type : type de resolution (min des dx, dy, dz ou 1 mm^3)
**	\param inter_type : type d'interpolation (0 plus proche, 1 lineaire, ....
**	\param save_type : type de sauvegarde (champs /transfo)
**	\param nomfichres : nom du fichier resultat
**	\param size_max : dimension maximum du 'cube' contenant l'image, 0 pour pas de taille max
**  \retval 1
**
******************************************************************************/
int	imx_zoom_isotrope_3d_p(grphic3d *imsrc, grphic3d *imdst, int res_type, int inter_type,int save_type, char* nomfichres, int size_max)
{
	transf3d *transfo=NULL;
	field3d  *champ=NULL;
	double   *param=NULL;
	int 	 wdth,hght,dpth,wold,hold,dold;
	float    dmin;
	char mesg[128];
  	int err=1;
    
	wold=imsrc->width;hold=imsrc->height;dold=imsrc->depth;
	
	if( res_type == 1 )
	{
		dmin = 1.f;
		printf(" ATTENTION : resolution finale a 1mm^3 !\n");
	}
	
	else // trouver le minimum de dx,dy,dz
	dmin = MINI(imsrc->dx, MINI(imsrc->dy, imsrc->dz));

	wdth  = (int)floor(imsrc->width*imsrc->dx/dmin);
	if (size_max && wdth>size_max)
	{
		sprintf(mesg,"ATTENTION recadrage suivant x  : wdth = %d > %d",wdth,size_max);
		PUT_WARN(mesg);
		wdth = size_max;
	}

	hght  = (int)floor(imsrc->height*imsrc->dy/dmin);
	if (size_max && hght>size_max)
	{
		sprintf(mesg,"ATTENTION recadrage suivant y  : hght = %d > %d",hght,size_max);
		PUT_WARN(mesg);
		hght = size_max;
	}

	dpth  = (int)floor(imsrc->depth*imsrc->dz/dmin);
	if (size_max && dpth>size_max)
	{
		sprintf(mesg,"ATTENTION recadrage suivant z  : dpth = %d > %d",dpth,size_max);
		PUT_WARN(mesg);
		dpth = size_max;
	}

  //creation de la transformation appropriee
	transfo=cr_transf3d(wdth,hght,dpth,NULL); if (!transfo) {fprintf(stderr, "erreur alloc memoire dans imx_zoom_isotrope_3d_p\n"); goto end_func; }
  	transfo->dx=dmin; transfo->dy=dmin; transfo->dz=dmin;
	transfo->nb_param=9;
	transfo->typetrans=RIGID3D;

	param=CALLOC(9,double); if (!param) {fprintf(stderr, "erreur alloc memoire dans imx_zoom_isotrope_3d_p\n"); goto end_func; }
	param[0]=param[1]=param[2]=0.0;

	param[3]=imsrc->dx*wold/2.0-dmin*wdth/2.0;
	param[4]=imsrc->dy*hold/2.0-dmin*hght/2.0;
	param[5]=imsrc->dz*dold/2.0-dmin*dpth/2.0;
	param[6]=dmin*wdth/2.0;
	param[7]=dmin*hght/2.0;
	param[8]=dmin*dpth/2.0;
	transfo->param=param;

 //application de la transfo
 imx_apply_transf_3d_p(imdst, imsrc, transfo, imx_choose_interpolation_fct(inter_type));
   
	//enregistrement du champ resultat
	if (nomfichres!=NULL)
	{
	  if (save_type==0)
	  {
     champ=transf_to_field_3d(transfo,NULL,NULL); if (!champ) {fprintf(stderr, "erreur alloc memoire dans imx_zoom_isotrope_3d_p\n"); goto end_func; } 
	   transfo=field_to_transf_3d(champ,imsrc,imdst);
	  }
	  save_transf_3d(transfo,nomfichres);
	}

 err=0;
  
end_func:

	if (champ) free_field3d(champ);
	if (transfo) free_transf3d(transfo);

	return err;
}

int	imx_zoom_isotrope_3d(int im_src, int im_dest, int res_type, int inter_type, int save_type, char * nomfichres)
{
	grphic3d * imsrc, * imdest;
  int err=0;
  
	imsrc = ptr_img_3d(im_src);
	imdest = ptr_img_3d(im_dest);

	imx_copie_param_3d_p(imsrc, imdest);

	err=imx_zoom_isotrope_3d_p(imsrc,imdest, res_type, inter_type, save_type, nomfichres, MAX_WIDTH_3D);

	imx_copie_visual_params_3d_p(imsrc, imdest);

	return err;
}


int quest_resolution_3d()
{
  int res_type=0;
  
  char *query[100]; 

	query[0]="resolution = min(dx, dy, dz)";
	query[1]="resolution = 1 mm";
	query[2]="\0";
	query[3]=NULL;

   res_type=GETV_QCM("Resolution ?",(char **)query);

   return res_type;
}

void zoom_isotrope_3d()
{
	int im_src,im_dest;
	int inter_type;
	char *	nomfichres=NULL;
	int save_type, res_type = 0;

	//images
	im_src  = GET_PLACE3D("image source");
	im_dest = GET_PLACE3D("image resultat");

	// zoom a dx=dy=dz=1 ou MIN(dx, dy, dz) ?
	res_type = quest_resolution_3d();

	//interpolation
	inter_type=imx_query_interpolation_type_3D(0);

	//question sur l'enregistrement du champ resultat
	nomfichres=quest_save_result_3d(&save_type);

	//zoom isotrope
	imx_zoom_isotrope_3d(im_src, im_dest, res_type, inter_type,save_type, nomfichres);

	show_picture_3d(im_dest);
  
	if (nomfichres) free(nomfichres);
}


/*******************************************************************************
********************************************************************************
************************** Zoom relatif ************************************
********************************************************************************
*******************************************************************************/

/*******************************************************************************
**        zoom_relatif_3d
**
**       zoom_relatif entre 2 image 3D
*******************************************************************************/
void zoom_relatif_3d(void)
{
 int im_ref,im_zoom,im_res;
 char /*quest[7],*/*nomfichres;
 int inter_type,save_type/*,i*/;


 /*question sur les images a recaler*/
 im_zoom=GET_PLACE3D(TEXT0027);
 im_ref=GET_PLACE3D(TEXT0031);
 im_res=GET_PLACE3D(TEXT0006);

 /*methode d'interpolation*/
 inter_type=imx_query_interpolation_type_3D(0);

 /*question sur l'enregistrement du champ resultat*/
 nomfichres=quest_save_result_3d(&save_type);

 imx_zoom_relatif_3d(im_zoom,im_ref,im_res,inter_type,save_type,nomfichres);

 if (nomfichres)
  free(nomfichres);
 show_picture_3d(im_res);

}

/*******************************************************************************
**        imx_zoom_relatif_3d
**
**       zoom_relatif entre 2 image 3D
*******************************************************************************/
int imx_zoom_relatif_3d(int im_zoom, int im_ref, int im_res, int inter_type, int save_type, char *nomfichres)
{
 grphic3d *imref,*imzoom,*imres;
 int err=0;

 imref=ptr_img_3d(im_ref);
 imzoom=ptr_img_3d(im_zoom);
 imres=ptr_img_3d(im_res);


 err=imx_zoom_relatif_3d_p(imzoom,imref,imres,inter_type,save_type,nomfichres);
 if (!err) imx_copie_visual_params_3d_p(imref, imres);

 return(0);

}

/*******************************************************************************
**        imx_zoom_relatif_3d_p
**
**       zoom_relatif entre 2 image 3D
*******************************************************************************/
int imx_zoom_relatif_3d_p(grphic3d *imzoom, grphic3d *imref, grphic3d *imres, int inter_type, int save_type, char *nomfichres)
{
 transf3d *transfo=NULL;
 field3d  *champ=NULL;
 double   *param=NULL;
 int	  wdth,hght,dpth,wold,hold,dold;

 // Controle s'il faut zoomer 
 if ( imzoom->dx == imref->dx && imzoom->dy == imref->dy && imzoom->dz == imref->dz )
   {
   if (imzoom != imres ) imx_copie_3d_p(imzoom,imres);
   return(1);  // Sortie si aucune modif 
   }

 wold=imzoom->width;hold=imzoom->height;dold=imzoom->depth;

 //on ne fait pas de subsampling => utiliser des fonctions specifique pour (cf les effets de recouvrement de spectre)
 if ((imzoom->dx<imref->dx)||(imzoom->dy<imref->dy)||(imzoom->dz<imref->dz)) { fprintf(stderr,"!!!attention risque de recouvrement de spectre\nutiliser plutot des fonctions appropriees\n"); }
 
 wdth=imref->width;hght=imref->height;dpth=imref->depth;

 transfo=cr_transf3d_p(imref,NULL);
 transfo->nb_param=9;
 transfo->typetrans=RIGID3D;

 param=CALLOC(9,double);
 param[0]=param[1]=param[2]=0.0;
 param[3]=imzoom->dx*wold/2.0-imref->dx*wdth/2.0;param[4]=imzoom->dy*hold/2.0-imref->dy*hght/2.0;param[5]=imzoom->dz*dold/2.0-imref->dz*dpth/2.0;
 param[6]=imref->dx*wdth/2.0;param[7]=imref->dy*hght/2.0;param[8]=imref->dz*dpth/2.0;
 transfo->param=param;

 imx_apply_transf_3d_p(imres, imzoom, transfo, imx_choose_interpolation_fct(inter_type));
  
 //enregistrement du champ resultat
 if (nomfichres!=NULL)
  {
   if (save_type==0)
    {
    champ=transf_to_field_3d(transfo,NULL,NULL);  
    transfo=field_to_transf_3d(champ,imref,imzoom);
    }
   save_transf_3d(transfo,nomfichres);
  }

 if (champ) free_field3d(champ);
 if (transfo) free_transf3d(transfo);

 return(0);

}

/*******************************************************************************
********************************************************************************
************************** Rotation de l'image ************************************
********************************************************************************
*******************************************************************************/

/*******************************************************************************
**        rot_3d
**
**       rotation de l'image 3D
*******************************************************************************/
void rot_3d()
{
 int inter_type, err=0;
 double ang1,ang2,ang3;
 int im_deb,im_res;

 //images a transformer et resultat
 im_deb=GET_PLACE3D(TEXT0001);
 im_res=GET_PLACE3D(TEXT0006);

 //angles de rotation
 ang1=(double) GET_FLOAT("Angle suivant X ?", 0, &err);
 ang2=(double) GET_FLOAT("Angle suivant Y ?", 0, &err);
 ang3=(double) GET_FLOAT("Angle suivant Z ?", 0, &err);

 //interpolation
 inter_type=imx_query_interpolation_type_3D(0);

 imx_rot_3d(im_deb,im_res,inter_type,ang1,ang2,ang3);

 show_picture_3d(im_res);
}

/*******************************************************************************
**        imx_rot_3d
**
**       rotation d'une image 3D
*******************************************************************************/
int imx_rot_3d(int im_ref, int im_res, int inter_type, double ang1, double ang2, double ang3)
{
 grphic3d *imref,*imres;

 imref=ptr_img_3d(im_ref);
 imres=ptr_img_3d(im_res);


 imx_rot_3d_p(imref,imres,inter_type,ang1,ang2,ang3);

 return(1);

}

/*******************************************************************************
**        imx_rot_3d_p
**
**       rotation d'une imageimage 3D
*******************************************************************************/
int imx_rot_3d_p(grphic3d *imref, grphic3d *imres, int inter_type, double ang1, double ang2, double ang3)
{
 grphic3d *imtemp,*imtres;
 transf3d *transfo;
 field3d  *champ;
 double   *param;
 int wdth,hght,dpth;


 /* Creation de imtemp, copie de im_ref dans imtemp*/
 imtemp=cr_grphic3d((grphic3d *)NULL);
 imx_copie_3d_p(imref,imtemp);
 imtres=cr_grphic3d(imtemp);

 wdth=imref->width;hght=imref->height;dpth=imref->depth;

 transfo=cr_transf3d(wdth,hght,dpth,NULL);
 transfo->nb_param=9;
 transfo->typetrans=RIGID3D;

 param=CALLOC(9,double);
 param[0]=ang1*PI/180.;
 param[1]=ang2*PI/180.;
 param[2]=ang3*PI/180.;
 param[3]=param[4]=param[5]=0.0;
 param[6]=wdth/2;param[7]=hght/2;param[8]=dpth/2;

 transfo->param=param;

 champ=transf_to_field_3d(transfo,NULL,NULL);

 (imx_choose_interpolation_fct(inter_type))(imtemp,champ,imtres);

 imx_copie_3d_p(imtres,imres);

 free_field3d(champ);
 free_transf3d(transfo);
 free_grphic3d(imtemp);
 free_grphic3d(imtres);

 return(1);

}

int transf_geom_to_field3d(transf3d *transfo, field3d *champ, grphic3d *imref, grphic3d *imtransf, field3d* points_calc)
{
  vector3d ***data=NULL;
  double param[15]={0.0};

  int i,j,k,wdth,hght,dpth;
  double x,y,z;
  double tx,ty,tz;
  double a11,a12,a13,a21,a22,a23,a31,a32,a33;
  double tpx11,tpx21,tpx31,tpy12,tpy22,tpy32;
  double dxRef, dyRef, dzRef;
  double dxReca, dyReca, dzReca;
  vector3d ***pts_calc=NULL;

  if ((!transfo)||(!champ))
  { fprintf(stderr, "manque de parametres dans transf_geom_to_field3d\n"); return 1; }
  
  if (points_calc) { wdth=(int)points_calc->width;hght=(int)points_calc->height;dpth=(int)points_calc->depth; }
  else
 {
  wdth=(int)transfo->width;hght=(int)transfo->height;dpth=(int)transfo->depth;
  //verification du champ
  if ((champ->width != wdth) || (champ->height != hght) || (champ->depth != dpth)) return 1;
 }

  data=champ->raw;

  //cas gere pour la compatibilte avec les anciens champs isotropes
  if ((!imref)&&(!imtransf))
  {
   dxRef=1.0; dyRef=1.0; dzRef=1.0;
   dxReca=1.0; dyReca=1.0; dzReca=1.0;
  }
  else
  {
   //si on a imref c'est lui qui a la priorite pour les parametres d'anisotropie
   if ((imref)&&(imref->dx)&&(imref->dy)&&(imref->dz))
   { dxRef=imref->dx; dyRef=imref->dy; dzRef=imref->dz; }
   else if ((transfo->dx)&&(transfo->dy)&&(transfo->dz))
   { dxRef=transfo->dx; dyRef=transfo->dy; dzRef=transfo->dz; }
   else { dxRef=1.0; dyRef=1.0; dzRef=1.0; }

   if ((imtransf)&&((imtransf->dx)&&(imtransf->dy)&&(imtransf->dz)))
   { dxReca=imtransf->dx; dyReca=imtransf->dy; dzReca=imtransf->dz; }
   else { dxReca=1.0; dyReca=1.0; dzReca=1.0; }
  }

  for (i=0;i<transfo->nb_param;i++) param[i]=transfo->param[i];

  //on se ramene a une transformation affine centree en (0;0;0)
  switch (transfo->typetrans)
  {
   case RIGID3D: rigid_to_rigid_global_zoom_3d(param);
   case RIGIDGLOBALZOOM3D: rigid_global_zoom_to_rigidz_3d(param);
   case RIGIDZOOM3D: rigidz_to_affine_decouple_3d(param);
   case AFFINEDECOUPLE3D: affine_decouple_to_affine_3d(param);
   case AFFINE3D : affine_to_affinesscg_3d(param); break;
   default: fprintf (stderr, "SHOULD NOT BE THERE IN transf_geom_to_field3d\n"); return 2;
  }

  a11=param[0];a12=param[1];a13=param[2];
  a21=param[3];a22=param[4];a23=param[5];
  a31=param[6];a32=param[7];a33=param[8];
  tx=param[9];ty=param[10];tz=param[11];

  tx/=dxReca; ty/=dyReca; tz/=dzReca;
  a11=a11*dxRef/dxReca; a12=a12*dyRef/dxReca; a13=a13*dzRef/dxReca;
  a21=a21*dxRef/dyReca; a22=a22*dyRef/dyReca; a23=a23*dzRef/dyReca;
  a31=a31*dxRef/dzReca; a32=a32*dyRef/dzReca; a33=a33*dzRef/dzReca;

  a11=a11-1.0; a22=a22-1.0; a33=a33-1.0;
  if (points_calc==NULL)
  {
   for (i=0;i<wdth;i++)
   {
    x=i;tpx11=a11*x;tpx21=a21*x;tpx31=a31*x;
    for (j=0;j<hght;j++)
    {
     y=j;tpy12=a12*y;tpy22=a22*y;tpy32=a32*y;
     for (k=0;k<dpth;k++)
     {
      z=k;
      data[i][j][k].x=(float)(tpx11+tpy12+a13*z+tx);
      data[i][j][k].y=(float)(tpx21+tpy22+a23*z+ty);
      data[i][j][k].z=(float)(tpx31+tpy32+a33*z+tz);
     }
    }
   }
  }
  else
  {
   pts_calc=points_calc->raw;
   for (i=0;i<wdth;i++)
   {
    for (j=0;j<hght;j++)
    {
     for (k=0;k<dpth;k++)
     {
      x=pts_calc[i][j][k].x;
      y=pts_calc[i][j][k].y;
      z=pts_calc[i][j][k].z;
      data[i][j][k].x=(float)(a11*x+a12*y+a13*z+tx);
      data[i][j][k].y=(float)(a21*x+a22*y+a23*z+ty);
      data[i][j][k].z=(float)(a31*x+a32*y+a33*z+tz);
     }
    }
   }
  }

 return 0;
}

