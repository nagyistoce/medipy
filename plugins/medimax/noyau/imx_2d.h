/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/**-----------------------------------------------------------------------
***
***	file:		imx_2d.h
***
***	project:	Imagix 2.01
***
***
***	description:	Header file for project imagix (2.01)
***
***
***	Copyright (c) 1993-2000, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***     Last user action by: Mr. ARMSPACH on Feb 2003
***
***---------------------------------------------------------------------*/

#ifndef __2D_INCL
#define __2D_INCL

// IMPORTANT !: When you use a graphic structure, always set pointers
//              to a NULL value.
struct s_graphic;

/*---- Definition structure s_graphic ----*/
typedef struct	s_graphic 
{      // For graphics...
  int		nr;               // From 0 to MAX_PICTURES-1

#ifdef __GTK_GUI
  // Private information:
  GdkGC		*gc;              // Appropriate graphic context.
//  guchar	*displayable_image_data;                         !!!!!unused??????
  t_CMapInfo	cmapinfo;        // Needed data for colormap
  unsigned char alpha_channel;  // % of transparency ( 100 -> full transp.)
  int		bar;              // indicates if palette's bar is shown
  TDisplayIndexPos indexes;   // indicates if palettes index is shown
#endif // __GTK_GUI

  TImageType	img_type;        // Type of picture
  TMaskType	mask_type;        // Type of mask
  int		bitppixel;        // Bit per pixel, Assumed 32

  struct s_graphic *mask;     // Masque de l'image

  // Public information:

  int		ref_count;        // number of references on this pic
  int		pos;              // Position on display 1 to MAXPICTURES
  int		zoom;             // Zoom factor
  int		zoom_x;           // x position start zoom (zoom_x=x*zoom)
  int		zoom_y;           // y position start zoom (zoom_y=y*zoom)
  int		x2dc;				// cf 3D;
  int		y2dc;				// cf 3D;
  TYPEMRI	cutoff_min;       // Valeur min de l'affichage
  TYPEMRI	cutoff_max;       // Valeur max de l'affichage

  char		type;             // b means BRUKER, v VISILOG.
  TYPEMRI	**mri;            // Picture data.
  TDimension	width, height;   // Width & Height.
  TYPEMRI	max_pixel;        // Max_pixel de mri.
  TYPEMRI	min_pixel;        // Min_pixel de mri.
  float		icomp;            // rcoeff=2.**icomp
  float		rcoeff;           // coefficient mult pour valeur reel
  float		dx, dy, dz;         // Resolution (mm) dans les 3D.

  Patient	patient;          // Indentification du patient

  int		mask_active;      // TRUE si le masque est actif
  int		Nbpts;            // Nombre de points formant le masque
  IMXPoint	tabcoord[IMX_MAX_DRAWROI_PIXELS];  // Tableau des points inscrits sur le masque

  int		file_number;      // Numero de l'image dans le fichier
  char		filename[FILE_LEN]; // Picture filename source if any.

}		grphic, *ptr_grphic;


/* ------------------------------------------------------------------------- */
/* ---- functions prototypes ----------------------------------------------- */
/* ------------------------------------------------------------------------- */

extern grphic	*ptr_img(int nr);		/* Give pointer image. */
extern grphic	*cr_grphic(grphic *temp);
extern grphic	*cr_grphic_array(int num, int width, int height);
extern grphic	*cr_grphic_modif(int width, int height, float icomp,
				 float rcoeff, int maxpixel);
extern TYPEMRI	**cr_mri(int width, int height);
extern void	free_mri(long int **mri);
extern void	free_grphic(grphic *temp);
extern void	free_grphic_array(grphic *array, int n);

extern int	ptr_has_mask(int nr);
extern int	ptr_has_mask_p(grphic *img);
extern grphic	*ptr_mask(int pos);
extern grphic	*ptr_mask_p(grphic* img);
extern void	free_mask (int pos);
extern void	free_mask_p(grphic *mask);
extern void	imx_info_mri_2d(int nr);
extern void	IMXPict2d_CleanPicture(void);
extern void	IMXPict2d_ResetPicture(int pos);
extern int	imx_brukermax(float amax, int *maxpixel, float *icomp, float *rcoeff);
extern int	imx_inimaxpixel(int im_1);
extern int	imx_inimaxpixel_p(grphic *im1);
extern int	imx_iniminpixel(int im_1);
extern int	imx_iniminpixel_p(grphic *im1);
extern int	imx_inimaxminpixel(int im_1);
extern int	imx_inimaxminpixel_p(grphic *im1);
extern int	imx_iniparam(int im_1, long int max_pixel, long int min_pixel,
			     float icomp, float rcoeff);
extern int	imx_iniparaimg(int im_1);
extern int	imx_iniparaimg_p(grphic *im1);
extern int	imx_copie_param(int im_1, int im_res);
extern int	imx_copie_param_p(grphic *im1, grphic *imres);


/* ------------------------------------------------------------------------- */
/* ---- Global variables --------------------------------------------------- */
/* ------------------------------------------------------------------------- */
extern	grphic	_imx_roi;		/* Region of interest.		*/
extern	grphic	*_imx_img;


#endif /* Fin __2D_INCL*/
