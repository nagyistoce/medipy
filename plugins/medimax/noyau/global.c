/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*!-----------------------------------------------------------------------
***	
***	\file:		global.c
***
***	project:	Imagix 2.01 
***			
***  \brief  Global variables for project.
***	
***	
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***---------------------------------------------------------------------*/

#include <config.h>
#include "imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_3d.h"

#ifdef __GTK_GUI
# include "noyau/imx_2d.h"
# include "recalage/imx_reca.h"
# include "serie/imx_seri.h"
# include "serie/seri_3d.h"
# include "outils/imx_stck.h"
# include "recalage/reca_3d.h"
# include "menu/menu.h"
#endif /*__GTK_GUI*/

#ifndef __IMAGIX_INCL_GLOBAL
#define __IMAGIX_INCL_GLOBAL

char    *_home;
char    *_imx_home;
char    *_name_header_file;
char    *_default_path;
char    *_tala_atlas_path;

char	_CurrentPath[256];

//#if defined (WIN32)
//int		errno;
//#endif

#endif

#ifndef IMAGIX_PARTIAL
#define IMAGIX_FULL
#endif

#ifdef IMAGIX_FULL
char	_is_bigEndian       = 1;    /* Big par defaut (HP) */
int   _is_filebigEndian   = FALSE; /*! Load File ...big or little endian*/
#endif // IMAGIX_FULL

UINT    _seuil_super_img    = 200;
float   _cutoff_max         = 0.;
int     _cutoff_img_nr      = 1;


#ifndef __IMAGIX_BEGIN__SET_GLOBALS
#define __IMAGIX_BEGIN__SET_GLOBALS

int	_automatic_mode	=ATM_MANUEL;
int	_prg_style	=ATM_PRG_STYLE__PATH;

int	_imagix_lst_itm=0;
int	_imagix_error=0;
int	_imagix_err_ext=0;

/*int	_full_reserved=ON;*/
int	_obj_color      = 4;            /* Initialisation couleur objet par default */
int	_txt_color      = 5;            /* Initialisation couleur texte par default */
int	_tala_color     = 6;           /* Couleur de la grille de Talairach. */
int	_gomme_color    = 0;            /* Initialisation couleur gomme par default */
int	_init_roi       = TRUE;         /* Initialisation RAZ ROI avant creation d'une ROI */
int	_read_format_3d = 2;            /* 1=LUKE, 2=IPB, 3=Rouffach  */

FILE   *_fp_auto =(FILE *)NULL;
FILE   *_fp_learn=(FILE *)NULL;
FILE   *_fp_lotus=(FILE *)NULL;

char	_atm_autofilename[FILE_LEN];
char    _lotus_filename  [FILE_LEN];


void  **_global_buff=NULL;

/* Bit 0 Reserved // ACTIVE OR NOT LOTUS FILE;
   Bit 2 ON   "   (...);
   Bit 4 ON   "   PUT_MSGWND(...);
   Bit 8 ON shows PUT_MESG(...);
   Bit10 ON   "   PUT_ERR(...);
   Bit11 ON   "   PUT_WARN(...);
   Bit12 ON   "   PUT_INFO(...);
   Bit13 ON   "   PUT_ERRFAT(...);
*/
unsigned long	_dialog_messages_status=0x00003f11l;

/*  */
#ifdef __GTK_GUI
//grphic	_imx_roi;
grphic	*_imx_img;
#endif /*__GTK_GUI*/

/* 4d image(s) List*/
t_ptrlist _img4d=NULL;
/*---*/

/* (2+1)d image(s) List*/
t_ptrlist _img3d=NULL;
/*---*/

/* Need in zone */
IMXPoint	_tabcoord[IMX_MAX_DRAWROI_PIXELS];
int		_Nbpts=0;


/* Variables fMRI. */
#ifdef __GTK_GUI
ptr_serie	_fmri = NULL;
ptr_serie_3d	_fmri_3d = NULL;
#endif /*__GTK_GUI*/

/* angles de rotation en angio*/
int _x_rotate = 0;
int _y_rotate = 0;
int _z_rotate = 0;

/* anciennes variables fMRI, encore utilisees pour cardio, etc... */
float	_paradigm[MAX_IMG_SERIE];
char	_serie_file[MAX_IMG_SERIE][FILE_LEN];
int	_serie_file_number[MAX_IMG_SERIE];
char	_serie_type[MAX_IMG_SERIE];
int	_nbpts_serie=0;

#ifdef TCLTK
int	_serie_file_creation_type=4;
#else
int _serie_file_creation_type=1;
#endif

#ifdef TCLTK
int	_serie_file_creation_type_3d=4;
#else
int _serie_file_creation_type_3d=1;
#endif

char	_path_paravision[FILE_LEN]; /* Contains path for paravision file */

/* variables pour statistiques */

int		_nbpts_stat=0;
double	_sx[MAX_WIDTH][MAX_HEIGHT];
double	_sy[MAX_WIDTH][MAX_HEIGHT];
double	_sx2[MAX_WIDTH][MAX_HEIGHT];
double	_sy2[MAX_WIDTH][MAX_HEIGHT];
double	_sxy[MAX_WIDTH][MAX_HEIGHT];
float	_x[MAX_IMG_SERIE];
float	_y[MAX_IMG_SERIE];
int		_nbpts_xy;                      
float	_x_histo[GRPHC_WDTH];
float	_y_histo[GRPHC_HGHT];
float	_xmin_histo=0.;
float	_xmax_histo=0.;
int	_nbpas_histo;

/* Initialisation histo cumule 3d */
float	*_x_histo_3d=NULL;
float	*_y_histo_3d=NULL;
float	_xmin_histo_3d=0.;
float	_xmax_histo_3d=0.;
int	_nbpas_histo_3d;

int		_mask_color = 5; 
int		_mask_color24b = 3134715;
int		_threshold_color = 1;
int 	_threshold_color24b = 16777215;

/*  Initialisation pour cardio  */
int		_speed_movie;
int		_mip_tag;
int		_movie_cardiowithroi=FALSE;
int		_movie_cardiocolorroi=5l;
int		_cardio_save_zoom=FALSE;

/*  Initialisation pour recalage  */
double	_sigma_robust=0.;



char	_messages[MESSAGES_LEN];	/* Contains message in mesg dialog */
char	_warning[WARNING_LEN];		/* Contains warning in mesg dialog */

long	 _lim_label=LIM_LABEL;  /* taille max du tableau pour imx_labelzi  */


#ifdef __GTK_GUI
GtkWidget	*_cmap;				/*Colors  	     */

/*  Initialisation pour X11  */
GdkGC _context;

/*  Initialisation pour objets  */
struct	stack  *_imagix_objects=NULL;
struct	stack  *_imagix_objects_3d=NULL;
struct	stack  *_imagix_texts  =NULL;
struct	stack  *_imagix_oft    =NULL;
struct	stack  *_loop_cmd__pos =NULL;
struct	stack  *_imagix_mire=NULL;

#endif /*__GTK_GUI*/




/*  Initialisation pour recalage */

#ifdef __GTK_GUI
double _limited_search_space[NB_PARAM_MT]= {
          			            0.5 , 0.5 , 0.5 , 0.2
			                   };

double _limited_search_space_recuit[NB_PARAM_MT]= {
			                           0.3 , 0.3 , 0.5 , 0.1
			                          };

double _limited_search_space_3D[NB_PARAM_MT_3D]= {
			                       0.5 , 0.5 , 0.5 ,
			                       0.5 , 0.5 , 0.5 ,
			                       0.8 , 0.8 , 0.8 
			                      };

double _limited_search_space_recuit_3D[NB_PARAM_MT_3D]= {
			                       0.3 , 0.3 , 0.3 ,
			                       0.5 , 0.5 , 0.5 ,
			                       0.1 , 0.1 , 0.1
			                      };
#endif /*__GTK_GUI*/

/* Pour la 3d */
grphic3d     *_imx_img3d;
grphic3d     _imx_roi3d;

/*Liste pour les images temporaires (dans les automates..)*/
t_ptrlist _img_tmp_3d=NULL;
/*---*/

T_drawing_3d _drawing_3d;

#ifdef __GTK_GUI
/* pour l'interfacage avec ImLib3D (repertoire cplusplus) */
# ifndef HAS_IMLIB3D
char *_imlib3d_improc_callback_name="imlib3d_improc_";
void imlib3d_add_menu_items(TMenuItem **pitems,TMenuItem *items0)
{
  printf("imlib3d_add_menu_items: no ImLib3D support compiled\n");
  *pitems = items0;
}
# endif

#endif /*__GTK_GUI*/

// Pour mtch_robust_3d.c //
 double _moy_imref_3d;
 double _moy_imreca_3d;

// double _ect_imref_3d;
// double _ect_imreca_3d;

char __get_value_buffer[256];


char* imx_int_to_str(int val) 
{ 
  sprintf(__get_value_buffer, "%d", val); 
  return (char*)__get_value_buffer;
}


char* imx_double_to_str(double val) 
{ 
  sprintf(__get_value_buffer, "%f", val); 
  return (char*)__get_value_buffer;
}

char* imx_float_to_str(float val) 
{ 
  sprintf(__get_value_buffer, "%f", val); 
  return (char*)__get_value_buffer;
}


#endif

