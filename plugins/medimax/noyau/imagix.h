/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/**-----------------------------------------------------------------------
***	
***	file:		imagix.h
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
***     Last user action by: Mr. ARMSPACH on Aug 29th 1994
***     Last user action by: Mr. ARMSPACH on Jan 29th 1995
***     Last user action by: Mr. ARMSPACH on Jan 29th 1996
***     Last user action by: Mr. ARMSPACH on Nov 21th 1996
***     Last user action by: Mr. BERST    on Feb 21th 2003
***
***---------------------------------------------------------------------*/

#ifndef __IMAGIX_INCL
#define __IMAGIX_INCL

#ifdef WIN32
#include <windows.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <math.h>

#ifdef __cplusplus
extern "C"{
#endif

#include "noyau/imx_types.h"
#include "noyau/imx_defines.h"
#include "outils/imx_list.h"
#include "outils/imx_stck.h"

/* ------------------------------------------------------------------------- */
/* ---- Constants ---------------------------------------------------------- */
/* ------------------------------------------------------------------------- */

/*  Definition des noms de programme externe */
#define XWPICK			"xwpick"
#define PNMTOTIFF		"pnmtotiff"
#define XWDTOPNM		"xwdtopnm"
#define W3				"netscape"

#define	MESSAGES_LEN	2048	/* Must be greater than 256	*/
#define	WARNING_LEN		2048	/* Must be greater than 256	*/

/* definitions globals */
#define	OFF			FALSE
#define	ON			!OFF

#define	LANG_FR		ON
#define	LANG_US		OFF
#ifdef LINUX
#define IMX_TMP_FOLDER    "/tmp"
#else 
#define IMX_TMP_FOLDER    "/var/tmp"
#endif /*LINUX*/
#define ATM_MAX_TMP_IMAGE 2000
#define ATM_MIN_TMP_IMAGE 100

/*  Definition des tailles de Medimax */
#define	DEFAULT_IMAGIX_WIDTH	(256)	/*   Taille initiale de	*/
#define	DEFAULT_IMAGIX_HEIGHT   (256)	/* la fenetre imagix.	*/

#define	MAX_WIDTH		((long)256) 	/* Largeur maxi d'1 img	multiple 4*/
#define	MAX_HEIGHT		((long)256) 	/* Hauteur maxi...	multiple 4*/
#define	MAX_WIDTH_3D	((long)256) 	/* Largeur maxi d'1 img	*/
#define	MAX_HEIGHT_3D	((long)256) 	/* Hauteur maxi...	*/
#define	MAX_DEPTH_3D	((long)256) 	/* Profondeur maxi...	*/

#define	DEFAULT_WIDTH_PICTURE	(MAX_WIDTH)
#define	DEFAULT_HEIGHT_PICTURE	(MAX_HEIGHT)

#define	MAX_PIC_IN_LINE		 (4)	/* Nombre maxi d'images dans une ligne.	*/
#define	MAX_PIC_IN_COLUMN 	 (4)	/* Nombre maxi d'images dans une colonne*/
#define	MAX_PICTURES		 (MAX_PIC_IN_LINE*MAX_PIC_IN_COLUMN)
#define	MAX_PIC_3D_IN_LINE 	 (2)	/* Nombre maxi d'images dans une ligne.	*/
#define	MAX_PIC_3D_IN_COLUMN     (2)	/* Nombre maxi d'images dans une colonne*/
#define	MAX_PICTURES_3D		 (MAX_PIC_3D_IN_LINE*MAX_PIC_3D_IN_COLUMN)
/*#define	IMX_GRAPHIC_WND_X	(1280-1024)*/
#define	IMX_GRAPHIC_WND_X	132	/* 165,133,100 When menu is VERTICAL */
#define	IMX_GRAPHIC_WND_Y	0		/* 38 	*/


#define	IMX_MAX_DRAWROI_PIXELS	4096

#define GRPHC_WDTH      256     /* Nb point x pour graphique (histo, ..) */
#define GRPHC_HGHT      256     /* Nb point y pour graphique (histo, ..) */

#define MAX_IMG_SERIE	128	/* Nombre max d'image pour les traitements en serie */
#define MAX_GRAPH		4	/* Nombre max de graphe a l'ecran */
#define MAX_TEXT_LEN	255	/* Nombre max de graphe a l'ecran */
 
/*  Determination of MACRO        */
#define	NUMBER(c)		(c>='0' && c<='9')
#define	MAJUSCULE(c)	(c>='A' && c<='Z')
#define	MINUSCULE(c)	(c>='a' && c<='z')
#define	SPECIAL(c)		(c=='*' || c=='.'  || c=='?' || c=='-' || c=='/')
#define	ASCII(c)		(MAJUSCULE(c) || MINUSCULE(c) || SPECIAL(c))
#define	NUMERIC(c)		(NUMBER(c))
#define	ALPHANUM(c)		(ASCII(c) || NUMBER(c))
#define	MINI(a,b)		((a>b)? (b):(a))
#define	MAXI(a,b)		((a>b)? (a):(b))
#define SIGN(a)         ((a)>=0.0 ? 1 : -1 )
#define DSIGN(a)        ((a)>0.0 ? 1.0 : -1.0 )
#define	SWAP(a,b)		{ a^=b; b^=a; a^=b; } /* swap a and b */
#define	LERP(a,l,h)		((l)+(((h)-(l))*(a))) /* linear inter from l to h */
#ifndef CLAMP
#define	CLAMP(v,l,h)	((v)<(l) ? (l) : (v) > (h) ? (h) : v)/* clamp v to the specified range */
#endif /* CLAMP */
#define GINT(n,d)		((div((int)n,(int)d).rem==0)? (div((int)n,(int)d).quot): (div((int)n,(int)d).quot+1))
#define	XOR(a,b)		(((a)||(b))&&(!((a)&&(b))))
#define	CALLOC(n,t)		((t*)calloc(n+1,sizeof(t)))
#define	MALLOC(n,t)		((t*)malloc(sizeof(t)*(n+1)))
#define	REALLOC(p,n,t)	((t*)realloc(p,sizeof(t)*(n+1)))
#define	FREE(p)			if(p) { free(p); p = NULL; }
#define CARRE(x)    ((x)*(x))


/* Imagix settings, controls */
#define	STATIC				OFF	/* Not change ! */
#define	DYNAMIC				!STATIC
#define	REFRESH				2 	/* See _refresh_status.	   	   */
#define	REFRESH_DRAWING_AREA 8 	/* See _refresh_status.		   */
#define	MESG_CLIENTS		OFF /* See or not client's mesg (0:Not)*/
#define	MESG_CONTROLER		OFF /* See or not evt.          (0:Not)*/ 
								/* (Doesn't control client's mesg  */
#define	REFRESH_CONTROLER	OFF	/* See or not refresh evt (0:Not)  */
#define	GV_EVT_CONTROLER	OFF	/* See or not getvalue evt. (0:Not)*/
#define	DRAW_CONTROLER		OFF /* See or not when draw action append*/
#define	MENU_CONTROLER		OFF	/* See or not text menu (0:Not)	   */
#define	STCK_CONTROLER		OFF /* See when pushing, pulling, etc  */
#define	MESG_CALLBACK		OFF /* See or not usr & itm I/O	   */
#define	MESG_DEBUG			OFF /* See or not debug messages	   */
#define PROG_DEBUG			OFF /* See or not debug messages specific to the programmer */
#define RECA_DEBUG			OFF	/* See or not debug mess registration */
#define MESG3D_DEBUG		ON  /* See or not debug mess registration */
#define TEST_DEBUG			ON  /* See or not debug mess test */

#define IMG_3D_USE			ON	/* See if 3D in use 	   */
#define	COLORbyCOLOR	  	ON  /* If you want to manage each color  */
								/* When ON, X11 settings takes times*/
#define TEST				OFF	/* For the test 	   */


#define	BRUKER_FORMAT	(0x112)
#define	VISILOG_FORMAT	(0x113)
#define	INTERFILE_FORMAT (0x114)
#define	PCX_FORMAT		(0x115)
#define	PCXBW_FORMAT	(0x116)
#define	TIF_FORMAT		(0x117)
#define	BMP_FORMAT		(0x11A)
#define	GIF_FORMAT		(0x120)
#define	PICT_FORMAT		(0x121)
#define	PPM_FORMAT		(0x122)
#define	RAW4_FORMAT		(0x123)
#define	PARAVISION_FORMAT (0x124)
#define IPB_FORMAT      (0x125)
#define IPB_2D_FORMAT   (0x126)
#define IPB_3D_FORMAT   (0x127)
#define SIEMENS_FORMAT  (0x128)
#define VANDERBILT_FORMAT (0x129)
#define SMIS_FORMAT     (0x12A)
#define ACRNEMA1L_FORMAT (0x12C) 
#define ACRNEMA1B_FORMAT (0x12D)
#define ACRNEMA2L_FORMAT (0x12E) 
#define ACRNEMA2B_FORMAT (0x12F)
#define ECAT7_FORMAT     (0x130)
#define DICMBIMP_FORMAT  (0x131) /*DICOM BigEndian Implicit */
#define DICMBEXP_FORMAT  (0x132) /*DICOM BigEndian Explicit */
#define DICMLIMP_FORMAT  (0x133) /*DICOM LittleEndian Implicit */
#define DICMLEXP_FORMAT  (0x134) /*DICOM LittleEndian Explicit */
#define DICM3DBIMP_FORMAT  (0x135) /*DICOM 3D BigEndian Implicit */
#define DICM3DBEXP_FORMAT  (0x136) /*DICOM 3D BigEndian Explicit */
#define DICM3DLIMP_FORMAT  (0x137) /*DICOM 3D LittleEndian Implicit */
#define DICM3DLEXP_FORMAT  (0x138) /*DICOM 3D LittleEndian Explicit */
#define ANALYZE_FORMAT_LITTLE_ENDIAN   (0x140) /*Format ANALYSE  little endian*/
#define ANALYZE_FORMAT_BIG_ENDIAN   (0x142) /*Format ANALYSE big endian*/
#define	HEADER_SIZE		((long)18432)	/* Value for Bruker format...	*/
#define PS_FORMAT		(0x139)
#define JPG_FORMAT		(0x140)
#define PSBW_FORMAT		(0x141)
#define IPB_WITHOUT_MASK_FORMAT		(0x142) /*utilise comme un sous format*/ 
#define IPB_WITH_MASK_FORMAT		(0x143) /*utilise comme un sous format*/ 

/* Loading & saving modes */
#define MODE_2D         0
#define MODE_3D         1

/* Imagix errors...	*/
#define	IMX_NOT_ENOUGH_CORE	1
#define	IMX_FILE_NOT_FOUND	2
#define	IMX_FILE_OPEN_ERROR	3
#define	IMX_PICTURE_NOT_FOUND	4
#define	IMX_NOT_A_PICTURE_FILE	5
#define	IMX_BAD_POSITION	6
#define	IMX_UNKNOWN_CMD		7
#define	IMX_UNKNOWN_CB		8
#define	IMX_UNKNOWN_TYPE	9

#define	IMX_NO_MenuCB		10
#define	IMX_NO_UserCB		11
#define	IMX_NO_ACCEPT		12
#define	IMX_NO_CANCEL		13
#define	IMX_NO_HELP			14



/* Imagix evenements...	*/
#define	IMX_EVT_END_INIT	108	/* Min. # of events to set imagix. */

#define	IMX_CLIENT_EVT		(ClientMessage)	/* See Xlib...	*/
#define	IMX_CALLBACK_EVT	(LASTEvent+1)
#define	IMX_LAST_EVT		(LASTEvent+2)

#define	IMX_END_GETVALUE	"Getvalue() has terminated"
#define	IMX_BEGIN_GETVALUE	"Getvalue() stops other evt"

#define	ATM_NEXT_CMD_MESG	"Ready__send_next_cmd"
#define	ATM_NO_MORE			"No_more_command"

#define	IMX_COLORMAP 		"imagix_colormap"


/* pour imx_menu2d.c*/
/* pour mani_3d.c */
#define MIRE_BCGRND 1
#define MIRE_FRGRND 2
#define MIRE_CIRCLE 3
#define MIRE_BOX 4
#define MIRE_DFT 25
#define SUP_CLMAX 5
#define SUP_DIMIMA 6
#define SUP_SYMIMA 7
#define SUP_EQUA 8
#define SUP_STEP 9
#define SUP_CLMAX_3D 10
#define SUP_EQUA_3D 11
#define SUP_EPI_3D 12
#define SUP_FMRI_3D 13
#define SUP_SEP_3D 14
#define SUP_EPI 15
#define SUP_FMRI 16
#define SUP_SEP 17
#define SUP_TRANSP_3D 18
/* Colors...	*/
#define	NB_RESERVED_COLORS	(7)
#define NCOLOR_FOR_TEXT		1
#define NCOLOR_FOR_OBJECT	5

/*   Definition pour l'algorithme de otsu */
#define GREY_LEVEL_NB		256
#define MAX_GREY_LEVEL		255
#define MAX_THRESHOLD_NB	5

/*    Constantes   */
#if defined(LINUX) || defined(__GNUC__) || defined(DEC)  || defined(HPUX) || defined(SGI) || defined(WIN32)
#define PI		3.141592654	/*   constante pi */
#endif
#define FLOATMAX 20000000
#define LONGMAX	2000000000
#define	PATH_LEN	512
#define	FILE_LEN	256
#define	NAME_LEN	 40	
#define	DATE_LENGTH	 35	
#define BAR_LENGTH    4

#define LIM_LABEL 30000  /* taille max du tableau pour imx_labelzi  */
#define LIM_TAB   2048  /* taille max d'un tableau   */

/*   For MouseButtonPressedCB Callback   */
#define MOUSEPOS	1
#define PTSBYPTS	2

/* Automate: */
#define	ATM_MANUEL	(1)
#define	ATM_LEARN	(2)
#define	ATM_AUTO	(4)
#define	ATM_REMOTE  (8)

#define	ATM_PRG_STYLE__PATH         (1)
#define	ATM_PRG_STYLE__REQUEST      (2)
#define	ATM_PRG_STYLE__PATHnREQUEST (4)
#define	ATM_PRG_STYLE__SWITCH_CODE  (8)

#define	ATM_PRG_STYLE__CURRENT      (8)

#define	ATM_SHOWALL     (OFF)   /* When you want to see each line read	*/
#define	ATM_SHOW_CMD	(OFF)   /* When you want to see commands executed */
#define	ATM_SHOW_LRN	(OFF)   /* When you want to see commands learnt */

#define	ATM_FILE__COMMENTARY    ('#')	 /* It's a char. */
#define	ATM_FILE__VARIABLE      ('@')	 /* It's a char.	*/
#define	ATM_FILE__CURRENT       ('$')	 /* It's a char.	*/
#define	ATM_FILE__WAIT_CMD      ("WAIT") /* Must be a string */
#define	ATM_FILE__LOOP_CMD      ("LOOP") /* Must be a string */
#define	ATM_FILE__END_CMD       ("END")	 /* Must be a string */

/* ------------------------------------------------------------------------- */
/* ---- Definition des structures ------------------------------------------ */
/* ------------------------------------------------------------------------- */



/*---- Definition structure s_open_file_table ----*/
typedef struct s_open_file_table { /* Which file are opened ??? */
  long position; 
  char file[256];
} oft;

/*---- Definition structure s_loop ----*/
typedef struct s_loop { /* When loop has been called, use this...	*/
  long  position; 
  long	variable;	/* Represents loop nr... 	*/
} t_loop;

/*---- Definition structure imx_colors ----*/
struct	imx_colors {
  int	nr;
  char  name[32];
  long	red;
  long	green;
  long  blue;
} ;


/* Structure pour la lecture d'une image au format PCX */
/*---- Definition structure s_pcx ----*/
typedef struct s_pcx
{
  int Version;
  int BitsPerPixel;
  int Width;
  int Height;
  int Planes;
  int BytesPerLine;

  unsigned char colormap[256][3];
  unsigned char *pixels;
} Pcx_Descr;




/* ------------------------------------------------------------------------- */
/* ---- functions prototypes ----------------------------------------------- */
/* ------------------------------------------------------------------------- */


	extern void		aff_cube(void);

	/* --------------------------------------------	imx_file.c: Very Useful	 */
	extern int	load_picture(int pos, char *file, int numfirst);		/* Load a picture */
	extern int	show_picture(int pos);		/* Show a picture */
	extern int	save_picture(int pos, char *file);		/* Save a picture */

	extern void	IMXExit(void);			/* Stop functions when error. */

	/* --------------------------------------------	imx_mouse.c: ----------- */
/*  #endif */


/* ------------------------------------------------------------------------- */
/* ---- macros ------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */
#define CMPSTR(str)	      		XmStringCreateLtoR (str, XmFONTLIST_DEFAULT_TAG)



/* ------------------------------------------------------------------------- */
/* ---- Global variables --------------------------------------------------- */
/* ------------------------------------------------------------------------- */

extern long	_lim_label;

extern  char    *_home;
extern  char    *_imx_home;
extern  char    *_name_header_file;
extern  char    *_default_path;
extern  char    *_tala_atlas_path;

extern	char	_is_bigEndian;
extern  int   _is_filebigEndian; 

extern	IMXPoint	_tabcoord[];	/* Buffer for points.		*/
extern  int		_Nbpts;			/*   which contains Nbpts points*/
extern  void	**_global_buff;	/* Contains values (for Getv).	*/

extern	int		_imagix_lst_itm;	/* Last item executed.		*/
extern	int		_imagix_error;		/* Internal imagix error...	*/
extern	int		_imagix_err_ext;	/* Content errno when err occurs*/

extern  int     _background_clr;
extern  int     _foreground_clr;

extern  int     _x_rotate;		/* rotate angles */
extern  int     _y_rotate;
extern  int     _z_rotate;

extern	int		_nbpts_stat;	/* Nb of pts for statistics */
extern	double	_sx[MAX_WIDTH][MAX_HEIGHT];  /* Images pour les stat serie */
extern	double	_sy[MAX_WIDTH][MAX_HEIGHT];  /* Images pour les stat serie */
extern	double	_sx2[MAX_WIDTH][MAX_HEIGHT];
extern	double	_sy2[MAX_WIDTH][MAX_HEIGHT];
extern	double	_sxy[MAX_WIDTH][MAX_HEIGHT];

extern	float	_x[MAX_IMG_SERIE];  /* X Area for serial treatment */
extern	float	_y[MAX_IMG_SERIE];  /* Y Area for serial treatment */
extern	int		_nbpts_xy;	    	/* Nb of points for the fonction _y=f(_x)*/
extern	float	_paradigm[MAX_IMG_SERIE]; /* Area for serial treatment */
extern	char	_serie_file[MAX_IMG_SERIE][FILE_LEN];/*Nom du file associated*/
extern	int		_serie_file_number[MAX_IMG_SERIE]; /*Numero img in the file*/
extern	char	_serie_type[MAX_IMG_SERIE];  /*define type (R) rest (A) Actif*/
extern	int		_nbpts_serie;				/* nombre d'images */
extern	int		_serie_file_creation_type; /* Can select la suite de fichier */
extern	int		_serie_file_creation_type_3d;/* Can select la suite de fichier*/

extern	float	_x_histo[GRPHC_WDTH];        /* X Area for serial treatment */
extern	float	_y_histo[GRPHC_HGHT];        /* Y Area for serial treatment */
extern	float	_xmin_histo;	/* Min de l'histo cumule */
extern	float	_xmax_histo;	/* Max de l'histo cumule */
extern	int		_nbpas_histo;	/* Nombre de pas de l'histo cumule */

extern  int     _mri_init;
extern  int     _mask_color;
extern  int     _mask_color24b;
extern  int		_threshold_color;
extern  int		_threshold_color24b;
extern	int		_speed_movie;
extern	int		_mip_tag;
extern	int		_movie_cardiowithroi;
extern	int		_movie_cardiocolorroi;
extern	int		_cardio_save_zoom;

/*  Variable globale pour recalage */
extern	double	_sigma_robust;

extern	char	_path_paravision[FILE_LEN]; /* Repertoire pour lecture fichier paravision */

extern	int		_prg_style;		/* Mode de programmation des atm*/
extern	int		_automatic_mode;	/* Mode fctmt de l'automate	*/
extern	FILE   *_fp_auto;		/*  Automatic file ptr.		*/
extern	FILE   *_fp_learn;		/*  Learning  file ptr.		*/
extern  FILE   *_fp_lotus;              /*  Lotus file ptr.             */
extern	char	_atm_autofilename[];	/* Current automatic filename.	*/
extern  char    _lotus_filename[];      /* Current lotus filename.      */

extern char _messages[];  
 /*char		_messages[MESSAGES_LEN]; *//* Contains message value	*/
extern char _warning[];  
/*char		_warning[WARNING _LEN]; *//* Contains message value	*/

extern	char	_CurrentPath[];	/* Last directory path.		*/
extern	int		_full_reserved;
extern  int		_obj_color;		/*  Choix de la couleur pour objets */
extern  int		_txt_color;		/*  Choix de la couleur pour objets */
extern	int		_tala_color;	/* Couleur de la grille de Talairach. */
extern  int		_gomme_color;	/*  Choix de la couleur pour la gomme */
extern  int		_init_roi;		/*Parametre d'initailisation de la roi
				                 lors de sa creation */


extern	int	_read_format_3d; /* Type de lecture pour garantir 
							une representation meca et radio */
extern	unsigned long		_dialog_messages_status;

extern struct	stack  *_imagix_objects;/* Contains imagix object on graphics*/ 
extern struct	stack  *_imagix_objects_3d;/*Contains imagix object3d on */ 
extern struct	stack  *_imagix_texts;  /* Contains imagix texts  on graphics*/ 
extern struct	stack  *_imagix_oft;	/* Contains open files in atm	*/
extern struct	stack  *_loop_cmd__pos;	/* Contains pos on loop cmd in
				           atm file.			*/
extern struct   stack  *_imagix_mire; /* Contains imagix mire 3d */

/* --------------------------- All functions --------------------------	*/
/* Global switch callbacks */
/* ----------------------- */
#if defined(LINUX) || defined(SGI) || defined(WIN32)      /* macro pour compenser un manque dans math.h */
#define log2(x) ( (double)( log(x)/log(2) ) )
#endif

#if defined(LINUX) || defined(DEC) || defined(SGI)
#define sgefa sgefa_
#define sgedi sgedi_
#endif

#if defined(WIN32)
//#define stat _stat
#define sleep _sleep
#endif

#if defined (LINUX) || defined (WIN32) 
#define FREAD_LONG(t,n,fp) (fread_long(t,n,fp))
#define FREAD_LONG2(t,n,fp) (fread_long2(t,n,fp))
#define FREAD_SHORT(t,n,fp) (fread_short(t,n,fp))
#define FREAD_SHORT2(t,n,fp) (fread_short2(t,n,fp))
#define FWRITE_LONG(t,n,fp) (fwrite_long(t,n,fp))
#else
#define FREAD_LONG(t,n,fp) (fread (t,sizeof(long),n,fp))
/* #define FREAD_SHORT(t,n,fp) (fread (t,sizeof(short),n,fp)) */
#define FREAD_SHORT(t,n,fp) (fread_short(t,n,fp))
#define FWRITE_LONG(t,n,fp) (fwrite (t,sizeof(long),n,fp))
#endif               

#ifdef DEBUG_MALLOC

extern void *malloc_debug(size_t size, int line, const char *file);
extern void *calloc_debug(size_t nelem, size_t elsize, int line, const char *file);
extern void *realloc_debug(void *ptr, size_t size, int line, const char *file);
extern void free_debug(void *ptr, int line, const char *file);

extern char *strdup_debug(const char *s, int line, const char *file);

// extern int _x_activate;

#define malloc(size) malloc_debug((size), __LINE__, __FILE__)
#define calloc(nelem, elsize) calloc_debug((nelem), (elsize), __LINE__, __FILE__)
#define realloc(ptr, size) realloc_debug((ptr), (size), __LINE__, __FILE__)
#define free(ptr) free_debug((ptr), __LINE__, __FILE__)

#define strdup(s) strdup_debug((s), __LINE__, __FILE__)

#endif

    
#ifdef __cplusplus
}
#endif


#endif  /* End of imagix.h 2.01 */					

