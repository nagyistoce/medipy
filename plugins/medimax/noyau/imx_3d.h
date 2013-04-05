/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

//
//  File    : imx_picture3d.h
//  Project : GTK-Imagix
//
//  Description:
//      ...
//
//  Copyright (c) 1993, ULP-IPB Strasbourg
//  All rights are reserved
//
//  >>  Created on ... ..th, .... by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//

#ifndef __3D_H
#define __3D_H

#include "noyau/imx_types.h"
#include "outils/imx_list.h"

#ifdef WIN32
#define DllExport   __declspec( dllexport )  // for windows
#else
#define DllExport  // Linux
#endif


/* ----------------------------- Macros -------------------------------	*/
#define TYPEMRI3D short   /* Indique le type (C) de donne de mri a modifier en meme temps que MAXMRI3D */

#define MAXMRI3D 32767   /* short -> 32767
			    unshort -> 	65535
			    unchar ->	255
			    char ->		128
			    int -> 	4194304	*/
#define MINMRI3D -32768   /* short -> -32768
			     unshort -> 0
			     unchar ->  0
			     char ->		-127
			     int -> 	-4194304	*/

/* Definir le Type pour fonction curve */
#define CURVE_LINE                 1
#define CURVE_MANUAL               2
#define CURVE_PARABOLE_SYM         3
#define CURVE_PARABOLE_QC          4
#define CURVE_POLYNOME	           5
#define CURVE_POLYLIGNE	           6
#define CURVE_DOMCONVEX	           7

#define MAX_CRAYON                 255
#define MAX_CTHETA                 180
#define MAX_CPHI                   180
#define MAX_CZOOM                  180
#define MIN_CRAYON                 0
#define MIN_CTHETA                 0
#define MIN_CPHI                   0
#define MIN_CZOOM                  0

#define CRAYON	                    1
#define CTHETA	                    2
#define CPHI	                    3
#define CZOOM	                    4

#define HORIZONTAL_SPACING         8
#define VERTICAL_SPACING           8
#define TOP_OFFSET                 (TDimension) 30
#define SCALE_HIGHLIGHT_THICKNESS  (TDimension) 2
#define SCALE_WIDTH                400
#define SCALE_HEIGHT               25

/* macro pour fMRI */
#define FMRI_MAX_ETATS     16
#define FMRI_SEUIL_DEFAUT  (20.0/100.0)

/* -- constantes et macros ------------------------------------------------- */
#define ERR_TYPING    1

#define REDUC_JUMP    0
#define REDUC_MEAN    1
#define REDUC_MEDIAN  2
#define REDUC_MIP     3

/* ----------------------------- New type ----------------------------- */
typedef struct Xvector_of_3D_points
{
  long x;
  long y;
  long z;
}XPoints_3d;

typedef struct Xitabvector_of_3D_points
{
  XPoints_3d *tab;
  long  nbpts;
}TabXPoints_3d;

typedef struct s_dpoint3d {
  double x,y,z;
} dpoint3d, *ptr_dpoint3d;

typedef struct s_dplan
{
  double a, b, c, d;
} dplan;

/******************************************************************************
 **  Definition du type t_curv pour la courbe curviligne
 **  Ce type regroupe :
 **  - les parametres de la courbe
 **  - le tableau des rayons (phi, th -> r)
 ******************************************************************************/
typedef enum
  {
    CURV_NONE,
    CURV_CURV,
    CURV_PLAN,
    CURV_DISPLAYED
  } TCurvState;


typedef struct s_curv
{
  TCurvState    state;                /* champs initialises ?        */

  double *vect;                /* coefficients de la courbe   */
  int     nb_params;           /* nombre de coefficients      */

  double theta_min, phi_min;   /* valeurs de debut des angles */
  double theta_max, phi_max;   /* valeurs de fin des angles   */
  double nb_theta,  nb_phi;    /* nombre de points            */

  float rayon;		     /* donnees de l'interface      */
  float oldrayon;
  float theta;
  float phi;
  float zoom;

  int   xg, yg, zg;            /* centre de gravite           */

  double *r;                   /* tableau des rayons          */
} t_curv;

/**************************************************
 **  Definition de la structure ref_talairach
 **  ce type regroupe les points de repere dans l'atlas de Talairach
 *****************************************************/
typedef struct s_label
{
  short area;
  short gyrus;

} label;

typedef struct s_labelBox
{ // On travaille en MNI, comme MNI est donn�par des �uations approximatives
  // la taille de la matrice a ��choisie de fa�n large afin d'�iter tout d�ordement
  // on a en gros : [coord max en x + coord max en -x], pareil pr y et z.
  label labelAt[80+80][80+110][90+60]; 
  
  // Tables de correspondance entre un entier (indice) et un label.
  char* areaLabelTable[29]; 
  char* gyriLabelTable[47];

} LabelBox;

typedef struct s_references_talairach
{

  /* Coordonnees des points de reference. */

  dpoint3d AC;	/* Anterior Commissure */
  dpoint3d PC;	/* Posterior Commissure */
  dpoint3d AP;	/* Anterior Point */
  dpoint3d PP;	/* Posterior Point */
  dpoint3d RM;	/* Right Medial */
  dpoint3d LM;	/* Left Medial */
  dpoint3d SP;	/* Superior Point */
  dpoint3d IP;	/* Inferior Point */
  dpoint3d IC;	/* Inferior Cerebellum (ne sert a rien) */


		/*Vecteurs de base de l'atlas de Talairach. */

  dpoint3d Tu;
  dpoint3d Tv;
  dpoint3d Tw;

  /* Coordonnees des points dans le repere (AC;Tu,Tv,Tw) */

  double x_AP, x_PC, x_PP;
  double z_IP, z_SP;
  double y_RM, y_LM;


  /* Etat des connaissances. */

  int etat;
  /* Au depart, lorsqu'aucune donnee n'a encore
     ete acquise, etat = 0.
     Quand on a tous les points et que
     le repere a ete determine, etat = 1. */


  /* Plan inter-hemispherique. */

  int p_sagi_ok;
  /* Si p_sagi_ok = 1, on connait le plan,
     sinon cela signifie qu'on ne connait pas le plan. */
  dplan p_sagi;
  /* Equation du plan. */



  /* Grille de Talairach. */

  struct stack * grid_id;
  /* Si la grille n'existe pas, grid_on = NULL,
     sinon, grid_on = l'adresse de l'objet definissant la grille.
     ! la grille existe n'entraine pas forcement qu'elle est tracee. */
  int grid_on;
  /*Si grid_on = 1, la grille est affichee,
    si grid_on = 0, elle n'est pas affichee. */

  LabelBox* labelbox;
  /* Contient les donn�s de labellisation qui sont charg�s �partir de fichiers
     �l'init du mode tala. */    

} ref_talairach;

/**************************************************
 **  Definition de la structure DTI
 **  ce type regroupe les informations caract�istiques 
 **  d'une s�uence d'images DTI
 *****************************************************/

typedef struct s_dti
{

int is_dti;		/* 1 si c'est une image appartenant a une sequence d'image dti,
                           2 si c'est un tenseur de diffusion,
                           3 si c'est une image brute recalee,
                           4 si c'est une image de diffusion calcul�e avec le voisinage
                           et  0 sinon */
int nb_directions;	/* Nombre de directions d'acquisition */
double gx;		/* gradient suivant x */
double gy;		/* gradient suivant y */
double gz;		/* gradient suivant z */
double b;		/* valeur du coefficient de diffusion b */
 
} t_dti;

/**************************************************
 **  Definition de la structure pta
 **  ce type regroupe les informations des points
 **  anatomiques d'une image 3d
 *****************************************************/
 
typedef struct s_pta
{
	TDimension coord[3];	// Coordonnées du point anatomique
	unsigned char ok;		// Indique si ce point anatomique a été désigné par l'utilisateur (il faut donc l'afficher)
	unsigned char type;		// Type : nez, oreille G. ou oreille D.
} t_pta;

/**************************************************
 **  Definition de la structure grphic3d
 **  ce type regroupe les informations des images 3D
 *****************************************************/

struct s_graphic3d;

typedef struct s_graphic3d
{
	// Private information:
	// For graphics...
	//  int        nr;               // ne sert pas (WIEST) From 0 to MAX_PICTURES_3D-1
	int        pos;              // Position on disp 1 to MAX_PICTURES_3D
	int        bitppixel;        // Bit per pixel, Assumed 32
	#ifdef  __GTK_GUI
	t_CMapInfo cmapinfo;         // Needed data for colormap
	#endif
	int        mask_active;      // TRUE si le masque est actif
	unsigned char alpha_channel; // % of transparency ( 100 -> full transp.)
	int        bar;              // indicates if palette's bar is shown
	
	// Plan de coupe visualise en multi planar
	int        x3dc, y3dc, z3dc;  // Position courante des plans de coupe
	int        zoom;		// Niveau de zoom de l'affichage
	int        zoom_x;		// Offset d'affichage x
	int        zoom_y;		//   "		      y
	int        zoom_z;		//   "		      z
	int        cutoff_min;       // Valeur affichage min
	int        cutoff_max;       // Valeur affichage max
	int img_coro,img_sagi;       // Position on display multi planar
	int img_tran,img_curv;       // Position on display multi planar
	TMaskType  mask_type;        // Type of mask
	
	// Public information:
	struct s_graphic3d *mask;
	TImageType img_type;              // Type of picture
	
	char       type;		      // b means BRUKER, v VISILOG.
	
	TYPEMRI3D  ***mri;  	      // Picture data.
	TDimension  mri_width, mri_height, mri_depth; // Width, Height and Depth de mri
	TDimension  width, height, depth; // Width, Height and Depth de l'image dans mri
	long       max_pixel;	      // Max_pixel.
	long       min_pixel;	      // Min_pixel.
	float      icomp;		      // rcoeff=2.**icomp
	float      rcoeff;  	      // coefficient mult pour valeur reelle
	float      dx, dy, dz;	      // Resolution 3D, taille des
						// voxels en mm
	
	Patient    patient; 	      // Indentification du patient
	char       filename[256];	      // Picture filename source if any.
	
	// Donnees relatives au referenciel de Talairach.
	ref_talairach  tala;
	
	// Donnees relatives a la courbe curviligne
	t_curv	   curv;
	
	// Donnees relatives aux images DTI
		t_dti dti;
	
	// annotations: liste d'informations.
	// Chaque information est attach� �une position donn�de l'image.
	// accesible uniquement en c++ en utilisant l'accesseur  get_annotations(...);
	//!!! Contient également les points anatomiques (nom commençant par "pta_" : à ignorer)
	void *annotations;
	
	// Gestionnaire de propri��
	void* properties;
	
	unsigned char hide_pta;		// Booléen : indique si l'on cache ou non les points anatomiques
	
} grphic3d, *ptr_grphic3d;



/**************************************************
 **  Definition de la structure imagecurve
 **  ce type regroupe  for image curve
 *****************************************************/
typedef struct s_imagecurve { /* For image curve */
  grphic3d         *image;
  float            rayon;
  float            oldrayon;
  float            theta;
  float            phi;
  float            zoom;
  double           *tabrayon;
  double           *tabcteta,*tabsteta,*tabcphi,*tabsphi;
  int              xg,yg,zg,nteta,nphi;
  double           teta1,teta2,phi1,phi2;
} imagecurve;


/*******************************************************
 **  Definition de la structure T_drawing_3d
 **  qui enregistre les voxels marques dans
 **  la visu3d
 ********************************************************/
typedef struct s_drawing_3d {
  unsigned int nb_obj; /* Nombre d'objets dans la derniere scene 3D */
  t_ptrlist *tab_obj;   /* Tableau des dessins sur les objets */
} T_drawing_3d;


#ifdef __cplusplus
extern "C" {
#endif
  extern DllExport grphic3d *ptr_img_3d(int nr);

  extern DllExport grphic3d *cr_grphic3d(grphic3d *temp);
  extern DllExport grphic3d *cr_grphic3d_modif(UINT width, UINT height, UINT depth, float icomp, float rcoeff, UINT maxpixel);
  extern DllExport void free_mri_3d(TYPEMRI3D ***mri);
  extern DllExport void free_grphic3d(grphic3d *temp);
  extern grphic3d *cr_grphic3d_array(UINT num, UINT width, UINT height, UINT depth);
  extern grphic3d **cr_ptr_grphic3d_array(UINT num, UINT width, UINT height, UINT depth);
  extern void free_grphic3d_array(grphic3d *array, int n);
  extern void free_ptr_grphic3d_array(grphic3d **array, int n);
  extern DllExport TYPEMRI3D ***cr_mri_3d(UINT width, UINT height, UINT depth);
  extern grphic3d *cr_mask_3d_ptr(grphic3d *img);
  extern grphic3d *init_mask_3d_ptr(grphic3d *img);
  extern int ptr_has_mask_3d(int pos);
  extern int ptr_has_mask_3d_p(grphic3d *img);
  extern void    free_mask_3d_p(grphic3d *mask) ;
  extern DllExport grphic3d *ptr_mask_3d(int pos3d);
  extern DllExport grphic3d *ptr_mask_3d_p(grphic3d *img);
  extern int resize_grphic3d_buffer(grphic3d *img, TDimension new_width, TDimension new_height, TDimension new_depth);
  extern int copie_image_param_adapt_buffer(grphic3d *imdeb, grphic3d *imres);

  extern int imx_inimaxminpixel_3d(int im_1);
  extern int imx_inimaxpixel_3d_p (grphic3d *im1);
  extern DllExport int imx_inimaxminpixel_3d_p (grphic3d *im1);
  extern int imx_iniminpixel_3d_p (grphic3d *im1);
  extern int imx_iniparaimg_3d_p (grphic3d *im1);
  extern int imx_iniparaimg_3d(int im_1);
  extern void imx_iniimg3d(int pos, int *img_coro, int *img_sagi, int *img_tran, int *img_curv);

  extern int  imx_copie_dti_params_3d_p(grphic3d *im1, grphic3d *im2);
  extern int imx_copie_param_3d(int im_1, int im_2) ;
  extern DllExport int imx_copie_param_3d_p(grphic3d *im1, grphic3d *im2) ;
  extern int imx_copie_visual_params_3d_p(grphic3d *im1, grphic3d *im2);
  extern DllExport int imx_copie_apparance_params_3d_p(grphic3d *im1, grphic3d *im2);
  extern int imx_copie_patient_params_3d_p(grphic3d *im1, grphic3d *im2);
  extern int imx_copie_dimension_params_3d_p(grphic3d *im1, grphic3d *im2);
  extern int imx_copie_image_params_3d_p(grphic3d *im1, grphic3d *im2);
  extern int imx_copie_dimension_param_3d_p(grphic3d *imsrc, grphic3d *imdst);
  extern int imx_iniparam_3d(int im_1, long int max_pixel, long int min_pixel,
			     float icomp, float rcoeff);

  extern DllExport int imx_brukermax_3d(float amax, float amin, grphic3d *graphic) ;
  extern void  imx_brukermax_update_3d(grphic3d *im);

  extern void IMXPict3d_CleanPicture( void );
  extern void IMXPict3d_ResetPicture( int pos3d );

  extern void SetImg3d(int pos, grphic3d *graphic);
  extern int  imx_pos_2dto3d(int pos) ;
  extern int  imx_ecran2dto3d(int pos2d, int x2d, int y2d, int *x3d, int *y3d, int *z3d) ;
  extern int  imx_ecran3dto2d(int pos2d, int *x2d, int *y2d, int x3d, int y3d, int z3d);


#ifdef __cplusplus
}
#endif


/* ------------------------- Global variables -------------------------	*/
extern grphic3d *_imx_img3d;
extern grphic3d _imx_roi3d;

extern float 	 *_x_histo_3d;    /* X Area for histo treatment */
extern float 	 *_y_histo_3d;    /* Y Area for histo treatment */
extern float    _xmin_histo_3d;  /* Min de l'histo cumule */
extern float    _xmax_histo_3d;  /* Max de l'histo cumule */
extern int      _nbpas_histo_3d; /* Nombre de pas de l'histo cumule */

extern T_drawing_3d _drawing_3d;

#endif  // #ifndef __3D_H
