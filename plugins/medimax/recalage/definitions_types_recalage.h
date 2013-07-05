/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*******************************************************************************/
/*!	   \file:		definitions_types_recalage.h
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

#ifndef __DEFINITIONS_TYPES_RECALAGE_H__
#define __DEFINITIONS_TYPES_RECALAGE_H__

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_3d.h"
#include "noyau/imx_2d.h"

//---TYPES LIES AUX TRANSFORMATIONS---//

/***************************** DEFINITION DU TYPE vector **********************/
typedef struct s_vector3d {
 float x,y,z;
} vector3d, *ptr_vector3d;

/***************************** DEFINITION DU TYPE dvector **********************/
typedef struct s_dvector3d {
 double x,y,z;
} dvector3d, *ptr_dvector3d;

/***************************** DEFINITION DU TYPE field ***********************/
typedef struct s_field3d {
  int		width,height,depth;
  float dx,dy,dz;              /*! dans le cas d'une transfo anisotropique*/
  vector3d	***raw;
} field3d, *ptr_field3d;


/***************************** DEFINITION DU TYPE transf3d ***********************/
typedef enum {CHAMPVOXEL3D,CHAMP3D,RIGID3D,RIGIDZOOM3D,AFFINE3D,BSPLINE3D,RIGIDGLOBALZOOM3D,AFFINEDECOUPLE3D} TRANSFO3D;
typedef struct s_transf3d {
  TDimension	width,height,depth;
  float dx,dy,dz;              /*! dans le cas d'une transfo anisotropique*/
  TRANSFO3D typetrans;
  int	nb_param;
  double *param;
  int	resol,degre; /*parametres utilises uniquement pour les ondelettes*/
  char  *trans_ini; /*nom de fichier de la transfo initiale si existe*/
  /*vector3d ***ch;*/	    /*pointeur sur les element du champ 3d si typetrans=CHAMP3D*/
  vector3d *ch;	/* fix by R Nabet */
} transf3d, *ptr_transf3d;

//prototype d'une fonction d'interpolation
typedef int (*InterpolationFct)(grphic3d *, field3d *, grphic3d *);

//---TYPES LIES A LA GESTION D'ERREUR---//

typedef enum
{
  NO_ERR,                      //pas d'erreur
  NB_BINS_NUL,                 //nb de bins nul dans un calcul d'histogramme
  ERREUR_ALLOC_MEMOIRE,        //allocation memoire qui a echoue
  PB_TAILLES_IMAGES,           //taille d'images incompatibles dans le recalage
  OVERLAP_NUL,                 //zone d'overlap des images ref et reca nulle dans le calcul de distance
  PB_TAILLE_CHAMP,             //taille du champ incompatible avec taille d'image a transformer
  PB_PARAMETRES,               //trop peu de parametres dans la methode simplexe
  TRANSFO_INADAPTEE            //la transformation passee en parametre est indaptee pour ce qu'on veut faire
} ERREUR_RECALAGE;


//---TYPES LIES AUX DISTANCES---//

//donnees liees a l'interpolation volume partiel
typedef struct s_VP_histos {
 int nb_voxels_recouvrement;
 int nb_cl1;
 int nb_cl2;
 int max1;
 int max2;
 double * histo_marginal1;
 double * histo_marginal2;
 double ** histo_conjoint;
 double **** histos_gradient;
 double entropie_conjointe;
 double entropie1;
 double entropie2;
 bool minimisation_par_gradient;
 } VP_histos, *ptr_VP_histos;

//definition des type de distance basees sur l'information mutuelle:
//ENTROPIE_CONJOINTE: l'entropie conjointe
//IM: l'information mutuelle
//IMNS: l'information mutuelle normalisee de Studholme
//IMNM1,IMNM2: les 2 informations mutuelles normalisees de Maes
typedef enum {ENTROPIE_CONJOINTE,IM,IMNS,IMNM1,IMNM2} IM3DTYPE;

//parametres a passer a une fonction distance
typedef struct
{

 grphic3d *imreca;
 grphic3d *imref;
 int nbclRef;
 int nbclReca;
 int maxRef;
 int maxReca;
 field3d *champ;
 VP_histos * donnees_VP;
 field3d *points_calcules;
 double ***weights;
 ERREUR_RECALAGE *err;
 transf3d *transfres;

} distance_param, *ptr_distance_param;

// definition du prototype d'une fonction distance
typedef double (*dist_func_t)(ptr_distance_param dist_param);
//et du prototype du gradient d'une fonction distance
typedef int (*gradient_func_t)(grphic3d *I2, grphic3d *Im, int nb_param, double *param, field3d *gradI2, double *grad, VP_histos *donnees_VP);

//---TYPES LIES A LA MINIMISATION---//

//prototype d'une fonction contrainte
typedef double (*contrainte_func_t)(int wdth, int hght, int dpth, int nb, double *param, field3d *champ_fin);

//parametres a passer a une methode de minimisation
typedef struct
{

 grphic3d *imref;
 grphic3d *imreca;
 grphic3d *imres;
 field3d *champ_ini;
 field3d *champ_fin;
 InterpolationFct inter;
 dist_func_t dist;
 contrainte_func_t contrainte;
 double *min_param;
 double *max_param;
 double *prec_param;
 ptr_distance_param dist_param;
 transf3d *transfres;
 ERREUR_RECALAGE *err;

} minimisation_param, *ptr_minimisation_param;

//protoype d'une fonction de minimisation
typedef double (*min_func_t)(ptr_minimisation_param);


//---TYPES LIES AU MATCHING---//

typedef int (*transf_func_t)(int nb_param, const double *param, field3d *champ,grphic3d* imref, grphic3d* imreca);

typedef int (*scal_func_t)(int pas, double **f, double **df);


/*DEFINITION DE LA STRUCTURE base qui definie une structure contenant
**la taille de l'image, le nombre d'element dans la base, la resolution courante
**et des pointeurs sur des tableaux contenant pour chaque indice les info sur la fonction
**equivalente dans la base*/
typedef struct s_base3d {
int	width,height,depth,nb_func,nbfx,nbfy,nbfz,resol; /*parametres divers*/
double	*px,*py,*pz;  /*valeur des deux parametres associe a cette fonction de la base*/
int     *idx,*idy,*idz; /*indice indiquant la position de la fonction dans la matrice de parametres*/
int	*x0,*y0,*z0;  /*point en haut a gauche du support la fonction*/
int	*x1,*y1,*z1;   /*point en bas a droite du support */
/*ATTENTION, la valeur des fonctions etant la meme pour chaque element de la base, il n'y a qu'un
**pointeur pour toute la base*/
double  *fx,*fy,*fz; /*pointeur sur les tableaux contenant la valeur des fonction 1D:									**F(i,j)=fi(i)*fj(j)*/
double  *dfx,*dfy,*dfz; /*pointeur sur les tableaux contenant les derivees de fi et fj suivant i et j*/

double  *filt; /*pointeur sur le tableau contenant le filtre de passage d'une resolution a l'autre*/
int	tfilt; /*taille du filtre*/
scal_func_t scal_func; /*pointeur sur la fonction d'echelle*/
} base3d, *ptr_base3d;

//intervalles de recherche possible pour un recalage
typedef enum {SMALL,MEDIUM,LARGE} eResearchInterv;
//precisions possibles pour un recalage
typedef enum {NORMAL,PRECIS,TRES_PRECIS,PRECISION_MAX} eMatchPrecision;



#endif /*__DEFINITIONS_TYPES_RECALAGE_H__*/
