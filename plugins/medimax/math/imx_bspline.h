/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/**-----------------------------------------------------------------------
***	
***	file:		imx_mtch.h
***
***	project:	Imagix 1.01
***			
***
***	description:    Fichier pour le matching des images 2D
***	
***	
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***
***---------------------------------------------------------------------*/

#ifndef _imx_bspline_h
#define _imx_bspline_h

/*DEFINITION DE LA STRUCTURE base qui definie une structure contenant
**la taille de l'image, le nombre d'element dans la base, la resolution courante
**et des pointeurs sur des tableaux contenant pour chaque indice les info sur la fonction
**equivalente dans la base*/
typedef struct s_base {
int	width,height,nb_func,nbfx,nbfy,resol; /*parametres divers*/
double	*px,*py;  /*valeur des deux parametres associe a cette fonction de la base*/
int     *idx,*idy; /*indice indiquant la position de la fonction dans la matrice de parametres*/
int	*x0,*y0;  /*point en haut a gauche du support la fonction*/
int	*x1,*y1;   /*point en bas a droite du support */
/*ATTENTION, la valeur des fonctions etant la meme pour chaque element de la base, il n'y a qu'un
**pointeur pour toute la base*/
double  *fx,*fy; /*pointeur sur les tableaux contenant la valeur des fonction 1D:									**F(i,j)=fi(i)*fj(j)*/
double  *dfx,*dfy; /*pointeur sur les tableaux contenant les derivees de fi et fj suivant i et j*/

double  *filt; /*pointeur sur le tableau contenant le filtre de passage d'une resolution a l'autre*/
int	tfilt; /*taille du filtre*/
int	(*scal_func)(); /*pointeur sur la fonction d'echelle*/
} base, *ptr_base;


/**************** GESTION DE LA BASE DE FONCTIONS *****************************/
int	Bsplined0(int pas, double **f, double **df);
int	Bsplined1(int pas, double **f, double **df);
int	Bsplined2(int pas, double **f, double **df);
int	Bsplineod1(int pas, double **f, double **df);
int	init_base(int wdth, int hght, int (*scal_func) (/* ??? */));
int	init_base_resol(int wdth, int hght, int (*scal_func) (/* ??? */), int resolution);
int	base_resol_up(double *param, int nb_param);
int	end_base(void);

extern base BASE;

#endif
