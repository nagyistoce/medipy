/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*******************************************************************************/
/*!	   \file:		minimisation_3d.c
***
***		project:	Imagix 2.01
***
***
***		\brief description: fonction de minimisation 
***
***
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***
********************************************************************************/

#include <config.h>
#include <stdio.h>
#include <stdlib.h> 	/* Needed when you use GET_FLOAT	*/
#include <math.h> 	/* Needed when you use math function	*/
#include <limits.h>
#include <float.h>
#include <time.h>

#include "noyau/imagix.h"                                                               
#include "noyau/gtkimagix.h"
#include "noyau/imx_3d.h"
#include "noyau/imx_lang.h"
#include "noyau/io/imx_head.h"
#include "recalage/chps_3d.h"
#include "recalage/distance_3d.h"
#include "recalage/mtch_3d.h"
#include "math/imx_matrix.h"
#include "outils/imx_sort.h"
#include "noyau/gui/imx_picture3d.h"
#include "noyau/mani_3d.h"
#include "traitement/norm_3d.h"
#include "recalage/mtch_robust_3d.h"
#include "outils/imx_misc.h"

#include "math/imx_bspline.h"	/* pour les fonctions d'interpolation */

#include "noyau/io/imx_log.h"
#include "math/oper_3d.h"
#include "traitement/trai_3d.h"

#include "recalage/gradient_distance_3d.h"
#include "recalage/transformations_3d.h"
#include "recalage/erreur_recalage.h"

#include "recalage/minimisation_3d.h"

#ifndef NB_ITER_MAX_LINEMIN_3D
#define NB_ITER_MAX_LINEMIN_3D 30 /*nombre d'iterations max pour la minimisation lineaire*/
#endif

/********* VARIABLES DEFINE*********/

/* constantes pour minimisation powel */
#ifndef POWELL_TINYRES
#define POWELL_TINYRES 0.5
#endif
#ifndef POWELL_TINY
#define POWELL_TINY 0.001
#endif

/*linemin*/
#define	NB_PAS_LINEMIN_3D	5      /*nb pas de decoupage pour minimisation linemin*/
/*recalage icm*/
#define	NB_ITER_MAX_ICM_3D	30 /*nb interations max pour minimisation icm*/
#define NB_REDUC_MAX_ICM_3D	3	/*nb de reduction maximum dans l'icm*/
#define	P_REDUC_3D		0.2	/*pourcentage mini de variation d'energie qui entraine*/
					/* une reduction de l'intervalle de recherche*/
#define C_REDUC_3D		0.5     /*coefficient de reduction de l'intervalle de recherche*/
/*recalage simplex*/
//#define NB_ITER_MAX_SIMPLEX_3D	500  	/*nb interations max pour minimisation simplex*/
#define NB_ITER_MAX_SIMPLEX_3D	200
#define P_STOP_SIMPLEX_3D	1	/*pourcentage mini de taille du simplex pour arret*/
/*recalage desc grad*/
//#define NB_ITER_MAX_DESC_GRAD_3D 50  	/*nb interations max pour minimisation descente de gradient*/
#define NB_ITER_MAX_DESC_GRAD_3D 50
/*recalage grad conj*/
//#define NB_ITER_MAX_GRAD_CONJ_3D 50  	/*nb interations max pour minimisation gradient conjugue*/
#define NB_ITER_MAX_GRAD_CONJ_3D 5
/*recalage desc grad mod*/
//#define NB_ITER_MAX_DESC_GRAD_MOD_3D 50  	/*nb interations max pour minimisation desc gradient modifiee*/
#define NB_ITER_MAX_DESC_GRAD_MOD_3D 5
/*recalage quasi newton mod*/
//#define NB_ITER_MAX_QUASI_NEWTON_MOD_3D 50  	/*nb interations max pour minimisation desc gradient modifiee*/
#define NB_ITER_MAX_QUASI_NEWTON_MOD_3D 5
/*recalage powell 2*/
#define NB_ITER_MAX_POWELL_3D 30  	/*nb interations max pour minimisation desc gradient modifiee*/

#ifndef NB_MINIMISATION_MAX_LINEMIN_3D
#define NB_MINIMISATION_MAX_LINEMIN_3D 1
#endif

#ifndef CGOLD
#define CGOLD 0.3819660
#endif
#ifndef GOLD
#define GOLD 1.618034
#endif
#ifndef ZEPS
#define ZEPS 1.0e-10 /*definition de epsilon*/
#endif
#ifndef EEPS
#define EEPS 1.0e-2
#endif
#ifndef MNBRAK_TINY
#define MNBRAK_TINY 1.0e-20
#endif
#ifndef GLIMIT
#define GLIMIT 100.0
#endif
#ifndef SGN
#define SGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#endif

#ifndef BRAIN_RADIUS
#define BRAIN_RADIUS 100.0
#endif

#ifndef EPS
#define EPS (1.0e-6)
#endif

/* seuil en-dessous duquel on arrete l'iteration de la minimisation locale :
un seuil plus eleve permet un recalage plus rapide, mais sans doute au prix
d'une moindre qualite du recalage */
#define SEUIL_FIN_LOCAL /*0.02*/0.01/*0.005*/
/* seuil en-dessous duquel on renonce a la minimisation locale, car on considere
que la matrice ne peut pas etre valablement inversee */
#define SEUIL_INVERSION_LOCAL 0.1/*1000000000000000000000000.*/

#define BRENT_GOLD 1.618034
#define BRENT_GLIMIT 100.0
#define BRENT_TINY 1.0e-20
#define BRENT_MAX(a,b) ((a) > (b) ? (a) : (b))
#define BRENT_SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define BRENT_SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

/*******************************************************************************
********************************************************************************
*************************** MINIMISATIONS **************************************
********************************************************************************
*******************************************************************************/
/*! \ingroup     MinimisationFonction  @{ */
/*******************************************************************************
**    fonc_transfo_3d(imref,imreca,imres,the_transf,inter,dist,champ_ini,champ_fin,nb,param)
*/
/*!
**    fonction de plusieurs variables a minimiser
**		\param imref : images ref.
**		\param imreca : images recal
**		\param imres : images resultat (E/S)
**      \param the_transf : pointeur sur la fonction de transformation
**                           int the_transf(nb_param,param,champ)
**      \param inter : pointeur sur la fonction d'interpolation
**                           int inter(imdeb,champ,imres)
**      \param dist : pointeur sur la fonction de distance
**                         double dist(im1,im2)
**		\param contrainte : fonction de calcul de la contrainte (NULL si pas de contrainte)
**      \param champ_ini : champ initial
**                              NULL = pas de champ initial
**      \param champ_fin : champ final resultat de la transfo
**		\param nb : nbre de parametres
**      \param param : resultats donnant la valeur des parametres au
**                         point definissant le minimum    (E/S)
**		\retval : le resultat du calcul de la distance (multiplie par le resultat de la contrainte si != NULL)
*******************************************************************************/
double fonc_transfo_3d(grphic3d *imref, grphic3d *imreca, grphic3d *imres,
			transf_func_t the_transf, InterpolationFct inter,
                        dist_func_t dist, contrainte_func_t contrainte,
                        field3d *champ_ini, field3d *champ_fin,
                        int nb, double *param, VP_histos *donnees_VP, field3d *points_calcules, double ***weights, ERREUR_RECALAGE *err)
{
 double y,c;
 int wdth,hght,dpth;
 ptr_distance_param dist_param=CALLOC(1, distance_param);
 transf3d *tmpTransf=NULL;

 wdth=imref->width;hght=imref->height;dpth=imref->depth;

 tmpTransf=cr_transf3d_p(imref, NULL);
 tmpTransf->nb_param=nb;
 tmpTransf->param=param;

 //on assume que si points_calcules!=NULL alors il s'agit d'un calcul de distance stochastique
 //un peu lourd : le calcul stochastique est independant de l'erreur erreur_IM_VP_stochastique_3d et n'est qu'une faon differente de calculer valable pour toutes les erreurs
 //ici on fait un cas a part pour ne pas rajouter une option dans le matching
/* if (points_calcules!=NULL)
 {
  if (the_transf==rigid_to_field_3d) rigid_to_field_champ_3d(nb,param,points_calcules,champ_fin,imref,imreca);
  else if (the_transf==rigidz_to_field_3d) rigidz_to_field_champ_3d(nb,param,points_calcules,champ_fin,imref,imreca);
  else if (the_transf==affine_to_field_3d) affine_to_field_champ_3d(nb,param,points_calcules,champ_fin,imref,imreca);
  else if (the_transf==rigid_global_zoom_to_field_3d) rigid_global_zoom_to_field_champ_3d(nb,param,points_calcules,champ_fin,imref,imreca);
  else if (the_transf==affine_decouple_to_field_3d) affine_decouple_to_field_champ_3d(nb,param,points_calcules,champ_fin,imref,imreca);
  else rigid_to_field_champ_3d(nb,param,points_calcules,champ_fin,imref,imreca);
 }
 else
 {
  //calcul de la transformation
  (*the_transf)(nb,param,champ_fin,imref, imreca);
 }*/
 if (the_transf==rigid_to_field_3d)
 {
  tmpTransf->typetrans=RIGID3D;
  transf_geom_to_field3d(tmpTransf, champ_fin, imref, imreca, points_calcules);
 }
 else if (the_transf==rigid_global_zoom_to_field_3d)
 {
  tmpTransf->typetrans=RIGIDGLOBALZOOM3D;
  transf_geom_to_field3d(tmpTransf, champ_fin, imref, imreca, points_calcules);
 }
 else if (the_transf==rigidz_to_field_3d)
 {
  tmpTransf->typetrans=RIGIDZOOM3D;
  transf_geom_to_field3d(tmpTransf, champ_fin, imref, imreca, points_calcules);
 }
 else if (the_transf==affine_decouple_to_field_3d)
 {
  tmpTransf->typetrans=AFFINEDECOUPLE3D;
  transf_geom_to_field3d(tmpTransf, champ_fin, imref, imreca, points_calcules);
 }
 else if (the_transf==affine_to_field_3d)
 {
  tmpTransf->typetrans=AFFINE3D;
  transf_geom_to_field3d(tmpTransf, champ_fin, imref, imreca, points_calcules);
 }
 else 
 {
  (*the_transf)(nb,param,champ_fin,imref, imreca);
 }

 //au champ final on rajoute le champ initial
 if (champ_ini!=NULL)
 {
  add_field_3d(champ_fin,champ_ini,champ_fin);
 }


 //calcul de l'erreur
 //cas general:
 if (inter!=inter_VP_3d)
 {
  //on cherche le champ qui donne la transformation inverse de imreca dans imref
  //donc le champ qui donne la transformation de imref vers imreca 
  (*inter)(imreca,champ_fin,imres);

  //on veut minimiser l'erreur entre imres et imref
  dist_param->imreca=imres; dist_param->imref=imref;
  dist_param->champ=champ_fin;
  dist_param->donnees_VP=donnees_VP;
  dist_param->points_calcules=points_calcules;
  dist_param->weights=weights;
  dist_param->err=err;
//  y=(*dist)(imres,imref,champ_fin,donnees_VP, points_calcules,weights, err);
  y=(*dist)(dist_param);
 }
 //cas de l'interpolation volume partiel:
 else
 {
  dist_param->imreca=imreca; dist_param->imref=imref;
  dist_param->champ=champ_fin;
  dist_param->donnees_VP=donnees_VP;
  dist_param->points_calcules=points_calcules;
  dist_param->weights=weights;
  dist_param->err=err;
//  y=(*dist)(imreca,imref,champ_fin,donnees_VP, points_calcules,weights, err);
  y=(*dist)(dist_param);
 }

 if (contrainte!=NULL)
   {
   c=(*contrainte)(wdth,hght,dpth,nb,param,champ_fin);
   y=y*c;
   }

 FREE(dist_param);
 tmpTransf->param=NULL; free_transf3d(tmpTransf);

 return(y);  
}
/*! @} */


double energie_transfo_3d(ptr_minimisation_param min_parameters)
{
 double y,c;
 int wdth,hght,dpth;
 ptr_distance_param dist_param=min_parameters->dist_param;
 grphic3d *imref=min_parameters->imref;
 grphic3d *imreca=min_parameters->imreca;
 grphic3d *imres=min_parameters->imres;
 field3d *champ_fin=min_parameters->champ_fin;
 field3d *champ_ini=min_parameters->champ_ini;
 field3d *points_calcules=dist_param->points_calcules;
 transf3d *transfo=min_parameters->transfres;
 dist_func_t dist=min_parameters->dist;
 InterpolationFct inter=min_parameters->inter;
 contrainte_func_t contrainte=min_parameters->contrainte;

 dist_param->transfres=min_parameters->transfres;

 wdth=imref->width;hght=imref->height;dpth=imref->depth;

 //calcul du champ de transformation
 if (points_calcules) //on fait un cas a part parce que tous les cas ne sont pas geres avec le calcul sur grille stochastique
 {
  if ((transfo->typetrans==RIGID3D)||(transfo->typetrans==RIGIDZOOM3D)||(transfo->typetrans==AFFINE3D)||(transfo->typetrans==RIGIDGLOBALZOOM3D)||(transfo->typetrans==AFFINEDECOUPLE3D))
  {
   transf_geom_to_field3d(transfo, champ_fin, imref, imreca, points_calcules);
  }
  else
  { fprintf (stderr, "le cas du recalage sur grille stochastique avec ce type de transformation n'est pas gere pour l'instant!!!"); return 1; }
 }
 else if (transfo->typetrans==BSPLINE3D) // dans ce cas precis les variables globales c'est penible:
                                         //on ne peut pas appeller transf_to_field_3d_noalloc
                                         //car elle reinitialise la base de BSpline qui est en global
 {
  base_to_field_3d(transfo->nb_param,transfo->param,champ_fin,imref, imreca);
 }
 else 
 {
  transf_to_field_3d_noalloc(transfo, champ_fin, imref, imreca);
 }

 //au champ final on rajoute le champ initial
 if (champ_ini!=NULL)
 {
  add_field_3d(champ_fin,champ_ini,champ_fin);
 }


 //calcul de l'erreur
 //cas general:
 if (!dist_utilisant_VP(dist))
 {
  //on cherche le champ qui donne la transformation inverse de imreca dans imref
  //donc le champ qui donne la transformation de imref vers imreca
  (*inter)(imreca,champ_fin,imres);

	if (imreca->mask)
	(*inter)(imreca->mask,champ_fin,imres->mask);
	

  y=(*dist)(dist_param);
 }
 //cas de l'interpolation volume partiel:
 else
 {
  y=(*dist)(dist_param);
 }

 if (contrainte!=NULL)
   {
   c=(*contrainte)(wdth,hght,dpth,transfo->nb_param,transfo->param,champ_fin);
   y=y*c;
   }

 return(y);
}

double eval_energie_direction(ptr_minimisation_param min_parameters, double *param_tmp, double *direction, double coef)
{
 int i;
 double Energie;
 int nb_param=min_parameters->transfres->nb_param;
 double *param=min_parameters->transfres->param;

 for (i=0;i<nb_param;i++) param_tmp[i]=param[i]+coef*direction[i];
 min_parameters->transfres->param=param_tmp;
 Energie=energie_transfo_3d(min_parameters);
 min_parameters->transfres->param=param;

 return Energie;
}


void mn_brak (double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, ptr_minimisation_param min_parameters, double *direction)  {

	double ulim, u, r, q, fu, dum;
 double *param_tmp=NULL;
 int nb_param=min_parameters->transfres->nb_param;
 ERREUR_RECALAGE *err=min_parameters->err;

 ///////allocations memoire///////
 param_tmp=CALLOC(nb_param, double);
 if ((!param_tmp))
 {
  RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans mnbrak\n");
 }

//	*fa = (*func) (*ax, OptType);
//	*fb = (*func) (*bx, OptType);

  *fa=eval_energie_direction(min_parameters, param_tmp, direction, *ax);
  *fb=eval_energie_direction(min_parameters, param_tmp, direction, *bx);

	if (*fb > *fa) {
		BRENT_SHFT (dum, *ax, *bx, dum)
		BRENT_SHFT (dum, *fb, *fa, dum)
	}

	*cx = (*bx) + BRENT_GOLD * (*bx - *ax);
//	*fc = (*func)(*cx, OptType);

  *fc=eval_energie_direction(min_parameters, param_tmp, direction, *cx);

	while (*fb > *fc)  {
		r = (*bx - *ax) * (*fb - *fc);
		q = (*bx - *cx) * (*fb - *fa);
		u = (*bx) - ((*bx - *cx) * q - (*bx - *ax) * r) /
			(2.0 * BRENT_SIGN (BRENT_MAX (fabs (q - r), BRENT_TINY), q - r));

		// u = minimum of the parabole going through ax, bx and cx
		ulim = (*bx) + BRENT_GLIMIT * (*cx - *bx);

		// Limit BRENT_MAX of bracket ...
		if ((*bx - u) * (u - *cx) > 0.0)  {
		  // first min. of the parabole between bx and cx
//			fu = (*func) (u, OptType);
			fu=eval_energie_direction(min_parameters, param_tmp, direction, u);

			if (fu < *fc)  {
				// the min is between bx and cx
				*ax = (*bx);
				*bx = u;
				*fa = (*fb);
				*fb = fu;
				return;

			} else if (fu > *fb) {
				// the min is between ax and bx
				*cx = u;
				*fc = fu;
				return;
			}

			// u est mal choisi , on a par ex
			// fa > fb > fu > fc
			// on va voir plus loin ...
			// bad u: we could have for ex. fa > fb > fu > fc
			// we're going further
			u  = (*cx) + BRENT_GOLD * (*cx - *bx);
//			fu = (*func) (u, OptType);
			fu=eval_energie_direction(min_parameters, param_tmp, direction, u);

		} else if ((*cx - u) * (u - ulim) > 0.0) {
		  // cx < u < ulim
			// if we decrease more (fu < fc)
			// we're going further
//			fu = (*func) (u, OptType);
			fu=eval_energie_direction(min_parameters, param_tmp, direction, u);

			if (fu < *fc) {
				BRENT_SHFT (*bx, *cx, u,*cx + BRENT_GOLD * (*cx - *bx))
//				BRENT_SHFT (*fb, *fc, fu, (*func) (u, OptType))
				BRENT_SHFT (*fb, *fc, fu, eval_energie_direction(min_parameters, param_tmp, direction, u))
			}

		} else if ((u - ulim) * (ulim - *cx) >= 0.0) {
			// we stumble to the limit
			u  = ulim;
//			fu = (*func) (u, OptType);
			fu=eval_energie_direction(min_parameters, param_tmp, direction, u);

			printf("\n Erreur Possible dans mnbrak ... boucle infinie");


			// Suggestion .... to be tested
			// BRENT_SHFT(*ax,*bx,*cx,u)
			// BRENT_SHFT(*fa,*fb,*fc,fu)
			// return ;

		} else {
			u  = (*cx) + BRENT_GOLD* (*cx - *bx);
//			fu = (*func) (u, OptType);
			fu=eval_energie_direction(min_parameters, param_tmp, direction, u);
		}


		// we're going further
		BRENT_SHFT (*ax, *bx, *cx, u)
		BRENT_SHFT (*fa, *fb, *fc, fu)
	}

end_func:

  FREE(param_tmp);

}

/*! \ingroup     MinimisationFonction  @{ */
/*******************************************************************************
**     linemin_3d(imref,imreca,imres,champ_ini,champ_fin,the_transf,inter
**                   ,dist,nb_param,min_param,max_param,prec_param,param,direc)
*/
/*!
**     Minimisation unidimensionnelle en utilisant la methode de Brent implementee
**     dans le "numerical recipes in C".
**
**     utilisation:
**     \param       imref : images ref
**     \param       imreca : images reca
**     \param       imres : images resultat (E/S)
**     \param       champ_ini : champ initial
**                              NULL = pas de champ initial
**     \param       champ_fin : champ final resultat de la transfo
**     \param       the_transf : pointeur sur la fonction de transformation
**                           int the_transf(nb_param,param,champ)
**     \param       inter : pointeur sur la fonction d'interpolation
**                           int inter(imdeb,champ,imres)
**     \param       dist : pointeur sur la fonction de distance
**                         double dist(im1,im2)
**     \param       nb_param: nombre de parametre a minimiser
**     \param       min_param: bornes inferieures des parametres
**     \param       max_param: bornes superieures des parametres
**     \param       prec_param: precision demandee pour chaque parametre
**     \param       param: resultats donnant la valeur des parametres au
**                         point definissant le minimum
**	   \param      direc: vecteur donnant la direction de minimisation
**     \retval     la fonction retourne un double contenant la valeur du minimum
**      Les parametres correspondant a ce minimum sont retournes dans param
*******************************************************************************/
double linemin_3d(ptr_minimisation_param min_parameters, double *direc)
{
 //les cf representent les coefs utilise pour se deplacer sur la droite passant par 'param' et ayant pour direction 'direc'
 //param2[i] = param[i] + coef*direc[i];
 double cfmin=0.0,cfmax=0.0, cfx, cfw, cfv, cfu, cfp, cfq, cfr, cfprec, cfm, cfetemp, cfd=0.0, cfopt, cfe=0.0, v1, v2;
 //parametres de tolerance tenant compte de la precision qu'on desire avoir sur le calcul des parametres optimaux 'prec_param'
 //et des erreurs de calculs dus a la precision machine
 double tol1, tol2;
 //energies calculees pour les parametres correspondant aux coefs: Ew -> cfw, EBorneMin -> cfmin ...
 double Ew, Ev, Ex, Eu, EParam, Eopt=DBL_MAX;
 //parametres correspondant aux coefs: paramx -> cfx, param_min -> cfmin ...
 double *paramu=NULL;
 int i, iter;
 int nb_minimisations=0;
 double Emin, Emax;

 int nb_param=min_parameters->transfres->nb_param;
 double *param=min_parameters->transfres->param;
 double *min_param=min_parameters->min_param;
 double *max_param=min_parameters->max_param;
 double *prec_param=min_parameters->prec_param;
 ERREUR_RECALAGE *err=min_parameters->err;

 ///////allocations memoire///////
 paramu=CALLOC(nb_param, double);
 if ((!paramu))
 {
  RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans linemin_3d\n");
 }

 ////////calcul de  cfmin et cfmax pour que param+cf*direc soit toujours compris entre min_param et max_param////////
 ////////calcul de  cfprec qui correspond a la precision a obtenir sur les cfopt pour avoir la precision prec_param sur les params////////

 //recherche de la premiere composante non nulle de la direction
// l=0;
// while (fabs(direc[l])<EPS) {l++;}
// //en cas de direction nulle on ne fait rien
// if (fabs(direc[l])<EPS)
// {
//   EParam=energie_transfo_3d(min_parameters);
//   if (erreur_critique(err)) goto end_func;
//   return EParam;
// }
// //initisalition de cfmin et cfmax
// if (direc[l]>0.0) {cfmin=(min_param[l]-param[l])/direc[l];cfmax=(max_param[l]-param[l])/direc[l];}
// if (direc[l]<0.0) {cfmin=(max_param[l]-param[l])/direc[l];cfmax=(min_param[l]-param[l])/direc[l];}
// cfprec=prec_param[l]/direc[l];
// //calcul de cfmin et cfmax
// for (i=0;i<nb_param;i++)
// {
//  if (direc[i]>0.0)
//  {
//   v1=(min_param[i]-param[i])/direc[i];v2=(max_param[i]-param[i])/direc[i];
//   if (v1<cfmin) cfmin=v1;
//   if (v2>cfmax) cfmax=v2;
//  }
//  if (direc[i]<0.0)
//  {
//   v1=(max_param[i]-param[i])/direc[i];v2=(min_param[i]-param[i])/direc[i];
//   if (v1<cfmin) cfmin=v1;
//   if (v2>cfmax) cfmax=v2;
//  }
//  if ((direc[i]!=0.0)&&(cfprec>(prec_param[i]/direc[i]))) cfprec=prec_param[i]/direc[i];
// }

 //calcul de cfmin, cfmax et cfprec
 cfmin=0.0; cfmax=0.0; cfprec=DBL_MAX;
 for (i=0;i<nb_param;i++)
 {
  if (direc[i]>EPS)
  {
   v1=(min_param[i]-param[i])/direc[i];v2=(max_param[i]-param[i])/direc[i];
   if (v1<cfmin) cfmin=v1;
   if (v2>cfmax) cfmax=v2;
  }
  if (direc[i]<-EPS)
  {
   v1=(max_param[i]-param[i])/direc[i];v2=(min_param[i]-param[i])/direc[i];
   if (v1<cfmin) cfmin=v1;
   if (v2>cfmax) cfmax=v2;
  }
  if ((fabs(direc[i])>EPS)&&(cfprec>(prec_param[i]/direc[i]))) cfprec=prec_param[i]/direc[i];
 }

 //calcul de l'energie correspondant aux parametres de depart
 EParam=energie_transfo_3d(min_parameters);
 if (erreur_critique(err)) goto end_func;

 //en cas de direction nulle on ne fait rien
 if (((cfmax-cfmin)<EPS)||(cfprec<EPS)) return EParam;
 cfu=0.0;
// for (i=0;i<nb_param;i++) paramu[i]=param[i]+cfmin*direc[i];
// min_parameters->transfres->param=paramu;
// Emin=energie_transfo_3d(min_parameters);
// for (i=0;i<nb_param;i++) paramu[i]=param[i]+cfmax*direc[i];
// Emax=energie_transfo_3d(min_parameters);
// min_parameters->transfres->param=param;

 mn_brak (&cfmin, &cfu, &cfmax, &Emin, &Eu, &Emax, min_parameters, direc);

 ///////initialisation///////
 cfx=cfw=cfv=cfu;
 Ew=Ev=Ex=Eu;

 ///////minimisation///////

 for (iter=0;iter<NB_ITER_MAX_LINEMIN_3D;iter++)
 {
  cfm=0.5*(cfmin+cfmax);
  tol1=cfprec*fabs(cfx)+ZEPS;
  tol2=2.0*tol1;
  //si on a atteint la precision voulue, on s'arrete
  if (fabs(cfx-cfm)<=(tol2-0.5*(cfmax-cfmin)))
  {
   break;
  }
  //l'intervalle dans lequel on va rechercher le nouveau point et plus grand que tol1:
  //on va tenter une interpolation parabolique avec les points correspondant a cfx, cfv et cfw
  if (fabs(cfe)>tol1)
  {
   cfr=(cfx-cfw)*(Ex-Ev);
   cfq=(cfx-cfv)*(Ex-Ew);
   cfp=(cfx-cfv)*cfq-(cfx-cfw)*cfr;
   cfq=2.0*(cfq-cfr);
   if (cfq>0.0) cfp=-cfp;
   cfq=fabs(cfq);
   cfetemp=cfe;
   cfe=cfd;
   //on teste si l'interpolation parabolique est pertinente
   //ie les 3 points ne doivent pas etre alignes
   if ((fabs(cfp)>=fabs(0.5*cfq*cfetemp))||(cfp<=cfq*(cfmin-cfx))||(cfp>=cfq*(cfmax-cfx)))
    //ca n'a pas marche: on prend la section du nombre d'or dans le plus grand segment
    cfd=CGOLD*(cfe=(cfx >= cfm ? cfmin-cfx : cfmax-cfx));
   else
   {
    //ca a marche: on prend le minimum de la parabole
    cfd=cfp/cfq;
    cfu=cfx+cfd;
    if (((cfu-cfmin)<tol2)||((cfmax-cfu)<tol2))
     cfd=SGN(tol1,cfm-cfx);
   }
  }
  //l'intervalle dans lequel on va rechercher le nouveau point est trop petit:
  //on choisit un autre intervalle
  else
  {
   cfe=(cfx>=cfm ? cfmin-cfx : cfmax-cfx);
   cfd=CGOLD*cfe;
  }

  //calcul du coef correspondant a l'optimum determine, des parametres et de l'energie associes
  cfu=(fabs(cfd) >= tol1 ? cfx+cfd : cfx + SGN(tol1,cfd));
  Eu=eval_energie_direction(min_parameters, paramu, direc, cfu);
  if (erreur_critique(err)) goto end_func;

  //si le point determine est meilleur que l'optimum precedent:
  if (Eu<=Ex)
  {
   //on remet les donnees a jour avec le nouvel optimum calcule
   if (cfu>=cfx) cfmin=cfx; else cfmax=cfx;
   BRENT_SHFT(cfv,cfw,cfx,cfu)
   BRENT_SHFT(Ev,Ew,Ex,Eu)
   nb_minimisations++;
   if (nb_minimisations>=NB_MINIMISATION_MAX_LINEMIN_3D) break;
  }
  //s'il n'est pas meilleur:
  else
  {
   //on se readapte en consequence
   if (cfu<cfx) cfmin=cfu; else cfmax=cfu;
   if ((Eu<=Ew)||(cfw==cfx))
   {
    cfv=cfw;
    cfw=cfu;
    Ev=Ew;
    Ew=Eu;
   }
   else if ((Eu<=Ev)||(cfv=cfx)||(cfv==cfw))
   {
    cfv=cfu;
    Ev=Eu;
   }
  }
 }

 //calcul et renvoie du point optimal
 cfopt=cfx;
 for (i=0;i<nb_param;i++) param[i]=param[i]+cfopt*direc[i];
 Eopt=Ex;

 if (err) *(err)=NO_ERR;

end_func:

 ///////liberations memoire///////
 if (paramu) FREE(paramu);

 return Eopt;
}

/*! @} */

/*! \ingroup     MinimisationFonction  @{ */
/*******************************************************************************
**     min_icm_3d(imref,imreca,imres,champ_ini,champ_fin,the_transf,inter
**                   ,dist,nb_param,min_param,max_param,prec_param,param)
*/
/*!     Minimisation par methode ICM
**
**     utilisation:
**     \param       imref,imreca,imres: images
**     \param       champ_ini : champ initial
**                              NULL = pas de champ initial
**     \param       champ_fin : champ final resultat de la transfo
**     \param       the_transf : pointeur sur la fonction de transformation
**                           int the_transf(nb_param,param,champ)
**     \param       inter : pointeur sur la fonction d'interpolation
**                           int inter(imdeb,champ,imres)
**     \param       dist : pointeur sur la fonction de distance
**                         double dist(im1,im2)
**     \param       nb_param: nombre de parametre a minimiser
**     \param       min_param: bornes inferieures des parametres
**     \param       max_param: bornes superieures des parametres
**     \param       prec_param: precision demandee pour chaque parametre
**     \param       param: resultats donnant la valeur des parametres au
**                         point definissant le minimum
**     \retval  la fonction retourne un double contenant la valeur du minimum
**      Les parametres correspondant a ce minimum sont retourne dans param
*******************************************************************************/
double min_icm_3d(ptr_minimisation_param min_parameters)
{
 int	nbitermax,nb,d,l;
 double *direc=NULL,Emin=DBL_MAX,Edeb,E1,E2;
 double min_par, max_par;
 double *old_param=NULL;

 int nb_param=min_parameters->transfres->nb_param;
 double *min_param=min_parameters->min_param;
 double *max_param=min_parameters->max_param;
 double *prec_param=min_parameters->prec_param;
 double *param=min_parameters->transfres->param;
 ERREUR_RECALAGE *err=min_parameters->err;

 ////initialisations et allocations memoires////
 nbitermax=NB_ITER_MAX_ICM_3D;

 old_param=CALLOC(nb_param,double);
 direc=CALLOC(nb_param,double);
 if ((!old_param)||(!direc))
 {
  RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans min_icm_3d\n");
 }
 
 //calcul de l'energie de depart
 Emin=energie_transfo_3d(min_parameters);
 if (erreur_critique(err)) goto end_func;
 Edeb=Emin;
 aff_log("Edebut: %.2f     nb_param: %d \n",Edeb,nb_param);

 ////minimisation////
 nb=0;
 E2=Emin;
 do
 {
  E1=E2;

  //recentrage de l'intervalle de recherche
  for (l=0;l<nb_param;l++)
  {
   min_par=min_param[l];
   max_par=max_param[l];
   min_param[l]=param[l]-(max_par-min_par)/2.0;
   max_param[l]=param[l]+(max_par-min_par)/2.0;
   old_param[l]=param[l];
  }
  
  for (d=0;d<nb_param;d++) //on regarde toutes les directions de l'espace
   if (prec_param[d]!=0.0)
    {//minimisation suivant le parametre d
     direc[d]=1.0;
     E2=linemin_3d(min_parameters, direc);
     if (erreur_critique(err)) goto end_func;
     direc[d]=0.0;
    }

  //si les parametres n'ont pas bouge de plus que la precision recherchee, on s'arrete
  if (arret_selon_parametres(param,old_param,prec_param,nb_param)) break;

  printf("nbiter:%d   Energie: %.2f  %%reduc: %.2f  \r",nb,E2,100.0*(E1-E2)/fabs(E1));
  
  nb++;
  
 } while(nb<nbitermax);

 Emin=E2;

 aff_log("nbiter:%d   Energie: %.2f  %%reduc: %.2f  \n",nb,Emin,100.0*(Edeb-Emin)/fabs(Edeb));

 if (err) *err=NO_ERR;

end_func:

 ////on termine////
 if (old_param) FREE(old_param);
 if (direc) FREE(direc);
 
 return(Emin);
}
/*! @} */

/*! \ingroup     MinimisationFonction  @{ */
/*******************************************************************************
**     min_simplex_3d(imref,imreca,imres,champ_ini,champ_fin,the_transf,inter
**                   ,dist,nb_param,min_param,max_param,prec_param,param)
**
*/
/*!    Minimisation par methode de minimisation simplex: implementation telle
**     decrite dans le "numerical recipes in C":
**     http://www.library.cornell.edu/nr/bookcpdf.html
**
**     utilisation:
**     \param       imref,imreca,imres: images
**     \param       champ_ini : champ initial
**                              NULL = pas de champ initial
**     \param       champ_fin : champ final resultat de la transfo
**     \param       the_transf : pointeur sur la fonction de transformation
**                           int the_transf(nb_param,param,champ)
**     \param       inter : pointeur sur la fonction d'interpolation
**                           int inter(imdeb,champ,imres)
**     \param       dist : pointeur sur la fonction de distance
**                         double dist(im1,im2)
**     \param       nb_param: nombre de parametre a minimiser
**     \param       min_param: bornes inferieures des parametres
**     \param       max_param: bornes superieures des parametres
**     \param       prec_param: precision demandee pour chaque parametre
**     \param       param: resultats donnant la valeur des parametres au
**                         point definissant le minimum
**     \retval       la fonction retourne un double contenant la valeur du minimum
**      Les parametres correspondant a ce minimum sont retourne dans param
*******************************************************************************/ 
double min_simplex_3d(ptr_minimisation_param min_parameters)
{
 int i, ihi=0, ilo=0, inhi=0, j, Nest, iter=0;
 double Esav, Ereflexion, Eopt=DBL_MAX, Edebut;
 double *psum=NULL, *energies=NULL, *p_reflexion=NULL;
 double **simplexe=NULL;

 int nb_param=min_parameters->transfres->nb_param;
 double *min_param=min_parameters->min_param;
 double *max_param=min_parameters->max_param;
 double *prec_param=min_parameters->prec_param;
 double *param=min_parameters->transfres->param;
 ERREUR_RECALAGE *err=min_parameters->err;

 Nest=1;
  for (i=0;i<nb_param;i++)
   if (prec_param[i]!=0.0) {Nest++;}

 if (Nest<4)
 {
  RAISE_ERR_TEXT(err, PB_PARAMETRES, "Trop peu de parametres a estimer pour la methode simplex\n");
 }

 //INITIALISATION
 simplexe=alloc_dmatrix(Nest,nb_param);
 psum=CALLOC(nb_param, double);
 energies=CALLOC(Nest, double);
 p_reflexion=CALLOC(nb_param, double);
 if ((!simplexe)||(!psum)||(!energies)||(!p_reflexion))
 {
  RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans min_simplex_3d\n");
 }

 //remplissage du simplex initial
 //mise a zero du simplex
 for (i=0;i<Nest;i++)
  for (j=0;j<nb_param;j++)
     simplexe[i][j]=0.0;
 //le point P0
 for (i=0;i<nb_param;i++)
  for (j=0;j<Nest;j++)
    simplexe[j][i]=param[i];

 min_parameters->transfres->param=simplexe[0];
 energies[0]=energie_transfo_3d(min_parameters);
 min_parameters->transfres->param=param;
 if (erreur_critique(err)) goto end_func;

 //les points 1 a N
 //on va prendre pour chaque point max_param ou min_param selon la direction de moindre energie
 j=1;
 for (i=0;i<nb_param;i++)
 {
  if (prec_param[i]!=0.0)
  {
   simplexe[j][i]=max_param[i];

   min_parameters->transfres->param=simplexe[j];
   energies[j]=energie_transfo_3d(min_parameters);
   min_parameters->transfres->param=param;

   if (erreur_critique(err)) goto end_func;
   Esav=energies[j];
   simplexe[j][i]=min_param[i];

   min_parameters->transfres->param=simplexe[j];
   energies[j]=energie_transfo_3d(min_parameters);
   min_parameters->transfres->param=param;

   if (erreur_critique(err)) goto end_func;
   if (energies[j]>Esav)
   {
    simplexe[j][i]=max_param[i];
    energies[j]=Esav;
   }

   j++;
  }
 }

 //calcul des sommes des valeurs des points du simplexe pour chaque parametre
 get_simplexe_psum(simplexe, psum, Nest, nb_param);

 //calcul et affichage de l'energie de depart
 Edebut=energie_transfo_3d(min_parameters);

 if (erreur_critique(err)) goto end_func;
 aff_log("Edebut: %.8f     nb_param: %d \n",Edebut,nb_param);

 //MINIMISATION
 while (iter<NB_ITER_MAX_SIMPLEX_3D)
 {
  ilo=0;
  //determination du plus mauvais, du deuxieme plus mauvais et du meilleur point
  if (energies[0]>energies[1])
   { ihi=0; inhi=1; }
  else
   { ihi=1; inhi=0; }
  for (i=0;i<Nest;i++)
  {
   if (energies[i]<energies[ilo]) ilo=i;
   if (energies[i]>energies[ihi])
    { inhi=ihi; ihi=i; }
   else if ((energies[i]>energies[inhi])&&(i!=ihi))
    inhi=i;
  }

  if (Edebut!=0.0) printf("iteration: %d  energie:%.8f  reduction : %.2f %% \r",iter,energies[ilo],(Edebut-energies[ilo])/fabs(Edebut)*100.0);
  else printf("iteration: %d  energie:%.8f  reduction : %.2f %% \r",iter,energies[ilo],0.0);

  //critere d'arret
  if (arret_selon_parametres(simplexe[ilo], simplexe[ihi], prec_param, nb_param)) break;
//  if (arret_selon_energie(energies[ilo], energies[ihi])) break;

  //nouvelle iteration
  min_parameters->transfres->param=p_reflexion;
  Ereflexion=reflexion_simplex_3d(min_parameters, simplexe, energies, psum, ihi, Nest, -1.0);
  min_parameters->transfres->param=param;
  if (erreur_critique(err)) goto end_func;
  //si la reflexion donne un meilleur resultat que le meilleur des points
  //on continue dans cette direction par un facteur 2
  if (Ereflexion<energies[ilo])
  {
   min_parameters->transfres->param=p_reflexion;
   Ereflexion=reflexion_simplex_3d(min_parameters, simplexe, energies, psum, ihi, Nest, 2.0);
   min_parameters->transfres->param=param;
   if (erreur_critique(err)) goto end_func;
  }
  //si on a un moins bon resultat que le deuxieme plus mauvais point
  //on teste une reflexion intermediaire
  else if (Ereflexion>=energies[inhi])
  {
   Esav=energies[ihi];

   min_parameters->transfres->param=p_reflexion;
   Ereflexion=reflexion_simplex_3d(min_parameters, simplexe, energies, psum, ihi, Nest, 0.5);
   min_parameters->transfres->param=param;

   if (erreur_critique(err)) goto end_func;
   //on n'a toujours pas ameliore le score du plus mauvais point:
   //on contracte le simplexe autour du meilleur point
   if (Ereflexion>=Esav)
   {
    for (i=0;i<Nest;i++)
    {
     if (i!=ilo)
     {
      for (j=0;j<nb_param;j++)
       simplexe[i][j]=psum[j]=0.5*(simplexe[i][j]+simplexe[ilo][j]);
      
      min_parameters->transfres->param=psum;
      energies[i]=energie_transfo_3d(min_parameters);
      min_parameters->transfres->param=param;

      if (erreur_critique(err)) goto end_func;
     }
    }
    get_simplexe_psum(simplexe, psum, Nest, nb_param);
   }

  }
  iter++;

 }

 //determination du meilleur point si on s'est arrete avant la convergence
 if (iter>=NB_ITER_MAX_SIMPLEX_3D)
 {
  ilo=0;
  for (i=0;i<Nest;i++)
   if (energies[i]<energies[ilo]) ilo=i;
 }

 for (i=0;i<nb_param;i++)  param[i]=simplexe[ilo][i];
 Eopt=energies[ilo];
 if (Edebut!=0.0) printf("iteration: %d  energie:%.8f  reduction : %.2f %% \n",iter,Eopt,(Edebut-Eopt)/fabs(Edebut)*100.0);
 else printf("iteration: %d  energie:%.8f  reduction : %.2f %% \n",iter,Eopt,0.0);

 if (err) *err=NO_ERR;

end_func:

 //liberation memoire
 if (p_reflexion) FREE(p_reflexion);
 if (energies) FREE(energies);
 if (psum) FREE(psum);
 if (simplexe) free_dmatrix(simplexe,Nest, nb_param);

 return Eopt;
}

/*! @} */

/*! \ingroup     MinimisationFonction  @{ */
/*******************************************************************************
**     min_desc_grad_3d(imref,imreca,imres,champ_ini,champ_fin,the_transf,inter
**                   ,dist,nb_param,min_param,max_param,prec_param,param)
*/
/*!    Minimisation par descente de gradient
**
**
**     utilisation:
**    \param        imref,imreca,imres: images
**    \param        champ_ini : champ initial
**                              NULL = pas de champ initial
**    \param        champ_fin : champ final resultat de la transfo
**    \param        the_transf : pointeur sur la fonction de transformation
**                           int the_transf(nb_param,param,champ)
**    \param        inter : pointeur sur la fonction d'interpolation
**                           int inter(imdeb,champ,imres)
**    \param        dist : pointeur sur la fonction de distance
**                         double dist(im1,im2)
**    \param        nb_param: nombre de parametre a minimiser
**    \param        min_param: bornes inferieures des parametres
**    \param        max_param: bornes superieures des parametres
**    \param        prec_param: precision demandee pour chaque parametre
**    \param        param: resultats donnant la valeur des parametres au
**                         point definissant le minimum
**    \retval       la fonction retourne un double contenant la valeur du minimum
**      Les parametres correspondant a ce minimum sont retourne dans param
*******************************************************************************/
double min_desc_grad_3d(ptr_minimisation_param min_parameters)
{
 int wdth,hght,dpth,l;
 field3d  *gradIm=NULL;
 double *grad=NULL;
 gradient_func_t gradient=NULL;
 double Edebut,Efin=DBL_MAX,E1,E2;
 int	iter;
 double min_par, max_par;
 double *old_param=NULL;

 grphic3d *imref=min_parameters->imref;
 grphic3d *imres=min_parameters->imres;
 dist_func_t dist=min_parameters->dist;
 int nb_param=min_parameters->transfres->nb_param;
 double *min_param=min_parameters->min_param;
 double *max_param=min_parameters->max_param;
 double *prec_param=min_parameters->prec_param;
 double *param=min_parameters->transfres->param;
 ptr_distance_param dist_param=min_parameters->dist_param;
 ERREUR_RECALAGE *err=min_parameters->err;


 wdth=imref->width;hght=imref->height;dpth=imref->depth;

 //mise en place des valeurs des fonctions en fonction de the_transf
 gradient=imx_choose_gradient_fct(min_parameters->transfres->typetrans,dist);

 //allocation memoire des variables
 old_param=CALLOC(nb_param,double);
 grad=alloc_dvector(nb_param);
 gradIm=cr_field3d(wdth,hght,dpth);
 if ((!old_param)||(!grad)||(!gradIm))
 {
  RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans min_desc_grad_3d\n");
 }

 //calcul de l'energie de depart
 Edebut=energie_transfo_3d(min_parameters);
 if (erreur_critique(err)) goto end_func;
 E2=Edebut;

 aff_log("Edebut: %.2f   nb_param: %d \n",Edebut,nb_param);

 iter=0;
 do
 {
  E1=E2;

  //recentrage de l'intervalle de recherche
  for (l=0;l<nb_param;l++)
  {
   min_par=min_param[l];
   max_par=max_param[l];
   min_param[l]=param[l]-(max_par-min_par)/2.0;
   max_param[l]=param[l]+(max_par-min_par)/2.0;
   old_param[l]=param[l];
  }
  
  //calcul du gradient au point defini par p
//  if (dist_utilisant_VP(dist)) calcul_histogrammes_gradient_VP(dist_param);
//  if (erreur_critique(err)) goto end_func;
  gradient(imref,imres,nb_param,param,gradIm,grad,dist_param->donnees_VP);

  //minimisation de l'energie suivant la direction donne par grad
  E2=linemin_3d(min_parameters, grad);
  if (erreur_critique(err)) goto end_func;

  //si les parametres n'ont pas bouge de plus que la precision recherchee, on s'arrete
  if (arret_selon_parametres(param,old_param,prec_param,nb_param)) break;

  printf("nbiter:%d   Energie: %.2f  %%reduc: %.2f  \r",iter,E2,100.0*(E1-E2)/fabs(E1));

  iter++;

 } while (iter<NB_ITER_MAX_DESC_GRAD_3D); 

  //calcul de l'energie finale
  Efin=energie_transfo_3d(min_parameters);
  if (erreur_critique(err)) goto end_func;
  aff_log("nbiter: %d    Efin: %.2f  reduction: %.2f %%  \n",iter,Efin,(Edebut-Efin)/fabs(Edebut)*100.0);

 if (err) *err=NO_ERR;

end_func:

 //liberation memoire des variables
 if (grad) free_dvector(grad,nb_param);
 if (gradIm) free_field3d(gradIm);
 if (old_param) FREE(old_param);

 return(Efin);
}
/*! @} */

/*! \ingroup     MinimisationFonction  @{ */
/*******************************************************************************
**     min_grad_conj_3d(imref,imreca,imres,champ_ini,champ_fin,the_transf,inter
**                   ,dist,nb_param,min_param,max_param,prec_param,param)
*/
/*!     Minimisation par la methode du gradient conjugue
**      methode de construction de direction de Polak-Ribiere
**
**     utilisation:
**     \param       imref,imreca,imres: images
**     \param		champ_ini : champ initial
**                              NULL = pas de champ initial
**     \param       champ_fin : champ final resultat de la transfo
**     \param       the_transf : pointeur sur la fonction de transformation
**                           int the_transf(nb_param,param,champ)
**     \param       inter : pointeur sur la fonction d'interpolation
**                           int inter(imdeb,champ,imres)
**     \param       dist : pointeur sur la fonction de distance
**                         double dist(im1,im2)
**     \param       nb_param: nombre de parametre a minimiser
**     \param       min_param: bornes inferieures des parametres
**     \param       max_param: bornes superieures des parametres
**     \param       prec_param: precision demandee pour chaque parametre
**     \param       param: resultats donnant la valeur des parametres au
**                         point definissant le minimum
**     \retval     la fonction retourne un double contenant la valeur du minimum
**      Les parametres correspondant a ce minimum sont retourne dans param
*******************************************************************************/
double min_grad_conj_3d(ptr_minimisation_param min_parameters)
{
 int wdth,hght,dpth,l;
 field3d  *gradIref=NULL;
 double *g=NULL,*h=NULL,*xi=NULL;
 gradient_func_t gradient=NULL;
 double Edebut,Efin=DBL_MAX,E1,E2;
 int	nb;
 double gg,dgg;
 double min_par, max_par;
 double *old_param=NULL;

 grphic3d *imref=min_parameters->imref;
 grphic3d *imres=min_parameters->imres;
 dist_func_t dist=min_parameters->dist;
 int nb_param=min_parameters->transfres->nb_param;
 double *min_param=min_parameters->min_param;
 double *max_param=min_parameters->max_param;
 double *prec_param=min_parameters->prec_param;
 double *param=min_parameters->transfres->param;
 ptr_distance_param dist_param=min_parameters->dist_param;
 ERREUR_RECALAGE *err=min_parameters->err;

  
 wdth=imref->width;hght=imref->height;dpth=imref->depth;

 //mise en place des valeurs des fonctions en fonction de the_transf
 gradient=imx_choose_gradient_fct(min_parameters->transfres->typetrans,dist);

 //allocation memoire des variables
 old_param=alloc_dvector(nb_param);
 g=alloc_dvector(nb_param);
 h=alloc_dvector(nb_param);
 xi=alloc_dvector(nb_param);
 gradIref=cr_field3d(wdth,hght,dpth);
 if ((!old_param)||(!g)||(!h)||(!xi)||(!gradIref))
 {
  RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans min_grad_conj_3d\n");
 }

 //calcul de l'energie de depart
  Edebut=energie_transfo_3d(min_parameters);
  if (erreur_critique(err)) goto end_func;
  E2=Edebut;

  //calcul du gradient au point defini par p
//  if (dist_utilisant_VP(dist)) calcul_histogrammes_gradient_VP(dist_param);
//  if (erreur_critique(err)) goto end_func;
  gradient(imref,imres,nb_param,param,gradIref,xi,dist_param->donnees_VP);
  for (l=0;l<nb_param;l++) {g[l]=-xi[l];xi[l]=h[l]=g[l];}

 aff_log("Edebut: %.2f   nb_param: %d \n",Edebut,nb_param);

 nb=0;
 do
 {
  E1=E2;

  //recentrage de l'intervalle de recherche
  for (l=0;l<nb_param;l++)
  {
   min_par=min_param[l];
   max_par=max_param[l];
   min_param[l]=param[l]-(max_par-min_par)/2.0;
   max_param[l]=param[l]+(max_par-min_par)/2.0;
   old_param[l]=param[l];
  }
  
  //minimisation de l'energie suivant la direction donne par gradk1
  E2=linemin_3d(min_parameters, xi);
  if (erreur_critique(err)) goto end_func;

  //si les parametres n'ont pas bouge de plus que la precision recherchee, on s'arrete
  if (arret_selon_parametres(param,old_param,prec_param,nb_param)) break;

  printf("nbiter:%d   Energie: %.2f  %%reduc: %.2f  \r",nb,E2,100.0*(E1-E2)/fabs(E1));
  
  //calcul du gradient au point defini par param
  //methode de Polak-Ribiere
//  if (dist_utilisant_VP(dist)) calcul_histogrammes_gradient_VP(dist_param);
//  if (erreur_critique(err)) goto end_func;
  gradient(imref,imres,nb_param,param,gradIref,xi,dist_param->donnees_VP);

  dgg=gg=0.0;
  for (l=0;l<nb_param;l++) {gg+=g[l]*g[l];dgg+=(xi[l]+g[l])*xi[l];}
  for (l=0;l<nb_param;l++) {g[l]=-xi[l];xi[l]=h[l]=g[l]+dgg/gg*h[l];}

  nb++;

 } while (gg!=0.0 && nb<NB_ITER_MAX_GRAD_CONJ_3D);

 //calcul de l'energie finale
 Efin=energie_transfo_3d(min_parameters);
 if (erreur_critique(err)) goto end_func;
 aff_log("nbiter: %d    Efin: %.2f  reduction: %.2f %%  \n",nb,Efin,(Edebut-Efin)/fabs(Edebut)*100.0);

 if (err) *err=NO_ERR;

end_func:

 //liberation memoire des variables
 if (gradIref) free_field3d(gradIref);
 if (g) free_dvector(g,nb_param);
 if (h) free_dvector(h,nb_param);
 if (xi) free_dvector(xi,nb_param);
 if (old_param) free_dvector(old_param,nb_param);

 return(Efin);
}
/*! @}  */

/*! \ingroup     MinimisationFonction  @{ */
/*******************************************************************************
**     min_desc_grad_mod_3d(imref,imreca,imres,champ_ini,champ_fin,the_transf,inter
**                   ,dist,nb_param,min_param,max_param,prec_param,param)
*/
/*!
**     Minimisation par la methode de descente de gradient modifiee
**
**     utilisation:
**     \param       imref,imreca,imres: images
**     \param       champ_ini : champ initial
**                              NULL = pas de champ initial
**     \param       champ_fin : champ final resultat de la transfo
**     \param       the_transf : pointeur sur la fonction de transformation
**                           int the_transf(nb_param,param,champ)
**     \param       inter : pointeur sur la fonction d'interpolation
**                           int inter(imdeb,champ,imres)
**     \param       dist : pointeur sur la fonction de distance
**                         double dist(im1,im2)
**     \param       nb_param: nombre de parametre a minimiser
**     \param       min_param: bornes inferieures des parametres
**     \param       max_param: bornes superieures des parametres
**     \param       prec_param: precision demandee pour chaque parametre
**     \param       param: resultats donnant la valeur des parametres au
**                         point definissant le minimum
**     \retval     la fonction retourne un double contenant la valeur du minimum
**      Les parametres correspondant a ce minimum sont retourne dans param
*******************************************************************************/
double min_desc_grad_mod_3d(ptr_minimisation_param min_parameters)
{
 int wdth,hght,dpth,i;
 field3d  *gradIref=NULL;
 double *g1=NULL,*g2=NULL,*p1=NULL,*p2=NULL;
 gradient_func_t gradient=NULL;
 double Edebut,Efin=DBL_MAX,E1,E2;
 int	nb,flag;
 double lambda,prstop;

 grphic3d *imref=min_parameters->imref;
 grphic3d *imres=min_parameters->imres;
 dist_func_t dist=min_parameters->dist;
 int nb_param=min_parameters->transfres->nb_param;
 double *min_param=min_parameters->min_param;
 double *max_param=min_parameters->max_param;
 double *prec_param=min_parameters->prec_param;
 double *param=min_parameters->transfres->param;
 ptr_distance_param dist_param=min_parameters->dist_param;
 ERREUR_RECALAGE *err=min_parameters->err;

 
 wdth=imref->width;hght=imref->height;dpth=imref->depth;

 //mise en place des valeurs des fonctions en fonction de the_transf
 gradient=imx_choose_gradient_fct(min_parameters->transfres->typetrans,dist);

 //allocation memoire des variables
 p1=alloc_dvector(nb_param);p2=alloc_dvector(nb_param);
 g1=alloc_dvector(nb_param);g2=alloc_dvector(nb_param);
 gradIref=cr_field3d(wdth,hght,dpth);
 if ((!p1)||(!p2)||(!g1)||(!g2)||(!gradIref))
 {
  RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans min_desc_grad_mod_3d\n");
 }

 //calcul de l'energie de depart 
  for (i=0;i<nb_param;i++) p2[i]=(min_param[i]+max_param[i])/2.0;
  min_parameters->transfres->param=p2;
  Edebut=energie_transfo_3d(min_parameters);
  min_parameters->transfres->param=param;
  if (erreur_critique(err)) goto end_func;
  E2=Edebut;

 //calcul du gradient au point defini par p2
//  if (dist_utilisant_VP(dist)) calcul_histogrammes_gradient_VP(dist_param);
//  if (erreur_critique(err)) goto end_func;
  gradient(imref,imres,nb_param,p2,gradIref,g2,dist_param->donnees_VP);

 //initialisation de lambda
 lambda=(max_param[0]-min_param[0])/fabs(g2[0]);
 for (i=1;i<nb_param;i++)
  if  ((max_param[i]-min_param[i])/fabs(g2[i])<lambda && prec_param[i]!=0.0) lambda=(max_param[i]-min_param[i])/fabs(g2[i]);

 aff_log("Edebut: %.2f   nb_param: %d \n",Edebut,nb_param);

 nb=flag=0;prstop=0.1/100.0;
 do
 {
  E1=E2;for (i=0;i<nb_param;i++) {p1[i]=p2[i];g1[i]=g2[i];}

  //calcul de p2 et de E2
  for (i=0;i<nb_param;i++) p2[i]=p1[i]-lambda*g1[i];
  min_parameters->transfres->param=p2;
  E2=energie_transfo_3d(min_parameters);
  min_parameters->transfres->param=param;
  if (erreur_critique(err)) goto end_func;

  //for (i=0;i<nb_param;i++) printf("%.2f ",p1[i]);printf("        %.2f \n",E1);
  //for (i=0;i<nb_param;i++) printf("%.2f ",p2[i]);printf("        %.2f \n",E2);

  if (E2<=E1)
   {//l'energie diminue
     lambda=lambda*2.0;
     if ((E1-E2)/E1<prstop*fabs(E1)) flag++; //l'energie a diminue de tres peu
    //calcul du gradient au point defini par p2
//     if (dist_utilisant_VP(dist)) calcul_histogrammes_gradient_VP(dist_param);
//     if (erreur_critique(err)) goto end_func;
     gradient(imref,imres,nb_param,p2,gradIref,g2,dist_param->donnees_VP);
   }
  else
   {//l'energie augmente
    lambda=lambda/2.0;
    //on retourne en p1
    for (i=0;i<nb_param;i++) {p2[i]=p1[i];g2[i]=g1[i];}
    E2=E1;
   }

  nb++;
  printf("nbiter: %d  energie: %.2f reduction: %.2f %% \r",nb,E2,100.0*(E1-E2)/fabs(E1));
 } while (nb<NB_ITER_MAX_DESC_GRAD_MOD_3D && flag<10);

  //calcul de l'energie finale
  min_parameters->transfres->param=p1;
  Efin=energie_transfo_3d(min_parameters);
  min_parameters->transfres->param=param;

  if (erreur_critique(err)) goto end_func;
  aff_log("nbiter: %d    Efin: %.2f  reduction: %.2f %%  \n",nb,Efin,(Edebut-Efin)/fabs(Edebut)*100.0);

 //resultat dans param
  for (i=0;i<nb_param;i++) param[i]=p1[i];

 if (err) *err=NO_ERR;

end_func:

 //liberation memoire des variables
 if (p1) free_dvector(p1,nb_param);
 if (p2) free_dvector(p2,nb_param);
 if (g1) free_dvector(g1,nb_param);
 if (g2) free_dvector(g2,nb_param);
 if (gradIref) free_field3d(gradIref);

 return(Efin);
}
/*! @}  */


/*! \ingroup     MinimisationFonction  @{ */
/*******************************************************************************
**     min_quasi_newton_mod_3d(imref,imreca,imres,champ_ini,champ_fin,the_transf,inter
**                   ,dist,nb_param,min_param,max_param,prec_param,param)
*/
/*!     Minimisation par la methode de quasi newton modifiee
**
**
**     utilisation:
**     \param       imref,imreca,imres: images
**     \param       champ_ini : champ initial
**                              NULL = pas de champ initial
**     \param       champ_fin : champ final resultat de la transfo
**     \param       the_transf : pointeur sur la fonction de transformation
**                           int the_transf(nb_param,param,champ)
**     \param       inter : pointeur sur la fonction d'interpolation
**                           int inter(imdeb,champ,imres)
**     \param       dist : pointeur sur la fonction de distance
**                         double dist(im1,im2)
**     \param       nb_param: nombre de parametre a minimiser
**     \param       min_param: bornes inferieures des parametres
**     \param       max_param: bornes superieures des parametres
**     \param       prec_param: precision demandee pour chaque parametre
**     \param       param: resultats donnant la valeur des parametres au
**                         point definissant le minimum
**     \retval       la fonction retourne un double contenant la valeur du minimum
**      Les parametres correspondant a ce minimum sont retourne dans param
*******************************************************************************/
double min_quasi_newton_mod_3d(ptr_minimisation_param min_parameters)
{
 int wdth,hght,dpth,i,j;
 field3d  *gradIref=NULL;
 double *g1=NULL,*g2=NULL,*p1=NULL,*p2=NULL,*popt=NULL;
 gradient_func_t gradient=NULL;
 double Edebut,Efin=DBL_MAX,E1,E2,Eopt;
 int	nb,flag;
 double **H=NULL;
 double lambda,prstop;

 grphic3d *imref=min_parameters->imref;
 grphic3d *imres=min_parameters->imres;
 dist_func_t dist=min_parameters->dist;
 int nb_param=min_parameters->transfres->nb_param;
 double *min_param=min_parameters->min_param;
 double *max_param=min_parameters->max_param;
 double *prec_param=min_parameters->prec_param;
 double *param=min_parameters->transfres->param;
 ptr_distance_param dist_param=min_parameters->dist_param;
 ERREUR_RECALAGE *err=min_parameters->err;

 
 wdth=imref->width;hght=imref->height;dpth=imref->depth;

 /*mise en place des valeurs des fonctions en fonction de the_transf*/
 gradient=imx_choose_gradient_fct(min_parameters->transfres->typetrans,dist);

 /*allocation memoire des variables*/
 p1=alloc_dvector(nb_param);p2=alloc_dvector(nb_param);popt=alloc_dvector(nb_param);
 g1=alloc_dvector(nb_param);g2=alloc_dvector(nb_param);
 gradIref=cr_field3d(wdth,hght,dpth);
  if ((! p1) || (! p2) || (! popt) || (! g1) || (! g2) || (! gradIref))
  {
    RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans min_quasi_newton_mod_3d\n");
  }

 /*calcul de l'energie de depart */
  for (i=0;i<nb_param;i++) p2[i]=popt[i]=(min_param[i]+max_param[i])/2.0;
  min_parameters->transfres->param=p2;
  Edebut=energie_transfo_3d(min_parameters);
  min_parameters->transfres->param=param;
  if (erreur_critique(err)) goto end_func;
  E2=Edebut;Eopt=Edebut;

 /*calcul du gradient au point defini par p2*/
//  if (dist_utilisant_VP(dist)) calcul_histogrammes_gradient_VP(dist_param);
//  if (erreur_critique(err)) goto end_func;
  gradient(imref,imres,nb_param,p2,gradIref,g2,dist_param->donnees_VP);

 /*initialisation de lambda*/
 lambda=(max_param[0]-min_param[0])/fabs(g2[0]);
 for (i=1;i<nb_param;i++)
  if  ((max_param[i]-min_param[i])/fabs(g2[i])<lambda && prec_param[i]!=0.0) lambda=(max_param[i]-min_param[i])/fabs(g2[i]);

/*allocation et initialisation des matrices Hessiennes*/
 H=alloc_dmatrix(nb_param,nb_param);
 for (i=0;i<nb_param;i++)
  for (j=0;j<nb_param;j++)
   if (i==j) H[i][j]=lambda/20.0; else H[i][j]=0.0;

 aff_log("Edebut: %.2f   nb_param: %d \n",Edebut,nb_param);

 nb=flag=0;prstop=0.1/100.0;
 do
 {
  E1=E2;for (i=0;i<nb_param;i++) {p1[i]=p2[i];g1[i]=g2[i];}
  /*for (i=0;i<nb_param;i++) for (j=0;j<nb_param;j++) H1[i][j]=H2[i][j];*/

  /*calcul de p2 et de E2*/
  for (i=0;i<nb_param;i++) p2[i]=p1[i];
  for (i=0;i<nb_param;i++)
   for (j=0;j<nb_param;j++) p2[i]=p2[i]-H[i][j]*g1[j];
  min_parameters->transfres->param=p2;
  E2=energie_transfo_3d(min_parameters);
  min_parameters->transfres->param=param;
  if (erreur_critique(err)) goto end_func;

  /*if (100.0*(E1-E2)/E1<0.5) flag++;
  if (100.0*(E2-E1)/E1>5)
  */
  if (100.0*(E1-E2)/fabs(E1)<0.5) flag++;
  if (100.0*(E2-E1)/fabs(E1)>5)
  {/*l'energie augmente trop*/
   for (i=0;i<nb_param;i++)
    for (j=0;j<nb_param;j++)
     if (i==j) H[i][j]=lambda; else H[i][j]=0.0;
   for (i=0;i<nb_param;i++) {p2[i]=p1[i];g2[i]=g1[i];}
   E2=E1;
  }
  else
  {
  /*calcul du gradient au point defini par p2*/
//   if (dist_utilisant_VP(dist)) calcul_histogrammes_gradient_VP(dist_param);
//   if (erreur_critique(err)) goto end_func;
   gradient(imref,imres,nb_param,p2,gradIref,g2,dist_param->donnees_VP);
   update_matrix_BFGS(H,H,p1,p2,g1,g2,nb_param);
  }

  if (E2<Eopt) {Eopt=E2;for (i=0;i<nb_param;i++) popt[i]=p2[i];}

  nb++;
  printf("nbiter: %d  energie: %.2f reduction: %.2f %% flag: %d \r",nb,E2,100.0*(E1-E2)/fabs(E1),flag);
 } while (nb<NB_ITER_MAX_QUASI_NEWTON_MOD_3D && flag<5);

  //resultat dans param
  for (i=0;i<nb_param;i++) param[i]=popt[i];

  //calcul de l'energie finale
  Efin=energie_transfo_3d(min_parameters);
  if (erreur_critique(err)) goto end_func;
  aff_log("nbiter: %d    Efin: %.2f  reduction: %.2f %%  \n",nb,Efin,(Edebut-Efin)/fabs(Edebut)*100.0);



 if (err) *err=NO_ERR;

end_func:

 /*liberation memoire des variables*/
 if (p1) free_dvector(p1,nb_param);
 if (p2) free_dvector(p2,nb_param);
 if (popt) free_dvector(popt,nb_param);
 if (g1) free_dvector(g1,nb_param);
 if (g2) free_dvector(g2,nb_param);
 if (gradIref) free_field3d(gradIref);
 if (H) free_dmatrix(H,nb_param,nb_param);

 return(Efin);
}
/*! @}  */


/*! \ingroup     MinimisationFonction  @{ */
/*******************************************************************************
**     min_flux_local_3d(imref,imreca,imres,champ_ini,champ_fin,the_transf,inter
**                   ,dist,nb_param,min_param,max_param,prec_param,param)
*/
/*!    Minimisation par la methode de quasi newton modifiee
**
**     utilisation:
**                \param imref,imreca,imres: images
**                \param   champ_ini : champ initial
**                              NULL = pas de champ initial
**                \param   champ_fin : champ final resultat de la transfo
**                \param   the_transf : pointeur sur la fonction de transformation
**                           int the_transf(nb_param,param,champ)
**                \param   inter : pointeur sur la fonction d'interpolation
**                           int inter(imdeb,champ,imres)
**                \param   dist : pointeur sur la fonction de distance
**                         double dist(im1,im2)
**                \param   nb_param: nombre de parametre a minimiser
**                \param   min_param: bornes inferieures des parametres
**                \param   max_param: bornes superieures des parametres
**                \param   prec_param: precision demandee pour chaque parametre
**                \param   param: resultats donnant la valeur des parametres au
**                         point definissant le minimum
**      		  \retval la fonction retourne un double contenant la valeur du minimum
**      					Les parametres correspondant a ce minimum sont retourne dans param
**
**   \remark   ATTENTION, cette methode de minimisation est basee sur la decomposition
**	en ondelettes du champ de deformation et ne peut-etre utilise dans aucun
**	autre cas. C'est egalement basee sur la minimisation locale de l'erreur
**	quadratique donc ca ne peut etre utilise que dans ce cas egalement
*******************************************************************************/
double min_flux_local_3d(ptr_minimisation_param min_parameters)
{
 double E1,E2,E3,Edeb,Emini,f,fdx,fdy,fdz,/*fd,*/f2,gx,gy,gz,diff,sumdiff,dpx,dpy,dpz;
 double	*px,*py,*pz,*fx,*fy,*fz,*dfx,*dfy,*dfz;
 int     *x0,*y0,*z0,*x1,*y1,*z1,*idx,*idy,*idz;
 int	i,j,k,l,x00,x11,y00,y11,z00,z11,im,jm,km,im0,jm0,km0,nbfx,nbfy,nbfz,nb;
 int	wdth,hght,dpth;
 double *p1,*p2;
 double a11,a12,a13,/*a21,*/a22,a23,/*a31,a32,*/a33,b1,b2,b3,det,ai11,ai12,ai13,ai22,ai23,ai33;
 int	support;
 field3d *grad;
 vector3d ***dg;
 /*double	ax,ay,az;*/
 int	***matparam;

 grphic3d *imref=min_parameters->imref;
 grphic3d *imres=min_parameters->imres;
 dist_func_t dist=min_parameters->dist;

 int nb_param=min_parameters->transfres->nb_param;
 double *param=min_parameters->transfres->param;

 wdth=imref->width;hght=imref->height;dpth=imref->depth;
 nb=0;

 /*test si on est dans le bon cas*/
// if (the_transf!=base_to_field_3d || dist!=erreur_quad_3d) {PUT_ERROR("function min_flux_local_3d aborted");return(0.0);}
 if (min_parameters->transfres->typetrans!=BSPLINE3D || dist!=erreur_quad_3d) {PUT_ERROR("function min_flux_local_3d aborted");return(0.0);}

 /*calcul du support de la fonction d'echelle*/
 {double *fx,*dfx;
 fx=dfx=NULL;
 support=(*(BASE3D.scal_func))(10,&fx,&dfx);
 free(fx);free(dfx);
 }

 /*on rempli les pointeurs avec les valeurs actuelles dans la base*/
 px=BASE3D.px;py=BASE3D.py;pz=BASE3D.pz;x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;
 idx=BASE3D.idx;idy=BASE3D.idy;idz=BASE3D.idz;
 nbfx=BASE3D.nbfx;nbfy=BASE3D.nbfy;nbfz=BASE3D.nbfz;

 /*allocation et remplissage de la matrice de correspondance position<=>numero parametre*/
 matparam=alloc_imatrix_3d(nbfx,nbfy,nbfz);
  /*mise a -1 de la matrice*/
   for (i=0;i<nbfx;i++) for (j=0;j<nbfy;j++) for (k=0;k<nbfz;k++) matparam[i][j][k]=-1;
  /*remplissage de la matrice*/
   for (l=0;l<nb_param/3;l++) matparam[idx[l]][idy[l]][idz[l]]=l;



 /*allocation memoire des variables*/
 p1=CALLOC(nb_param,double);p2=CALLOC(nb_param,double);
 for (i=0;i<nb_param;i++) p1[i]=p2[i]=param[i];
 grad=cr_field3d(wdth,hght,dpth);dg=grad->raw;


 /*calcul de l'energie de depart*/
 min_parameters->transfres->param=p2;
 Edeb=E2=energie_transfo_3d(min_parameters);
 min_parameters->transfres->param=param;
 aff_log("Energie de depat: %f \n",E2);

 do{
  E1=E2;
  for (i=0;i<nb_param;i++) p1[i]=p2[i];

  for (im0=0;im0<support;im0++)
   for (jm0=0;jm0<support;jm0++)
    for (km0=0;km0<support;km0++)
    {
     /*calcul du gradient*/
     imx_gradient_3d_p(imres,grad,2,4);

     for (im=im0;im<nbfx;im+=support)
      for (jm=jm0;jm<nbfy;jm+=support)
       for (km=km0;km<nbfz;km+=support)
       {
        l=matparam[im][jm][km];
	if (l!=-1)
	{
	 x00=x0[l];x11=x1[l];y00=y0[l];y11=y1[l];z00=z0[l];z11=z1[l];
	 a11=a12=a13=a22=a23=a33=b1=b2=b3=0.0;sumdiff=0.0;
	  /*remplissage*/
 	  for (i=x00;i<x11;i++)
   	   for (j=y00;j<y11;j++)
    	    for (k=z00;k<z11;k++)
   	     {
     	      diff=imref->mri[i][j][k]-imres->mri[i][j][k];sumdiff+=diff;
     	      gx=dg[i][j][k].x;gy=dg[i][j][k].y;gz=dg[i][j][k].z;
    	      f=fx[i-x00]*fy[j-y00]*fz[k-z00];f2=f*f;
	      fdx=dfx[i-x00]*fy[j-y00]*fz[k-z00];fdy=fx[i-x00]*dfy[j-y00]*fz[k-z00];fdz=fx[i-x00]*fy[j-y00]*dfz[k-z00];
	      a11+=f2*gx*gx;a12+=f2*gx*gy;a13+=f2*gx*gz;
             	            a22+=f2*gy*gy;a23+=f2*gy*gz;
					  a33+=f2*gz*gz;
	      b1+=f*gx*diff;b2+=f*gy*diff;b3+=f*gz*diff;
      	     }
          /*calcul du determinant de la matrice*/
           det=a11*(a22*a33-a23*a23)-a12*(a12*a33-a23*a13)+a13*(a12*a23-a22*a13);
           if (fabs(det)>SEUIL_INVERSION_LOCAL)
           {/*calcul de la matrice inverse*/
            ai11=(a22*a33-a23*a23)/det;ai12=(a23*a13-a12*a33)/det;ai13=(a12*a23-a13*a22)/det;
      		      			ai22=(a11*a33-a13*a13)/det;ai23=(a12*a13-a11*a23)/det;
	   							  ai33=(a11*a22-a12*a12)/det;
	    /*calcul de la variation des coefficients*/
  	    dpx=b1*ai11+b2*ai12+b3*ai13;dpy=b1*ai12+b2*ai22+b3*ai23;dpz=b1*ai13+b2*ai23+b3*ai33;
	    if ((dpx*dpx+dpy*dpy+dpz*dpz)<=12)
	    {p2[3*l]=p1[3*l]+dpx;p2[3*l+1]=p1[3*l+1]+dpy;p2[3*l+2]=p1[3*l+2]+dpz;}

	   }
	}
       }

      min_parameters->transfres->param=p2;
      E3=energie_transfo_3d(min_parameters);
      min_parameters->transfres->param=param;
      printf("nbiter: %d  energie: %.2f reduction: %.2f %% \r",nb,E3,100.0*(E1-E3)/fabs(E1));
    }

  /*calcul de E2*/
   min_parameters->transfres->param=p2;
   E2=energie_transfo_3d(min_parameters);
   min_parameters->transfres->param=param;
   nb++;printf("nbiter: %d  energie: %.2f reduction: %.2f %% \r",nb,E2,100.0*(E1-E2)/fabs(E1));
 }while ((E1-E2)/fabs(E1)>SEUIL_FIN_LOCAL);

 if (E2<E1) for (i=0;i<nb_param;i++) p1[i]=p2[i];
 min_parameters->transfres->param=p1;
 Emini=energie_transfo_3d(min_parameters);
 min_parameters->transfres->param=param;

 aff_log("nbiter: %d    Efin: %.2f  reduction: %.2f %%  \n",nb,Emini,(Edeb-Emini)/fabs(Edeb)*100.0);
 for (i=0;i<nb_param;i++) param[i]=p1[i];
 aff_log("Energie finale %f \n",Emini);
 free_field3d(grad);free(p1);free(p2);
 free_imatrix_3d(matparam);
 
 return(Emini);
}

/*! @} */

/*! \ingroup     MinimisationFonction  @{ */
// ================================================================================================

/*! \brief Powell method

*   \param nNbParam       Number of parameters to estimate.
*
*   \param dTabParam      Table of parameter vectors. \a dTabParam is a 2D table defined as:
*                         \a dTabParam [(\a nNbParam + 1)][\a nNbParam ]
*
*   \param dTabDirection  Table of direction vectors. \a dTabDirection is a 2D table defined as:
*                         \a dTabDirection [\a nNbParam ][\a nNbParam ]
*
*   \param dTolerance     Tolerance threshold used by the linear optimization function.
*                         Needed by \a pLinMinFct
*
*   \param OptType        Way of optimization: minimization or maximization.
*                         Needed by \a pLinMinFct.
*
*   \param dDeltaMax      Half interval around the initial position in which the linear optimization
*                         will look for the optimum.
*                         Needed by \a pLinMinFct.
*
*   \param pCostFct       Pointer to the cost function.
*
*   \param pLinMinFct     Pointer to the linear optimization function.
*
*   Powell is a method for the minimization or maximization of a function \a pCostFct in a \a
*   nNbParam - dimension space. This method uses \a pLinMinFct to optimze \a pCostFct according to
*   one variable. Powell method will find a direction set for \a pLinMinFct: \a pCostFct will be
*   optimized along each direction, without interfering with the optimizations along
*   previous directions.
*/
// ================================================================================================
double min_powel_3d(ptr_minimisation_param min_parameters)
{
 double retValue = DBL_MAX;
 int i,j,l;

 double **dTabParam=NULL;
 double **dTabDirection=NULL;

 int nb_param=min_parameters->transfres->nb_param;
 double *param=min_parameters->transfres->param;
 ERREUR_RECALAGE *err=min_parameters->err;
 int nbParamUtiles=0;
 double *prec_param=min_parameters->prec_param;


 for (i=0;i<nb_param;i++)
 {
  if (prec_param[i]!=0.0) nbParamUtiles++;
 }

 //Allocation des tableaux de parametres
 dTabParam = CALLOC(nbParamUtiles+1, double*);
 dTabDirection = CALLOC(nbParamUtiles, double*);
 if ((!dTabParam)||(!dTabDirection))
 {
  RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans min_powel_3d\n");
 }
 for (i=0;i<nbParamUtiles+1;i++)
 {
  dTabParam[i] = CALLOC(nb_param, double);
  if (!dTabParam[i])
  {
   RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans min_powel_3d\n");
  }
 }
 for (i=0;i<nbParamUtiles;i++)
 {
  dTabDirection[i] = CALLOC(nb_param, double);
  if (!dTabDirection[i])
  {
   RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans min_powel_3d\n");
  }
 }

 //initialisation des tableaux de parametres
 for (i=0;i<nbParamUtiles+1;i++)
  for (j=0;j<nb_param;j++)
   dTabParam[i][j]=param[j];

 l=0;
 for(i=0;i<nb_param;i++)
  if (prec_param[i]!=0.0)
  {
   for (j=0;j<nb_param;j++)
   {
    dTabDirection[l][j] = (i==j) ? 1.0 :0.0; 
   }
  l++;
  }

 //minimisation par la methode de powel
 retValue = PowellOptimization(min_parameters,dTabParam,dTabDirection);
 if (erreur_critique(err)) goto end_func;

 //modification de parametres initiaux pour le retour
 for (i=0;i<nb_param;i++) param[i]=dTabParam[nbParamUtiles][i];

 retValue=energie_transfo_3d(min_parameters);
 if (erreur_critique(err)) goto end_func;

 if (err) *err=NO_ERR;

end_func:

 for (i=0;i<nbParamUtiles;i++)
 {
  if (dTabParam[i]) FREE(dTabParam[i]);
  if (dTabDirection[i]) FREE(dTabDirection[i]);
 }
 if (dTabParam[nbParamUtiles+1]) FREE(dTabParam[nbParamUtiles+1]);
 if (dTabParam) FREE(dTabParam);
 if (dTabDirection) FREE(dTabDirection);

 return retValue;
}
/*! @} */

/*! \ingroup     MinimisationFonction  @{ */
/*******************************************************************************
**     min_powel2_3d(imref,imreca,imres,champ_ini,champ_fin,the_transf,inter
**                   ,dist,nb_param,min_param,max_param,prec_param,param)
*/
/*!     Minimisation par methode Powell telle que decrite dans
**      "Numerical Recipes in C"
**      http://www.library.cornell.edu/nr/bookcpdf.html
**
**     utilisation:
**     \param       imref,imreca,imres: images
**     \param       champ_ini : champ initial
**                              NULL = pas de champ initial
**     \param       champ_fin : champ final resultat de la transfo
**     \param       the_transf : pointeur sur la fonction de transformation
**                           int the_transf(nb_param,param,champ)
**     \param       inter : pointeur sur la fonction d'interpolation
**                           int inter(imdeb,champ,imres)
**     \param       dist : pointeur sur la fonction de distance
**                         double dist(im1,im2)
**     \param       nb_param: nombre de parametre a minimiser
**     \param       min_param: bornes inferieures des parametres
**     \param       max_param: bornes superieures des parametres
**     \param       prec_param: precision demandee pour chaque parametre
**     \param       param: resultats donnant la valeur des parametres au
**                         point definissant le minimum
**     \retval  la fonction retourne un double contenant la valeur du minimum
**      Les parametres correspondant a ce minimum sont retourne dans param
*******************************************************************************/
double min_powel2_3d(ptr_minimisation_param min_parameters)
{
  int i, ibig, j, iter=0, n_max=0;
  double del, t, *pt=NULL, *ptt=NULL, *xit=NULL;
  double energie1=0.0, energie2=0.0, energie_min=DBL_MAX, energie_deb=0.0;
  int nb_iter_max=NB_ITER_MAX_POWELL_3D;
  double **xi=NULL;
  double minparam, maxparam;
  double *old_param=NULL;

 int nb_param=min_parameters->transfres->nb_param;
 double *min_param=min_parameters->min_param;
 double *max_param=min_parameters->max_param;
 double *prec_param=min_parameters->prec_param;
 double *param=min_parameters->transfres->param;
 ERREUR_RECALAGE *err=min_parameters->err;

  ////allocations memoire et initialisations////
  pt=CALLOC(nb_param,double); 
  ptt=CALLOC(nb_param,double);
  xit=CALLOC(nb_param,double);
  old_param=CALLOC(nb_param,double);
  xi=alloc_dmatrix(nb_param,nb_param);
	if ((!pt)||(!ptt)||(!xit)||(!old_param)||(!xi))
	{
		RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans min_powel2_3d\n");
	}

  for (i=0;i<nb_param;i++) pt[i]=param[i];

  //determination de l'indice max des parametres qu'on peut faire varier
  for (i=0;i<nb_param;i++)
   if (prec_param[i]!=0.0) n_max=i;

  //construction des directions initiales qui sont ici la base de l'espace des parametres
  //si la precision d'un parametre est mise a 0, on prend pour ce parametre sa valeur dans le vecteur param initial
  for (i=0;i<nb_param;i++)
   for (j=0;j<nb_param;j++)
   {
    if (prec_param[i]!=0.0)
    {
     if (i==j)
      xi[i][j]=1.0;
     else
      xi[i][j]=0.0;
    }
    else
     xi[i][j]=param[j];
   }

  //calcul de l'energie initiale
  min_parameters->transfres->param=pt;
  energie_deb=energie1=energie_transfo_3d(min_parameters);
  min_parameters->transfres->param=param;
  if (erreur_critique(err)) goto end_func;
  aff_log("Edebut: %.8f     nb_param: %d \n",energie1,nb_param);

  do
  {
    energie_min=energie1;
    ibig=0;
    del=0.0;

    //on recentre l'intervalle de recherche
    for (j=0;j<nb_param;j++)
    {
     if (prec_param[j]!=0.0)
     {
      minparam=min_param[j];
      maxparam=max_param[j];
      min_param[j]=param[j]-(maxparam-minparam)/2.0;
      max_param[j]=param[j]+(maxparam-minparam)/2.0;
      old_param[j]=param[j];
     }
    }

    //minimisation lineaire selon chaque parametre
    for (i=0;i<nb_param;i++)
     if (prec_param[i]!=0.0)
     {//minimisation suivant le parametre d
      for (j=0;j<nb_param;j++) xit[j]=xi[i][j];
      energie2=energie1;
      energie1=linemin_3d(min_parameters,xit);
      if (erreur_critique(err)) goto end_func;
      //on determine la direction pour laquelle la decroissance de l'energie est maximale
      if ((energie2-energie1)>del)
      { del=energie2-energie1; ibig=i; }
     }

    printf("nbiter:%d   Energie: %.8f  %%reduc: %.2f  \r",iter,energie1,100.0*(energie_deb-energie1)/fabs(energie_deb));

    //si les parametres n'ont pas bouge de plus que la precision recherchee, on s'arrete
    if (arret_selon_parametres(param,old_param,prec_param,nb_param)) break;

    //calcul de la direction moyenne et de la direction extrapolee
    for (j=0;j<nb_param;j++)
    {
     //direction moyenne dans laquelle on s'est deplace lors de la minimisation
     xit[j]=param[j]-pt[j];
     //direction extrapolee en se deplacant de 2 fois la direction moyenne
     ptt[j]=2.0*param[j]-pt[j];
     //sauvegarde des parametres calcules lors de la minimisation
     pt[j]=param[j];
    }

    //calcul de l'energie obtenue pour la direction extrapolee
    min_parameters->transfres->param=ptt;
    energie2=energie_transfo_3d(min_parameters);
    min_parameters->transfres->param=param;
    if (erreur_critique(err)) goto end_func;

    if (energie2<energie_min)
    {
     t=2.0*(energie_min-2.0*energie1+energie2)*sqrt(energie_min-energie1-del)-del*sqrt(energie_min-energie2);
     if (t<0.0)
     {
      //on cherche le minimum dans la nouvelle direction
      energie1=linemin_3d(min_parameters, xit);
      if (erreur_critique(err)) goto end_func;
      //et on la sauvegarde
      for (j=0;j<nb_param;j++)
      {
       xi[ibig][j]=xi[n_max][j];
       xi[n_max][j]=xit[j];
       pt[j]=param[j];
      }
     }
    }

    iter ++;

  } while (iter<nb_iter_max);

  energie_min=energie1;
  aff_log("nbiter:%d   Energie: %.8f  %%reduc: %.2f  \n",iter,energie_min,100.0*(energie_deb-energie_min)/fabs(energie_deb));

 if (err) *err=NO_ERR;

end_func:

  ////allocations memoire////
  if (old_param) FREE(old_param);
  if (pt) FREE(pt);
  if (ptt) FREE(ptt);
  if (xit) FREE(xit);
  if (xi) free_dmatrix(xi,nb_param,nb_param);

  return energie_min;
}
/*! @} */

/*! \ingroup     MinimisationFonction  @{ */
/*******************************************************************************
**     min_simplex_robuste_3d(imref,imreca,imres,champ_ini,champ_fin,the_transf,inter
**                   ,dist,nb_param,min_param,max_param,prec_param,param)
**
*/
/*!    Minimisation par methode simplex iterative robuste: implementation simplex telle
**     decrite dans le "numerical recipes in C". plusieurs minimisation simplex sont faites
**     avec a chaque fois une mise a jour des poids utilises pour le calcul de la distance robuste
**
**     utilisation:
**     \param       imref,imreca,imres: images
**     \param       champ_ini : champ initial
**                              NULL = pas de champ initial
**     \param       champ_fin : champ final resultat de la transfo
**     \param       the_transf : pointeur sur la fonction de transformation
**                           int the_transf(nb_param,param,champ)
**     \param       inter : pointeur sur la fonction d'interpolation
**                           int inter(imdeb,champ,imres)
**     \param       dist : pointeur sur la fonction de distance
**                         double dist(im1,im2)
**     \param       nb_param: nombre de parametre a minimiser
**     \param       min_param: bornes inferieures des parametres
**     \param       max_param: bornes superieures des parametres
**     \param       prec_param: precision demandee pour chaque parametre
**     \param       param: resultats donnant la valeur des parametres au
**                         point definissant le minimum
**     \retval       la fonction retourne un double contenant la valeur du minimum
**      Les parametres correspondant a ce minimum sont retourne dans param
*******************************************************************************/ 
double min_simplex_robuste_3d(ptr_minimisation_param min_parameters)
{
 int i, ihi=0, ilo=0, inhi=0, j, Nest, iter=0;
 double /*rtol,*/ Esav, Ereflexion, Eopt=DBL_MAX, Edebut;
 double *psum=NULL, *energies=NULL, *p_reflexion=NULL;
 double **simplexe=NULL;
 double minparam, maxparam;
 double *p_deb=NULL;
 int wdth, hght, dpth;

 grphic3d *imref=min_parameters->imref;
 dist_func_t dist=min_parameters->dist;
 int nb_param=min_parameters->transfres->nb_param;
 double *min_param=min_parameters->min_param;
 double *max_param=min_parameters->max_param;
 double *prec_param=min_parameters->prec_param;
 double *param=min_parameters->transfres->param;
 ERREUR_RECALAGE *err=min_parameters->err;

 wdth=imref->width; hght=imref->height; dpth=imref->depth;
   
 Nest=1;
  for (i=0;i<nb_param;i++)
   if (prec_param[i]!=0.0) {Nest++;}

 if (Nest<4) 
 {
  RAISE_ERR_TEXT(err, PB_PARAMETRES, "Trop peu de parametres a estimer pour la methode simplex robuste\n");
 }

 //INITIALISATION
 //allocation memoire des vecteurs et matrice
 simplexe=alloc_dmatrix(Nest,nb_param);
 psum=CALLOC(nb_param, double);
 energies=CALLOC(Nest, double);
 p_reflexion=CALLOC(nb_param, double);
 p_deb=CALLOC(nb_param, double);
	if ((!simplexe)||(!psum)||(!energies)||(!p_reflexion)||(!p_deb))
	{
		RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans min_simplex_robuste_3d\n");
	}

 //on fait n minimisations simplexe avec un update des poids du robuste apres chaque minimisation
 do {

 iter=0;

 //on recentre l'intervalle de recherche
 for (j=0;j<nb_param;j++)
    {
     if (prec_param[j]!=0.0)
     {
      minparam=min_param[j];
      maxparam=max_param[j];
      min_param[j]=param[j]-(maxparam-minparam)/2.0;
      max_param[j]=param[j]+(maxparam-minparam)/2.0;
      p_deb[j]=param[j];
     }
    }
      
 //remplissage du simplex initial
 //mise a zero du simplex
 for (i=0;i<Nest;i++)
  for (j=0;j<nb_param;j++)
     simplexe[i][j]=0.0;
 //le point P0
 for (i=0;i<nb_param;i++)
  for (j=0;j<Nest;j++)
    simplexe[j][i]=param[i];

 min_parameters->transfres->param=simplexe[0];
 energies[0]=energie_transfo_3d(min_parameters);
 min_parameters->transfres->param=param;

 //les points 1 a N
 //on va prendre pour chaque point max_param ou min_param selon la direction de moindre energie
 j=1;
 for (i=0;i<nb_param;i++)
 {
  if (prec_param[i]!=0.0)
  {
   simplexe[j][i]=max_param[i];
   min_parameters->transfres->param=simplexe[j];
   energies[j]=energie_transfo_3d(min_parameters);
   min_parameters->transfres->param=param;
   if (erreur_critique(err)) goto end_func;
   Esav=energies[j];
   simplexe[j][i]=min_param[i];
   min_parameters->transfres->param=simplexe[j];
   energies[j]=energie_transfo_3d(min_parameters);
   min_parameters->transfres->param=param;
   if (erreur_critique(err)) goto end_func;
   if (energies[j]>Esav)
   {
    simplexe[j][i]=max_param[i];
    energies[j]=Esav;
   }

   j++;
  }
 }

 //calcul des sommes des valeurs des points du simplexe pour chaque parametre
 get_simplexe_psum(simplexe, psum, Nest, nb_param);

 //calcul et affichage de l'energie de depart
 Edebut=energie_transfo_3d(min_parameters);
 if (erreur_critique(err)) goto end_func;
 aff_log("Edebut: %.8f     nb_param: %d \n",Edebut,nb_param);

 //MINIMISATION
 while (iter<NB_ITER_MAX_SIMPLEX_3D)
 {
  ilo=0;
  //determination du plus mauvais, du deuxieme plus mauvais et du meilleur point
  if (energies[0]>energies[1])
   { ihi=0; inhi=1; }
  else
   { ihi=1; inhi=0; }
  for (i=0;i<Nest;i++)
  {
   if (energies[i]<energies[ilo]) ilo=i;
   if (energies[i]>energies[ihi])
    { inhi=ihi; ihi=i; }
   else if ((energies[i]>energies[inhi])&&(i!=ihi))
    inhi=i;
  }

  printf("iteration: %d  energie:%.8f  reduction : %.2f %% \r",iter,energies[ilo],(Edebut-energies[ilo])/fabs(Edebut)*100.0);

  //critere d'arret
  if (arret_selon_parametres(simplexe[ilo], simplexe[ihi], prec_param, nb_param)) break;
//  if (arret_selon_energie(energies[ilo], energies[ihi])) break;

  //nouvelle iteration
  min_parameters->transfres->param=p_reflexion;
  Ereflexion=reflexion_simplex_3d(min_parameters, simplexe, energies, psum, ihi, Nest, -1.0);
  min_parameters->transfres->param=param;
  if (erreur_critique(err)) goto end_func;

  //si la reflexion donne un meilleur resultat que le meilleur des points
  //on continue dans cette direction par un facteur 2
  if (Ereflexion<energies[ilo])
  {
   min_parameters->transfres->param=p_reflexion;
   Ereflexion=reflexion_simplex_3d(min_parameters, simplexe, energies, psum, ihi, Nest, 2.0);
   min_parameters->transfres->param=param;
   if (erreur_critique(err)) goto end_func;
  }
  //si on a un moins bon resultat que le deuxieme plus mauvais point
  //on teste une reflexion intermediaire
  else if (Ereflexion>=energies[inhi])
  {
   Esav=energies[ihi];
   min_parameters->transfres->param=p_reflexion;
   Ereflexion=reflexion_simplex_3d(min_parameters, simplexe, energies, psum, ihi, Nest, 0.5);
   min_parameters->transfres->param=param;
   if (erreur_critique(err)) goto end_func;
   //on n'a toujours pas ameliore le score du plus mauvais point:
   //on contracte le simplexe autour du meilleur point
   if (Ereflexion>=Esav)
   {
    for (i=0;i<Nest;i++)
    {
     if (i!=ilo)
     {
      for (j=0;j<nb_param;j++)
       simplexe[i][j]=psum[j]=0.5*(simplexe[i][j]+simplexe[ilo][j]);

      min_parameters->transfres->param=psum;
      energies[i]=energie_transfo_3d(min_parameters);
      min_parameters->transfres->param=param;
      if (erreur_critique(err)) goto end_func;
     }
    }
    get_simplexe_psum(simplexe, psum, Nest, nb_param);
   }

  }
  iter++;

 }

 //determination du meilleur point si on s'est arrete avant la convergence
 if (iter>=NB_ITER_MAX_SIMPLEX_3D)
 {
  ilo=0;
  for (i=0;i<Nest;i++)
   if (energies[i]<energies[ilo]) ilo=i;
 }

 for (i=0;i<nb_param;i++)  param[i]=simplexe[ilo][i];
 Eopt=energies[ilo]; 
 printf("iteration: %d  energie:%.8f  reduction : %.2f %% \n",iter,Eopt,(Edebut-Eopt)/fabs(Edebut)*100.0);

 //si la minimisation n'a pas fait avancer les choses, on s'arrete
 if (arret_selon_parametres(param, p_deb, prec_param, nb_param)) break;
 
 imx_aff_param(RIGID3D, param);

 //on recalcule le sigma dans le cas du quadratique, du woods et du ratio de correlation robustes, on le met a une valeur arbitraire dans le cas contraire
 if (dist==erreur_quad_robust2_3d) compute_robust_weights_quad(min_parameters);
 else if ((dist==erreur_woods_robust2_3d)||(dist==erreur_CR_robust_3d)) compute_robust_weights_woods(min_parameters);
 else _sigma_robust=80.0;
 
 } while (TRUE);

 if (err) *err=NO_ERR;

end_func:

 //liberation memoire
 if (p_reflexion) FREE(p_reflexion);
 if (energies) FREE(energies);
 if (psum) FREE(psum);
 if (simplexe) free_dmatrix(simplexe,Nest, nb_param);
 if (p_deb) FREE(p_deb);

 return Eopt;
}
/*! @} */

// ================================================================================================

/*! \brief Powell method

*   \param imref         : image de reference
*   \param imreca        : image a recaler
*   \param imres         : image resultat
*   \param champ_ini     :
*   \param champ_fin
*   \param the_transf    : fonction de transformation
*   \param inter         : fonction d'interpolation
*   \param dist          : fonction de cout
*   \param contrainte    : contrainte (NULL)
*   \param nNbParam      : nombre de params a minimiser
*   \param min_param     : bornes inf des parametres
*   \param max_param     : bornes sup des parametres
*   \param prec_param    : precision sur les parametres
*   \param dTabParam     : valeurs des parametres
*   \param dTabDirection : vecteur de direction de recherche initial
*   \param donnees_VP    : pour l'interpolation VP?
*/
// ================================================================================================
double PowellOptimization (ptr_minimisation_param min_parameters, double **dTabParam,double **dTabDirection)
{
  int      i=0, nCurrentParam=0, nPowellIter=0;
  int      s_bas = 0;
  bool     bParamChange;
  double * dPreviousParam = NULL;
  double * dLambda = NULL;
  int nb=0;
  double ** pdVecProduct = NULL;
  double   dNorme = 0.0;
  double *p;
  double Eini=0.0,E1=0.0,E2=0.0;
  double max_p, min_p;

  int nNbParam=min_parameters->transfres->nb_param;
  double *param=min_parameters->transfres->param;
  double *max_param=min_parameters->max_param;
  double *min_param=min_parameters->min_param;
  double *prec_param=min_parameters->prec_param;
  ERREUR_RECALAGE *err=min_parameters->err;
  int nbParamUtiles=0;

  for (i=0;i<nNbParam;i++)
  {
   if (prec_param[i]!=0.0) nbParamUtiles++;
  }

  bParamChange = TRUE;

  //allocations memoire
  dPreviousParam = CALLOC(nNbParam, double);
  pdVecProduct = CALLOC((nbParamUtiles-1),double*);
  dLambda = CALLOC(nbParamUtiles, double);
  p=CALLOC(nNbParam,double);

  if ((dPreviousParam == NULL)||(pdVecProduct == NULL)||(dLambda == NULL)||(p == NULL))
  {
		RAISE_ERR_TEXT(err, ERREUR_ALLOC_MEMOIRE, "pb d'allocation memoire dans min_powel2_3d\n");
  }

  for (i=0; i<nNbParam; i++) dPreviousParam[i] = 0.0;

  //calcul de l'energie de depart
  min_parameters->transfres->param=dTabParam[0];
  Eini=energie_transfo_3d(min_parameters);
  min_parameters->transfres->param=param;
  E2=Eini;
  aff_log("Edebut: %.2f     nb_param: %d \n",Eini,nNbParam);

// Powell optimization start
// =========================

  for (nPowellIter = 0; (nPowellIter<nNbParam) && bParamChange; nPowellIter++)
  {
    E1=E2;
    s_bas = 0;
    for (nCurrentParam = 0; nCurrentParam < nbParamUtiles; nCurrentParam++)
    {
     //recentrage de l'intervalle de recherche
     for (i=0;i<nNbParam;i++)
     {
      min_p=min_param[i];
      max_p=max_param[i];
      min_param[i]=p[i]-(max_p-min_p)/2.0;
      max_param[i]=p[i]+(max_p-min_p)/2.0;
      dTabParam [nCurrentParam + 1][i]=dTabParam[nCurrentParam][i];
     }

     //optimisation dasn la direction consideree
     nb++;
     min_parameters->transfres->param=dTabParam [nCurrentParam + 1];
     E2 = linemin_3d(min_parameters, dTabDirection [nCurrentParam]);
     min_parameters->transfres->param=param;

     printf("nbiter:%d   Energie: %.2f  %%reduc: %.2f   param : %d\r",nb,E2,100.0*(E1-E2)/E1,nCurrentParam);

    if (!arret_selon_parametres(dTabParam [nCurrentParam + 1],dTabParam [nCurrentParam],prec_param,nNbParam))
       s_bas ++ ;
     else
       ; // do nothing    

    } //end for nCurrentParam

    if (!s_bas)
    {
      break;
    }
    else
    {
      ; // do nothing
    }

    // Shift the search directions for the next iteration
    // ==================================================

    for (nCurrentParam = 0; nCurrentParam < (nbParamUtiles - 1); nCurrentParam++)
      for (i=0; i<nNbParam; i++)
        dTabDirection [nCurrentParam][i] = dTabDirection [nCurrentParam + 1][i];

    // New search direction
    // ====================

    if ((s_bas == 1) && (arret_selon_parametres(dTabParam [nbParamUtiles],dTabParam [nbParamUtiles-1],prec_param,nNbParam)))
    {
      printf( "PowellOptimization --> Bad directions\n");

      for (i=0; i<(nbParamUtiles - 1); i++) pdVecProduct [i] = dTabDirection [i];

      PowellVectorProduct (pdVecProduct, (nbParamUtiles - 1), dTabDirection [nbParamUtiles - 1]);
      // end bad directions
    }
    else
    {
      for (i=0; i<nNbParam; i++)
      {
        dTabDirection [nbParamUtiles - 1][i] = dTabParam [nbParamUtiles][i] - dTabParam [0][i];
      }
    }// end computation of the new search direction



    // normalization of the direction
    // ==============================

    dNorme = 0.0;

    for (i=0; i<nNbParam; i++)
    {
      dNorme += (dTabDirection [nbParamUtiles - 1][i] * dTabDirection [nbParamUtiles - 1][i]);
    }

    dNorme = sqrt (dNorme);

    bParamChange = FALSE;

    if (dNorme > POWELL_TINY)
    {
      for (i=0; i<nNbParam; i++)
      {
        dTabDirection [nbParamUtiles - 1][i] = dTabDirection [nbParamUtiles -1][i] / dNorme;
      }


      // optimization along the new direction
      // ------------------------------------

      //recentrage de l'intervalle de recherche
      for (i=0;i<nNbParam;i++)
      {
       min_p=min_param[i];
       max_p=max_param[i];
       min_param[i]=p[i]-(max_p-min_p)/2.0;
       max_param[i]=p[i]+(max_p-min_p)/2.0;
       dTabParam [0][i] = dTabParam [nbParamUtiles][i];
      }
      min_parameters->transfres->param=dTabParam [0];
      E2 = linemin_3d(min_parameters, dTabDirection [nbParamUtiles - 1]);
      min_parameters->transfres->param=param;

      for (i=0; (i<nNbParam); i++)
      {
        bParamChange     |= (dTabParam [0][i] != dPreviousParam[i]);
        dPreviousParam[i] = dTabParam [0][i];
      }

    } // end if (dNorme > POWELL_TINY)
    else
    {
      // Norm too small. We stop Powell. bParamChange is already FALSE.
      // --------------------------------------------------------------

      printf( "\tNorme = 0.001 => exit...\n");

      bParamChange = FALSE;
    } // end dNorme == 0

  } // end for nPowellIter

  printf("nbiter:%d   Energie: %.2f  %%reduc: %.2f   param : %d\n",nb,E2,100.0*(E1-E2)/E1,nCurrentParam);

end_func :

 if (dPreviousParam != NULL)
  {
    free(dPreviousParam);
    dPreviousParam = NULL;
  }

  if (pdVecProduct != NULL)
  {
    free(pdVecProduct);
    pdVecProduct = NULL;
  }

  if (dLambda) free(dLambda);
  if (p) free(p);

  return 0.0;
}//end PowellOptimization

// ================================================================================================

/*! \brief Vector product for Powell method.

*   \param dVectSrc     Array of source vectors.
*
*   \param nNbVect      Number of vectors. These vectors are assumed to be \a (nNbVect+1) dimension.
*
*   \param dVectResult  Pointer to a \(nNbVect+1) dimension vector result.
*                       The memory should be already allocated
*
*   This function will process a vector product of \a nNbVect vectors of \a (nNbVect+1) dimension.
*   The vector result will be stored in \a dVectResult.
*/
// ================================================================================================


void PowellVectorProduct (double ** dVectSrc, int nNbVect, double * dVectResult)  {
	double ** dSubmatrix  = NULL;
	int    * nSignTab    = NULL;
	int		   i, j, k;
	int      nNumVect;
	int      nVectDimension = nNbVect + 1;

	dSubmatrix  = (double**)malloc((nNbVect)*sizeof(double*));
	for (i=0;i<nNbVect;i++)
		dSubmatrix[i]=(double*)malloc((nNbVect)*sizeof(double));
	nSignTab    = (int*)malloc(nVectDimension*sizeof(int));


	// Initialize nSignTab
	for (i=0; i<nVectDimension; i++)  {
		if (i & 0x0001)  {
			nSignTab [i] = -1;
		} else {
			nSignTab [i] = 1;
		}
	}


	for (i=0; i<nVectDimension; i++)
	{
		k = 0;

		for (j=0; j<nVectDimension; j++)
		{
			if (j != i)
			{
				for (nNumVect=0; nNumVect<nNbVect; nNumVect++)
				{
					dSubmatrix [nNumVect][k]  = dVectSrc[nNumVect][j];
				}
				k++;
			}
		}
		dVectResult[i] = Determinant(dSubmatrix,nNbVect);
		dVectResult[i] = nSignTab[i] * dVectResult[i];
	 }


	if (dSubmatrix != NULL)  {
		for (i=0;i<nNbVect;i++)
			free(dSubmatrix[i]);
		free(dSubmatrix);
		dSubmatrix = NULL;
	}

	if (nSignTab != NULL)  {
		free(nSignTab);
		nSignTab = NULL;
	}
}

/*******************************************************************************
**     get_simplexe_psum(simplexe, psum, Nest, nb_param)
*/
/*!    calcule les sommes des valeurs de chaque point pour tous les parametres
**     du simplexe
**
**
**     utilisation:
**     \param       simplexe: le tableau contenant le simplexe
**     \param       psum : le vecteur alloue ou seront stockees les sommes
**     \param       Nest : le nombre de points dans le simplexe
**     \param       nb_param : le nombre de parametres dans le type de transformation utilise
*******************************************************************************/
void get_simplexe_psum(double **simplexe, double *psum, int Nest, int nb_param)
{
 int i,j;

 for (i=0;i<nb_param;i++)
 {
  psum[i]=0.0;
  for (j=0;j<Nest;j++)
   psum[i]+=simplexe[j][i];
 }

 return;
}

/*******************************************************************************
**     reflexion_simplex_3d(imref, imreca, imres, champ_ini, champ_fin,
**                        the_transf, inter, dist, contrainte, nb_param,
**                        prec_param, donnees_VP, simplexe, energies, psum, ihi, Nest, fac)
*/
/*!    calcul le point reflechi du plus mauvais point par rapport a l'hyperplan forme
**     par les autres points
**
**
**     utilisation:
**     \param       imref,imreca,imres: images
**     \param       champ_ini : champ initial
**                              NULL = pas de champ initial
**     \param       champ_fin : champ final resultat de la transfo
**     \param       the_transf : pointeur sur la fonction de transformation
**     \param       inter : pointeur sur la fonction d'interpolation
**     \param       dist : pointeur sur la fonction de distance
**     \param       contrainte : pointeur sur la fonction de contrainte
**     \param       prec_param : precision demandee pour chaque parametre
**     \param       donnees_VP : structure utilisee pour les distances se
**                               basant sur le volume partiel
**     \param       simplexe : le tableau represant le simplexe courant
**     \param       energies : vecteur contenant l'energie associe a chaque point du simplexe
**     \param       psum : vecteur contenant les sommes de valeurs du simplexe pour chaque parametre
**     \param       ihi : indice associe au plus mauvais point du simplexe
**     \param       Nest : nombre de points du simplexe
**     \param       fac : facteur de reflexion
**     \retval      retourne l'energie associee au point reflechi
*******************************************************************************/
double reflexion_simplex_3d(ptr_minimisation_param min_parameters,
                        double ** simplexe, double *energies,
                        double *psum, int ihi, int Nest, double fac)
{
 int j;
 double fac1, fac2, Ereflexion;
 double *p_reflexion=NULL;
 int nb_param=min_parameters->transfres->nb_param;

 double *prec_param=min_parameters->prec_param;

 p_reflexion=min_parameters->transfres->param;

 //calcul du parametre reflechi  du plus mauvais point du simplexe
 //au travers de l'hyperplan forme par les autres points
 //et de l'energie associee
 fac1=(1.0-fac)/(Nest-1);
 fac2=fac1-fac;
 for (j=0;j<nb_param;j++)
 {
  if (prec_param[j]!=0.0)
   p_reflexion[j]=psum[j]*fac1-simplexe[ihi][j]*fac2;
  else
   p_reflexion[j]=simplexe[ihi][j];
 }

 Ereflexion=energie_transfo_3d(min_parameters);

 //si on a un meilleur point, on remplace le plus mauvais
 if (Ereflexion<energies[ihi])
 {
  energies[ihi]=Ereflexion;
  for (j=0;j<nb_param;j++)
   if (prec_param[j]!=0.0)
   {
    psum[j]+=p_reflexion[j]-simplexe[ihi][j];
    simplexe[ihi][j]=p_reflexion[j];
   }
 }

 return Ereflexion;
}

/*******************************************************************************
**     compute_robust_weights_quad(imref, imreca, imres, transf, inter, champ_ini, champ_fin, nb, param, weights
*/
/*!    calcul les poids associes au calcul robuste de l'erreur quadratique
**     selon la methode de geman-mc lure
**
**     utilisation:
**     \param       imref,imreca,imres: images ref et reca et reca transformee
**     \param       transf, inter, champ_ini, champ_fin, nb, param: parametres servant a calculer l'image transformee
**     \param       weights : tableau des poids calcules alloue
*******************************************************************************/
void compute_robust_weights_quad(ptr_minimisation_param min_parameters)
{
 int i,j,k;
 int wdth,hght,dpth;
 double diff, moy=0.0, ect=0.0, nb_voxels;
 TYPEMRI3D m1, m2;
 double sigma, diff2, sigma2;
 TYPEMRI3D ***imrefMRI=NULL, ***imresMRI=NULL;

 grphic3d *imref=min_parameters->imref;
 grphic3d *imreca=min_parameters->imreca;
 grphic3d *imres=min_parameters->imres;
 field3d *champ_ini=min_parameters->champ_ini;
 field3d *champ_fin=min_parameters->champ_fin;
 InterpolationFct inter=min_parameters->inter;
 double ***weights=min_parameters->dist_param->weights;

 wdth=imref->width; hght=imref->height; dpth=imref->depth;

 nb_voxels=(double)(wdth*hght*dpth);

 transf_to_field_3d_noalloc(min_parameters->transfres, champ_fin, imref, imreca);
 
 //au champ final on rajoute le champ initial
 if (champ_ini!=NULL)
 {
  add_field_3d(champ_fin,champ_ini,champ_fin);
 }

 (*inter)(imreca,champ_fin,imres);

 imresMRI=imres->mri; imrefMRI=imref->mri;
 //on calcule la variance de l'image difference: imreca transformee moins imref
 //cette variance est assimilee a la variance de l'erreur residuelle et servira
 //de valeur de rejection dans le calul des poids
 for (i=0;i<wdth;i++)
 {
  for (j=0;j<hght;j++)
  {
   for (k=0;k<dpth;k++)
   {
    diff=fabs(imresMRI[i][j][k]-imrefMRI[i][j][k]);
    moy += diff;
    ect += diff*diff;
   }
  }
 }

 moy=moy/nb_voxels;
 ect=ect/nb_voxels-moy*moy;

 //calcul des poids robustes par la formule de geman-mc lure
 sigma=ect;

 for (i=0;i<wdth;i++)
 {
  for (j=0;j<hght;j++)
  {
   for (k=0;k<dpth;k++)
   {
    m1=imrefMRI[i][j][k]; m2=imresMRI[i][j][k];
    diff=(double)(m1-m2); diff2=diff*diff;
    sigma2=sigma*sigma;
    weights[i][j][k]=(float)2.0*sigma2/((sigma2+diff2)*(sigma2+diff2));
   }
  }
 }

}

/*******************************************************************************
**     compute_robust_weights_woods(imref, imreca, imres, transf, inter, champ_ini, champ_fin, nb, param, weights
*/
/*!    calcul les poids associes au calcul robuste du critere de woods (sert aussi pour le ratio de correlation)
**     selon la methode de geman-mc lure
**
**     utilisation:
**     \param       imref,imreca,imres: images ref et reca et reca transformee
**     \param       transf, inter, champ_ini, champ_fin, nb, param: parametres servant a calculer l'image transformee
**     \param       weights : tableau des poids calcules alloue
*******************************************************************************/
void compute_robust_weights_woods(ptr_minimisation_param min_parameters)
{
 int i,j,k,wdth,hght,dpth;
 int nbcl1,nbcl2,max1,max2,cl1,cl2;
 int nbVoxelsCommuns=0;
 double *moyennesParClasse, *variancesParClasse;
 int *nbVoxelsParClasse;
 double **region;
 double sigma2, diff, diff2, w;

 grphic3d *imref=min_parameters->imref;
 grphic3d *imreca=min_parameters->imreca;
 grphic3d *imres=min_parameters->imres;
 field3d *champ_ini=min_parameters->champ_ini;
 field3d *champ_fin=min_parameters->champ_fin;
 InterpolationFct inter=min_parameters->inter;
 double ***weights=min_parameters->dist_param->weights;

 /////calcul de l'image transformee
 transf_to_field_3d_noalloc(min_parameters->transfres, champ_fin, imref, imreca);

 //au champ final on rajoute le champ initial
 if (champ_ini!=NULL)
 {
  add_field_3d(champ_fin,champ_ini,champ_fin);
 }

 (*inter)(imreca,champ_fin,imres);

 wdth=imres->width;hght=imres->height;dpth=imres->depth;
 //calcul des maxpixel dans les deux images
 //et du nombre de classes de valeurs de voxels qu'on va utiliser
// nbcl1=imx_calculer_nombre_de_classes_3d_p(imres, &max1);
// nbcl2=imx_calculer_nombre_de_classes_3d_p(imref, &max2);

 nbcl1=min_parameters->dist_param->nbclReca; nbcl2=min_parameters->dist_param->nbclRef;
 max1=min_parameters->dist_param->maxReca; max2=min_parameters->dist_param->maxRef;

 region=CALLOC(nbcl2,double *);
 for (i=0;i<nbcl2;i++) region[i]=alloc_dvector(1);
 
 //allocations memoires et initialisation
 nbVoxelsParClasse = CALLOC(nbcl2,int);
 moyennesParClasse = CALLOC(nbcl2,double);
 variancesParClasse = CALLOC(nbcl2,double);

 //partition de imres en region qui correspondent a des intensite homogenes dans imref
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    //on fait une classe a part pour les intensites nulles qui ne sont porteuses d'aucune information
    if (imref->mri[i][j][k]<=0) cl2=0;
    else cl2=MINI(((nbcl2-1)*imref->mri[i][j][k])/(max2+1)+1, (nbcl2-1));
    if (imres->mri[i][j][k]<=0) cl1=0;
    else cl1=imres->mri[i][j][k];
    //on n'utilise que les informations provenant de la zone de recouvrement des images
    if ((cl2>0)&&(cl1>0))
    {
     region[cl2][nbVoxelsParClasse[cl2]]=(double) cl1;
     nbVoxelsParClasse[cl2]++;
     nbVoxelsCommuns++;
     region[cl2]=realloc_dvector(region[cl2],(nbVoxelsParClasse[cl2]+1));
    }
   }

 //si il n'y a pas de recouvrement des 2 images on ne fait rien
 if (nbVoxelsCommuns==0) return ;

 //calcul des moyennes et des variances par classe
 //on ne fait le calcul que sur la zone de recouvrement des images
  for (i=1;i<nbcl2;i++)
  {
   if (nbVoxelsParClasse[i]>0)
   {
    //calul de la moyenne robuste (mediane)
    moyennesParClasse[i]=quick_select(nbVoxelsParClasse[i]/2,region[i],nbVoxelsParClasse[i]);
    for (j=0;j<nbVoxelsParClasse[i];j++)
    {
     diff=region[i][j]-moyennesParClasse[i];
     w=diff*diff;
     variancesParClasse[i]+=w;
    }
   }
  }

  for (i=1;i<nbcl2;i++)
  {
   if (nbVoxelsParClasse[i]>0)
   {
    variancesParClasse[i]=variancesParClasse[i]/(double)nbVoxelsParClasse[i];
   }
  }

  //calcul des poids
  for (i=0;i<wdth;i++)
   for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
    {
     //on fait une classe a part pour les intensites nulles qui ne sont porteuses d'aucune information
     //if ((im1->mri[i][j][k]<=0)||(im2->mri[i][j][k]<=0)) cl2=0;
     if (imref->mri[i][j][k]<=0) cl2=0;
     else cl2=MINI(((nbcl2-1)*imref->mri[i][j][k])/(max2+1)+1, (nbcl2-1));
     if (imres->mri[i][j][k]<=0) cl1=0;
     else cl1=imres->mri[i][j][k];
     if (cl2>0)
     {
      diff=(double)(cl1-moyennesParClasse[cl2]);
      diff2=diff*diff;
      //calcul du facteur de rejection:
      //ici ce facteur, qui est la variance des erreurs residuelles, est assimile a la variance correspondant a chaque region
      sigma2=variancesParClasse[cl2]*variancesParClasse[cl2];
      if (sigma2>0)
      {
       weights[i][j][k]=(float)2.0*sigma2/((sigma2+diff2)*(sigma2+diff2));
      }
     }
     else weights[i][j][k]=0.0;
    }
  
 //liberation memoire
 for (i=0;i<nbcl2;i++) free_dvector(region[i],(nbVoxelsParClasse[i]+1));
 FREE(region);
 FREE(nbVoxelsParClasse);
 FREE(moyennesParClasse);
 FREE(variancesParClasse);
}

/*******************************************************************************
**     tri_double_index(double *arr, int taille, int *index)
*/
/*!    permet de trier un tableau en ordre croissant a partir d'un tableau d'indexes
**
**     utilisation:
**     \param       arr: tableau a trier
**     \param       taille : taille du tableau
**     \param       index : tableau d'indexes
*******************************************************************************/
void tri_double_index(double *arr, int taille, int *index)
{
  int i,j,k;
  int temp;

  for (i=0;i<taille;i++)
  {
   for (j=0;j<i;j++)
   {
    if (arr[index[i]]<arr[index[j]])
    {
     temp=index[i];
     for (k=i;k>j;k--) index[k]=index[k-1];
     index[j]=temp;
     break;
    }
   }
  }

  return;
}

/*******************************************************************************
**        imx_query_minimisation_type_3D
*/
/*!
**       \brief choix de la fonction de minimisation
**
**		ouvre une boite de dialogue pour le choix de la fonction de minimisation
** 		\retval int : le numero de la fonction de cout
**			- ICM -> 0
**			- simplex -> 1
**			- Powel -> 2
**			- etc ...
**
*******************************************************************************/
int imx_query_minimisation_type_3D(int dist_type)
{
	char *query[100];
	query[0]="ICM";
	query[1]="simplex";
	query[2]="Powel";
  query[3]="Powel2";
    	if ((dist_type==0)||((dist_type>=11)&&(dist_type<=15)))
    	{
		query[2]="Descente gradient";
		query[3]="Gradient conjugue";
		query[4]="Descente gradient modifiee";
		query[5]="Quasi Newton modifiee";
		query[6]="Powel";
    query[7]="Powel2";
		query[8]="\0";
		query[9]=NULL;
	}
  else if ((dist_type==17)||(dist_type==18)||(dist_type==20))
    	{
		query[2]="Powel";
    query[3]="Powel2";
    query[4]="simplex robuste";
		query[5]="\0";
		query[6]=NULL;
	}
	else
	{
		query[4]="\0";
		query[5]=NULL;
	}

	return GETV_QCM("Minimisation",(char **)query);
}

/*------------------------------------------------*/
/*!
**	\brief determine la fonction de minimisation
**
**	\param  min_type : le type de minimisation choisie
**	\param  dist_type : le type de distance choisie
**	\retval min_func_t : un pointeur sur la fonction de minimisation
*/
/*------------------------------------------------*/
min_func_t
imx_choose_minimisation_fct(int min_type,int dist_type)
{
  min_func_t minimisation;
  if ((dist_type==0)||((dist_type>=11)&&(dist_type<=15)))
  {
    switch (min_type)
    {
      case 0: minimisation=min_icm_3d;aff_log("ICM ");break;
	  case 1: minimisation=min_simplex_3d;aff_log("SIMPLEX ");break;
	  case 2: minimisation=min_desc_grad_3d;aff_log("DESCGRAD ");break;
	  case 3: minimisation=min_grad_conj_3d;aff_log("GRADCONJ ");break;
	  case 4: minimisation=min_desc_grad_mod_3d;aff_log("DESCGRADMOD ");break;
	  case 5: minimisation=min_quasi_newton_mod_3d;aff_log("QUASINEWTONMOD ");break;
	  case 6: minimisation=min_powel_3d;aff_log("POWEL ");break;
    case 7: minimisation=min_powel2_3d;aff_log("POWEL2 ");break;
	  default:minimisation=min_simplex_3d;aff_log("SIMPLEX ");break;
    }
  }
  else if ((dist_type==17)||(dist_type==18)||(dist_type==20))
  {
    switch (min_type)
    {
    case 0: minimisation=min_icm_3d;aff_log("ICM ");break;
	  case 1: minimisation=min_simplex_3d;aff_log("SIMPLEX ");break;
	  case 2: minimisation=min_powel_3d;aff_log("POWEL ");break;
    case 3: minimisation=min_powel2_3d;aff_log("POWEL2 ");break;
    case 4: minimisation=min_simplex_robuste_3d;aff_log("SIMPLEX ROBUSTE ");break;
	  default:minimisation=min_simplex_3d;aff_log("SIMPLEX ");break;
    }
  }
  else
  {
    switch (min_type)
    {
      case 0: minimisation=min_icm_3d;aff_log("ICM ");break;
	  case 1: minimisation=min_simplex_3d;aff_log("SIMPLEX ");break;
	  case 2: minimisation=min_powel_3d;aff_log("POWEL ");break;
    case 3: minimisation=min_powel2_3d;aff_log("POWEL2 ");break;
	  default:minimisation=min_simplex_3d;aff_log("SIMPLEX ");break;
    }
  }
  return minimisation;
}

/*******************************************************************************
**     arret_selon_parametres(new_param, old_param, prec_param, nb_param)
*/
/*!    determine si on a atteint la precision sur l'evolution des parametres
**
**
**     utilisation:
**     \param       new_param, old_param: parametres avant et apres l'iteration de minimisation
**     \param       prec_param : precisions voulues sur les parametres
**     \param       nb_param : nb de parametres
**     \retval      determine si on peut arreter la minimisation
*******************************************************************************/
bool arret_selon_parametres(double *new_param, double *old_param, double *prec_param, int nb_param)
{
 int i;
 bool arret=TRUE;

 for (i=0;i<nb_param;i++)
 {
  if (prec_param[i]!=0.0) arret = arret&&(fabs(new_param[i]-old_param[i])<=prec_param[i]);
 }

 return arret;
}

/*******************************************************************************
**     arret_selon_energie(E1, E2)
*/
/*!    determine si on a atteint la precision sur l'evolution de l'energie
**
**
**     utilisation:
**     \param       E1, E2: energies avant et apres l'iteration de minimisation
**     \retval      determine si on peut arreter la minimisation
*******************************************************************************/
bool arret_selon_energie(double E1, double E2)
{
 bool arret;
 double rtol;

 rtol=2.0*fabs(E1-E2)/(fabs(E1)+fabs(E2)+EEPS);
 arret=(rtol<EEPS);

 return arret;
}

