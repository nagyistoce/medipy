/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

//
/*!  \file    : imx_list.h
//  Project : Imagix
//
//  Description:
//      Lists in C
//
//  Copyright (c) 1993, ULP-IPB Strasbourg
//  All rights are reserved
//
//  >>  Created on ... ..th, .... by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//
*/

#ifndef __IMX_LIST_H
#define __IMX_LIST_H

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief type element de liste            */
typedef void *ptrlistitem;				

/*! \brief type maillon de chaine de liste  */
typedef struct s_listnode				
{
	struct s_listnode *prev;  /*!< pointeur sur l'element precedent */
	struct s_listnode *next;  /*!< pointeur sur l'element suivant */
	ptrlistitem item;		  /*!< pointeur sur les donnees*/
} t_listnode, *t_ptrlistnode;

/*! \brief type liste                       */
typedef struct s_list					
{
	t_ptrlistnode first;    /*!< pointeur sur le premier element */
	t_ptrlistnode curs;		/*!< pointeur sur l'element courant */
	int           nbitems;  /*!< nbre d'element de la liste*/
} t_list , *t_ptrlist;

#ifndef __DOXYGEN_SHOULD_SKIP_THIS__
int list_init(t_ptrlist *plist);                           /* initialisation d'une liste      */
int list_erase(t_ptrlist *plist);                          /* suppression de la liste         */
int list_clean(t_ptrlist *plist);                          /* nettoyage de la liste           */

int list_movefirst(t_ptrlist *plist);                      /* curseur sur premier element     */
int list_movelast(t_ptrlist *plist);                       /* curseur sur dernier element     */
int list_movenext(t_ptrlist *plist);                       /* curseur sur element suivant     */
int list_moveprev(t_ptrlist *plist);                       /* curseur sur element precedent   */
int list_movetop(t_ptrlist *plist);                        /* curseur sur haut de pile        */
int list_movedeep(t_ptrlist *plist);                       /* element empile suivant          */
int list_moveto(t_ptrlist *plist, int n);                         /* curseur sur l'element n         */

int list_delete(t_ptrlist *plist);                         /* supprime l'element courant      */
int list_insert(t_ptrlist *plist, ptrlistitem item);                         /* insere avant l'element courant  */
int list_push(t_ptrlist *plist, ptrlistitem item);                           /* ajoute en tete de liste         */
ptrlistitem list_pop(t_ptrlist *plist);                    /* supprime le premier element     */

ptrlistitem list_item(t_ptrlist *plist);                   /* element courant                 */
ptrlistitem list_top(t_ptrlist *plist);                    /* premier element                 */

int list_isempty(t_ptrlist *plist);                        /* la liste est vide ?             */
int list_nbitems(t_ptrlist *plist);                        /* nombre d'elements de la liste ? */
int list_pos(t_ptrlist *plist);                            /* position du curseur ?           */
int list_isatend(t_ptrlist *plist);                        /* curseur en fin de liste ?       */
int list_isatbeg(t_ptrlist *plist);                        /* curseur en debut de liste ?     */

char *list_tostring(t_ptrlist *plist, int size);                     /* convertit la liste en chaine    */
int list_sortbystring(t_ptrlist *plist);                   /* tri sur chaine                  */
#endif /*__DOXYGEN_SHOULD_SKIP_THIS__*/

#ifdef __cplusplus
}
#endif

#endif  // #ifndef __IMX_LIST_H
