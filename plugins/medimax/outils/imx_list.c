/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*!  \file    : imx_list.c 
**   \brief     Lists in C
**
**  Copyright (c) 1993, ULP-IPB Strasbourg
**  All rights are reserved
**
**  >>  Created on ... ..th, .... by F. BOULEAU
**  >>  Last user action on ... ..th, 2001 by ...
******************************************************************************
**                         LISTES CHAINEES BILATERES 
*******************************************************************************
**
** >> Derniere revision par F. BOULEAU le 17 Aout 1999
**
******************************************************************************
**
*/

#include <config.h>
#include <stdio.h>
#include <string.h>

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "outils/imx_list.h"


/******************************************************************************
** -- list_init() -------------------------------------------------------------
**
** >> Derniere revision par F. BOULEAU le 17 Aout 1999
*/
/*!  \brief Initialisation d'une nouvelle liste
**	 \param plist  : la liste a initialiser
**	 \retval 0
******************************************************************************/

int list_init(t_ptrlist *plist)
{
	*plist            = CALLOC(1, t_list);
	(*plist)->nbitems = 0;
	(*plist)->first   = (*plist)->curs = (t_ptrlistnode) NULL;

	return 0;
}

/******************************************************************************
** -- list_erase() ------------------------------------------------------------
**
** >> Derniere revision par F. BOULEAU le 18 Aout 1999
*/
/*! \brief Suppression d'une liste
**	\param plist : la liste a supprimer
**	\retval 0
******************************************************************************/

int list_erase(t_ptrlist *plist)
{
	t_ptrlist list = *plist;

	if(list == NULL)
		return 0;

	while(list->first != (t_ptrlistnode) NULL)
	{
		list->curs  = list->first;
		list->first = list->first->next;

		if (list->curs->item != NULL)
			FREE(list->curs->item);

		FREE(list->curs);
	}

	FREE(list);
	*plist = NULL;
	
	return 0;
}

/******************************************************************************
** -- list_clean() ------------------------------------------------------------
**
** >> Derniere revision par F. BOULEAU le 17 Aout 1999
*/
/*!
**  \brief Suppression des elements d'une liste
**	\param plist : la liste a vider
**	\retval  0
******************************************************************************/

int list_clean(t_ptrlist *plist)
{
	t_ptrlist list = *plist;

    if(!list)
        return 0;

	while(list->first != (t_ptrlistnode) NULL)
	{
		list->curs  = list->first;
		list->first = list->first->next;

		if (list->curs->item != NULL)
			FREE(list->curs->item);

		FREE(list->curs);
	}

	list->curs    = list->first;
	list->nbitems = 0;

	return 0;
}

/******************************************************************************
** -- list_movefirst() --------------------------------------------------------
**
** >> Derniere revision par F. BOULEAU le 17 Aout 1999
*/
/*!  \brief Deplace le curseur sur le premier element de la liste
**	 \param plist : la liste
**	 \retval 0
**	 \remark plist est modifiee
**
******************************************************************************/

int list_movefirst(t_ptrlist *plist)
{
	(*plist)->curs = (*plist)->first;
	return 0;
}

/******************************************************************************
** -- list_movelast() ---------------------------------------------------------
**
** >> Derniere revision par F. BOULEAU le 17 Aout 1999
*/
/*!	\brief Deplace le curseur sur le dernier element de la liste
**	\param plist : la liste
**	\retval 0
**
******************************************************************************/

int list_movelast(t_ptrlist *plist)
{
	(*plist)->curs = (t_ptrlistnode) NULL;
	return 0;
}

/******************************************************************************
** -- list_movenext() ---------------------------------------------------------
**
** >> Derniere revision par F. BOULEAU le 17 Aout 1999
*/
/*!	\brief Deplace le curseur sur l'element suivant
**	\param plist : la liste
**	\retval 0
**
******************************************************************************/

int list_movenext(t_ptrlist *plist)
{
	if ((*plist)->curs != (t_ptrlistnode) NULL)
		(*plist)->curs = (*plist)->curs->next;
	return 0;
}

/******************************************************************************
** -- list_moveprev() ---------------------------------------------------------
**
** >> Derniere revision par F. BOULEAU le 17 Aout 1999
*/
/*!  \brief Deplace le curseur sur l'element precedent
**	\param plist : la liste
**	\retval 0
**
******************************************************************************/

int list_moveprev(t_ptrlist *plist)
{
	t_ptrlist list = *plist;
	
	if (list->curs != (t_ptrlistnode) NULL)
	{
		if (list->curs->prev != (t_ptrlistnode) NULL)
			list->curs = list->curs->prev;
	}
	else
	{
		list->curs = list->first;

		if (list->curs != (t_ptrlistnode) NULL)
			while(list->curs->next != (t_ptrlistnode) NULL)
				list->curs = list->curs->next;
	}
		
	return 0;
}

/******************************************************************************
** -- list_movetop() ----------------------------------------------------------
**
** >> Derniere revision par F. BOULEAU le  5 Oct. 1999
*/
/*! \brief  Deplace le curseur sur l'element en haut de pile
**	\param plist : la liste
**	\retval  1 si element courant de list est NULL, 0 sinon
**
******************************************************************************/

int list_movetop(t_ptrlist *plist)
{
	t_ptrlist list = *plist;
	list->curs = list->first;
	return (list->curs == (t_ptrlistnode) NULL);
}

/******************************************************************************
** -- list_movedeep() ---------------------------------------------------------
**
** >> Derniere revision par F. BOULEAU le  5 Oct. 1999
*/
/*! \brief Deplace le curseur sur l'element suivant de la pile
**	\param plist : la liste
**	\retval  1 si element courant de list est NULL, 0 sinon
**
******************************************************************************/

int list_movedeep(t_ptrlist *plist)
{
	t_ptrlist list = *plist;
	if (list->curs != (t_ptrlistnode) NULL)
		list->curs = list->curs->next;
	return (list->curs == (t_ptrlistnode) NULL);
}

/******************************************************************************
** -- list_moveto() -----------------------------------------------------------
**
** >> Derniere revision par F. BOULEAU le 17 Aout 1999
*/
/*! \brief  Deplace le curseur a une position determinee
**	\param plist : la liste
**	\param n : numero du noeud
**	\retval n
**
******************************************************************************/

int list_moveto(t_ptrlist *plist, int n)
{
	t_ptrlist list = *plist;

	list->curs = list->first;

	for( ; n > 1 && list->curs != (t_ptrlistnode) NULL ; n--)
		list->curs = list->curs->next;

	return n;
}

/******************************************************************************
** -- list_delete() -----------------------------------------------------------
**
** >> Derniere revision par F. BOULEAU le 17 Aout 1999
*/
/*! \brief Supprime l'element courant
**	\param plist : la liste
**	\retval 0
**
******************************************************************************/

int list_delete(t_ptrlist *plist)
{
	t_ptrlist     list = *plist;
	t_ptrlistnode tmp;

	if (list->curs != (t_ptrlistnode) NULL)
	{
		tmp = list->curs;
		list->curs = list->curs->next;

		if (tmp->next != (t_ptrlistnode) NULL)
			tmp->next->prev = tmp->prev;

		if (tmp->prev == (t_ptrlistnode) NULL)
			list->first = list->curs;
		else
			tmp->prev->next = tmp->next;

		list->nbitems--;

		if (tmp->item != NULL)
			FREE(tmp->item);

		FREE(tmp);
	}

	return 0;
}

/******************************************************************************
** -- list_insert() -----------------------------------------------------------
**
** >> Derniere revision par F. BOULEAU le 17 Aout 1999
*/
/*! \brief Insere un element a la position du curseur
**	\param plist : la liste
**	\param item : l'item a inserer
**	\retval 0
**
******************************************************************************/

int list_insert(t_ptrlist *plist, ptrlistitem item)
{
	t_ptrlist     list = *plist;
	t_ptrlistnode tmp;
	t_ptrlistnode curs;

	/* creee le nouvel element */

	tmp       = CALLOC(1, t_listnode);
	tmp->item = item;

	curs = list->curs;

	/* insertion en queue de liste */
	if (curs == (t_ptrlistnode) NULL)
	{
		curs = list->first;

		/* cherche le dernier element de la liste */
		while(curs          != (t_ptrlistnode) NULL 
			  && curs->next != (t_ptrlistnode) NULL)
			curs = curs->next;

		/* la liste est vide */
		if (curs == (t_ptrlistnode) NULL)
		{
			list->first = tmp;
			tmp->next   = tmp->prev = (t_ptrlistnode) NULL;
		}
		/* insere nouvel element */
		else
		{
			curs->next = tmp;
			tmp->prev  = curs;
			tmp->next  = (t_ptrlistnode) NULL;
		}
	}
	/* insertion entre deux elements */
	else
	{
		/* chainage des elements */

		tmp->next = curs;
		tmp->prev = curs->prev;

		curs->prev = tmp;

		if (tmp->prev != (t_ptrlistnode) NULL)
			tmp->prev->next = tmp;

		if (list->first == curs)
			list->first = tmp;

		list->curs = tmp;
	}

	list->nbitems++;

	return 0;
}

/******************************************************************************
** -- list_push() -------------------------------------------------------------
**
** >> Derniere revision par F. BOULEAU le 17 Aout 1999
**
*/
/*! \brief Insere un element en tete de liste
**	\param plist : la liste
**	\param item : l'item a inserer
**	\retval 0
**
******************************************************************************/

int list_push(t_ptrlist *plist, ptrlistitem item)
{
	t_ptrlist     list = *plist;
	t_ptrlistnode tmp  = CALLOC(1, t_listnode);

	tmp->prev = (t_ptrlistnode) NULL;
	tmp->next = list->first;
	tmp->item = item;

	if (list->first != (t_ptrlistnode) NULL)
		list->first->prev = tmp;

	list->first = tmp;
	list->nbitems++;

	return 0;
}

/******************************************************************************
** -- list_pop() --------------------------------------------------------------
**
** >> Derniere revision par F. BOULEAU le 17 Aout 1999
*/
/*! \brief Supprime l'element en tete de liste
**	\param plist : la liste
**	\retval l'item supprime
**
******************************************************************************/

ptrlistitem list_pop(t_ptrlist *plist)
{
	t_ptrlist     list = *plist;
	t_ptrlistnode tmp  = list->first;
    ptrlistitem   item = NULL;

	if (tmp != (t_ptrlistnode) NULL)
	{
		list->first = tmp->next;

		if (tmp->next != (t_ptrlistnode) NULL) 
		{
			tmp->next->prev = (t_ptrlistnode) NULL;
			if(list->curs == list->first)
				list->curs = list->curs->next;
		}
		else
			list->curs = NULL;

        item = tmp->item;

		FREE(tmp);

		list->nbitems--;
	}

	return item;
}

/******************************************************************************
** -- list_item() -------------------------------------------------------------
**
** >> Derniere revision par F. BOULEAU le 17 Aout 1999
*/
/*! \brief Renvoie un pointeur sur l'element courant
**	\param plist : la liste
**	\retval l'element courant
**  \remark la liste doit etre valide
******************************************************************************/

ptrlistitem list_item(t_ptrlist *plist)
{
	return ((*plist)->curs ? (*plist)->curs->item : (ptrlistitem) NULL);
}

/******************************************************************************
** -- list_top() --------------------------------------------------------------
**
** >> Derniere revision par F. BOULEAU le 17 Aout 1999
*/
/*! \brief Renvoie un pointeur sur l'element en tete de liste
**	\param plist : la liste
**	\retval l'element en tete de liste
**
******************************************************************************/

ptrlistitem list_top(t_ptrlist *plist)
{
	return ((*plist)->first ? (*plist)->first->item : (ptrlistitem) NULL);
}

/******************************************************************************
** -- list_isempty() ----------------------------------------------------------
**
** >> Derniere revision par F. BOULEAU le 17 Aout 1999
*/
/*! \brief Indique si la liste est vide
**	\param plist : la liste
**	\retval 1 si la liste est vide, 0 sinon
**
******************************************************************************/

int list_isempty(t_ptrlist *plist)
{
	return ((*plist)->nbitems == 0 ? 1 : 0);
}

/******************************************************************************
** -- list_nbitems() ----------------------------------------------------------
**
** >> Derniere revision par F. BOULEAU le 17 Aout 1999
*/
/*! \brief Nombre d'elements de la liste
**	\param plist : la liste
**	\retval nbre d'element de la liste
**	\remark la liste passee en parametre doit etre valide
******************************************************************************/

int list_nbitems(t_ptrlist *plist)
{
	return (*plist)->nbitems;
}

/******************************************************************************
** -- list_pos() --------------------------------------------------------------
**
** >> Derniere revision par F. BOULEAU le 17 Aout 1999
*/
/*! \brief Position du curseur
**	\param plist : la liste
**	\retval position du curseur dans la liste
**
******************************************************************************/

int list_pos(t_ptrlist *plist)
{
	t_ptrlist     list = *plist;
	t_ptrlistnode tmp  = list->first;
	int           i=0;

	if (list->curs == NULL)
		i = list->nbitems + 1;
	else
		for(i = 1 ; tmp != list->curs ; i++)
			tmp = tmp->next;
	
	return i;
}

/******************************************************************************
** -- list_tostring() ---------------------------------------------------------
**
** >> Derniere revision par F. BOULEAU le 17 Aout 1999
*/
/*! \brief Renvoie la liste sous forme de chaine de caracteres
**	\param plist : la liste
**	\param size : la taille de la chaine qui sera alloue
**	\retval char *  : une chaine de caracteres contenant les elements caste en (char*)
**
******************************************************************************/

char *list_tostring(t_ptrlist *plist, int size)
{
	t_ptrlist      list = *plist;
	char          *s    = CALLOC(size, char);
	t_ptrlistnode  tmp  = list->first;

	strcpy(s, "");

	while(tmp != NULL)
	{
		if (tmp->item != (ptrlistitem) NULL)
		{
			if (tmp != list->first)
				strcat(s, ", ");
			strcat(s, (char *) tmp->item);
		}

		tmp = tmp->next;
	}

	return s;
}

/******************************************************************************
** -- list_sortbystring() -----------------------------------------------------
**
** >> Derniere revision par F. BOULEAU le 17 Aout 1999
**
*/
/*! \brief Trie une liste de chaines de caracteres
**  On recherche le plus petit element parmi les n - i derniers elements de la
**  la liste, puis on le place en position i et on recommence avec i + 1 
**  (avec i = 0 au depart).
**	\param plist : la liste
**	\retval 0
**
******************************************************************************/

int list_sortbystring(t_ptrlist *plist)
{
	t_ptrlist     list     = *plist;
	t_ptrlistnode tmp_depl = list->first;
	t_ptrlistnode tmp_curs;
	t_ptrlistnode tmp_new;
	ptrlistitem   tmp_old;

	while(tmp_depl != (t_ptrlistnode) NULL)
	{
		tmp_curs = tmp_depl->next;
		tmp_new  = tmp_depl;

		/* cherche le plus petit element de la liste commencant par tmp_depl */

		while(tmp_depl->item != (ptrlistitem)   NULL 
			  && tmp_curs    != (t_ptrlistnode) NULL)
		{
			if (tmp_curs->item == NULL
			    || strcmp((char *) tmp_curs->item, 
				          (char *) tmp_new->item) < 0)
				tmp_new = tmp_curs;

			tmp_curs = tmp_curs->next;
		}

		/* le curseur reste a la meme position */

		if (list->curs == tmp_new)
			list->curs = tmp_depl;
		else if (list->curs == tmp_depl)
			list->curs = tmp_new;

		/* permutation des elements : le plus petit element de la sous-liste */
		/* vient en tete de cette sous-liste                                 */

		if (tmp_new != tmp_depl)
		{
			tmp_old        = tmp_depl->item;
			tmp_depl->item = tmp_new->item;
			tmp_new->item  = tmp_old;
		}

		tmp_depl = tmp_depl->next;
	}

	return 0;
}

/******************************************************************************
** -- list_isatend() ----------------------------------------------------------
**
** >> Derniere revision par F. BOULEAU le 19 Aout 1999
*/
/*! \brief Indique si le curseur est en fin de liste
**	\param plist : la liste
**	\retval 1 si la fin : 0 sinon
**
******************************************************************************/

int list_isatend(t_ptrlist *plist)
{
	return ((*plist)->curs == NULL);
}

/******************************************************************************
** -- list_isatbeg() ----------------------------------------------------------
**
** >> Derniere revision par F. BOULEAU le 19 Aout 1999
*/
/*! \brief  Indique si le curseur est en debut de liste
**	\param plist : la liste
**	\retval 1 si le debut : 0 sinon
**
******************************************************************************/

int list_isatbeg(t_ptrlist *plist)
{
	return ((*plist)->first == (*plist)->curs);
}
