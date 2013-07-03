/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*******************************************************************************/
/*!	   \file:		erreur_recalage.c
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
#include "recalage/erreur_recalage.h"


bool erreur_critique(ERREUR_RECALAGE *err)
{
 if (!err) return FALSE;

 if (((*err)!=NO_ERR)&&((*err)!=NB_BINS_NUL)&&((*err)!=OVERLAP_NUL)) return TRUE;

 return FALSE;
}
