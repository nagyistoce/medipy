/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*******************************************************************************/
/*!	   \file:		erreur_recalage.h
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

#ifndef __ERREUR_RECALAGE_H__
#define __ERREUR_RECALAGE_H__

#include "recalage/definitions_types_recalage.h"

#define RAISE_ERR_TEXT(err_ptr, err_type, texte) { aff_log(texte); if (err_ptr) *err_ptr=err_type; goto end_func; }

extern bool erreur_critique(ERREUR_RECALAGE *err);

#endif /*__ERREUR_RECALAGE_H__*/
