/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*!
//  File    : imx_tmp_grphic3d.h
//  Project : GTK-Imagix
//
//  Description: Gestion des images 3D temporaires
//      ...
//
//  Copyright (c) 1993, ULP-IPB Strasbourg
//  All rights are reserved
//
//
*/

#ifndef _IMX_TMP_GRPHIC3D_H
#define __IMX_TMP_GRPHIC3D_H

/*function declaration*/
extern grphic3d *ptr_img_tmp_3d(int nr);
extern grphic3d *cr_img_tmp_3d(t_ptrlist liste,int nr);
extern int free_img_tmp_3d(t_ptrlist list,int nr);

/* global variable*/
extern t_ptrlist _img_tmp_3d;


#endif /*__IMX_TMP_GRPHIC3D_H*/
