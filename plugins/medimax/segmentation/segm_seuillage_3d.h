/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*!************************************************************************
***	
***	\file:		segm_seuillage_3d.h
***
***	project:	Imagix 2.01
***			
***
***	\brief description:    Fichier source pour segmentation 3D
***	
***	
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	Copyright (c) 1997, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***
*************************************************************************/

#ifndef __SEGM_SEUILLAGE_3D_H
#define __SEGM_SEUILLAGE_3D_H


extern void    threh_cor_3d(void) ;
extern int    imx_threh_cor_3d(int im_1, int im_res, float bas, float haut) ;
extern void    threhcluster_cor_3d(void) ;
extern int    imx_threhcluster_cor_3d(int im_1, int im_res, float bas, float haut) ;
extern int    imx_threhcluster_cor_3d_p(grphic3d *im1, grphic3d *imres, float bas, float haut) ;

#endif /*__SEGM_SEUILLAGE_3D_H*/

