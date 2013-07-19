/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*******************************************************************************/
/*!	   \file:		transformations_3d.h
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

#ifndef __TRANSFORMATIONS_3D_H__
#define __TRANSFORMATIONS_3D_H__

//---zoom isotrope---//
extern void zoom_isotrope_3d();
extern int  imx_zoom_isotrope_3d(int im_src, int im_dest, int res_type, int inter_type, int save_type, char * nomfichres);
extern int  imx_zoom_isotrope_3d_p(grphic3d *imsrc, grphic3d *imdst, int res_type, int inter_type,int save_type, char* nomfichres, int size_max);

//--- Zoom relatif ---//
extern void zoom_relatif_3d(void);
extern int imx_zoom_relatif_3d(int im_zoom, int im_ref, int im_res, int inter_type, int save_type, char *nomfichres);
extern int imx_zoom_relatif_3d_p(grphic3d *imzoom, grphic3d *imref, grphic3d *imres, int inter_type, int save_type, char *nomfichres);

//--- Rotation de l'image ---//
extern void rot_3d();
extern int imx_rot_3d(int im_ref, int im_res, int inter_type, double ang1, double ang2, double ang3);
extern int imx_rot_3d_p(grphic3d *imref, grphic3d *imres, int inter_type, double ang1, double ang2, double ang3);

//---calcul du champ corrrespondant a une transformation---//
extern int transf_geom_to_field3d(transf3d *transfo, field3d *champ, grphic3d *imref, grphic3d *imtransf, field3d* points_calc);



#endif /*__TRANSFORMATIONS_3D_H__*/
