/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

//
//  File    : imx_picture3d.h
//  Project : GTK-Imagix
//
//  Description:
//      ...
//
//  Copyright (c) 1993, ULP-IPB Strasbourg
//  All rights are reserved
//
//  >>  Created on ... ..th, .... by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//

#ifndef __IMX_PICTURE3D_H
#define __IMX_PICTURE3D_H

#ifdef WIN32
#define DllExport   __declspec( dllexport )  // for windows
#else
#define DllExport  // Linux
#endif

#include "noyau/imx_3d.h"
#include "noyau/imx_2d.h"

extern int DllExport show_picture_3d (int pos);
extern int picture_3dto2d(grphic3d *img) ;
extern int picture_3dto2d_param(grphic3d *img, grphic *imt, grphic *ims, grphic *imc, grphic *imcurv);
extern int picture_3dto2d_copy(grphic3d *img, grphic *imt, grphic *ims, grphic *imc, grphic *imcurv);
extern int picture_mask_3dto2d(grphic3d *img);
extern void Refresh_3d(void);
extern void Zoom_Vi_3d(void) ;
extern void Zoom_Vin_3d(void) ;
extern void Zoom_Vim_3d(void) ;
extern void Zoom_Vi_3d_pos(void) ;
extern void Zoom_Vi_3d_pos_xyz(int x, int y, int z, int pos ,int nzoom);
extern int imx_reinitialise_visual_params_3d_p(grphic3d *imres);




#endif  // #ifndef __IMX_PICTURE2D_H
