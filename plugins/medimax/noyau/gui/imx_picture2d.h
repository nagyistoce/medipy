/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

//
//  File    : imx_picture2d.h
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

#ifndef __IMX_PICTURE2D_H
#define __IMX_PICTURE2D_H

#include "noyau/imx_types.h"
#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum
{
    GENERIC = 0,
    ZOOM_1,
    ZOOM_2,
    ZOOM_4
} TPict2dDisplayFuncId;

typedef void (*TPict2dDisplayFunc)(ptr_grphic, TMaskType, UCHAR *);

#ifdef __GTK_GUI
typedef struct
{
    TMaskType               masktype;
    TPict2dDisplayFuncId    id;
    t_CMapType              type;
    t_CMapDisplayMode       mode;
    TPict2dDisplayFunc      func;
} TPict2dDisplayFuncItem;
#endif /* __GTK_GUI */

// extern UCHAR *sp_build_displayable_image_data(grphic *, int, int, int, int, int, int, int, int, int, UCHAR *);

extern void IMXPict2d_BuildDisplayableImageData(ptr_grphic pImage, UCHAR *pBuff, TMaskType mask_type);
extern void IMXPict2d_Initialize();

#ifdef __GTK_GUI
extern void (*IMXPict2d_ShowBar)(grphic *pImage, UCHAR *pDest);
extern void IMXPict2d_ShowIndexes(int noCMap, ptrColormap cmap, TDisplayIndexPos nIndexPos);
extern int	unzoom_img(int pos);
//extern void Refresh_Order(void)
#endif /* __GTK_GUI */

#ifdef __cplusplus
}
#endif

#endif  // #ifndef __IMX_PICTURE2D_H
