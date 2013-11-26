/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

//
//  File    : zone_3d.h
//  Project : MEDIMAX
//
//  Description:
//      ...
//
//  Copyright (c) 1993, ULP-IPB Strasbourg
//  All rights are reserved
//
//  >>  Created on ... ..th, 2001 by ...
//  >>  Last user action on ... ..th, 2001 by ...
//

#ifndef __ZONE_3D_H
#define __ZONE_3D_H

#ifdef __cplusplus 
extern "C" 
{
#endif /*__cplusplus*/ 

// Definition de la structure de donnees pour les fichier ROI 3D

#define COMPRESSED 0x80000000L;
#define UNCOMPRESSED 0x0L;

#define BUFFER_SIZE   16384

typedef struct {
    char sign[8];
    int width, height, depth;
    int bitperpixel;
    int icomp;
    int min_pixel, max_pixel;
    int cutoff_min, cutoff_max;
    int x3dc, y3dc, z3dc;
    int zoom, zoom_x, zoom_y, zoom_z;
    int pos_x, pos_y, pos_z;
    float dx, dy, dz;
    struct s_patient patient;
} ROI_3D;

extern int  show_roi_3d(int pos3d) ;
extern void Affichage_roi_3d(int pos3d) ;
extern void free_display_mask_3d_p(grphic3d *im);
extern int  draw_roibythreshold_3d(grphic3d *roi) ;
extern void combine_data(short int *mri, long int *buffer, int num_items, int shift, int mode) ;
extern void Aff_Roi_3d(void) ;
extern void Aff_nRoi_3d(void) ;
extern void Unaff_Roi_3d(void) ;
extern void Unaff_nRoi_3d(void) ;
extern void Eff_Mask_3d(void) ;
extern void Eff_nMask_3d(void) ;
extern void Sup_Mask_3d(void) ;
extern void Sup_nMask_3d(void) ;
extern void ptr_mask_activate_3d(int pos3d, BOOL bActive);
extern void ptr_mask_activate_3d_ptr(grphic3d *img, BOOL bActive);
extern void op_log_mask_3d();
extern void position_ds_talairach (grphic3d *im, dpoint3d point, char *reponse);
extern void Iz_DrawClosed_3d(void);
extern void Iz_DrawOpened_3d(void);
extern void fill_mask_3d(void);
extern void imx_getMaskStatus_3d_p(grphic3d* img,char *retString);
extern int imx_change_mask_type_3d_p(grphic3d *img);
extern void Iz_threh_3d();
extern int  imx_Iz_threh_3d_p(grphic3d *im1, grphic3d *imres, float bas, float haut);
extern void Iz_threh_man_3d(void);


#ifdef __cplusplus 
}
#endif /*__cplusplus*/ 

#endif
