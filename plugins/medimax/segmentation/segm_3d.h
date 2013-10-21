/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "noyau/imx_lang.h"

/* AUTOMATICALLY EXTRACTED PROTOTYPES (MB 6/12/2000) */
#ifdef __cplusplus 
extern "C" 
{
#endif /*__cplusplus*/ 
extern void threh_cluster_3d(void) ;
extern void opt_threh_3d(void) ;
extern int  imx_opt_threh_3d(int im_1, int im_res) ;
extern int  seg_cerv_grow_dis_3d(void) ;
extern int  imx_seg_cerv_grow_dis_3d(int im_deb, int im_res, int x0, int y0, int z0) ;
extern void imx_moyenne_3d_p(grphic3d *imdeb, grphic3d *imres) ;
extern int  imx_seg_cerv_grow_dis_3d_p(grphic3d *imdeb, grphic3d *imres, int x0, int y0, int z0) ;
extern int  grow_spher_cerv_3d(int option) ;
extern int  imx_grow_spher_cerv_3d(int im_deb, int im_res, int x0, int y0, int z0, int option) ;
extern int  cerveau_to_ROI_3d(int option) ;
extern int  imx_cerveau_to_ROI_3d(int im_deb, int im_res, int x0, int y0, int z0, int option) ;
extern int  imx_cerveau_to_ROI_3d_p(grphic3d *imdeb, grphic3d *imres, int x0, int y0, int z0, int option) ;
extern void otsu_head_3d(void) ;
extern void imx_otsu_head_3d(int im_deb, int im_res);
extern void imx_otsu_head_3d_p(grphic3d *imdeb, grphic3d *imres);
extern DllExport int imx_threh_img_3d_p (grphic3d *im1, grphic3d *imres, float bas, float haut);
extern int	thre_cerv_3d(void);
extern int	imx_threh_img_3d(int im_1, int im_res, float bas, float haut);
extern void thresh_img_with_mask_3d();

extern int imx_threh_img_outliers_3d_p (grphic3d *im1, grphic3d *imres, float bas, float haut);

#ifdef __cplusplus 
}
#endif /*__cplusplus*/ 
