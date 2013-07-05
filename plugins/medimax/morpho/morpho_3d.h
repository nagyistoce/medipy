/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

 
/* AUTOMATICALLY EXTRACTED PROTOTYPES (MB 6/12/2000) */
#ifdef __cplusplus 
extern "C" 
{
#endif /*__cplusplus*/ 
extern void     eros_3d(void) ;
extern int     imx_eros_3d(int im_1, int im_res, int numas, int niter) ;
extern int     imx_eros_3d_p(grphic3d *im1, grphic3d *imres, int numas, int niter) ;
extern void    dilat_3d(void) ;
extern int     imx_dilat_3d(int im_1, int im_res, int numas, int niter) ;
extern int     imx_dilat_3d_p(grphic3d *im1, grphic3d *imres, int numas, int niter) ;
extern void    GradM_3d(void) ;
extern int     imx_GradM_3d(int im_1, int im_res, int numas, int numgrad, int niter) ;
extern int     imx_GradM_3d_p(grphic3d *im1, grphic3d *imres, int numas, int numgrad, int niter) ;
extern void    open_3d(void) ;
extern int     imx_open_3d(int im_1, int im_res, int numas, int niter) ;
extern int     imx_open_3d_p(grphic3d *im1, grphic3d *imres, int numas, int niter) ;
extern void    close_3d(void) ;
extern int     imx_close_3d(int im_1, int im_res, int numas, int niter) ;
extern int     imx_close_3d_p(grphic3d *im1, grphic3d *imres, int numas, int niter) ;
extern void    tophat_3d(void) ;
extern int     imx_tophat_3d(int im_1, int im_res, int numas, int black_white_tophat, int niter) ;
extern int     imx_tophat_3d_p(grphic3d *im1, grphic3d *imres, int numas, int black_white_tophat, int niter) ;
extern void reconstruction_geod_dilat_3d(void) ;
extern int     imx_reconstruction_geod_dilat_3d(int im_1, int im_2, int im_res, int answer_nr) ;
extern int     imx_reconstruction_geod_dilat_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres, int numas) ;
extern void     imx_reconstruction_geod_dilat_bin_3d_p(grphic3d *imtemp2, grphic3d *imres, int x3d, int y3d, int z3d, int im_res, int numas) ;
extern void reconstruction_geod_eros_3d(void) ;
extern int     imx_reconstruction_geod_eros_3d(int im_1, int im_2, int im_res, int answer_nr) ;
extern int     imx_reconstruction_geod_eros_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres, int numas) ;
extern void    sharpen_morpho_3d(void) ;
extern int     imx_sharpen_morpho_3d(int im_1, int im_res, int numas, int niter) ;
extern int     imx_sharpen_morpho_3d_p(grphic3d *im1, grphic3d *imres, int numas, int niter) ;
extern void    contrast_morpho_3d(void) ;
extern int     imx_contrast_morpho_3d(int im_1, int im_res, int numas, int niter) ;
extern int     imx_contrast_morpho_3d_p(grphic3d *im1, grphic3d *imres, int numas, int niter) ;
extern int     eros_bin_3d(void) ;
extern void    FAS_3d(void) ;
extern void    FASS_3d(void) ;
extern XPoints_3d  *masque_3d(int numas) ;
extern int     imx_eros_bin_3d(int im_1, int im_res, int numas, int niter) ;
extern int     imx_eros_bin_3d_p(grphic3d *im1, grphic3d *imres, int numas) ;


#ifdef __cplusplus 
}
#endif /*__cplusplus*/ 
