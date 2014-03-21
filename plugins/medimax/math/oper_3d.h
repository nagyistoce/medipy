/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _oper_3d_h
#define _oper_3d_h

#ifdef __cplusplus
extern "C"{
#endif

extern int imx_add_3d (int im_1, int im_2, int im_res);
extern int imx_sub_3d (int im_1, int im_2, int im_res);
extern int imx_mul_3d (int im_1, int im_2, int im_res);
extern int imx_div_3d (int im_1, int im_2, int im_res);
extern int imx_abs_3d (int im_1, int im_res);
extern int imx_sqr_3d (int im_1, int im_res);
extern int imx_log_3d (int im_1, int im_res);
extern int imx_exponent_3d (int im_1, int im_res);
extern int imx_max_3d (int im_1, int im_2, int im_res);
extern int imx_min_3d (int im_1, int im_2, int im_res);
extern void imx_add_coe_3d (int im_1, float coefficient, int im_res);
extern void imx_sub_coe_3d (int im_1, float coefficient, int im_res);
extern void imx_mul_coe_3d (int im_1, float coefficient, int im_res);
extern void imx_div_coe_3d (int im_1, float coefficient, int im_res);
extern int imx_and_3d (int im_1, int im_2, int im_res);
extern int imx_or_3d (int im_1, int im_2, int im_res);
extern int imx_eg_3d (int im_1, int im_2, int im_res);
extern int imx_inv_3d (int im_1, int im_res);
extern int imx_and2img_3d (int im_1, int im_2);
extern int imx_add_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern int imx_max_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern int imx_min_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern int imx_sub_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern int imx_subabs_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern int imx_mul_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern int imx_div_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern int imx_sqr_3d_p (grphic3d *im1, grphic3d *imres);
extern int imx_log_3d_p (grphic3d *im1, grphic3d *imres);
extern int imx_exponent_3d_p (grphic3d *im1, grphic3d *imres);
extern int imx_abs_3d_p (grphic3d *im1, grphic3d *imres);
extern int imx_mul_coe_3d_p (grphic3d *im1, float coefficient, grphic3d *imres);
extern int imx_div_coe_3d_p (grphic3d *im1, float coefficient, grphic3d *imres);
extern int imx_add_coe_3d_p (grphic3d *im1, float coefficient, grphic3d *imres);
extern int imx_sub_coe_3d_p (grphic3d *im1, float coefficient, grphic3d *imres);
extern int imx_and_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern int imx_or_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern int imx_inv_3d_p (grphic3d *im1, grphic3d *imres);
extern int imx_and2img_3d_p (grphic3d *im1, grphic3d *im2);
extern int imx_eg_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern int imx_subabs_3d (int im_1, int im_2, int im_res);

extern int IMXColor_BindColormap_3d_ptr(ptr_grphic3d ptrPic, int nNewCMap);

#ifdef __cplusplus
}
#endif

#endif /* _oper_3d_h */
 
/* AUTOMATICALLY EXTRACTED PROTOTYPES (MB 6/12/2000) */
#ifdef __cplusplus 
extern "C" 
{
#endif /*__cplusplus*/ 
extern int    imx_inv_3d_p(grphic3d *im1, grphic3d *imres) ;
#ifdef __cplusplus 
}
#endif /*__cplusplus*/ 
