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

extern int imx_mul_p(grphic *im1, grphic *im2, grphic *imres) ;
extern int imx_mul_coe_p(grphic *im1, float coefficient, grphic *imres) ;

extern int imx_inimaxminpixel_p(grphic *im1);

extern int imx_copie_param_p(grphic *im1, grphic *imres);

extern void op_math(void);
extern void	op_math_coe(void);
extern void	op_log(void);

extern int	imx_iniparaimg(int im_1);
extern int	imx_iniparam(int im_1, long int max_pixel, long int min_pixel, float icomp, float rcoeff);
extern int	imx_inimaxminpixel(int im_1);
extern int	imx_and(int im_1, int im_2, int im_res);
extern int	imx_add_p(grphic *im1, grphic *im2, grphic *imres);
extern int	imx_copie_param(int im_1, int im_res);
extern int	imx_iniparaimg_p(grphic *im1);
extern int	imx_mul(int im_1, int im_2, int im_res);
extern int	imx_add_coe(int im_1, float coefficient, int im_res);
extern int	imx_brukermax(float amax, int *maxpixel, float *icomp, float *rcoeff);
#ifdef __cplusplus 
}
#endif /*__cplusplus*/ 
