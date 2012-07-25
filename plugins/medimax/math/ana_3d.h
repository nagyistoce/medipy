/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/


//extern int imx_histogram_3d_p(grphic3d *im1, grphic3d *roi, float xmin, float xmax, long int Nbpas, float *x, float *y);

 
/* AUTOMATICALLY EXTRACTED PROTOTYPES (MB 6/12/2000) */
#ifdef __cplusplus 
extern "C" 
{
#endif /*__cplusplus*/ 
	
	extern void   select_cell_3d(void);
  extern void   select_multi_cell_3d(void);
	extern int	  imx_select_cell_3d(int im_1, int im_res,int choix);
	extern int 	  imx_select_cell_3d_p(grphic3d *im1, grphic3d *imres,int choix);
	extern int    imx_select_cell_value_3d(int im_1, int im_res, long int valeur) ;
	extern int    imx_select_cell_value_3d_p(grphic3d *im1, grphic3d *imres, long int valeur) ;
	extern long    *fct_volume_cell_3d_p(grphic3d *im1) ;
	extern int	statistique_3d(int pos3d, float amaxh, float aminh);
	extern int  imx_statistique_3d(int im_1, int roi_, float aminh, float amaxh, int *inb, float *moy, float *ecty, float *skew, float *kurt, float *mediane, float *per5quantile, float *per95quantile) ;
	extern int  imx_statistique_3d_p(grphic3d* im_1, grphic3d *mask3d, float aminh, float amaxh, int *inb, float *moy, float *ecty, float *skew, float *kurt, float *mediane, float *per5quantile, float *per95quantile) ;

	int    imx_select_cell_value_3d_p(grphic3d *im1, grphic3d *imres, long int valeur);
	extern int	imx_histogram_3d_p(grphic3d *im1, grphic3d *roi, float xmin, float xmax, long int Nbpas, float *x, float *y);
  
  extern int	  imx_select_multi_cell_3d(int im_1, int im_res,int choix);
  extern int 	  imx_select_multi_cell_3d_p(grphic3d *im1, grphic3d *imres,int choix);
#ifdef __cplusplus 
}
#endif /*__cplusplus*/ 
