/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*-----------------------------------------------------------------------
***
***     file:           reca_3d.h
***
***     project:        Imagix 1.01
***
***
***     description:    Recalage d' images 3D 
***
***
***     Copyright (c) 1993, ULP-IPB Strasbourg.
***     All rights are reserved.
***
***   Last action by user  C.NIKOU  on 10 May 1996
***
***---------------------------------------------------------------------*/


#ifndef CBIZARRE
#define CBIZARRE


#ifdef LINUX
#define MAX_RAND_3D                	(pow(2.0,31.0)-1.0)

#else 
#define MAX_RAND_3D                  	(32767.0) 
#endif

#define EPSILON_3D                (10e-9)  
#define MAX_ITER_3D                (100)    
#define TOO_SMALL_3D               (10e-20)



 /****** tx,ty,tz,sx,sy,sz,thetax,thetay,thetaz  ******/

#define NB_PARAM_MT_3D               (9)

#define TX_3D            (0)
#define TY_3D            (1)
#define TZ_3D            (2)
#define ZOOM_X_3D        (3)
#define ZOOM_Y_3D        (4)
#define ZOOM_Z_3D        (5)
#define THETA_X_3D       (6)
#define THETA_Y_3D       (7)
#define THETA_Z_3D       (8)

#define MARGE_T_X_3D                   (5.0)
#define MARGE_T_Y_3D                   (5.0)
#define MARGE_T_Z_3D                   (5.0)
#define MARGE_THETA_X_3D             (PI/30.0)  
#define MARGE_THETA_Y_3D             (PI/30.0)  
#define MARGE_THETA_Z_3D             (PI/30.0)  
#define MARGE_SCALE_X_3D               (0.05)
#define MARGE_SCALE_Y_3D               (0.05)
#define MARGE_SCALE_Z_3D               (0.05)

 /** Nombre des iterations pour le calcul des hyperparametres ***/

#define MAX_REC_HYPER_3D           (50) 
#define MAX_ICM_HYPER_3D           (50)     

 /** Nombre max d' iterations pour le recuit simule(icm) simple ***/

#define MAX_RECUIT_ITER_MT_3D       (25) 
#define MAX_RECUIT_ITER_B_3D       (2*NB_VEC_PRO)


#define NB_CONFIG_MAX_3D                 (200) 
#define NB_CONFIG_T_3D                   (100) 
#define NB_CONFIG_ZOOM_3D                (100) 
#define NB_CONFIG_THETA_3D               (100) 


#define DIMIN_ENERG_3D            (0.05)

#define TEMPERATURE_INITIALE_3D                (1000.0)
#define COEFF_DE_TEMPERATURE_INITIALE_3D                (1.0)

#define TEMP_COEFF_BITMAP_3D                             (0.05)
#define TEMP_COEFF_LEAST_SQUARES_3D                      (0.1)
#define TEMP_COEFF_RATIO_3D                              (0.05)
#define TEMP_COEFF_CORRELATION_3D                        (0.1)
#define TEMP_COEFF_GEMAN_ROBUST_3D                       (0.01)
#define TEMP_COEFF_TUKEY_BIWEIGHT_ROBUST_3D              (0.01)
#define TEMP_COEFF_TRUNCATED_LEAST_SQUARES_ROBUST_3D     (0.1)

#define TAUX_DE_REFROIDISSEMENT_RECUIT_MT_3D               (0.5)
#define TAUX_DE_REFROIDISSEMENT_RECUIT_B_3D                (0.6)


#define SEED_3D                        (0) 
#define MAX_EXPONENT_3D             (1000.0)
#define MIN_EXPONENT_3D             (-1000.0)
#define ENERGIE_MAX_3D             (256.0*256.0)
#define MARGE_3D                         (10)


 /*** Configuration building *****/ 

#define REGULAR_3D		        (0)
#define IRREGULAR_3D	        (1)

 
#define ZOOM_OFF_3D			0
#define ZOOM_ON_3D			1

 
 
 /** Energy function type ***/

#define BITMAP_ENERGY_3D             		(0)
#define SQUARES_DIFF_ENERGY_3D       		(1)
#define WOODS_RATIO_ENERGY_3D        		(2)
#define WOODS_ROBUST_RATIO_ENERGY_3D        	(3)
#define ENTROPY_ENERGY_3D			(4)
#define MUTUAL_INFO_3D				(5)

#define LAMBDA_ENTRO_3D				(0.5)

#define PARAM_GRAD_3D           (3)

#define CONS_NORMA_3D        ((RIGHT_T-LEFT_T)/(RIGHT_THETA-LEFT_THETA))


#define MAX_ITER_GRADIENT_3D                    (10)
#define CRITERE_GRADIENT_CONJUGEE           (1.0e-6) 


 /******* Type de rotation ********/


#define TRILINEAR_ROTATION_3D                  (0) 
#define NEAREST_NEIGHBOUR_ROTATION_3D          (1)
#define THREE_STEP_ROTATION_3D                 (2)
#define SINC_ROTATION_3D                       (3)


#define SINC_WINDOW_3D                         (3)


 /***** Type d l'interpolation pour la rotation a 3 pas (3 step)  *******/

#define SINC_INTERPOLATION_3D               (1)
#define KEY_CUBIC_INTERPOLATION_3D          (2)
#define TRIANGULAR_INTERPOLATION_3D         (3)




 /** Low-pass filter dimension for pyramidal decomposition - always odd for computing purposes **/

#define ALPHA_GAUSS_PARAM_3D               (0.4) 
#define ALPHA_LAPLACE_PARAM_3D             (0.5) 
#define PYRAMID_FILTER_DIM_3D              (5) 
#define MAX_PYRAMID_LEVELS_3D              (3) /* total levels */
#define WITH_INTERPOLATION_3D              (1) 
#define WITHOUT_INTERPOLATION_3D           (2) 

/***** Pour raffiner et limiter la recherche *******/
#define LIMITED_SEARCH_SPACE_3D           (30.0/100.0) 
#define PYRAMID_LIMITED_SEARCH_SPACE_3D           (20.0/100.0) 



/***** Thresholds used to create bitmap image   *****/

#define LOW_THRESHOLD_3D                  (40)
#define HIGH_THRESHOLD_3D                (250)

#define THE_IMAGES_ARE_THRESHOLDED_3D		(0)
#define THE_IMAGES_ARE_NOT_THRESHOLDED_3D      (1)

#define MAX_VOLUME_SAMPLING_STEP          (81) /* pow(3,4+1) */
#define MAX_PYRAMID LEVELS                (5)   /* 4+1 */ 


#define NB_ENERGY_MODIF_3D                             2
#define NB_FINAL_ENERGY_MODIF_3D                       2
#define NB_PYRAMID_ENERGY_MODIF_3D                     1


#define NB_PARAM_ZOOM_3D                  (3)


#define WOODS_MAX_GREY_3D			(256)
#define ENTROPY_MAX_GREY_3D			(256)
#define RUN_OF_LENGTH_3D			(2)

#ifndef WOODS_SIGMA_ROBUST_3D
#define WOODS_SIGMA_ROBUST_3D			(5)
#endif

typedef struct vector_of_3D_points 
 {
  double x;
  double y;
  double z;
 }Points_3d; 






 double  _const_norma_3d; 

 int    _x_image_center_3d; 
 int    _y_image_center_3d; 
 int    _z_image_center_3d; 


 int    _nb_config_tab_3d[NB_PARAM_MT_3D];



extern double _limited_search_space_3D[NB_PARAM_MT_3D];
extern double _limited_search_space_recuit_3D[NB_PARAM_MT_3D];

 
/* AUTOMATICALLY EXTRACTED PROTOTYPES (MB 6/12/2000) */
#ifdef __cplusplus 
extern "C" 
{
#endif /*__cplusplus*/ 

extern double find_angle(double cosa, double sina) ;
extern void Rec_Recu_3d(void) ;
extern void Rec_Recu_3d_Z(void) ;
extern void Rec_Icm_3d(void) ;
extern void Rec_Icm_3d_Z(void) ;
extern void zoom_image_3d(int interpollation);
extern void imx_zoom_image_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres,int interpolation);
extern void modify_image_resolution_3d (grphic3d *image, int max);
extern void transformation_rigide_3d(grphic3d *image, grphic3d *imres, double tx, double ty, double tz,
									 double sx, double sy, double sz, double thetax, double thetay, double thetaz,
									 int xg, int yg, int zg, int TYPE);
extern void calcul_barycentre_de_roi_3d_p(grphic3d *roi, double *x, double *y, double *z);
extern void imx_compute_uniformity_3d (int im_1, int im_2, int im_histo, char *hfname, int maxg);
extern void imx_compute_uniformity_3d_p (grphic3d *im1, grphic3d *im2, grphic *imhisto, char *hfname, int maxg);
extern void compute_uniformity_3d(void);
extern int rotation_translation_scaling_3d (grphic3d *image, grphic3d *imres, double *param_mt, int xg, int yg, int zg, int TYPE, int sampling_step);
extern void initialiser_param_mt_3d (grphic3d *imref, grphic3d *imreca, double *param, int pr_ax, int *xg, int *yg, int *zg);
extern int rotation_translation_scaling_3d_nikou (grphic3d *image, grphic3d *imres, double *param_mt, XPoints_3d cgdeb, XPoints_3d cgres, int TYPE, int sampling_step);
extern void imx_rec_pr_ax_3d_p (grphic3d *imref, grphic3d *imreca, grphic3d *imres, int *xg, int *yg, int *zg, double *tx, double *ty, double *tz, double *thetax, double *thetay, double *thetaz);
extern void subtract_roi_3d (grphic3d *im1, grphic3d *im2, grphic3d *imres);
extern void imx_rec_pr_ax_3d (int im_ref, int im_reca, int im_res);
extern void imx_zoom_image_3d(int im_1, int im_2, int im_res,int interpolation);
extern void enregis_param_reca_mt_3d(double *param_mt, int xg, int yg, int zg, char *file_param);


#ifdef __cplusplus 
}
#endif /*__cplusplus*/ 

#endif
