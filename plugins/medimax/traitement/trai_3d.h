/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/


#ifndef __TRAI_3D_H__
#define __TRAI_3D_H__

#ifdef WIN32
#define DllExport   __declspec( dllexport )  // for windows
#else
#define DllExport  // Linux
#endif

/* AUTOMATICALLY EXTRACTED PROTOTYPES (MB 6/12/2000) */
#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/
  extern void ctr_labele_3d(void);
  extern void ctr_labele_rempl_3d(void);
  extern void d_euclidean_mapping_3d(void) ;
  extern void sync_labels_3d(void);
  extern int imx_sync_labels_3d(int im_src, int im_ref, int im_res);
  extern int imx_sync_labels_3d_p(grphic3d *im_src, grphic3d *im_ref, grphic3d *im_res);
  extern void sync_labels_sep_3d(void);
  extern int imx_sync_labels_sep_3d(int im_src, int im_ref);
  extern int imx_sync_labels_sep_3d_p(grphic3d *im_src, grphic3d *im_ref);
  extern int imx_labelcroix_3d(int im_1, int im_res);
  extern int imx_labelcroix_3d_p(grphic3d *im1, grphic3d *imres) ;
  extern DllExport int imx_labelcube_3d_p(grphic3d *im1, grphic3d *imres) ;
  extern int imx_ordonner_labels_3d(int im_1, int im_res);
  extern int imx_ordonner_labels_3d_p(grphic3d *im1, grphic3d *imres) ;
  extern int imx_ctr_binaire_3d (int im_1, int im_res);
  extern int imx_ctr_binaire_3d_p(grphic3d *im1, grphic3d *imres) ;
  extern int imx_ctr_labele_3d (int im_1, int im_res);
  extern int imx_ctr_labele_3d_p(grphic3d *im1, grphic3d *imres) ;
  extern int imx_ctr_labele_rempl_3d (int im_1, int im_res);
  extern int imx_ctr_labele_rempl_3d_p(grphic3d *im1, grphic3d *imres) ;
  extern void sub_3d(void) ;
  extern void subabs_3d(void) ;
  extern void sub_normamax2_3d(void) ;
  extern int imx_sub_normamax2_3d(int im_1, int im_2, int im_res);
  extern int imx_sub_normamax2_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres);
  extern void sub_normpixel2_3d(void) ;
  extern int imx_sub_normpixel2_3d (int im_1, int im_2, int im_res);
  extern void sub_normpixelmean2_3d(void) ;
  extern void pix_compte_3d(void) ;
  extern void shell_sort(int n, float *data) ;
  extern float kolmogorov_proba(float lambda) ;
  extern float imx_sigma_from_VOI_p(grphic3d *imref, int start_x, int start_y, int start_z, int size) ;
  extern int labelcube_3d(void);
  extern int imx_labelcube_3d (int im_1, int im_res);
  extern void imx_fusion_moy_3d(int pos_img_sup, int pos_img_fond);
  extern void imx_fusion_mix_3d(int pos_img_sup, int pos_img_fond);
  extern void imx_sigma_median_3d_p (grphic3d *imref, float *median);
  extern void imx_sigma_median_3d(int im_ref, float *median);
  extern int	fct_grow_label_3d(int i0, int j0, int k0, int i, int j, int k, grphic3d *im1, float *tab, int n);
  extern int	fct_grow_con_pos_3d(int i0, int j0, int k0, int i, int j, int k, grphic3d *im1, float *tab, int n);
  extern int	imx_sub_normpixelmean2_3d_p(grphic3d * im_1, grphic3d * im_2, grphic3d * im_res);
  extern int	imx_sub_normpixelmean2_3d(int im_1, int im_2, int im_res);
  extern void remplit_objet(grphic3d *im1, int obj);

  extern void imx_d_euclidean_mapping_3d (int im_deb, int im_res_x, int im_res_y, int im_res_z, int im_res_m);
  extern void imx_d_euclidean_mapping_3d_p (grphic3d *imdeb, grphic3d *imresx, grphic3d *imresy, grphic3d *imresz, grphic3d *imresm);
  extern void border_processing_3d (grphic3d *imresx, grphic3d *imresy, grphic3d *imresz, grphic3d *imresm);

//---------------------
  extern void standardize_noise_3d(void);
  extern int imx_standardize_noise_3d(int im_ref, int im_res, float sigma);
  extern void imx_standardize_noise_3d_p(grphic3d *imref, grphic3d *imres, float sigma);

//---methodes de vraisemblance statistique pour eliminer le bruit---//
  extern int compmanuel_3d(void) ;
  extern int imx_compmanuel_3d(int im_1, int im_2, int im_res, int N, double seuil) ;
  extern int imx_compmanuel_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres, int N, double seuil) ;
  extern int compauto_3d(void) ;
  extern int imx_compauto_3d(int im_1, int im_2, int im_res, int N, int answer_nr, float K) ;
  extern int imx_compauto_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres, int N, int answer_nr, float K) ;

//---test de kolmogorov---//
  extern void kolmogorov_3d(void) ;
  extern int imx_kolmogorov_3d(int im_1, int im_2, int im_res, int N, float fausse_alarme, int type_remplissage) ;
  extern int imx_kolmogorov_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres, int N, float fausse_alarme, int type_remplissage) ;

//---test du maximum de vraisemblance---//
  extern void vraisemblance_3d(void) ;
  extern int imx_vraisemblance_3d(int im_1, int im_2, int im_res, int N, double fausse_alarme, int type_remplissage) ;
  extern int imx_vraisemblance_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres, int N, double fausse_alarme, int type_remplissage) ;

//--- Comparaison par Rang ---//
  extern void comparaison_par_rang(void);
  extern void imx_comparaison_par_rang(int im_ref, int im_cmp, int im_res, int methode);
  extern void imx_comparaison_par_rang_p(grphic3d *imref, grphic3d *imcmp, grphic3d *imres);
  extern void imx_comparaison_par_int_dis_p(grphic3d *imref, grphic3d *imcmp, grphic3d *imres);
  //extern void imx_rechelonage_p(grphic3d *imref, grphic3d *imres);

//---test du maximumde vraisemblance---//
  extern void matched_filter_3d(void) ;
  extern void unit_matched_filter_3d(void) ;
  extern int imx_matched_filter_3d(int im_1, int im_2, int im_res, int N, double fausse_alarme, int type_remplissage, double sigma) ;
  extern int imx_matched_filter_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres, int N, double fausse_alarme, int type_remplissage, double sigma) ;

//---extraction du cerveau---//
  extern int imx_extract_cerv_3d_p (grphic3d *im1, grphic3d *imres);
  extern int imx_extract_cerv_3d (int im_1, int im_res);

//---segmentation par growing---//
  extern int	growing_3d_p(int x, int y, int z, grphic3d *im1, grphic3d *imres, int (*fct_grow_3d) (/* ??? */), float *tableau, int ntab, long int valeur, int typevois);
  extern int imx_growing_cerv_3d_p(grphic3d *im1, grphic3d *imres, int x0, int y0, int z0) ;
  extern int growing_cerv_fonction_3d_p(int x, int y, int z, grphic3d *im1, grphic3d *imres, float pourcent);
  int imx_growing_cerv_3d (int im_1, int im_res, int x0, int y0, int z0);

//---test de student---//
  extern void student_3d();
  extern int imx_student_3d(int im_1, int im_2, int im_res, int N, float fausse_alarme, int type_remplissage);
  extern int imx_student_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres, int N, float fausse_alarme, int type_remplissage);

//---distance de chamfer---//
  extern void chamfer_distance_3d();
  extern void imx_chamfer_distance_3d (int im_deb, int im_res);
  extern void imx_chamfer_distance_3d_p (grphic3d *imdeb, grphic3d *imres);

//---remplissage de trous---//
  extern int imx_hole_fill_3d (int im_1, int im_res);
  extern int imx_hole_fill_3d_dir (int im_1, int im_res);
  extern int imx_hole_fill_3d_cpc (int im_1, int im_res, int type);
  extern void imx_hole_fill_3d_p (grphic3d *im1, grphic3d *imres);
  extern void imx_hole_fill_3d_dir_p (grphic3d *imdeb, grphic3d *imres);
  extern void imx_hole_fill_3d_cpc_p (grphic3d *imdeb, grphic3d *imres, int type);
  extern int hole_fill_3d_cpc(int type);
  extern void hole_fill_3d();
  extern void hole_fill_3d_dir();

//---seuillage par hysteresis---//
  extern void imx_hysteresis_thresholding_3d_p (grphic3d *imdeb, grphic3d *imres, int cc_length);
  extern void imx_hysteresis_thresholding_3d (int im_deb, int im_res, int cc_length);
  extern void hysteresis_thresholding_3d();

//-------------------FILTRAGE----------------------//

//---canny deriche---//
  double filter_d(int x, double alpha, double c);
  double filter_l(int x, double alpha, double s);
  void imx_canny_deriche_3d (int im_deb, int im_res, double alpha, int window, int cc_length);
  void imx_canny_deriche_3d_p (grphic3d *imdeb, grphic3d *imres, double alpha, int window, int cc_length);
  void compute_gx_3d (grphic3d *imdeb, grphic3d *gx, double alpha, int wnd);
  void compute_gy_3d (grphic3d *imdeb, grphic3d *gy, double alpha, int wnd);
  void compute_gz_3d (grphic3d *imdeb, grphic3d *gz, double alpha, int wnd);
  void compute_module_3d (grphic3d *gx, grphic3d *gy, grphic3d *gz, grphic3d *imres);
  void gradient_maxima_3d (grphic3d *gx, grphic3d *gy, grphic3d *gz, grphic3d *imres);
  void initialize_gradient_3d (grphic3d *gx, grphic3d *gy, grphic3d *gz);
  int interpol_gradient_3d (grphic3d *gx, grphic3d *gy, grphic3d *gz, double x, double y, double z);
  extern void  canny_deriche_3d();

//----norme quadratique carree du gardient d'un tableau de double
  extern void norme_sqr_quad_gradient();
  extern int imx_norme_sqr_quad_gradient_3d(int im_src, int im_dst);
  extern int imx_norme_sqr_quad_gradient_3d_p(double *image, double *gradient, int width, int height, int depth);

//----filtrage moyenneur
  extern void moyen_filter_3d(void);
  extern int imx_moyen_filter_3d (int im_deb, int im_res, int wnd);
  extern int imx_moyen_filter_3d_p(grphic3d *imdeb, grphic3d *imres, int wnd);


//----filtrage median
  extern void median_3d(void);
  extern int imx_median_3d (int im_deb, int im_res, int wnd);
  extern int imx_median_3d_p(grphic3d *imdeb, grphic3d *imres, int wnd);

//----filtrage gaussien
  extern void gaussian_filter_3d(void);
  extern void imx_gaussian_filter_3d (int im_deb, int im_res, int wnd, double sigma);
  extern void imx_gaussian_filter_3d_p(grphic3d *imdeb, grphic3d *imres, int wnd, double sigma) ;
  extern int imx_filtre_gaussien_anisotrope_3d_p(grphic3d *imdeb, grphic3d *imres, unsigned int wnd, double sigmax, double sigmay, double sigmaz);
  extern void imx_gaussian_filterModif_3d(int im_res, int wndX,int wndY,int wndZ,double,double,double) ;
  extern void imx_gaussian_filterModif_3d_p(grphic3d *imres, int wndX,int wndY,int wndZ,double,double,double) ;
  extern void sc_gaussian_filterModif(float *y, long int nbpts, int wnd, double sigma,double * filter);

//----filtre approximant une equation de diffusion
  extern void diffusion_filter_3d(void) ;
  extern void imx_diffusion_filter_3d(int im_deb, int im_res, unsigned int nb_it);
  extern void imx_diffusion_filter_3d_p(grphic3d *imdeb, grphic3d *imres, unsigned int nb_it);

//----filtre moyenneur
  extern void filtre_moyenneur_3d();
  extern int imx_filtre_moyenneur_3d(int im_src, int im_dst, int taille_filtre);
  extern int imx_filtre_moyenneur_3d_p(grphic3d *imsrc, grphic3d *imdst, int taille_filtre);

//----filtre de diffusion anisotrope
  extern void filtre_diffusion_anisotrope_3d();
  extern int imx_filtre_diffusion_anisotrope_3d(int im_src, int im_dst, int nb_iter);
  extern int imx_filtre_diffusion_anisotrope_3d_p(grphic3d *imsrc, grphic3d *imdst, int nb_iter);

  extern int convoluer_filtre_1D(double *ligneIn, double *ligneOut, int tailleLigne, double *filtre, int tailleFiltre);
  extern int convoluer_filtre_1D_sans_bord(double *ligneIn, double *ligneOut, int tailleLigne, double *filtre, int tailleFiltre);

//---misc---//

extern int image_pixel_3d (int x, int y, int z, grphic3d *im1resc);
void find_sign_3d (grphic3d *imdeb, grphic3d *imresx, grphic3d *imresy, grphic3d *imresz);

extern int imx_query_filter_dimension();

void histogramme_cumule_dec_p(int *histo, grphic3d *image, int size, int *min, int *amplitude);

#ifdef __cplusplus
}
#endif /*__cplusplus*/


#endif //__TRAI_3D_H__
