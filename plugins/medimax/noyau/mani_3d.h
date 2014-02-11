/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/


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
extern int  imx_shift_bou_3d(void) ;
extern int  imx_shift_bou_3d_p(grphic3d *im1, grphic3d *imres, int shft_x, int shft_y, int shft_z) ;
extern int  imx_copie_3d(int im_deb, int im_res) ;
extern int  imx_copie_img_mask_3d(int im_deb, int im_res) ;
extern int  imx_copie_positif_3d(int im_deb, int im_res) ;
extern int  imx_copie_positif_min_3d(int im_deb, int im_res, int wnd) ;
extern int DllExport imx_copie_3d_p(grphic3d *imdeb, grphic3d *imres) ;
extern int  imx_copie_positif_3d_p(grphic3d *imdeb, grphic3d *imres) ;
extern int  imx_copie_positif_min_3d_p(grphic3d *imdeb, grphic3d *imres,int wnd) ;
extern int  imx_copie_byte_3d (int im_deb, int im_res);
extern int  imx_copie_byte_maxmin_3d (int im_deb, int im_res);
extern int  imx_maxrcoeff_multicoupe_p(grphic3d *imres, int nbslice, double *tabrcoeff3d, float max, int read_format_3d) ;
extern void modify_pixsize_3d(void) ;
extern void modify_norrcoeff_3d(void) ;
extern void imx_modify_norrcoeff_3d_p(grphic3d* imchoix,int valeur);
extern int  rot_droite_3d(int im_deb, int im_res) ;
extern int  rot_gauche_3d(int im_deb, int im_res) ;
extern int  rot_bas_3d(int im_deb, int im_res) ;
extern int  rot_haut_3d(int im_deb, int im_res) ;
extern int  rot_pro_dr_3d(int im_deb, int im_res) ;
extern int  rot_pro_gau_3d(int im_deb, int im_res) ;

//---repositionnement de l'image---//
extern void repositionner_image_3d();
extern int imx_repositionner_image_3d(int im_deb, int im_res, int vue_sagittale, int pos_nez, int pos_sommet_crane);

extern grphic3d *imx_adapte_taille_3d(int im_1, grphic3d *imres, int wdth, int hght, int dpth) ;
extern grphic3d *imx_adapte_taille_3d_p(grphic3d *im1, grphic3d *imres, int wdth, int hght, int dpth);

extern int Miroir_ttaxe_3d(int im_deb,int im_res,int axe); // modif hennequin 07/04/2005
extern int Miroir_ttaxe_3d_p(grphic3d *im_deb,grphic3d *im_res,int axe); // ajout hennequin 07/04/2005
extern int Miroir_3d_x(int im_deb,int im_res);
extern int Miroir_3d_y(int im_deb,int im_res);
extern int Miroir_3d_z(int im_deb,int im_res);
extern void Miroir_3d(void);

extern int imx_histogram_equalization_3d_p (grphic3d *im1, grphic3d *imres, grphic3d *roi, float *hist, long int Nbpas);
extern void imx_cut_head_3d (int im_deb, int im_res, double hd_hght);
extern void imx_cut_head_3d_p (grphic3d *imdeb, grphic3d *imres, double hd_hght);

extern void Rot_Dr_3d(void) ;
extern void Rot_Ga_3d(void) ;
extern void Rot_Ba_3d(void) ;
extern void Rot_Ha_3d(void) ;
extern void Rot_Pro_Dr_3d(void) ;
extern void Rot_Pro_Gau_3d(void) ;
//extern void Rot_3d(void) ;
extern void Move_3d_Inc(void) ;
extern void Copie_3d(void) ;
extern void Copie_positif_3d(void) ;
extern void Copie_positif_min_3d(void) ;
extern void Copie_F1_3d(void) ;
extern void Copie_imgmask_3d(void) ;
extern int  imx_Copie_imgmask_3d(int im_deb, int im_res);
extern int  imx_Copie_imgmask_3d_p(grphic3d *imdeb, grphic3d *im_res);
extern int	imx_copie_f1_3d(int im_deb, int im_res);
extern int	imx_copie_f1_3d_p(grphic3d *imdeb,grphic3d *imres);
extern void cent_3d(void) ;
extern int  imx_cent_3d(int im_deb, int im_res, long int *decal_x, long int *decal_y, long int *decal_z ) ;
extern int  imx_cent_3d_p(grphic3d *imdeb, grphic3d *imres, long int *decal_x, long int *decal_y, long int *decal_z ) ;
extern void bary_3d(void) ;
extern int  imx_bary_3d(int im_deb, int im_res, long int *bary_x, long int *bary_y, long int *bary_z) ;
extern int    imx_bary_3d_p(grphic3d *imdeb, grphic3d * imres, long int *bary_x, long int *bary_y, long int *bary_z);
extern void calcul_barycentre_de_img_3d_p(grphic3d *img, long int *x, long int *y, long int *z);
extern void clean_board_3d(void) ;
extern int  imx_clean_board_3d(int im_deb, int im_res, int nbpts) ;
extern int  imx_clean_board_3d_p(grphic3d *imdeb, grphic3d *imres, int nbpts) ;
extern void imx_equalize_slices_3d(int im_1, int im_2);
extern void imx_equalize_slices_3d_p(grphic3d *im1, grphic3d *im2);
extern void  equalize_dxdydz_3d(void);
extern void imx_equalize_cutoff_3d(void);
extern void set_size_3d(void);
extern void set_size_cg_3d(void);
extern void set_size_cg_trf_3d(void);
extern void set_x3dcy3dcz3dc_3d(void);
extern void imx_set_size_3d(int im_1,int nw,int nh,int nd);
extern void imx_set_size_3d_p(grphic3d *im1,int nw,int nh,int nd);
extern void imx_set_size_cg_3d(int im_1,int nw,int nh,int nd);
extern void imx_set_size_cg_3d_p(grphic3d *im1,unsigned int nw,unsigned int nh,unsigned int nd);
extern void imx_set_size_cg_trf_3d(int im_1,int nw,int nh,int nd,char *filename);
extern void imx_set_x3dcy3dcz3dc_3d_3d_p(grphic3d *im1);
extern void imx_energie_intra_image (int im_1);
extern void imx_energie_inter_images_3d (int im_1, int im_2);
extern void imx_energie_intra_image_3d (int im_1);
extern double imx_energie_inter_images_3d_p(grphic3d *im1, grphic3d *im2);
extern double imx_energie_intra_image_3d_p(grphic3d *im1);
extern void energie_inter_images_3d ( void);
extern void energie_intra_image_3d ( void);
extern void imx_calcul_barycentre_de_img_3d_p(grphic3d *im1, long int *x, long int *y, long int *z);
extern void imx_barycentre_pondere_3d_p(grphic3d *img, int seuil, float *x,float *y, float *z);
extern int       decalage_3d(int im_dep, int im_res) ;
extern int   imx_decalage_3d(int im_deb, int im_res, int deci, int decj, int deck) ;
extern int    imx_decalage_3d_p(grphic3d *imdeb, grphic3d *imres, int deci, int decj, int deck) ;
extern void rescal3d(grphic3d *im1resc) ;


#ifdef __cplusplus 
}
#endif /*__cplusplus*/ 
