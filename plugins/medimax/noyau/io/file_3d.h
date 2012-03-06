/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/**-----------------------------------------------------------------------
***	
***	file:		imx_3d.h
***
***	project:	Imagix 1.01
***			
***
***	description:	Header file for project imagix (1.01)
***	
***	
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***     Last user action by: Mr. ARMSPACH on Juny  1995
***
***---------------------------------------------------------------------*/

#ifndef __FILE_3D_H
#define __FILE_3D_H

#ifdef __GTK_GUI
#include "outils/imx_stck.h" /* _drawing_3d utilise les listes */
#endif /* __GTK_GUI */

#ifdef __cplusplus 
extern "C" 
{
#endif /*__cplusplus*/  


  /* --------------------------- All functions --------------------------	*/
#ifdef __GTK_GUI
  extern int    load_picture_3d_cst(int pos, const char *file, int max, int step, int numfirst) ;
  extern int    load_picture_3d(int pos, char *file, int max, int step, int numfirst) ;
  extern int    load_picture_raw_3d(const char *file, int numfirst, int wdth, int hght, int dpth, int pos, int taille) ;
  extern grphic3d *load_mri_bruker_3d(const char *file, grphic3d *graphic, int dpth, int numfirst) ;
#endif /* __GTK_GUI */
  extern int load_mri_ipb_3d(const char *file, grphic3d *image, int number) ;
  extern int save_mri_ipb_3d_p(char *file, grphic3d *image) ;
  extern int    load_mri_3d(const char *file, grphic3d *graphic, int numfirst, int max, int step) ;
 
#ifdef __GTK_GUI
  extern int save_pictures_3d(char *file, int start, int num_images, int step) ;
  extern int save_pictures_3d_std(char *file, int start, int num_images, int step) ;
  extern int save_pictures_3d_erase(char *file, int start, int num_images, int step) ;
  extern int save_mri_ipb_3d(char *file, int start, int num_images, int step) ;
  extern int save_mri_img_3d(const char *file, int start, int num_images, int step);
  extern int save_mri_img_3d_p(const char *file_img, grphic3d *image);
  extern int add_mask_Ipb_3d(char * file, int pos);
  extern int update_tala(int pos, int index, const char *file) ;
  extern long    buff3pic(long int *buff, grphic3d *graphic, int dpth, float anc_rcoef) ;
  extern int read_raw_3d(const char *file, int numfirst, int filewdth, int filehght, int filedpth, grphic3d *img, int format, int typemri) ;
  extern float    pix_value_3d_p(grphic3d *im1, float *X, int TYPE) ;
  extern int     interp_3d(int im_deb, int im_res, int TYPINTERP) ;
  extern int     imx_interp_3d(int im_deb, int im_res, int TYPINTERP) ;
  extern float sin_car(float nb) ;
  extern int     interp_tri_3d(int im_deb, int im_res) ;
  extern int     imx_interp_tri_3d(int im_deb, int im_res) ;
  extern int imx_spline_tab(IMXPoint *tabcoord, int *Nbpts, IMXPoint *tabres, int *Nbptsres);
  extern void Open_3d(void) ;
  extern void nOpen_3d(void) ;
  extern void Save_3d_xml(void) ;
  extern void Save_3d_Ipb(void) ;
  extern void Save_3d_Ipb_std(void) ;
  extern void Save_3d_Ipb_erase(void) ;
  extern void Save_3d_Ipb_mask_erase(void) ;
  extern void Save_3d_analyze(void) ;
  extern void Save_Ipbr_3d(void) ;
  extern void Save_3d_Ipb_mask(void) ;
  extern void nSave_3d_Ipb(void) ;
  extern void Int_Lin_3d(void) ;
  extern void Int_Tri_3d(void) ;
  extern void Int_Sinc_3d(void) ;
  extern void Open_raw_3d(void) ;
  extern void nOpen_raw_3d(void) ;
  extern int scaleimg_3d_p(grphic3d *imdeb, grphic3d *imres) ;
  extern int imx_scaleimg_3d_p(grphic3d *imdeb, grphic3d *imres, int meth, int factx, int facty, int factz) ;
  extern long getPosFile(const char *file, int numfirst, int filewdth, int filehght, int filedpth, int format, int typemri) ;
  extern void essai_spline(void);
  extern int  check_ipb_format_3d(char * file);
//  extern int  check_ipb_format_3d(const char * file);

extern void targz_Ipb_file(char *file);
extern void Save_3d_Ipb_targz(void);
extern void nSave_3d_Ipb_targz(void);
extern void Save_3d_Ipb_mask_targz(void);
extern void delete_ipb_file(char* file);
extern int is_targz(char* file);
extern int untargz_Ipb_file(char *file);

  /*
  ** read_raw_float_3d
  */
  /*! insere le contenu du fichier file dans l'image img.
    \param file : fichier a lire
    \param img : image a remplir
    \param endianess : parametre d'endianess du fichier (ne sert pas
  pour l'instant)  
   */
  extern void read_raw_float_3d(const char* file, grphic3d* img, int endianess);

#endif /* __GTK_GUI */

#ifdef __cplusplus 
}
#endif /*__cplusplus*/ 
#endif /* __FILE_3D_H */
