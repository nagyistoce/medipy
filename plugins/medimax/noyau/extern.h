/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef __EXTERN_H
#define __EXTERN_H

/*  #ifdef __cplusplus */
/*  extern "C" { */
/*  int    is_mri(const char *, const char *); */
/*  } */
/*  #else */
/*protoize:???*/ /* pour la definitionde HEADER */
#include "noyau/io/imx_head.h"
int    is_mri        (char *path, char *file);/* Determine if a file is a mri.    */
long   *getmri        (char *file, int numfirst);/* Load in memory a bruker picture.    */
long    getmri_index3D  (char *file, int index2D);/* Converts a 2D mri index in a 3D mri index */
long    getmri_numberimg(char *file);/* Gets the number of image in a file*/
long    getmri_width     (char *file, int number, int mode);/* Calculate mri width.        */
long    getmri_height    (char *file, int number, int mode);/* Calculate mri height.        */
long    getmri_depth    (char *file, int number, int mode);/* Calculate mri depth.        */
long    getmri_minpixel    (char *file, int number, int mode);/* Calculate mri minpixel.        */
long    getmri_maxpixel    (char *file, int number, int mode);/* Calculate mri maxpixel.        */
long    getmri_cutoff_min (char *file, int number, int mode);/* Calculate mri cutoff min.    */
long    getmri_cutoff_max (char *file, int number, int mode);/* Calculate mri cutoff max.    */
float    getmri_icomp    (char *file, int number, int mode);/* Calculate mri icomp.        */
float    getmri_dx(char *file, int number, int mode); 
float    getmri_dy(char *file, int number, int mode); 
float    getmri_dz(char *file, int number, int mode); 
int      getmri_mask_type(char *file, int number, int mode); 
long    getmri_dattype  (char *file, int number, int mode);/* Calculate mri icomp.        */
char   *getmri_name (char *file);
char   *getmri_dbirth(char *file);
char   *getmri_dexamen(char *file);
char   *getmri_medecin(char *file);
char   *getfilename_d3proc(char *file);
char   *getfilename_2dseq(char *file);
char   *getfilename_procs(char *file);
char   *getfilename_reco(char *file);
char   *getfilename_subject(char *file);
char   *getfilename_imnd(char *file);
char   *getfilename_acqp(char *file);
char   *getheader_interfile(char *file, char *stringbis, int type, int number);
long    putheader_interfile(char *file, char *string, int type, void *itemadr, int pos);
int     put_header      (HEADER *header, char *string, int type, void *itemadr, int pos);
int     load_header     (char *file, HEADER *header);/* For putheader_interfile function */
int     save_header     (HEADER *header, char *file);/* For putheader_interfile function */
int    sp_cpymri_picture(long int *dest, long int *src, long int wdth, long int hght);
long    putmri_icomp    (char *src_file, char *dst_file, long int icomp);
long    putmri_maxpix    (char *src_file, char *dst_file, long int maxpix);
long    putmri_minpix    (char *src_file, char *dst_file, long int minpix);
long    putmri_width    (char *src_file, char *dst_file, long int width);
long    putmri_height    (char *src_file, char *dst_file, long int height);
long    buff2pic    (long int *buff, long int **pic, grphic *image);
grphic *load_mri    (char *file, grphic *graphic, int numfirst);/* Load pictures.            */
void    refresh_image_display(grphic *ima, int x, int y, int w, int h); 
int     check_file_2d    (char *file);/* Check file for save */

void     imx_info_mri_2d    (int nr);
void     imx_info_mri_3d    (int nr);
void    ChangeGC    (int pos);
/* Build an image that can be used by XCreateImage from a grphic image */ 
int    sp_cpymri_picture(long int *dest, long int *src, long int wdth, long int hght);
int    cpymri_picture    (long int *dest, long int *src, long int wdth, long int hght);
void    charge_petite_picture(int num_image, char *file, int numfirst, int n, int posdsimg, int norm);
long    fread_long      (long int *t, int nitem, FILE *fp);
long    fread_long2      (long int *t, int nitem, FILE *fp);
long    fread_short     (short *t, int nitem, FILE *fp);
long    fread_short2     (short *t, int nitem, FILE *fp);
long    fwrite_long    (long int *t, int nitem, FILE *fp);

/*protoize:???*/ /* pour la definition de : byte*/
#include"noyau/io/imx_file.h"
void put_byte           (byte b, FILE *hand);
void put_word           (int i, FILE *hand);
void put_line           (byte *line, FILE *hand);
FILE* p_openr        (char *name);

/*  #endif */

#endif /*  __EXTERN_H */
