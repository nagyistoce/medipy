/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/* -----------------------------------------------------------------------
**    
**     File:       img_file.c
**
**     Project:    Imagix 1.01
**
**     Description:
**     Manages bruker files and pictures in memory
**
**     Copyright (c) 1993, ULP-IPB Strasbourg.
**     All rights are reserved.
**
**     >>  Last user action by: Mr. ETTER      on          1994
**     >>  Last user action by: Mr. ARMSPACH   on Apr 09th 1995
**     >>  Last user action by: Mr. GOUVERNEUR on Nov 24th 1998
**     >>  Last user action by: Mr. BOULEAU    on Nov 18th 2000
**         - moved #includes, #defines, type declarations, etc. from
**           imx_file.c to imx_file.h
**
**---------------------------------------------------------------------*/

#ifndef __IMX_FILE_H
#define __IMX_FILE_H

#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <errno.h>

#define     NOT_SAVE    1

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_2d.h"	   /* Global header for project...    */
#include "noyau/imx_3d.h"	   /* Global header for project...    */
#include "noyau/io/imx_head.h"    /* About header...	*/
#include "noyau/imx_lang.h"    /* */
typedef unsigned char byte;
#include "noyau/extern.h"	   /* */
#ifdef __GTK_GUI
#include "noyau/gui/gtkmmtools.h"
#endif /* __GTK_GUI */

/* Fonction declared in this file:
** -------------------------------
*/


#define PCX_MAGIC           0x0a    /* PCX magic number                       */

#define TAILLE_BLOC_HEADER          4096
#define TAILLE_LIGNE_MIN_HEADER     256
#define INCREMENT_LIGNE_HEADER      100

/*  #ifdef __cplusplus */ 
/*  extern "C" { */
/*  #endif */

/*  #if defined(ANSIC) || defined(__cplusplus) */

extern void ReadGIF(FILE *fd, int imageNumber, int pos, int reduc, char *file);
extern void loadpcx_p(char *file, grphic *graphic, int reduc);
extern int GetByte ( FILE *fp );
extern int GetWord ( FILE *fp );
/*  extern long GetDword(FILE *fp); */
extern char    **Make_files(char *file, int step, int max);
extern long    getmri_readformat(char *file);
extern int     getmri_paraversion(char *file);
extern char* Ecat7_getFromMatrix(char *file, unsigned short ndMatrix, unsigned short offset, unsigned short len, char *data, char *type);
extern char *decode_pos(unsigned char *header, long unsigned int pos, int type);
extern char *Acrnema_getGroup(char *file, unsigned short g1, unsigned short g2, char *data, int format);
extern char *Acrnema_getGroupHeader(char *file, short unsigned int g1, short unsigned int g2, char *data, int format);
extern long	 Acrnema_getGroupSizeandOffset(char *file, short unsigned int g1, short unsigned int g2, long *offset, int format);
extern char* getchar_header(char *fichier,  int offset, int taille);

/*  #else */

extern unsigned long GetDword  (FILE *fp);
extern void    loadpcx_p       (char *file, grphic *graphic, int reduc);
extern int     Acrnema_Dump	();
extern int     Affiche         (void);
extern int     load_picture    (int pos, char *file, int numfirst);
extern int     show_picture    (int pos);     /* Show picture number -1 to MAX_PICTURES */
extern int     load_pictures   (int first_pos, int nb_img, char **pictures);
extern int     save_picture    (int pos, char *file);
extern int     save_pictures   (int first_pos, int nb_img, char **pictures);
extern void     IMXFile_Capture (int type);

/*
**  Ecat7_getFromMatrix : lecture dans les fichiers de type ECAT 7
*/

extern grphic 	*ptr_img            (int nr);
extern void    	Make_AtmFile   	    (char *file, int nb_img, char **pictures);
extern void    	Make_AtmFileII      (char *file, char *wildcard_file);
extern char 	*Acrnema_getGroup	(char *file, short unsigned int g1, short unsigned int g2, char *data, int format);
extern char    *Ecat7_getFromMatrix(char *file, short unsigned int indMatrix, short unsigned int offset, short unsigned int len, char *data, char *type);

/* only to compile program */

extern grphic	        *save_mri_bruker    (char *file, grphic *graphic);
extern int             save_mri_ipb_2d     (char *file, int start, int num_images, int step);
extern long	        *getmri_convert     (long int *buff, void *buffer, int data_size, int format, int num_items);
extern int             load_mri_ipb_2d     (char *file, grphic *image, int pos);
extern char *          getchar_header      (char *fichier, int offset, int taille);
extern double          getdouble_header    (char *fichier, int offset);
extern float           getfloat_header     (char *fichier, int offset);
extern unsigned int    getuint_header      (char *fichier, int offset);
extern long 			getmri_readformat   (char *file);
extern float			mri_point_3d        (long int pixelscreen, long int im_1);
/*  #endif */  /* ANSIC */ 

/*  #ifdef __cplusplus */
/*  } */
/*  #endif */

#endif  /* #ifndef __IMX_FILE_H */
 
#ifdef __cplusplus 
extern "C" 
{
#endif /*__cplusplus*/ 

  extern char   **get_rep(char *rep, int *max);
  extern char   *incfic(char *name, int step);
  extern int incfic_img(char *name, int step, int *numfirst);
  extern char   *incfic_I00X(char *name, int step);
  extern char   *incname_img(char *name);

  extern int     is_mri(char *path, char *file) ;
  extern long   *getmri(char *file, int numfirst) ;
  extern long    getmri_index3D(char *file, int index2D) ;
  extern long    getmri_width(char *file, int number, int mode) ;
#ifdef __GTK_GUI
  extern t_CMapDisplayMode getmri_dspStyle(char *file, int pos, int mode) ;
  extern t_CMapOverBlanckingMode getmri_dspOverBlancking(char *file, int pos, int mode) ;
#endif /* __GTK_GUI */
  extern int getmri_nCMap(char *file, int pos, int mode) ;
  extern int getmri_nCMap2(char *file, int pos, int mode) ;
  extern long    getmri_height(char *file, int number, int mode) ;
  extern long    getmri_depth(char *file, int number, int mode) ;
  extern long    getmri_numberimg(char *file) ;
  extern long    getmri_minpixel(char *file, int number, int mode) ;
  extern long    getmri_maxpixel(char *file, int number, int mode) ;
  extern long    getmri_cutoff_min(char *file, int number, int mode) ;
  extern long    getmri_cutoff_max(char *file, int number, int mode) ;
  extern float    getmri_icomp(char *file, int number, int mode) ;
  extern float    getmri_fov(char *file) ;
  extern float    getmri_dx(char *file, int number, int mode) ;
  extern float    getmri_dy(char *file, int number, int mode) ;
  extern float    getmri_dz(char *file, int number, int mode) ;
  extern int      getmri_mask_type(char *file, int number, int mode) ;
  extern int     getmri_endian(char *file);
  extern long    getmri_dattype(char *file, int number, int mode) ;
  extern long    getmri_readformat(char *file) ;
  extern char*   getmri_name(char *file) ;
  extern char*   getmri_medecin(char *file) ;
  extern char*   getmri_dbirth(char *file) ;
  extern char*   getmri_dexamen(char *file) ;
  extern char    *getfilename_d3proc(char *file) ;
  extern char    *getfilename_procs(char *file) ;
  extern char    *getfilename_2dseq(char *file) ;
  extern char    *getfilename_img(char *file) ;
  extern char    *getfilename_reco(char *file) ;
  extern char    *getfilename_imnd(char *file) ;
  extern char    *getfilename_acqp(char *file) ;
  extern char    *getfilename_subject(char *file) ;
  extern char     *getheader_interfile(char *file, char *stringbis, int type, int number) ;
  extern int     load_header(char *file, HEADER *header) ;
  extern int     put_header(HEADER *header, char *string, int type, void *itemadr, int pos) ;
  extern int     save_header(HEADER *header, char *file) ;
  extern long    putheader_interfile(char *file, char *string, int type, void *itemadr, int pos) ;
  extern long    putmri_icomp(char *src_file, char *dst_file, long int icomp) ;
  extern long    putmri_maxpix(char *src_file, char *dst_file, long int maxpix) ;
  extern long    putmri_minpix(char *src_file, char *dst_file, long int minpix) ;
  extern long    putmri_width(char *src_file, char *dst_file, long int width) ;
  extern long    putmri_height(char *src_file, char *dst_file, long int height) ;
  extern long    buff2pic(long int *buff, long int **pic, grphic *image) ;
  extern grphic *load_mri(char *file, grphic *graphic, int numfirst) ;
  extern void refresh_image_display(grphic *ima, int x, int y, int w, int h) ;
  extern void erase_display(int pos, int x, int y, int w, int h) ;
  extern int check_file_2d(char *file) ;
  extern int check_file_3d(char *file) ;
  extern int check_for_mask(char *file) ;
  extern grphic *save_mri_bruker(char *file, grphic *graphic) ;
  extern int save_pictures_2d(char *file, int start, int num_images, int step, int format) ;
  extern int     save_mri_ipb_2d(char *file, int start, int num_images, int step) ;
  extern int     save_mri_ipb_2d_p(char *file, grphic *image) ;
  extern void save_mri_raw32 (int pos, char *file) ;
  extern void save_mri_raw8 (int pos, char *file) ;
  extern int save_mri_pcx (int pos, char *file) ;
  extern void put_line(byte *line, FILE *hand) ;
  extern void put_word(int i, FILE *hand) ;
  extern void put_byte(byte b, FILE *hand) ;
  extern int     save_mri_nraw8(int first_pos, int nb_img, char **pictures) ;
  extern int     save_mri_nraw32(int first_pos, int nb_img, char **pictures) ;
  extern int     save_mri_npcx(int first_pos, int nb_img, char **pictures) ;
  extern void imx_info_mri_2d(int nr) ;
  extern void imx_info_mri_3d(int nr) ;
  extern void    ChangeGC(int pos) ;
  extern float mri_point_3d(long int pixelscreen, long int im_1) ;
  extern int    sp_cpymri_picture(long int *dest, long int *src, long int wdth, long int hght) ;
  extern int    cpymri_picture(long int *dest, long int *src, long int wdth, long int hght) ;
  extern int    Affiche(void) ;
  extern int    load_picture(int pos, char *file, int numfirst) ;
  extern int show_picture(int pos) ;
  extern int    load_pictures(int first_pos, int nb_img, char **pictures) ;
  extern int    save_picture(int pos, char *file) ;
  extern int    save_pictures(int first_pos, int nb_img, char **pictures) ;
  extern int     capture_picture(int pos, char *file, int type, int r) ;
  extern int    capture_pictures(int first_pos, int nb_img, char **pictures, int type, int r) ;
  extern void    IMXFile_Capture(int type) ;
  extern grphic *GetMire(int type, int width, int height, grphic *mire) ;
  extern int    AffMire(int pos, int type, int width, int height) ;
  extern grphic *ptr_img(int nr) ;
  extern void    Make_AtmFile(char *file, int nb_img, char **pictures) ;
  extern void    Make_AtmFileII(char *file, char *wildcard_file) ;
  extern FILE *p_openr(char *name) ;
  /*BADMATCH extern p_close(FILE *f) ; */
  /*BADMATCH extern allocrow(int cols, int size) ; */
  /*BADMATCH extern allocarray(int cols, int rows, int size) ; */
  /*BADMATCH extern loadpcx(char *file, int pos, int reduc) ; */
  /*BADMATCH extern loadpcx_p(char *file, grphic *graphic, int reduc) ; */
  extern unsigned char read_pcx_file(char *filename, Pcx_Descr *descr, char **status) ;
  /*BADMATCH extern read_pcx_image(FILE *fp, unsigned char *buf, int BytesPerLine, int Planes, int Height) ; */
  /*BADMATCH extern pcx_planes_to_pixels(unsigned char *pixels, unsigned char *bitplanes, int bytesperline, int planes, int bitsperpixel) ; */
  /*BADMATCH extern pcx_unpack_pixels(unsigned char *pixels, unsigned char *bitplanes, int bytesperline, int planes, int bitsperpixel) ; */
  extern int GetByte(FILE *fp) ;
  extern int GetWord(FILE *fp) ;
  extern unsigned long GetDword(FILE *fp) ;
  /*BADMATCH extern p_keymatch(char *str, char *keyword, int minchars) ; */
  /*BADMATCH extern loadgif(char *file, int pos, int reduc) ; */
  /*BADMATCH extern ReadGIF(FILE *fd, int imageNumber, int pos, int reduc, char *file) ; */
  /*BADMATCH extern ReadColormap(FILE *fd, int number, unsigned char (*buffer)[256]) ; */
  /*BADMATCH extern DoExtension(FILE *fd, int label) ; */
  /*BADMATCH extern GetDataBlock(FILE *fd, unsigned char *buf) ; */
  /*BADMATCH extern GetCode(FILE *fd, int code_size, int flag) ; */
  /*BADMATCH extern LWZReadByte(FILE *fd, int flag, int input_code_size) ; */
  /*BADMATCH extern ReadImage(FILE *fd, int len, int height, unsigned char (*cmap)[256], int interlace, int ignore, int n, int reduc, char *file) ; */
  extern int	imx_saveonfloppy_sa(void);
  extern int	imx_saveonfloppy(void);
  extern int	compressave_serial_sa(void);
  extern int	compressave_serial(void);
  extern int Lecture_nimg_ds_nWins(int n, int norm) ;
  extern int Lecture_1img_ds_1Win(int n, int norm) ;
  extern void charge_petite_picture (int num_image, char *file, int numfirst, int n, int posdsimg, int norm) ;
  extern long fread_long (long int *t, int nitem, FILE *fp) ;
  extern long fread_short (short *t, int nitem, FILE *fp) ;
  extern long fwrite_long (long int *t, int nitem, FILE *fp) ;
  extern int read_raw (char *file, int dimx, int dimy, int pos, int taille) ;
  extern int     load_mri_ipb_2d(char *file, grphic *image, int pos) ;
  extern char* getchar_header(char *fichier, int offset, int taille) ;
  extern unsigned int getuint_header(char *fichier, int offset) ;
  extern double getdouble_header(char *fichier, int offset) ;
  extern float getfloat_header(char *fichier, int offset) ;
  extern void Open_Cmd(void) ;
  extern void nOpen_Cmd(void) ;
  extern void Open_Norm(void) ;
  extern void nOpen_Norm(void) ;
  extern void n2Open_Cmd(void) ;
  extern void Open_PCX(void) ;
  extern void Open_GIF(void) ;
  extern void Open_Raw(void) ;
  extern void Open_nRaw(void) ;
  extern void Save_Bru(void) ;
  extern void Save_Raw8(void) ;
  extern void Save_Raw32(void) ;
  extern void Save_Pcx(void) ;
  extern void Save_2d_Ipb(void) ;
  extern void Save_2d_mask(void) ;
  extern int  add_mask_Ipb_2d(char * file, int pos);
  extern void Save_Ipbr(void) ;
  extern void Save_nBru(void) ;
  extern void Save_nRaw8(void) ;
  extern void Save_nRaw32(void) ;
  extern void Save_nPcx(void) ;
  extern void Save_nIpb2d(void) ;
  extern void Refresh_Order(void) ;
  extern void Zoom_Visu(void) ;
  extern void Zoom_Visun(void) ;
  extern int scaleimg_p(grphic *imdeb, grphic *imres) ;
  extern int imx_scaleimg_p(grphic *imdeb, grphic *imres, int meth, int factx, int facty) ;
  extern char *Acrnema_getGroup(char *file, short unsigned int g1, short unsigned int g2, char *data, int format) ;
  extern long *getmri_convert(long int *buff, void *buffer, int data_size, int format, int num_items) ;
  extern char* Ecat7_getFromMatrix(char *file, short unsigned int indMatrix, short unsigned int offset, short unsigned int len, char *data, char *type) ;
  extern void read_filetxt_3d(void);
  extern int read_filetxt(char *file, grphic3d *roi);
  extern void    IMXFile_Capture_linux(int type);
#ifdef __cplusplus 
}
#endif /*__cplusplus*/ 
