/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/**************************************************************************/
/*!    
***     \file:       imx_file.c
*** 
***     Project:    Imagix 1.01
***
***     \brief Description:
***     Manages bruker files and pictures in memory
***
***     Copyright (c) 1993, ULP-IPB Strasbourg.
***     All rights are reserved.
***
***     >>  Last user action by: Mr. ETTER      on          1994
***     >>  Last user action by: Mr. ARMSPACH   on Apr 09th 1995
***     >>  Last user action by: Mr. GOUVERNEUR on Nov 24th 1998
***     >>  Last user action by: Mr. BOULEAU    on Nov 18th 2000
#include <config.h>
***         - moved #includes, #defines, type declarations, etc.
***           from imx_file.c to noyau/io/imx_file.h
***     >>  Last user action by: Mr. CHEVRIER   on Oct 25th 2001
***         - added prototypes and headers.
*** 
**************************************************************************/

#include <config.h>

#ifdef WIN32
# include <windows.h>
# include <io.h>
#endif 

#ifndef WIN32
#include <unistd.h>
#else
#ifndef S_ISDIR
#define S_ISDIR(m) (((m) & _S_IFDIR) == _S_IFDIR)
#endif
#endif

#include <limits.h>

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/io/imx_file.h"
#include "noyau/extern.h"
#include "outils/imx_misc.h"

#ifdef __GTK_GUI
#include "serie/imx_seri.h"
#include "noyau/gui/interface.h"
#include "noyau/gui/imx_picture2d.h"
#include "noyau/gui/imx_picture_common.h"
#include "noyau/io/imx__atm.h"
#include "outils/imx_tool_2d.h"
#include "noyau/gui/imx_logical_cmap.h"
#include "traitement/defor_3d.h"
#include "noyau/gui/imxmouse.h"
#include "math/operation.h"
#include "noyau/imx_zone.h"
#include "noyau/gui/imx_screenshot_gtk.h"
#endif /* __GTK_GUI */


/* -----------------------------------------------------------------*/
/* -----------------------------------------------------------------*/
/*------------------------------     Declaration pour PCX lecture  */
extern void loadpcx(char *file, int pos, int reduc);

extern FILE* p_openr(char *name);
extern void p_close(FILE *f);
extern void pcx_unpack_pixels(unsigned char *pixels, unsigned char *bitplanes, int bytesperline, int planes, int bitsperpixel);
extern void loadgif(char *file, int pos, int reduc);
extern void ReadGIF(FILE *fd, int imageNumber, int pos, int reduc, char *file);
extern int ReadColormap(FILE *fd, int number, unsigned char (*buffer)[256]);
extern int DoExtension(FILE *fd, int label);
extern int GetDataBlock(FILE *fd, unsigned char *buf);
extern int GetCode(FILE *fd, int code_size, int flag);
extern int LWZReadByte(FILE *fd, int flag, int input_code_size);
extern void ReadImage(FILE *fd, int len, int height, unsigned char (*cmap)[256], int interlace, int ignore, int n, int reduc, char *file);
extern int Lecture_nimg_ds_nWins(int n, int norm);
extern int Lecture_1img_ds_1Win(int n, int norm);
extern void charge_petite_picture (int num_image, char *file, int numfirst, int n, int posdsimg, int norm);
extern int read_raw (char *file, int dimx, int dimy, int pos, int taille);
extern int     load_mri_ipb_2d(char *file, grphic *image, int pos);
extern char* allocrow(int cols, int size);
extern char** allocarray(int cols, int rows, int size);
extern void read_pcx_image(FILE *fp, unsigned char *buf, int BytesPerLine, int Planes, int Height);
extern void pcx_planes_to_pixels(unsigned char *pixels, unsigned char *bitplanes, int bytesperline, int planes, int bitsperpixel);

// Ajout de prototypes. Remy CHEVRIER. Le 25 octobre 2001.
extern int IMXColor_BindColormap_ptr(ptr_grphic ptrPic, int nNewCMap);
extern int imx_iniparaimg_p (grphic *im1);


// Fonctions d�lar�s dans le fichier.
int	imx_nonoise_mask(grphic *mask, grphic *im, int res_pos);

typedef unsigned char gray;

#define ARGS(alist) ()

typedef gray pixval;

typedef struct 
{
  pixval r, g, b;
} pixel;

#define P_GETR(p) ((p).r)
#define P_GETG(p) ((p).g)
#define P_GETB(p) ((p).b)

#define P_ASSIGN(p,red,grn,blu) do { (p).r = (red); (p).g = (grn); (p).b = (blu); } while ( 0 )

/* Luminance macro. */

#define P_LUMIN(p) ( 0.299 * P_GETR(p) + 0.587 * P_GETG(p) + 0.114 * P_GETB(p) )

#define pcx_allocarray( cols, rows ) ((pixel**) allocarray( cols, rows, sizeof(pixel) ))



/*----------------------------------------------------------------------------*/
#define PCX_HDR_SIZE    128             /* size of PCX header           */
#define PCX_256_COLORS  0x0c            /* magic number for 256 colors  */

#define MAXCOLORS       256

/*       Declaration pour lecture fichier GIF     ---------------------------*/
#define    MAXCOLORMAPSIZE        256

#define CM_RED        0
#define CM_GREEN    1
#define CM_BLUE        2

#define    MAX_LWZ_BITS        12

#define INTERLACE        0x40
#define LOCALCOLORMAP    0x80
#define BitSet(byte, bit)    (((byte) & (bit)) == (bit))

#define    ReadOK(file,buffer,len)    (fread(buffer, len, 1, file) != 0)

#define LM_to_uint(a,b)            (((b)<<8)|(a))

struct {
  unsigned int    Width;
  unsigned int    Height;
  unsigned char    Colormap[3][MAXCOLORMAPSIZE];
  unsigned int    BitPixel;
  unsigned int    ColorResolution;
  unsigned int    Background;
  unsigned int    AspectRatio;
} GifScreen;

struct {
  int    transparent;
  int    delayTime;
  int    inputFlag;
  int    disposal;
} Gif89 = { 
  -1, -1, -1, 0 };

pixel    *Image = NULL;
int    verbose;
int    showComment;



/* ---is_mri()----------------------------------------------------------*/
/*    Is a file specified is a magnetic resonance image ? This fct      */
/*  tries to answer to this question...                                 */
/*                                                                      */
/*    Call:    int    is_mri(path,file);                                */
/*            path: where to find picture file,                         */
/*            file: picture filename.                                   */
/*    Return: 0< :Err appends...                                        */
/*        0  :Probably no header in this file...                        */
/*        1  :One header...                                             */
/*        1+2:One header and data :means it might be mri                */
/*                                                                      */
/*    NOTE: We never be sure that file specified is a picture... But,   */
/*          we check that file specified might be one ! To check if     */
/* it is a picture we first look 1st 24bits. If these bits are:         */
/*        0x4784f5: it might be a BRUKER picture...                     */
/*        0x6931  : it might be a VISILOG picture...                    */
/*        0x0x21494e54  : it might be a INTERFILE picture...            */
/*  For Bruker picture                                                  */
/* We also look if file contains enough bits to be a picture. A file    */
/* with an 256x256 picture definition must contain at least one header  */
/* of 18423bits and 256*256bits :total: 283..bits.                      */
/* If it is 128x128 picture resolution, file will be compose with 83968 */
/* (divided in two parts: header 18432, data 128*128 ).                 */
/* Finally, if it is only header file, the lenght will be 18432bits.    */
/*     So, to be sure that file might be a picture, file will contain   */
/* 0x4784f5, and be 18432 longer or x resolution * y resolution + 18432 */
/* longer.                                                              */
/*                                                                      */
/*    VISILOG Picture cannot be use with this application.              */
/*                                                                      */
/* ---------------------------------------------------------------------*/
int    
is_mri(char *path, char *file)
{
  long st_three_sixbit,t[4],size=0l;
  char cur_dir[255],str[80],msg_err[80];
  FILE *fp;
  char siemens[8];  /* for siemens */
  int repo;         /* for Siemens */
  int len;
  short endian=0;  /* for ACR-NEMA */
  char *data, *data2;
  struct stat	buf; 

  int acrnema_version;
  unsigned long lgrp;
  char a,b,c,d;
  char header[5];
  int indice_of_data=0; // indice du tableau data

  strcpy(cur_dir,path);
  if (file!=NULL)
    strcat(cur_dir,file);

  stat(cur_dir, &buf);
  /* First: Check first bits values...    */
  if((fp=fopen(cur_dir,"rb"))==NULL || S_ISDIR(buf.st_mode)) {
    if(_dialog_messages_status&0x00000002l)
      {
	char e[256];
	sprintf(e, "Open file error on >%s<", cur_dir);
	PUT_ERROR( e);
      }		
    perror(NULL);
    return(-1);
  }
  /*fread(t,sizeof(long),3,fp);*/
  /* MB: ca l'air d'etre un gros melange la dedans */
  /* pour ipb il ne faut il lire un long ou une chaine au debut?*/
  FREAD_LONG(t,3,fp);
  fseek(fp,0l,SEEK_END);
  size=(long)ftell(fp);
  fseek(fp,0l,SEEK_SET);
  fread(header,4,1,fp);
  header[4]=0;
  fclose(fp);
  st_three_sixbit=t[0];


  if((st_three_sixbit&0xffffffff)==0x21494e54) {
    return(INTERFILE_FORMAT);
  }

  if((st_three_sixbit&0xffffffff)==0x23235449) {
    return(PARAVISION_FORMAT);
  }

  if((st_three_sixbit&0xffffffff)==0x49504232) {
    return(IPB_2D_FORMAT);
  }

  if((st_three_sixbit&0xffffffff)==0x49504233 || !strcmp(header,"IPB3")) {
    return(IPB_3D_FORMAT);
  }

  if((st_three_sixbit&0xffffffff)==0x23495042) {
    return(IPB_FORMAT);
  }

  if((st_three_sixbit&0xffffffff)==0x0a050108)  {
    /* PCX_FORMAT  */
    return(PCX_FORMAT);
  }

  if((st_three_sixbit&0xffffff00)==0x4784f500) {
    /* BRUKER:    */
    /* NOTE: 4: Means number of eight bits to get one pixel.    */
    if(size==(long)(4l*256l*256l+18432l) ||
       size==(long)(4l*512l*512l+18432l) ||
       size==(long)(4l*64l*64l+18432l) ||
       size==(long)(4l*128l*128l+18432l)/* || Only when header
					   size==(long)(18432l)*/) {
      /* Bruker file... */
      return(BRUKER_FORMAT);
    }
    strcpy(msg_err,"Err, not a picture,\nOr not for this program.\n");
    PUT_MESG(msg_err);
    if(_dialog_messages_status&0x00000002l)
      { PUT_ERROR(msg_err); 
      }
  }

  if((st_three_sixbit&0xffff0000)==0x69310000) {
    /*VISILOG:    */
    sprintf(str," visilog file\n");
    if(_dialog_messages_status&0x00000002l) {
      PUT_ERROR("Visilog file !\nCannot read.");
      return(-5);    /* ERR: Program cannot manage VISILOG... */
    }
    return(VISILOG_FORMAT);
  }

  if((st_three_sixbit&0xffffffff)==0x0000015c)
    {
      /*ANALIZE format:    */
      return(ANALYZE_FORMAT_BIG_ENDIAN);
    }
	
  if((st_three_sixbit&0xffffffff)==0x5c010000) 
    {
      /*ANALIZE format:    */
      return(ANALYZE_FORMAT_LITTLE_ENDIAN);
    }

  if(size==(long)(4l*256l*256l)
     || size==(long)(4l*512l*512l)
     || size==(long)(4l*64l*64l)
     || size==(long)(4l*128l*128l)) {
    return(RAW4_FORMAT);
  }

  /* Format SMIS */
  if ( ((len = strlen(path)) >= 4)  
       && (!strcmp(path + len -4, ".sur") || !strcmp(path + len -4, ".SUR"))) {
    return(SMIS_FORMAT);
  }

  /* debut SIEMENS Magnetom Vision format :
     format tres specifique, les 1ers bits etant des dates, nous irons chercher
     le mot cle: SIEMENS , normalement situe a l'offset 96 (sur 7 caracteres).
     Implemente par cyril feruglio le 23-04-1998*/    
  strcpy(siemens, getchar_header(cur_dir, 0x60, 7));
  repo = (strcmp(siemens,"SIEMENS"));
  if(repo==0){
    return(SIEMENS_FORMAT);
  }
  /* fin Siemens Magnetom format */
	
  /* debut Siemens Vanderbilt Format Feruglio Cyril Mai 1998 */
  strcpy(siemens, getchar_header(cur_dir, 0, 5));
  repo = (strcmp(siemens,"Modal"));
  if(repo==0){
    return(VANDERBILT_FORMAT);
  }
  /* fin Siemens Vanderbilt format */

  /* Format ECAT 7 */
  if(!strcmp("MATRIX",getchar_header(cur_dir,0,6)))return ECAT7_FORMAT;
  /* fin format ECAT7 */

  /* Format Dicom */
  if(!strcmp("DICM",getchar_header(cur_dir,0x80,4))) {
    unsigned char *data=(unsigned char *)NULL;
    unsigned char *data2=(unsigned char *)NULL;
    unsigned char *data3=(unsigned char *)NULL;
	  
    /* Recherche si le fichier est Littlendian ou bigendian */
    data=(unsigned char *)Acrnema_getGroupHeader(cur_dir, 0x02, 0x10, (char *)data, DICMLEXP_FORMAT);
    /* data TransfertSyntaxUId */
    data2=(unsigned char *)Acrnema_getGroupHeader(cur_dir, 0x02, 0x02, (char *)data2, DICMLEXP_FORMAT);
    /* data2 MediaStorageSOPClassUID */
    data3=(unsigned char *)Acrnema_getGroupHeader(cur_dir, 0x0008, 0x0016, (char *)data2, DICMLEXP_FORMAT);
    /* data3 SOPClassUID */
	
    /* en fonction du fichier dont nous disposons, les champs sont remplis avec plus ou moins de precision
       on utilise le fait qu'une NMImage soit 3D.
        Deux tags nous donnent la valeur qui nous permet de l'identifier...Si l'un d'eux n'est 
       pas rempli, on utilise le second.
       Data nous donne l'Endian du fichier : Big ou Little
    */			
//    printf("valeur de		data : %s\tvaleur de data2 : %s\n", data, data2);
	  	  
    if ( (!strcmp ((const char *) data2,"1.2.840.10008.5.1.4.1.1.20")) 
	 || (data3 && !strcmp ((const char *) data3,"1.2.840.10008.5.1.4.1.1.20"))) 
      {
	if ( strcmp ((const char *) data,"1.2.840.10008.1.2.1")==0) return DICM3DLEXP_FORMAT;
	if ( strcmp ((const char *) data,"1.2.840.10008.1.2.2")==0) return DICM3DBEXP_FORMAT;
	if ( strcmp ((const char *) data,"1.2.840.10008.1.2")  ==0) return DICM3DLIMP_FORMAT;
	return DICM3DLIMP_FORMAT;
      }else{
	if ( strcmp ((const char *) data,"1.2.840.10008.1.2.1")==0) return DICMLEXP_FORMAT;
	if ( strcmp ((const char *) data,"1.2.840.10008.1.2.2")==0) return DICMBEXP_FORMAT;
	if ( strcmp ((const char *) data,"1.2.840.10008.1.2")  ==0) return DICMLIMP_FORMAT;
	return DICMLIMP_FORMAT;
      }
  }
  /* fin format dicom */

  /* Format ACR-NEMA  Burger Laurent      novembre 1999*/
  fp=fopen(cur_dir,"rb");
  if (fp==NULL) 
    {
      char e[256];
      sprintf(e, "Open file error on >%s<", file);
      PUT_ERROR( e);	
      return (0);
    }
  lgrp=0;
  a=b=c=d=0;
  fread (&d, 1, 1, fp);
  fread (&c, 1, 1, fp);
  fread (&b, 1, 1, fp);
  fread (&a, 1, 1, fp);
  lgrp=(d<<24)+(c<<16)+(b<<8)+a;
  while (!feof(fp) && lgrp!=0x08001000 && lgrp!=0x00080010)
    {	
      d=c;
      c=b;
      b=a;
      a=0;
      fread(&a, 1, 1, fp);
      lgrp=(d<<24)+(c<<16)+(b<<8)+a;
    }
  /* une petite manipe car : 
     on recherche 08 00 10 00 ou 00 08 00 10, mais rien ne nous dit que l'on n'a pas trouv�
     00 08 00 10 00, alors on verifie */
  d=c;
  c=b;
  b=a;
  a=0;
  fread(&a, 1, 1, fp);
  if (((d<<24)+(c<<16)+(b<<8)+a) == 0x08001000) lgrp=(d<<24)+(c<<16)+(b<<8)+a;
  if (!feof(fp))
    {
      if (lgrp==0x08001000 ){ endian=0;}
      else if (lgrp==0x00080010) { endian=1;}
      data=(char *)malloc(257);
      data2=data;
      fread(data, 256, 1, fp);
      fclose (fp);
	
      while ( (indice_of_data<256) && *data!='A' ) 
	{ indice_of_data++; data++; } 

      if (strstr(data, "ACR-NEMA 1.0")!=NULL) acrnema_version=1;
      else if (strstr(data, "ACR-NEMA 2.0")!=NULL) acrnema_version=2;
      else acrnema_version=-1;

      switch (acrnema_version)
	{
	case 1:
	  if (endian==0) return ACRNEMA1L_FORMAT;
	  else return ACRNEMA1B_FORMAT;
	  break;

	case 2:
	  if (endian==0) return ACRNEMA2L_FORMAT;
	  else return ACRNEMA2B_FORMAT;
	  break;
	
	}

      free (data2);
    }
  else fclose (fp); 	
  /*** fin de la d�ection du format ACR-NEMA 2.0 ***/

  strcpy(msg_err,"Err, not a picture,\nOr not for this program.\n");
  PUT_MESG(msg_err);
  if(_dialog_messages_status&0x00000002l)
    { PUT_ERROR(msg_err); 
    }

  return(0);
}

/* --getmri() ---------------------------------------------------------    
**
**    Gets a bruker picture and returns location in memory.     
**
**    Usage:    getmri(file)
**        file: char*: Represents a file witch can contains path 
**
**    Return: This fct returns a long pointer on memory
**        location of bruker file picture contain.
**                     
**                 
**    Example: void proc()
**              { long *ptr, ...
**                ...
**
**            ... ; 
**            if((ptr=getmri("/root/sub_directory/PICTURE")==NULL)
**                         PUT_WARN("Cannot read picture... ");
**                      else printf("picture loaded... ");
**                      ...
*/
long   *getmri(char *file, int numfirst)
{
  /* #Parametre: fichier - pointeur o�se trouve l'image...       */
  /* 1 Est-ce un fichier ?                                       */
  /* 2 Est-ce une image BRUKER ???                               */
  /* 3 Lecture et mise �jour d'un buffer                        */
  /*    deb correspond a l'offset lecture dans le fichier image   */

  char str[80], *file_mri;
  long *buff = NULL, size;
  FILE *fp;
  long format;
  char datafile[PATH_LEN];      /* For SIEMENS Vanderbilt format  */
  int num_items = 0, data_size = 0;
  void *buffer = NULL;
  long offset = 0;
  struct stat status;
  char v[16];
  long image, number_img, slice;
  long width, height, depth;

  short readFormat;

  errno = 0;

  if((fp = fopen(file,"rb"))==NULL) {
    _imagix_error = IMX_FILE_NOT_FOUND;        /* Open error    */
    _imagix_err_ext = errno;
    return(NULL);
  }

  format = is_mri(file, NULL);
	
  if(format != BRUKER_FORMAT
     && format != INTERFILE_FORMAT
     && format != IPB_2D_FORMAT
     && format != IPB_3D_FORMAT
     && format != RAW4_FORMAT
     && format != PARAVISION_FORMAT
     && format != PCX_FORMAT
     && format != SIEMENS_FORMAT
     && format != VANDERBILT_FORMAT
     && format != SMIS_FORMAT
     && format != ACRNEMA2L_FORMAT
     && format != ACRNEMA2B_FORMAT
     && format != ACRNEMA1L_FORMAT
     && format != ACRNEMA1B_FORMAT
     && format != ECAT7_FORMAT
     && format != DICMBIMP_FORMAT
     && format != DICMBEXP_FORMAT
     && format != DICMLIMP_FORMAT
     && format != DICMLEXP_FORMAT
     ){
    PUT_WARN( "Warning, Bad file...");
    _imagix_error = IMX_NOT_A_PICTURE_FILE;    /* Bad file    */
    _imagix_err_ext = format;
    fclose(fp);
    return(NULL);
  }

  if(format == BRUKER_FORMAT){
    /*   Affichage dans la boite de dialogue   */
    sprintf(str,">%s< bruker file\n", file);
    PUT_MESG( str);

    data_size = 4;
    fseek(fp, 0, SEEK_END);
    offset = 18432;
    num_items = (ftell(fp) - offset)/4;
  }

#ifdef __GTK_GUI
  if(format == PARAVISION_FORMAT) {
    fclose(fp);

    /*   Affichage dans la boite de dialogue   */
    sprintf(str,">%s< paravision file\n", file);
    PUT_MESG( str);
    file_mri = getfilename_2dseq(file);

    if((fp = fopen(file_mri,"rb")) == NULL) {
      _imagix_error = IMX_FILE_NOT_FOUND;        /* Open error    */
      _imagix_err_ext = errno;
      return(NULL);
    }

    data_size = getmri_dattype(file, 0, MODE_2D);
    num_items = getmri_width(file, 0, MODE_2D) * getmri_height(file, 0, MODE_2D);
    offset = (numfirst-1) * num_items * data_size;

    FREE(file_mri);
  }
#endif /* __GTK_GUI */

  if(format == INTERFILE_FORMAT) {
    fclose(fp);

    /*   Affichage dans la boite de dialogue   */
    sprintf(str,">%s< interfile file\n", file);
    PUT_MESG( str);

    /* gets the name of the data file */
    if(strlen(file) < 3) {
      char e[256];
      sprintf(e, "ERREUR: nom de fichier invalide >%s<", file);
      PUT_ERROR(e);
      return(NULL);
    }

    strncpy(datafile, file, PATH_LEN);
    strcpy(datafile+strlen(datafile)-3, "IMG");

    if((fp = fopen(datafile, "rb")) == NULL) {
      _imagix_error = IMX_FILE_NOT_FOUND;        /* Open error    */
      _imagix_err_ext = errno;
      return(NULL);
    }

    num_items = getmri_width(file, 0, MODE_2D) * getmri_height(file, 0, MODE_2D);
    data_size = getmri_dattype(file, 0, MODE_2D);   /*sizeof(char);*/
    offset = (numfirst-1) * num_items;
  }

  if (format == IPB_2D_FORMAT || format == IPB_3D_FORMAT) {
    fclose(fp);

    /* gets the name of the data file */
    if(strlen(file) < 3) {
      char e[256];
      sprintf(e, "ERREUR: nom de fichier invalide >%s<", file);
      PUT_ERROR(e);
      return(NULL);
    }

    strncpy(datafile, file, PATH_LEN);
    strcpy(datafile+strlen(datafile)-3, "img");

    /* shows picture type */
    if(format == IPB_2D_FORMAT)
      sprintf(str, ">%s< ipb 2d file\n", file);
    else
      sprintf(str, ">%s< ipb 3d file\n", file);
    PUT_MESG( str);

    /* gets data offset in file */
    if(format == IPB_2D_FORMAT) {
      num_items = getmri_width(file, numfirst, MODE_2D) * getmri_height(file, numfirst, MODE_2D);
      data_size = getmri_dattype(file, numfirst, MODE_2D);

      strcpy(v, getheader_interfile(file, "offset in raw=", ENTIER_IPB_2D, numfirst));
      offset = atol(v);
    }
    else {
      number_img = getmri_numberimg(file);
      image = 1;

      width = getmri_width(file, image, MODE_3D);
      height = getmri_height(file, image, MODE_3D);
      depth = getmri_depth(file, image, MODE_3D);
      slice = depth;
      data_size = getmri_dattype(file, image, MODE_3D);

      size = width * height * depth * data_size * number_img;
      if(stat(datafile, &status) != 0) {
	char e[256];
	sprintf(e, "[getmri] Error: cannot get file size from >%s<", datafile);
	PUT_ERROR(e);
	return(NULL);
      }
      if(size != status.st_size) {
	PUT_WNDMSG( "File size doesn't match !\nIf this is an old file please convert it !");
	return(NULL);
      }

      offset = width * height * depth * data_size;

      while(image < number_img && slice < numfirst) {
	image ++;
	width = getmri_width(file, image, MODE_3D);
	height = getmri_height(file, image, MODE_3D);
	depth = getmri_depth(file, image, MODE_3D);
	slice += depth;
	data_size = getmri_dattype(file, image, MODE_3D);
	offset += width * height * depth * data_size;
      }

      if(slice < numfirst) {
	/* erreur */
	sprintf(str, "[getmri] only %ld slices in file \"%s\"\n", slice, file);
	PUT_MESG( str);
	_imagix_error = IMX_PICTURE_NOT_FOUND;
	_imagix_err_ext = 0;
	return(NULL);
      }

      offset -= (slice - (numfirst - 1)) * width * height * data_size;
      num_items = width * height;
    }

    if((fp = fopen(datafile, "rb")) == NULL) {
      _imagix_error = IMX_FILE_NOT_FOUND;        /* open error    */
      _imagix_err_ext = errno;
      return(NULL);
    }
  }

  if(format == RAW4_FORMAT) {
    /*   Affichage dans la boite de dialogue   */
    sprintf(str,">%s< raw file\n", file);
    PUT_MESG( str);

    data_size = sizeof(long);
    fseek(fp, 0, SEEK_END);
    num_items = ftell(fp)/4;
    offset = 0;
  }

  /* ----Format Siemens Magnetom Vision, Feruglio Cyril Mai 1998------ */
  if(format == SIEMENS_FORMAT){
    /*   Affichage dans la boite de dialogue   */
    sprintf(str,">%s< Siemens file\n", file);
    PUT_MESG( str);

    /*    Calcul de la taille...    */
    fseek(fp, 0, SEEK_END);

    offset = 6144;
    data_size = sizeof(short);
    num_items = (ftell(fp) - 6144)/2;
  }

  /* ----Siemens Vanderbilt Format, Feruglio Cyril Mai 1998------ */
  if(format == VANDERBILT_FORMAT){
    fclose(fp);

    /* gets the name of the data file */
    if(strlen(file) < 12) {
      char e[256];
      sprintf(e, "ERREUR: nom de fichier invalide >%s<", file);
      PUT_ERROR(e);
      return(NULL);
    }

    strncpy(datafile, file, PATH_LEN);
    strcpy(datafile+strlen(datafile)-12, "image.bin");
    if(MESG_DEBUG) printf("file:%s", datafile);

    if((fp = fopen(datafile, "rb")) == NULL) {
      _imagix_error = IMX_FILE_NOT_FOUND;        /* Open error    */
      _imagix_err_ext = errno;
      return(NULL);
    }
    /*   Affichage dans la boite de dialogue   */
    sprintf(str, ">%s< Vanderbilt Siemens file\n", datafile);
    PUT_MESG( str);

    /*    Calcul de la taille...    */
    fseek(fp, 0, SEEK_END);
    data_size = sizeof(short);
    num_items = (ftell(fp) - getmri_depth(file, 0, MODE_2D)) / 2;
    offset = 0;
  }

  if(format == SMIS_FORMAT) {
    data_size = sizeof(short);
    num_items = getmri_width(file, 0, MODE_2D) * getmri_height(file, 0, MODE_2D);
    offset = 512;
  }

  /* for ECAT7 format */
  if(format == ECAT7_FORMAT){
    num_items=getmri_width(file, 1, MODE_2D)*getmri_height(file, 1, MODE_2D);
    data_size=getmri_dattype(file,0,MODE_2D);
    offset=0x600+num_items*data_size*(numfirst-1);
  }

  /* For ACR-NEMA 2.0 and ACR-NEMA 1.0 and DICOM Format */
  if (format==ACRNEMA2L_FORMAT || format==ACRNEMA2B_FORMAT ||
      format==ACRNEMA1L_FORMAT || format==ACRNEMA1B_FORMAT ||
      format==DICMBIMP_FORMAT  || format==DICMBEXP_FORMAT  ||
      format==DICMLIMP_FORMAT  || format==DICMLEXP_FORMAT)
    {
      unsigned int x,y;
      unsigned int sizeTotale;
      unsigned char *data=(unsigned char *)NULL;
      /* Le fichier est rest�ouvert pour les op�ations restantes qui ne sont
	 pas applicables pour le format ACR-NEMA ,donc on le ferme */
      fclose (fp);
      /* lecture des data images pour fichier acr-nema2 et DICOM*/
      if (  format==ACRNEMA2L_FORMAT || format==ACRNEMA2B_FORMAT
	    ||format==DICMLIMP_FORMAT  || format==DICMBIMP_FORMAT
	    ||format==DICMLEXP_FORMAT  || format==DICMBEXP_FORMAT ){
	buffer=Acrnema_getGroup(file, 0x7fe0, 0x10, (char *)buffer,format);
	/* lecture de la taille d'un pixel */
	data=(unsigned char *)Acrnema_getGroup(file, 0x28, 0x100, (char *)data, format);
	data_size=0;
	memcpy(&data_size, data, 2);
	free (data);
	data_size=data_size/8;
	if (_is_bigEndian) data_size=longEndianConversion(data_size);
	if(data_size != sizeof(long)) {
	  x=getmri_width(file, 1, MODE_2D);
	  y=getmri_height(file, 1, MODE_2D);
	  num_items=x*y;
	  if((buff = CALLOC(num_items, long)) == NULL) {
	    PUT_ERROR("[getmri] memory allocation error (1)");
	    FREE(buffer);
	    return(NULL);
	  }
	}
      }
      else /* lecture des data images pour fichier acr-nema1 */
	{
	  /* Calcul taille image */
	  x=getmri_width(file, 1, MODE_2D);
	  y=getmri_height(file, 1, MODE_2D);
	  /* lecture de la taille d'un pixel */
	  data=(unsigned char *)Acrnema_getGroup(file, 0x28, 0x100, (char *)data, format);
	  data_size=0;
	  memcpy(&data_size, data, 2);
	  free (data);
	  data_size=data_size/8;
	  if (_is_bigEndian) data_size=longEndianConversion(data_size);                                        
	  /*  taille totale*/
	  sizeTotale=x*y*data_size;	
	  num_items=x*y;
	  fp=fopen(file, "rb");
	  buffer=(unsigned char *)malloc(sizeTotale);
	  fseek(fp, -sizeTotale, SEEK_END);
	  fread(buffer, sizeTotale, 1, fp);
	  fclose (fp);				
	  /*buffer=Acrnema1_getRaw(file, buffer, format, sizeTotale);*/
	  if(data_size != sizeof(long)) 
	    {
	      if((buff = CALLOC(num_items, long)) == NULL) 
		{
		  PUT_ERROR( "[getmri] memory allocation error (1)");
		  FREE(buffer);
		  return(NULL);
		}
	    }
	}
    }
  /* lecture des data image pas ACRNEMA et DICOM */

  else 
    {
      if((buffer = calloc(num_items, data_size)) == NULL) {
	PUT_ERROR( "[getmri] memory allocation error (1)");
	fclose(fp);
	return(NULL);
      }

      if(data_size != sizeof(long)) {
	if((buff = CALLOC(num_items, long)) == NULL) {
	  PUT_ERROR( "[getmri] memory allocation error (1)");
	  FREE(buffer);
	  fclose(fp);
	  return(NULL);
	}
      }

      fseek(fp, offset, SEEK_SET);
      if(fread(buffer, data_size, num_items, fp) != num_items) {
	char e[256];
	sprintf(e, "[getmri] unable to read picture data from file \"%s\"\n", file);
	PUT_ERROR(e);
        FREE(buffer);
        if(data_size != sizeof(long))
	  FREE(buff);
	fclose(fp);
        return(NULL);
      }
      fclose(fp);
    }

  readFormat=getmri_readformat(file);
  buff = getmri_convert(buff, buffer, data_size, format, num_items);
  return(buff);
}

/* --getmri_index3D() ---------------------------------------------------
**
**    Returns a 3D mri index from a 2D mri index
**
**    Usage:    getmri_index3D(file, index2D)
**        file: char*: file name with path
**        index2D: int: number of the slice (or 2D mri number)
**
**    Return: This function returns a 3D mri number
**                     
*/
long    getmri_index3D(char *file, int index2D)
{
  char msg[64];
  long image, number_img, slice;

  if(is_mri(file, NULL) != IPB_3D_FORMAT) {
    return(index2D);
  }

  else {
    number_img = getmri_numberimg(file);

    image = 1;
    slice = getmri_depth(file, image, MODE_3D);

    while(image < number_img && slice < index2D) {
      image ++;
      slice += getmri_depth(file, image, MODE_3D);
    }

    if(slice < index2D) {
      /* error */
      sprintf(msg, "[getmri_index3D] only %ld slices in file \"%s\"\n", slice, file);
      PUT_MESG( msg);
      _imagix_error = IMX_PICTURE_NOT_FOUND;
      _imagix_err_ext = 0;
      return(0);
    }

    return(image);
  }
}

/* --getmri_width() ---------------------------------------------------
**
**    Gets width for a bruker picture     // Assumed pixel is 32 bits
**              and interfile format
**
**    Usage:    getmri_width(file,number,number)
**        file: char*: Represents a file witch can contains path 
**        number: int: number of the picture in the file
**        mode: int: mode to read the pictures in the file
**
**    Return: This fct returns a long witch indicates width
**                     
*/
long    getmri_width(char *file, int number, int mode)
{
  long size,width;
  double x;
  FILE *fp;
  char v[16],*file_d3proc;
  int Version, BitsPerPixel;
  int  Xmin, Ymin, Xmax, Ymax;
  char mess[80];
  long format;
  unsigned long value;
	
  format = is_mri(file,NULL);

  /*    For Brucker file format   */
  if(format == BRUKER_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0l);
      }

      /* Calcul de la taille...    */
      fseek(fp,0l,SEEK_END);
      size=ftell(fp);

      size-=HEADER_SIZE;
      size/=4;            /* Cause 32bits par points...    */

      fclose(fp);

      switch(size) {
      case 128l*128l : /* Image 128*128... */
	return((long)128);
      case 256l*256l : /* Image 256*256... */
	return((long)256);
      default: /* Unknow size */
	x=(double)size;
	width=(long)sqrt(x);
	if((long)width*(long)width==size)
	  return((long)width);
	else
	  return(0l);
      }
    }

#ifdef __GTK_GUI
  /*  For PARAVISION format   */
  if(format == PARAVISION_FORMAT)
    {
      file_d3proc=getfilename_d3proc(file);
      strcpy(v, getheader_interfile(file_d3proc, "##$IM_SIX=", ENTIER, 0));
      FREE(file_d3proc);
      return(atol(v));
    }
#endif /* __GTK_GUI */

  /*  For Interfile format   */
  if(format == INTERFILE_FORMAT)
    {
      strcpy(v, getheader_interfile(file, "!matrix size [1]:=", ENTIER, 0));
      return(atol(v));
    }

  /*    For raw4 file format   */
  if(format == RAW4_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	fclose(fp);
	return(0l);
      }

      /* Calcul de la taille...    */
      fseek(fp,0l,SEEK_END);
      size=ftell(fp);

      size/=4;            /* Cause 32bits par points...    */

      fclose(fp);

      switch(size) {
      case 64l*64l : /* Image 64*64... */
	return((long)64);
      case 128l*128l : /* Image 128*128... */
	return((long)128);
      case 256l*256l : /* Image 256*256... */
	return((long)256);
      case 512l*512l : /* Image 512*512... */
	return((long)512);
      default: /* Unknow size */
	x=(double)size;
	width=(long)sqrt(x);
	if((long)width*(long)width==size)
	  return((long)width);
	else
	  return(0l);
      }
    }

  /*   For Siemens Magnetom Vision format, Cyril Feruglio mai 1998*/
  if(format == SIEMENS_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0l);
      }

      /* Calcul de la taille ...    */
      fseek(fp,0l,SEEK_END);
      size = ftell(fp);
      size -= 6144;           /* Car les donnees commence au bloc 6144 */
      size /= 2;            /* Car 16 bits par points...    */

      fclose(fp);

      switch(size) {
      case 64l*64l : /* Image 64*64... */
	return((long)64);
      case 128l*128l : /* Image 128*128... */
	return((long)128);
      case 256l*256l : /* Image 256*256... */
	return((long)256);
      case 512l*512l : /* Image 512*512... */
	return((long)512);
      default: /* taille inconnue */
	x=(double)size;
	width=(long)sqrt(x);
	if((long)width*(long)width==size)
	  return((long)width);
	else
	  return(0);
      }
    }

  /* For SIEMENS Vanderbilt FORMAT Feruglio Cyril Mai 1998 */
  if(format == VANDERBILT_FORMAT)
    {
      strcpy(v, getheader_interfile(file, "Rows := ", ENTIER, 0));
      printf("width= %s",v);
      return(atol(v));
    }

  /*   For IPB 2D file format */
  if(format == IPB_2D_FORMAT)
    {
      /* Gets only first parameter in 3D mode */
      if(mode == MODE_3D)
	number = 0;

      strcpy(v, getheader_interfile(file, "width=", ENTIER_IPB_2D, number));
      return(atol(v));
    }

  /*   For IPB 3D file format */
  if(format == IPB_3D_FORMAT)
    {
      /* Gets an image number corresponding to the slice given by 'number' */
      if(mode == MODE_2D)
	number = getmri_index3D(file, number);

      strcpy(v, getheader_interfile(file, "width=", ENTIER_IPB_2D, number));
      return(atol(v));
    }

  /*    For pcx file format   */
    
  if(format == PCX_FORMAT)
    {
      fp = p_openr(file);
      //
      // read the PCX header
      //
      if (GetByte(fp) != PCX_MAGIC)
	return(0);

      Version     = GetByte(fp);  // get version #

      if (GetByte(fp) != 1)       // check for PCX run length encoding
        {
	  sprintf(mess,"%s has unknown encoding scheme \n", file );
	  PUT_MESG(mess);
	  return(0);
        }
      BitsPerPixel= GetByte(fp); printf( "Bpp = %d\n", BitsPerPixel);
      Xmin        = GetWord(fp); printf( "Xmin = %d\n", Xmin);
      Ymin        = GetWord(fp); printf( "Ymin = %d\n", Ymin);
      Xmax        = GetWord(fp); printf( "Xmax = %d\n", Xmax);
      Ymax        = GetWord(fp); printf( "Ymax = %d\n", Ymax);
      width       = (Xmax - Xmin) + 1;
      //Height      = (Ymax - Ymin) + 1;
      fclose(fp);
      return((long)width);
    }

  /* For smis format */
  if(format == SMIS_FORMAT) {
    if((fp = fopen(file, "rb")) == NULL) {
      return(0);
    }

    value = GetDword(fp);

    fclose(fp);
    if (_is_bigEndian) return (longEndianConversion(value));
    return((unsigned long) value);
  }
	
	
  /* for ACR-NEMA Format    Burger Laurent  novembre 1999 */
  if (format==ACRNEMA2B_FORMAT ||format==ACRNEMA2L_FORMAT || format==ACRNEMA1B_FORMAT ||format==ACRNEMA1L_FORMAT) 
    {
	
      char *data=(char *)NULL;
      unsigned short width=0;
		

      data=Acrnema_getGroup(file, 0x28, 0x10, data,format);
      if (data!=(char *)NULL)	memcpy(&width, data, 2);
      if (((format==ACRNEMA2B_FORMAT || format==ACRNEMA1B_FORMAT) && !_is_bigEndian) 
	  || ((format==ACRNEMA2L_FORMAT || format==ACRNEMA1L_FORMAT) && _is_bigEndian)) 
	width=shEndianConversion(width);
      free (data);
      return (width);
    }

  /* for dicom 2D format */
  if(format == DICMBIMP_FORMAT || format == DICMLIMP_FORMAT ||
     format == DICMBEXP_FORMAT || format == DICMLEXP_FORMAT)
    {
      char *data=(char *)NULL;
      unsigned short width=0;
		

      data=Acrnema_getGroup(file, 0x28, 0x10, data,format);
      if (data!=(char *)NULL)	memcpy(&width, data, 2);
      if ( _is_bigEndian) /* On suppose endian dans is_mri =0 */ 
	width=shEndianConversion(width);
      free (data);
      return (width);
    }
	
  /* for dicom 3D format */
  if(format == DICM3DBIMP_FORMAT || format == DICM3DLIMP_FORMAT ||
     format == DICM3DBEXP_FORMAT || format == DICM3DLEXP_FORMAT)
    {
      char *data=(char *)NULL;
      unsigned short width=0;
		

      data=Acrnema_getGroup(file, 0x28, 0x10, data,format);
      if (data!=(char *)NULL)	memcpy(&width, data, 2);
      if ( _is_bigEndian) 
	width=shEndianConversion(width);
      free (data);
      return (width);
    }
	
  /* for ECAT7 format */
  if (format == ECAT7_FORMAT){
    char *data=(char*)malloc(3);
    long ret;
    Ecat7_getFromMatrix(file,number,0x4,2,data,"int");
    ret=atoi(data);
    FREE(data);
    return ret;
  }

  return(0l);
}

/* --getmri_height() --------------------------------------------------
**
**    Gets height for a bruker picture    // Assumed pixel is 32 bits
**               and for interfile format  
**
**    Usage:    getmri_height(file, number,number)
**        file: char*: Represents a file witch can contains path 
**        number: int: number of the picture in the file
**        mode: int: mode to read the pictures in the file
**
**    Return: This fct returns a long witch indicates height
**                     
*/
long    getmri_height(char *file, int number, int mode)
{
  long size,height;
  double x;
  FILE *fp;
  char v[16],*file_d3proc;
  int  Version;
  int  Xmin, Ymin, Xmax, Ymax;
  int  BitsPerPixel;
  char mess[80];
  long format;
  unsigned long value;

  format = is_mri(file,NULL);

  if(format == BRUKER_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0l);
      }

      /* Calcul de la taille...    */
      fseek(fp,0l,SEEK_END);
      size=ftell(fp);

      size-=HEADER_SIZE;
      size/=4;            /* Cause 32bits par points...    */

      fclose(fp);

      switch(size) {
      case (long)128l*128l : /* Image 128*128... */
	return((long)128);
      case (long)256l*256l : /* Image 256*256... */
	return((long)256);
      default: /* Unknow size */
	x=(double)size;
	height=(long)sqrt(x);
	if((long)height*(long)height==size)
	  return((long)height);
	else
	  return(0);
      }
    }

#ifdef __GTK_GUI
  /*  For PARAVISION format   */
  if(format == PARAVISION_FORMAT)
    {
      file_d3proc=getfilename_d3proc(file);
      strcpy(v, getheader_interfile(file_d3proc, "##$IM_SIY=", ENTIER, 0));
      FREE(file_d3proc);
      return(atol(v));
    }
#endif /* __GTK_GUI */

  /*  For Interfile format   */
  if(format == INTERFILE_FORMAT)
    {
      strcpy(v, getheader_interfile(file, "!matrix size [2]:=", ENTIER, 0));
      return(atol(v));
    }

  /*    For Siemens Magnetom Vision format, Cyril Feruglio mai 1998  */
  if(format == SIEMENS_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0l);
      }

      /*  Calcul de la taille ... */

      fseek(fp,0l,SEEK_END);
      size = ftell(fp);
      size -= 6144;             /* Car les donnees commence au bloc 6144 */
      size /= 2;        /* Car 16 bits par points...    */

      fclose(fp);

      switch(size) {
      case 64l*64l : /* Image 64*64... */
	return((long)64);
      case 128l*128l : /* Image 128*128... */
	return((long)128);
      case 256l*256l : /* Image 256*256... */
	return((long)256);
      case 512l*512l : /* Image 512*512... */
	return((long)512);
      default: /* taille inconnue */
	x=(double)size;
	height=(long)sqrt(x);
	if((long)height*(long)height==size)
	  return((long)height);
	else
	  return(0l);
      }
    }

  /* For SIEMENS Vanderbilt FORMAT Feruglio Cyril Mai 1998 */
  if(format == VANDERBILT_FORMAT)
    {
      strcpy(v, getheader_interfile(file, "Columns := ", ENTIER, 0));
      printf("height= %s",v);
      return(atol(v));
    }

  /*    For raw4 file format   */
  if(format == RAW4_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0l);
      }

      /* Calcul de la taille...    */
      fseek(fp,0l,SEEK_END);
      size=ftell(fp);

      size/=4;            /* Cause 32bits par points...    */

      fclose(fp);

      switch(size) {
      case 64l*64l : /* Image 64*64... */
	return((long)64);
      case 128l*128l : /* Image 128*128... */
	return((long)128);
      case 256l*256l : /* Image 256*256... */
	return((long)256);
      case 512l*512l : /* Image 512*512... */
	return((long)512);
      default: /* Unknow size */
	x=(double)size;
	height=(long)sqrt(x);
	if((long)height*(long)height==size)
	  return((long)height);
	else
	  return(0l);
      }
    }

  /*   For IPB 2D file format */
  if(format == IPB_2D_FORMAT)
    {
      /* Gets only first parameter in 3D mode */
      if(mode == MODE_3D)
	number = 0;

      strcpy(v, getheader_interfile(file, "height=", ENTIER_IPB_2D, number));
      return(atol(v));
    }

  /*   For IPB 3D file format */
  if(format == IPB_3D_FORMAT)
    {
      /* Gets an image number corresponding to the slice given by 'number' */
      if(mode == MODE_2D)
	number = getmri_index3D(file, number);

      strcpy(v, getheader_interfile(file, "height=", ENTIER_IPB_2D, number));
      return(atol(v));
    }

  /*    For pcx file format   */
    
  if(format == PCX_FORMAT)
    {
      fp = p_openr(file);
      //
      // read the PCX header
      //
      if (GetByte(fp) != PCX_MAGIC)
        {
	  return(0);
        }
      Version     = GetByte(fp);  // get version #

      if (GetByte(fp) != 1)       // check for PCX run length encoding
        {
	  sprintf(mess,"%s has unknown encoding scheme \n", file );
	  PUT_MESG(mess);
	  return(0);
        }
      BitsPerPixel= GetByte(fp);
      Xmin        = GetWord(fp);
      Ymin        = GetWord(fp);
      Xmax        = GetWord(fp);
      Ymax        = GetWord(fp);
      height      = (Ymax - Ymin) + 1;
      fclose(fp);
      return((long)height);
    }
    

  /* For smis format */
  if(format == SMIS_FORMAT) {
    if((fp = fopen(file, "rb")) == NULL) {
      return(0);
    }

    value = GetDword(fp);

    fclose(fp);
    if (_is_bigEndian) return (longEndianConversion(value));
    return((unsigned long) value);
  }
	
  /* for Acr-NEMA format  Burger Laurent novembre 1999 */
  if (format==ACRNEMA2L_FORMAT || format==ACRNEMA2B_FORMAT || format==ACRNEMA1B_FORMAT || format==ACRNEMA1L_FORMAT)
    {
      void *data=(char *)NULL;
      unsigned short height=0;
		

      data=Acrnema_getGroup(file, 0x28, 0x11, (char *)data,format);
      if (data!=(char *)NULL) memcpy(&height, data, 2);
      if (((format==ACRNEMA2B_FORMAT || format==ACRNEMA1B_FORMAT) && !_is_bigEndian) 
	  || ((format==ACRNEMA2L_FORMAT || format==ACRNEMA1L_FORMAT) && _is_bigEndian)) 
	height=shEndianConversion(height);	
      free (data);
      return (height);
    }

  /* for dicom 2D format */
  if(format == DICMBIMP_FORMAT || format == DICMLIMP_FORMAT ||
     format == DICMBEXP_FORMAT || format == DICMLEXP_FORMAT)
    {
      void *data=(char *)NULL;
      unsigned short height=0;


      data=Acrnema_getGroup(file, 0x28, 0x11, (char *)data,format);
      if (data!=(char *)NULL) memcpy(&height, data, 2);
      if ( _is_bigEndian) 
	height=shEndianConversion(height);	
      free (data);
      return (height);
    }
	
  /* for dicom 3D format */
  if(format == DICM3DBIMP_FORMAT || format == DICM3DLIMP_FORMAT ||
     format == DICM3DBEXP_FORMAT || format == DICM3DLEXP_FORMAT)
    {
      void *data=(char *)NULL;
      unsigned short height=0;


      data=Acrnema_getGroup(file, 0x28, 0x11, (char *)data,format);
      if (data!=(char *)NULL) memcpy(&height, data, 2);
      if ( _is_bigEndian)  
	height=shEndianConversion(height);	
      free (data);
      return (height);
    }
  /* for ECAT7 format */
  if (format == ECAT7_FORMAT){
    char *data=(char*)malloc(3);
    long ret;
    Ecat7_getFromMatrix(file,number,0x6,2,data,"int");
    ret=atoi(data);
    FREE(data);
    return ret;
  }

  return(0);
}

/* --getmri_depth() ---------------------------------------------------
**
**    Gets depth for a bruker picture     // Assumed pixel is 32 bits
**              anf interfile format
**
**    Usage:    getmri_depth(file,number,number)
**        file: char*: Represents a file witch can contains path 
**        number: int: number of the picture in the file
**        mode: int: mode to read the pictures in the file
**
**    Return: This fct returns a long witch indicates width
**                     
*/
long    getmri_depth(char *file, int number, int mode)
{
  long width;
  FILE *fp;
  char v[80],*file_acqp,*file_reco, *file_mri;
  int  Version;
  int  Xmin, Ymin, Xmax, Ymax;
  int  BitsPerPixel;
  char mess[80];
  int  dim_acq;
  long format;

  format = is_mri(file,NULL);

  /*    For Brucker file format   */
  if(format == BRUKER_FORMAT)
    {
      return(1);
    }

#ifdef __GTK_GUI
  /*  For PARAVISION format  */ 
  if(format == PARAVISION_FORMAT)
    {
      file_acqp=getfilename_acqp(file);
      strcpy(v, getheader_interfile(file_acqp, "##$ACQ_dim=", ENTIER, 0));
      dim_acq=atol(v);
      FREE(file_acqp);

      switch (dim_acq) {
      case 2:   

	file_mri = getfilename_imnd(file);
	strcpy(v, getheader_interfile(file_mri,
				      "##$IMND_slice_scheme=", ASCII8B, 0));
	if (strcmp(v, "Angio"))
	  { 
	    strcpy(v, getheader_interfile(file_mri,
					  "##$IMND_n_slices=", ENTIER, 0));
	  }
	else
	  { 
	      
	    strcpy(v, getheader_interfile(file_mri,
					  "##$IMND_slicepack_n_slices=(", TAB_PARA_ENTIER, 0));
	  }
	FREE(file_mri);
	return atol(v);
      case 3: 
	file_reco=getfilename_reco(file);
	strcpy(v, getheader_interfile(file_reco,
				      "##$RECO_size=(", TAB_PARA_ENTIER, 2));
	FREE(file_reco);
	return(atol(v));
      default:
	return(0l);
      }
    }


#endif /* __GTK_GUI */
  /*  For Interfile format   */
  if(format == INTERFILE_FORMAT)
    {
      strcpy(v, getheader_interfile(file, "!total number of images:=", ENTIER, 0));
      return(atol(v));
    }

  /*    For raw4 file format   */
  if(format == RAW4_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0l);
      }

      fclose(fp);
      return(1l);
    }

  /*   For Siemens Magnetom Vision format, Cyril Feruglio Mai 1998  */
  if(format == SIEMENS_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0l);
      }

      fclose(fp);
      return(1l);
    }

  /*   For IPB 2D file format */
  if(format == IPB_2D_FORMAT)
    {
      if(mode == MODE_3D)
	return(getmri_numberimg(file));

      /* there's no depth for 2D IPB pictures */
      return(0);
    }

  /*   For IPB 3D file format */
  if(format == IPB_3D_FORMAT)
    {
      if(mode == MODE_2D)
	return(0);

      strcpy(v, getheader_interfile(file, "depth=", ENTIER_IPB_2D, number));
      return(atol(v));
    }

  /*    For pcx file format   */
    
  if(format == PCX_FORMAT)
    {
      fp = p_openr(file);
      //
      // read the PCX header
      //
      if (GetByte(fp) != PCX_MAGIC)
	return(0);

      Version     = GetByte(fp);  // get version #

      if (GetByte(fp) != 1)       // check for PCX run length encoding
        {
	  sprintf(mess,"%s has unknown encoding scheme \n", file );
	  PUT_MESG(mess);
	  return(0);
        }
      BitsPerPixel= GetByte(fp);
      Xmin        = GetWord(fp);
      Ymin        = GetWord(fp);
      Xmax        = GetWord(fp);
      Ymax        = GetWord(fp);
      width       = (Xmax - Xmin) + 1;
      //Height      = (Ymax - Ymin) + 1;
      fclose(fp);
      return(1);
    }
    

  /* For SIEMENS Vanderbilt FORMAT Feruglio Cyril Mai 1998 */
  if(format == VANDERBILT_FORMAT)
    {
      strcpy(v, getheader_interfile(file, "Slices := ", ENTIER, 0));
      printf("depth= %s",v);
      return(atol(v));
    }

  if(format == SMIS_FORMAT)
    {
      char tmp[16];
      int depth;
      sprintf(tmp,":IM_TOTAL_SLICES");
      strcpy(v, getheader_interfile(file, tmp, ENTIER, 0));
      depth = atoi(v);
      return(depth);
    }
	
	
  if (format == ACRNEMA2L_FORMAT || format == ACRNEMA2B_FORMAT 
      || format==ACRNEMA1B_FORMAT || format==ACRNEMA1L_FORMAT) 
    return 0;

  /* Pour format dicom 2D */
  if(format == DICMBIMP_FORMAT || format == DICMLIMP_FORMAT ||
     format == DICMBEXP_FORMAT || format == DICMLEXP_FORMAT)
    {   
      char *nbrslice=NULL;
      long depth;

      nbrslice=Acrnema_getGroup(file, 0x28, 0x08, (char *)nbrslice, format);
	
      if ( strlen(nbrslice)==0 ){
	return 1;
      }else{
	depth=atol(nbrslice);	     
	free (nbrslice);
	return (depth);
      }
    }
    
  /* Pour format dicom 3D */
  if(format == DICM3DBIMP_FORMAT || format == DICM3DLIMP_FORMAT ||
     format == DICM3DBEXP_FORMAT || format == DICM3DLEXP_FORMAT)
    {   
      char *nbrslice=NULL;
      long depth;

      /* NumberOfFrames : depth donnne le nombre d'images dans le fichier 
	 ATTENTION : On n'utilise pas le tag NumberInAcquisition */
      nbrslice=Acrnema_getGroup(file, 0x28, 0x08, (char *)nbrslice, format);
	
      if ( strlen(nbrslice)==0 ){
	return 0;
      }else{
	depth=atol(nbrslice);	     
	free (nbrslice);
	return (depth);
      }
    }
	 
  /* for ECAT7 format */
  if (format == ECAT7_FORMAT){
    char *data=(char*)malloc(3);
    long ret;
    Ecat7_getFromMatrix(file,number,0x8,2,data,"int");
    ret=atoi(data);
    FREE(data);
    return ret;
  }

  return(0);
}

/* --getmri_numberimg() ---------------------------------------------------
**
**    Gets number image in the file    
**
**    Usage:    getmri_numberimg(file)
**        file: char*: Represents a file witch can contains path 
**
**    Return: This fct returns a long witch indicates width
**                     
*/
long    getmri_numberimg(char *file)
{
  long width;
  FILE *fp;
  char v[16],*file_d3proc;
  int  Version;
  int  Xmin, Ymin, Xmax, Ymax;
  int  BitsPerPixel;
  char mess[80];
  long format;

  format = is_mri(file,NULL);

  /*    For Brucker file format   */
  if(format == BRUKER_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0l);
      }

      fclose(fp);
      return(1l);
    }

#ifdef __GTK_GUI
  /*  For PARAVISION format   */
  if(format == PARAVISION_FORMAT)
    {
      file_d3proc=getfilename_d3proc(file);
      strcpy(v, getheader_interfile(file_d3proc, "##$IM_SIZ=", ENTIER, 0));
      FREE(file_d3proc);
      return(atol(v));
    }
#endif /* __GTK_GUI */

  /*  For Interfile format   */
  if(format == INTERFILE_FORMAT)
    {
      strcpy(v, getheader_interfile(file, "!total number of images:=", ENTIER, 0));
      return(atol(v));
    }

  /*    For raw4 file format   */
  if(format == RAW4_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0l);
      }

      fclose(fp);
      return(1l);
    }

  /*   For Siemens Magnetom Vision format, Cyril Feruglio Mai 1998  */
  if(format == SIEMENS_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0l);
      }

      fclose(fp);
      return(1l);
    }

  /*   For Siemens Vanderbilt format, Cyril Feruglio Juin 1998  */
  if(format == VANDERBILT_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0l);
      }

      fclose(fp);
      return(1l);
    }

  /*   For IPB 2D & 3D file formats*/
  if(format == IPB_2D_FORMAT || format == IPB_3D_FORMAT)
    {
      strcpy(v, getheader_interfile(file, "number of images=", ENTIER, 0));
      return(atol(v));
    }

  /*    For pcx file format   */
    
  if(format == PCX_FORMAT)
    {
      fp = p_openr(file);
      //
      // read the PCX header
      //
      if (GetByte(fp) != PCX_MAGIC)
	return(0);

      Version     = GetByte(fp);  // get version #

      if (GetByte(fp) != 1)       // check for PCX run length encoding
        {
	  sprintf(mess,"%s has unknown encoding scheme \n", file );
	  PUT_MESG(mess);
	  return(0);
        }
      BitsPerPixel= GetByte(fp);
      Xmin        = GetWord(fp);
      Ymin        = GetWord(fp);
      Xmax        = GetWord(fp);
      Ymax        = GetWord(fp);
      width       = (Xmax - Xmin) + 1;
      // Height      = (Ymax - Ymin) + 1;
      fclose(fp);
      return((long)1);
    }
    

  if(format == SMIS_FORMAT) 
    return(1);
	
  if (format == ACRNEMA2L_FORMAT || format == ACRNEMA2B_FORMAT 
      || format==ACRNEMA1B_FORMAT || format==ACRNEMA1L_FORMAT) 
    return 1;

  if(format == DICMBIMP_FORMAT || format == DICMLIMP_FORMAT ||
     format == DICMBEXP_FORMAT || format == DICMLEXP_FORMAT)
    return 1;

  if(format == DICM3DBIMP_FORMAT || format == DICM3DLIMP_FORMAT ||
     format == DICM3DBEXP_FORMAT || format == DICM3DLEXP_FORMAT)
    return 1;

  if(format==ECAT7_FORMAT){
    int count=0,tmp,indRecord=1;
    char *data=(char*)malloc(10);
    int offset=0x200;
    char temp;

    fp=fopen(file,"rb");
    fseek(fp,offset+0xC,SEEK_SET);
    fread(data,2,1,fp);
    if(!_is_bigEndian){
      temp=data[0];
      data[0]=data[1];
      data[1]=temp;
    }
    tmp=*(unsigned short*)data;
    count=tmp;
    while(tmp==31){
      count+=31;
      /* nombre de records occupes */
      fseek(fp,0x200*indRecord+0x1FC,SEEK_SET);
      fread(data,2,1,fp);
      if(!_is_bigEndian){
	temp=data[0];
	data[0]=data[1];
	data[1]=temp;
      }
      tmp=*(unsigned short*)data;
      /* indice de la prochaine table de records */
      fseek(fp,0x200*indRecord+0x1F4,SEEK_CUR);
      fread(data,2,1,fp);
      if(!_is_bigEndian){
	temp=data[0];
	data[0]=data[1];
	data[1]=temp;
      }
      if(*(unsigned short*)data)indRecord=*(unsigned short*)data;
    }
    fseek(fp,0x200*indRecord+0xC,SEEK_SET);
    fread(data,2,1,fp);
    if(!_is_bigEndian){
      temp=data[0];
      data[0]=data[1];
      data[1]=temp;
    }
    count+=*(unsigned short*)data;
    FREE(data);
    return count;
  }


  return(0);
}

/* --getmri_minpixel() ------------------------------------------------    
**
**    Gets minpixel for a bruker picture
**
**    Usage:    getmri_minpixel(bruker_file,number,number)
**        file: char*: Represents a file. (May contains path )
**        number: int: number of the picture in the file
**        mode: int: mode to read the pictures in the file
**
**    Return: This fct returns a long which indicates minpixel.
**                     
*/
long    getmri_minpixel(char *file, int number, int mode)
     /* Assumed: Bruker file */
           
         

{
  char header[HEADER_SIZE+1],v[16],*file_procs,*file_reco,t[16];
  FILE *fp;
  long format;

  format = is_mri(file,NULL);

  /*    For Brucker file format   */
  if(format == BRUKER_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0l);
      }

      fseek(fp,0l,SEEK_SET);
      fread(header,sizeof(char),sizeof(header),fp);
      strcpy(v,decode_pos((unsigned char *)header,MVAR_OFF+3*(2+6),ENTIER));
      fclose(fp);
      return(atol(v));
    }

#ifdef __GTK_GUI
  /*  For PARAVISION format   */
  if(format == PARAVISION_FORMAT)
    {
      file_procs=getfilename_procs(file);
      strcpy(v, getheader_interfile(file_procs, "##$YMIN_p=", ENTIER, 0));
      FREE(file_procs);

      /*cas ou c'est une image de phase*/
      file_reco=getfilename_reco(file);
      strcpy(t, getheader_interfile(file_reco, "##$RECO_image_type=", ASCII8B, 0));
      if (strcmp(t,"PHASE_IMAGE")==0)
        {
	  strcpy(v,"-2147483647");
        }
      FREE(file_reco);

      return(atol(v));
    }
#endif /* __GTK_GUI */

  /*    For Interfile file format   */
  if(format == INTERFILE_FORMAT)
    {
      return(0);
    }

  /*    For raw4 file format   */
  if(format == RAW4_FORMAT)
    {
      return(0);
    }

  /*    For Siemens Magnetom format, Cyril Feruglio, Mai 1998   */
  if(format == SIEMENS_FORMAT)
    {
      return(0);
    }

  /* For SIEMENS Vanderbilt FORMAT Feruglio Cyril Mai 1998 */
  if(format == VANDERBILT_FORMAT)
    {
      strcpy(v, getheader_interfile(file, "Smallest pixel value := ", REEL, 0));
      if(MESG_DEBUG) printf("MinPixel='%s', ",v);
      /* return(atof(v));*/
      return(0);
    }

  /* For IPB 2D file format */
  if(format == IPB_2D_FORMAT)
    {
      /* Gets only first parameter in 3D mode */
      if(mode == MODE_3D)
	number = 0;

      strcpy(v, getheader_interfile(file, "min pixel=", ENTIER_IPB_2D, number));
      return(atol(v));
    }

  /* For IPB 3D file formats */
  if(format == IPB_3D_FORMAT)
    {
      /* Gets an image number corresponding to the slice given by 'number' */
      if(mode == MODE_2D)
	number = getmri_index3D(file, number);

      strcpy(v, getheader_interfile(file, "min pixel=", ENTIER_IPB_2D, number));
      return(atol(v));
    }

  /* For PCX file format */
  if(format == PCX_FORMAT)
    {
      return(0);
    }

  /* For SMIS file format */
  if(format == SMIS_FORMAT)
    return(0);
    

  /* For ACRNEMA 2.0 format */
  if (format == ACRNEMA2L_FORMAT || format == ACRNEMA2B_FORMAT 
      || format==ACRNEMA1B_FORMAT || format==ACRNEMA1L_FORMAT) 
    return 0;

  /* for dicom format */
  if(format == DICMBIMP_FORMAT || format == DICMLIMP_FORMAT ||
     format == DICMBEXP_FORMAT || format == DICMLEXP_FORMAT)
    {
      char *data=(char *)NULL;	
      int min;
      short buffer,buffer2;
      short minpixel;
				
      data=Acrnema_getGroup(file, 0x28, 0x0106, data, format);
      if (!strlen(data)==0) /* Ce groupe existe */ 
	{
	  memcpy(&minpixel, data, 2);
	  if ( _is_bigEndian) minpixel=shEndianConversion(minpixel);
	  free(data);
	  return(minpixel);
	}
      else /* Calcul max a partir de  bits stored */
	{
	  data=Acrnema_getGroup(file, 0x28, 0x0101, data, format);
	  if (data!=(char *)NULL)	memcpy(&buffer, data, 2);
	  if ( _is_bigEndian) /* On suppose endian dans is_mri =0 */ 
	    buffer=shEndianConversion(buffer);
	  FREE(data);
	  data=Acrnema_getGroup(file, 0x28, 0x0103, data, format);
	  if (data!=(char *)NULL)	memcpy(&buffer2, data, 2);
	  if ( _is_bigEndian) /* On suppose endian dans is_mri =0 */ 
	    buffer2=shEndianConversion(buffer2);
	  FREE(data);
	  if (buffer2==0) { /* Pas de bit de signe */
	    return (0);
	  }
	  else /* avec bit de signe */
	    {
	      if(buffer==32)
		min=-(int)(pow((double)2,(double)31)-1)-1;
	      else 
		min=-(int)pow((double)2,(double)(buffer)-1)-1;
	    }
	  return (min);
	}
    }
		
  /* for 3D dicom format */
  if(format == DICM3DBIMP_FORMAT || format == DICM3DLIMP_FORMAT ||
     format == DICM3DBEXP_FORMAT || format == DICM3DLEXP_FORMAT)
    {
      char *data=(char *)NULL;	
      short buffer,buffer2;
      int min;
      short minpixel;
				
      data=Acrnema_getGroup(file, 0x28, 0x0108, data, format);
      if (!strlen(data)==0) /* Ce groupe existe */ 
	{
	  memcpy(&minpixel, data, 2);
	  if ( _is_bigEndian) minpixel=shEndianConversion(minpixel);
	  free(data);
	  return(minpixel);
	}
      else /* calcul du min a partir du SmallestImagePixelValue */
	{
	  data=Acrnema_getGroup(file, 0x28, 0x0106, data, format);
	  if (!strlen(data)==0) /* Ce groupe existe et le groupe 0x0108 n'existe pas*/ 
	    {
	      memcpy(&minpixel, data, 2);
	      if ( _is_bigEndian) minpixel=shEndianConversion(minpixel);
	      free(data);
	      return(minpixel);
	    }
	  else /* Calcul max a partir de  bits stored */
	    { 
	      data=Acrnema_getGroup(file, 0x28, 0x0101, data, format);
	      if (data!=(char *)NULL)	memcpy(&buffer, data, 2);
	      if ( _is_bigEndian) /* On suppose endian dans is_mri =0 */ 
		buffer=shEndianConversion(buffer);
	      FREE(data);
	      data=Acrnema_getGroup(file, 0x28, 0x0103, data, format);
	      if (data!=(char *)NULL)	memcpy(&buffer2, data, 2);
	      if ( _is_bigEndian) /* On suppose endian dans is_mri =0 */ 
		buffer2=shEndianConversion(buffer2);
	      FREE(data);
	      if (buffer2==0) { /* Pas de bit de signe */
		return (0);
	      }
	      else /* avec bit de signe */
		{
		  if(buffer==32)
		    min=-(int)(pow((double)2,(double)31)-1)-1;
		  else 
		    min=-(int)pow((double)2,(double)(buffer)-1)-1;
		}
	      return (min);
	    } 
	}
    }
  /* end of 3D dicom format*/

  /* For ECAT7 format */
  if(format==ECAT7_FORMAT) return 0;

    
  return(0);
}

/* --getmri_maxpixel() ------------------------------------------------    
**
**    Gets maxpixel for a bruker picture
**
**    Usage:    getmri_maxpixel(bruker_file,number,number)
**        file: char*: Represents a file. (May contains path )
**        number: int: number of the picture in the file
**        mode: int: mode to read the pictures in the file
**
**    Return: This fct returns a long which indicates maxpixel.
**                     
*/
long    getmri_maxpixel(char *file, int number, int mode)
     /* Assumed: Bruker file */
           
         
{
  char header[HEADER_SIZE+1],v[16],t[16];
  char *file_procs;
  char *file_reco;
  FILE *fp;
  int format;

  format = is_mri(file,NULL);

  /*    For Brucker file format   */
  if(format == BRUKER_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0);
      }

      fseek(fp,0l,SEEK_SET);
      fread(header,sizeof(char),sizeof(header),fp);

      strcpy(v,decode_pos((unsigned char *)header,MVAR_OFF+3*(2+5),ENTIER));
      fclose(fp);
      return(atol(v));
    }

#ifdef __GTK_GUI
  /*  For PARAVISION format   */
  if(format == PARAVISION_FORMAT)
    {
      file_procs=getfilename_procs(file);
      strcpy(v, getheader_interfile(file_procs, "##$YMAX_p=", ENTIER, 0));
      FREE(file_procs);

      /* cas ou c'est une image de phase */
      file_reco=getfilename_reco(file);
      strcpy(t, getheader_interfile(file_reco, "##$RECO_image_type=", ASCII8B, 0));
      if (strcmp(t,"PHASE_IMAGE")==0)
        {
	  strcpy(v,"2147483647");
        }
      FREE(file_reco);
      return(atol(v));
    }
#endif /* __GTK_GUI */

  /*    For Interfile file format   */
  if(format == INTERFILE_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0);
      }
      fclose(fp);
      return(255);
    }

  /*    For raw4 file format   */
  if(format == RAW4_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0);
      }
      fclose(fp);
      return(100000);
    }

  /* For Siemens Magnetom format, Feruglio Cyril 04-05-1998 */
  if(format == SIEMENS_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0);
      }
      fclose(fp);
      return(65536);
    }

  /* For SIEMENS Vanderbilt FORMAT Feruglio Cyril Mai 1998 */
  if(format == VANDERBILT_FORMAT)
    {
      strcpy(v, getheader_interfile(file, "Largest pixel value := ", REEL, 0));
      if(MESG_DEBUG) printf("MinPixel='%s', ",v);
      return(atol(v)*2048);
    }

  /*   For IPB 2D file format */
  if(format == IPB_2D_FORMAT)
    {
      /* Gets only first parameter in 3D mode */
      if(mode == MODE_3D)
	number = 0;

      strcpy(v, getheader_interfile(file, "max pixel=", ENTIER_IPB_2D, 0));
      return(atol(v));
    }

  /*   For IPB 3D file format */
  if(format == IPB_3D_FORMAT)
    {
      /* Gets an image number corresponding to the slice given by 'number' */
      if(mode == MODE_2D)
	number = getmri_index3D(file, number);

      strcpy(v, getheader_interfile(file, "max pixel=", ENTIER_IPB_2D, 0));
      return(atol(v));
    }

  /*    For pcx file format   */
  if(format == PCX_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0);
      }
      fclose(fp);
      return(255);
    }

  /* For SMIS file format */
  if(format == SMIS_FORMAT) {
    /* 12 bits acquisition */
    return(4095);
  }

	
  /* For ACRNEMA 2.0 format */
  if (format == ACRNEMA2L_FORMAT || format == ACRNEMA2B_FORMAT 
      || format==ACRNEMA1B_FORMAT || format==ACRNEMA1L_FORMAT)
    {
      char *data=(char *)NULL;	
      int x,y,max;
		
      x=getmri_width(file, 1, MODE_2D);
      y=getmri_height(file, 1, MODE_2D);
		
      /* on d�ecte la taille de la zone de donn�s */
      data=Acrnema_getGroup(file, 0x28, 0x0101, data, format);
      max=(int)(pow((double)2,(int)((atol(data)/(x*y))*8)));   
      free (data);
      return (max);
    }		
	
  /* for dicom 2D format */
  if(format == DICMBIMP_FORMAT || format == DICMLIMP_FORMAT ||
     format == DICMBEXP_FORMAT || format == DICMLEXP_FORMAT)
    {
      char *data=(char *)NULL;	
      int max;
      short buffer,buffer2;
      short maxpixel;
				
      data=Acrnema_getGroup(file, 0x28, 0x0107, data, format);
      if (!strlen(data)==0) /* Ce groupe existe */ 
	{
	  memcpy(&maxpixel, data, 2);
	  if ( _is_bigEndian) maxpixel=shEndianConversion(maxpixel);
	  free(data);
	  return(maxpixel);
	}
      else /* Calcul max a partir de  bits stored */
	{
				
	  data=Acrnema_getGroup(file, 0x28, 0x0101, data, format);
	  if (data!=(char *)NULL)	memcpy(&buffer, data, 2);
	  if ( _is_bigEndian) /* On suppose endian dans is_mri =0 */ 
	    buffer=shEndianConversion(buffer);
	  FREE(data);
	  data=Acrnema_getGroup(file, 0x28, 0x0103, data, format);
	  if (data!=(char *)NULL)	memcpy(&buffer2, data, 2);
	  if ( _is_bigEndian) /* On suppose endian dans is_mri =0 */ 
	    buffer2=shEndianConversion(buffer2);
	  FREE(data);
	  if (buffer2==0) { /* Pas de bit de signe */
	    max=(int)pow((double)2,(double)(buffer))-1;
	  }
	  else
	    {
	      if(buffer==32)
		max=(int)(pow((double)2,(double)31)-1)-1;
	      else 
		max=(int)pow((double)2,(double)(buffer)-1)-1;
	    }
	  return (max);
	}
    }	   		
		
  /* for dicom 3D format */
  if(format == DICM3DBIMP_FORMAT || format == DICM3DLIMP_FORMAT ||
     format == DICM3DBEXP_FORMAT || format == DICM3DLEXP_FORMAT)
    {
      char *data=(char *)NULL;	
      short maxpixel;
				
      data=Acrnema_getGroup(file, 0x28, 0x0109, data, format);
      if (!strlen(data)==0) /* Ce groupe existe */ 
	{
	  memcpy(&maxpixel, data, 2);
	  if ( _is_bigEndian) maxpixel=shEndianConversion(maxpixel);
	  free(data);
	  return(maxpixel);
	}
      else /* calcul de maxpixel a partir de LargestImagePixelValue */
	{
	  data=Acrnema_getGroup(file, 0x28, 0x0107, data, format);
	  if (!strlen(data)==0) /* Ce groupe existe et le group 0x109 n'existe pas */ 
	    {
	      memcpy(&maxpixel, data, 2);
	      if ( _is_bigEndian) maxpixel=shEndianConversion(maxpixel);
	      free(data);
	      return(maxpixel);
	    }
	  else /* Calcul max a partir de  bits stored */
	    {		  
	      /* finalement on lui fait retourne une valeur 0 pour recalcul
		 du maximum dans imx_3d */
	      return (0);
	      /*
		data=Acrnema_getGroup(file, 0x28, 0x0101, data, format);
		if (data!=(char *)NULL)	memcpy(&buffer, data, 2);
		if ( _is_bigEndian) // On suppose endian dans is_mri =0
		buffer=shEndianConversion(buffer);
		FREE(data);
		data=Acrnema_getGroup(file, 0x28, 0x0103, data, format);
		if (data!=(char *)NULL)	memcpy(&buffer2, data, 2);
		if ( _is_bigEndian) // On suppose endian dans is_mri =0  
		buffer2=shEndianConversion(buffer2);
		FREE(data);
		if (buffer2==0) { // Pas de bit de signe 
		max=(int)pow((double)2,(double)(buffer))-1;
		}
		else
		{
		if(buffer==32)
		max=(int)(pow((double)2,(double)31)-1)-1;
		else 
		max=(int)pow((double)2,(double)(buffer)-1)-1;
		}
		return (max);
			
	      */
			
	    }
	}
    }		
		
  /* For ECAT7 format */
  if (format == ECAT7_FORMAT){
    char* data=(char*)malloc(10);
    int datType=atoi(Ecat7_getFromMatrix(file,number,0,2,data,"int"));
    FREE(data);
    switch(datType-1){
    case 0:
      return (long int)pow((double)2,(double)8);
      break;
    case 1:
    case 5:
      return (long int)pow((double)2,(double)15);
      break;
    case 2:
    case 3:
    case 4:
      return (long int)pow((double)2,(double)31);
      break;
    }	
  }



  return(0);
}

/* --getmri_cutoff_min() ------------------------------------------------    
**
**    Gets cutoff_min for a picture
**
**    Usage:    getmri_cutoff_min(file,number,number)
**        file: char*: Represents a file. (May contains path )
**        number: int: number of the picture in the file
**        mode: int: mode to read the pictures in the file
**
**    Return: This fct returns a long which indicates minpixel.
**                     
*/
long    getmri_cutoff_min(char *file, int number, int mode)
     /* Assumed: Bruker file */
           
         

{
  char header[HEADER_SIZE+1],v[16],*file_procs,*file_reco,t[16];
  FILE *fp;
  long format;
  long value;

  format = is_mri(file,NULL);

  /*    For Brucker file format   */
  if(format == BRUKER_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0l);
      }

      fseek(fp,0l,SEEK_SET);
      fread(header,sizeof(char),sizeof(header),fp);
      strcpy(v,decode_pos((unsigned char *)header,MVAR_OFF+3*(2+6),ENTIER));
      fclose(fp);
      return(atol(v));
    }

#ifdef __GTK_GUI
  /*  For PARAVISION format   */
  if(format == PARAVISION_FORMAT)
    {
      file_procs=getfilename_procs(file);
      strcpy(v, getheader_interfile(file_procs, "##$YMIN_p=", ENTIER, 0));
      FREE(file_procs);

      /*cas ou c'est une image de phase*/
      file_reco=getfilename_reco(file);
      strcpy(t, getheader_interfile(file_reco, "##$RECO_image_type=", ASCII8B, 0));
      if (strcmp(t,"PHASE_IMAGE")==0)
        {
	  strcpy(v,"-2147483647");
        }
      FREE(file_reco);
      return(atol(v));
    }
#endif /* __GTK_GUI */

  /*    For Interfile file format   */
  if(format == INTERFILE_FORMAT)
    {
      return(0);
    }

  /*    For raw4 file format   */
  if(format == RAW4_FORMAT)
    {
      return(0);
    }

  /*    For Siemens Magnetom format, Cyril Feruglio, Mai 1998   */
  if(format == SIEMENS_FORMAT)
    {
      return(0);
    }

  /* For SIEMENS Vanderbilt FORMAT Feruglio Cyril Mai 1998 */
  if(format == VANDERBILT_FORMAT)
    {
      strcpy(v, getheader_interfile(file, "Smallest pixel value := ", REEL, 0));
      if(MESG_DEBUG) printf("MinPixel='%s', ",v);
      /* return(atof(v));*/
      return(0);
    }

  /*   For IPB 2D file formats */
  if(format == IPB_2D_FORMAT)
    {
      /* Gets only first parameter in 3D mode */
      if(mode == MODE_3D)
	number = 0;

      strcpy(v, getheader_interfile(file, "cutoff min=", ENTIER_IPB_2D, number));
      if((value = atol(v)) == 0)
	return(getmri_minpixel(file, number, MODE_2D));
      else
	return(value);
    }

  /*   For IPB 3D file formats */
  if(format == IPB_3D_FORMAT)
    {
      /* Gets an image number corresponding to the slice given by 'number' */
      if(mode == MODE_2D)
	number = getmri_index3D(file, number);

      strcpy(v, getheader_interfile(file, "cutoff min=", ENTIER_IPB_2D, number));
      if((value = atol(v)) == 0)
	return(getmri_minpixel(file, number, MODE_2D));
      else
	return(value);
    }

  /*    For PCX file format   */
  if(format == PCX_FORMAT)
    {
      return(0);
    }

  /* For SMIS file format */
  if(format == SMIS_FORMAT) 
    return(0);

  /* For ACRNEMA file format */
  if (format == ACRNEMA2L_FORMAT || format == ACRNEMA2B_FORMAT 
      || format==ACRNEMA1B_FORMAT || format==ACRNEMA1L_FORMAT) 
    return 0;

  /* for dicom 2Dformat */
  if(format == DICMBIMP_FORMAT || format == DICMLIMP_FORMAT ||
     format == DICMBEXP_FORMAT || format == DICMLEXP_FORMAT)
    return(getmri_minpixel(file, 1, MODE_3D));

  /* for dicom 3D format */
  if(format == DICM3DBIMP_FORMAT || format == DICM3DLIMP_FORMAT ||
     format == DICM3DBEXP_FORMAT || format == DICM3DLEXP_FORMAT)
    return(getmri_minpixel(file, 1, MODE_3D));

  /* For ECAT7 format */
  if(format==ECAT7_FORMAT) return 0;


  return(0);
}

/* --getmri_cutoff_max() ------------------------------------------------    
**
**    Gets cutoff_max for a picture
**
**    Usage:    getmri_cutoff_max(bruker_file,number,number)
**        file: char*: Represents a file. (May contains path )
**        number: int: number of the picture in the file
**        mode: int: mode to read the pictures in the file
**
**    Return: This fct returns a long which indicates maxpixel.
**                     
*/
long    getmri_cutoff_max(char *file, int number, int mode)
     /* Assumed: Bruker file */
           
         
{
  char header[HEADER_SIZE+1],v[16],t[16];
  char *file_procs;
  char *file_reco;
  FILE *fp;
  long format;
  long value;

  format = is_mri(file,NULL);

  /*    For Brucker file format   */
  if(format == BRUKER_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0l);
      }

      fseek(fp,0l,SEEK_SET);
      fread(header,sizeof(char),sizeof(header),fp);

      strcpy(v,decode_pos((unsigned char *)header,MVAR_OFF+3*(2+5),ENTIER));
      fclose(fp);
      return(atol(v));
    }

#ifdef __GTK_GUI
  /*  For PARAVISION format   */
  if(format == PARAVISION_FORMAT)
    {
      file_procs=getfilename_procs(file);
      strcpy(v, getheader_interfile(file_procs, "##$YMAX_p=", ENTIER, 0));
      FREE(file_procs);

      /*cas ou c'est une image de phase*/
      file_reco=getfilename_reco(file);
      strcpy(t, getheader_interfile(file_reco, "##$RECO_image_type=", ASCII8B, 0));
      if (strcmp(t,"PHASE_IMAGE")==0)
        {
	  strcpy(v,"2147483647");
        }
      FREE(file_reco);
      return(atol(v));
    }
#endif /* __GTK_GUI */

  /*    For Interfile file format   */
  if(format == INTERFILE_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0l);
      }
      fclose(fp);
      return(255l);
    }

  /*    For raw4 file format   */
  if(format == RAW4_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0l);
      }
      fclose(fp);
      return(100000l);
    }

  /* For Siemens Magnetom format, Feruglio Cyril 04-05-1998 */
  if(format == SIEMENS_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0l);
      }
      fclose(fp);
      return(65536l);
    }

  /* For SIEMENS Vanderbilt FORMAT Feruglio Cyril Mai 1998 */
  if(format == VANDERBILT_FORMAT)
    {
      strcpy(v, getheader_interfile(file, "Largest pixel value := ", REEL, 0));
      if(MESG_DEBUG) printf("MinPixel='%s', ",v);
      return(atol(v)*2048);
    }

  /*   For IPB 2D file format */
  if(format == IPB_2D_FORMAT)
    {
      /* Gets only first parameter in 3D mode */
      if(mode == MODE_3D)
	number = 0;

      strcpy(v, getheader_interfile(file, "cutoff max=", ENTIER_IPB_2D, number));
      if((value = atol(v)) == 0)
	return(getmri_maxpixel(file, number, MODE_2D));
      else
	return(value);
    }

  /*   For IPB 3D file format */
  if(format == IPB_3D_FORMAT)
    {
      /* Gets an image number corresponding to the slice given by 'number' */
      if(mode == MODE_2D)
	number = getmri_index3D(file, number);

      strcpy(v, getheader_interfile(file, "cutoff max=", ENTIER_IPB_2D, number));
      if((value = atol(v)) == 0)
	return(getmri_maxpixel(file, number, MODE_2D));
      else
	return(value);
    }

  /*    For pcx file format   */
  if(format == PCX_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0);
      }
      fclose(fp);
      return(255);
    }

  if(format == SMIS_FORMAT) {
    /* 12 bits acquisition */
    return(4095);
  }

  if (format == ACRNEMA2L_FORMAT || format == ACRNEMA2B_FORMAT 
      || format==ACRNEMA1B_FORMAT || format==ACRNEMA1L_FORMAT)
    {
      char *data=(char *)NULL;	
      int x,y,max;

      x=getmri_width(file, 1, MODE_2D);
      y=getmri_height(file, 1, MODE_2D);
      /* on d�ecte la taille de la zone de donn�s */
      data=Acrnema_getGroup(file, 0x28, 0x0101, data, format);
      max=(int)pow((double)2,(int)(atol(data)/(x*y)*8));	/* tailleTotale/(x*y)*8 */
      free (data);
      return (max);
    }		
		
  if(format == DICMBIMP_FORMAT || format == DICMLIMP_FORMAT ||
     format == DICMBEXP_FORMAT || format == DICMLEXP_FORMAT)
    {
      return(getmri_maxpixel(file, 1, MODE_3D));
    }	
	
	
  /* for dicom 3D format */
  if(format == DICM3DBIMP_FORMAT || format == DICM3DLIMP_FORMAT ||
     format == DICM3DBEXP_FORMAT || format == DICM3DLEXP_FORMAT)
    {
      return(getmri_maxpixel(file, 1, MODE_3D));
    }		

  /* For ECAT7 format */
  if (format == ECAT7_FORMAT){
    char* data=(char*)malloc(10);
    int datType=atoi(Ecat7_getFromMatrix(file,number,0,2,data,"int"));
    FREE(data);
    switch(datType-1){
    case 0:
      return (long int)pow((double)2,(double)8);
      break;
    case 1:
    case 5:
      return (long int)pow((double)2,(double)15);
      break;
    case 2:
    case 3:
    case 4:
      return (long int)pow((double)2,(double)31);
      break;
    }	
  }

  return(0);
}

/* --getmri_icomp() ------------------------------------------------    
**
**    Gets icomp for a bruker picture
**
**    Usage:    getmri_icomp(bruker_file,number,number)
**        file: char*: Represents a file. (May contains path )
**        number: int: number of the picture in the file
**        mode: int: mode to read the pictures in the file
**
**    Return: This fct returns a long which indicates icomp.
**                     
*/
float    getmri_icomp(char *file, int number, int mode)
     /* Assumed: Bruker file */
           
         
{
  char header[HEADER_SIZE+1],v[16],t[16];
  char *file_procs;
  char *file_reco;
  FILE *fp;
  long format;

  format = is_mri(file,NULL);

  /*    For Brucker file format   */
  if(format == BRUKER_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0l);
      }

      fseek(fp,0l,SEEK_SET);
      fread(header,sizeof(char),sizeof(header),fp);
      strcpy(v,decode_pos((unsigned char *)header,MVAR_OFF+3*(2+7),ENTIER));
      fclose(fp);
      return((float)(atof(v)));
    }

#ifdef __GTK_GUI
  /*  For PARAVISION format   */
  if(format == PARAVISION_FORMAT)
    {
      file_procs=getfilename_procs(file);
      strcpy(v, getheader_interfile(file_procs, "##$NC_proc=", ENTIER, 0));
      FREE(file_procs);

      /*cas ou c'est une image de phase*/
      file_reco=getfilename_reco(file);
      strcpy(t, getheader_interfile(file_reco, "##$RECO_image_type=", ASCII8B, 0));
      if (strcmp(t,"PHASE_IMAGE")==0)
        {
	  strcpy(v,"-29.34850387");
        }
      FREE(file_reco);

      return((float)(atof(v)));
    }
#endif /* __GTK_GUI */

  if(format == INTERFILE_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0l);
      }
      fclose(fp);
      return(0l);
    }

  if(format == RAW4_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0l);
      }
      fclose(fp);
      return(0l);
    }

  /*   For Siemens Magnetom Vision format, Cyril Feruglio Mai 1998 */
  if(format == SIEMENS_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0l);
      }
      fclose(fp);
      return(0l);
    }

  /*   For Siemens Vanderbilt format, Cyril Feruglio Mai 1998 */
  if(format == VANDERBILT_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0l);
      }
      fclose(fp);
      return(0l);
    }

  /*   For IPB 2D file format */
  if(format == IPB_2D_FORMAT)
    {
      /* Gets only first parameter in 3D mode */
      if(mode == MODE_3D)
	number = 0;

      strcpy(v, getheader_interfile(file, "icomp=", REEL_IPB_2D, number));
      return((float)(atof(v)));
    }

  /*   For IPB 3D file format */
  if(format == IPB_3D_FORMAT)
    {
      /* Gets an image number corresponding to the slice given by 'number' */
      if(mode == MODE_2D)
	number = getmri_index3D(file, number);

      strcpy(v, getheader_interfile(file, "icomp=", REEL_IPB_2D, number));
      return((float)(atof(v)));
    }

  if(format == PCX_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0);
      }
      fclose(fp);
      return(0);
    }

  /* For SMIS file format */
  if(format == SMIS_FORMAT) 
    return(0);

  /* For ACRNEMA 2.0 format */
  if (format == ACRNEMA2L_FORMAT || format == ACRNEMA2B_FORMAT 
      || format==ACRNEMA1B_FORMAT || format==ACRNEMA1L_FORMAT) 
    return 0;
		
  /* for dicom format */
  if(format == DICMBIMP_FORMAT || format == DICMLIMP_FORMAT ||
     format == DICMBEXP_FORMAT || format == DICMLEXP_FORMAT)
    return 0;
	
  /* for dicom 3D format */
  if(format == DICM3DBIMP_FORMAT || format == DICM3DLIMP_FORMAT ||
     format == DICM3DBEXP_FORMAT || format == DICM3DLEXP_FORMAT)
    return 0;
	   
  /* For ECAT 7 format */
  if (format == ECAT7_FORMAT) return 0;

  return(0);
}

/* --getmri_fov() ------------------------------------------------    
**
**    Gets fov (Field of view) for a picture
**
**    Usage:    getmri_fov(bruker_file)
**        file: char*: Represents a file. (May contains path )
**
**    Return: This fct returns a long which indicates maxpixel.
**                     
*/
float    getmri_fov(char *file)
     /* Assumed: Bruker file */
{
  char header[HEADER_SIZE+1],v[16];
  char *file_d3proc;
  FILE *fp;
  double fov;
  long format;

  format = is_mri(file,NULL);

  /*    For Brucker file format   */
  if(format == BRUKER_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) /* Open file error... */
	return(0);

      fseek(fp, 0, SEEK_SET);
      fread(header,sizeof(char),sizeof(header),fp);

      strcpy(v,decode_pos((unsigned char *)header,MVAR_OFF+3*(2+39),REEL));
      fclose(fp);
      return(atof(v));
    }

#ifdef __GTK_GUI
  /*  For PARAVISION format   */
  if(format == PARAVISION_FORMAT)
    {
      file_d3proc=getfilename_d3proc(file);
      strcpy(v, getheader_interfile(file_d3proc, "##$PR_STA=(", TAB_PARA_REEL, 0));
      FREE(file_d3proc);
      return(atof(v));
    }
#endif /* __GTK_GUI */

  /*    For Interfile file format   */
  if(format == INTERFILE_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) /* Open file error... */
	return(0);

      fclose(fp);
      return(1);
    }

  /*    For raw4 file format   */
  if(format == RAW4_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) /* Open file error... */
	return(0);

      fclose(fp);
      return(1);
    }

  /*    For Siemens Magnetom Vision format, Cyril Feruglio Mai 1998 */

  if(format == SIEMENS_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) /* Open file error... */
	return(0);

      /* offset 3744: FOV Row, offset 3752: FOV Column*/
      /*protoize:???*/ /*ok?*/
      fov = getdouble_header(file/*fp*/, 3744);
      return((float)fov);
    }

  /*    For SIEMENS Vanderbilt format   */
  if(format == VANDERBILT_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) /* Open file error... */
	return(0);

      fclose(fp);
      return(1);
    }

  /*    For IPB file formats   */
  if(format == IPB_2D_FORMAT || format == IPB_3D_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) /* Open file error... */
	return(0);

      fclose(fp);
      return(1);
    }

  /*    For PCX file format   */
  if(format == PCX_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) /* Open file error... */
	return(0);

      fclose(fp);
      return(1);
    }

  /* For smis format */
  if(format == SMIS_FORMAT) 
    return(1);

  /* For ACRNEMA 2.0 format */
  if (format == ACRNEMA2L_FORMAT || format == ACRNEMA2B_FORMAT 
      || format==ACRNEMA1B_FORMAT || format==ACRNEMA1L_FORMAT) 
    return 1;
		
  /* for dicom format */
  if(format == DICMBIMP_FORMAT || format == DICMLIMP_FORMAT ||
     format == DICMBEXP_FORMAT || format == DICMLEXP_FORMAT)
    return 1;

  /* for dicom 3D format */
  if(format == DICM3DBIMP_FORMAT || format == DICM3DLIMP_FORMAT ||
     format == DICM3DBEXP_FORMAT || format == DICM3DLEXP_FORMAT)
    return 1;
  /* end of dicom 3D format */

  /* For ECAT7 format */
  if (format == ECAT7_FORMAT)return 1;

  return(0);
}

/* --getmri_dx() ------------------------------------------------    
**
**    Gets dx for a picture (dimension x in mm)
**
**    Usage:    getmri_dx(bruker_file,number,number)
**        file: char*: Represents a file. (May contains path )
**        number: int: number of the picture in the file
**        mode: int: mode to read the pictures in the file
**
**    Return: This fct returns a long which indicates maxpixel.
**                     
*/
float    getmri_dx(char *file, int number, int mode)
     /* Assumed: Bruker file */
           
         
{
  float fov;
  int wdth;
  FILE *fp;
  char v[16],*file_d3proc;
  double dx;
  long format;

  format = is_mri(file,NULL);

  /*    For Brucker file format   */
  if(format == BRUKER_FORMAT)
    {
      fov=getmri_fov(file);
      wdth=getmri_width(file,1, MODE_2D);
      return(10.*fov/wdth);
    }

#ifdef __GTK_GUI
  /*  For PARAVISION format   */
  if(format == PARAVISION_FORMAT)
    {
      file_d3proc=getfilename_d3proc(file);
      strcpy(v, getheader_interfile(file_d3proc, "##$PR_STA=", REEL_IPB_2D, 1));
      if (atof(v)==0.0)
	strcpy(v, getheader_interfile(file_d3proc, "##$PR_STA=(", TAB_PARA_REEL, 0));

      FREE(file_d3proc);
      wdth=getmri_width(file,1, MODE_2D);
      fov=10.*atof(v)/wdth;
      return(fov);
    }
#endif /* __GTK_GUI */

  /*    For Interfile file format   */
  if(format == INTERFILE_FORMAT)
    {
      strcpy(v, getheader_interfile(file, "!scaling factor (mm/pixel) [1]:=", REEL, 0));
      return(atof(v));
    }

  /*    For raw4 file format   */
  if(format == RAW4_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) /* Open file error... */
	return(0);

      fclose(fp);
      return(1);
    }

  /*   For Siemens Magnetom Vision format, Cyril Feruglio Mai 1998  */
  if(format == SIEMENS_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) /* Open file error... */
	return(0);

      dx = getdouble_header(file, 5000);
      fclose(fp);
      return(dx);
    }

  /* For SIEMENS Vanderbilt FORMAT Feruglio Cyril Mai 1998 */
  if(format == VANDERBILT_FORMAT)
    {
      strcpy(v, getheader_interfile(file, "Pixel size := ", REEL, 0));
      printf("dx='%s', ",v);
      return(atof(v));
    }

  /*   For IPB 2D file formats */
  if(format == IPB_2D_FORMAT)
    {
      /* Gets only first parameter in 3D mode */
      if(mode == MODE_3D)
	number = 0;

      strcpy(v, getheader_interfile(file, "dx=", REEL_IPB_2D, number));
      return(atof(v));
    }

  /*   For IPB 3D file formats */
  if(format == IPB_3D_FORMAT)
    {
      /* Gets an image number corresponding to the slice given by 'number' */
      if(mode == MODE_2D)
	number = getmri_index3D(file, number);

      strcpy(v, getheader_interfile(file, "dx=", REEL_IPB_2D, number));
      return(atof(v));
    }

  /*    For PCX file format   */
  if(format == PCX_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) /* Open file error... */
	return(0.);

      fclose(fp);
      return(1.);
    }

  /* For SMIS file format */
  if(format == SMIS_FORMAT) {
    float FOVx;
    float RESx;
    char tmp[64];

    sprintf(tmp,":IM_FOV");
    strcpy(v, getheader_interfile(file, tmp, REEL, 0));
    FOVx = atof(v);

    sprintf(tmp,":IM_RESOLUTION");
    strcpy(v, getheader_interfile(file, tmp, REEL, 0));
    RESx = atof(v);

    return(FOVx/RESx);
  }

			
  /* For ACRNEMA 2.0 format */
  if (format == ACRNEMA2L_FORMAT || format == ACRNEMA2B_FORMAT 
      || format==ACRNEMA1B_FORMAT || format==ACRNEMA1L_FORMAT)
    {
      char *st=(char *)NULL;
      int y;

      st=Acrnema_getGroup(file, 0x18, 0x50, st,format);
      y=atol(st);
      free (st);
      return (y);
    }
	
  /* for dicom format -> no dx in documentation */
  if(format == DICMBIMP_FORMAT || format == DICMLIMP_FORMAT ||
     format == DICMBEXP_FORMAT || format == DICMLEXP_FORMAT)
    { char *data=NULL,*data2;
    float dx;
    data=Acrnema_getGroup(file, 0x28, 0x30, data, format);
	if (!data || !strcmp(data, "")) {
		PUT_ERR("Erreur de format du fichier DICOM (2D)");
		return 0; }
		
    data2=(char*)strtok(data,"\\");
    data2=(char*)strtok(NULL,"\\");
    dx=atof(data2);
    FREE(data);
    return(dx);
    }
	
  /* for dicom 3D format  */
  if(format == DICM3DBIMP_FORMAT || format == DICM3DLIMP_FORMAT ||
     format == DICM3DBEXP_FORMAT || format == DICM3DLEXP_FORMAT)
    { char *data=NULL,*data2;
    float dx;
	data=Acrnema_getGroup(file, 0x28, 0x30, data, format);
	if (!data || !strcmp(data, "")) {
		PUT_ERR("Erreur de format du fichier DICOM (3D)");
		return 0; }
    data2=(char*)strtok(data,"\\");
    data2=(char*)strtok(NULL,"\\");
    dx=atof(data2);
    FREE(data);
    return(dx);
    }
  /* end of dicom 3D format */
	
  /* For ECAT7 Format */
  if (format == ECAT7_FORMAT){
    char* data=(char*)malloc(10);
    float ret=10*atof(Ecat7_getFromMatrix(file,number,0x22,4,data,"float"));
    FREE(data);
    return ret;
  }

  return(0);
}

/* --getmri_dy() ------------------------------------------------    
**
**    Gets dy for a picture (dimension y in mm)
**
**    Usage:    getmri_dy(bruker_file,number,number)
**        file: char*: Represents a file. (May contains path )
**        number: int: number of the picture in the file
**        mode: int: mode to read the pictures in the file
**
**    Return: This fct returns a long which indicates maxpixel.
**                     
*/
float    getmri_dy(char *file, int number, int mode)
     /* Assumed: Bruker file */
           
         
{
  float fov;
  int hght;
  FILE *fp;
  char dx[16];
  char ligne[26];
  char v[16],*file_d3proc;
  double dy;
  long format;

  format = is_mri(file,NULL);

  /*    For Brucker file format   */
  if(format == BRUKER_FORMAT)
    {
      fov=getmri_fov(file);
      hght=getmri_height(file,1, MODE_2D);
      return(10.*fov/hght);
    }

#ifdef __GTK_GUI
  /*  For PARAVISION format   */
  if(format == PARAVISION_FORMAT)
    {
      file_d3proc=getfilename_d3proc(file);
      strcpy(v, getheader_interfile(file_d3proc, "##$PR_STA=", REEL_IPB_2D, 2));
      if (atof(v)==0.0)          
	strcpy(v, getheader_interfile(file_d3proc, "##$PR_STA=(", TAB_PARA_REEL, 1));
	
      FREE(file_d3proc);
      hght=getmri_height(file,1, MODE_2D);
      fov=10.*atof(v)/hght;
      return(fov);
    }
#endif /* __GTK_GUI */

  /*    For Interfile file format   */
  if(format == INTERFILE_FORMAT)
    {
      strcpy(v, getheader_interfile(file, "!scaling factor (mm/pixel) [2]:=", REEL, 0));
      return((float) atof(v));
    }

  /*    For raw4 file format   */
  if(format == RAW4_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0.);
      }
      fclose(fp);
      return(1.);
    }

  /*   For Siemens Magnetom Vision format, Cyril Feruglio Mai 1998  */
  if(format == SIEMENS_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) /* Open file error... */
	return(0);

      fclose(fp);
      dy = getdouble_header(file, 5008);
      return(dy);
    }

  /* For SIEMENS Vanderbilt FORMAT Feruglio Cyril Mai 1998 */
  if(format == VANDERBILT_FORMAT)
    {
      /*protoize:???*/
      sprintf(dx,"%f", getmri_dx(file, 0/*???*/,MODE_2D));
      strcpy(ligne, "Pixel size := ");
      strcat(ligne, dx);
      strcat(ligne, " : ");
      strcpy(v, getheader_interfile(file, ligne, REEL, 0));
      printf("\n ligne: '%s' \n",ligne);
      printf(" dy='%s', ",v);
      return(atof(v));
    }

  /*   For IPB 2D file format */
  if(format == IPB_2D_FORMAT)
    {
      /* Gets only first parameter in 3D mode */
      if(mode == MODE_3D)
	number = 0;

      strcpy(v, getheader_interfile(file, "dy=", REEL_IPB_2D, number));
      return(atof(v));
    }

  /*   For IPB 3D file format */
  if(format == IPB_3D_FORMAT)
    {
      /* Gets an image number corresponding to the slice given by 'number' */
      if(mode == MODE_2D)
	number = getmri_index3D(file, number);

      strcpy(v, getheader_interfile(file, "dy=", REEL_IPB_2D, number));
      return(atof(v));
    }

  /*    For PCX file format   */
  if(format == PCX_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) /* Open file error... */
	return(0);

      fclose(fp);
      return(1);
    }

  /* For SMIS file format */
  if(format == SMIS_FORMAT) {
    float FOVx,FOVy;
    float RESx,RESy;
    char tmp[64];

    sprintf(tmp,":IM_FOV");
    strcpy(v, getheader_interfile(file, tmp, REEL, 0));
    FOVx = atof(v);
    sprintf(tmp,"%s %.2f",tmp,atof(v));
    strcpy(v, getheader_interfile(file, tmp, REEL, 0));
    FOVy = atof(v);

    sprintf(tmp,":IM_RESOLUTION");
    strcpy(v, getheader_interfile(file, tmp, REEL, 0));
    RESx = atof(v);
    //strcpy(v, getheader_interfile(file, ":IM_RESOLUTION", REEL_IPB_2D, 2)); cette commande ne donne pas le bon resultat  d'ou l'astuce
    sprintf(tmp,"%s %.2f,",tmp,atof(v));
    strcpy(v, getheader_interfile(file, tmp, REEL, 0));
    RESy = atof(v);

    sprintf(tmp,":IM_RESOLUTION");

    return(FOVy/RESy);
  }

	
  /* For ACRNEMA 2.0 format */
  if (format == ACRNEMA2L_FORMAT || format == ACRNEMA2B_FORMAT 
      || format==ACRNEMA1B_FORMAT || format==ACRNEMA1L_FORMAT) 
    return 1;

  /* for dicom 2D format */
  if(format == DICMBIMP_FORMAT || format == DICMLIMP_FORMAT ||
     format == DICMBEXP_FORMAT || format == DICMLEXP_FORMAT)
    { char *data=NULL,*data2;
    float dy;
    data=Acrnema_getGroup(file, 0x28, 0x30, data, format);
	data=Acrnema_getGroup(file, 0x28, 0x30, data, format);
	if (!data || !strcmp(data, "")) {
		PUT_ERR("Erreur de format du fichier DICOM (2D)");
		return 0; }
    data2=(char*)strtok(data,"\\");
    dy=atof(data2);
    data2=(char*)strtok(NULL,"\\");
    FREE(data);
    return(dy);
    }

  /* for dicom 3D format */
  if(format == DICM3DBIMP_FORMAT || format == DICM3DLIMP_FORMAT ||
     format == DICM3DBEXP_FORMAT || format == DICM3DLEXP_FORMAT)
    { char *data=NULL,*data2;
    float dy;
    data=Acrnema_getGroup(file, 0x28, 0x30, data, format);
	data=Acrnema_getGroup(file, 0x28, 0x30, data, format);
	if (!data || !strcmp(data, "")) {
		PUT_ERR("Erreur de format du fichier DICOM (3D)");
		return 0; }
    data2=(char*)strtok(data,"\\");
    dy=atof(data2);
    data2=(char*)strtok(NULL,"\\");
    FREE(data);
    return(dy);
    }
  /* end of dicom 3D format */

  /* For ECAT7 Format */
  if (format == ECAT7_FORMAT){
    char* data=(char*)malloc(10);
    float ret;
    if(mode==MODE_3D)
      ret=10*atof(Ecat7_getFromMatrix(file,number,0x2A,4,data,"float"));
    else
      ret=10*atof(Ecat7_getFromMatrix(file,number,0x26,4,data,"float"));
		
    FREE(data);
    return ret;
  }


  return(0);
}

/* --getmri_dz() ------------------------------------------------    
**
**    Gets dz for a picture (slice thickness in mm)
**
**    Usage:    getmri_dz(bruker_file,number,number)
**        file: char*: Represents a file. (May contains path )
**        number: int: number of the picture in the file
**        mode: int: mode to read the pictures in the file
**
**    Return: This fct returns a long which indicates maxpixel.
**                     
*/
float    getmri_dz(char *file, int number, int mode)
     /* Assumed: Bruker file */
           
         
{
  char header[HEADER_SIZE+1],v[16],*file_d3proc,*file_acqp;
  FILE *fp;
  int dpth,dim_acq,i;
  float fov;
  long format;

  format = is_mri(file,NULL);

  /*    For Brucker file format   */
  if(format == BRUKER_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0l);
      }

      fseek(fp,0l,SEEK_SET);
      fread(header,sizeof(char),sizeof(header),fp);

      strcpy(v,decode_pos((unsigned char *)header,MVAR_OFF+3*(2+41),REEL));
      fclose(fp);
      return(atof(v));
    }

#ifdef __GTK_GUI
  /*  For PARAVISION format   */
  if(format == PARAVISION_FORMAT)
    {
      file_acqp=getfilename_acqp(file);
      strcpy(v, getheader_interfile(file_acqp, "##$ACQ_dim=", ENTIER, 0));
      dim_acq=atol(v);
      FREE(file_acqp);

      switch (dim_acq) {
      case 2:   /* acquisition image 2D */
	file_acqp=getfilename_acqp(file);    
	strcpy(v, getheader_interfile(file_acqp, "##$ACQ_slice_sepn=(", TAB_PARA_REEL, 2));

	//            strcpy(v, getheader_interfile(file_acqp, "##$ACQ_slice_thick=", REEL, 0));
	FREE(file_acqp);
	return(atof(v));
      case 3:   /* acquisition image 3D */

	file_d3proc=getfilename_d3proc(file);
	strcpy(v, getheader_interfile(file_d3proc, "##$PR_STA=", REEL_IPB_2D, 3));
	if (atof(v)==0.0)            
	  strcpy(v, getheader_interfile(file_d3proc, "##$PR_STA=(", TAB_PARA_REEL, 2));
	
	FREE(file_d3proc);
	dpth=getmri_depth(file,1, MODE_2D);
	fov=10.*atof(v)/dpth;
	return(fov);
      }
    }
#endif /* __GTK_GUI */

  /*    For Interfile file format   */
  if(format == INTERFILE_FORMAT)
    {
      strcpy(v, getheader_interfile(file, "slice thickness (pixels):=", ENTIER, 0));
      i=atol(v);
      if (i==0) i=1;
      strcpy(v, getheader_interfile(file, "!scaling factor (mm/pixel) [2]:=", REEL, 0));
      /*  Pour les images de Elscint de l'IPB */
      if (atof(v)==1.7) strcpy(v,"4.5"); 
      if (atof(v)==3.3||atof(v)==3.32) strcpy(v,"5.7"); 
      return((float) atof(v)*i);
    }

  /*    For raw4 file format   */
  if(format == RAW4_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0.);
      }
      fclose(fp);
      return(1.);
    }

  /*    For Siemens magnetom vision format, Cyril Feruglio Mai 1998   */
  if(format == SIEMENS_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0.);
      }
      fclose(fp);
      return(1.);
    }

  /* For SIEMENS Vanderbilt FORMAT Feruglio Cyril Mai 1998 */
  if(format == VANDERBILT_FORMAT)
    {
      strcpy(v, getheader_interfile(file, "Slice thickness := ", REEL, 0));
      printf("dz='%s', ",v);
      return(atof(v));
    }

  /*   For IPB 2D file format */
  if(format == IPB_2D_FORMAT)
    {
      /* Gets only first parameter in 3D mode */
      if(mode == MODE_3D)
	number = 0;

      strcpy(v, getheader_interfile(file, "dz=", REEL_IPB_2D, number));
      return(atof(v));
    }

  /*   For IPB 3D file formats */
  if(format == IPB_3D_FORMAT)
    {
      /* Gets an image number corresponding to the slice given by 'number' */
      if(mode == MODE_2D)
	number = getmri_index3D(file, number);

      strcpy(v, getheader_interfile(file, "dz=", REEL_IPB_2D, number));
      return(atof(v));
    }

  /*    For PCX file format   */
  if(format == PCX_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	return(0.);
      }
      fclose(fp);
      return(1.);
    }

  /* For SMIS file format */
  if(format == SMIS_FORMAT) {
    float FOVx,FOVy,FOVz;
    //		float RESz;
    char tmp[128];

    sprintf(tmp,":IM_FOV");
    strcpy(v, getheader_interfile(file, tmp, REEL, 0));
    FOVx = atof(v);
    sprintf(tmp,"%s %.2f",tmp,atof(v));
    strcpy(v, getheader_interfile(file, tmp, REEL, 0));
    FOVy = atof(v);
    sprintf(tmp,"%s %.2f",tmp,atof(v));
    strcpy(v, getheader_interfile(file, tmp, REEL, 0));
    FOVz = atof(v);

    /*		sprintf(tmp,":IM_TOTAL_SLICES");
		strcpy(v, getheader_interfile(file, tmp, ENTIER, 0));
		RESz = atof(v);
    */

    return(FOVz); // on a pas besoin de faire de calcul puisque FOVz concerne qu'une seule coupe
  }

  /* For ACRNEMA 2.0 format */
  if (format == ACRNEMA2L_FORMAT || format == ACRNEMA2B_FORMAT 
      || format==ACRNEMA1B_FORMAT || format==ACRNEMA1L_FORMAT)
    return 1;
		 
  /* for dicom format */
  if(format == DICMBIMP_FORMAT || format == DICMLIMP_FORMAT ||
     format == DICMBEXP_FORMAT || format == DICMLEXP_FORMAT)
    { char* data=NULL;
    float dz;
    data=Acrnema_getGroup(file, 0x18, 0x88, data, format);
//	data=Acrnema_getGroup(file, 0x28, 0x30, data, format);
	if (!data || !strcmp(data, "")) {
		PUT_ERR("Erreur de format du fichier DICOM (2D)");
		return 0; }
    dz=atof(data);
    FREE(data);
    return(dz);
    }

  /* for dicom 3D format */
  if(format == DICM3DBIMP_FORMAT || format == DICM3DLIMP_FORMAT ||
     format == DICM3DBEXP_FORMAT || format == DICM3DLEXP_FORMAT)
    { char* data=NULL;
    float dz;
    data=Acrnema_getGroup(file, 0x18, 0x88, data, format);
//	data=Acrnema_getGroup(file, 0x28, 0x30, data, format);
	if (!data || !strcmp(data, "")) {
		PUT_ERR("Erreur de format du fichier DICOM (3D)");
		return 0; }
    dz=atof(data);
    FREE(data);
    return(dz);
    }
  /* end of dicom 3D format */

  /* For ECAT7 Format */
  if (format == ECAT7_FORMAT){
    char* data=(char*)malloc(10);
    float ret;
    if(mode==MODE_3D)
      ret=10*atof(Ecat7_getFromMatrix(file,number,0x26,4,data,"float"));
    else
      ret=10*atof(Ecat7_getFromMatrix(file,number,0x2A,4,data,"float"));
		
    FREE(data);
    return ret;
  }

  return(0.);
}

/* --getmri_mask_type() ------------------------------------------------    
**
**    Gets mask_type for a picture (Binary ....)
**
**    Usage:    getmri_mask_type(bruker_file,number,number)
**        file: char*: Represents a file. (May contains path )
**        number: int: number of the picture in the file
**        mode: int: mode to read the pictures in the file
**
**    Return: This fct returns a int which indicates mask_type.
**                     
*/

int  getmri_mask_type(char *file, int number, int mode)
     /* Assumed: Bruker file */
{
  char v[16];
  long format;

  format = is_mri(file,NULL);

  if(format == IPB_3D_FORMAT || format == IPB_2D_FORMAT)
    {
      /* Gets a mask type in file */
      strcpy(v, getheader_interfile(file, "MASK_TYPE=", ENTIER, number));
      return(atoi(v));
    }
  return 0;
}

/* --getmri_endian() ---------------------------------------------------
**
**    Gets endianity  of a bruker file
**              anf interfile format
**
**    Usage:    getmri_endian(file)
**        file: char*: Represents a file witch can contains path
**
**    Return: This fct returns a TRUE if bigendian
**
*/
int     getmri_endian(char *file)
{
  int format;
  int is_filebigEndian;
  char v[16];
  char *file_reco;
       
  format = is_mri(file, NULL);
    
  /*  Recherche type d'ecriture sur disque du fichier data en fonction du format */
  switch (format)
    {

    case DICM3DBIMP_FORMAT :
    case DICM3DLIMP_FORMAT :
    case DICM3DBEXP_FORMAT :
    case DICM3DLEXP_FORMAT :
    case ANALYZE_FORMAT_LITTLE_ENDIAN :
      is_filebigEndian=FALSE;
      break;
    case PARAVISION_FORMAT :
      file_reco= getfilename_reco(file);
      strcpy(v, getheader_interfile(file_reco, "##$RECO_byte_order=", STRING, 0)); 
      if (!strcmp(v,"bigEndian"))
	is_filebigEndian = TRUE;
      else
	is_filebigEndian = FALSE;
      free(file_reco);
      break;
    case IPB_3D_FORMAT :
      strcpy(v,getheader_interfile(file,"endian=",STRING,0));
      if (!strncmp(v,"littleEndian",12))
	is_filebigEndian = FALSE;
      else
	is_filebigEndian = TRUE;
      break;
    default:
      is_filebigEndian=TRUE;
      break;
    }

    
  return is_filebigEndian;
}


/* --getmri_dattype() ---------------------------------------------------
**
**    Gets type of data  a bruker file 
**              anf interfile format
**
**    Usage:    getmri_dattype(file,number,number)
**        file: char*: Represents a file witch can contains path 
**        number: int: number of the picture in the file
**        mode: int: mode to read the pictures in the file
**
**    Return: This fct returns a long which indicates the type of data
**                     
*/
long    getmri_dattype(char *file, int number, int mode)
{
  long dattype;
  FILE *fp;
  char v[16],*file_d3proc;
  int  Version;
  int  Xmin, Ymin, Xmax, Ymax;
  int  BitsPerPixel;
  char mess[80];
  long format;

  format = is_mri(file,NULL);

  /*    For Brucker file format   */
  if(format == BRUKER_FORMAT)
    {
      return(4l);
    }

#ifdef __GTK_GUI                                  
  /*  For PARAVISION format   */
  if(format == PARAVISION_FORMAT)
    {
      file_d3proc=getfilename_d3proc(file);
      strcpy(v, getheader_interfile(file_d3proc, "##$DATTYPE=", STRING, 0));
      if (!strcmp(v,"ip_short"))
	dattype = 2;
      else
        {
          strcpy(v, getheader_interfile(file_d3proc, "##$DATTYPE=", ENTIER, 0));
          FREE(file_d3proc);
          dattype=atol(v)-1;
        }
      return(dattype);
    }
#endif /* __GTK_GUI */

  /*  For Interfile format  */
  if(format == INTERFILE_FORMAT)
    {
      strcpy(v, getheader_interfile(file, "!number of bytes per pixel:=", ENTIER, 0));
      return(atol(v));
    }

  /*    For raw4 file format  A FAIRE  */
  if(format == RAW4_FORMAT)
    {
      return(4l);
    }

  /*    For Siemens Magnetom Vision format, Cyril Feruglio Mai 1998   */
  if(format == SIEMENS_FORMAT)
    {
      return(4l);
    }

  /* For SIEMENS Vanderbilt FORMAT Feruglio Cyril Mai 1998 */
  if(format == VANDERBILT_FORMAT)
    {
      strcpy(v, getheader_interfile(file, "Bits allocated := ", ENTIER, 0));
      printf("Bits Allocated='%s'\n", v);
      return((atol(v))/8);
    }

  /*   For IPB 2D file format */
  if(format == IPB_2D_FORMAT)
    {
      /* Gets only first parameter in 3D mode */
      if(mode == MODE_3D)
	number = 0;

      strcpy(v, getheader_interfile(file, "bits per pixel=", ENTIER_IPB_2D, number));
      return((long)(atol(v)/8));
    }

  /*   For IPB 3D file format */
  if(format == IPB_3D_FORMAT)
    {
      strcpy(v, getheader_interfile(file, "bits per pixel=", ENTIER_IPB_2D, number));
      return((long)(atol(v)/8));
    }

  /*    For pcx file format A FAIRE  */
    
  if(format == PCX_FORMAT)
    {
      fp = p_openr(file);
      //
      // read the PCX header
      //
      if (GetByte(fp) != PCX_MAGIC)
        {
	  return(0);
        }
      Version     = GetByte(fp);  // get version #

      if (GetByte(fp) != 1)       // check for PCX run length encoding
        {
	  sprintf(mess,"%s has unknown encoding scheme \n", file );
	  PUT_MESG(mess);
	  return(0);
        }

      BitsPerPixel= GetByte(fp);
      Xmin        = GetWord(fp);
      Ymin        = GetWord(fp);
      Xmax        = GetWord(fp);
      Ymax        = GetWord(fp);
      //width       = (Xmax - Xmin) + 1;
      //Height      = (Ymax - Ymin) + 1;
      fclose(fp);
      dattype=BitsPerPixel/8;
      return((long)dattype);
    }
    

  /* For SMIS file format */
  if(format == SMIS_FORMAT) 
    return(2);

  /* For ACRNEMA 2.0 format */
  if (format == ACRNEMA2L_FORMAT || format == ACRNEMA2B_FORMAT 
      || format==ACRNEMA1B_FORMAT || format==ACRNEMA1L_FORMAT) {
    char *data=(char *)NULL;	
    unsigned short dataSize;

    data=Acrnema_getGroup(file, 0x28, 0x0101, data, format);
    memcpy(&dataSize, data, 2);
    free (data);
    return (dataSize/8);
  }	

  /* for dicom 2Dformat */
  if(format == DICMBIMP_FORMAT || format == DICMLIMP_FORMAT ||
     format == DICMBEXP_FORMAT || format == DICMLEXP_FORMAT)
    {
      char *data=(char *)NULL;    
      int max;
      short buffer;
	    	   
      data=Acrnema_getGroup(file, 0x28, 0x0101, data, format);
      buffer=atoi(data);
      FREE(data);
      if(buffer==32)max=(int)(pow((double)2,(double)31)-1);
      else max=(int)pow((double)2,(double)(buffer)-1);   
      return (abs(max));
    }	


  /* for dicom 3D format */
  if(format == DICM3DBIMP_FORMAT || format == DICM3DLIMP_FORMAT ||
     format == DICM3DBEXP_FORMAT || format == DICM3DLEXP_FORMAT)
    {
      char *data=(char *)NULL;	

      short data_size;
				
      data=Acrnema_getGroup(file, 0x28, 0x0100, data, format);
      if (data!=(char *)NULL)	memcpy(&data_size, data, 2);
      if ( _is_bigEndian)  
	data_size=shEndianConversion(data_size);
      free (data);
      return (data_size/8);
    }	

  /* For ECAT7 format */
  if (format == ECAT7_FORMAT){
    char* data=(char*)malloc(10);
    int datType=atoi(Ecat7_getFromMatrix(file,number,0,2,data,"int"));
    FREE(data);
    switch(datType+1){
    case 0:
      return 1l;
      break;
    case 1:
    case 5:
      return 2l;
      break;
    case 2:
    case 3:
    case 4:
      return 4l;
      break;
    }	
  }


  return(0);
}

#ifdef __GTK_GUI
/* --getmri_paraversion() ---------------------------------------------------
**
**    Gets paravision's version 

**    Usage:    getmri_paraversion(file)
**        file: char*: Represents a file witch can contains path 
**
**    Return: This fct returns a long witch indicates readformat
**                     
*/
int getmri_paraversion(char *file)
{
  char *file_acqp;
  int version=0;
  char v[32];
  file_acqp = getfilename_acqp(file);

  if (file_acqp)
    {
      strcpy(v, getheader_interfile(file_acqp,"##$ACQ_sw_version=", TAB_PARA_STRING, 0));
      if (strstr(v,"<PV 2."))
	{// on est en paravision 2
	  version = 2;
	}
      else 
	{// on est en paravision 1
	  version = 1;
	}
      FREE(file_acqp);
    }
	
  return version;
}
#endif // __GTK_GUI

/* --getmri_readformat() ---------------------------------------------------
**
**    Gets type of readformat for 3D image 
**    Est utile pour connaitre comment afficher une image.
**    Attention retour utilise dans read_raw_3d 
**             et dans lecture suite de fichier 2D          
**
**    Usage:    getmri_readformat(file)
**        file: char*: Represents a file witch can contains path 
**
**    Return: This fct returns a long witch indicates readformat
**                     
*/
long    getmri_readformat(char *file)
{
  char *file_mri;
  char chaine[50],chaine1[50];
  long format;
  char *data=(char *)NULL;
  long typeHFS;
  float val1,val2,val3,val4,val5,val6;
  int ival1,ival2,ival3,ival4,ival5,ival6;
	
  format = is_mri(file,NULL);

  /*    For Brucker file format   */
  if(format == BRUKER_FORMAT)
    {
      return(0);
    }

#ifdef __GTK_GUI
  /*  For PARAVISION format   */
  if(format == PARAVISION_FORMAT)
    {
      file_mri=getfilename_imnd(file);
      strcpy(chaine, getheader_interfile(file_mri,
					 "##$IMND_slice_orient=", ASCII8B, 0));
      if (strcmp(chaine,"Sagittal_Head_Foot")==0) {
	printf("ATTENTION le format Sagittal a ete pris pour cette image \n");
	return(2l);
      }
      else if (strcmp(chaine,"Coronal_Head_Foot")==0) {
	long retVal=0;
	if (getmri_paraversion(file)>1)
	  {// on est en paravision 2
	    printf("ATTENTION format Coronal_Head_Foot Paravision 2. \n");
	    retVal=11;
	  } 
	else
	  {// pas en paravision 2
	    printf("ATTENTION format a rechercher pour cette image \n");
	    retVal=0;
	  }
	return(retVal);
      }
      else if (strcmp(chaine,"Transverse_Posterior_Anterior")==0) {
	return(6l);
      }
      else if (strcmp(chaine,"Arbitrary_Oblique")==0) {
	return(6l);
      }
      else if (strcmp(chaine,"Transverse_Left_Right")==0) {
	strcpy(chaine, getheader_interfile(file_mri,
					   "##$IMND_dummy_method=", ASCII8B, 0));
	strcpy(chaine1, getheader_interfile(file_mri,
					    "##$IMND_method=", ASCII8BNL, 0));
	if (strcmp(chaine,"GE3D")==0||strcmp(chaine1,"<GE3D>\012")==0) { 
	  /*Ancienne version epilepsie Rouffach et IPB */
	  file_mri=getfilename_reco(file);
	  strcpy(chaine, getheader_interfile(file_mri
					     , "##$RECO_transpose_dim="
					     ,ENTIER, 0));
	  if (atol(chaine)==0) {
	    /*format rouffach 1*/
	    return(3l);
	  }
	  else {
	    /*format rouffach 2*/
	    return(0l);
	  }
	}
	else /* Nouvelle version Transverse_Left_Right Mai 1998i
		On considere toujours commencer par le bas les coupes
		Pas de verification de ACQ_slice_offset*/
	  return(6l);
      }
      else {  /*   Valeur par defaut */
	/*format IPB par defaut*/
	printf("ATTENTION le format IPB (Transverse) a ete pris par defaut pour cette image \n");
	return(2l);
      }

      FREE(file_mri);
    }
#endif /* __GTK_GUI */

  /*  For Interfile format  */
  if(format == INTERFILE_FORMAT)
    {
      return(1);
    }

  /*    For raw4 file format  A FAIRE  */
  if(format == RAW4_FORMAT)
    {
      return(0);
    }

  /* For Siemens Vanderbilt format */
  if(format == VANDERBILT_FORMAT)
    {
      return(0);
    }

  /*   For IPB 2D & 3D file format */
  if(format == IPB_2D_FORMAT || format == IPB_3D_FORMAT)
    {
      return(_read_format_3d); /*  Pas tres juste */
    }

  /*    For pcx file format A FAIRE  */
  if(format == PCX_FORMAT)
    {
      return(0);
    }

  /* For SMIS file format */
  if(format == SMIS_FORMAT) 
    {
      int x_angle,y_angle,z_angle;
      char tmp[64];
      char v[16];
      sprintf(tmp,":IM_ORIENTATION");
      strcpy(v, getheader_interfile(file, tmp, ENTIER, 0));
      x_angle = atoi(v);

      sprintf(tmp,"%s %d",tmp,atoi(v));
      strcpy(v, getheader_interfile(file, tmp, ENTIER, 0));
      y_angle = atoi(v);

      sprintf(tmp,"%s %d",tmp,atoi(v));
      strcpy(v, getheader_interfile(file, tmp, ENTIER, 0));
      z_angle = atoi(v);
		
      //		printf("Orientation prise par defaut :0\n");
      return(0);
    }
	
  /* For ACRNEMA 2.0 or 1.0 format */
  if (format == ACRNEMA2L_FORMAT || format == ACRNEMA2B_FORMAT 
      || format==ACRNEMA1B_FORMAT || format==ACRNEMA1L_FORMAT) 
    {
      data=Acrnema_getGroup(file,0x18 , 0x5100, data, format);
      if (strstr(data, "HFS")) typeHFS=9;
      else typeHFS=9;
      free (data);
      return typeHFS;
    }

  /* for dicom 2D and 3D format */
  if(format == DICMBIMP_FORMAT || format == DICMLIMP_FORMAT ||
     format == DICMBEXP_FORMAT || format == DICMLEXP_FORMAT  ||
     format == DICM3DBIMP_FORMAT || format == DICM3DLIMP_FORMAT ||
     format == DICM3DBEXP_FORMAT || format == DICM3DLEXP_FORMAT)
    {
      char *data=NULL,*data1,*data2,*data3,*data4,*data5,*data6;

      /* Recherche par Image orientation patient */
      data=Acrnema_getGroup(file, 0x20, 0x37, data,format);
      if (!strlen(data)==0) /* Ce groupe existe */ 
	{
	  data1=(char*)strtok(data,"\\");
	  data2=(char*)strtok(NULL,"\\");
	  data3=(char*)strtok(NULL,"\\");
	  data4=(char*)strtok(NULL,"\\");
	  data5=(char*)strtok(NULL,"\\");
	  data6=(char*)strtok(NULL,"\\");
	  val1=atof(data1); ival1=(int)(val1 < 0 ? val1 - 0.5 : val1 + 0.5);
	  val2=atof(data2); ival2=(int)(val2 < 0 ? val2 - 0.5 : val2 + 0.5);
	  val3=atof(data3); ival3=(int)(val3 < 0 ? val3 - 0.5 : val3 + 0.5);
	  val4=atof(data4); ival4=(int)(val4 < 0 ? val4 - 0.5 : val4 + 0.5);
	  val5=atof(data5); ival5=(int)(val5 < 0 ? val5 - 0.5 : val5 + 0.5);
	  val6=atof(data6); ival6=(int)(val6 < 0 ? val6 - 0.5 : val6 + 0.5);
	  if(ival1==-1&&ival2==0&&ival3==0&&
	     ival4==0&&ival5==0&&ival6==-1) 
	    {
	      free(data);
	      printf("ATTENTION read_format dans le cas : \"-1/0/0/0/0/-1\"\n");
	      return(0);
	    } /* -1/0/0/0/0/-1 */
	  else if(ival1==0&&ival2==1&&ival3==0&&
		  ival4==0&&ival5==0&&ival6==-1) 
	    {  
	      free(data);
	      printf("ATTENTION read_format dans le cas : \"0/1/0/0/0/-1\"\n");
	      return(10);
	    }  /* 0/1/0/0/0/-1 */
	  else if(ival1==1&&ival2==0&&ival3==0&&
		  ival4==0&&ival5==0&&ival6==0) 
	    {  
	      free(data);
	      printf("ATTENTION read_format dans le cas : \"1/0/0/0/0/0\"\n");
	      return(9);
	    }  /* 1/0/0/0/0/0 */
	  else 
	    {
	      free(data); 
	      printf("ATTENTION read_format 1 pris par defaut\n");
	      return(9);
	    }
	}

      /* recherche par type d'examen */
      data=Acrnema_getGroup(file, 0x18, 0x24, data,format);
      if (!strlen(data)==0) /* Ce groupe existe */ 
	{
	  if (!strcmp(data,"mpr ")) 
	    {  free(data); return 10; } 
	  else if ( !strcmp(data,"se1 "))
	    {  free(data); 
	    data=Acrnema_getGroup(file, 0x8, 0x103e, data,format);
	    if( !strcmp(data,"IRM-FONCTIONNELLE/T2-SE-TE"))
	      { free(data); return 9; }
	    else 
	      { free(data);  return 10; }
	    }				 
	  else if ( !strcmp(data,"tirm1 "))
	    {  free(data); return 9;  }
	  else if ( !strcmp(data,"tse2-5"))
	    {  free(data); return 9;  }
	  else if ( !strcmp(data,"ep2d_z"))
	    {  free(data); return 9;  }
	  else /* default */
	    free(data); return 10;	
	}
      /* Recherche par Image position patient */
      data=Acrnema_getGroup(file, 0x20, 0x32, data,format);
      if (!strlen(data)==0) /* Ce groupe existe */ 
	{
	  data1=(char*)strtok(data,"\\");
	  data2=(char*)strtok(NULL,"\\");
	  data3=(char*)strtok(NULL,"\\");
	  if(atoi(data1)!=0) {free(data);return(4);} /* x/0/0 */
	  else if (atoi(data2)!=0) {free(data);return(0);} /* 0/x/0 */
	  else {free(data);return(1);}  /* 0/0/0 */
	}
      printf("ATTENTION read_format 2 pris par defaut\n");
      if (format == DICMBIMP_FORMAT || format == DICMLIMP_FORMAT ||
	  format == DICMBEXP_FORMAT || format == DICMLEXP_FORMAT) return 10;
      return 6; /* Par defaut il semblerait que 6 soit bien adapte, a verifier!!*/	  
    }
	
  /* For ECAT7 format */
  if(format==ECAT7_FORMAT){
    return 7;
  }

  return(0);
}

/* --getmri_name() ---------------------------------------------------
**
**      Return the name of the patient (paravision file only)
**
**
**      Usage:  getmri_name(file)
**              file: char*: Represents a file witch can contains path
**
**      Return: This fct returns a char*
**
*/
char*   getmri_name(char *file)
{
  char *v=(char*)malloc(50),*file_subject=(char *)NULL;
  char tmp[50];
  long format;

  format = is_mri(file,NULL);

  /*    For Brucker file format   */
  if(format == BRUKER_FORMAT)
    {
      free(v);
      v=NULL;
      return(v);
    }

#ifdef __GTK_GUI
  /*  For PARAVISION format   */
  if(format == PARAVISION_FORMAT)
    {
      file_subject=getfilename_subject(file);
      strcpy(tmp, getheader_interfile(file_subject, "##$SUBJECT_name_string=", ASCII8BNL, 0));
      FREE(file_subject);
      strncpy(v,tmp+1,strlen(tmp)-3);
      *(v+strlen(tmp)-3) = '\0';
      return(v);
    }
#endif /* __GTK_GUI */

  /*  For Interfile format   */
  if(format ==INTERFILE_FORMAT)
    {
      strcpy(v,getheader_interfile(file, "patient name:=", ASCII8B, 0));
      return((char*) v);
    }

  /*    For raw4 file format   */
  if(format == RAW4_FORMAT)
    {
      free(v);
      v=NULL;
      return(v);
    }

  /*   For Siemens Magnetom Vision format, Cyril Feruglio Mai 1998 */
  if(format == SIEMENS_FORMAT)
    {
      strcpy(v,getchar_header(file, 0x300, 25));
      return(v);
    }

  /* For SIEMENS Vanderbilt FORMAT Feruglio Cyril Mai 1998 */
  if(format == VANDERBILT_FORMAT)
    {
      strcpy(v, getheader_interfile(file, "Other Patient ID's :=", ASCII8B, 0));
      return(v);
    }

  /*   For IPB file formats */
  if(format == IPB_2D_FORMAT || format == IPB_3D_FORMAT)
    {
      strcpy(v, getheader_interfile(file, "patient name=", STRING, 0));
      return(v);
    }

  /*    For pcx file format   */
  if(format == PCX_FORMAT)
    {
      free(v);
      v=NULL;
      return(v);
    }

  /* For SMIS file format */
  if(format == SMIS_FORMAT) {
    free(v);
    v=NULL;
    return(v);
  }
	
  /* For Acrnema format   Burger Laurent   novembre 1999*/
  if (format == ACRNEMA2L_FORMAT || format==ACRNEMA2B_FORMAT 
      || format==ACRNEMA1B_FORMAT || format==ACRNEMA1L_FORMAT)
    {
      file_subject=Acrnema_getGroup(file, 0x10, 0x10, file_subject,format);	
      sprintf(v, "%s", file_subject);
      free (file_subject);
      return v;
    }
	
  /* for dicom format */
  if(format == DICMBIMP_FORMAT || format == DICMLIMP_FORMAT ||
     format == DICMBEXP_FORMAT || format == DICMLEXP_FORMAT)
    {
      file_subject=Acrnema_getGroup(file, 0x10, 0x10, file_subject,format);	
      sprintf(v, "%s", file_subject);
      free (file_subject);
      return v;
    }

  /* for dicom 3D format */
  if(format == DICM3DBIMP_FORMAT || format == DICM3DLIMP_FORMAT ||
     format == DICM3DBEXP_FORMAT || format == DICM3DLEXP_FORMAT)
    {
      file_subject=Acrnema_getGroup(file, 0x10, 0x10, file_subject,format);
      /* group = PatientsName */	
      sprintf(v, "%s", file_subject);
      free (file_subject);
      return v;
    }
  /* end of dicom 3D format */

  /* For ECAT7 format */
  if(format==ECAT7_FORMAT){
    strcpy(v,getchar_header(file, 182, 32));
    return(v);		
  }

  free(v);
  v=NULL;
  return(v);
}

/* --getmri_medecin() ---------------------------------------------------
**
**      Return the name of the doctor (paravision file only)
**
**
**      Usage:  getmri_medecin(file)
**              file: char*: Represents a file witch can contains path
**
**      Return: This fct returns a char*
**
*/
char*   getmri_medecin(char *file)
{
  char *v=(char*)malloc(50);
  long format;

  format = is_mri(file,NULL);

  /*    For Brucker file format   */
  if(format == BRUKER_FORMAT)
    {
      free(v);
      v=NULL;
      return(v);
    }

  /*  For PARAVISION format   */
  if(format == PARAVISION_FORMAT)
    {
      sprintf(v,"Dr Namer");
      return((char*)v);
    }

  /*  For Interfile format   */
  if(format == INTERFILE_FORMAT)
    {
      free(v);
      v=NULL;
      return(v);
    }

  /*    For raw4 file format   */
  if(format == RAW4_FORMAT)
    {
      free(v);
      v=NULL;
      return(v);
    }

  /*   For Siemens Magnetom Vision format, Cyril Feruglio Mai 1998 */
  if(format == SIEMENS_FORMAT)
    {
      strcpy(v,getchar_header(file, 213, 25));
      return(v);
    }

  /*   For IPB file formats */
  if(format == IPB_2D_FORMAT || format == IPB_3D_FORMAT)
    {
      strcpy(v, getheader_interfile(file, "doctor name=", STRING, 0));
      return(v);
    }

  /*    For pcx file format   */
  if(format == PCX_FORMAT)
    {
      free(v);
      v=NULL;
      return(v);
    }

  /* For SMIS file format */
  if(format == SMIS_FORMAT) {
    free(v);
    v=NULL;
    return(v);
  }

  /* For Acrnema format */
  if (format == ACRNEMA2L_FORMAT || format==ACRNEMA2B_FORMAT 
      || format==ACRNEMA1B_FORMAT || format==ACRNEMA1L_FORMAT)
    {
      char *medecin=(char *)NULL;

      medecin=Acrnema_getGroup(file, 0x8, 0x80, medecin,format);
      sprintf(v, "%s", medecin);
      free (medecin);
      return v;
    }

  /* for dicom format */
  if(format == DICMBIMP_FORMAT || format == DICMLIMP_FORMAT ||
     format == DICMBEXP_FORMAT || format == DICMLEXP_FORMAT)
    {
      char *medecin=(char *)NULL;

      medecin=Acrnema_getGroup(file, 0x8, 0x80, medecin,format);
      sprintf(v, "%s", medecin);
      free (medecin);
      return v;
    }

  /* for dicom 3D format */
  if(format == DICM3DBIMP_FORMAT || format == DICM3DLIMP_FORMAT ||
     format == DICM3DBEXP_FORMAT || format == DICM3DLEXP_FORMAT)
    {
      char *medecin=(char *)NULL;

      medecin=Acrnema_getGroup(file, 0x8, 0x80, medecin,format);
      /* group = InstitutionName */
      sprintf(v, "%s", medecin);
      free (medecin);
      return v;
    }
  /* end of dicom 3D format */

  /* For ECAT7 format */
  if(format==ECAT7_FORMAT){
    strcpy(v,getchar_header(file, 264, 32));
    return(v);		
  }


  free(v);
  v=NULL;
  return(v);
}

/* --getmri_dbirth() ---------------------------------------------------
**
**      Return the date of birth of the patient (paravision file only)
**
**
**      Usage:  getmri_dbirth(file)
**              file: char*: Represents a file witch can contains path
**
**      Return: This fct returns a char*
**
*/
char*   getmri_dbirth(char *file)
{
  char *v=(char*)malloc(50),*file_subject=(char *)NULL;
  char tmp[50];
  long format;

  format = is_mri(file,NULL);

  /*    For Brucker file format   */
  if(format == BRUKER_FORMAT)
    {
      free(v);
      v=NULL;
      return(v);
    }

#ifdef __GTK_GUI
  /*  For PARAVISION format   */
  if(format == PARAVISION_FORMAT)
    {
      file_subject=getfilename_subject(file);
      strcpy(tmp, getheader_interfile(file_subject, "##$SUBJECT_dbirth=", ASCII8BNL, 0));
      FREE(file_subject);
      strncpy(v,tmp+1,strlen(tmp)-3);
      *(v+strlen(tmp)-3) = '\0';
      return(v);
    }
#endif /* __GTK_GUI */

  /*  For Interfile format   */
  if(format == INTERFILE_FORMAT)
    {
      free(v);
      v=NULL;
      return(v);
    }

  /*    For raw4 file format   */
  if(format == RAW4_FORMAT)
    {
      free(v);
      v=NULL;
      return(v);
    }

  /*   For Siemens Magnetom Vision format, Cyril Feruglio Mai 1998 */
  if(format == SIEMENS_FORMAT)
    {
      strcpy(v,getchar_header(file, 5956, 12));
      return(v);
    }

  /*   For IPB file format */
  if(format == IPB_2D_FORMAT || format == IPB_3D_FORMAT)
    {
      strcpy(v, getheader_interfile(file, "date of birth=", STRING, 0));
      return(v);
    }

  /*    For pcx file format   */
  if(format == PCX_FORMAT)
    {
      free(v);
      v=NULL;
      return(v);
    }

  /* For SMIS file format */
  if(format == SMIS_FORMAT) {
    free(v);
    v=NULL;
    return(v);
  }

  /* For Acrnema format */
  if (format == ACRNEMA2L_FORMAT || format==ACRNEMA2B_FORMAT 
      || format==ACRNEMA1B_FORMAT || format==ACRNEMA1L_FORMAT)
    {
      file_subject=Acrnema_getGroup(file, 0x10, 0x30, file_subject,format);
      sprintf(v, "%s", file_subject);
      free (file_subject);
      return v;
    }

  /* for dicom format */
  if(format == DICMBIMP_FORMAT || format == DICMLIMP_FORMAT ||
     format == DICMBEXP_FORMAT || format == DICMLEXP_FORMAT)
    {
      file_subject=Acrnema_getGroup(file, 0x10, 0x30, file_subject,format);
      sprintf(v, "%s", file_subject);
      free (file_subject);
      return v;
    }

  /* for dicom 3D format */
  if(format == DICM3DBIMP_FORMAT || format == DICM3DLIMP_FORMAT ||
     format == DICM3DBEXP_FORMAT || format == DICM3DLEXP_FORMAT)
    {
      file_subject=Acrnema_getGroup(file, 0x10, 0x30, file_subject,format);
      /* group = PatientsBirthDate */
      sprintf(v, "%s", file_subject);
      free (file_subject);
      return v;
    }
  /* end of dicom 3d format */

  /* For ECAT 7 format */
  if(format == ECAT7_FORMAT){
    v[0]='\0';
    return v;
  }
  free(v);
  v=NULL;
  return(v);
}

/* --getmri_dexamen() ---------------------------------------------------
**
**      Return the date of the exam of the patient (paravision file only)
**
**
**      Usage:  getmri_dexamen(file)
**              file: char*: Represents a file witch can contains path
**
**      Return: This fct returns a char*
**
*/
char*   getmri_dexamen(char *file)
{
  char *v=(char*)malloc(50),*file_subject=(char *)NULL;
  char tmp[50];
  long format;

  format = is_mri(file,NULL);

  /*    For Brucker file format   */
  if(format == BRUKER_FORMAT)
    {
      free(v);
      v=NULL;
      return(v);
    }

#ifdef __GTK_GUI
  /*  For PARAVISION format   */
  if(format == PARAVISION_FORMAT)
    {
      file_subject=getfilename_subject(file);
      strcpy(tmp,getheader_interfile(file_subject,"##$SUBJECT_date=",ASCII8BNL,0));
      FREE(file_subject);
      strncpy(v,tmp+9,strlen(tmp)-11);
      *(v+strlen(tmp)-11) = '\0';
      return(v);
    }
#endif /* __GTK_GUI */

  /*  For Interfile format   */
  if(format == INTERFILE_FORMAT)
    {
      free(v);
      v=NULL;
      return(v);
    }

  /*    For raw4 file format   */
  if(format == RAW4_FORMAT)
    {
      free(v);
      v=NULL;
      return(v);
    }

  /*   For Siemens Magnetom Vision format, Cyril Feruglio Mai 1998 */
  if(format == SIEMENS_FORMAT)
    {
      strcpy(v,getchar_header(file, 5559, 12));
      return(v);
    }

  /*   For IPB file formats */
  if(format == IPB_2D_FORMAT || format == IPB_3D_FORMAT)
    {
      strcpy(v, getheader_interfile(file, "date of exam=", STRING, 0));
      return(v);
    }

  /*    For pcx file format   */
  if(format == PCX_FORMAT)
    {
      free(v);
      v=NULL;
      return(v);
    }

  /* For SMIS file format */
  if(format == SMIS_FORMAT) {
    free(v);
    v=NULL;
    return(v);
  }

  /* For Acrnema format */
  if (format == ACRNEMA2L_FORMAT || format==ACRNEMA2B_FORMAT 
      || format==ACRNEMA1B_FORMAT || format==ACRNEMA1L_FORMAT)
    {
      file_subject=Acrnema_getGroup(file, 0x8, 0x20, file_subject, format);
      sprintf(v, "%s", file_subject);
      free(file_subject);
      return v;
    }

  /* for dicom format */
  if(format == DICMBIMP_FORMAT || format == DICMLIMP_FORMAT ||
     format == DICMBEXP_FORMAT || format == DICMLEXP_FORMAT)
    {
      file_subject=Acrnema_getGroup(file, 0x8, 0x20, file_subject, format);
      sprintf(v, "%s", file_subject);
      free(file_subject);
      return v;
    }

  /* for dicom 3D format */
  if(format == DICM3DBIMP_FORMAT || format == DICM3DLIMP_FORMAT ||
     format == DICM3DBEXP_FORMAT || format == DICM3DLEXP_FORMAT)
    {
      file_subject=Acrnema_getGroup(file, 0x8, 0x20, file_subject, format);
      /* group = StudyDate */
      sprintf(v, "%s", file_subject);
      free(file_subject);
      return v;
    }
  /* end of dicom 3d format */

  /* For ECAT7 format */
  if(format==ECAT7_FORMAT){
    v[0]='\0';
    return v;
  }
	
  free(v);
  v=NULL;
  return(v);
}

#ifdef __GTK_GUI
/*
**  -- getmri_dspStyle() ------------------------------------------------------
**
**  Gets display style
**
**  Return value:
**  - display style of picture
**  - DSP_POSITIVE otherwise
**
**  >>  Created by Mr BOULEAU on Oct, 31st 2000
**
*/

t_CMapDisplayMode getmri_dspStyle(char *file, int pos, int mode)
{
  t_CMapDisplayMode res = DSP_POSITIVE;
  long format;

  format = is_mri(file, NULL);

  if (mode == MODE_2D && format == IPB_3D_FORMAT)
    pos = getmri_index3D(file, pos);

  if (format == IPB_3D_FORMAT || format == IPB_2D_FORMAT)
    {
      res = atoi(getheader_interfile(file, "display style=", ENTIER_IPB_2D,
				     pos));

      if(res == 0)
	res = DSP_POSITIVE;
    }
    
  return res;
}

/*
**  -- getmri_dspOverBlancking() ------------------------------------------------------
**
**  Gets display overblancking style
**
**  Return value:
**  - display style of picture
**  - OVERBLANCKING_NONE otherwise
**
**  >>  Created by Mr BOULEAU on Oct, 31st 2000
**
*/
t_CMapOverBlanckingMode getmri_dspOverBlancking(char *file, int pos, int mode)
{
  t_CMapOverBlanckingMode res = OVERBLANCKING_NONE;
  long format;

  format = is_mri(file, NULL);

  if (mode == MODE_2D && format == IPB_3D_FORMAT)
    pos = getmri_index3D(file, pos);

  if (format == IPB_3D_FORMAT || format == IPB_2D_FORMAT)
    {
      res = atoi(getheader_interfile(file, "overblancking=", ENTIER_IPB_2D,
				     pos));

      if(res == 0)
	res = OVERBLANCKING_NONE;
    }

  return res;
}





#endif /* __GTK_GUI */

/*
**  -- getmri_nCMap() ---------------------------------------------------------
**
**  Gets first colormap id of pictures. If colormap doesn't exist, creates it.
**
**  Return value:
**  - id of colormap to create
**  - 0 otherwise
**
**  >>  Created by Mr BOULEAU on Oct, 31st 2000
**
*/

int getmri_nCMap(char *file, int pos, int mode)
{
  int res = 1;
  long format;

  format = is_mri(file, NULL);

  if(mode == MODE_2D && format == IPB_3D_FORMAT)
    pos = getmri_index3D(file, pos);

  if(format == IPB_3D_FORMAT || format == IPB_2D_FORMAT)
    res = atol(getheader_interfile(file, "colormap=", ENTIER_IPB_2D, pos));
    
  return res;
}

#ifdef __GTK_GUI
/* --getfilename_d3proc() ------------------------------------------------    
**
**    Gets filename of d3proc in structure PARAVISION
**
**    Usage:    getfilename_d3proc(paravision_subjeuct_file)
**        paravision_subjeuct_file: char*: Represents a file.
**
**    Return: This fct returns a char with indicates filename.
**                     
*/
char    *getfilename_d3proc(char *file)
     /* Assumed: Bruker file */
{
  FILE *fp;
  char *filename;
  int k;

  /*    For Brucker file format   */
  if((is_mri(file,NULL))==PARAVISION_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	fclose(fp);
	return((char *) NULL);
      }

      filename=CALLOC(FILE_LEN,char);
      for(k=strlen(file);file[k]!='/';k--)  ;
      strncpy(filename,file,++k);
      strcat(filename,_path_paravision);
      strcat(filename,"/pdata/1/d3proc");
      fclose(fp);
      return((char *) filename);
    }

  return((char *) NULL);
}

/* --getfilename_procs() ------------------------------------------------    
**
**    Gets filename of procs in structure PARAVISION
**
**    Usage:    getfilename_procs(paravision_subjeuct_file)
**        paravision_subjeuct_file: char*: Represents a file.
**
**    Return: This fct returns a char with indicates filename.
**                     
*/
char    *getfilename_procs(char *file)
     /* Assumed: Bruker file */
{
  FILE *fp;
  char *filename;
  int k;
  struct stat status;

  /*    For Brucker file format   */
  if((is_mri(file,NULL))==PARAVISION_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	fclose(fp);
	return((char *) NULL);
      }

      filename=CALLOC(FILE_LEN,char);
      for(k=strlen(file);file[k]!='/';k--)  ;
      strncpy(filename,file,++k);
      strcat(filename,_path_paravision);
      strcat(filename,"/pdata/1/procs");    
      if(stat(filename, &status) == -1) 
	{
	  /* file doesn't exists */
	  strncpy(filename,file,k);
	  filename[k]='\0';
	  strcat(filename,_path_paravision);
	  strcat(filename,"/pdata/1/d3proc");    
	}

      fclose(fp);
      return((char *) filename);
    }

  return((char *) NULL);
}

/* --getfilename_2dseq() ------------------------------------------------    
**
**    Gets filename of 2dseq in structure PARAVISION
**
**    Usage:    getfilename_2dseq(paravision_subjeuct_file)
**        paravision_subjeuct_file: char*: Represents a file.
**
**    Return: This fct returns a char with indicates filename.
**                     
*/
char    *getfilename_2dseq(char *file)
     /* Assumed: Bruker file */
{
  FILE *fp;
  char *filename;
  int k;

  /*    For Brucker PARAVISION file format   */
  if((is_mri(file,NULL))==PARAVISION_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	fclose(fp);
	return((char *) NULL);
      }

      filename=CALLOC(FILE_LEN,char);
      for(k=strlen(file);file[k]!='/';k--)  ;
      strncpy(filename,file,++k);
      strcat(filename,_path_paravision);
      strcat(filename,"/pdata/1/2dseq");
      fclose(fp);
      return((char *) filename);
    }

  return((char *) NULL);
}

#endif /* __GTK_GUI */

/* --getfilename_img() ------------------------------------------------    
**
**    Gets filename of img in structure INTERFILE, IPB
** and Siemens Vanderbilt
**    Usage:    getfilename_2dseq(paravision_subjeuct_file)
**        paravision_subjeuct_file: char*: Represents a file.
**
**    Return: This fct returns a char with indicates filename.
**                     
*/
char    *getfilename_img(char *file)
     /* Assumed: Bruker file */
{
  FILE *fp;
  char *filename;

  /*    For INTERFILE file format  */
  if((is_mri(file,NULL))==INTERFILE_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	fclose(fp);
	return((char *) NULL);
      }
      filename=CALLOC(strlen(file),char);
      strcpy(filename,file);
      strcpy(filename+strlen(file)-3,"IMG");
      fclose(fp);
      return((char *) filename);
    }

  /*    For IPB file format   */
  if((is_mri(file,NULL))==IPB_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	fclose(fp);
	return((char *) NULL);
      }

      filename=CALLOC(strlen(file),char);
      strcpy(filename,file);
      strcpy(filename+strlen(file)-3,"img");
      fclose(fp);
      return((char *) filename);
    }


  /*   For IPB 2D file format. Ne prend que la premiere valeur si plusieurs images dans le fichier */
  if((is_mri(file,NULL))==IPB_2D_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	fclose(fp);
	return((char *) NULL);
      }

      filename=CALLOC(strlen(file),char);
      strcpy(filename,file);
      strcpy(filename+strlen(file)-3,"img");
      fclose(fp);
      return((char *) filename);
    }


  /*    For Siemens vanderbilt file format   */
  if((is_mri(file,NULL))==VANDERBILT_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	_imagix_error=IMX_FILE_NOT_FOUND;        /* Open error    */
	_imagix_err_ext=errno;
	fclose(fp);
	return((char *) NULL);
      }

      filename=CALLOC(strlen(file),char);
      strcpy(filename,file);
      strcpy(filename+strlen(file)-12,"image.bin");
      if(MESG_DEBUG) printf("file:%s",filename);
      fclose(fp);
      return((char *) filename);
    }


  return((char *) NULL);
}


/* --getfilename_reco() ------------------------------------------------    
**
**    Gets filename of reco in structure PARAVISION
**
**    Usage:    getfilename_2dseq(paravision_subjeuct_file)
**        paravision_subjeuct_file: char*: Represents a filegetfilename_reco.
**
**    Return: This fct returns a char with indicates filename.
**                     
*/
char    *getfilename_reco(char *file)
     /* Assumed: Bruker file */
{
  FILE *fp;
  char *filename;
  int k;

  /*    For Brucker PARAVISION file format   */
  if((is_mri(file,NULL))==PARAVISION_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	fclose(fp);
	return((char *) NULL);
      }

      filename=CALLOC(FILE_LEN,char);
      for(k=strlen(file);file[k]!='/';k--)  ;
      strncpy(filename,file,++k);
      strcat(filename,_path_paravision);
      strcat(filename,"/pdata/1/reco");
      fclose(fp);
      return((char *) filename);
    }

  return((char *) NULL);
}

/* --getfilename_imnd() ------------------------------------------------    
**
**    Gets filename of imnd in structure PARAVISION
**
**    Usage:    getfilename_2dseq(paravision_subjeuct_file)
**        paravision_subjeuct_file: char*: Represents a file.
**
**    Return: This fct returns a char with indicates filename.
**                     
*/
char    *getfilename_imnd(char *file)
     /* Assumed: Bruker file */
{
  FILE *fp;
  char *filename;
  int k;

  /*    For Brucker PARAVISION file format   */
  if((is_mri(file,NULL))==PARAVISION_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	fclose(fp);
	return((char *) NULL);
      }

      filename=CALLOC(FILE_LEN,char);
      for(k=strlen(file);file[k]!='/';k--)  ;
      strncpy(filename,file,++k);
      strcat(filename,_path_paravision);
      strcat(filename,"/imnd");
      fclose(fp);
      return((char *) filename);
    }

  return((char *) NULL);
}

/* --getfilename_acqp() ------------------------------------------------    
**
**    Gets filename of acqp in structure PARAVISION
**
**    Usage:    getfilename_acqp(paravision_subjeuct_file)
**        paravision_subjeuct_file: char*: Represents a file.
**
**    Return: This fct returns a char with indicates filename.
**                     
*/
char    *getfilename_acqp(char *file)
     /* Assumed: Bruker file */
{
  FILE *fp;
  char *filename;
  int k;

  /*    For Brucker PARAVISION file format   */
  if((is_mri(file,NULL))==PARAVISION_FORMAT)
    {
      if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
	fclose(fp);
	return((char *) NULL);
      }

      filename=CALLOC(FILE_LEN,char);
      for(k=strlen(file);file[k]!='/';k--)  ;
      strncpy(filename,file,++k);
      strcat(filename,_path_paravision);
      strcat(filename,"/acqp");
      fclose(fp);
      return((char *) filename);
    }

  return((char *) NULL);
}

/* --getfilename_subject() ------------------------------------------------
**
**      Gets filename of subject in structure PARAVISION
**
**      Usage:  getfilename_subject(paravision_subject_file)
**              paravision_subject_file: char*: Represents a file.
**
**      Return: This fct returns a char with indicates filename.
**
*/
char    *getfilename_subject(char *file)
     /* Assumed: Bruker file */
{
  FILE *fp;
  char *filename;
  int k;

  /*    For Brucker file format   */
  /*  if((is_mri(file,NULL))==PARAVISION_FORMAT)         DOMTMP */
  {
    if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
      fclose(fp);
      return((char *) NULL);
    }

    filename=CALLOC(FILE_LEN,char);
    for(k=strlen(file);file[k]!='/';k--)  ;
    strncpy(filename,file,++k);
    strcat(filename,"/subject");
    fclose(fp);
    return((char *) filename);
  }
  return((char *) NULL);
}

/* --getheader_interfile() -------------------------------------------------*/
/*!
**   \brief  Gets information stored in the header of the interfile
**      format
**
**    Usage:    getheader_interfile(file,stringbis,type,number)
**    \param  file: char*: file name 
**    \param  stringbis : search line
**    \param  type : type of data to read
**      - ENTIER : recupere un int derriere le string 
**      - REEL   : recupere un float derriere le string 
**      - TAB_PARA_RELLE   :  reel ligne suivant pour paravision  
**      - TAB_PARA_ENTIER   :  entier ligne suivant pour paravision  
**      - TAB_PARA_STRING   :  string ligne suivant pour paravision  
**      - TAB_PARA   :  
**      - ASCII8B   : recupere un *char derriere le string 
**      - ENTIER_IPB   : recupere un int indexe par numero derriere le string 
**                     expl : width[numero]:int
**                            -----   <-- string (width)
**      -REEL_IPB  :recupere un float indexe par numero derriere le string
**
**    \param  number : dans le cas d'un vecteur???
**
**
**    \retval  This fct returns a char recupered, NULL en cas d'erreur
**                     
*/
char     *getheader_interfile(char *file, char *stringbis, int type, int number)
{
  char *fmt1,*fmt2;
  FILE *fp;
  int i,k,d;
  float f=0.0;
  static char s[8192], string[100];  

  if((fp=fopen(file,"rb"))==NULL) { 
    /* printf("File open error: %s \n", file);*/
    return(0l);    /*Open file error... */
  }
  fmt1=CALLOC(8192,char);
  fmt2=CALLOC(8192,char);

  *s = '\0';

  if (type==ENTIER_IPB || type==REEL_IPB || (type==STRING && number)) {
    sprintf(string,"%s[%d]=",stringbis,number);
  }
  else {
    strcpy(string,stringbis);
  }
  while ((fgets(fmt1,8192,fp))!=NULL)
    {
      if (strstr(fmt1,string))
        {
	  switch(type) {
	  case ENTIER:
	    sprintf(fmt2,"%s %%d",string);
	    sscanf(fmt1,fmt2,&i);
	    sprintf(s,"%d",i);
	    break;
	  case REEL:
	    sprintf(fmt2,"%s %%f",string);
	    sscanf(fmt1,fmt2,&f);
	    sprintf(s,"%f",f);
	    break;
	  case TAB_PARA:
	    sprintf(fmt2,"%s %%d",string);
	    sscanf(fmt1,fmt2,&i);
	    fgets(fmt1,8192,fp);
	    if (number==0) {  /* 1 element de la suite */
	      strcpy(fmt2,"%d");
	      sscanf(fmt1,fmt2,&i);
	      sprintf(s,"%d",i);
	      break;
	    }
	    strcpy(fmt2,"%*d"); /*   sinon  */
	    for (k=1;k<i;k++) {
	      if (number==k) {
		strcat(fmt2,"%d");
		sscanf(fmt1,fmt2,&i);
		sprintf(s,"%d",i);
		break;
	      }
	      else
		strcat(fmt2,"%*d");
	    }
	    if(k==i) sprintf(s,"%d",0);
	    break;
	  case TAB_PARA_REEL:
	    sprintf(fmt2,"%s %%d",string);
	    sscanf(fmt1,fmt2,&i);
	    fgets(fmt1,8192,fp);
	    if (number==0) {  /* 1 element de la suite */
	      strcpy(fmt2,"%f");
	      sscanf(fmt1,fmt2,&f);
	      sprintf(s,"%f",f);
	      break;
	    }
	    strcpy(fmt2,"%*f"); /*   sinon  */
	    for (k=1;k<i;k++) {
	      if (number==k) {
		strcat(fmt2,"%f");
		sscanf(fmt1,fmt2,&f);
		sprintf(s,"%f",f);
		break;
	      }
	      else
		strcat(fmt2,"%*f");
	    }
	    if(k==i) sprintf(s,"%f",0.);
	    break;
	  case TAB_PARA_ENTIER:
	    sprintf(fmt2,"%s %%d",string);
	    sscanf(fmt1,fmt2,&i);
	    fgets(fmt1,8192,fp);
	    if (number==0) {  /* 1 element de la suite */
	      strcpy(fmt2,"%d");
	      sscanf(fmt1,fmt2,&d);
	      sprintf(s,"%d",d);
	      break;
	    }
	    strcpy(fmt2,"%*d"); /*   sinon  */
	    for (k=1;k<i;k++) {
	      if (number==k) {
		strcat(fmt2,"%d");
		sscanf(fmt1,fmt2,&d);
		sprintf(s,"%d",d);
		break;
	      }
	      else
		strcat(fmt2,"%*d");
	    }
	    if(k==i) sprintf(s,"%d",0);
	    break;
	  case TAB_PARA_STRING:
	    fgets(fmt1,8192,fp);
	    sprintf(s,"%s",fmt1);
	    break;
	  case ASCII8B:
	    sprintf(fmt2,"%s %%s",string);
	    sscanf(fmt1,fmt2,s);
	    break;
	  case STRING:
	    sprintf(fmt2,"%s%%[^\n]",string);
	    sscanf(fmt1,fmt2,s);
	    break;
	  case ASCII8BNL:
	    if (fgets(fmt1,8192,fp) == NULL) {
	      fclose(fp);
	      return((char*)NULL);
	    }
	    strcpy(s,fmt1);
	    break;
	  case ENTIER_IPB:  /* lit le premier entier qui suit la chaine */
	    sprintf(fmt2,"%s %%d",string/*,number*/);
	    sscanf(fmt1,fmt2,&i);
	    sprintf(s,"%d",i);
	    break;
	  case ENTIER_IPB_2D:  /* lit un entier a une position donnee */
	    sprintf(fmt2,"%s",string);
	    for (k=1;k<number;k++) {
	      strcat(fmt2,"%*d");
	    }
	    strcat(fmt2,"%d");
	    sscanf(fmt1,fmt2,&i);
	    sprintf(s,"%d",i);
	    break;
	  case REEL_IPB:  /* lit le premier reel qui suit la chaine */
	    sprintf(fmt2,"%s %%f",string/*,(float)number*/);
	    sscanf(fmt1,fmt2,&f);
	    sprintf(s,"%f",f);
	    break;
	  case REEL_IPB_2D:  /* lit un reel a une position donnee */
	    sprintf(fmt2,"%s",string);
	    for (k=1;k<number;k++) {
	      strcat(fmt2,"%*f");
	    }
	    strcat(fmt2,"%f");
	    sscanf(fmt1,fmt2,&f);
	    sprintf(s,"%f",f);
	    break;
	  default: 
	    if(MESG_DEBUG) printf("Default case... */");
	  }

	  FREE(fmt1);
	  FREE(fmt2);
	  fclose(fp);
	  return((char*)s);
        }
    }

  FREE(fmt1);
  FREE(fmt2);
  fclose(fp);
  *s=0;
  return(s);
}

/* --load_header() -------------------------------------------------------------
**
**    Loads and split in lines a header
**
**    Usage:    load_header(file, header)
**        file: char*: file name 
**              header : pointer to struct HEADER
**
**    Return: 1 - succesful
**              0 - error
**                     
*/
int     load_header(char *file, HEADER *header)
{
  struct stat status;
  FILE *fp;
  long size;
  char *buffer;
  int maxlines;
  char **lines;
  int i, j;
  int state;

  if(stat(file, &status) == -1) { /* file doesn't exists */
    lines = (char **)malloc(1024 * sizeof(char *));
    if(lines == NULL) {
      PUT_ERROR("[load_header] memory allocation error (1)!\n");
      return(0);
    }

    header->exists = 0;
    header->maxlines = 1024;
    header->numlines = 0;
    header->lines = lines;
  }
  else { /* file exists */
    if((fp = fopen(file, "rb")) == NULL) {
      PUT_ERROR("[load_header] cannot open file \n");
      return(0);
    }

    /* gets file size */
    size = status.st_size;

    if((buffer = (char *)malloc(size+1)) == NULL) {
      fclose(fp);
      PUT_ERROR("[load_header] memory allocation error (2)!\n");
      return(0);
    }

    /* loads header */
    if(fread(buffer, sizeof(char), size, fp) != size) {
      fclose(fp);
      FREE(buffer);
      PUT_ERROR("[load_header] cannot read file \n" );
      return(0);
    }

    buffer[size] = '\0';

    fclose(fp);

    /* allocates ((size + 1)/2 + 1) char* so we are sure that we can add a new line */
    maxlines = (size + 1)/2 + 1;

    /*protoize:???*/ /* verifier si chaque element de lines est bien alloue */
    if((lines = (char **)malloc(maxlines * sizeof(char *))) == NULL) {
      FREE(buffer);
      PUT_ERROR("[load_header] memory allocation error (3)!\n");
      return(0);
    }

    /* split buffer in lines */
    state = 0;

    for(i=0, j=0; i<size; i++) {
      if(buffer[i] == '\n') {
	buffer[i] = '\0';
	state = 0;
      }
      else {
	if(state == 0) {
	  lines[j++] = buffer + i;
	  state = 1;
	}
      }
    }

    header->exists = 1;
    header->maxlines = maxlines;
    header->numlines = j;
    header->lines = lines;
  }

  return(1);
}

/* --put_header() --------------------------------------------------------------
 */
/*!    Puts a info in a header
**
**    Usage:    put_header(header, string, type, itemadr, pos)
**    \param    header : pointer to struct HEADER
**    \param    string : string to match
**    \param	type : type of data (ENTIER, REEL, STRING)
**    \param    itemadr : address of data item
**    \param    pos : position in the line where to add the info
**
**    \retval  1 successful, 0  error
**              
**                     
*/
int     put_header(HEADER *header, char *string, int type, void *itemadr, int pos)
{
  char *newline;
  int maxlines, numlines;
  char **lines;
  int i;
  char *p = NULL;
  int state;
  int numwords;
  int length;
  char value[64];

  maxlines = header->maxlines;
  numlines = header->numlines;
  lines = header->lines;

  if(lines == NULL) {
    PUT_ERROR("[put_header] Invalid header !\n");
    return(0);
  }

  if((newline = (char *)malloc(256 * sizeof(char))) == NULL) {
    PUT_ERROR("[put_header] memory allocation error (1)!\n");
    return(0);
  }

  if(numlines > 0) { /* Search the line */
    for(i=0; i<numlines; i++) {
      if(strstr(lines[i], string) != NULL)
	{
	  /* we've found the string */
	  break;
	}
    }
  }
  else {
    i = numlines = 0;
  }

  /* Warning: don't modify 'i' it's used later */

  if(i == numlines) { /* case where the string was not found */
    strcpy(newline, string);

    numwords = 1;
    while(numwords++ < pos)
      strcat(newline, "0 ");
  }
  else {
    p = lines[i];

    if(strlen(p) > 192) {
      if((newline = (char *)realloc(newline, strlen(p) + 256)) == NULL) {
	PUT_ERROR("[put_header] reallocation failed !\n");
	return(0);
      }
    }   

    /* searches the '=' sign */
    while(*p != '\0' && *p != '=')
      p ++;

    if(!*p) { /* we have reached the end of the string */
      /* Supprime car erreur en ecrivant IPB3 
	 printf( "[put_header] no '=' in string\n"); */
      return(0);
    }

    state = 0;
    numwords = 0;

    /* skip the '=' sign */
    p ++;

    while(numwords < pos && *p != '\0') {
      if(*p == ' ')
	{
	  state = 0;
	  p++;
	}
      else {
	if(state == 0)
	  {
	    state = 1;
	    numwords ++;
	  }
	else
	  p ++;
      }

    }

    length = p - lines[i];

    strncpy(newline, lines[i], length);
    newline[length] = '\0';

    /* if there's not enough values, fill with zeroes */
    while(numwords++ < pos - 1)
      strcat(newline, "0 ");
  }

  switch(type)
    {
    case STRING:
      if(strlen((char *)itemadr) < 256)
	sprintf(value, "%s ", (char *)itemadr);
      else
	{
	  PUT_WARN("put_header : STRING too long!!!\n");
	  value[0] = '\0';
	}
      break;

    case ENTIER:
      //marat	sprintf(value, "%d ", *(long *)itemadr);
      sprintf(value, "%ld ", *(long *)itemadr);
        
      break;

    case REEL:
      sprintf(value, "%f ", *(double *)itemadr);
      break;

    default:
      return(0);
    }

  /* adds the value to the string */
  strcat(newline, value);

  /* we have found the line */
  if(i < numlines) {
    /* skip the value */
    if (pos!= 0) {
      while(*p != '\0' && *p != ' ')
	p ++;
      while(*p != '\0' && *p == ' ')
	p ++;

      /* copies the rest of the line */
      if(*p)
	strcat(newline, p);
    }
  }
  /* we must add a new line */
  else {
    /* verifier qu'on ne depasse pas header->maxlines !!! */
    if(numlines >= maxlines) {
      maxlines += 1024;
      /*protoize:???*/ /* idem*/
      if((lines = (char **)/*???*/realloc(lines, maxlines * sizeof(char *))) == NULL) {
	PUT_ERROR("[put_header] memory reallocation error(1)!\n");
	return(0);
      }
    }
    numlines ++;
  }

  /* replaces or add the line */
  lines[i] = newline;

  /* updates header structure */
  header->maxlines = maxlines;
  header->numlines = numlines;
  header->lines = lines;

  return(1);
}

/* --save_header() -------------------------------------------------------------
 */
/*!    Saves a header
**
**    Usage:    save_header(header, file)
**    \param    header : pointer to struct HEADER
**    \param    file : header file
**  
**    \retval 1 succesful, 0  error
**                     
*/
int     save_header(HEADER *header, char *file)
{
  char fileold[256];
  char filenew[256];
  FILE *fp;
  int i;

  if(header->lines == NULL) {
    PUT_ERROR("[save_header] Invalid header !\n");
    return(0);
  }

  strcpy(filenew, file);

  if(header->exists) { /* file exists */
    strcat(filenew, ".tmp");
  }

  /* writes the new header */
  if((fp = fopen(filenew, "w")) == NULL) {
    char e[256];
    sprintf(e, "[save_header] cannot create file \"%s\"\n", filenew);
    PUT_ERROR(e);
    return(0);
  }

  for(i=0; i<header->numlines; i++) {
    fprintf(fp, "%s\n", header->lines[i]);
  }

  fclose(fp);

  if(header->exists) {
    strcpy(fileold, file);
    /* strcat(fileold, ".ipb"); */

    if(remove(fileold) == -1) {
      char e[256];
      sprintf(e, "[save_header] unable to remove file \"%s\"\n", fileold);
      PUT_ERROR(e);
      return(0);
    }
    if(rename(filenew, fileold) == -1) {
      char e[256];
      sprintf(e, "[save_header] unable to rename file \"%s\" in \"%s\"\n", filenew, fileold);
      PUT_ERROR(e);
      return(0);
    }
  }

  header->exists = 1;
  header->maxlines = 0;
  header->numlines = 0;
  FREE(header->lines);
  header->lines = NULL;

  return(1);
}

/* --putheader_interfile() --------------------------------------------------
**
**    Puts information in the header of the file
**
**    Usage:    putheader_interfile(file, str, type, number)
**        file: char*: file name 
**              string : search line
**              type : type of data to read
**              number : 
**
**    Return: 1 - succesful
**              0 - error
**                     
*/
long    putheader_interfile(char *file, char *string, int type, void *itemadr, int pos)
{
  printf("[putheader_interfile] function no longer supported !\n\
Please use functions \"load_header/put_header/save_header\" !\n");

  return(0);
}

/* --putmri_icomp() ------------------------------------------------    
**    Puts icomp in a bruker picture
**
**    Usage:    putmri_icomp(src_file,dst_file,icomp)
**        src_file: char*: Represents a file. (May contains path )
**        dst_file: char*: Represents a file. (May contains path )
**        icomp: long: New icomp value.
**
**    Return: This fct returns a long.
**        0 when error occured.
*/
long    putmri_icomp(char *src_file, char *dst_file, long int icomp)
     /* Assumed: Bruker file */
     /* Assumed: Bruker file */
           
{ 
  char header[HEADER_SIZE],v[16];
  int ret,fd;
  FILE *src_fp;

  if((src_fp=fopen(src_file,"rb"))==NULL) { /* Open file error... */
    fclose(src_fp);
    return(0l);
  }
  fseek(src_fp,0l,SEEK_SET);
  fread(header,sizeof(char),sizeof(header),src_fp);
  fclose(src_fp);

  /*protoize:???*/ /*confusion open fopen??*/
#ifdef WIN32
  if((fd=_open(dst_file,_O_RDWR))==0) { /* Open file error... */
    close(fd);
    return(0l);
  }
#else
  if((fd=/*f*/open(dst_file,O_RDWR))==0) { /* Open file error... */
    close(fd);
    return(0l);
  }
#endif
  sprintf(v,"%ld",icomp);
  ret=codify_pos((unsigned char *)header,MVAR_OFF+3*(2+7),ENTIER,v);
  lseek(fd,0l,SEEK_SET);
  write(fd,header,sizeof(header));
  close(fd);

  return(ret);
}

/* --putmri_maxpix() ------------------------------------------------    
**    Puts maxpixel in a bruker picture
**
**    Usage:    putmri_maxpix(src_file,dst_file,maxpix)
**        src_file: char*: Represents a file. (May contains path )
**        dst_file: char*: Represents a file. (May contains path )
**        maxpix:long: Maxpixel.
**
**    Return: This fct returns a long.
**        0 when error occured.
*/
long    putmri_maxpix(char *src_file, char *dst_file, long int maxpix)
     /* Assumed: Bruker file */
     /* Assumed: Bruker file */
            
{ 
  char header[HEADER_SIZE],v[16];
  int ret,fd;
  FILE *src_fp;

  if((src_fp=fopen(src_file,"rb"))==NULL) { /* Open file error... */
    fclose(src_fp);
    return(0l);
  }
  fseek(src_fp,0l,SEEK_SET);
  fread(header,sizeof(char),sizeof(header),src_fp);
  fclose(src_fp);

  /*protoize:???*/ /*confusion oen fopen?*/
  if((fd=/*f*/open(dst_file,O_RDWR))==0) { /* Open file error... */
    close(fd);
    return(0l);
  }
  sprintf(v,"%ld",maxpix);
  ret=codify_pos((unsigned char *)header,MVAR_OFF+3*(2+5),ENTIER,v);
  lseek(fd,0l,SEEK_SET);
  write(fd,header,sizeof(header));
  close(fd);

  return(ret);
}

/* --putmri_minpix() ------------------------------------------------    
**    Puts minpixel in a bruker picture
**
**    Usage:    putmri_maxpix(src_file,dst_file,minpix)
**        src_file: char*: Represents a file. (May contains path )
**        dst_file: char*: Represents a file. (May contains path )
**        minpix:long: Maxpixel.
**
**    Return: This fct returns a long.
**        0 when error occured.
*/
long    putmri_minpix(char *src_file, char *dst_file, long int minpix)
     /* Assumed: Bruker file */
     /* Assumed: Bruker file */
             
{ 
  char header[HEADER_SIZE],v[16];
  int ret,fd;
  FILE *src_fp;

  if((src_fp=fopen(src_file,"rb"))==NULL) { /* Open file error... */
    fclose(src_fp);
    return(0l);
  }
  fseek(src_fp,0l,SEEK_SET);
  fread(header,sizeof(char),sizeof(header),src_fp);
  fclose(src_fp);

  /*protoize:???*/ /*confusion open fopen??*/
  if((fd=/*f*/open(dst_file,O_RDWR))==0) { /* Open file error... */
    close(fd);
    return(0l);
  }
  sprintf(v,"%ld",minpix);
  ret=codify_pos((unsigned char *)header,MVAR_OFF+3*(2+6),ENTIER,v);
  lseek(fd,0l,SEEK_SET);
  write(fd,header,sizeof(header));
  close(fd);

  return(ret);
}

/* --putmri_width() ------------------------------------------------    
**
**    Puts width in a bruker picture
**
**    Usage:    putmri_width(src_file,dst_file,width)
**        src_file: char*: Represents a file. (May contains path )
**        dst_file: char*: Represents a file. (May contains path )
**        width:long:  Width.
**
**    Return: This fct returns a long.
**        0 when error occured.
*/
long    putmri_width(char *src_file, char *dst_file, long int width)
     /* Assumed: Bruker file */
     /* Assumed: Bruker file */
            
{ 
  char header[HEADER_SIZE],v[16];
  int ret,fd;
  FILE *src_fp;

  if((src_fp=fopen(src_file,"rb"))==NULL) { /* Open file error... */
    fclose(src_fp);
    return(0l);
  }
  fseek(src_fp,0l,SEEK_SET);
  fread(header,sizeof(char),sizeof(header),src_fp);
  fclose(src_fp);

  if((fd=open(dst_file,O_RDWR))==0) { /* Open file error... */
    close(fd);
    return(0l);
  }
  sprintf(v,"%ld",width);
  ret=codify_pos((unsigned char *)header,MVAR_OFF+3*(2+4),ENTIER,v);
  lseek(fd,0l,SEEK_SET);
  write(fd,header,sizeof(header));
  close(fd);

  return(ret);
}

/* --putmri_height() -----------------------------------------------    
**
**    Puts height in a bruker picture
**
**    Usage:    putmri_height(src_file,dst_file,height)
**        src_file: char*: Represents a file. (May contains path )
**        dst_file: char*: Represents a file. (May contains path )
**        height:long:  Height.
**
**    Return: This fct returns a long.
**        0 when error occured.
*/
long    putmri_height(char *src_file, char *dst_file, long int height)
     /* Assumed: Bruker file */
     /* Assumed: Bruker file */
             
{ 
  char header[HEADER_SIZE],v[16];
  int ret,fd;
  FILE *src_fp;

  if((src_fp=fopen(src_file,"rb"))==NULL) { /* Open file error... */
    fclose(src_fp);
    return(0l);
  }
  fseek(src_fp,0l,SEEK_SET);
  fread(header,sizeof(char),sizeof(header),src_fp);
  fclose(src_fp);

  if((fd=open(dst_file,O_RDWR))==0) { /* Open file error... */
    close(fd);
    return(0l);
  }
  sprintf(v,"%ld",height);
  ret=codify_pos((unsigned char *)header,MVAR_OFF+3*(2+3),ENTIER,v);
  lseek(fd,0l,SEEK_SET);
  write(fd,&header,sizeof(header));
  close(fd);

  return(ret);
}

/* --buff2pic() -------------------------------------------------------    
**    Transforms a buffer to a picture[MAX_WIDTH][MAX_HEIGHT].
**
*/
long    buff2pic(long int *buff, long int **pic, grphic *image)
     /* long pic[][MAX_HEIGHT]; */
              
           
               
{ 
  long i,j,k;
  int w_step,h_step,step;
  char ss[100];
	
  w_step = (image->width + MAX_WIDTH - 1) / MAX_WIDTH;
  h_step = (image->height + MAX_HEIGHT - 1) / MAX_HEIGHT;
  step = (w_step > h_step)? w_step: h_step;
  if(step > 1) {
    sprintf(ss,"Warning : image has been reduced by a factor %d\n",step);
    PUT_WNDMSG(ss);
    image->width    = image->width/step;
    image->height   = image->height/step;
    image->dx = image->dx * step;
    image->dy = image->dy * step;
  }
	
  k=0;
  for(j=0;j<image->height;j++) {
    for(i=0;i<image->width;i++) {
      pic[i][j]=buff[k];
      k=k+step;
    }
    k=k+image->width*step*(step-1);	
  }


  return(k);
}

/* --load_mri() -------------------------------------------------------    
**
**    This function loads a mri.
**
**    Usage:    load_mri(file,graphic,numfirst)
**        file: char *: Designs file to load.
**        graphic: grphic:Where to stock file.
**              numfirst : numero de l'image dans le fichier
**
*/
grphic *load_mri(char *file, grphic *graphic, int numfirst)
{
  int format;
  FILE *fp;
  FILE *fp_bin;
  long *buff, max, i, min;
  long wdth, hght, size;
  int number_img;
  char *tmp;
  int n,idCMap;
  char file_bin[PATH_LEN];      /*For SIEMENS Vanderbilt format */



  format = is_mri(file,NULL);

  /* For interfile, IPB_2D, Paravision format */
  if(format == INTERFILE_FORMAT
     || format == IPB_2D_FORMAT
     || format == PARAVISION_FORMAT)
    {
      if ((number_img=getmri_numberimg(file))<numfirst)
        {
	  char s[256];
	  sprintf(s,"Sorry, cannot load image number %d only %d images \n"
		  ,numfirst,number_img);
	  PUT_MESG(s);
	  _imagix_error=IMX_PICTURE_NOT_FOUND;
	  _imagix_err_ext=errno;
	  return(NULL);
        }
    }

  /* for interfile, ipb, bruker, paravision, raw4, SMIS formats */
  if(format == INTERFILE_FORMAT
     || format == IPB_2D_FORMAT
     || format == IPB_3D_FORMAT
     || format == BRUKER_FORMAT
     || format == PARAVISION_FORMAT
     || format == RAW4_FORMAT
     || format == SMIS_FORMAT) {

    errno = 0;
    strncpy(graphic->filename, file, 256);
    graphic->file_number = numfirst;

    if((buff = getmri(file, numfirst)) == NULL) {
      char s[256];
      sprintf(s, "Sorry, cannot load %s...(*errno:%d) \n", file, errno);
      PUT_MESG( s);
      _imagix_error = IMX_PICTURE_NOT_FOUND;
      _imagix_err_ext = errno;
      return(NULL);
    }
		
    wdth = getmri_width(file, numfirst, MODE_2D);
    hght = getmri_height(file, numfirst, MODE_2D);
    graphic->width = wdth;
    graphic->height= hght;
    graphic->dx = getmri_dx(file, numfirst, MODE_2D);
    graphic->dy = getmri_dy(file, numfirst, MODE_2D);
    graphic->dz = getmri_dz(file, numfirst, MODE_2D);
      
    buff2pic(buff, graphic->mri, graphic);

    graphic->min_pixel = getmri_minpixel(file, numfirst, MODE_2D);
    graphic->max_pixel = getmri_maxpixel(file, numfirst, MODE_2D);
    graphic->cutoff_min = getmri_cutoff_min(file, numfirst, MODE_2D);
    graphic->cutoff_max = getmri_cutoff_max(file, numfirst, MODE_2D);
    graphic->icomp = getmri_icomp(file, numfirst, MODE_2D);
    graphic->rcoeff = (double)pow((double)2, (double)graphic->icomp);

    tmp = getmri_name(file);
    if (tmp != (char*)NULL) strcpy(graphic->patient.name, tmp);
    else memset(graphic->patient.name, 0, sizeof(graphic->patient.name));
    FREE(tmp);

    tmp = getmri_dbirth(file);
    if (tmp != (char*)NULL) strcpy(graphic->patient.d_birth, tmp);
    else memset(graphic->patient.d_birth, 0, sizeof(graphic->patient.d_birth));
    FREE(tmp);

    tmp = getmri_dexamen(file);
    if (tmp != (char*)NULL) strcpy(graphic->patient.d_examen, tmp);
    else memset(graphic->patient.d_examen, 0, sizeof(graphic->patient.d_examen));
    FREE(tmp);

    tmp = getmri_medecin(file);
    if (tmp != (char*)NULL) strcpy(graphic->patient.medecin, tmp);
    else memset(graphic->patient.medecin,0, sizeof(graphic->patient.medecin));
    FREE(tmp);

    graphic->pos =-1;        /* Not determined yet ! */
    graphic->img_type = IMAGE_NORMAL;
    graphic->type = 'b';
    graphic->bitppixel = 32;
    graphic->zoom = 1;
    graphic->zoom_x = 0;
    graphic->zoom_y = 0;

#ifdef __GTK_GUI
    graphic->cmapinfo.dspStyle = getmri_dspStyle(file, numfirst, MODE_2D);
    graphic->cmapinfo.dspOverBlancking = getmri_dspOverBlancking(file, numfirst, MODE_2D);

    idCMap = getmri_nCMap(file, numfirst, MODE_2D);
    n = (idCMap == 0 ? 1 : IMXColor_NewColormap(idCMap, 0, NULL, 0));

    IMXColor_BindColormap_ptr(graphic, n);
    IMXColor_UpdateColormap();
#endif /* __GTK_GUI */
        

    /* For Raw4 format */
    if(format == RAW4_FORMAT
       || format == SMIS_FORMAT)
      imx_iniparaimg_p(graphic);

    FREE(buff);     /* Free width*height*4 + 18432 */
    return((grphic *)&graphic);
  }

#ifdef __GTK_GUI
  /* For PCX format */
  if(format == PCX_FORMAT) {
    strncpy(graphic->filename,file,256);
    graphic->file_number=numfirst;
    wdth=getmri_width(file,0, MODE_2D);
    hght=getmri_height(file,0, MODE_2D);
    graphic->width =MINI(wdth,MAX_WIDTH);
    graphic->height=MINI(hght,MAX_HEIGHT);
    graphic->min_pixel=getmri_minpixel(file,0, MODE_2D);
    graphic->max_pixel=getmri_maxpixel(file,0, MODE_2D);
    graphic->cutoff_min=getmri_cutoff_min(file,0, MODE_2D);
    graphic->cutoff_max=getmri_cutoff_max(file,0, MODE_2D);
    graphic->icomp =getmri_icomp(file,0, MODE_2D);
    graphic->rcoeff=(double)pow((double)2,(double)graphic->icomp);
    /*protoize:???*/
    graphic->dx=getmri_dx(file, 0, MODE_2D);
    graphic->dy=getmri_dy(file, 0, MODE_2D);
    graphic->dz=getmri_dz(file, 0, MODE_2D);

    loadpcx_p(file,graphic,1);

    graphic->pos=-1;        /* Not determine yet ! */
    graphic->img_type =IMAGE_NORMAL;
    graphic->type='b';
    graphic->bitppixel=32;
    graphic->zoom = 1;
    graphic->zoom_x = 0;
    graphic->zoom_y = 0;
    return((grphic *)&graphic);
  }
#endif /* __GTK_GUI */

  /*---- For Siemens Magnetom Vision format, Feruglio Cyril 24-04-1998----*/

  if(format == SIEMENS_FORMAT) {
    strncpy(graphic->filename,file,256);
    printf("fichier: <%s>\n",file);
    graphic->file_number=numfirst;
    /*protoize:???*/
    if((buff=getmri(file,numfirst))==NULL) {
      char s[256];
      printf("loadmri error !");
      sprintf(s,"Sorry, cannot load %s...(errno:%d from imx_file) \n",file,errno);
      PUT_MESG(s);
      _imagix_error=IMX_PICTURE_NOT_FOUND;
      _imagix_err_ext=errno;
      return(NULL);
    }

    wdth=getmri_width(file,1, MODE_2D);
    hght=getmri_height(file,1, MODE_2D);
    graphic->width =wdth;
    graphic->height=hght;
    graphic->icomp =getmri_icomp(file,1, MODE_2D);
    graphic->rcoeff=(double)pow((double)2, (double)graphic->icomp);;
    graphic->dx=getmri_dx(file, 0, MODE_2D);
    graphic->dy=getmri_dy(file, 0, MODE_2D);
    graphic->dz=getmri_dz(file, 0, MODE_2D);
    graphic->pos=-1;        /* Not determine yet ! */
    graphic->img_type =IMAGE_NORMAL;
    graphic->type='b';
    graphic->bitppixel=16;
    buff2pic(buff,graphic->mri,graphic);
    max = 0;
    min = 65536;

    fp = fopen(file,"r");

    fseek(fp,0l,SEEK_END);
    size=(long)ftell(fp);
    fclose(fp);
    size=(size-6144)/2;
    for(i=0; i<size; i++){
      if(buff[i]> max)
	max =(long)buff[i];
      if(buff[i]< min)
	min = (long)buff[i];
    }
    graphic->min_pixel=min;
    graphic->max_pixel=max;
    graphic->cutoff_min=min;
    graphic->cutoff_max=max;
    tmp = getmri_name(file);
    if (tmp != (char*)NULL) strcpy(graphic->patient.name,tmp);
    else memset(graphic->patient.name, 0, sizeof(graphic->patient.name));
    FREE(tmp);
    tmp = getmri_dbirth(file);
    if (tmp != (char*)NULL) strcpy(graphic->patient.d_birth,tmp);
    else memset(graphic->patient.d_birth, 0, sizeof(graphic->patient.d_birth));
    FREE(tmp);
    tmp = getmri_dexamen(file);
    if (tmp != (char*)NULL) strcpy(graphic->patient.d_examen,tmp);
    else memset(graphic->patient.d_examen, 0, sizeof(graphic->patient.d_examen));
    FREE(tmp);
    tmp = getmri_medecin(file);
    if (tmp != (char*)NULL) strcpy(graphic->patient.medecin,tmp);
    else memset(graphic->patient.medecin, 0, sizeof(graphic->patient.medecin));
    FREE(tmp);
    FREE(buff);

    graphic->zoom = 1;
    graphic->zoom_x = 0;
    graphic->zoom_y = 0;

    return((grphic *)&graphic);

  }

  /*  For SIEMENS Vanderbilt FORMAT, par Feruglio Cyril 1998 */
  if(format == VANDERBILT_FORMAT) {
    strcpy(file_bin,file);
    sprintf(file_bin + strlen(file_bin)-12,"image.bin");
    if(MESG_DEBUG) printf("file: %s",file_bin);
    strncpy(graphic->filename,file_bin,256);
    graphic->file_number=numfirst;
    /*protoize:???*/ /*ok?*/
    if((buff=/*  (long) */getmri(file,numfirst))==NULL) {
      char s[256];
      printf("loadmri error !");
      sprintf(s,"Sorry, cannot load %s...(errno:%d from imx_file) \n",file,errno);
      PUT_MESG(s);
      _imagix_error=IMX_PICTURE_NOT_FOUND;
      _imagix_err_ext=errno;
      return(NULL);
    }
    wdth=getmri_width(file,1, MODE_2D);
    hght=getmri_height(file,1, MODE_2D);
    graphic->width =wdth;
    graphic->height=hght;
    if(MESG_DEBUG) printf("width='%ld', height='%ld'\n",wdth, hght);
    graphic->icomp =getmri_icomp(file,1, MODE_2D);
    graphic->rcoeff=1;
    graphic->dx=getmri_dx(file, 0, MODE_2D);
    graphic->dy=getmri_dy(file, 0, MODE_2D);
    graphic->dz=getmri_dz(file, 0, MODE_2D);
    graphic->pos=-1;        /* Not determine yet ! */
    graphic->img_type =IMAGE_NORMAL;
    graphic->type='b';
    graphic->bitppixel=getmri_dattype(file, 0, MODE_2D);
    if(MESG_DEBUG) printf("entree dans buff2pic");
    buff2pic(buff,graphic->mri,graphic);
    max = 0;
    min = 65536;
    fp_bin = fopen(file_bin,"r");
    fseek(fp_bin,0l,SEEK_END);
    size=(long)ftell(fp_bin);
    fclose(fp_bin);
    size-=6144;
    for(i=0; i<size; i++){
      if(buff[i]> max)
	max =(long)buff[i];
      if(buff[i]< min)
	min = (long)buff[i];
    }
    graphic->min_pixel=min;
    graphic->max_pixel=max;
    graphic->cutoff_min=min;
    graphic->cutoff_max=max;
    if(MESG_DEBUG) printf("min(%ld) & max(%ld) oki", min, max);
    /*  ???????? FREE(tmp); */
    FREE(buff);

    graphic->zoom = 1;
    graphic->zoom_x = 0;
    graphic->zoom_y = 0;

    return((grphic *)&graphic);
  }

  /* For ACR-NEMA 2.0 format Burger Laurent, novembre 1999 */
  if (format==ACRNEMA2L_FORMAT || format==ACRNEMA2B_FORMAT 
      || format==ACRNEMA1B_FORMAT || format==ACRNEMA1L_FORMAT 
      || format == DICMBIMP_FORMAT || format == DICMLIMP_FORMAT
      || format == DICMBEXP_FORMAT || format == DICMLEXP_FORMAT)
    {
      graphic->file_number = numfirst;

      if((buff = getmri(file, numfirst)) == NULL) {
	char s[256];
	sprintf(s, "Sorry, cannot load %s...(*errno:%d) \n", file, errno);
	PUT_MESG( s);
	_imagix_error = IMX_PICTURE_NOT_FOUND;
	_imagix_err_ext = errno;
	return(NULL);
      }
	
      wdth=getmri_width(file,1, MODE_2D);
      hght=getmri_height(file,1, MODE_2D);
      graphic->width =wdth;
      graphic->height=hght;
      graphic->icomp =getmri_icomp(file,1, MODE_2D);
      graphic->rcoeff=1;
      graphic->dx=getmri_dx(file, 0, MODE_2D);
      graphic->dy=getmri_dy(file, 0, MODE_2D);
      graphic->dz=getmri_dz(file, 0, MODE_2D);
      graphic->pos=-1;       
      graphic->img_type =IMAGE_NORMAL;
      graphic->type='b';
      graphic->bitppixel=getmri_dattype(file, 0, MODE_2D);
	
        
      buff2pic(buff,graphic->mri,graphic);
      /* recherche du maximum et du maximum dans la partie donn�s */
      max = 0;
      min = 65536;
      for(i=0; i<graphic->width*graphic->height; i++){
	if(buff[i]> max)
	  max =(long)buff[i];
	if(buff[i]< min)
	  min = (long)buff[i];
      }
		
      graphic->min_pixel=min;
      graphic->max_pixel=max;
      graphic->cutoff_min=min;
      graphic->cutoff_max=max;
        
      tmp = getmri_name(file);
      if (tmp != (char*)NULL) strcpy(graphic->patient.name,tmp);
      else memset(graphic->patient.name, 0, sizeof(graphic->patient.name));
      FREE(tmp);

      tmp = getmri_dbirth(file);
      if (tmp != (char*)NULL) strcpy(graphic->patient.d_birth,tmp);
      else memset(graphic->patient.d_birth, 0, sizeof(graphic->patient.d_birth));
      FREE(tmp);
      tmp = getmri_dexamen(file);
      if (tmp != (char*)NULL) strcpy(graphic->patient.d_examen,tmp);
      else memset(graphic->patient.d_examen, 0, sizeof(graphic->patient.d_examen));
      FREE(tmp);
      tmp = getmri_medecin(file);
      if (tmp != (char*)NULL) strcpy(graphic->patient.medecin,tmp);
      else memset(graphic->patient.medecin, 0, sizeof(graphic->patient.medecin));
      FREE(tmp);
      FREE(buff);

      graphic->zoom = 1;
      graphic->zoom_x = 0;
      graphic->zoom_y = 0;

      return((grphic *)&graphic);

    }

  /* For ECAT7 format */
  if (format==ECAT7_FORMAT)
    {
      graphic->file_number = 1;

      if((buff = getmri(file, numfirst)) == NULL) {
	char s[256];
	sprintf(s, "Sorry, cannot load %s...(*errno:%d) \n", file, errno);
	PUT_MESG( s);
	_imagix_error = IMX_PICTURE_NOT_FOUND;
	_imagix_err_ext = errno;
	return(NULL);
      }
	
      wdth=getmri_width(file,1, MODE_2D);
      hght=getmri_height(file,1, MODE_2D);
      graphic->width =wdth;
      graphic->height=hght;
      graphic->icomp =1;
      graphic->rcoeff=1;
      graphic->dx=getmri_dx(file, 0, MODE_2D);
      graphic->dy=getmri_dy(file, 0, MODE_2D);
      graphic->dz=getmri_dz(file, 0, MODE_2D);
      graphic->pos=-1;       
      graphic->img_type =IMAGE_NORMAL;
      graphic->type='b';
      graphic->bitppixel=getmri_dattype(file, 0, MODE_2D);
	
      for(i=0; i<graphic->width*graphic->height; i++)
	if((long)buff[i]>60000)buff[i]=0l;

      buff2pic(buff,graphic->mri,graphic);
      /* recherche du maximum et du maximum dans la partie donn�s */
      max = 0;
      min = 65536;
      for(i=0; i<graphic->width*graphic->height; i++){
	if(buff[i]> max)
	  max =(long)buff[i];
	if(buff[i]< min)
	  min = (long)buff[i];
      }
		
      graphic->min_pixel=min;
      graphic->max_pixel=max;
      graphic->cutoff_min=min;
      graphic->cutoff_max=max;

		

      tmp=(char*)malloc(50);
      tmp = getmri_name(file);
      if (tmp != (char*)NULL) strcpy(graphic->patient.name,tmp);
      else memset(graphic->patient.name, 0, sizeof(graphic->patient.name));
        

      FREE(tmp);
      FREE(buff);

      graphic->zoom = 1;
      graphic->zoom_x = 0;
      graphic->zoom_y = 0;

      return((grphic *)&graphic);
    }

  /* 96/08/29  HMO  
     " La fonction doit retourner NULL si la lecture n'est pas faite " */

  return (NULL);
}


/****  check_file_2d()  ********************************************
 **
 **  >> Last user action by: F. BOULEAU on Aug. 10th, 1999
 **
 **  Controle l'autorisation de creer ou modifier le fichier file
 **
 **  Si file existe alors on demande 
 **    si on ajoute au fichier
 **  Si la reponse est non on demande si on supprime le fichier
 **  Si reponse oui le fichier est vide mais existe.
 **
 **  retourne 0 en cas de succes, c.a.d on peu creer ou ajouter au fichier
 **  Retourne 1 si impossibilite de sauvegarder dans  le fichier
 **             ou si on ne le souhaite pas ecraser le fichier.
 **
 ************/
int check_file_2d(char *file)
{
  char headerfile[PATH_LEN];
  char datafile[PATH_LEN];
  char crvfile[PATH_LEN];
  char s[64];
  struct stat statbuff;
  int erreur;


  strcpy(headerfile, file);
  strcpy(datafile, file);
  strcpy(crvfile, file);
  strcat(headerfile, ".ipb");
  strcat(datafile, ".img");
  strcat(crvfile, ".crv");

  if(stat(headerfile, &statbuff) == 0) 
    {
      if(GET_EXIST(TEXT0166, &erreur) != 1) 
	{
	  if(GET_EXIST(TEXT0169, &erreur) != 1) 
	    {
	      return(1);  /* On ne touche pas au fichier */
            }
			
	  /* On supprime les fichiers .ipb, .crv et .img */

	  if(remove(headerfile) == -1) 
	    {
	      sprintf(s, "[save_header] unable to remove file \"%s\"\n",
		      headerfile);    
	      PUT_ERROR( s);
	      return(1); /* Erreur */
	    }
          	
	  if(remove(datafile) == -1) 
	    {
	      sprintf(s, "[save_header] unable to remove file \"%s\"\n",
		      datafile);
	      PUT_ERROR( s);
	      return(1); /* Erreur */
	    }

	  if(stat(crvfile, &statbuff) == 0) 
	    {
	      /* le fichier existe ? */
	      if(remove(crvfile) == -1) 	/* on le supprime */
		{ 
		  sprintf(s, "[save_header] unable to remove file \"%s\"\n",
			  crvfile);
		  PUT_ERROR( s);
		  return(1); 				/* Erreur */
		}
	      return(0);
	    }
	}
    }

  return(0); /* Le fichier n'existe pas, on peut le creer */
}
/****  check_for_mask()  ********************************************
 **
 **  >> Last user action by: T. BERST on Aug. 10th, 1999
 **
 **  Controle l'existence d'un mask associe a une image
 *	return 0 si pas de mask, !=0 sinon
 *
 *********************************************************************/
int check_for_mask(char *file)
{
  char answer[80];
  strcpy(answer,getheader_interfile(file,"MASK_FILE=",STRING,(int)NULL));
  return strlen(answer);
}

/****  check_file_3d()  ********************************************
 **
 **  >> Last user action by: F. BOULEAU on Aug. 10th, 1999
 **
 **  Controle l'autorisation de creer ou modifier le fichier file
 **
 **  Si file existe alors on demande 
 **    si on ajoute au fichier
 **  Si la reponse est non on demande si on supprime le fichier
 **  Si reponse oui le fichier est vide mais existe.
 **
 **  retourne 0 en cas de succes, c.a.d on peu creer ou ajouter au fichier
 **  Retourne 1 si impossibilite de sauvegarder dans  le fichier
 **             ou si on ne le souhaite pas ecraser le fichier.
 **
 ************/
int check_file_3d(char *file)
{
  char file_without_ext[PATH_LEN];
  char headerfile[PATH_LEN];
  char datafile[PATH_LEN];
  char crvfile[PATH_LEN];
  char s[64];
  struct stat statbuff;
  int erreur;

  remove_file_extension(file,file_without_ext);
  strcpy(headerfile, file_without_ext);
  strcpy(datafile, file_without_ext);
  strcpy(crvfile, file_without_ext);
  strcat(headerfile, ".ipb");
  strcat(datafile, ".img");
  strcat(crvfile, ".crv");

  if(stat(headerfile, &statbuff) == 0) 
    {
      if(GET_EXIST(TEXT0166, &erreur) != 1) 
	{
	  if(GET_EXIST(TEXT0169, &erreur) != 1) 
	    return(1);  /* On ne touche pas au fichier */
	 	 
	  /* On supprime les fichiers .ipb, .crv et .img */

	  if(remove(headerfile) == -1) 
	    {
	      sprintf(s, "[save_header] unable to remove file \"%s\"\n",
		      headerfile);	 
	      PUT_ERROR( s);
	      return(1); /* Erreur */
	    }  
         
	  if(remove(datafile) == -1) 
	    {
	      sprintf(s,"[save_header] unable to remove file \"%s\"\n",
		      datafile);
	      PUT_ERROR( s);
	      return(1); /* Erreur */
	    }
	 	 
	  if(stat(crvfile, &statbuff) == 0) 
	    { 
	      /* le fichier existe ? */
              if(remove(crvfile) == -1)  /* on le supprime */
	 	{ 
		  sprintf(s, "[save_header] unable to remove file \"%s\"\n",
			  crvfile);
		  PUT_ERROR( s);
		  return(1);				 /* Erreur */
		}
              return(0);
	    }
	}
    }
	
  return(0); /* Le fichier n'existe pas, on peut le creer */
}

#ifdef __GTK_GUI
/* --save_mri_bruker() --------------------------------------------------    
**
**    save_mri_bruker() save a picture.
**
**    Usage ...
*/
grphic *save_mri_bruker(char *file, grphic *graphic)
{ 
  long *buff,header[HEADER_SIZE],i,j,wdth,hght;
  long k=0l;
  FILE *fp,*fp_src;
  char str[80];
  int err;


  /* Save mri */
  /* 1st -Check if a mri filename already exists !
  ** 2nd -Save file if possible
  */

  /* Ist  PART: Check if mri filename already exists. */
  errno=0;
  if((fp=fopen(file,"rb"))==NULL) {
    if(errno==2) {  /* Not such file or directory */

      /* -> Save ! */         }
    else { /* Error happened ! */
      char s[256];
      sprintf(s,"Error occured:\n  Cannot save mri in file >%s<\n",file);
      PUT_ERROR(s);
      fclose(fp);
      return(NULL);
    } 
  }
  else { /* File exists */
    char s[256];

    sprintf(s,"Delete old file ?\n  (%s)\n",file);
    if(XGET_EXIST(s, &err)==0) {

      sprintf(s,"Cannot save mri in file >%s<\n",file);
      PUT_ERROR(s);
      fclose(fp);
      return(NULL);
      fclose(fp);
    } 
  }

  /* IInd PART: Save ROI */
  /* 1st: Save an header... */
  if((fp=fopen(file,"wb"))==NULL) { 
    char s[256];
    sprintf(s,"Cannot write mri in file>%s<.\n",file);
    PUT_ERROR(s);
    fclose(fp);
    return((grphic *)NULL);
  }

  sprintf(str,graphic->filename);
  if((strlen(graphic->filename)<=0)
     || (is_mri("",graphic->filename)!=BRUKER_FORMAT) )
    {
      sprintf(str,getenv("NAME_HEADER_FILE"));
      if((fp_src=fopen(str,"rb"))==NULL) {
	PUT_ERROR("Cannot save\n  Cannot open an header !");
	fclose(fp);
	return((grphic *)NULL);
      }
    }
  else
    if((fp_src=fopen(str,"rb"))==NULL) {
      PUT_ERROR("Cannot save\n  Cannot open for an header !");
      fclose(fp);
      return((grphic *)NULL);
    }
  fseek(fp_src,0l,SEEK_SET);
  /*fread(header,sizeof(long),sizeof(header)/4,fp_src);*/
  FREAD_LONG(header,sizeof(header)/4,fp_src);

  wdth=graphic->width;
  hght=graphic->height;

  buff=CALLOC(wdth*hght,long);

  for(j=0;j<hght;j++)
    for(i=0;i<wdth;i++) {
      buff[k++]=graphic->mri[i][j];
    }

  sprintf(str,"%d",graphic->height);
  /*protoize:???*/
  codify_pos((unsigned char*)header/*???*/,MVAR_OFF+3*(2+3),ENTIER,str);    /* Heigth */

  sprintf(str,"%d",graphic->width);
  codify_pos((unsigned char*)header/*???*/,MVAR_OFF+3*(2+4),ENTIER,str);    /* Width  */
  sprintf(str,"%ld",graphic->max_pixel);
  codify_pos((unsigned char*)header/*???*/,MVAR_OFF+3*(2+5),ENTIER,str);      /* Max PIX */
  sprintf(str,"%ld",graphic->min_pixel);
  codify_pos((unsigned char*)header/*???*/,MVAR_OFF+3*(2+6),ENTIER,str);      /* min PIX */
  sprintf(str,"%f",graphic->icomp);
  codify_pos((unsigned char*)header/*???*/,MVAR_OFF+3*(2+7),ENTIER,str);      /* Global exp */

  fseek(fp,0l,SEEK_SET);
  /*fwrite(header,sizeof(long),sizeof(header)/4,fp);*/
  FWRITE_LONG(header,sizeof(header)/4,fp);

  fseek(fp,HEADER_SIZE,SEEK_SET);
  /*fwrite(buff,sizeof(long),k,fp);*/
  FWRITE_LONG(buff,k,fp);

  fclose(fp);
  fclose(fp_src);
  FREE(buff);

  /* Unavailable ! 
     putmri_icomp (file,file,graphic->icomp);
     putmri_maxpix(file,file,graphic->max_pixel);
     putmri_minpix(file,file,graphic->min_pixel);
     putmri_width (file,file,graphic->width);
     putmri_height(file,file,graphic->height);
  */

  PUT_INFO("Mri saved.\n");
  sprintf(graphic->filename,file);
  return((grphic *)NULL);
}

/* --save_pictures_2d() ---------------------------------------------------------
**
**    save_pictures_2d() saves more pictures in IPB 2D file format.
**
**    Usage ...
**   file :       nom du fichier
**   start :      numero de la premier image a sauver 
**   num_images : nombre d'images a sauver 
**   step :       pas si plusieurs images
**   format : format image a sauver (defut IPB_FORMAT)
**
*/
int save_pictures_2d(char *file, int start, int num_images, int step, int format)
{
  int err;
  
  switch(format) {
  
  case IPB_2D_FORMAT:
    if (check_file_2d(file)!=0 ) return (0);
    err=save_mri_ipb_2d(file, start, num_images, step);
    break;


  default:
    if (check_file_2d(file)!=0 ) return (0);
    err=save_mri_ipb_2d(file, start, num_images, step);
    
    break;
  }

  return(err);
}


/* --save_mri_ipb_2d() ---------------------------------------------------------
**
**    save_mri_ipb_2d() saves a picture in IPB 2D file format.
**
**    Usage ...
**   file :       nom du fichier
**   start :      numero de la premier image a sauver 
**   num_images : nombre d'images a sauver 
**   step :       pas si plusieurs images
**
**   return : 1 if Ok  !=1 if error
*/
int     save_mri_ipb_2d(char *file, int start, int num_images, int step)
{
  grphic *image;
  int end;
  char s[256];
    
  end = (num_images - 1) * step + start;
  if(end > MAX_PICTURES) {
    end = MAX_PICTURES;
    PUT_WARN("[save_mri_ipb_d] warning: cannot save all images !\n");
  }
  
  while(start <= end) {
    if((image = ptr_img(start)) == NULL) {
      sprintf(s,"[save_mri_ipb_2d] unable to get image %d info!\n", start);
      PUT_ERROR(s);
      return(0);
    }
  
    if ( save_mri_ipb_2d_p(file, image) == 0) 
      return(0) ;
    start += step;
  }

  return(1);
}

/* --save_mri_ipb_2d_p() ---------------------------------------------------------
**
**    save_mri_ipb_2d_p() saves a picture in IPB 2D file format.
**
**    Usage ...
**   file :       nom du fichier
**   image : pointeur sur structure grphic 
**
**   return : 1 if Ok  !=1 if error
*/
int     save_mri_ipb_2d_p(char *file, grphic *image)
{
  struct stat statbuff;
  char s[256];
  int  index;
  long entier;
  double reel;
  long offset;
  char *buffchar = NULL;
  short *buffshort = NULL;
  long *bufflong = NULL;
  void *buffer;
  int i, j, num_items;
  FILE *fp;
  char headerfile[256], datafile[256];
  HEADER header;

  strcpy(headerfile, file);
  strcat(headerfile, ".ipb");
  strcpy(datafile, file);
  strcat(datafile, ".img");

  index = 0;
  offset = 0;

  if(stat(headerfile, &statbuff) == 0) { /* File exists */
    if((index = getmri_numberimg(headerfile)) == 0) {
      PUT_ERROR("[save_mri_ipb_2d] Error reading header (1)!\n");
      return(0);
    }

    for(i=0; i<index; i++) {
      int width, height, data_size;

      width = getmri_width(headerfile, i+1, MODE_2D);
      height = getmri_height(headerfile, i+1, MODE_2D);
      data_size = getmri_dattype(headerfile, i+1, MODE_2D);
      offset += width * height * data_size;
    }    		
  }

  index ++;

  if((fp = fopen(datafile, "ab")) == NULL) {
    sprintf(s,"Unable to save file %s !\n", file);
    PUT_WNDMSG(s);
    return(0);
  }

  if(load_header(headerfile, &header) == 0) {
    PUT_ERROR("[save_mri_ipb_2d] Error reading header (3)!\n");
    fclose(fp);
    return(0);
  }

  put_header(&header, "IPB2", STRING, (void *)"", 0);
  put_header(&header, "patient name=",  STRING, image->patient.name,    0);
  put_header(&header, "date of birth=", STRING, image->patient.d_birth, 0);
  put_header(&header, "doctor name=",   STRING, image->patient.medecin, 0);
  put_header(&header, "date of exam=",  STRING, image->patient.d_examen,0);
  put_header(&header, "number of images=", ENTIER, &index, 0);

  entier = (image->cmapinfo.cmap ? image->cmapinfo.cmap->noMode : 0);
  put_header(&header, "colormap=", ENTIER, &entier, index);
  put_header(&header, "display style=", ENTIER, &image->cmapinfo.dspStyle, 
	     index);
  put_header(&header, "overblancking=", ENTIER, &image->cmapinfo.dspOverBlancking,
	     index);

  entier = image->width;	put_header(&header, "width=", ENTIER, &entier, 
					   index);
  entier = image->height;	put_header(&header, "height=", ENTIER, &entier, 
					   index);
  reel   = image->rcoeff;	put_header(&header, "rcoeff=", REEL, &reel, index);
  reel   = image->icomp;	put_header(&header, "icomp=", REEL, &reel, index);
  entier = image->max_pixel;  put_header(&header, "max pixel=", ENTIER, 
					 &entier, index);
  entier = image->min_pixel;  put_header(&header, "min pixel=", ENTIER, 
					 &entier, index);
  entier = image->cutoff_max; put_header(&header, "cutoff max=", ENTIER, 
					 &entier, index);
  entier = image->cutoff_min; put_header(&header, "cutoff min=", ENTIER, 
					 &entier, index);
  entier = image->bitppixel;  put_header(&header, "bits per pixel=", ENTIER, 
					 &entier, index);
  reel   = image->dx; 	put_header(&header, "dx=", REEL, &reel, index);
  reel   = image->dy; 	put_header(&header, "dy=", REEL, &reel, index);
  reel   = image->dz; 	put_header(&header, "dz=", REEL, &reel, index);
  entier = offset;		put_header(&header, "offset in raw=", ENTIER, &entier, index);

  offset += image->width * image->height * (image->bitppixel/8);

  num_items = image->width * image->height;

  switch(image->bitppixel/8)
    {
    case sizeof(char):
      buffchar = (char *)malloc(num_items);
      if(buffchar == NULL) {
	PUT_ERROR("[save_mri_ipb_2d] memory allocation error (1)\n");
	return(0);
      }
      buffer = buffchar;
      for(j=0; j<image->height; j++)
	for(i=0; i<image->width; i++)
	  *buffchar ++ = image->mri[i][j];
      break;

    case sizeof(short):
      buffshort = (short *)malloc(num_items * sizeof(short));
      if(buffshort == NULL) {
	PUT_ERROR("[save_mri_ipb_2d] memory allocation error (1)\n");
	return(0);
      }
      buffer = buffshort;
      for(j=0; j<image->height; j++)
	for(i=0; i<image->width; i++)
	  if (_is_bigEndian) *buffshort ++ = image->mri[i][j];
	  else	*buffshort ++ = shEndianConversion(image->mri[i][j]);
      break;

    case sizeof(long):
      bufflong = (long *)malloc(num_items * sizeof(long));
      if(bufflong == NULL) {
	PUT_ERROR("[save_mri_ipb_2d] memory allocation error (1)\n");
	return(0);
      }
      buffer = (void*) bufflong;
      for(j=0; j<image->height; j++)
	for(i=0; i<image->width; i++)
	  if (_is_bigEndian) *bufflong ++ = image->mri[i][j];
	  else	*bufflong ++ = longEndianConversion(image->mri[i][j]);
      /* *bufflong ++ = image->mri[i][j]; */
      break;

    default:
      PUT_ERROR("[save_mri_ipb_2d] invalid data size\n");
      return(0);
    }

  if(fwrite(buffer, image->bitppixel/8, num_items, fp) != num_items) {
    PUT_ERROR("[save_mri_ipb_2d] error writing data !\n");
    fclose(fp);
    return(0);
  }

  FREE(buffer);

  index ++;

  entier = index - 1;
  put_header(&header, "number of images=", ENTIER, &entier, 0);

  if(save_header(&header, headerfile) == 0) {
    sprintf(s,"Unable to save file header %s !\n", file);
    PUT_ERROR(s);
    fclose(fp);
    return(0);
  }

  fclose(fp);
  return(1);
}

/* --save_mri_raw32() -------------------------------------------------------    
**
**   Permet de sauver une image dans un fichier
**   en raw data 32 bits
**
******************************************************/
void save_mri_raw32 (int pos, char *file)
     /* 32 bits */
        
           
{
  int i,j,k,err,*buff;
  grphic *img_prt = NULL;
  FILE *hand;

  if((hand=fopen(file,"rb"))==NULL)
    {
      if(!errno==2)              /* else file does not exist ! */
        { /* Error appended ! */

	  char s[256];
	  _imagix_err_ext=errno;
	  sprintf(s,"Error occured:\n  Cannot save mri in file >%s<\n",file);
	  if(_dialog_messages_status&0x00000002l)
	    { PUT_ERROR(s);
	    }
	  fclose(hand);
	  return;
        }
    }
  else { /* File exists */

    char s[256];
    sprintf(s,"Delete old file ?\n  (%s)\n",file);
    if(XGET_EXIST(s, &err)==0)
      {
	sprintf(s,"Cannot save mri in file >%s<\n",file);
	PUT_ERROR(s);
	fclose(hand);
	return;
      }
  }


  if(pos>0 && pos<=MAX_PICTURES)
    {
      if((hand=fopen(file,"w"))==NULL) { 
	char s[256];
	sprintf(s,"Cannot write mri in file>%s<.\n",file);
	PUT_ERROR(s);
	fclose(hand);
	return;
      }
      img_prt=ptr_img(pos);
      /* Preparation+Sauvegarde du tableau mri en int ('256'*'256' * 4)  */
      k=0;
      buff=CALLOC(MAX_WIDTH*MAX_HEIGHT,int);
      for (j=0;j<MAX_HEIGHT;j++)
	for (i=0;i<MAX_WIDTH;i++)
	  buff[k++]=(int)img_prt->mri[i][j];
      /*fwrite( buff ,sizeof(int), k, hand);*/
      FWRITE_LONG((void *) buff,k,hand);
      fclose(hand);
      /*    Enregistrement fin             */
    }
  sprintf(img_prt->filename,file);
}

/**************************************************
 **   Permet de sauver une image dans un fichier
 **   en raw data 8 bits
 **
 ******************************************************/

void save_mri_raw8 (int pos, char *file)
     /* 8 bits */
        
           
{
  int i,j,k,err;
  int *buff;
  unsigned char *bufftmp;
  grphic *img_prt = NULL;
  FILE *hand;

  if((hand=fopen(file,"rb"))==NULL)
    {
      if(!errno==2)              /* else file does not exist ! */
        { /* Error appended ! */

	  char s[256];
	  _imagix_err_ext=errno;
	  sprintf(s,"Error occured:\n  Cannot save mri in file >%s<\n",file);
	  if(_dialog_messages_status&0x00000002l)
            { PUT_ERROR(s);
            }
	  fclose(hand);
	  return;
        }
    }
  else { /* File exists */

    char s[256];
    sprintf(s,"Delete old file ?\n  (%s)\n",file);
    if(XGET_EXIST(s, &err)==0)
      {
	sprintf(s,"Cannot save mri in file >%s<\n",file);
	PUT_ERROR(s);
	fclose(hand);
	return;
      }
  }



  if(pos>0 && pos<=MAX_PICTURES)
    {
      if((hand=fopen(file,"w"))==NULL) { 
	char s[256];
	sprintf(s,"Cannot write mri in file>%s<.\n",file);
	PUT_ERROR(s);
	fclose(hand);
	return;
      }
      img_prt=ptr_img(pos);
      /* Preparation+Sauvegarde du tableau mri en unsigned char ('256'*'256' * 1)  */
      k=0;
      buff=CALLOC(MAX_WIDTH*MAX_HEIGHT,int);
      bufftmp=CALLOC(MAX_WIDTH*MAX_HEIGHT,unsigned char);
      for (j=0;j<MAX_HEIGHT;j++)
        {
	  for (i=0;i<MAX_WIDTH;i++)
            {
	      buff[k++]=(int) ((float)img_prt->mri[i][j]*255./(float)img_prt->max_pixel);
	      bufftmp[k-1]=buff[k-1];
	      /*printf ("bufftmp =%c = rapport(%f)  *  buff =%d\n",bufftmp[k-1],255./(float)img_prt->max_pixel,img_prt->mri[i][j]);
	       */
            }
        }
      fwrite( (unsigned char*) bufftmp, sizeof(unsigned char), k, hand);
      /* printf ( "result= %d\n",fwrite( (int) buff,  sizeof(int),  k, hand)); */
      fclose(hand);
      /*    Enregistrement fin             */
    }
  sprintf(img_prt->filename,file);
}
#endif /* __GTK_GUI */


/* ---- put_line() ---------------------------------------
**
**
**-----------------------------------------*/
void put_line(byte *line, FILE *hand)
{
  /*protoize:???*/ /* conversions int <-> byte * !?*/
  byte *ptr,*end;
  int current,previous,count;

  ptr=line;
  end=ptr+255;
  previous=*ptr++;
  count=1;

  while (ptr < end)
    {
      current=*ptr++;
      if(current==previous && count < 63)
        {
	  count++;
	  continue;
        }
      if (count>1||(previous & 0xc0)==0xc0)
        {
	  count |= 0xc0;
	  put_byte(count,hand);
        }
      put_byte(previous,hand);
      previous=current;
      count=1;
    }
  if (count>1||(previous & 0xc0)==0xc0)
    {
      count |= 0xc0;
      put_byte(count,hand);
    }
  put_byte(previous,hand);
}

/* ---- put_word() ---------------------------------------
**
**  Ecrit un mot dans le fichier decrit par hand
**
**-----------------------------------------*/
void put_word(int i, FILE *hand)
{	    
  if (_is_bigEndian)
    {
      put_byte(i & 0xff,hand);
      put_byte((i>>8) & 0xff,hand);
    }
  else
    {
      put_byte((i>>8)& 0xff,hand);
      put_byte(i & 0xff,hand);
    }

}


/* ---- put_byte() ---------------------------------------
**
**  Ecrit un byte dans un fichier decrit par hanh (pointeur)
**
**-----------------------------------------*/
void put_byte(byte b, FILE *hand)
{
  fputc(((byte)(b)),hand);
}


#ifdef __GTK_GUI
/**************************************************
 **   Permet de sauver plusieurs images dans des fichiers
 **   en raw data 8 bits
 **
 ******************************************************/

int     save_mri_nraw8(int first_pos, int nb_img, char **pictures)
{ 
  int i,nb_loaded=0;
  if(first_pos+nb_img>(MAX_PICTURES+1))
    return(-1);
  if(first_pos<1 || nb_img<1) return(-1);
  for(i=0;i<nb_img;i++) {
    if(pictures[i]!=NULL && pictures!=NULL) {
      save_mri_raw8(first_pos+i,pictures[i]);
      nb_loaded++;
    }
  }
  return(nb_loaded);
}

/**************************************************
 **   Permet de sauver plusieurs images dans des fichiers
 **   en raw data 32 bits
 **
 ******************************************************/

int     save_mri_nraw32(int first_pos, int nb_img, char **pictures)
{ 
  int i,nb_loaded=0;
  if(first_pos+nb_img>(MAX_PICTURES+1))
    return(-1);
  if(first_pos<1 || nb_img<1) return(-1);
  for(i=0;i<nb_img;i++) {
    if(pictures[i]!=NULL && pictures!=NULL) {
      save_mri_raw32(first_pos+i,pictures[i]);
      nb_loaded++;
    }
  }
  return(nb_loaded);
}

/* --- save_mri_npcx() -------------------------------------
**
**   Permet de sauver plusieurs images dans des fichiers
**   en format PCX
**
******************************************************/
int     save_mri_npcx(int first_pos, int nb_img, char **pictures)
{ 
  int i,nb_loaded=0;
  if(first_pos+nb_img>(MAX_PICTURES+1))
    return(-1);
  if(first_pos<1 || nb_img<1) return(-1);
  for(i=0;i<nb_img;i++) {
    if(pictures[i]!=NULL && pictures!=NULL) {
      save_mri_pcx(first_pos+i,pictures[i]);
      nb_loaded++;
    }
  }
  return(nb_loaded);
}




/* Public functions:
** ----------------- */

int    Affiche(void)
{ 
  int i;

  for(i=0;i<MAX_PICTURES;i++)
    load_mri("picture",(grphic *)&_imx_img[i],0);

  if(MESG_DEBUG) printf("Loaded\n");

  for(i=0;i<(MAX_PICTURES);i++)
    IMXWidget_ShowMRI(i+1,(grphic *)&_imx_img[i]);

  return(MAX_PICTURES);
}


int    load_picture(int pos, char *file, int numfirst)
{

  char *buff;
  int n;

  if(pos>MAX_PICTURES || pos<-MAX_PICTURES ) {
    PUT_WARN("Error Loading Image |position|>MAX_PICTURES_3D\n"); 
    return(-1);
  }

  if(_bIsParalActivated)
    {
      // TODO:PARAL
      n    = strlen(file);
      buff = CALLOC(n + 3, char);
      strcpy(buff, file);
        
      if(buff[n - 4] == '.')
	buff[n - 4] = '\0';
       
      strcat(buff, ".*");

      file = IMXAtm_RCP(buff);

      FREE(buff);
    }

  if(load_mri(file,(grphic *)&_imx_img[pos-1],numfirst)==NULL) {
    _imagix_error=IMX_BAD_POSITION;
    _imagix_err_ext=errno;
    return(-1);
  }
  if (is_mri(file,NULL) == IPB_2D_FORMAT && check_for_mask(file) )
    {
      char answer[80];
      strcpy(answer,getheader_interfile(file,"MASK_FILE=",STRING,(int)NULL));
      answer[strlen(answer)-1]='\0';   // il y a toujours un espace en trop a la fin
      if(load_mri(answer,ptr_img(-pos),numfirst)==NULL) {
	_imagix_error=IMX_BAD_POSITION;
	_imagix_err_ext=errno;
	return(-1);
      }
      else
	{
	  ptr_mask_activate(pos,TRUE);
	}
      //	((grphic *)&_imx_img[pos-1])->mask_active=TRUE;
    }

  return(IMXWidget_ShowMRI(pos,(grphic *)&_imx_img[pos-1]));
}

int    load_pictures(int first_pos, int nb_img, char **pictures)
{ 
  int i,nb_loaded=0;

  if(MESG_DEBUG) printf("--->>II %d %d \n",first_pos,nb_img);
  if(first_pos+nb_img>(MAX_PICTURES+1))
    return(-1);
  if(first_pos<1 || nb_img<1) return(-1);

  if(MESG_DEBUG) printf("Aventi. II frst %d  nb %d\n",first_pos,nb_img);
  for(i=0;i<nb_img;i++) {
    if(MESG_DEBUG) printf("Aventi i%d... >%s<\n",i,pictures[i]);
    if(pictures[i]!=NULL && pictures!=NULL) {
      load_picture(first_pos+i,pictures[i],1);
      nb_loaded++;
    }
  }

  return(nb_loaded);
}

int    save_picture(int pos, char *file)
{
  if(pos<=0 || pos>MAX_PICTURES) {
    _imagix_error=IMX_BAD_POSITION;
    _imagix_err_ext=pos;
    PUT_ERROR("Cannot save picture (save_pictiure)");
    return(-1);
  }

  if(MESG_DEBUG) printf("Save> %d\n",pos);
  if(pos>0 && pos<=MAX_PICTURES) {
    if(MESG_DEBUG) printf("Save picture from %d\n",pos);
    save_mri_bruker(file,(grphic *)&_imx_img[pos-1]);
  }

  return(pos);
}


int    save_pictures(int first_pos, int nb_img, char **pictures)
{
  int i, nb_loaded = 0;

  if(first_pos+nb_img>(MAX_PICTURES+1))
    return(-1);

  if(first_pos<1 || nb_img<1) return(-1);

  for(i=0;i<nb_img;i++) {
    if(pictures[i]!=NULL && pictures!=NULL) {
      save_picture(first_pos+i, pictures[i]);
      nb_loaded++;
    }
  }

  return(nb_loaded);
}

/******************************************************************************
 ** -- Open_Cmd() --------------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July 22nd, 1999
 **
 **  Opens a 2d picture.
 **
 ******************************************************************************/

void Open_Cmd(void)
{ /* Call file selection box and:            */
  /* Open Bruker picture -  Lire une image bruker.    */
  char *ptr;
  int pos,numfirst,e;

  ptr=GET_FILE_IMG(NULL, &numfirst, &e);
  if (e) return;
  /*
    if(MESG_DEBUG) printf("File is >%s<\n",ptr);
  */
  if(ptr==NULL || strlen(ptr)<=0)
    {
      PUT_WARN("Cannot read picture\n");
    }
  else
    {
      switch(is_mri(ptr,NULL))
	{
	case BRUKER_FORMAT :
	  pos=GET_WND(NULL);
	  load_picture(pos,ptr,1);
	  break;

	case PARAVISION_FORMAT :
	  pos=GET_WND(NULL);
	  load_picture(pos,ptr,numfirst);
	  break;

	case INTERFILE_FORMAT :
	  pos=GET_WND(NULL);
	  load_picture(pos,ptr,numfirst);
	  break;

	case RAW4_FORMAT :
	  pos=GET_WND(NULL);
	  load_picture(pos,ptr,1);
	  break;

	case SIEMENS_FORMAT :
	  pos=GET_WND(NULL);
	  load_picture(pos,ptr,1);
	  break;

	case VANDERBILT_FORMAT :
	  pos=GET_WND(NULL);
	  load_picture(pos,ptr,1);
	  break;

	case PCX_FORMAT :
	  pos=GET_WND(NULL);
	  load_picture(pos,ptr,1);
	  break;
            
	case IPB_3D_FORMAT :
	case IPB_2D_FORMAT :
	  pos=GET_WND(NULL);
	  load_picture(pos,ptr,numfirst);
	  break;
            
	case SMIS_FORMAT:
	  pos = GET_WND( NULL);
	  load_picture(pos, ptr, 1);
	  break;
			
	case ACRNEMA2L_FORMAT:
	  pos = GET_WND( NULL);
	  load_picture(pos, ptr, 1);
	  break;

	case ACRNEMA2B_FORMAT:
	  pos = GET_WND( NULL);
	  load_picture(pos, ptr, 1);
	  break;            

	case ACRNEMA1B_FORMAT :
	  pos = GET_WND( NULL);
	  load_picture(pos, ptr, 1);
	  break;            

	case ACRNEMA1L_FORMAT :
	  pos = GET_WND( NULL);
	  load_picture(pos, ptr, 1);
	  break;            

	case DICMBIMP_FORMAT :
	case DICMLIMP_FORMAT :
	case DICMBEXP_FORMAT :
	case DICMLEXP_FORMAT :
	  pos = GET_WND( NULL);
	  load_picture(pos, ptr, 1);
	  break;

		
	case ECAT7_FORMAT :
	  numfirst= GET_INT(TEXT0183, 0, &e);
	  if (e) break;
	  pos = GET_WND( NULL);
	  load_picture(pos, ptr, numfirst);
	  break;
	default: 
	  break;
	}
    }
}

#endif /* __GTK_GUI */
/****************************************************************/
/*   Fonctions de lecture dans un fichier binaire:        */
/*    lis du char de l'u.int du float ou du double selon la   */
/*   fonction appellee, les fonction prennent un offset en      */
/*  argument qui permet de dire ou l'on commence la lecture,    */
/*  dans le cas des chars on doit aussi preciser (2ieme arg.)   */
/* la taille de la chaine de caracteres que l'on veut lire      */
/*                                */
/*              Par Feruglio Cyril le 11 Mai 1998        */
/****************************************************************/



/*------------------ getchar_header -----------------------*/

char* getchar_header(char *fichier, int offset, int taille)
{
  FILE *fp;
  char *lu=(char *)calloc(taille+1, sizeof(char));
  fp = fopen(fichier, "r");
  if(fp == NULL){
    char e[256];
    sprintf(e,"getchar_header erreur: Impossible d'ouvrir le fichier <%s>.\n", fichier);
    PUT_ERROR(e);
    fclose(fp);
    return(NULL);
  }
  /*fflush(stdin);
    fflush(stdout);*/
  rewind(fp);
  fseek(fp, offset, SEEK_SET);
  fread(lu, sizeof(char), taille, fp);
  lu[taille]='\0';
  fclose(fp);
  return(lu);
  FREE(lu);   /* on nettoie... */
}


/*------------------ getuint_header ---------------------*/

unsigned int getuint_header(char *fichier, int offset)
{

  FILE *fp;
  unsigned int lu = 0;

  fp = fopen(fichier, "r");
  if(fp == NULL){
    char e[256];
    sprintf(e, "getuint_header erreur: Impossible d'ouvrir le fichier <%s>.\n", fichier);
    PUT_ERROR(e);
    fclose(fp);
    return(0);
  }

  /*fflush(stdin);
    fflush(stdout);*/
  rewind(fp);
  fseek(fp, offset, SEEK_SET);
  /*protoize:???*/ /* !!!!! */
  fread((void *)&lu/*  lu */, sizeof(unsigned int), 1, fp);
  fclose(fp);
  printf( "getuint_header: %d",lu);
  return(lu);
}


/*------------------- getdouble_header -------------------*/

double getdouble_header(char *fichier, int offset)
{
  FILE *fp;
  double lu;
  fp = fopen(fichier, "r");
  if(fp == NULL){
    char e[256];
    sprintf(e, "getdouble_header erreur: Impossible d'ouvrir le fichier <%s>.\n", fichier);
    PUT_ERROR(e);
    fclose(fp);
    return(0.);
  }

  /*fflush(stdin);
    fflush(stdout);*/
  fseek(fp, offset, SEEK_SET);
  fread((void*)&lu, sizeof(double), 1, fp);
  fclose(fp);
  puts("fin getdouble_header");
  return(lu);
}

/*------------------- getfloat_header -------------------*/

float getfloat_header(char *fichier, int offset)
{

  FILE *fp;
  float lu;

  fp = fopen(fichier, "r");
  if(fp == NULL){
    char e[256];
    sprintf(e, "getfloat_header erreur: Impossible d'ouvrir le fichier <%s>.\n",fichier);
    PUT_ERROR(e);
    fclose(fp);
    return(0);
  }

  /*fflush(stdin);
    fflush(stdout);*/
  rewind(fp);
  fseek(fp, offset, 0);
  fread((void *)&lu, sizeof(float), 1, fp);

  fclose(fp);
  return(lu);
}


int GetByte(FILE *fp)
{
  int c;

  if ((c = getc(fp)) == EOF)
    {
      if(MESG_DEBUG) printf( "unexpected end of file");
      return(0);
    }

  return c;
}

int GetWord(FILE *fp)
{
  int c;

  c  = GetByte(fp);
  c |= (GetByte(fp) << 8);

  return c;
}

unsigned long GetDword(FILE *fp)
{
  unsigned long value;

  if(fread(&value, 4, 1, fp) != 1) {
    if(MESG_DEBUG) printf( "[GetDword] read error !\n");
    return(0);
  }

  return(value);
}

/*------------------------------------------------------------------
 * Ecat7_getFromMatrix()
 *
 * args : char * file   :   nom du fichier ECAT
 *        US indMatrix  :   numero de la matrice dans fichier
 *        US offset     :   deplacement relativement au debut
 *                          de la matrice 
 *        US len        :   longueur de la donnee a recuperer
 *        char * data   :   donnee ou stocker le resultat
 *        char * type   :   type du resultat
 *
 * ret:   char *            retourne la valeur placee dans data
 *
 * res : cette fonction lit dans file, dans le header de la matrice
 *       indMatrix une donne de type type au deplacement offset 
 *		 et la stocke dans data avant de la renvoyer
 *
 * Peloutier jean charles 04/2000
 *-----------------------------------------------------------------*/
char* Ecat7_getFromMatrix(char *file, short unsigned int indMatrix, short unsigned int offset, short unsigned int len, char *data, char *type)
{
  FILE *fp;
  int i,dst,isOk=1;
  char *buff=(char*)malloc(10*sizeof(char));
  char *m_int=(char*)malloc(2);
  char *m_float=(char*)malloc(4);
  char *m_long=(char*)malloc(4);
  char tmp;
  fp=fopen(file,"rb");


  if(strcmp(type,"int") && strcmp(type,"long") && strcmp(type,"float") && strcmp(type,"char")) isOk=0;
  /* calcul du deplacement en fonction de l'indice de matrice 
     apres verification de la validite de ce dernier */
  if(isOk){
    if(indMatrix<32){
      fseek(fp,0x206+indMatrix*0x10,SEEK_SET);
      fread(m_int,2,1,fp);
      if(!_is_bigEndian){
	tmp=m_int[0];
	m_int[0]=m_int[1];
	m_int[1]=tmp;
      }
      if(!(dst=*(unsigned short*)m_int))isOk=0;
    }else{
      dst=1;
      for(i=1;i<(indMatrix/32);++i){
	fseek(fp,0x1F6+0x200*dst,SEEK_SET);
	fread(m_int,2,1,fp);
	if(!_is_bigEndian){
	  tmp=m_int[0];
	  m_int[0]=m_int[1];
	  m_int[1]=tmp;
	}
	if(!(dst=*(unsigned short*)m_int))isOk=0;
      }
      fseek(fp,(indMatrix%32)*0x10,SEEK_CUR);
      fread(m_int,2,1,fp);
      if(!_is_bigEndian){
	tmp=m_int[0];
	m_int[0]=m_int[1];
	m_int[1]=tmp;
      }
      if(!(dst=*(unsigned short*)m_int))isOk=0;
    }
    fseek(fp,(dst-1)*0x200+offset,SEEK_SET);
  }
  /* recuperation des donnees */
  if(isOk)
    {
      if(!strcmp(type,"int")){
	fread(m_int,2,1,fp);	
	if(!_is_bigEndian){
	  tmp=m_int[0];
	  m_int[0]=m_int[1];
	  m_int[1]=tmp;
	}
	sprintf(data,"%d",*(unsigned short*)m_int);
      }else if(!strcmp(type,"long")){
	fread(m_long,4,1,fp);	
	if(!_is_bigEndian){
	  tmp=m_int[0];
	  m_long[0]=m_long[1];
	  m_long[1]=tmp;
	  tmp=m_int[2];
	  m_long[2]=m_long[3];
	  m_long[3]=tmp;
	}
	sprintf(data,"%ld",*(long*)m_int);
      }else if(!strcmp(type,"float")){
	fread(m_float,4,1,fp);
	if(!_is_bigEndian){	
	  tmp=m_int[0];
	  m_float[0]=m_float[2];
	  m_float[2]=tmp;
	  tmp=m_float[1];
	  m_float[1]=m_float[3];
	  m_float[3]=tmp;
	}
	sprintf(data,"%f",*(float*)m_float);
      }
    }
  else data=(char*)NULL;

  FREE(m_long);
  FREE(m_float);
  FREE(buff);
  FREE(m_int);
  fclose(fp);
  return data;
}

/****** Acrnema_getGroup() ********************
 *  Renvoie les donn�s correspondant a un groupe
 *  et son sous-groupe dans un fichier de format
 *  ACR-NEMA 2.0 et DICOM
 *  Si le groupe est 7fe0 0010 retourne la totalit�*  du contenu de ce groupe qui est le groupe image.
 *  Burger Laurent    novembre 1999
 *  JP Armspach		Janvier 2001
 ***********************************************/
char *Acrnema_getGroup(char *file, short unsigned int g1, short unsigned int g2, char *data, int format)
{
  FILE *fp;





  long taille,offset,taille_fichier;





  /* Recherche taille et offset des donne du MetaElementGroup */
  taille=Acrnema_getGroupSizeandOffset(file, g1, g2, &offset, format);

  /* debut */
  fp=fopen(file, "rb");
  fseek(fp, 0, SEEK_END);
  taille_fichier=ftell(fp);
  fseek(fp, offset, SEEK_SET);
  
  /* Controle si fin de fichier ; MetaElementGroup n'existe pas*/
  if (offset==taille_fichier) {
    data=(char *)malloc(2*sizeof(char));
    strcpy(data,"\0");
  }   
  else /* Lecture du contenu du MetaElementGroup */
    {
      data=(char *)malloc(taille+1);
      fread(data, taille, 1, fp);
      data[taille]='\0';
    }
  fclose (fp);	

  return (data);
}

/****** Acrnema_getGroupSeq() ********************
 *  Renvoie les donn�s correspondant a un groupe
 *  et son sous-groupe dans un fichier de format
 *  ACR-NEMA 2.0 et DICOM
 *  Si le groupe est 7fe0 0010 retourne la totalit�*  du contenu de ce groupe qui est le groupe image.
 *  Burger Laurent    novembre 1999
 *  JP Armspach		Janvier 2001
 ***********************************************/
char *Acrnema_getGroupSeq(char *file, short unsigned int g1, short unsigned int g2, short unsigned int g12, short unsigned int g22,char *data, int format)
{
  FILE *fp;



  long taille,offset,taille_fichier;






  /* Recherche taille et offset des donne du MetaElementGroup */
  taille=Acrnema_getGroupSizeandOffset(file, g1, g2, &offset, format);

  /* debut */
  fp=fopen(file, "rb");
  fseek(fp, 0, SEEK_END);
  taille_fichier=ftell(fp);
  fseek(fp, offset, SEEK_SET);
  
  /* Controle si fin de fichier ; MetaElementGroup n'existe pas*/
  if (offset==taille_fichier) {
    data=(char *)malloc(2*sizeof(char));
    strcpy(data,"\0");
  }   
  else /* Lecture du contenu du MetaElementGroup */
    {
      data=(char *)malloc(taille+1);
      fread(data, taille, 1, fp);
      data[taille]='\0';
    }
  fclose (fp);	

  return (data);
}

/****** Acrnema_getGroupSizeandOffset() ********************
 *  Renvoie la taille et la position dans le fichiers des donn�s 
 *  correspondant a un groupe
 *  et son sous-groupe dans un fichier de format
 *  ACR-NEMA 2.0 et DICOM
 *  retourne une taille=0 si le groupe n'existe pas ou est vide.
 *  offset est egal a la taille dufichier si le groupe n'existe pas.
 *
 *  JP Armspach		Janvier 2001
 ***********************************************/
long Acrnema_getGroupSizeandOffset(char *file, short unsigned int g1, short unsigned int g2, long *offset, int format)
{
  FILE *fp;

  unsigned short Grp1,Grp2, GrpTmp;
  unsigned int taille;
  int mustConvert;
  char a,b,c,d;
  unsigned long lgrp;

  /* Controle fichier bon format */
  if (format!=ACRNEMA2L_FORMAT && format!=ACRNEMA2B_FORMAT 
      && format!=ACRNEMA1B_FORMAT && format!=ACRNEMA1L_FORMAT
      && format!=DICMBIMP_FORMAT  && format!=DICMLIMP_FORMAT
      && format!=DICMBEXP_FORMAT  && format!=DICMLEXP_FORMAT
      && format!=DICM3DBIMP_FORMAT  && format!=DICM3DLIMP_FORMAT
      && format!=DICM3DBEXP_FORMAT  && format!=DICM3DLEXP_FORMAT)
    {
      printf("Error, bad file format in Acrnema_getGroup\n");
      taille=0; (*offset)=0;
      return (taille);
    }
    	  
  /* Choix de conversion */		  
  if (   ((format==ACRNEMA2B_FORMAT || format==ACRNEMA1B_FORMAT) && !_is_bigEndian) 
	 || ((format==ACRNEMA2L_FORMAT || format==ACRNEMA1L_FORMAT) && _is_bigEndian) 
	 || ((format==DICMBIMP_FORMAT  || format==DICMBEXP_FORMAT)  && !_is_bigEndian) 
	 || ((format==DICMLIMP_FORMAT  || format==DICMLEXP_FORMAT)  && _is_bigEndian) 
	 || ((format==DICM3DBIMP_FORMAT  || format==DICM3DBEXP_FORMAT)  && !_is_bigEndian) 
	 || ((format==DICM3DLIMP_FORMAT  || format==DICM3DLEXP_FORMAT)  && _is_bigEndian)) 
    mustConvert=1;
  else 
    mustConvert=0;

  fp=fopen(file, "rb");
      
  /* Recherche de la taille du MetaElementGroupLength puis saut*/
  lgrp=0;
  a=b=c=d=0;
  fread (&d, 1, 1, fp);
  fread (&c, 1, 1, fp);
  fread (&b, 1, 1, fp);
  fread (&a, 1, 1, fp);
  lgrp=(d<<24)+(c<<16)+(b<<8)+a;
  while (!feof(fp) && lgrp!=0x02000000)
    {
      d=c;
      c=b;
      b=a;
      a=0;
      fread(&a, 1, 1, fp);
      lgrp=(d<<24)+(c<<16)+(b<<8)+a;
    }
  if (  format==ACRNEMA1B_FORMAT || format==ACRNEMA1L_FORMAT
	|| format==DICMBEXP_FORMAT  || format==DICMLEXP_FORMAT 
	|| format==DICMBIMP_FORMAT  || format==DICMLIMP_FORMAT 
	|| format==DICM3DBEXP_FORMAT  || format==DICM3DLEXP_FORMAT 
	|| format==DICM3DBIMP_FORMAT  || format==DICM3DLIMP_FORMAT) 
    {
      taille = 0;
      fseek (fp, 2, SEEK_CUR);
      fread(&taille, 2, 1, fp);
    }
  else fread(&taille, 4, 1, fp);
  if (mustConvert)  taille=longEndianConversion(taille);
  if (taille ==4) fread(&taille, 4, 1, fp);
  if (mustConvert)  taille=longEndianConversion(taille);
  fseek(fp, taille, SEEK_CUR);

  /* Lecture des infos */
  Grp1=Grp2=0;taille=0; 					  
  while (!feof(fp) && (Grp1!=g1 ||Grp2!=g2))
    {
      if (  format==ACRNEMA1B_FORMAT || format==ACRNEMA1L_FORMAT
	    || format==DICMBEXP_FORMAT  || format==DICMLEXP_FORMAT 
	    || format==DICMBIMP_FORMAT  || format==DICMLIMP_FORMAT 
	    || format==DICM3DBEXP_FORMAT  || format==DICM3DLEXP_FORMAT 
	    || format==DICM3DBIMP_FORMAT  || format==DICM3DLIMP_FORMAT) 
	{
	  fseek(fp, taille, SEEK_CUR);
	  fread(&Grp1, 2, 1, fp);
	  if (mustConvert) Grp1=shEndianConversion(Grp1);
	  fread(&Grp2, 2, 1, fp);
	  if (mustConvert) Grp2=shEndianConversion(Grp2);
    			  
	  taille=0;
	  /* les 14 lignes de codes suivantes ont ete realise dans un souci de performance du programme */
	  if ( (Grp1<g1) && (Grp2==0x0000) && (Grp1!=0x7fe0)) 
	    {
	      if (format==DICMBEXP_FORMAT  || format==DICMLEXP_FORMAT
		  || format==DICM3DBEXP_FORMAT  || format==DICM3DLEXP_FORMAT)
		{
		  fread(&taille, 2, 1, fp);		    
		  fread(&taille, 2, 1, fp);
		}
	      else
		{
		  fread(&taille, 4, 1, fp);
		} 
	      if (mustConvert) taille=longEndianConversion(taille);
	      if (taille==4) fread(&taille, 4, 1, fp);
	      /*	Chaque tag dont le sous_groupe est 0x0000 est un GroupLength, c'est a dire
			un tag qui defini la longueur du group concerne ; ex :
			- (0x0010,0x0000) 64   4  PatientGroupLength
			donc si on cherche un tag dont le groupe est superieur a celui rencontre,
			on ne traite pas chaque sous groupe du GroupLength (puisque le tag
			recherche ne s'y trouve pas) mais on saute le groupe en utilisant la longueur 
			donner par le GroupLength.
			Ce traitement ne doit pas etre appliquer sur le PixelDataGroupLength afin
			de pouvoir avoir acces au data images sans probleme.
			
			il y a deux possibilite, soit le format est implicit, soit explicit : 
			
			1 - En explicit le premier fread nous permet d'obtenir la VR du DataElement (toujours UL 
			dont la longueur est de 2 bytes fixed et la valeur 4.
			Le second fread permet d'obtenir la longueur en bytes du champs qui contient la 
			valeur rechercher. Si la taille est bien egale a 4, on lit le fichier sur 4 bytes 
			pour obtenir la longueur du GroupLength grace a l'endian conversion realisee plus bas.
			
			2 - En implicit, il n'y a pas de VR, on lit directement la taille du champ qui contient la valeur 
			et qui doit etre de 4 bytes. si la taille est bien egale a quatre, on avance toujours de 4 bytes dans 
			le fichier et de la meme maniere quand explicit, on converti la valeur quelques lignes plus bas
	      */	
	    }
	  else
	    {
	      /* Pour le format Dicom en Explicit */
	      if ( format ==DICMBEXP_FORMAT || format == DICMLEXP_FORMAT ||
		   format ==DICM3DBEXP_FORMAT || format == DICM3DLEXP_FORMAT) 
		{    
		  /*    les 3 tag suivants (groupe = Grp1, sous_groupe = Grp2) posent un probleme en explicit
			puisqu'ils n'ont pas de VR. Ils representent respectivement 'Item Delimitation Item' , 
			'Sequence Delimitation Item' et 'Data Element Item'.
			Lorsque l'on rencontre ces tag, ils doivent etre traites comme en implicit.
			On lit donc juste le champ ou se trouve la taille de l'element, valeur qui sera convertit
			et testee plus bas. (la taille finale renvoye devra etre de 0, c'est le cas pour les deux 
			derniers tag alors que celle du premier est de : 4294967295... c'est a dire la longueur
			d'une 'sequence with undefined length' et qui sera mise a 0 plus bas).
		  */  
		  if ( (Grp1 == 0xfffe && Grp2 == 0xe000) || (Grp1 == 0xfffe && Grp2 == 0xe00d) 
		       || (Grp1 == 0xfffe && Grp2 == 0xe0dd)  ) {
		    fread(&taille, 4, 1, fp);
		  }
		  else
		    {
		      fread (&GrpTmp,2 ,1, fp);
		      /* OB, OW, SQ, UN sont des tags speciaux et qui doivent etre traites differement
			 des autres tags. Voir documentation DICOM.
		      */
		      if (   GrpTmp == 0x4f42 || GrpTmp == 0x424f  /* OB */
			     || GrpTmp == 0x4f57 || GrpTmp == 0x574f  /* OW */ 
			     || GrpTmp == 0x5351 || GrpTmp == 0x5153  /* SQ */ 
			     || GrpTmp == 0x554e || GrpTmp == 0x4e55) /* UN */
			{ 
    		   	  fread(&taille, 2, 1, fp);
			  fread(&taille, 4, 1, fp);
			}
		      else
			fread(&taille, 2, 1, fp);  
		    } 
		}   /* entre autre pour DICMBIMP_FORMAT et DICMLIMP_FORMAT et les implicit Dicom 3D FORMAT */
	      else 
		fread(&taille, 4, 1, fp);
	    }	
	  if (mustConvert) taille=longEndianConversion(taille);
	  if (taille == UINT_MAX)
	    taille=0;
	}
    }

  /* Controle si fin de fichier ; MetaElementGroup n'existe pas*/
  if (feof(fp)) {
    taille=0;
    (*offset)=ftell(fp);
  }
  else
    (*offset)=ftell(fp);
		
  fclose (fp);	

  return (taille);
}

/****** Acrnema_getGroupHeader() ********************
 *Renvoie les donn�s correspondant a un groupe
 * et son sous-groupe du groupe header dans un fichier de format
 * ACR-NEMA 2.0 et DICOM
 *  Burger Laurent       novembre 1999
 *  jean-Paul Armspach	Janvier 2001
 ***********************************************/
char *Acrnema_getGroupHeader(char *file, short unsigned int g1, short unsigned int g2, char *data, int format)
{
  FILE *fp;

  unsigned short Grp1,Grp2, GrpTmp;
  unsigned int taille;
  int mustConvert;
  char a,b,c,d;
  unsigned long lgrp;

  if (format!=ACRNEMA2L_FORMAT && format!=ACRNEMA2B_FORMAT 
      && format!=ACRNEMA1B_FORMAT && format!=ACRNEMA1L_FORMAT
      && format!=DICMBIMP_FORMAT  && format!=DICMLIMP_FORMAT
      && format!=DICMBEXP_FORMAT  && format!=DICMLEXP_FORMAT
      && format!=DICM3DBIMP_FORMAT  && format!=DICM3DLIMP_FORMAT
      && format!=DICM3DBEXP_FORMAT  && format!=DICM3DLEXP_FORMAT)
    {
      printf("Error, bad file format \n");
      data=(char *)NULL;
      return (data);
    }
				
		
  if (   ((format==ACRNEMA2B_FORMAT || format==ACRNEMA1B_FORMAT) && !_is_bigEndian) 
	 || ((format==ACRNEMA2L_FORMAT || format==ACRNEMA1L_FORMAT) && _is_bigEndian) 
	 || ((format==DICMBIMP_FORMAT  || format==DICMBEXP_FORMAT)  && !_is_bigEndian) 
	 || ((format==DICMLIMP_FORMAT  || format==DICMLEXP_FORMAT)  && _is_bigEndian) 
	 || ((format==DICM3DBIMP_FORMAT  || format==DICM3DBEXP_FORMAT)  && !_is_bigEndian) 
	 || ((format==DICM3DLIMP_FORMAT  || format==DICM3DLEXP_FORMAT)  && _is_bigEndian)) 




    mustConvert=1;
  else 
    mustConvert=0;


  fp=fopen(file, "rb");			
  /* Recherche de la taille du MetaElementGroupLength puis saut*/
  lgrp=0;
  a=b=c=d=0;
  fread (&d, 1, 1, fp);
  fread (&c, 1, 1, fp);
  fread (&b, 1, 1, fp);
  fread (&a, 1, 1, fp);
  lgrp=(d<<24)+(c<<16)+(b<<8)+a;
  while (!feof(fp) && lgrp!=0x02000000/*0x08001000*/ )
    {
      d=c;
      c=b;
      b=a;
      a=0;
      fread(&a, 1, 1, fp);
      lgrp=(d<<24)+(c<<16)+(b<<8)+a;
    }
  taille = 0;
  fseek (fp, 2, SEEK_CUR);
  fread(&taille, 2, 1, fp);
  if (mustConvert)  taille=longEndianConversion(taille);
 		
  /* Lecture des infos */
  Grp1=Grp2=0;						
  while (!feof(fp) && (Grp1!=g1 ||Grp2!=g2))
    {
      if (  format==ACRNEMA1B_FORMAT || format==ACRNEMA1L_FORMAT
	    || format==DICMBEXP_FORMAT  || format==DICMLEXP_FORMAT 
	    || format==DICMBIMP_FORMAT  || format==DICMLIMP_FORMAT 
	    || format==DICM3DBEXP_FORMAT  || format==DICM3DLEXP_FORMAT 
	    || format==DICM3DBIMP_FORMAT  || format==DICM3DLIMP_FORMAT) 
	{
	  fseek(fp, taille, SEEK_CUR);
	  fread(&Grp1, 2, 1, fp);
	  if (mustConvert) Grp1=shEndianConversion(Grp1);
	  fread(&Grp2, 2, 1, fp);
	  if (mustConvert) Grp2=shEndianConversion(Grp2);
						
	  taille=0;
	  if ( format ==DICMBEXP_FORMAT || format == DICMLEXP_FORMAT ||
	       format ==DICM3DBEXP_FORMAT || format == DICM3DLEXP_FORMAT) { 
	    fread (&GrpTmp,2 ,1, fp);
	    /*printf (" %2s \n",&GrpTmp);*/
	    if (   GrpTmp == 0x4f42 || GrpTmp == 0x424f  /* OB */
		   || GrpTmp == 0x4f57 || GrpTmp == 0x574f  /* OW */ 
		   || GrpTmp == 0x5351 || GrpTmp == 0x5153  /* SQ */ 
		   || GrpTmp == 0x554e || GrpTmp == 0x4e55) /* UN */{ 
	      fread(&taille, 2, 1, fp);
	      fread(&taille, 4, 1, fp);
	    }
	    else fread(&taille, 2, 1, fp);  
	  } /*entre autre pour DICMBIMP_FORMAT et DICMLIMP_FORMAT*/
	  else fread(&taille, 4, 1, fp);
	  if (mustConvert) taille=longEndianConversion(taille);
	}
    }

  if (feof(fp)) data=(char *)NULL;
		
  else
    {
      data=(char *)malloc((taille+1) * sizeof(char));
      fread(data, taille, 1, fp);
      data[taille]='\0';
      /* puts("valeur lue :");
	 puts(data);*/
    }
  fclose (fp);	

  return (data);

}
#ifdef __GTK_GUI

/*******************************************
 ** --  get_rep() ----------------------------
 **
 **  Permet de rechercher le contenu d'un repertoire (rep)
 **  retour dans le tableau **quest le resultat.
 **  max : nombre de ligne dans quest
 **  ATTENTION : Ce tableau a ete cree par CALLOC donc free
 **  
 ********************************************/
char   **get_rep(char *rep, int *max)
{
  char **quest;
  int num_entries,len=80;
  int i,k=0;
  char **namelist, **list;

  quest=/*CALLOC(len,char)*/CALLOC(len/*/4*/,char *);
  if ((num_entries = gtkmmscandir(rep,&namelist,F_SORT)) <0) {
    quest[0]=CALLOC(len,char);
    strcpy(quest[0],NULL);
    return((char **) quest);
  }
  printf("Number of entries is %d\n",num_entries);
  if (num_entries) {
    if (MESG_DEBUG) printf("Entries are :");
    for (i=0,list=namelist; i<num_entries ;i++) 
      {
	if (MESG_DEBUG) 
	  printf(" %s\n", *list);
	if (!(!strcmp(*list,".")
	      ||!strcmp(*list,"..") 
	      ||!strcmp(*list,"subject")) ) 
	  {
	    quest[k]=CALLOC(len,char);
	    strcpy(quest[k++],*list);
	  }
	FREE(*list);
	*list++;
      }
    quest[k]=CALLOC(len,char);
    strcpy(quest[k++],"\0");
    *max=k;
    FREE(namelist);
  }
  return((char **) quest);
}

/*******************************************
 ** --  incfic_img(name,step,numfirst) ----------------------------
 **
 **  Permet d'incrementer l'extension d'un  nom de fichier
 **  name : nom du fichier
 **  step : pas de l'increment
 **  numfirst : si fichier a image multiple, increment de numfirst
 ** 
 **  Retourne le numero de l'image dans le fichier
 ********************************************/

int incfic_img(char *name, int step, int *numfirst)
{
  int i, format;

  format = is_mri(name,NULL);

  if(format == PCX_FORMAT
     || format == SMIS_FORMAT
     || format == ACRNEMA2B_FORMAT
     || format == ACRNEMA2L_FORMAT) {
    for(i=0;i<step;i++)	incname_img(name);
  }
  else {
    if(format ==PARAVISION_FORMAT
       || format ==INTERFILE_FORMAT
       || format ==IPB_2D_FORMAT
       || format ==IPB_3D_FORMAT)
      *numfirst = *numfirst+step;
    else incfic_I00X(name,step);
  }
  return(*numfirst);
}
/*******************************************
 ** --  incname_img(name) ----------------------------
 **
 **  Permet d'incrementer le nom du fichier de l'image
 **  (partie precedant l'extension)
 ********************************************/
char   *incname_img(char *name)
{  char *form,*nametemp;
 char *tab[10];
 char arg1[FILE_LEN];
 char *delims=".";
 int count,len,i;


 strcpy(arg1,name);
 nametemp=arg1;
 len=strlen(nametemp);
 if(len<=0) return(NULL);
 for (i=0;i<10;i++)
   tab[i]=/*CALLOC(FILE_LEN,char *)*/CALLOC(FILE_LEN*4,char);

 count=0;
 form = strtok(nametemp,delims); 
 strcpy(tab[count],form);
 while ((form=strtok((char *)NULL,delims)) != NULL)
   {
     count=count+1;
     if (count<=11)
       strcpy(tab[count],form);
   }
 if ((count > 11)&&(MESG_DEBUG))
   printf(" ERROR   dans incname_img depassement tableau");
  
 incfic(tab[count-1],1);
 strcpy(name,tab[0]);
 for (i=1;i<=count;i++)
   {
     strcat(name,delims);
     strcat(name,tab[i]);
   }

 for (i=0;i<10;i++) FREE(tab[i]);
 return(NULL);
}

#endif /* __GTK_GUI */

/* --fread_long() ---------------------------------------------------------    
**
**   Permet la lecture des images independemment du type de machine
**
**/

long fread_long (long int *t, int nitem, FILE *fp)
{
  int i;
  long size_tot=0,size=1;
  BYTE *buff;
  buff = (BYTE *) malloc (4);

  for ( i=0; (i<nitem) && (size==1) ; i++ )
  {
      size=fread ( buff, (size_t) 4, (size_t) 1, fp );
      t[i]=( (((unsigned long) buff[0])<<24) + (((unsigned long) buff[1])<<16) + (((unsigned long) buff[2])<<8) +
	     buff[3] );    
      size_tot+=size;
    }
  free (buff);
  return(size_tot);
}

long fread_long2 (long int *t, int nitem, FILE *fp)
{
//  int i;
  long size_tot=0;//,size=1;
//  BYTE *buff;
//  buff = (BYTE *) malloc (4);

/*  for ( i=0; (i<nitem) && (size==1) ; i++ )
  {
      size=fread ( buff, (size_t) 4, (size_t) 1, fp );
      t[i]=( (((unsigned long) buff[0])<<24) + (((unsigned long) buff[1])<<16) + (((unsigned long) buff[2])<<8) +
	     buff[3] );
      size_tot+=size;
    }*/
  size_tot=fread (t,(size_t) sizeof(long), (size_t) nitem, fp);
//  free (buff);
  return(size_tot);
}

/* --fread_short() ---------------------------------------------------------    
**
**  Permet la lecture des images independemment du type de machine
**
**/
long fread_short (short *t, int nitem, FILE *fp)
{
  int i;
  long size_tot=0,size=1;
  BYTE *buff;
  buff = (BYTE *) malloc (2);

  for ( i=0; (i<nitem) && (size==1) ; i++ )
    {
      size=fread ( buff, (size_t) 2, (size_t) 1, fp );
      t[i]=(((( unsigned short) buff[0])<<8) + buff[1] );
      size_tot+=size;
    }
  free (buff);
  return(size_tot);
}

long fread_short2 (short *t, int nitem, FILE *fp)
{
//  int i;
  long size_tot=0;//,size=1;
//  BYTE *buff;
//  buff = (BYTE *) malloc (2);

/*  for ( i=0; (i<nitem) && (size==1) ; i++ )
    {
      size=fread ( buff, (size_t) 2, (size_t) 1, fp );
      t[i]=(((( unsigned short) buff[0])<<8) + buff[1] );
      size_tot+=size;
    }*/
  size_tot=fread (t,(size_t) sizeof(short), (size_t) nitem, fp);

//  free (buff);
  return(size_tot);
}


/* --fwrite_long() ---------------------------------------------------------    
**
**i  Permet l'ecriture des images independemment du type de machine
**
**/
long fwrite_long (long int *t, int nitem, FILE *fp)
{
  int i;
  long size_tot=0,size=4;
  char buffstr[2];
  union trans 
  {
    BYTE buff[4];
    long tt;
  }tr;
  for ( i=0; (i<nitem) && (size==4) ; i++ )
    {
      tr.tt=t[i];
      buffstr[1]='\0';
      buffstr[0]=tr.buff[3];
      size =fwrite ( buffstr, (size_t) 1, (size_t) 1, fp );
      buffstr[0]=tr.buff[2];
      size+=fwrite ( buffstr, (size_t) 1, (size_t) 1, fp );
      buffstr[0]=tr.buff[1];
      size+=fwrite ( buffstr, (size_t) 1, (size_t) 1, fp );
      buffstr[0]=tr.buff[0];
      size+=fwrite ( buffstr, (size_t) 1, (size_t) 1, fp );
      size_tot+=size;
    }
  return(size_tot);
}

/* -- Make_AtmFile() --------------------------------------------------
**
*/
void    Make_AtmFile(char *file, int nb_img, char **pictures)
{ 
  FILE *fp;
  int i,errno;

  /* Ist  PART: Check if roi exists. */
  errno=0;
  if((fp=fopen(file,"r"))==NULL) {
    if(errno==2 || errno==0) {  /* Not such file or directory */

      /* -> Save */         }
    else { /* Error appended ! */
      char s[256];
      _imagix_error=IMX_FILE_OPEN_ERROR;
      _imagix_err_ext=errno;
      sprintf(s,"Error occured:\n  Cannot create atm parameter file >%s<\n"
	      "  Errno:%d\n",file,errno);
      PUT_ERROR(s);
      fclose(fp);
      return;
    } 
  }
  else { /* File exists */
    int err; 
    char s[256];
    sprintf(s,"Delete old file ?\n  (%s)\n",file);
    if(XGET_EXIST(s, &err)==0) {

      sprintf(s,"Cannot create atm parameter file >%s<\n",file);
      PUT_ERROR(s);
      fclose(fp);
      return;
    } 
  }
  fclose(fp);

  if((fp=fopen(file,"w"))==NULL) { 
    char s[256],tmp[256];
    _imagix_error=IMX_FILE_OPEN_ERROR;
    _imagix_err_ext=errno;
    sprintf(tmp,"Error occured:\n* Cannot create atm parameter file >%s<\n"
            ,file);
    sprintf(s,"%s  errno=%d\n\n",tmp,errno);
    PUT_ERROR(s);
    fclose(fp);
    return;
  }

  for(i=0;i<nb_img;i++) {
    if(pictures[i]!=NULL && pictures!=NULL)
      fprintf(fp,"%s\n",pictures[i]);
  }

  fclose(fp);
}

/* -- Make_AtmFileII() ------------------------------------------------
**
*/
void    Make_AtmFileII(char *file, char *wildcard_file)
{ 
  FILE *fp;
  char cmd[FILE_LEN*2+255];

  /* Ist  PART: Check if roi exists. */
  errno=0;
  if((fp=fopen(file,"r"))==NULL) {
    if(errno==2 || errno==0) {  /* Not such file or directory */

      /* -> Save */         }
    else { /* Error appended ! */
      char s[256];
      _imagix_error=IMX_FILE_OPEN_ERROR;
      _imagix_err_ext=errno;
      sprintf(s,"Error occured:\n  Cannot create atm parameter file >%s<\n"
	      "  Errno:%d\n",file,errno);
      PUT_ERROR(s);
      fclose(fp);
      return;
    } 
  }
  else { /* File exists */
    int err; 
    char s[256];
    sprintf(s,"Delete old file ?\n  (%s)\n",file);
    if(XGET_EXIST(s, &err)==0) {

      sprintf(s,"Cannot create atm parameter file >%s<\n",file);
      PUT_ERROR(s);
      fclose(fp);
      return;
    } 
  }
  fclose(fp);

  if((fp=fopen(file,"w"))==NULL) { 
    char s[256],tmp[256];
    _imagix_error=IMX_FILE_OPEN_ERROR;
    _imagix_err_ext=errno;
    sprintf(tmp,"Error occured:\n* Cannot create atm parameter file >%s<\n"
            ,file);
    sprintf(s,"%s  errno=%d\n\n",tmp,errno);
    PUT_ERROR(s);
    fclose(fp);
    return;
  }
  fclose(fp);

  sprintf(cmd,"ls -1 %s >> %s",wildcard_file,file);
  system(cmd);
}

/*******************************************
 ** --  OpenTempfile() ----------------------------
 **
 **  Ouverture d'un fichier temporaire
 **  Retourne l'identification du graphe.
 ********************************************/
int     OpenTempfile(char *nametemp, FILE **fp)
{
#ifndef WIN32
  int id;
  long pid;

  id=10;
  pid=getpid();
  sprintf(nametemp,"%s%ld%s%d","/tmp/imagix_",pid,"_",id);
  if((*fp=fopen(nametemp,"w"))==NULL)
    {
#ifdef __GTK_GUI
      if(_dialog_messages_status&0x00000002l)
	PUT_ERROR("Open file error");
#endif /* __GTK_GUI */
      return(-1);
    }
#endif // WIN32
  return(0);
}


/*******************************************
 ** --  CloseTempfile() ----------------------------
 **
 **  Ferme le fichier temporaire
 ********************************************/
int     CloseTempfile(FILE *fp)
{

  fclose (fp); 

  return(0);
}

/*******************************************
 ** --  WriteTempfile() ----------------------------
 **
 **  Permet d'ecrire dans le fichier temporaire.
 ********************************************/
int     WriteTempfile(FILE *fp, int type, char *str, int i, int f)
{


  switch(type) {

  case '1' : /* float */
    sprintf(str,"%f",(float)f);
    break;

  case '2' : /* int */
    sprintf(str,"%d",f);
    break;

  }

  fprintf(fp,"%s\n",str);

  return(0);
}

#ifdef __GTK_GUI
/******************************************************************************
 ** -- nOpen_Cmd() -------------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July 22nd, 1999
 **
 **  Opens multiple 2d pictures.
 **
 ******************************************************************************/

void nOpen_Cmd(void)
{/* Load several files.             */
  /*   1st, Call file selection box and:    */
  /*   2nd, Make severals files.         */
  char *file,**files;
            
  char name[FILE_LEN];
	    
  int pos,i,max,step=0,numfirst,e,num;

  file=GET_FILE_IMG(NULL, &numfirst, &e);
  if (e) return;
  if(file==NULL || strlen(file)<=0)
    PUT_WARN("Sorry, cannot read picture !\n");
  else {
    switch(is_mri(file,NULL))
      {
      case BRUKER_FORMAT :
	pos=GET_WND(NULL);
	max= GET_INT(TEXT0004, 0, &e);
	if (!e)
	  step= GET_INT(TEXT0024, 0, &e);
	if (e) break;
	/*  traitement cas particulier
	    Si max(nombre)=0 alors affichage d'une seule image
	    qui sera image+pas
	*/
	if (max==0)
	  {
	    strcpy(name,file);
	    incfic(name,step);
	    strcpy(file,name);
	    max  = 1;
	    step = 1;
	  }
	files=Make_files(file,step,max);
	if(load_pictures(pos,MINI(max,MAX_PICTURES-pos+1)
			 ,files)==-1)
	  PUT_ERROR("Cannot read picture.\n");
	for(i=0;i<max;i++)
	  FREE(files[i]);
	FREE(files);
	break;

      case PARAVISION_FORMAT :
	pos=GET_WND(NULL);
	max= GET_INT(TEXT0004, 0, &e);
	if (e) break;
	step= GET_INT(TEXT0024, 0, &e);
	if (e) break;
	for (i=0;i<max;i++)
	  {
	    load_picture(pos+i,file,numfirst);
	    numfirst=numfirst+step;
	  }
	break;
      case SIEMENS_FORMAT :
	pos=GET_WND(NULL);
	max= GET_INT(TEXT0004, 0, &e);
	if (!e)
	  step= GET_INT(TEXT0024, 0, &e);
	if (e) break;
	/*  traitement cas particulier
	    Si max(nombre)=0 alors affichage d'une seule image
	    qui sera image+pas
	*/
	if (max==0)
	  {
	    strcpy(name,file);
	    incfic(name,step);
	    strcpy(file,name);
	    max=1;
	    step=1;
	  }
	files=Make_files(file,step,max);
	if(load_pictures(pos,MINI(max,MAX_PICTURES-pos+1)
			 ,files)==-1)
	  PUT_ERROR("Cannot read picture.\n");
	for(i=0;i<max;i++)
	  FREE(files[i]);
	FREE(files);
	break;

      case VANDERBILT_FORMAT :
	pos=GET_WND(NULL);
	max= GET_INT(TEXT0004, 0, &e);
	if (!e)
	  step= GET_INT(TEXT0024, 0, &e);
	if (e) break;
	/*  traitement cas particulier
	    Si max(nombre)=0 alors affichage d'une seule image
	    qui sera image+pas
	*/
	if (max==0)
	  {
	    strcpy(name,file);
	    incfic(name,step);
	    strcpy(file,name);
	    max=1;
	    step=1;
	  }
	files=Make_files(file,step,max);
	if(load_pictures(pos,MINI(max,MAX_PICTURES-pos+1)
			 ,files)==-1)
	  PUT_ERROR("Cannot read picture.\n");
	for(i=0;i<max;i++)
	  FREE(files[i]);
	FREE(files);
	break;


      case INTERFILE_FORMAT :
	pos=GET_WND(NULL);
	max= GET_INT(TEXT0004, 0, &e);
	if (e) break;
	step= GET_INT(TEXT0024, 0, &e);
	if (e) break;
	/*  numfirst=GET_INT(
	    "Numero first image ?",NULL,&e);
	    if (e) break;
	*/   for (i=0;i<max;i++)
	  {
	    load_picture(pos+i,file,numfirst);
	    numfirst=numfirst+step;
	  }
	break;

      case RAW4_FORMAT :
	pos=GET_WND(NULL);
	max= GET_INT(TEXT0004, 0, &e);
	if (e) break;
	step= GET_INT(TEXT0024, 0, &e);
	if (e) break;
	/*  traitement cas particulier
	    Si max(nombre)=0 alors affichage d'une seule image
	    qui sera image+pas
	*/
	if (max==0)
	  {
	    strcpy(name,file);
	    incfic(name,step);
	    strcpy(file,name);
	    max=1;
	    step=1;
	  }
	files=Make_files(file,step,max);
	if(load_pictures(pos,MINI(max,MAX_PICTURES-pos+1)
			 ,files)==-1)
	  PUT_ERROR("Cannot read picture.\n");
	for(i=0;i<max;i++)
	  FREE(files[i]);
	FREE(files);
	break;
      case IPB_2D_FORMAT:
      case IPB_3D_FORMAT:
	pos=GET_WND(NULL);
	max= GET_INT(TEXT0004, 0, &e);
	if (e) break;
	step= GET_INT(TEXT0024, 0, &e);
	if (e) break;
	for (i=0;i<max;i++)
	  {
	    load_picture(pos+i,file,numfirst);
	    numfirst=numfirst+step;
	  }
	break;
      case PCX_FORMAT :
	pos=GET_WND(NULL);
	max= GET_INT(TEXT0004, 0, &e);
	if (e) break;
	step= GET_INT(TEXT0024, 0, &e);
	if (e) break;
	/*  traitement cas particulier
	    Si max(nombre)=0 alors affichage d'une seule image
	    qui sera image+pas
	*/
	if (max==0)
	  {
	    strcpy(name,file);
	    /*protoize:???*/
	    incname_img(name/*,step, numfirst*/);
	    strcpy(file,name);
	    max=1;
	    step=1;
	  }
	files=Make_files(file,step,max);
	if(load_pictures(pos,MINI(max,MAX_PICTURES-pos+1)
			 ,files)==-1)
	  PUT_ERROR("Cannot read picture.\n");
	for(i=0;i<max;i++)
	  FREE(files[i]);
	FREE(files);
	break;

      case SMIS_FORMAT:
	pos = GET_WND( NULL);
	max =  GET_INT( TEXT0004, 0, &e);
	if (e) break;
	step = GET_INT( TEXT0024, 0, &e);
	if (e) break;

	if (max == 0) {
	  strcpy(name, file);
	  incname_img(name);
	  strcpy(file, name);
	  max = 1;
	  step = 1;
	}

	for(i=0; i<max; i++) {
	  if(load_picture(pos+i, file, 1) == -1) {
	    PUT_ERROR("Cannot read picture.\n");
	    break;
	  }
	  strcpy(name, file);
	  incfic_img(name, step,&numfirst);
	  strcpy(file, name);
	}
	break;

      case ACRNEMA2L_FORMAT:
	pos = GET_WND( NULL);
	max =  GET_INT( TEXT0004, 0, &e);
	if (e) break;
	step = GET_INT( TEXT0024, 0, &e);
	if (e) break;

	strcpy(name, file);
                   
	for(i=0; i<max; i++) {
	  if(load_picture(pos+i, name, 1) == -1) {
	    PUT_ERROR("Cannot read picture.\n");
	    break;
	  }
	  incfic_img(name, step, &numfirst);
	}
	break;

      case ACRNEMA2B_FORMAT:
	pos = GET_WND( NULL);
	max =  GET_INT( TEXT0004, 0, &e);
	if (e) break;
	step = GET_INT( TEXT0024, 0, &e);
	if (e) break;

	strcpy(name, file);
                   
	for(i=0; i<max; i++) {
	  if(load_picture(pos+i, name, 1) == -1) {
	    PUT_ERROR("Cannot read picture.\n");
	    break;
	  }
	  incfic_img(name, step, &numfirst);
	  /*incfic(name, step);
	    strcpy(file, name);*/
	}
				
	break;            

      case ACRNEMA1B_FORMAT:
	pos = GET_WND( NULL);
	max =  GET_INT( TEXT0004, 0, &e);
	if (e) break;
	step = GET_INT( TEXT0024, 0, &e);
	if (e) break;

	strcpy(name, file);
                   
	for(i=0; i<max; i++) {
	  if(load_picture(pos+i, name, 1) == -1) {
	    PUT_ERROR("Cannot read picture.\n");
	    break;
	  }
	  incfic_img(name, step, &numfirst);
	  /*incfic(name, step);
	    strcpy(file, name);*/
	}
				
	break;      
					
      case ACRNEMA1L_FORMAT:
	pos = GET_WND( NULL);
	max =  GET_INT( TEXT0004, 0, &e);
	if (e) break;
	step = GET_INT( TEXT0024, 0, &e);
	if (e) break;

	strcpy(name, file);
                   
	for(i=0; i<max; i++) {
	  if(load_picture(pos+i, name, 1) == -1) {
	    PUT_ERROR("Cannot read picture.\n");
	    break;
	  }
	  incfic_img(name, step, &numfirst);
	  /*incfic(name, step);
	    strcpy(file, name);*/
	}
				
	break;            
      case DICMBIMP_FORMAT :
      case DICMBEXP_FORMAT :
      case DICMLIMP_FORMAT :
      case DICMLEXP_FORMAT :
	pos = GET_WND( NULL);
	max =  GET_INT( TEXT0004, 0, &e);
	if (e) break;
	step = GET_INT( TEXT0024, 0, &e);
	if (e) break;

	/*                    strcpy(name, file);
			      i=strlen(name);
			      memcpy(ext,name+i-4,4);
			      ext[4]='\0';
			      i-=5;
			      a=0;
			      while(name[i]>='0' && name[i]<='9'){
			      a++;i--;
			      }
			      i++;
			      for(j=i;j<a+i;j++)tcur[j-i]=name[j];
			      tcur[a]='\0';
			      cur=atoi(tcur);
			      memcpy (pref,name,i);
			      pref[i]='\0';
			      sprintf(m_name,"%s%d%s",pref,cur,ext);

			      for(i=0; i<max; i++) {
			      if(load_picture(pos+i, m_name, 1) == -1) {
			      PUT_ERROR("Cannot read picture.\n");
			      break;
			      }
			      cur+=step;
			      sprintf(m_name,"%s%d%s",pref,cur,ext);
	*/
	strcpy(name, file);
                   
	for(i=0; i<max; i++) {
	  if(load_picture(pos+i, name, 1) == -1) {
	    PUT_ERROR("Cannot read picture.\n");
	    break;
	  }
	  incfic_img(name, step, &numfirst);
	}
	break;			
					

      case ECAT7_FORMAT:
	pos = GET_WND( NULL);
	num= GET_INT(TEXT0183, 0, &e);
	if (e) break;
	max =  GET_INT( TEXT0004, 0, &e);
	if (e) break;
	step = GET_INT( TEXT0024, 0, &e);
	if (e) break;

	for(i=0;i<max;++i){
	  if(load_picture(pos+i, file, num) == -1) {
	    PUT_ERROR("Cannot read picture.\n");
	    break;
	  }
	  num+=step;
	}
	break;
      default : 
	break;
      }
  }
}

/******************************************************************************
 ** -- Open_Norm() -------------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July 22nd, 1999
 **
 **  Opens a 2d picture either with or without normalisation.
 **
 ******************************************************************************/

void Open_Norm(void)
{
  int n,e,norm;
  char *qcm[] = {"2x2", "3x3", "4x4","5x5", "6x6", "More ...", "" };
  char *qcm2[] = { "Without normalisation", "With normalisation", "" };

  n = 2 + GETV_QCM( TEXT0043,(char **)qcm);
  if (n==7)
    {
      n = (GET_INT(TEXT0118, 0,&e));
      if (e || n<1) return;
    }
  
  norm=GETV_QCM(TEXT0043,(char **)qcm2);
  
  Lecture_1img_ds_1Win(n,norm);
}

/******************************************************************************
 ** -- nOpen_Norm() ------------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July 22nd, 1999
 **
 **  Opens multiple 2d pictures either with or without normalisation.
 **
 ******************************************************************************/

void nOpen_Norm(void)
{
  char *quest[7];
  int i,n,e,norm;

  for(i=0;i<7;i++)
    /*protoize:???*/
    quest[i]=CALLOC(80,char /*  * */);
  strcpy(quest[0],"2x2");
  strcpy(quest[1],"3x3");
  strcpy(quest[2],"4x4");
  strcpy(quest[3],"5x5");
  strcpy(quest[4],"6x6");
  strcpy(quest[5],"More ...");
  strcpy(quest[6],"\0");

  n=2+GETV_QCM(TEXT0043,(char **)quest);
  if (n==7)
    {
      n = GET_INT(TEXT0118, 0,&e);
      if (e || n<1) return;
    }

  for(i=0;i<3;i++)
    /*protoize:???*/
    quest[i]=CALLOC(80,char /*  * */);
  strcpy(quest[0],"Without normalisation");
  strcpy(quest[1],"With normalisation");
  strcpy(quest[2],"\0");
  norm=GETV_QCM(
		TEXT0043,(char **)quest);

  Lecture_nimg_ds_nWins(n,norm);
}

/******************************************************************************
 ** -- n2Open_Norm() -----------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July 22nd, 1999
 **
 **  Opens multiple 2d pictures with offset
 **
 ******************************************************************************/

void n2Open_Cmd(void)
{
  /* Load several files with offset        */
  /*   1st, Call file selection box and:    */
  /*   2nd, Make severals files.         */
  char *file,**files;
  char name[FILE_LEN];
  int pos,i,max,step,offset,e;

  file=GET_FILE(NULL, &e);
  if (e) return;
  if(file==NULL || strlen(file)<=0)
    PUT_WARN("Sorry, cannot read picture !\n");
  else {
    pos=GET_WND(NULL);
    offset= GET_INT("Offset ?", 0, &e);
    if (e) return;
    max= GET_INT(TEXT0004, 0, &e);
    if (e) return;
    step= GET_INT(TEXT0024, 0, &e);
    if (e) return;
    if (offset!=0)
      {
	strcpy(name,file);
	incfic(name,offset);
	strcpy(file,name);
      }
    files=Make_files(file,step,max);
    /*protoize:???*/
    if(load_pictures(pos,MINI(max,MAX_PICTURES-pos+1),files/*,1*/)
       ==-1)
      PUT_ERROR("Cannot read picture.\n");
    for(i=0;i<max;i++)
      FREE(files[i]);
    FREE(files);
  }

}

/* --read_filetxt_3d() --------------------------------------------------
**
**    permet de lire un fichier contenant des coordonnees de point en 3D
**
**    Usage:    getheader_interfile(file,roi)
**        file: char*: file name 
**              roi : affichage du resultat
**
**	  On peut dans cette fonction commencer a la ligne : nu_ligne                     
*/
void
read_filetxt_3d(void)
{
  char *fmt1,*fmt2;
  char s[128];
  FILE *fp;
  int i,j,k,i1,j1,k1,rien;
  int nu_ligne,count,type;
  int radius,err;

  char *file;
  grphic3d *roi;
 
  roi=ptr_img_3d(0);
  file=GET_FILE("*.txt", &err);
  /*	strcpy(file,"/home/armspach/toto.txt"); */
	
  if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
    return;
  }
  fmt1=CALLOC(4096,char);
  fmt2=CALLOC(4096,char);

  count=0;
  nu_ligne=3;
  type=2;
  radius=3;

  while ((fgets(fmt1,4096,fp))!=NULL)
    {
      count++;
      if (count >nu_ligne-1 ) {  
        sscanf(fmt1,"%d %d %d %s %d %d %d %d",&i,&j,&k,s,&rien,&i1,&j1,&k1);
	printf("I=%d, j=%d, k=%d  s=%s I1=%d, j1=%d, k1=%d\n",i,j,k,s,i1,j1,k1);
        switch (type) {
	case 1:
	  roi->mri[i][j][k]=1;
	  break;

	case 2:
	  /*	create_sphere_3d(roi,radius,10,&nx,&ny);
		printf("nx=%d ny=%d\n",nx,ny); */ 
	  //creation_sphere_3d(roi,radius,i,j,k,1); 
	  UnAvailable();
	  /*	imx_midpoint_sphere_3d_p(roi,radius); */
	  break;
	} 
      }
    }

  FREE(fmt1);
  FREE(fmt2);
  fclose(fp);
  /*protoize:???*/
  return;
}

/******************************************************************************
 ** -- Open_PCX() --------------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July 22nd, 1999
 **
 **  Opens a PCX file.
 **
 ******************************************************************************/

void Open_PCX(void)
{/* Call file selection box and:            */
  /* Open PCX picture -  Lire une image PCX.    */
  char *ptr,*quest[16];
  int pos,reduc,e,i;

  ptr=GET_FILE(NULL, &e);
  if (e) return;
  if(ptr==NULL || strlen(ptr)<=0) {
    PUT_WARN("Cannot read picture\n");
  }
  else {
    pos=GET_PLACE("Image Depart de l'image :");
    for(i=0;i<16;i++)
      /*protoize:???*/
      quest[i]=CALLOC(80,char /*  * */);

    strcpy(quest[0],"No reduction");
    strcpy(quest[1],"Reduction 2");
    strcpy(quest[2],"Reduction 4");
    strcpy(quest[3],"Reduction 8");
    strcpy(quest[4],"\0");
    reduc=GETV_QCM(
		   TEXT0043,(char **)quest);
    switch (reduc)
      {
      case 0 : reduc = 1; break;
      case 1 : reduc = 2; break;
      case 2 : reduc = 4; break;
      case 3 : reduc = 8; break;
      default : //do nothin
	break;
      }
    loadpcx(ptr,pos,reduc);

    for(i=0;i<16;i++)
      FREE(quest[i]);
  }
} 

/******************************************************************************
 ** -- Open_GIF() --------------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July 22nd, 1999
 **
 **  Opens a GIF file.
 **
 ******************************************************************************/

void Open_GIF(void)
{/* Call file selection box and:            */
  /* Open GIF picture -  Lire une image GIF.    */
  char *ptr,*quest[16];
  int pos,reduc,e,i;

  ptr=GET_FILE(NULL, &e);
  if (e) return;
  if(ptr==NULL || strlen(ptr)<=0) {
    PUT_WARN("Cannot read picture\n");
  }
  else {
    pos=GET_PLACE("Image Depart de l'image :");
    for(i=0;i<16;i++)
      /*protoize:???*/
      quest[i]=CALLOC(80,char /*  * */);

    strcpy(quest[0],"No reduction");
    strcpy(quest[1],"Reduction 2");
    strcpy(quest[2],"Reduction 4");
    strcpy(quest[3],"Reduction 8");
    strcpy(quest[4],"\0");
    reduc=GETV_QCM(
		   TEXT0043,(char **)quest);
    switch (reduc)
      {
      case 0 : reduc = 1; break;
      case 1 : reduc = 2; break;
      case 2 : reduc = 4; break;
      case 3 : reduc = 8; break;
      default : //do nothin
	break;
      }
    loadgif(ptr,pos,reduc);

    for(i=0;i<16;i++)
      FREE(quest[i]);
  }
} 

/******************************************************************************
 ** -- Open_Raw() --------------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July 22nd, 1999
 **
 **  Opens a raw file.
 **
 ******************************************************************************/

void Open_Raw(void)
{ /* Lecture d'une image raw x.y. nb octets */
  int pos,taille,dimx,dimy,i,e;
  char *file,*quest[16];
  for(i=0;i<6;i++)
    /*protoize:???*/
    quest[i]=CALLOC(80,char /*  * */);
  file=GET_FILE(NULL, &e);
  if (e) return;

  strcpy(quest[0],"8  bits");
  strcpy(quest[1],"16 bits");
  strcpy(quest[2],"       ");
  strcpy(quest[3],"32 bits");
  strcpy(quest[4],"Cancel");
  strcpy(quest[5],"\0");
  taille=GETV_QCM(
		  TEXT0043,(char **)quest);
  for(i=0;i<6;i++)
    FREE(quest[i]);
  if (taille==4) return;

  dimx= GET_INT(TEXT0125, 0, &e);
  if (e) return;
  dimy= GET_INT(TEXT0126, 0, &e);
  if (e) return;
  if (dimx > MAX_WIDTH || dimy > MAX_HEIGHT ) {
    char s[256];
    sprintf(s,"Sorry, cannot load image size greater than %ldx%ld \n",MAX_WIDTH,MAX_HEIGHT);
    PUT_MESG(s);
    return;
  }
  pos=GET_WND(NULL);
  read_raw(file,dimx,dimy,pos,taille+1);
}

/******************************************************************************
 ** -- Open_nRaw() -------------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July 22nd, 1999
 **
 **  Opens a raw file.
 **
 ******************************************************************************/

void Open_nRaw(void)
{ /* Lecture d'une image raw x.y. nb octets */
  int pos,taille,dimx,dimy,max,step=0,i,e;
  char *file,*quest[16],**files;
  for(i=0;i<6;i++)
    /*protoize:???*/quest[i]=CALLOC(80,char /* * */);
  file=GET_FILE(NULL, &e);
  if (e) return;
  max= GET_INT(TEXT0004, 0, &e);
  if (!e)
    step= GET_INT(TEXT0024, 0, &e);
  if (e) return;
  files=Make_files(file,step,max);

  strcpy(quest[0],"8  bits");
  strcpy(quest[1],"16 bits");
  strcpy(quest[2],"        ");
  strcpy(quest[3],"32 bits");
  strcpy(quest[4],"Cancel");
  strcpy(quest[5],"\0");
  taille=GETV_QCM(
		  TEXT0043,(char **)quest);
  for(i=0;i<6;i++)
    FREE(quest[i]);
  if (taille==4) return;

  dimx= GET_INT(TEXT0125, 0, &e);
  if (e) return;
  dimy= GET_INT(TEXT0126, 0, &e);
  if (e) return;
  if (dimx > MAX_WIDTH || dimy > MAX_HEIGHT ) {
    char s[256];
    sprintf(s,"Sorry, cannot load image size greater than %ldx%ld \n",MAX_WIDTH,MAX_HEIGHT);
    PUT_MESG(s);
    return;
  }
  pos=GET_WND(NULL);
  for (i=0;i<max;i++)
    {
      if (i+pos > MAX_PICTURES) return ;
      read_raw(files[i],dimx,dimy,i+pos,taille+1);
    }
  for(i=0;i<max;i++)
    FREE(files[i]);
  FREE(files);
}

/******************************************************************************
 ** -- Save_Bru() --------------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July 22nd, 1999
 **
 ******************************************************************************/

void Save_Bru(void)
{
  char *str;
  int   e;
  
  str=(char *) malloc(PATH_LEN);
  sprintf(str,"File ? in %s",_CurrentPath);
  save_picture(GET_WND(NULL), 
	       SAVE_FILE(str,NULL, &e));
  FREE(str);
}

/******************************************************************************
 ** -- Save_Raw8() --------------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July 22nd, 1999
 **
 ******************************************************************************/

void Save_Raw8(void)
{
  char *str;
  int   e;

  str=(char *) malloc(PATH_LEN);
  sprintf(str,"File ? in %s",_CurrentPath);
  save_mri_raw8(GET_WND(NULL),
		SAVE_FILE(str,NULL, &e));
  FREE(str);
}

/******************************************************************************
 ** -- Save_Raw32() -------------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July 22nd, 1999
 **
 ******************************************************************************/

void Save_Raw32(void)
{
  char *str;
  int   e;

  str=(char *) malloc(PATH_LEN);
  sprintf(str,"File ? in %s",_CurrentPath);
  save_mri_raw32(GET_WND(NULL),
		 SAVE_FILE(str,NULL, &e));
  FREE(str);
}

/******************************************************************************
 ** -- Save_Pcx() --------------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July 22nd, 1999
 **
 ******************************************************************************/

void Save_Pcx(void)
{
  char *str;
  int   e;

  str=(char *) malloc(PATH_LEN);
  sprintf(str,"File ? in %s",_CurrentPath);
  save_mri_pcx(GET_WND(NULL),
	       SAVE_FILE(str,NULL, &e));
  FREE(str);
}

/******************************************************************************
 ** -- Save_2d_Ipb() -----------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July 22nd, 1999
 **
 ******************************************************************************/

void Save_2d_Ipb(void)
{
  char *file;
  int pos,e;
  char *str;

  str=(char *) malloc(PATH_LEN);
  sprintf(str,"File ? in %s",_CurrentPath);
  file=SAVE_FILE(str,NULL, &e);
  FREE(str);
  if(file == NULL || strlen(file) <= 0)
    PUT_WARN("Sorry, cannot save!\n");
  else {
    pos=GET_WND(NULL);
    if (save_pictures_2d(file,pos,1,1,IPB_2D_FORMAT) == 0)
      printf("Error save format IPB 2D\n");
  }
}

/******************************************************************************
 ** -- Save_2d_mask() -----------------------------------------------------------
 **
 ** >> Last user action by: T. BERST on July 8th, 2002
 **
 ******************************************************************************/

void Save_2d_mask(void)
{
  char *file;
  int pos,e;
  char *str;

  str=(char *) malloc(PATH_LEN);
  sprintf(str,"File ? in %s",_CurrentPath);
  file=SAVE_FILE(str,NULL, &e);
  FREE(str);
  if(file == NULL || strlen(file) <= 0)
    PUT_WARN("Sorry, cannot save!\n");
  else 
    {
      char mask_file[PATH_LEN];
      strcpy(mask_file,file);
      strcat(mask_file,"_mask");
      pos=GET_WND(NULL);

      if (save_pictures_2d(file,pos,1,1,IPB_2D_FORMAT) == 0)
	{	printf("Error save format IPB 2D\n");
	return;
	}
      if (save_pictures_2d(mask_file,-pos,1,1,IPB_2D_FORMAT) == 0)
	{	printf("Error save format IPB 2D\n");
	return;
	}
      //update the header file
      if (!add_mask_Ipb_2d(file,pos))
	{
	  PUT_ERROR(TEXT0422);
	  return ;
	}	
    }
}
/* --add_mask_Ipb_3d_3d() ---------------------------------------------------------
**
**    add_mask_Ipb_3d() adds mask item in the *.ipb.
**
**    Usage ...
**   file :       nom du fichier
**   pos :        position de l'image dont on veut sauver le masque 
**   
**	return 0 si echec 1 sinon 
*/

int add_mask_Ipb_2d(char * file, int pos)
{ 
  HEADER header;
  char s[512];
  char mask_file[512];
	
  strcpy(mask_file,file);
  strcat(file,".ipb");
  strcat(mask_file,"_mask.ipb");

  if(load_header(file, &header) == 0) {
    PUT_ERROR( "[save_mri_ipb_2d] Error reading header (3)!\n");
    return(0);
  }

  put_header(&header, "MASK_FILE=", STRING, mask_file, 0);

  if(save_header(&header, file) == 0) {
    sprintf(s,"Unable to save file header %s !\n", file);
    PUT_ERROR(s);
    return(0);
  }
  return(1);
}
/******************************************************************************
 ** -- Save_2d_Ipbr() ----------------------------------------------------------
 **
 ** >> Derniere revision par F. BOULEAU le 10 Aout 1999
 **
 **  Enregistrement et reduction d'une image 2d au format IPB
 **
 ******************************************************************************/

void Save_Ipbr(void)
{
  char *file;
  int   pos,e;
  char *str;
  grphic *tmp;

  str = (char *) malloc(PATH_LEN);
  sprintf(str, "File ? in %s", _CurrentPath);
  file = SAVE_FILE( str, NULL, &e);
  FREE(str);
	
  if(file == NULL || strlen(file) <= 0)
    PUT_WARN("Sorry, cannot save!\n");
  else 
    {
      pos = GET_WND( NULL);

      tmp = cr_grphic(ptr_img(pos));
      imx_copie_p(ptr_img(pos), tmp);
      scaleimg_p(tmp, ptr_img(pos));

      if (save_pictures_2d(file, pos, 1, 1, IPB_2D_FORMAT) == 0)
	printf( "Error save format IPB 2D\n");

      imx_copie_p(tmp, ptr_img(pos));
      free_grphic(tmp);
    }
}

/******************************************************************************
 ** -- Save_nBru() -------------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July 22nd, 1999
 **
 ******************************************************************************/

void Save_nBru(void)
{/* Save several files.             */
  /*   1st, Get primary file            */
  /*   2nd, Make severals files.         */
  char *file,**files;
  int pos,i,max,step,e;
  char *str;

  str=(char *) malloc(PATH_LEN);
  sprintf(str,"File ? in %s",_CurrentPath);
  file=SAVE_FILE(str,NULL, &e);
  FREE(str);
  if(file==NULL || strlen(file)<=0)
    PUT_WARN("Sorry, cannot save!\n");
  else {
    pos=GET_WND(NULL);
    max= GET_INT(TEXT0004, 0, &e);
    if (e) return;
    step= GET_INT(TEXT0024, 0, &e);
    if (e) return;
    files=Make_files(file,step,max);
    if(save_pictures(pos,max,files)==-1);
    PUT_ERROR("Cannot some pictures.\n");
    for(i=0;i<max;i++)
      FREE(files[i]);
    FREE(files);
  }

}

/******************************************************************************
 ** -- Save_nRaw8() ------------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July 22nd, 1999
 **
 ******************************************************************************/

void Save_nRaw8(void)
{/* Save several filesi raw data 8.             */
  /*   1st, Get primary file            */
  /*   2nd, Make severals files.         */
  char *file,**files;
  int pos,i,max,step,e;
  char *str;

  str=(char *) malloc(PATH_LEN);
  sprintf(str,"File ? in %s",_CurrentPath);
  file=SAVE_FILE(str,NULL, &e);
  FREE(str);
  if(file==NULL || strlen(file)<=0)
    PUT_WARN("Sorry, cannot save!\n");
  else {
    pos=GET_WND(NULL);
    max= GET_INT(TEXT0004, 0, &e);
    if (e) return;
    step= GET_INT(TEXT0024, 0, &e);
    if (e) return;
    files=Make_files(file,step,max);
    if(save_mri_nraw8(pos,max,files)==-1);
    PUT_ERROR("Cannot some pictures.\n");
    for(i=0;i<max;i++)
      FREE(files[i]);
    FREE(files);
  }

}

/******************************************************************************
 ** -- Save_nRaw32() -----------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July 22nd, 1999
 **
 ******************************************************************************/

void Save_nRaw32(void)
{/* Save several files in raw data 32    */
  /*   1st, Get primary file            */
  /*   2nd, Make severals files.         */
  char *file,**files;
  int pos,i,max,step,e;
  char *str;

  str=(char *) malloc(PATH_LEN);
  sprintf(str,"File ? in %s",_CurrentPath);
  file=SAVE_FILE(str,NULL, &e);
  FREE(str);
  if(file==NULL || strlen(file)<=0)
    PUT_WARN("Sorry, cannot save!\n");
  else {
    pos=GET_WND(NULL);
    max= GET_INT(TEXT0004, 0, &e);
    if (e) return;
    step= GET_INT(TEXT0024, 0, &e);
    if (e) return;
    files=Make_files(file,step,max);
    if(save_mri_nraw32(pos,max,files)==-1);
    PUT_ERROR("Cannot some pictures.\n");
    for(i=0;i<max;i++)
      FREE(files[i]);
    FREE(files);
  }

}

/******************************************************************************
 ** -- Save_nPcx() -------------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July 22nd, 1999
 **
 ******************************************************************************/

void Save_nPcx(void)
{/* Save several files in PCX format            */
  /*   1st, Get primary file            */
  /*   2nd, Make severals files.         */
  char *file,**files;
  int pos,i,max,step,e;
  char *str;

  str=(char *) malloc(PATH_LEN);
  sprintf(str,"File ? in %s",_CurrentPath);
  file=SAVE_FILE(str,NULL, &e);
  FREE(str);
  if(file==NULL || strlen(file)<=0)
    PUT_WARN("Sorry, cannot save!\n");
  else {
    pos=GET_WND(NULL);
    max= GET_INT(TEXT0004, 0, &e);
    if (e) return;
    step= GET_INT(TEXT0024, 0, &e);
    if (e) return;
    files=Make_files(file,step,max);
    if(save_mri_npcx(pos,max,files)==-1);
    PUT_ERROR("Cannot some pictures.\n");
    for(i=0;i<max;i++)
      FREE(files[i]);
    FREE(files);
  }

}

/******************************************************************************
 ** -- Save_nIpb2d() -----------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July 22nd, 1999
 **
 ******************************************************************************/

void Save_nIpb2d(void)
{
  char *file;
  int pos,max,step,e;
  char *str;

  str=(char *) malloc(PATH_LEN);
  sprintf(str,"File ? in %s",_CurrentPath);
  file=SAVE_FILE(str,NULL, &e);
  FREE(str);
  if(file == NULL || strlen(file) <= 0)
    PUT_WARN("Sorry, cannot save!\n");
  else {
    pos=GET_WND(NULL);
    max= GET_INT(TEXT0004, 0, &e);
    if (e) return;
    step= GET_INT(TEXT0024, 0, &e);
    if (e) return;
    if (save_pictures_2d(file,pos,max,step,IPB_2D_FORMAT) == 0)
      printf("error save format IPB 2D\n");
    return;
  }
}

/******************************************************************************
 ** -- Zoom_Visu() -------------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July 22nd, 1999
 **
 ******************************************************************************/

void Zoom_Visu(void)
{ /* Zoom for visualisation... */

  char *quest[16];
  int pos;
  long nzoom;

  pos=GET_WND(NULL);

  quest[0]="No reduction";
  quest[1]="Reduction half";
  quest[2]="Zoom 2";
  quest[3]="Zoom 4";
  quest[4]="\0";
  nzoom=GETV_QCM(
		 TEXT0043,(char **)quest);
  nzoom = (nzoom==0)? 1 : (nzoom==1) ? -2 : (nzoom==3)? 4:2;

  if ( nzoom==1) { /* Pas de deplacement de l'image */
    ptr_img(pos)->zoom_x=0;
    ptr_img(pos)->zoom_y=0;
    ptr_img(pos)->zoom=nzoom;
    show_picture(pos);
  }
  else /* Zoom avec deplacement   */ 
    {
      ptr_img(pos)->zoom=nzoom;
      show_picture(pos);
      IMXMouse_MoveImage(pos);
    }
}

/******************************************************************************
 ** -- Zoom_Visun() -------------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July 22nd, 1999
 **
 ******************************************************************************/

void Zoom_Visun(void)
{ /* Zoom for visualisation... */

  int pos,max,i,e;
  char *quest[16];
  long nzoom;
  pos=GET_WND(NULL);
  max= GET_INT(TEXT0004, 0, &e);

  if (!e)
    {
      quest[0]="No reduction";
      quest[1]="Reduction half";
      quest[2]="Zoom 2";
      quest[3]="Zoom 4";
      quest[4]="\0";
      nzoom=GETV_QCM(TEXT0043,(char **)quest);
      nzoom = (nzoom==0)? 1 : (nzoom==1) ? -2 : (nzoom==3)? 4:2;

      if(MESG_DEBUG) printf("  nzoom=%ld max=%lud %d \n",nzoom
			    ,(MAX_WIDTH/ptr_img(pos)->width)
			    ,-2>2 );

      for (i=0;i<max;i++)
	{
	  if(-2>2l)
	    {
	      if(MESG_DEBUG) printf("  nzoom=%ld\n",nzoom);
	      ptr_img(pos+i)->zoom=MAX_WIDTH/ptr_img(pos+i)->width;
	    }
	  else
	    ptr_img(pos+i)->zoom=nzoom;

	  show_picture(pos+i);
	}
    }
}

/******************************************************************************
 ** -- Zoom_Vi_2d_pos_xyz() -------------------------------------------------------------
 **
 ** Zoom une image 2d avec gestion d'un offset (zoom_x, zoom_y)
 **
 ******************************************************************************/

void Zoom_Vi_2d_pos_xy(int x, int y,int pos ,int nzoom)
{
  int wdth,hght;
  int zoom_x,zoom_y;
  int dflt_size=DEFAULT_IMAGIX_WIDTH;
	
  wdth = ptr_img(pos)->width;
  hght = ptr_img(pos)->height;
	
  ptr_img(pos)->zoom	=	nzoom;
	
  // calcul de la position du coin superieur gauche

  zoom_x=x-dflt_size/(nzoom *2);

  if ( (zoom_x+wdth/nzoom)>wdth )
    zoom_x = wdth-wdth/nzoom;
  ptr_img(pos)->zoom_x	=	MAXI(0,zoom_x); 
	
  zoom_y=y-dflt_size/(nzoom *2);

  if ( (zoom_y+hght/nzoom)>hght )
    zoom_y = hght-hght/nzoom;
  ptr_img(pos)->zoom_y	=	MAXI(0,zoom_y); 
	
  //---------------------------------//
  //-mise a jour du masque si existe //
  if (ptr_img(pos)->mask)
    {
      ptr_img(pos)->mask->zoom		= nzoom;
      ptr_img(pos)->mask->zoom_x	= MAXI(0,zoom_x);
      ptr_img(pos)->mask->zoom_y	= MAXI(0,zoom_y);
    }
	
  //----------------------------------------
  ptr_img(pos)->zoom_x *= nzoom;
  ptr_img(pos)->zoom_y *= nzoom;
  if (ptr_img(pos)->mask)
    {
      ptr_img(pos)->mask->zoom_x *= nzoom;
      ptr_img(pos)->mask->zoom_y *= nzoom;
    }
  /*
    printf("zoomx%d zoomy%d zoom%d\n",ptr_img(pos)->zoom_x,ptr_img(pos)->zoom_y,nzoom);
    printf("%d %d %d\n",x,y,z);
  */
  show_picture(pos);
}


/******************************************************************************
 ** -- scaleimg_p() ------------------------------------------------------------
 **
 ** >> Derniere revision par F. BOULEAU le  4 Aout 1999
 **
 **  Redimensionne une image 2d selon une methode parmi plusieurs.
 **
 **  Entree : imres : image resultat
 **
 **  Sortie : imres est modifie
 **
 ******************************************************************************/

int scaleimg_p(grphic *imdeb, grphic *imres)
{
  int		 meth,						/* methode de redimensionnement      */
    factx, facty,				/* facteurs de reduction             */
    rep;						/* reponse au qcm                    */
  int	 err;						/* erreur qcm                        */
  char	*qcm[] = {"1:2","1:4","Personnalise","" };
  char *qcm2[] = {"Saut de pixels", "Moyenne des pixels",
		  "Mediane des pixels", "MIP", ""};

  rep = GETV_QCM( "Reduction ?", (char **) qcm);

  switch(rep)
    {
    case 0:
      factx = facty = 2;
      break;

    case 1:		/* 1:2 ou 1:4 */
      factx = facty = 4;
      break;

    case 2:		/* Personnalise */
      factx =  GET_INT("X", 0, &err);

      if(!err)
	{
	  PUT_ERROR("Typing error!\n");
	  return ERR_TYPING;
	} 

      facty =  GET_INT("Y", 0, &err);

      if(!err)
	{
	  PUT_ERROR("Typing error!\n");
	  return ERR_TYPING;
	} 

      break;

    default:
      PUT_ERROR("Typing error!\n");
      return ERR_TYPING;
    }

 

  meth = GETV_QCM("Type de reduction ?", (char **) qcm2);


  return imx_scaleimg_p(imdeb, imres, meth, factx, facty);
}

/******************************************************************************
 ** -- imx_scaleimg_p() --------------------------------------------------------
 **
 ** >> Derniere revision par F. BOULEAU le  4 Aout 1999
 **
 **  Redimensionne une image 2d selon une methode parmi plusieurs.
 **
 **  Entree : imdeb : image de reference
 **           imres : image resultat
 **           meth  : methode a utiliser
 **           factx : facteur de reduction en x
 **           facty : facteur de reduction en y
 **
 **  Sortie : imres est modifie
 **
 ******************************************************************************/

int imx_scaleimg_p(grphic *imdeb, grphic *imres, int meth, int factx, int facty)
     /* image de reference                */
     /* image resultat                    */
     /* methode de redimensionnement      */
     /* facteurs de reduction             */
{
  int  i, j, k;						/* indices de parcours de l'image    */
  int  v;								/* valeur du pixel calcule           */
  int  resi, resj;					/* indices des pixels de l'image res */
  int *tri;							/* valeurs des pixels triees         */
  int  n;								/* parcours du tableau "tri"         */
  int  tmp;							/* stockage temporaire des pixels    */
  int  mini=0, maxi=0, once = 0;			/* cutoff_min, cutoff_max            */

  /* -- calcul des dimensions de la nouvelle image ----------------------- */

  imx_copie_param_p(imdeb, imres);

  imres->width = imdeb->width / factx 
    + (imdeb->width % factx > 0 ? 1 : 0);
  imres->height = imdeb->height / facty 
    + (imdeb->height % facty > 0 ? 1 : 0);

  imres->dx = imres->dx * factx;
  imres->dy = imres->dy * facty;

  /* -- calcul des pixels ------------------------------------------------ */

  tri = CALLOC(factx * facty, int);

  for(i = resi = 0 ; i < imdeb->width ; i += factx, resi++)
    for(j = resj = 0 ; j < imdeb->height ; j += facty, resj++)
      {
	/* traitement des pixels de la "fenetre" de parcours de l'image  */
	/* ("fenetre" : rectangle de factx sur facty pixels)             */

	switch(meth)
	  {
	  case REDUC_JUMP:	/* un pixel sur n */
	    imres->mri[resi][resj] = imdeb->mri[i][j];
	    break;
			
	  case REDUC_MEAN:	/* moyenne des pixels */
	    for(v = k = 0 ; k < factx * facty ; k++)
	      v += imdeb->mri[i + k / factx][j + k % facty];
	    imres->mri[resi][resj] = v / (factx * facty);
	    break;
			
	  case REDUC_MEDIAN:	/* valeur mediane des pixels */
	    for(k = 0 ; k < factx * facty ; k++)
	      {
		v = imdeb->mri[i + k / factx][j + k % facty];

		/* insertion de la valeur dans le tableau "tri" */

		for(n = 0 ; n < k ; n++)
		  {
		    if (tri[n] > v)
		      {
			tmp = v;
			v = tri[n];
			tri[n] = tmp;
		      }
		  }

		tri[n] = v;
	      }

	    /* recuperation de la valeur mediane */
				
	    imres->mri[resi][resj] = tri[factx * facty / 2];
	    break;
				
	  case REDUC_MIP:		/* valeur maximale des pixels */
	    for(v = k = 0 ; k < factx * facty ; k++)
	      v = MAXI(v, imdeb->mri[i + k / factx][j + k % facty]);
	    imres->mri[resi][resj] = v;
	    break;

	    if (once == 0)
	      {
		mini = maxi = imres->mri[resi][resj];
		once = 1;
	      }
	    else
	      {
		mini = MINI(mini, imres->mri[resi][resj]);
		maxi = MAXI(maxi, imres->mri[resi][resj]);
	      }
	  }
      }			/* fin boucle i et j */

  imres->cutoff_min = mini;
  imres->cutoff_max = maxi;

  FREE(tri);

  return 0;
}


/* --read_filetxt() --------------------------------------------------
**
**    permet de lire un fichier contenant des coordonnees de point en 3D
**
**    Usage:    getheader_interfile(file,roi)
**        file: char*: file name 
**              roi : affichage du resultat
**
**	  On peut dans cette fonction commencer a la ligne : nu_ligne                     
*/
int
read_filetxt(char *file, grphic3d *roi)
{
  char *fmt1,*fmt2;
  FILE *fp;
  int i,j,k,i1,j1,k1,rien;
  int nu_ligne,count;
  int err;
 
  file=GET_FILE("*.txt", &err);
  /*	strcpy(file,"/home/armspach/toto.txt"); */
	
  if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
    return(0l); 
  }
  fmt1=CALLOC(4096,char);
  fmt2=CALLOC(4096,char);

  count=0;
  nu_ligne=3;

  while ((fgets(fmt1,4096,fp))!=NULL)
    {
      count++;
      if (count >nu_ligne-1 ) {  
	sscanf(fmt1,"%d %d %d %d %d %d %d %d",&i,&j,&k,&rien,&rien,&i1,&j1,&k1);
	printf("I=%d, j=%d, k=%d  I1=%d, j1=%d, k1=%d\n",i,j,k,i1,j1,k1);
	/*      if (strstr(fmt1,string))
		{
		sprintf(fmt2,"%s %%d",string);
		sscanf(fmt1,fmt2,&i);
		sprintf(s,"%d",i);
		}
	*/	
      }
    }

  FREE(fmt1);
  FREE(fmt2);
  fclose(fp);
  /*protoize:???*/
  return(0/*???*/);
}

/* --read_filetxt_sep_3d() --------------------------------------------------
**
**    permet de lire un fichier contenant des coordonnees de point en 3D
**    puis interprette le 4 champ pour decider de la couleur
**         + : rouge
**         - : vert
**         0 : rien
**    Usage:    read_filetxt_sep_3d()
**
**	  On peut dans cette fonction commencer a la ligne : nu_ligne                     
*/
void
read_filetxt_sep_3d(void)
{
  char *fmt1,*fmt2;
  char s[128];
  FILE *fp;
  int i,j,k,i1,j1,k1,rien;
  int nu_ligne,count,type;
  int radius,err;

  char *file;
  grphic3d *roi;
 
  roi=ptr_img_3d(0);
  file=GET_FILE("*.txt", &err);
  /*	strcpy(file,"/home/armspach/toto.txt"); */
	
  if((fp=fopen(file,"rb"))==NULL) { /* Open file error... */
    /*protoize:???*/
    return;
    /*          return(0l); */
  }
  fmt1=CALLOC(4096,char);
  fmt2=CALLOC(4096,char);

  count=0;
  nu_ligne=3;
  type=2;
  radius=3;

  while ((fgets(fmt1,4096,fp))!=NULL)
    {
      count++;
      if (count >nu_ligne-1 ) {  
	sscanf(fmt1,"%d %d %d %s %d %d %d %d",&i,&j,&k,s,&rien,&i1,&j1,&k1);
	printf("I=%d, j=%d, k=%d  s=%s I1=%d, j1=%d, k1=%d\n",i,j,k,s,i1,j1,k1);
	switch (type) {
	case 1:
	  if (strstr(s,"+")) 
	    roi->mri[i][j][k]=1;
	  if (strstr(s,"-")) 
	    roi->mri[i][j][k]=-1;
	  break;

	case 2:
	  /*	create_sphere_3d(roi,radius,10,&nx,&ny);
		printf("nx=%d ny=%d\n",nx,ny); */ 
	  if (strstr(s,"+")) 
	    //creation_sphere_3d(roi,radius,i,j,k,1); 
	    UnAvailable();
	  if (strstr(s,"-")) 
	    //creation_sphere_3d(roi,radius,i,j,k,-1); 
	    UnAvailable();
	  /*	imx_midpoint_sphere_3d_p(roi,radius); */
	  break;
	} 
      }
    }

  FREE(fmt1);
  FREE(fmt2);
  fclose(fp);
  /*      return; */
}
 
 
/* ------- imx_saveonfloppy() ------------------------------------------
***
***
***
***
***
***
***
*** ---------------------------------------------------------------------*/

int	imx_saveonfloppy(void)
{
  int pos,max,step,i,n=0,e,numfirst;
  char *file;
  char filename[FILE_LEN],filefloppy[FILE_LEN],namefile[FILE_LEN];
  char nom[5];
  grphic *image;


  PUT_WNDMSG(TEXT0135);
  file=GET_FILE_IMG(NULL, &numfirst, &e);
  CLEAR_WNDMSG();
  if (e) {	
    PUT_WARN("Sorry, cannot read picture !\n");
    return(-1);
  }
  for (i=(strlen(file));file[i]!='/';i--)
    n=i;
  for (i=n;file[i]!='.'&&i < n+6 ;i++) {
    namefile[i-n]=file[i];
    //namefile[i-n+1]=NULL;   <- n'importe quoi... RC le 25.10.2001
    namefile[i-n+1]=0;
  }
  strcpy(filename,file);

  /*pos=GET_WND(NULL);
    max= GET_INT(TEXT0004, 0, &e);
    if (e) return(-1);
    step= GET_INT(TEXT0024, 0, &e);
    if (e) return(-1);*/

  max= GET_INT(TEXT0004, 0, &e);
  if (e) return(-1);
  pos=1;
  step=1;

  strcpy(filefloppy,_home);
  strcat(filefloppy,"/floppy/");
  strcat(filefloppy,namefile);
  strcat(filefloppy,"01.PCX");

  /* paravision*/

  if (is_mri(file,NULL)==PARAVISION_FORMAT) 	
    {
      load_picture(pos,filename,numfirst);
      image=ptr_img(pos);
      for(i=0;i<4;i++) nom[i]=image->patient.name[i];
      strcpy(filefloppy,_home);
      strcat(filefloppy,"/floppy/");
      strcat(filefloppy,nom);
      strcat(filefloppy,_path_paravision);
      strcat(filefloppy,"01.PCX");
    }
  /* Deleted Wiest 2004 
     if(_cardio_save_zoom==TRUE)	init_zoom2_movie(pos,file,numfirst);
  */
  _movie_cardiowithroi=FALSE;

  // Deleted Wiest 2004 cf deleted/cardio.c
  // IMXFile_Animate(pos,max,step,file,1,NULL);

  XGET_OKCANCEL(&e);
  if(e == 1) // OK
    {
      if(_cardio_save_zoom==TRUE)	unzoom_img(pos);	
      return(-1);
    }
	
  for (i=0;i<max;i++)
    {
      load_picture(pos,filename,numfirst);
      save_mri_pcx(pos,filefloppy);
      incfic_img(filename,step,&numfirst);
      incfic(filefloppy,step);
    }
  if(_cardio_save_zoom==TRUE)	unzoom_img(pos);	
  return(1);

}

/*--------  compressave_serial_sa()   ----------------------------------
***
***	Compression et sauvegarde des images
***		
***	Compression d'une s�ie d'images et sauvegarde sur disquette:
***	- nom de la s�ie d'images
***	- sans v�ification par animation
***	- compression de la s�ie d'images par imx_compres_serial dans imx_zone.c
***	- sauvegarde avec l'extension .C001 etc sur disquette
***	et dans le r�ertoire temp de coeur.
***
------------------------------------------------------------------*/
int	compressave_serial_sa(void)
{
  int pos1=0,pos2,number,step,i,j,n=0,e,numfirst;
  char *file,img_name[FILE_LEN];
  char filename[FILE_LEN],filerest[FILE_LEN],cmd[FILE_LEN];
  char nom[5];
  grphic *mask,*image;

  PUT_WNDMSG(TEXT0135);
  file=GET_FILE_IMG(NULL, &numfirst, &e);
  CLEAR_WNDMSG();
  if (e) return(-1);
  if(file==NULL || strlen(file)<=0)
    {
      PUT_WARN("Sorry, cannot read picture !\n");
      return(-1);
    }

  if (is_mri(file,NULL) !=PARAVISION_FORMAT)
    {
      /*   recherche image num�o 1 de la s�ie   */

      for(i=0;file[i]!='.';i++)
	{
	  filename[i]=file[i];
	  //filename[i+1]=NULL;  <- faux... RC.
	  filename[i+1]=0;
	}
      strcpy(file,filename);
      strcat(file,".I001");

      /*   extraction du nom de l'image    */

      for (i=(strlen(file));file[i]!='/';i--)
	n=i;
      for (i=n;file[i]!='.'&&i < n+6 ;i++)
	{
	  img_name[i-n]=file[i];
	  //img_name[i-n+1]=NULL; <- faux. RC.
	  img_name[i-n+1]=0;
	}

      /*   fichiers des images r�ultat   */

      strcpy(filerest,_home);
      strcat(filerest,"/");
      strcat(filerest,img_name);
      strcat(filerest,"01.PCX");
    }
  else
    {
      image=ptr_img(pos1);
      load_mri(file,image,numfirst);
      for (i=0;i<4;i++) nom[i]=image->patient.name[i];
      nom[4]='\0';
      strcpy(filerest,_home);
      strcat(filerest,"/");
      strcat(filerest,nom);
      strcat(filerest,_path_paravision);
      strcat(filerest,"01.PCX");
    } 


  /*   choix du nombre d'images, du pas et des positions  */

  /*	PUT_WNDMSG(TEXT0066);
	pos1=GET_WND(NULL);*/
  pos1=1;
  /*	PUT_WNDMSG(TEXT0006);
	pos2=GET_WND(NULL);*/
  pos2=2;
  /*	CLEAR_WNDMSG();*/
  number= GET_INT(TEXT0112, 0, &e);
  if (e) return(-1);
  /*	step= GET_INT(TEXT0024, 0, &e);
	if (e) return(-1);*/
  step=1;
	
  /*   suppression du bruit   */

  /*   Initialisation du masque */

  mask=ptr_img(pos2);
  imx_iniparaimg(pos2);
  mask->max_pixel=1;
  mask->min_pixel=0;
  mask->cutoff_max=1;
  mask->cutoff_min=0;
  for(j=0;j<_imx_img->height;j++)
    for(i=0;i<_imx_img->width;i++)
      mask->mri[i][j]=0;
		
  /*   visualisation image de depart  */ 

  image=ptr_img(pos1);
  if(load_mri(file,image,numfirst)==NULL)
    return(-1);
  show_picture(pos1);


  /*   cr�tion du masque   */

  if(imx_nonoise_mask(mask,image,pos2)!=1)	return(-1);

  /*   visualisation du masque pour confirmation   */

  show_picture(pos2);


  for (i=0;i<number;i++)
    {
      load_picture(pos1,file,numfirst);
      imx_and(pos1,pos2,pos1);
      save_mri_pcx(pos1,filerest);
      incfic_img(file,step,&numfirst);
      incfic(filerest,step);
    }

  XGET_OKCANCEL(&e);
  if(e==1)
    {
      strcpy(cmd,"rm ");
      strcat(cmd,_home);
      strcat(cmd,"/");
      if (is_mri(file,NULL)!=PARAVISION_FORMAT)
	strcat(cmd,img_name); 
      else 
	{  
	  strcat(cmd,nom);
	  strcat(cmd,_path_paravision);
	} 
      strcat(cmd,"*.PCX ");
      system(cmd);

      return(-1);
    }
  /*   sauvegarde sur disquette  */
  strcpy(cmd,"mv ");
  strcat(cmd,_home);
  strcat(cmd,"/");
  if (is_mri(file,NULL)!=PARAVISION_FORMAT)
    strcat(cmd,img_name); 
  else 
    { 
      strcat(cmd,nom);
      strcat(cmd,_path_paravision);
    } 
  strcat(cmd,"*.PCX ");
  strcat(cmd,_home);
  strcat(cmd,"/floppy/");
  system(cmd);
  return(1);

}

/* ------- imx_saveonfloppy_sa() ------------------------------------------
***
***
***
***		sans animation
***
***
***
*** ---------------------------------------------------------------------*/

int	imx_saveonfloppy_sa(void)
{
  int pos,max,step,i,n=0,e,numfirst;
  char *file;
  char filename[FILE_LEN],filefloppy[FILE_LEN],namefile[FILE_LEN];
  char nom[5];
  char separe;
  grphic *image;

#if defined (WIN32)
  separe='\\';
#else
  separe='/';
#endif		

  PUT_WNDMSG(TEXT0135);
  file=GET_FILE_IMG(NULL, &numfirst, &e);
  CLEAR_WNDMSG();
  if (e)
    {	
      PUT_WARN("Sorry, cannot read picture !\n");
      return(-1);
    }
  for (i=(strlen(file));file[i]!=separe;i--)
    n=i;
  for (i=n;file[i]!='.'&&i < n+6 ;i++)
    {
      namefile[i-n]=file[i];
      //namefile[i-n+1]=NULL;
      namefile[i-n+1]=0; // RC, le 25.10.2001
    }
  strcpy(filename,file);
  /*pos=GET_WND(NULL);
    max= GET_INT(TEXT0004, 0, &e);
    if (e) return(-1);
    step= GET_INT(TEXT0024, 0, &e);
    if (e) return(-1);*/
  pos=1;
  step=1;
  max= GET_INT(TEXT0004, 0, &e);
  if (e) return(-1);

  strcpy(filefloppy,_home);
  strcat(filefloppy,"/floppy/");
  strcat(filefloppy,namefile);
  strcat(filefloppy,"01.PCX");
	
  /* paravision : modification du nom du fichier de sauvegarde */

  if (is_mri(file,NULL)==PARAVISION_FORMAT) 	
    {
      load_picture(pos,filename,numfirst);
      image=ptr_img(pos);
      for(i=0;i<4;i++) nom[i]=image->patient.name[i];
      strcpy(filefloppy,_home);
      strcat(filefloppy,"/floppy/");
      strcat(filefloppy,nom);
      strcat(filefloppy,_path_paravision);
      strcat(filefloppy,"01.PCX");
    }

  for (i=0;i<max;i++)
    {
      load_picture(pos,filename,numfirst);
      save_mri_pcx(pos,filefloppy);
      incfic_img(filename,1,&numfirst);
      incfic(filefloppy,step);
    }

  return(1);

}

/*******************************************
 ** --  loadpcx() ----------------------------
 **
 ** Lecture fichier PCX
 ********************************************/

void 
loadpcx(char *file, int pos, int reduc)
{
  register int        i;
  FILE                *ifp;
  int                 Version;
  int                 Xmin, Ymin, Xmax, Ymax;
  int                 Width, Height;
  register int        x, y;
  int                 Planes;
  int                 BitsPerPixel;
  int                 BytesPerLine;
  unsigned char       Red[MAXCOLORS], Green[MAXCOLORS], Blue[MAXCOLORS];
  unsigned char       *pcximage;
  unsigned char       *pcxplanes;
  unsigned char       *pcxpixels;
  int im_res[2*MAX_PICTURES];
  int j,k,l;
  int pasx,pasy,pasposx,pasposy;
  grphic *imres[2*MAX_PICTURES];
  char tab[80];
  char mess[80];
  pixel               **pixels;

  /*-----------------------------------------*/
  /*-----------------------------------------*/


  ifp = p_openr(file);
  /*
   * read the PCX header
   */
  if (GetByte(ifp) != PCX_MAGIC)
    {
      sprintf(mess,"%s is not a PCX file \n", file );
      PUT_MESG(mess);
      return;
    }
  else if(MESG_DEBUG) printf("%s est un fichier PCX\n", file );
  Version     = GetByte(ifp);  /* get version #                       */

  if (GetByte(ifp) != 1)       /* check for PCX run length encoding   */
    {
      sprintf(mess,"%s has unknown encoding scheme \n", file );
      PUT_MESG(mess);
      return;
    }
  BitsPerPixel= GetByte(ifp);
  Xmin        = GetWord(ifp);
  Ymin        = GetWord(ifp);
  Xmax        = GetWord(ifp);
  Ymax        = GetWord(ifp);

  Width       = (Xmax - Xmin) + 1;
  Height      = (Ymax - Ymin) + 1;


  /*-----------------------------------------*/
  sprintf(tab,">%s< PCX\ninfo : Taille de l'image PCX %dx%d \n",file,Width,Height);
  PUT_MESG(tab);

  im_res[0]=pos;
  pasx= GINT(GINT(Width,MAX_WIDTH),reduc);
  pasy= GINT(GINT(Height,MAX_HEIGHT),reduc);
  if (pasx*pasy> 2*MAX_PICTURES)
    {
      PUT_MESG("Erreur : Taille de l'image trop grande\n");
      return;
    }

  for (j=0;j<pasy;j++)
    {
      l=pos+j*MAX_PIC_IN_LINE;
      for (i=0;i<pasx;i++)
	im_res[i+j*pasx]=(l>((GINT(pos,MAX_PIC_IN_LINE)+j)*MAX_PIC_IN_LINE)||l>MAX_PICTURES) ? 0 :l++;
    }

  for (j=0;j<pasy;j++)
    for (i=0;i<pasx;i++)
      if(im_res[i+j*pasx]!=0)
	imres[i+j*pasx]=ptr_img(im_res[i+j*pasx]);
  /*-----------------------------------------*/


  (void) GetWord(ifp);                /* ignore horizontal resolution */
  (void) GetWord(ifp);                /* ignore vertical resolution   */

  /*
   * get the 16-color color map
   */
  for (i = 0; i < 16; i++)
    {
      Red[i]    = GetByte(ifp);
      Green[i]  = GetByte(ifp);
      Blue[i]   = GetByte(ifp);
    }

  (void) GetByte(ifp);                /* skip reserved byte    */
  Planes      = GetByte(ifp);         /* # of color planes     */
  BytesPerLine= GetWord(ifp);         /* # of bytes per line   */
  (void) GetWord(ifp);                /* ignore palette info   */

  /*
   * check that we can handle this image format
   */
  switch (BitsPerPixel)
    {
    case 1:
      if (Planes > 4)
	{
	  PUT_MESG("can't handle image with more than 4 planes\n");
	  return;
	}
      break;

    case 2:
    case 4:
    case 8:
      if (Planes == 1)
	break;
    default:
      {
	sprintf(mess,"can't handle %d bits per pixel image with %d planes\n",
		BitsPerPixel,Planes);
	PUT_MESG(mess);
	return;
      }
    }


  /*
   * read the pcx format image
   */
  fseek(ifp, (long)PCX_HDR_SIZE, 0);
  pcximage    = (unsigned char *)allocrow(BytesPerLine * Planes, Height);
  read_pcx_image(ifp, pcximage, BytesPerLine, Planes, Height);

  /*
   * 256 color images have their color map at the end of the file
   * preceeded by a magic byte
   */
  if (BitsPerPixel == 8)
    {
      if (GetByte(ifp) != PCX_256_COLORS)
	{
	  PUT_MESG("bad color map signature\n");
	  return;
	}
      for (i = 0; i < MAXCOLORS; i++)
	{
	  Red[i]      = GetByte(ifp);
	  Green[i]    = GetByte(ifp);
	  Blue[i]     = GetByte(ifp);
	}
    }

  pixels      = pcx_allocarray(Width, Height);
  pcxpixels   = (unsigned char *)allocrow(Width+7, 1);

  /*
   * convert the image
   */
  for (y = 0; y < Height; y++)
    {
      pcxplanes = pcximage + (y * BytesPerLine * Planes);

      if (Planes == 1)
	{
	  pcx_unpack_pixels(pcxpixels, pcxplanes,
			    BytesPerLine, Planes, BitsPerPixel);
	}
      else
	{
	  pcx_planes_to_pixels(pcxpixels, pcxplanes,
			       BytesPerLine, Planes, BitsPerPixel);
	}

      for (x = 0; x < Width; x++)
	{
	  i = pcxpixels[x];
	  P_ASSIGN(pixels[y][x], Red[i], Green[i], Blue[i]);
	  /*==========================>*/



	  pasposx= GINT(GINT(x+1,MAX_WIDTH),reduc);
	  pasposy= GINT(GINT(y+1,MAX_HEIGHT),reduc);
	  pasposx--;
	  pasposy--;
	  k=pasposx+(pasx*pasposy);
	  i=(x-pasposx*MAX_WIDTH*reduc)/reduc;
	  j=(y-pasposy*MAX_HEIGHT*reduc)/reduc;
	  if(im_res[k]!=0)
	    imres[k]->mri[i][j]=
	      (gray)(P_LUMIN(pixels[y][x])+0.5);

	  /*==========================>*/


	}
    }

  /*==========================>*/
  for (y=0; y<(reduc*pasy*MAX_HEIGHT);y++)
    for (x=0; x< (reduc*pasx*MAX_WIDTH); x++)
      {
	if(x>= Width || y >= Height)
	  {
	    pasposx= GINT(GINT(x+1, MAX_WIDTH),reduc);
	    pasposy= GINT(GINT(y+1, MAX_HEIGHT),reduc);
	    pasposx--;
	    pasposy--;
	    k=pasposx+(pasx*pasposy);
	    i=(x-pasposx*MAX_WIDTH*reduc)/reduc;
	    j=(y-pasposy*MAX_HEIGHT*reduc)/reduc;
	    if(im_res[k]!=0)
	      imres[k]->mri[i][j]=0;
	  }
      }

  for (j=0;j<pasy;j++)
    for (i=0;i<pasx;i++)
      {
	imx_iniparam(im_res[i+j*pasx],255,0,0.,1.);
	show_picture(im_res[i+j*pasx]);
      }

  /*==========================>*/


  p_close(ifp);

  FREE(pixels);
  FREE(pcxpixels);
  FREE(pcximage);


}

/*******************************************
 ** --  loadpcx_p() ----------------------------
 **
 ** Lecture fichier PCX
 ********************************************/

void 
loadpcx_p(char *file, grphic *graphic, int reduc)
{
  register int        i;
  FILE                *ifp;
  int                 Version;
  int                 Xmin, Ymin, Xmax, Ymax;
  int                 Width, Height;
  register int        x, y;
  int                 Planes;
  int                 BitsPerPixel;
  int                 BytesPerLine;
  unsigned char       Red[MAXCOLORS], Green[MAXCOLORS], Blue[MAXCOLORS];
  unsigned char       *pcximage;
  unsigned char       *pcxplanes;
  unsigned char       *pcxpixels;
  pixel               **pixels;

  /*-----------------------------------------*/
  int im_res[2*MAX_PICTURES];
  int j,k,l;
  int pasx,pasy,pasposx,pasposy;
  grphic *imres[2*MAX_PICTURES];
  char tab[80];
  char mess[80];
  /*-----------------------------------------*/

  int pos=1;


  ifp = p_openr(file);
  /*
   * read the PCX header
   */
  if (GetByte(ifp) != PCX_MAGIC)
    {
      sprintf(mess,"%s is not a PCX file \n", file );
      PUT_MESG(mess);
      return;
    }
  else if(MESG_DEBUG) printf("%s est un fichier PCX\n", file );
  Version     = GetByte(ifp);  /* get version #                       */

  if (GetByte(ifp) != 1)       /* check for PCX run length encoding   */
    {
      sprintf(mess,"%s has unknown encoding scheme \n", file );
      PUT_MESG(mess);
      return;
    }
  BitsPerPixel= GetByte(ifp);
  Xmin        = GetWord(ifp);
  Ymin        = GetWord(ifp);
  Xmax        = GetWord(ifp);
  Ymax        = GetWord(ifp);

  Width       = (Xmax - Xmin) + 1;
  Height      = (Ymax - Ymin) + 1;


  /*-----------------------------------------*/
  sprintf(tab,">%s< PCX\ninfo : Taille de l'image PCX %dx%d \n",file,Width,Height);
  PUT_MESG(tab);

  im_res[0]=pos;
  pasx= GINT(GINT(Width,MAX_WIDTH),reduc);
  pasy= GINT(GINT(Height,MAX_HEIGHT),reduc);
  if (pasx*pasy> 2*MAX_PICTURES)
    {
      PUT_MESG("Erreur : Taille de l'image trop grande\n");
      return;
    }

  for (j=0;j<pasy;j++)
    {
      l=pos+j*MAX_PIC_IN_LINE;
      for (i=0;i<pasx;i++)
	im_res[i+j*pasx]=(l>((GINT(pos,MAX_PIC_IN_LINE)+j)*MAX_PIC_IN_LINE)||l>MAX_PICTURES) ? 0 :l++;
    }

  if (pasx!=1 || pasy!=1) {
    return;
  }

  imres[0]=graphic;
  /*   for (j=0;j<pasy;j++)
       for (i=0;i<pasx;i++)
       if(im_res[i+j*pasx]!=0)
       imres[i+j*pasx]=ptr_img(im_res[i+j*pasx]);
  */
  /*-----------------------------------------*/


  (void) GetWord(ifp);                /* ignore horizontal resolution */
  (void) GetWord(ifp);                /* ignore vertical resolution   */

  /*
   * get the 16-color color map
   */
  for (i = 0; i < 16; i++)
    {
      Red[i]    = GetByte(ifp);
      Green[i]  = GetByte(ifp);
      Blue[i]   = GetByte(ifp);
    }

  (void) GetByte(ifp);                /* skip reserved byte    */
  Planes      = GetByte(ifp);         /* # of color planes     */
  BytesPerLine= GetWord(ifp);         /* # of bytes per line   */
  (void) GetWord(ifp);                /* ignore palette info   */

  /*
   * check that we can handle this image format
   */
  switch (BitsPerPixel)
    {
    case 1:
      if (Planes > 4)
	{
	  PUT_MESG("can't handle image with more than 4 planes\n");
	  return;
	}
      break;

    case 2:
    case 4:
    case 8:
      if (Planes == 1)
	break;
    default:
      {
	sprintf(mess,"can't handle %d bits per pixel image with %d planes\n",
		BitsPerPixel,Planes);
	PUT_MESG(mess);
	return;
      }
    }


  /*
   * read the pcx format image
   */
  fseek(ifp, (long)PCX_HDR_SIZE, 0);
  pcximage    = (unsigned char *)allocrow(BytesPerLine * Planes, Height);
  read_pcx_image(ifp, pcximage, BytesPerLine, Planes, Height);

  /*
   * 256 color images have their color map at the end of the file
   * preceeded by a magic byte
   */
  if (BitsPerPixel == 8)
    {
      if (GetByte(ifp) != PCX_256_COLORS)
	{
	  PUT_MESG("bad color map signature\n");
	  return;
	}
      for (i = 0; i < MAXCOLORS; i++)
	{
	  Red[i]      = GetByte(ifp);
	  Green[i]    = GetByte(ifp);
	  Blue[i]     = GetByte(ifp);
	}
    }

  pixels      = pcx_allocarray(Width, Height);
  pcxpixels   = (unsigned char *)allocrow(Width+7, 1);

  /*
   * convert the image
   */
  for (y = 0; y < Height; y++)
    {
      pcxplanes = pcximage + (y * BytesPerLine * Planes);

      if (Planes == 1)
	{
	  pcx_unpack_pixels(pcxpixels, pcxplanes,
			    BytesPerLine, Planes, BitsPerPixel);
	}
      else
	{
	  pcx_planes_to_pixels(pcxpixels, pcxplanes,
			       BytesPerLine, Planes, BitsPerPixel);
	}

      for (x = 0; x < Width; x++)
	{
	  i = pcxpixels[x];
	  P_ASSIGN(pixels[y][x], Red[i], Green[i], Blue[i]);
	  /*==========================>*/



	  pasposx= GINT(GINT(x+1,MAX_WIDTH),reduc);
	  pasposy= GINT(GINT(y+1,MAX_HEIGHT),reduc);
	  pasposx--;
	  pasposy--;
	  k=pasposx+(pasx*pasposy);
	  i=(x-pasposx*MAX_WIDTH*reduc)/reduc;
	  j=(y-pasposy*MAX_HEIGHT*reduc)/reduc;
	  if(im_res[k]!=0)
	    imres[k]->mri[i][j]=
	      (gray)(P_LUMIN(pixels[y][x])+0.5);

	  /*==========================>*/


	}
    }

  /*==========================>*/
  for (y=0; y<(reduc*pasy*MAX_HEIGHT);y++)
    for (x=0; x< (reduc*pasx*MAX_WIDTH); x++)
      {
	if(x>= Width || y >= Height)
	  {
	    pasposx= GINT(GINT(x+1,imres[0]->width),reduc);
	    pasposy= GINT(GINT(y+1,imres[0]->height),reduc);
	    pasposx--;
	    pasposy--;
	    k=pasposx+(pasx*pasposy);
	    i=(x-pasposx*MAX_WIDTH*reduc)/reduc;
	    j=(y-pasposy*MAX_HEIGHT*reduc)/reduc;
	    if(im_res[k]!=0)
	      imres[k]->mri[i][j]=0;
	  }
      }


  /*==========================>*/


  p_close(ifp);

  FREE(pixels);
  FREE(pcxpixels);
  FREE(pcximage);


}

/* --save_mri_pcx() --------------------------------------    
**
**   Permet de sauver une image dans un fichier PCX
**   en 8 bits
**
******************************************************/
int save_mri_pcx (int pos, char *file)
{
  int i,j,err;
  byte line[MAX_WIDTH];
  grphic *img_prt;
  FILE *hand;
  IMXColor colors[256];

  double I,Red,Green,Blue;
  int    ncolors=256;


  /* Couleurs     */

  for(i=0;i<ncolors;i++) {
    I=(double)i;
    Red=255*(1-exp(-I/100))/(1-exp(-1.90));
    if(Red>255) Red  =255;
    if(i  <100) Green=0;  
    else Green=255*I/(255-100)-100*255/(255-100);
    if(i  <190) Blue =0;  
    else Blue =198*I/(255-190)-190*198/(255-190);

    colors[i].red  =(int)Red  ;
    colors[i].green=(int)Green;
    colors[i].blue =(int)Blue ;
  }


  /*  Controle du fichier   */

  if((hand=fopen(file,"rb"))==NULL)
    {
      if(!errno==2)              /* else file does not exist ! */
	{ /* Error appended ! */

	  char s[256];
	  _imagix_err_ext=errno;
	  sprintf(s,"Error occured:\n  Cannot save mri in file >%s<\n",file);
	  if(_dialog_messages_status&0x00000002l)
	    { PUT_ERROR(s);
	    }
	  fclose(hand);
	  return(0);
	}
    }
  else { /* File exists */

    char s[256];
    sprintf(s,"Delete old file ?\n  (%s)\n",file);
    if(XGET_EXIST(s, &err)==0)
      {
	sprintf(s,"Cannot save mri in file >%s<\n",file);
	PUT_ERROR(s);
	fclose(hand);
	return(0);
      }
  }



  if(pos>0 && pos<=MAX_PICTURES)
    {


      img_prt=ptr_img(pos);
      if((hand=fopen(file,"wb"))==NULL) { 
	char s[256];
	sprintf(s,"Cannot write mri in file>%s<.\n",file);
	PUT_ERROR(s);
	fclose(hand);
	return(0);
      }

      /*  Header     */
      fseek(hand, SEEK_SET, 0);
      /* Paintbrush */

      put_byte(0x0a,hand);     /* 10 == PC-Paintbrush version  */
      put_byte(0x03,hand);     /* 5 == version 3.0 with palette info */
      put_byte(0x01,hand);     /* 1 == .PCX run-length encoding  */
      put_byte(8,hand);        /* 8 bits per pixel  */

      /* Window */

      put_word(0,hand);        /* Min x */ 
      put_word(0,hand);        /* Min y*/ 
      if (!_is_bigEndian) 
	{
	  put_word(shEndianConversion(img_prt->width)-1,hand);      /* Max x */ 
	  put_word(shEndianConversion(img_prt->height)-1,hand);      /* Max y */ 
	} 
      else 
	{ 
	  put_word(img_prt->width-1,hand);      /* Max x */
	  put_word(img_prt->height-1,hand);      /* Max y */ 
	} /* Colors */
      if (!_is_bigEndian) 
	{ 
	  put_word(shEndianConversion(256),hand); /* X dpi */ 
	  put_word(shEndianConversion(256),hand);      /* Y dpi */
	} 
      else 
	{ 
	  put_word(256,hand); 
	  put_word(256,hand); 
	}

      for (i=0;i<16;++i) 
	{
	  put_byte(colors[i].red,hand);
	  put_byte(colors[i].green,hand);
	  put_byte(colors[i].blue,hand);
	}
      for (;i<16;++i)  
	{
	  put_byte(255,hand);
	  put_byte(255,hand);
	  put_byte(255,hand);
	}
      put_byte(0,hand);        /* reserved */
      put_byte(1,hand);        /* Number of color planes */
      if (!_is_bigEndian) 
	{
	  put_word(shEndianConversion(255),hand);      /* Max x */
	  put_word(shEndianConversion(1),hand);      /* Max x */	  	
	}
      else 
	{
	  put_word(255,hand);      /* bytes per line */
	  put_word(1,hand);        /* palette info */
	}
      /* fill the end of header */

      for(i=0;i<58;++i) {  
	put_byte(0,hand); 
      }

      /* Write image */
      printf("cutoff_min:%ld,%f cutoff_max:%ld,%f\n",img_prt->cutoff_min,img_prt->cutoff_min*img_prt->rcoeff,img_prt->cutoff_max,img_prt->cutoff_max*img_prt->rcoeff);
      for (j=0;j<img_prt->height;j++)
	{ 
	  for (i=0;i<img_prt->width;i++)
	    { 
	      line[i]=(int) (((img_prt->mri[i][j]-img_prt->cutoff_min)/(float)(img_prt->cutoff_max-img_prt->cutoff_min)) *255);
	    }
	  put_line(line,hand);
	}

      /* Write palette */

      put_byte(0x0c,hand);
      for (i=0;i<256;++i) {
	put_byte(colors[i].red,hand);
	put_byte(colors[i].green,hand);
	put_byte(colors[i].blue,hand);

      }
      for (;i<256;++i)  {
	put_byte(255,hand);
	put_byte(255,hand);
	put_byte(255,hand);
      }
      fclose(hand);
    }
  return(0);


}
#endif /* __GTK_GUI */

FILE*
p_openr(char *name)
{
  FILE* f;

#if defined (WIN32)
  if ( strcmp( name, "-" ) == 0 )
    f = stdin;
  else
    {
      f = fopen( name, "rb" );
      if ( f == NULL )
	{
	  if(MESG_DEBUG) printf(" name ");
	  return(NULL);
	}
    }
  return f;   
#else
  if ( strcmp( name, "-" ) == 0 )
    f = stdin;
  else
    {
      f = fopen( name, "r" );
      if ( f == NULL )
	{
	  if(MESG_DEBUG) printf(" name ");
	  return(NULL);
	}
    }
  return f;
#endif
}

/*---------------------------------------------------------*/
void
p_close(FILE *f)
{
  fflush( f );
  if ( ferror( f ) )
    printf( "a file read or write error occurred at some point" );
  if ( f != stdin )
    if ( fclose( f ) != 0 )
      printf( "fclose" );
}

#ifdef __GTK_GUI
/*
 * convert packed pixel format into 1 pixel per byte
 */
void
pcx_unpack_pixels(unsigned char *pixels, unsigned char *bitplanes, int bytesperline, int planes, int bitsperpixel)
{
  register int        bits;

  if (planes != 1)
    {
      PUT_MESG("can't handle packed pixels with more than 1 plane\n");
      return;
    }
  if (bitsperpixel == 8)
    {
      while (--bytesperline >= 0)
	*pixels++ = *bitplanes++;
    }
  else if (bitsperpixel == 4)
    {
      while (--bytesperline >= 0)
	{
	  bits        = *bitplanes++;
	  *pixels++   = (bits >> 4) & 0x0f;
	  *pixels++   = (bits     ) & 0x0f;
	}
    }
  else if (bitsperpixel == 2)
    {
      while (--bytesperline >= 0)
	{
	  bits        = *bitplanes++;
	  *pixels++   = (bits >> 6) & 0x03;
	  *pixels++   = (bits >> 4) & 0x03;
	  *pixels++   = (bits >> 2) & 0x03;
	  *pixels++   = (bits     ) & 0x03;
	}
    }
  else if (bitsperpixel == 1)
    {
      while (--bytesperline >= 0)
	{
	  bits        = *bitplanes++;
	  *pixels++   = ((bits & 0x80) != 0);
	  *pixels++   = ((bits & 0x40) != 0);
	  *pixels++   = ((bits & 0x20) != 0);
	  *pixels++   = ((bits & 0x10) != 0);
	  *pixels++   = ((bits & 0x08) != 0);
	  *pixels++   = ((bits & 0x04) != 0);
	  *pixels++   = ((bits & 0x02) != 0);
	  *pixels++   = ((bits & 0x01) != 0);
	}
    }
}

/*******************************************
 ** --  loadgif() ----------------------------
 **
 ** Lecture fichier GIF
 ********************************************/
void
loadgif(char *file, int pos, int reduc)
{

  FILE        *in;
  int        imageNumber;

  imageNumber = 1; /* Image number */
  verbose = FALSE; /*  verbose */
  showComment = FALSE; /* -comments */

  in = p_openr(file);
  ReadGIF(in, imageNumber, pos, reduc, file);
  p_close(in);

}

void
ReadGIF(FILE *fd, int imageNumber, int pos, int reduc, char *file)
{
  unsigned char    buf[16];
  unsigned char    c;
  unsigned char    localColormap[3][MAXCOLORMAPSIZE];
  int        useGlobalColormap;
  int        bitPixel;
  int        imageCount = 0;
  char        version[4];
  char            mess[80];

  if (! ReadOK(fd,buf,6))
    {
      PUT_MESG("error reading magic number\n");
      return;
    }

  if (strncmp((char *)buf,"GIF",3) != 0)
    {
      PUT_MESG("not a GIF file\n");
      return;
    }

  strncpy(version, (char *)buf + 3, 3);
  version[3] = '\0';

  if ((strcmp(version, "87a") != 0) && (strcmp(version, "89a") != 0))
    {
      PUT_MESG("bad version number, not '87a' or '89a'\n");
      return;
    }

  if (! ReadOK(fd,buf,7))
    {
      PUT_MESG("failed to read screen descriptor\n");
      return;
    }

  GifScreen.Width           = LM_to_uint(buf[0],buf[1]);
  GifScreen.Height          = LM_to_uint(buf[2],buf[3]);
  GifScreen.BitPixel        = 2<<(buf[4]&0x07);
  GifScreen.ColorResolution = (((buf[4]&0x70)>>3)+1);
  GifScreen.Background      = buf[5];
  GifScreen.AspectRatio     = buf[6];

  if (BitSet(buf[4], LOCALCOLORMAP)) {    /* Global Colormap */
    if (ReadColormap(fd,GifScreen.BitPixel,GifScreen.Colormap))
      {
	PUT_MESG("error reading global colormap\n");
	return;
      }
  }

  if (GifScreen.AspectRatio != 0 && GifScreen.AspectRatio != 49) {
    float    r;
    r = ( (float) GifScreen.AspectRatio + 15.0 ) / 64.0;
    if(MESG_DEBUG) printf("warning - non-square pixels; to fix do a 'pnmscale -cscale g'");
  }

  for (;;) {
    if (! ReadOK(fd,&c,1))
      {
	PUT_MESG("EOF / read error on image data\n");
	return;
      }

    if (c == ';') {        /* GIF terminator */
      if (imageCount < imageNumber)
	{
	  sprintf(mess,"only %d image%s found in file \n",
		  imageCount, imageCount>1?"s":"" );
	  PUT_MESG(mess);
	  return;
	}
      return;
    }

    if (c == '!') {     /* Extension */
      if (! ReadOK(fd,&c,1))
	{
	  PUT_MESG("EOF / read error on extention function code\n");
	  return;
	}
      DoExtension(fd, c);
      continue;
    }

    if (c != ',') {        /* Not a valid start character */
      if(MESG_DEBUG) printf("bogus character 0x%02x, ignoring",c);
      continue;
    }

    ++imageCount;

    if (! ReadOK(fd,buf,9))
      {
	PUT_MESG("couldn't read left/top/width/height\n");
	return;
      }

    useGlobalColormap = ! BitSet(buf[8], LOCALCOLORMAP);

    bitPixel = 1<<((buf[8]&0x07)+1);

    if (! useGlobalColormap) {
      if (ReadColormap(fd, bitPixel, localColormap))
	{
	  PUT_MESG("error reading local colormap\n");
	  return;
	}
      /*protoize:???*/
      ReadImage(fd, LM_to_uint(buf[4],buf[5]),
		LM_to_uint(buf[6],buf[7]), localColormap,
		BitSet(buf[8], INTERLACE), imageCount != imageNumber, pos, reduc ,NULL /*???*/);
    } else {
      ReadImage(fd, LM_to_uint(buf[4],buf[5]),
		LM_to_uint(buf[6],buf[7]), GifScreen.Colormap,
		BitSet(buf[8], INTERLACE), imageCount != imageNumber, pos, reduc, file);
    }

  }
}

int
ReadColormap(FILE *fd, int number, unsigned char (*buffer)[256])
{
  int        i;
  unsigned char    rgb[3];

  for (i = 0; i < number; ++i) {
    if (! ReadOK(fd, rgb, sizeof(rgb)))
      {
	PUT_MESG("bad colormap\n");
	return(TRUE);
      }

    buffer[CM_RED][i] = rgb[0] ;
    buffer[CM_GREEN][i] = rgb[1] ;
    buffer[CM_BLUE][i] = rgb[2] ;
  }
  return FALSE;
}

int
DoExtension(FILE *fd, int label)
{
  static char    buf[256];
  char        *str;

  switch (label) {
  case 0x01:        /* Plain Text Extension */
    str = "Plain Text Extension";
#ifdef notdef
    if (GetDataBlock(fd, (unsigned char*) buf) == 0)
      ;

    lpos   = LM_to_uint(buf[0], buf[1]);
    tpos   = LM_to_uint(buf[2], buf[3]);
    width  = LM_to_uint(buf[4], buf[5]);
    height = LM_to_uint(buf[6], buf[7]);
    cellw  = buf[8];
    cellh  = buf[9];
    foreground = buf[10];
    background = buf[11];

    while (GetDataBlock(fd, (unsigned char*) buf) != 0) {
      P_ASSIGN(image[ypos][xpos],
	       cmap[CM_RED][v],
	       cmap[CM_GREEN][v],
	       cmap[CM_BLUE][v]);
      ++index;
    }

    return FALSE;
#else
    break;
#endif
  case 0xff:        /* Application Extension */
    str = "Application Extension";
    break;
  case 0xfe:        /* Comment Extension */
    str = "Comment Extension";
    while (GetDataBlock(fd, (unsigned char*) buf) != 0) {
      if ((showComment)
	  &&(MESG_DEBUG)) printf("gif comment: %s", buf );
    }
    return FALSE;
  case 0xf9:        /* Graphic Control Extension */
    str = "Graphic Control Extension";
    (void) GetDataBlock(fd, (unsigned char*) buf);
    Gif89.disposal    = (buf[0] >> 2) & 0x7;
    Gif89.inputFlag   = (buf[0] >> 1) & 0x1;
    Gif89.delayTime   = LM_to_uint(buf[1],buf[2]);
    if ((buf[0] & 0x1) != 0)
      Gif89.transparent = buf[3];

    while (GetDataBlock(fd, (unsigned char*) buf) != 0)
      ;
    return FALSE;
  default:
    str = buf;
    sprintf(buf, "UNKNOWN (0x%02x)", label);
    break;
  }

  if(MESG_DEBUG) printf("got a '%s' extension", str );

  while (GetDataBlock(fd, (unsigned char*) buf) != 0)
    ;

  return FALSE;
}

int    ZeroDataBlock = FALSE;

int
GetDataBlock(FILE *fd, unsigned char *buf)
{
  unsigned char    count;

  if (! ReadOK(fd,&count,1)) {
    if(MESG_DEBUG) printf("error in getting DataBlock size" );
    return -1;
  }

  ZeroDataBlock = count == 0;

  if ((count != 0) && (! ReadOK(fd, buf, count))) {
    if(MESG_DEBUG) printf("error in reading DataBlock" );
    return -1;
  }

  return count;
}

int
GetCode(FILE *fd, int code_size, int flag)
{
  static unsigned char    buf[280];
  static int        curbit, lastbit, done, last_byte;
  int            i, j, ret;
  unsigned char        count;

  if (flag) {
    curbit = 0;
    lastbit = 0;
    done = FALSE;
    return 0;
  }

  if ( (curbit+code_size) >= lastbit) {
    if (done) {
      if (curbit >= lastbit)
	{
	  PUT_MESG("ran off the end of my bits\n");
	  return(-1);
	}
      return -1;
    }
    buf[0] = buf[last_byte-2];
    buf[1] = buf[last_byte-1];

    if ((count = GetDataBlock(fd, &buf[2])) == 0)
      done = TRUE;

    last_byte = 2 + count;
    curbit = (curbit - lastbit) + 16;
    lastbit = (2+count)*8 ;
  }

  ret = 0;
  for (i = curbit, j = 0; j < code_size; ++i, ++j)
    ret |= ((buf[ i / 8 ] & (1 << (i % 8))) != 0) << j;

  curbit += code_size;

  return ret;
}

int
LWZReadByte(FILE *fd, int flag, int input_code_size)
{
  static int    fresh = FALSE;
  int        code, incode;
  static int    code_size, set_code_size;
  static int    max_code, max_code_size;
  static int    firstcode, oldcode;
  static int    clear_code, end_code;
  static int    table[2][(1<< MAX_LWZ_BITS)];
  static int    stack[(1<<(MAX_LWZ_BITS))*2], *sp;
  register int    i;

  if (flag) {
    set_code_size = input_code_size;
    code_size = set_code_size+1;
    clear_code = 1 << set_code_size ;
    end_code = clear_code + 1;
    max_code_size = 2*clear_code;
    max_code = clear_code+2;

    GetCode(fd, 0, TRUE);

    fresh = TRUE;

    for (i = 0; i < clear_code; ++i) {
      table[0][i] = 0;
      table[1][i] = i;
    }
    for (; i < (1<<MAX_LWZ_BITS); ++i)
      table[0][i] = table[1][0] = 0;

    sp = stack;

    return 0;
  } else if (fresh) {
    fresh = FALSE;
    do {
      firstcode = oldcode =
	GetCode(fd, code_size, FALSE);
    } while (firstcode == clear_code);
    return firstcode;
  }

  if (sp > stack)
    return *--sp;

  while ((code = GetCode(fd, code_size, FALSE)) >= 0) {
    if (code == clear_code) {
      for (i = 0; i < clear_code; ++i) {
	table[0][i] = 0;
	table[1][i] = i;
      }
      for (; i < (1<<MAX_LWZ_BITS); ++i)
	table[0][i] = table[1][i] = 0;
      code_size = set_code_size+1;
      max_code_size = 2*clear_code;
      max_code = clear_code+2;
      sp = stack;
      firstcode = oldcode =
	GetCode(fd, code_size, FALSE);
      return firstcode;
    } else if (code == end_code) {
      int        count;
      unsigned char    buf[260];

      if (ZeroDataBlock)
	return -2;

      while ((count = GetDataBlock(fd, buf)) > 0)
	;

      if (count != 0)
	if(MESG_DEBUG) printf("missing EOD in data stream (common occurence)");
      return -2;
    }

    incode = code;

    if (code >= max_code) {
      *sp++ = firstcode;
      code = oldcode;
    }

    while (code >= clear_code) {
      *sp++ = table[1][code];
      if (code == table[0][code])
	{
	  PUT_MESG("circular table entry BIG ERROR\n");
	  return(-1);
	}
      code = table[0][code];
    }

    *sp++ = firstcode = table[1][code];

    if ((code = max_code) <(1<<MAX_LWZ_BITS)) {
      table[0][code] = oldcode;
      table[1][code] = firstcode;
      ++max_code;
      if ((max_code >= max_code_size) &&
	  (max_code_size < (1<<MAX_LWZ_BITS))) {
	max_code_size *= 2;
	++code_size;
      }
    }

    oldcode = incode;

    if (sp > stack)
      return *--sp;
  }
  return code;
}

void
ReadImage(FILE *fd, int len, int height, unsigned char (*cmap)[256], int interlace, int ignore, int n, int reduc, char *file)
{
  unsigned char    c;
  int        v;
  int        xpos = 0, ypos = 0, pass = 0;
  /*-----------------------------------------*/
  int im_res[2*MAX_PICTURES];
  int i,j,k,l;
  pixel val;
  int pasx,pasy,pasposx,pasposy;
  grphic *imres[2*MAX_PICTURES];
  char tab[80];

  sprintf(tab,">%s< GIF\ninfo : Taille de l'image GIF %dx%d \n",file,len,height);
  PUT_MESG(tab);

  im_res[0]=n;
  pasx= GINT(GINT(len,MAX_WIDTH),reduc);
  pasy= GINT(GINT(height,MAX_HEIGHT),reduc);
  if (pasx*pasy> 2*MAX_PICTURES)
    {
      PUT_MESG("Erreur : Taille de l'image trop grande\n");
      return;
    }

  for (j=0;j<pasy;j++)
    {
      l=n+j*MAX_PIC_IN_LINE;
      for (i=0;i<pasx;i++)
	im_res[i+j*pasx]=(l>((GINT(n,MAX_PIC_IN_LINE)+j)*MAX_PIC_IN_LINE)||l>MAX_PICTURES) ? 0 :l++;
    }

  for (j=0;j<pasy;j++)
    for (i=0;i<pasx;i++)
      if(im_res[i+j*pasx]!=0)
	imres[i+j*pasx]=ptr_img(im_res[i+j*pasx]);
  /*-----------------------------------------*/

  /*
  **  Initialize the Compression routines
  */
  if (! ReadOK(fd,&c,1))
    {
      PUT_MESG("EOF / read error on image data\n");
      return;
    }

  if (LWZReadByte(fd, TRUE, c) < 0)
    {
      PUT_MESG("error reading image\n");
      return;
    }

  /*
  **  If this is an "uninteresting picture" ignore it.
  */
  if (ignore) {
    if ((verbose)&&(MESG_DEBUG))
      printf("skipping image..." );

    while (LWZReadByte(fd, FALSE, c) >= 0)
      ;
    return;
  }


  if ((verbose)&&(MESG_DEBUG))
    printf("reading %d by %d%s GIF image",
	   len, height, interlace ? " interlaced" : "" );

  while ((v = LWZReadByte(fd,FALSE,c)) >= 0 ) {
    P_ASSIGN(val, cmap[CM_RED][v],
	     cmap[CM_GREEN][v], cmap[CM_BLUE][v]);
    /*==========================>*/



    pasposx= GINT(GINT(xpos+1,MAX_WIDTH),reduc);
    pasposy= GINT(GINT(ypos+1,MAX_HEIGHT),reduc);
    pasposx--;
    pasposy--;
    k=pasposx+(pasx*pasposy);
    i=(xpos-pasposx*MAX_WIDTH*reduc)/reduc;
    j=(ypos-pasposy*MAX_HEIGHT*reduc)/reduc;
    if(im_res[k]!=0)
      imres[k]->mri[i][j]=
	(gray)(P_LUMIN(val)+0.5);

    /*==========================>*/

    ++xpos;
    if (xpos == len) {
      xpos = 0;
      if (interlace) {
	switch (pass) {
	case 0:
	case 1:
	  ypos += 8; 
	  break;
	case 2:
	  ypos += 4; 
	  break;
	case 3:
	  ypos += 2; 
	  break;
	}

	if (ypos >= height) {
	  ++pass;
	  switch (pass) {
	  case 1:
	    ypos = 4; 
	    break;
	  case 2:
	    ypos = 2; 
	    break;
	  case 3:
	    ypos = 1; 
	    break;
	  default:
	    goto fini;
	  }
	}
      } else {
	++ypos;
      }
    }
    if (ypos >= height)
      break;

  }

  /*===========>*/

 fini:
  if (LWZReadByte(fd,FALSE,c)>=0)
    if(MESG_DEBUG) printf("too much input data, ignoring extra...");

  if ((verbose)&&(MESG_DEBUG))
    printf("Fichier GIF...\n");

  /*==========================>*/
  for (ypos=0; ypos<(reduc*pasy*MAX_HEIGHT);ypos++)
    for (xpos=0; xpos< (reduc*pasx*MAX_WIDTH); xpos++)
      {
	if(xpos>= len || ypos >= height)
	  {
	    pasposx= GINT(GINT(xpos+1,imres[0]->width),reduc);
	    pasposy= GINT(GINT(ypos+1,imres[0]->height),reduc);
	    pasposx--;
	    pasposy--;
	    k=pasposx+(pasx*pasposy);
	    i=(xpos-pasposx*MAX_WIDTH*reduc)/reduc;
	    j=(ypos-pasposy*MAX_HEIGHT*reduc)/reduc;
	    if(im_res[k]!=0)
	      imres[k]->mri[i][j]=0;
	  }
      }

  for (j=0;j<pasy;j++)
    for (i=0;i<pasx;i++)
      {
	imx_iniparam(im_res[i+j*pasx],255,0,0.,1.);
	show_picture(im_res[i+j*pasx]);
      }

  /*==========================>*/

}




/***********************************************************************
 **
 **     Procedure  Lecture_nimg_ds_nWins(n,norm)         HMO 4/3/96
 **    int n    : Dimension de la fenetre souhaite 
 **           (4 pour une fen. avec 16 images)
 **      int norm : si != 0 alors il y normalisation ( 65535 )
 **
 **     Description:
 **     Permet de charger & visualiser n images dans plusieurs fenetres
 **    La taille des images dependra de n & leurs disposition du nombre
 **    total images a visualiser.
 **
 **     Remarque:
 **    La procedure fait appel a charge_petite_picture(...)
 **
 ***********************************************************************/

int Lecture_nimg_ds_nWins(int n, int norm)
{
  int max,pos,nbimage,nbligne,i,j,k,l=1,c=0;
  grphic *img_tmp;

  if( !imx_serie_file_creation(3) )   return(0);
  pos=GET_WND(NULL);

  max=_nbpts_serie;
  nbimage=max/(n*n);
  if (max-nbimage*n*n) nbimage++;
  if (nbimage+pos > 12)
    {
      nbimage =12-pos+1 ;
      if (max > nbimage*n*n) max=nbimage*n*n;
    }
  nbligne=max/(nbimage*n);
  if (max-nbligne*nbimage*n) nbligne++;
  if (PROG_DEBUG) printf ("NBLIGNE=%d    NBIMAGE=%d \n",nbligne,nbimage);

  for (k=pos;k<(pos+nbimage);k++)
    {
      img_tmp=ptr_img(k);
      for (j=0;j<MAX_HEIGHT;j++)
	for (i=0;i<MAX_WIDTH;i++)
	  img_tmp->mri[i][j]=(long)0;
    }

  do
    {
      for ( i=pos ; i<(pos+nbimage) && (c<max) ; i++ )
	for ( j=(1+(l-1)*n) ; (j<=n+(l-1)*n) && (c<max) ; j++ )
	  {
	    charge_petite_picture(i,_serie_file[c],_serie_file_number[c],n,j,norm);
	    c++;
	    if (PROG_DEBUG) printf ("compteur= %d    pos= %d    posdsimg= %d\n",c,i,j);
	  }
      l++;
    }    while ( (l<=nbligne) && (c<max) );

  for (i=pos;i<(pos+nbimage);i++)
    {
      imx_inimaxminpixel(i);
      show_picture(i);
    }
  return(1);
}


/***********************************************************************
 **
 **     Procedure  Lecture_1img_ds_1Win(n,norm)         HMO 4/3/96
 **    int n    : Dimension de la fenetre souhaite 
 **           (4 pour une fen. avec 16 images)
 **      int norm : si != 0 alors il y normalisation ( 65535 )
 **
 **     Description:
 **     Permet de charger & visualiser 1 image dans une fenetres
 **    La taille de l'image dependra de n & sa position est demander en
 **      indiquant la fourchette
 **
 **     Remarque:
 **    La procedure fait appel a charge_petite_picture(...)
 **
 ***********************************************************************/

int Lecture_1img_ds_1Win(int n, int norm)
{
  char *file,str[20];
  int e,pos,posdsimg,i,j,numfirst;
  grphic *img_tmp;

  file=GET_FILE_IMG(NULL, &numfirst, &e);
  if (e) return (0);
  pos=GET_WND(NULL);
  sprintf (str,"De 1 a %d",n*n);
  do {
    posdsimg= GET_INT(str, 0, &e);
  }
  while ( posdsimg>(n*n) || posdsimg < 1 );
  if (e)  return (0);


  /* Normalisation prealable */
  if (norm)
    {
      img_tmp=ptr_img(pos);
      if (img_tmp->max_pixel)
	for (j=0;j<MAX_HEIGHT;j++)
	  for (i=0;i<MAX_WIDTH;i++)
	    img_tmp->mri[i][j]=(TYPEMRI)(img_tmp->mri[i][j]*65535./img_tmp->max_pixel);
    }

  /*charge_petite_picture (pos,file,n,posdsimg,norm); */
  charge_petite_picture(pos,file,numfirst,n,posdsimg,norm);
  imx_inimaxminpixel(pos);
  show_picture(pos);
  return(1);
}

/***********************************************************************
 **
 **     Procedure  charge_petite_picture(...)         HMO 4/3/96
 **    int num_image : entier qui correspond a la fenetre
 **    char *file : nom complet du fichier image
 **    int n    : Dimension de la fenetre souhaite 
 **           (4 pour une fen. avec 16 images)
 **      int posdsimg : entier qui defini la position de l'image dans la fenetre 
 **                La numerotation se fait ligne apres ligne
 **               (expl pour n=2):          *********
 **                        * 1 * 2 *
 **                        *********
 **                        * 3 * 4 *
 **                        *********    
 **      int norm : si != 0 alors il y normalisation ( 65535 )
 **
 **     Description:
 **     Permet de charger 1 image dans une fenetres
 **    La taille de l'image dependra de n & sa position est demander en
 **      indiquant la fourchette
 **
 ***********************************************************************/

void charge_petite_picture (int num_image, char *file, int numfirst, int n, int posdsimg, int norm)
{
  grphic *im_tmp,*im1;
  int i,j,reste,l,c;

  im_tmp=cr_grphic(NULL);
  if(load_mri(file,(grphic *)im_tmp,numfirst)!=NULL)
    {
      im1=ptr_img(num_image);
      im1->rcoeff=1;
      im1->icomp=0;
      im1->dx=0;
      im1->dy=0;
      im1->dz=0;
      sprintf(im1->filename,"%c",'\0');


      reste=posdsimg-(posdsimg/n)*n;
      if (reste) l=(posdsimg+(n-reste))/n;
      else l=posdsimg/n;
      c=posdsimg-(l-1)*n;

      for (j=0;j<(MAX_HEIGHT/n);j++)
	for (i=0;i<(MAX_WIDTH/n);i++)
	  {
	    if (norm)
	      im1->mri[j+(c-1)*(MAX_HEIGHT/n)][i+(l-1)*(MAX_WIDTH/n)]=(TYPEMRI)(im_tmp->mri[j*n][i*n]*65535./im_tmp->max_pixel);
	    else 
	      im1->mri[j+(c-1)*(MAX_HEIGHT/n)][i+(l-1)*(MAX_WIDTH/n)]=im_tmp->mri[j*n][i*n];
	  }
      free_grphic(im_tmp);
    }
}




/**********************************************************
 *** int read_raw ()            HMO 310596
 ***
 ***    Permet de lire une image rawdata quelque soit 
 ***     taille de celle ci
 ***
 ***     PARAMETRES
 ***        file : nom de l'image
 ***        pos  : numero de la fenetre pour le resultat        
 ***        dimx,dimy : dimension de l'image
 ***        taille : nombre d'octets par point (1byte, 2bytes etc)
 ***             taille doit etre inferieure a 5 !!!
 ***    REMARQUE
 ***        La fonction teste la compatibilite entre 
 ***        dimx,dimy et taille d'une part et la taille
 ***        en octet de l'image d'autre part.
 ***
 ***********************************************************/

int read_raw (char *file, int dimx, int dimy, int pos, int taille)
{
  int i,j,l,size;
  unsigned char *buff8b;
  long *buff32b;
  short *buff16b;
  grphic *img;
  FILE *fp;

  if ( (fp=fopen(file,"rb"))==NULL ) {
    printf("Le fichier %s n'existe pas \n",file);
  }
  fseek(fp,0l,SEEK_END);
  size=ftell(fp);
  fseek(fp,0l,SEEK_SET);

  if ( fp==NULL || size!=(dimx*dimy*taille) )
    {
      PUT_MESG("Reading error or not in the right format\n");
      return (-1);
    }

  img=ptr_img(pos);

  switch (taille) {
  case 1:  /*On lit un fichier 8 bits*/
    /*allocation du buffer a la bonne taille*/
    if((buff8b=CALLOC(dimx*dimy,unsigned char))==NULL) {
      /* Cannot alloc... */
      fclose(fp);
      return(0);
    }

    size=fread((unsigned char*)buff8b,sizeof(unsigned char)
	       ,dimx*dimy,fp);
    l=0;
    for (j=0;j<dimy;j++)
      for (i=0;i<dimx;i++) {
	img->mri[i][j]=buff8b[l++];
      }
    FREE(buff8b);
    break;

  case 2:  /*On lit un fichier 16 bits*/
    /*allocation du buffer a la bonne taille*/
    if((buff16b=CALLOC(dimx*dimy,short))==NULL) {
      /* Cannot alloc... */
      fclose(fp);
      return(0);
    }

    size=FREAD_SHORT((short*)buff16b,dimx*dimy,fp);
    l=0;
    for (j=0;j<dimy;j++)
      for (i=0;i<dimx;i++) {
	img->mri[i][j]=buff16b[l++];
      }
    FREE(buff16b);
    break;

  case 4:  /*On lit un fichier 32 bits*/
    /*allocation du buffer a la bonne taille*/
    if((buff32b=CALLOC(dimx*dimy,long))==NULL) {
      /* Cannot alloc... */
      fclose(fp);
      return(0);
    }

    size=FREAD_LONG((long*)buff32b,dimx*dimy,fp);
    l=0;
    for (j=0;j<dimy;j++)
      for (i=0;i<dimx;i++) {
	img->mri[i][j]=buff32b[l++];
      }
    FREE(buff32b);
    break;

  }


  img->rcoeff=1;
  imx_inimaxminpixel (pos);
  img->width=dimx;
  img->height=dimy;
  strcpy (img->filename,file);
  show_picture(pos);
  return(0);
}

int     load_mri_ipb_2d(char *file, grphic *image, int pos)
{
  printf( "[load_mri_ipb_2d] this function is no longer defined !!!\n\
Please use \"load_mri\"\n");

  return(0);
}

/*--------  compressave_serial()   ----------------------------------
***
***	Compression et sauvegarde des images
***		
***	Compression d'une s�ie d'images et sauvegarde sur disquette:
***	- nom de la s�ie d'images
***	- v�ification par animation
***	- compression de la s�ie d'images par imx_compres_serial dans imx_zone.c
***	- sauvegarde avec l'extension .C001 etc sur disquette
***	et dans le r�ertoire temp de coeur.
***
------------------------------------------------------------------*/
int	compressave_serial(void)
{
  int pos1=0,pos2,number,step,i,j,n=0,e,numfirst;
  char nom[5],*file,img_name[FILE_LEN],cmd[FILE_LEN];
  char filename[FILE_LEN],filerest[FILE_LEN];
  grphic *mask,*image;

  PUT_WNDMSG(TEXT0135);
  file=GET_FILE_IMG(NULL, &numfirst, &e);
  CLEAR_WNDMSG();
  if (e) return(-1);
  if(file==NULL || strlen(file)<=0)
    {
      PUT_WARN("Sorry, cannot read picture !\n");
      return(-1);
    }

  if (is_mri(file,NULL) != PARAVISION_FORMAT)
    {
      /*   recherche image num�o 1 de la s�ie   */

      for(i=0;file[i]!='.';i++)
	{
	  filename[i]=file[i];
	  //filename[i+1]=NULL;  <- correction, 0 et non NULL.  RC octobre 2001
	  filename[i+1]=0;
	}
      strcpy(file,filename);
      strcat(file,".I001");

      /*   extraction du nom de l'image    */

      for (i=(strlen(file));file[i]!='/';i--)
	n=i;
      for (i=n;file[i]!='.'&&i < n+6 ;i++)
	{
	  img_name[i-n]=file[i];
	  //img_name[i-n+1]=NULL;  <- correction, 0 et non NULL.  R. Chevrier octobre 2001
	  img_name[i-n+1]=0;
	}

      /*   fichiers des images r�ultat   */

      strcpy(filerest,_home);
      strcat(filerest,"/");
      strcat(filerest,img_name);
      strcat(filerest,"01.PCX");
    }
  else
    {
      image=ptr_img(pos1);
      load_mri(file,image,numfirst);
      for (i=0;i<4;i++) nom[i]=image->patient.name[i];
      nom[4]='\0';
      strcpy(filerest,_home);
      strcat(filerest,"/");
      strcat(filerest,nom);
      strcat(filerest,_path_paravision);
      strcat(filerest,"01.PCX");
    } 

  /*   choix du nombre d'images, du pas et des positions  */

  /*	PUT_WNDMSG(TEXT0066);
	pos1=GET_WND(NULL);*/
  pos1=1;
  /*	PUT_WNDMSG(TEXT0006);
	pos2=GET_WND(NULL);*/
  pos2=2;
  /*	CLEAR_WNDMSG();*/
  number= GET_INT(TEXT0112, 0, &e);
  if (e) return(-1);
  /*	step= GET_INT(TEXT0024, 0, &e);
	if (e) return(-1);*/
  step=1;
	

  /*   suppression du bruit */  
  /*   Initialisation du masque */

  mask=ptr_img(pos2);
  imx_iniparaimg(pos2);
  mask->max_pixel=1;
  mask->min_pixel=0;
  mask->cutoff_max=1;
  mask->cutoff_min=0;
  for(j=0;j<_imx_img->height;j++)
    for(i=0;i<_imx_img->width;i++)
      mask->mri[i][j]=0;
		
  /*   visualisation image de depart  */ 

  image=ptr_img(pos1);
  if(load_mri(file,image,numfirst)==NULL)
    return(-1);
  show_picture(pos1);


  /*   cr�tion du masque   */

  if(imx_nonoise_mask(mask,image,pos2)!=1)	return(-1);

  /*   visualisation du masque pour confirmation   */

  show_picture(pos2);

  for (i=0;i<number;i++)
    {
      load_picture(pos1,file,numfirst);
      imx_and(pos1,pos2,pos1);
      save_mri_pcx(pos1,filerest);
      incfic_img(file,step,&numfirst);
      incfic(filerest,step);
    }

  strcpy(filerest,_home);
  strcat(filerest,"/");
  if (is_mri(file,NULL)!=PARAVISION_FORMAT)
    strcat(filerest,img_name); 
  else 
    { 
      strcat(filerest,nom);
      strcat(filerest,_path_paravision);
    } 
  strcat(filerest,"01.PCX");

  /*   animation pour confirmation   */
  /* Deleted Wiest 2004 cf deleted/cardio.c
     if(_cardio_save_zoom==TRUE)  init_zoom2_movie(pos1, filerest, 1);
  */
	
  _movie_cardiowithroi=FALSE;
  /* Deleted Wiest 2004 cf deleted/cardio.c

  IMXFile_Animate(pos1,number,step,filerest,
  1,NULL);
  */
  XGET_OKCANCEL(&e);
  if(e==1)
    {
      if(_cardio_save_zoom==TRUE)	unzoom_img(pos1);	
      strcpy(cmd,"rm ");
      strcat(cmd,_home);
      strcat(cmd,"/");
      if (is_mri(file,NULL)!=PARAVISION_FORMAT)
	strcat(cmd,img_name); 
      else 
	{ 
	  strcat(cmd,nom);
	  strcat(cmd,_path_paravision);
	} 
      strcat(cmd,"*.PCX ");
      system(cmd);

      return(-1);
    }

  /*   sauvegarde sur disquette  */

  strcpy(cmd,"mv ");
  strcat(cmd,_home);
  strcat(cmd,"/");
  if (is_mri(file,NULL)!=PARAVISION_FORMAT)
    strcat(cmd,img_name); 
  else 
    { 
      strcat(cmd,nom);
      strcat(cmd,_path_paravision);
    } 
  strcat(cmd,"*.PCX ");
  strcat(cmd,_home);
  strcat(cmd,"/floppy/");
  system(cmd);
  return(1);

}


/* ---	imx_nonoise_mask() -------------------------------
***
***
***	cr� un masque pour �iminer le bruit.
***
***	- l'image support pour la cr�tion du masque est _serie_file[0]
***	- le masque est stock�graphiquement dans mask 
***
--------------------------------------------------------------*/
int	imx_nonoise_mask(grphic *mask, grphic *im, int res_pos)
{
  int d,i=0,j,k,l,x,y,n,m,p,q,s,seuil;
  int hght,wdth,ntabgrow;
  float val,valbas,valhaut;
  float tabgrow[3];

  hght=im->height;
  wdth=im->width;

  /*   calcul de la valeur haute de seuillage   */

  s=0;
  q=0;
  for(j=0;j<hght/20;j++)
    for(i=0;i<wdth/20;i++)
      {	
	s=s+im->mri[i][j];
	q=q+1;
      }
  valhaut=3.5*s/q;
  valbas=0.;
  val=valhaut/2;

  for(k=0;k<4;k++)
    {
      if(k==0){
	j=hght/2;
	i=2;}
      if(k==1){
	j=2;
	i=wdth/2;}
      if(k==2){
	j=hght/2;
	i=wdth-2;}
      if(k==3){
	j=hght-2;
	i=wdth/2;}

      mask->mri[i][j]=1;
      tabgrow[0]=val;
      tabgrow[1]=valbas;
      tabgrow[2]=valhaut;
      ntabgrow=3;
      growing(i,i,im,mask,fct_grow_pixel,tabgrow,ntabgrow);
    }

  /*   suppression totale du bruit par zones   */

  d=32;
  n=wdth/d;
  m=hght/d;
  seuil=n*m*2/5;

  for(q=0;q<d;q++)
    for(p=0;p<d;p++)
      {
	s=0;
	x=p*n;
	y=q*m;
	for(j=y;j<(y+m);j++)
	  for(i=x;i<(x+n);i++)
	    if(mask->mri[i][j]==0)	s++;
	if(s<=seuil)
	  {	/*   mise �0 de la zone  */
	    for(j=y;j<(y+m);j++)
	      for(i=x;i<(x+n);i++)
		mask->mri[i][j]=0;
	  }
	else
	  {	/*   mise �1 de la zone   */
	    for(j=y;j<(y+m);j++)
	      for(i=x;i<(x+n);i++)
		mask->mri[i][j]=1;
	  }
      }

  /*   Suppression des zones noires enclav�s   */
	
  for(q=0;q<d;q++)
    for(p=0;p<d;p++)
      {
	s=0;
	if(p==0||p==(d-1)||q==0||q==(d-1))
	  seuil=3;
	else
	  seuil=5;

	if(p==0)
	  {
	    k=0;
	    l=1;
	  }
	else
	  {
	    if(p==(d-1))
	      {
		k=-1;
		l=0;
	      }
	    else
	      {
		k=-1;
		l=1;
	      }
	  }
	if(q==0)
	  {
	    x=0;
	    y=1;
	  }
	else
	  {
	    if(q==(d-1))
	      {
		x=-1;
		y=0;
	      }
	    else
	      {
		x=-1;
		y=1;
	      }
	  }

	if(mask->mri[p*n][q*m]==0)
	  {
	    for(i=k;i<(l+1);i++)
	      for(j=x;j<(y+1);j++)
		if(mask->mri[(p+i)*n][(q+j)*m]==1)	s++;
	    if(s>=seuil)
	      {
		for(j=q*m;j<(q*m+m);j++)
		  for(i=p*n;i<(p*n+n);i++)
		    mask->mri[i][j]=1;
	      }
	  }
	else
	  {
	    for(i=k;i<(l+1);i++)
	      for(j=x;j<(y+1);j++)
		if(mask->mri[(p+i)*n][(q+j)*m]==0)	s++;
	    if(s>=seuil+1)
	      {
		for(j=q*m;j<(q*m+m);j++)
		  for(i=p*n;i<(p*n+n);i++)
		    mask->mri[i][j]=0;
	      }
	  }
      }

  return(1);
}

/*---------------------------------------------------------*/
/* Variable-sized arrays. */

char*
allocrow(int cols, int size)
{
  register char* itrow;

  itrow = (char*) malloc( cols * size );
  if ( itrow == (char*) 0 )
    printf( "out of memory allocating a row" );
  return itrow;
}


/*---------------------------------------------------------*/
char**
allocarray(int cols, int rows, int size)
{
  char** its;
  int i;

  its = (char**) malloc( rows * sizeof(char*) );
  if ( its == (char**) 0 )
    printf( "out of memory allocating an array" );
  its[0] = (char*) malloc( rows * cols * size );
  if ( its[0] == (char*) 0 )
    printf( "out of memory allocating an array" );
  for ( i = 1; i < rows; ++i )
    its[i] = &(its[0][i * cols * size]);
  return its;
}


/*-----------------------------------------------------------------*/
void
read_pcx_image(FILE *fp, unsigned char *buf, int BytesPerLine, int Planes, int Height)
{
  int         c;
  int         nbytes;
  int         count;
  char        mess[80];

  nbytes      = BytesPerLine * Planes * Height;

  while (nbytes > 0)
    {
      c    = GetByte(fp);
      if ((c & 0xc0) != 0xc0)
	{
	  *buf++    = c;
	  --nbytes;
	  continue;
	}

      count    = c & 0x3f;
      c    = GetByte(fp);
      if (count > nbytes)
	{
	  sprintf(mess,"repeat count spans end of image, count = %d, nbytes = %d \n", count, nbytes);
	  PUT_MESG(mess);
	  return;
	}
      nbytes    -= count;
      while (--count >= 0)
	*buf++ = c;
    }
}

/* -- read_pcx_file() -------------------------------------------------------
**
**    Lit un fichier de type PCX et renvoie une structure contenant
**    les informations du header ainsi qu'un tableau contenant
**    la table de couleur et bien sur les pixels de l'image.
**
**    Une valeur de retour non nulle indique un probleme decrit par la
**    chaine status.
**
**    Usage: read_pcx_file(filename, descr, status)
**       filename : char*      : le nom du fichier pcx a lire
**       descr    : Pcx_Descr* : le descripteur du fichier
**       status   : char**     : pointeur sur une chaine decrivant une
**                               erreur eventuelle
**
**    SVP: Si *status n'est pas egal a NULL au retour, pensez a en faire
**         un free apres usage...
**
*/
unsigned char read_pcx_file(char *filename, Pcx_Descr *descr, char **status)
{
  FILE *ifp;
  char string[100];
  int  Xmin, Ymin, Xmax, Ymax;
  register int i;
  register int x, y;
  unsigned char *pcximage;
  unsigned char *pcxplanes;
  unsigned char *pcxpixels;
  unsigned char *tmppix;


  /* On essaye d'ouvrir le fichier */
  ifp = p_openr(filename);
  if (ifp == NULL){
    sprintf(string, "Unable to open file : %s\n", filename);
    *status = strdup(string);  
    return 1;
  }

  /* On lit le magic number (pr etre sur que c'est bien du PCX) */
  if (GetByte(ifp) != PCX_MAGIC){
    sprintf(string, "%s is not a PCX file \n", filename);
    *status = strdup(string);  
    return 1;
  }

  /* get version # */
  descr->Version = GetByte(ifp); 

  /* check for PCX run length encoding   */
  if (GetByte(ifp) != 1){
    sprintf(string, "%s has unknown encoding scheme \n", filename);
    *status = strdup(string);  
    return 1;
  }

  /* Read some header info */
  descr->BitsPerPixel= GetByte(ifp);
  Xmin = GetWord(ifp);
  Ymin = GetWord(ifp);
  Xmax = GetWord(ifp);
  Ymax = GetWord(ifp);

  descr->Width  = (Xmax - Xmin) + 1;
  descr->Height = (Ymax - Ymin) + 1;

  GetWord(ifp); /* ignore horizontal DPI resolution */
  GetWord(ifp); /* ignore vertical DPI resolution   */


  if (MESG_DEBUG){
    sprintf(string, ">%s< PCX\ninfo : Taille de l'image PCX %dx%d \n", filename, descr->Width, descr->Height);
    PUT_MESG( string );
  }


  /* On note la 16-color colormap (au cas ou...) */
  for (i=0; i<16; i++){
    descr->colormap[i][0] = GetByte(ifp);
    descr->colormap[i][1] = GetByte(ifp);
    descr->colormap[i][2] = GetByte(ifp);
  }

  GetByte(ifp);                      /* skip reserved byte  */
  descr->Planes = GetByte(ifp);       /* # of color planes   */
  descr->BytesPerLine = GetWord(ifp); /* # of bytes per line */
  GetWord(ifp);                      /* ignore palette info */

  
  /* Check that we can handle this image format */
  switch (descr->BitsPerPixel)
    {
    case 1:
      if (descr->Planes > 4){
	sprintf(string, "Can't handle image with more than 4 planes\n");
	*status = strdup(string);
	return 1;
      }
      break;

    case 2:
    case 4:
    case 8:
      if (descr->Planes == 1) break;
    default:
      {
	sprintf(string, "Can't handle %d bits per pixel image with %d planes\n", descr->BitsPerPixel, descr->Planes);
	*status = strdup(string);
	return 1;
      }
    }


  /* read the pcx format image */
  fseek(ifp, (long)PCX_HDR_SIZE, SEEK_SET);
  pcximage = (unsigned char *)allocrow(descr->BytesPerLine * descr->Planes, descr->Height);
  read_pcx_image(ifp, pcximage, descr->BytesPerLine, descr->Planes, descr->Height);


  pcxpixels = (unsigned char *)allocrow(descr->Width+7, 1);
  descr->pixels = (unsigned char *)calloc((long)descr->Width * (long)descr->Height, sizeof(unsigned char));
  tmppix = descr->pixels;

  /* convert the image */
  for (y=0; y<descr->Height; y++){
    pcxplanes = pcximage + (y * descr->BytesPerLine * descr->Planes);

    if (descr->Planes == 1){
      pcx_unpack_pixels(pcxpixels, pcxplanes, descr->BytesPerLine,
			descr->Planes, descr->BitsPerPixel);
    } else {
      pcx_planes_to_pixels(pcxpixels, pcxplanes, descr->BytesPerLine,
			   descr->Planes, descr->BitsPerPixel);
    }
    for (x=0; x<descr->Width; x++) *(tmppix++) = pcxpixels[x];
  }


  /* 256 color images have their color map at the end of the file
   * preceeded by a magic byte */
  if (descr->BitsPerPixel == 8) {
    if (GetByte(ifp) != PCX_256_COLORS){
      FREE(pcximage);
      sprintf(string, "Bad color map signature\n");
      *status = strdup(string);
      return 1;
    }
    for (i=0; i<256; i++){
      descr->colormap[i][0] = GetByte(ifp);
      descr->colormap[i][1] = GetByte(ifp);
      descr->colormap[i][2] = GetByte(ifp);
    }
  }


  p_close(ifp);
  FREE(pcxpixels);
  FREE(pcximage);

  *status = NULL; 
  return 0;
}



/*
 * convert multi-plane format into 1 pixel per byte
 */
void
pcx_planes_to_pixels(unsigned char *pixels, unsigned char *bitplanes, int bytesperline, int planes, int bitsperpixel)
{
  int  i, j;
  int  npixels;
  unsigned char    *p;

  if (planes > 4)
    {
      PUT_MESG("can't handle more than 4 planes\n");
      return;
    }
  if (bitsperpixel != 1)
    {
      PUT_MESG("can't handle more than 1 bit per pixel\n");
      return;
    }
  /*
   * clear the pixel buffer
   */
  npixels = (bytesperline * 8) / bitsperpixel;
  p    = pixels;
  while (--npixels >= 0)
    *p++ = 0;

  /*
   * do the format conversion
   */
  for (i = 0; i < planes; i++)
    {
      int pixbit, bits, mask;

      p    = pixels;
      pixbit    = (1 << i);
      for (j = 0; j < bytesperline; j++)
	{
	  bits = *bitplanes++;
	  for (mask = 0x80; mask != 0; mask >>= 1, p++)
	    if (bits & mask)
	      *p |= pixbit;
	}
    }
}

/*******************************************
 ** --  IMXFile_Capture ----------------------------
 **
 **  Capture d'une serie d'images.
 **  type : type de sauvegarde (PCX,TIF,...)
 ********************************************/

void    IMXFile_Capture(int type)
{ 
  char *quest[16],*file=NULL,**files;
  int e,pos=0,r,max,step;

  /* 1st: What do you want to save ? 
  ** -------------------------------- */
  quest[0] = "Cancel";
  quest[1] = "Print a picture ?";
  quest[2] = "Print <n> pictures in <n> files ?";
  quest[3] = "All pictures in one file?";
  quest[4] = "A window ?";
  quest[5] = "Print a part of screen ?";
  quest[6] = "";

  r=GETV_QCM( TEXT0043,(char **)quest);

  /* => Always in a or more file... */
  if (r!=0) /* Ajout Helgo OHLENBUSCH 06/12/95 */
    file=SAVE_FILE("File ?",NULL, &e);
  if(file==NULL || strlen(file)<=0) {
    PUT_ERROR("Cannot print!\n");
    return;
  }

  /*   Sauvegarde proprement dite...   */
  if (r==1||r==2)
    pos=GET_WND(NULL);
  if (r==2)
    {
      max= GET_INT(TEXT0004, 0, &e);
      step= GET_INT(TEXT0024, 0, &e);
      files=Make_files(file,step,max);
      capture_pictures(pos,max,files,type,r);
    }
  else
    capture_picture(pos,file,type,r);

}

/* --capture_picture() --------------------------------------------------
**
**  pos:    pos of the image to print if demanded
**  r:          What to print
**  type:     which type of image
**  file :     name of file for temporary or file save 
**
*/
int     capture_picture(int pos, char *file, int type, int r)
{
  UINT    wnd_id;
  char cmd[255];
  char str[255];
  int e;

  if (!((type==PCX_FORMAT)||(type==PICT_FORMAT)
	||(type==TIF_FORMAT)||(type==GIF_FORMAT)
	||(type==PPM_FORMAT)||(type==PCXBW_FORMAT)))
    {
      sprintf(str,"Erreur:  Format d'image inconnu pour export\n");
      PUT_MESG(str);
      return(-1);
    }

  /*   Numero d'identification de la window a sauver  */
  switch(r) {
  case 0: /* Cancel */
    break;
  case 1: /* Save a picture */
    wnd_id=(UINT)GDK_WINDOW_XWINDOW(GTK_WIDGET(IMXWidget_PosToWidget(pos)->window));
    sprintf(str,"-window 0x%x",wnd_id);
    break;
  case 2: /* Print <n> pictures in <n> files */
    wnd_id=(UINT)GDK_WINDOW_XWINDOW(GTK_WIDGET(IMXWidget_PosToWidget(pos)->window));
    sprintf(str,"-window 0x%x",wnd_id);
    break;
  case 3: /*  all pictures. */
    wnd_id=(UINT)GDK_WINDOW_XWINDOW(GTK_WIDGET(_wnd_pictures_area->window));
    sprintf(str,"-window 0x%x",wnd_id);
    break;
  case 4: /*  A window. */
    sprintf(str,"-local");
    break;
  case 5: /*  Print a part of screen. */
    sprintf(str," ");
    break;
  }


  /*      Sauvegarde   */

  switch(type) {
  case PCX_FORMAT: /* PCX file*/
    sprintf(cmd,XWPICK " -format pcx %s %s",str,file);
    e=system(cmd);
    break;
  case PCXBW_FORMAT: /* PCX file in gray */
    sprintf(cmd,XWPICK " -gray -format pcx %s %s",str,file);
    e=system(cmd);
    break;
  case GIF_FORMAT: /* GIF file*/
    sprintf(cmd,XWPICK " -format gif %s %s",str,file);
    e=system(cmd);
    break;
  case PICT_FORMAT: /* PICT file*/
    sprintf(cmd,XWPICK " -format pict %s %s",str,file);
    e=system(cmd);
    break;
  case PPM_FORMAT: /* PPM file*/
    sprintf(cmd,XWPICK " -format ppm %s %s",str,file);
    e=system(cmd);
    break;
  case TIF_FORMAT: /* TIF   file*/
    sprintf(cmd,
	    XWPICK " -format ppm %s | /usr/local/bin/"PNMTOTIFF" -none >%s"
	    ,str,file);
    e=system(cmd);
    break;
  }

  printf("comd=%s\n",cmd);
  return(1);
}


/*********************************************
 **
 **   Permet de capturer des images dans un format donne
 **
 *************************************************/
int    capture_pictures(int first_pos, int nb_img, char **pictures, int type, int r)
{ 
  int i,nb_loaded=0;

  if(first_pos+nb_img>(MAX_PICTURES+1))
    return(-1);
  if(first_pos<1 || nb_img<1) return(-1);

  for(i=0;i<nb_img;i++) {
    if(pictures[i]!=NULL && pictures!=NULL) {
      capture_picture(first_pos+i,pictures[i],type,r);
      nb_loaded++;
    }
  }

  return(nb_loaded);
}

/*******************************************
 ** --  IMXFile_Capture_linux ----------------------------
 **
 **  Capture d'une serie d'images.
 **  type : type de sauvegarde  
 ********************************************/

void    IMXFile_Capture_linux(int type)
{ 
  char *quest[16],tmp[256],*file=NULL;
  char finalFname[256];
  int e=0,pos=0,r=0,flg=0;
  /* hide unused windows*/
  HIDE_MESG() ; HIDE_WNDMSG() ; HIDE_WARN() ; HIDE_ERRFAT(); 
  HIDE_ERR() ; HIDE_INFO();    
  /* 1st: What do you want to save ? 
  ** -------------------------------- */
  for (e=0;e<16;e++)
    quest[e] = CALLOC(64,char);
		
  sprintf(quest[0],"Cancel");
  sprintf(quest[1],"Print a picture 2D?");
  sprintf(quest[2],"Print a picture 3D?");
  sprintf(quest[3],"3D image slices in <n> files?");
  sprintf(quest[4],"3D image slices in 1 file?");
  sprintf(quest[5],"all pictures 2D in <n> files?");
  sprintf(quest[6],"all pictures 2D in 1 file?");
  sprintf(quest[7],"%c",'\0');

  r=GETV_QCM( TEXT0043,(char **)quest);
  if (!r) {PUT_ERROR("Cannot print!\n"); 	for (e=0;e<16;e++) FREE(quest[e]); return;}
    
  file=GET_FILE("File ?", &e); 
  if (e) {PUT_ERROR("Cannot print!\n");	for (e=0;e<16;e++) FREE(quest[e]); return;}

  sprintf(finalFname,"%s",file);
  /*   Sauvegarde proprement dite...   */
  sprintf(quest[0],"Cancel");
  sprintf(quest[1],"garder le fond noir?");
  sprintf(quest[2],"supprimer le fond noir?");
  sprintf(quest[3],"%c",'\0');
  switch (r)
    {
    case 1 :{
      pos = GET_WND("");	
      r=GETV_QCM( TEXT0043,(char **)quest);
      if (!r) {PUT_ERROR("Cannot print!\n");for (e=0;e<16;e++) FREE(quest[e]);return;}
      if (r==1) flg = BACKGROUND; else flg = NO_BACKGROUND ; 
      flg = IMX_ScreenShot_2D(pos,finalFname,flg,type);
    }
      break;
    case 2 :{
      pos = GET_WND3D("");
      r=GETV_QCM( TEXT0043,(char **)quest);
      if (!r) {PUT_ERROR("Cannot print!\n");for (e=0;e<16;e++) FREE(quest[e]);return;}
      if (r==1) flg = BACKGROUND; else flg = NO_BACKGROUND ; 
      flg = IMX_ScreenShot_3D(pos,finalFname,flg,type); 
    }
      break;
    case 3 :{
      pos = GET_WND3D("");
      flg = IMX_ScreenShot_n3DinN(pos,finalFname,BACKGROUND,type);
    }
      break;
    case 4 :{
      int img_by_line,step;			   
      pos = GET_WND3D("");
      r=GETV_QCM( TEXT0043,(char **)quest);
      if (!r) {PUT_ERROR("Cannot print!\n");for (e=0;e<16;e++) FREE(quest[e]);return;}
      if (r==1) flg = BACKGROUND; else flg = NO_BACKGROUND ; 
      img_by_line =  GET_INT("nbre d'images par ligne", 0, &e);
      if (e) {PUT_ERROR("Canceled\n");for (e=0;e<16;e++) FREE(quest[e]);return;}
      step =  GET_INT( "pas", 1,  &e);
      if (e) {PUT_ERROR("Canceled\n");for (e=0;e<16;e++) FREE(quest[e]);return;}
      flg = IMX_ScreenShot_n3Din1(pos,finalFname,flg,type,img_by_line,step); 
    }
      break;
    case 5 :flg = IMX_ScreenShot_n2DinN(finalFname,BACKGROUND,type); 
      break;
    case 6 :flg = IMX_ScreenShot_n2Din1(finalFname,BACKGROUND,type); 
      break;

    default:break;
    }
  if (!flg)
    {
      sprintf(tmp,"error saving picture \n");
      PUT_MESG(tmp);
      for (e=0;e<16;e++) FREE(quest[e]);
      return;	
    }

  for (e=0;e<16;e++)
    FREE(quest[e]);

}

#endif /* __GTK_GUI */
