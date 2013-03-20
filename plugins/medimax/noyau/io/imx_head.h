/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef imx_head_h
#define imx_head_h
/**-----------------------------------------------------------------------
***	
***	file:		imx_header.h
***
***	project:	Imagix 1.01 
***			
***
***	description:	Fct about bruker header
***	
***	
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***---------------------------------------------------------------------*/


#ifndef __BRUKER_Header
#define __BRUKER_Header

/* -------------------------- Macros --------------------------	*/

	/* Bruker's pictures contain two parts:		*/
	/*	- header section,       : informations on picture...	*/
        /*               informations stored with 24bits...		*/
	/*	- data section, picture : value for pixel intensity...	*/ 
        /*               informations stored with 32bits...		*/


#define	HEXA	    (1)
#define	ASC6B	    (2)
#define	ASCII8B	    (4)
#define ASCII8BNL   (5)
#define STRING      (6)
#define	ENTIER	    (8)
#define	REEL	   (16)
#define	L_HEXA	   (32)
#define	TAB_PARA   (64)
#define	TAB_PARA_REEL (65)
#define	TAB_PARA_ENTIER (66)
#define ENTIER_IPB (67)
#define REEL_IPB   (68)
#define ASCII_IPB  (69)
#define ENTIER_IPB_2D (80)
#define REEL_IPB_2D   (81)
#define	TAB_PARA_STRING (82)
 
#define	PDAT_OFF	0x1000
#define	MVAR_OFF	0x0800

#define CREATE      0
#define APPEND      1

#endif

#ifndef __PICTURES_DB
#define __PICTURES_DB 

#endif

/* structure for putheader_interfile */
typedef struct {
    int exists;         /* Set to 1 if file exists */
    int maxlines;
    int numlines;
    char **lines;
} HEADER;

#endif /*imx_head_h*/
 
/* AUTOMATICALLY EXTRACTED PROTOTYPES (MB 6/12/2000) */
#ifdef __cplusplus 
extern "C" 
{
#endif /*__cplusplus*/ 
extern char   *decode_pos(unsigned char *header, long unsigned int pos, int type) ;
extern int     codify_pos(unsigned char *header, long unsigned int pos, int type, char *value) ;
#ifdef __cplusplus 
}
#endif /*__cplusplus*/ 
