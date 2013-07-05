/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

//
//  File    : imx_tools.h
//  Project : GTK-Imagix
//
//  Description:
//      Miscellaneous and useful functions
//
//  Copyright (c) 1993, ULP-IPB Strasbourg
//  All rights are reserved
//
//  >>  Created on May 03rd, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, ... by ...
//

#ifndef __IMX_TOOLS_H
#define __IMX_TOOLS_H

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_FGETS_BUFFER 1024

extern char *imx_trim(char *szExpr);
extern char *imx_spctotab(const char *str, int col);

extern void endianDetection(void);
extern unsigned short shEndianConversion(short unsigned int);
extern unsigned long longEndianConversion (long unsigned int);
extern char **Make_files(char *file, int step, int max);
extern char *incfic(char *name, int step);
extern int incfic_img(char *name, int step, int *numfirst);
extern char *incname_img(char *name);
extern char *incfic_I00X(char *name, int step);
extern void UnAvailable(void);
extern char *imx_GetEnvironment(char *lbl);
extern int imx_IsNumericString(char *str);
extern int imx_decfilename(const char *absol, char *base, char *name);
extern void UnAvailable();
extern char *imx_ltoa(long n);
extern char *ftoa(double n);
extern char *imx_fgets(FILE *fp);
extern void imx_fwrite(void *tblItems, int nbItems, int nSizeItem, FILE *fp);
extern void imx_fread(void *tblItems, int nbItems, int nSizeItem, FILE *fp);
extern int imx_system(char *szCommand, BOOL bWait);
extern int put_file_extension(char * filename, char *ext,char * res);
extern int remove_file_extension(char * filename, char *res);

#ifdef __cplusplus
}
#endif

#endif  // #ifndef __IMX_TOOLS_H
