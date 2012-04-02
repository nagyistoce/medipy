/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

// -----------------------------------------------------------------------
// 	
//  File:       imx_types.h
//  Project:    Imagix 2.01
//  
//  Description: Header file for project imagix (2.01)
//  
//  Copyright (c) 1993-2000, ULP-IPB Strasbourg.
//  All rights are reserved.
//  
//      Created by: Fabien BOULEAU on Jun 21st, 2001
//      Last user action by: ... on ... ..th, ....
//  
// ------------------------------------------------------------------------

#ifndef __TYPES_H
#define __TYPES_H


#ifndef FALSE
#define FALSE 0
#endif /* FALSE */
#ifndef TRUE
#define TRUE 1
#endif /* TRUE */

//  project-dependant #include forbidden here !!! 
//  (standard system libs might be accepted)

#ifndef __cplusplus
typedef int bool;
#endif

typedef unsigned int UINT;
typedef unsigned short USHORT;
typedef unsigned long ULONG;
#ifndef WIN32
typedef bool BOOL;
#endif
typedef unsigned char UCHAR;
typedef long int LONG;
typedef short int SHORT;
typedef void * LPVOID;

typedef UCHAR BYTE;

typedef UINT TDimension;
typedef UINT TMenuId;

typedef struct
{
    ULONG  pixel;
    USHORT red;
    USHORT green;
    USHORT blue;
} IMXColor;

typedef struct
{
    IMXColor *colors;
    UINT     size;
} IMXColormap;

typedef struct
{
    SHORT x;
    SHORT y;
} IMXPoint;

typedef enum
{
    TEXT_NONE,
    TEXT_ALIGN_RIGHT,
    TEXT_ALIGN_LEFT
} TDisplayIndexPos;

/* -- Definition des types d'images dans s_graphic 2D et 3D */
typedef enum
{
    IMAGE_NORMAL
} TImageType;

typedef enum
{
    MASK_BINARY,
    MASK_OPAQUE,
    MASK_SUPERPOSED,
    MASK_TRANSPARENCY
} TMaskType;


typedef struct s_patient
{
    char  name[256];
    char  d_birth[35];
    char  medecin[256];
    char  d_examen[35];
    char  d_stim[35];
    float percent_pem;
} Patient;    
#endif  /* #ifndef __TYPES_H*/
