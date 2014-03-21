/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

//
//  File    : imx_export_file.h
//  Project : MEDIMAX
//
//  Description:
//      ...
//
//  Copyright (c) 1993, ULP-IPB Strasbourg
//  All rights are reserved
//
//  >>  Created on ... ..th, 2001 by ...
//  >>  Last user action on ... ..th, 2001 by ...
//

#ifndef __IMX_EXPORT_FILE_H
#define __IMX_EXPORT_FILE_H

#ifdef __cplusplus 
extern "C" 
{
#endif


  extern void	spI_xprintf(char *s_src1, int i, char *s_src2);
  extern void	spF_xprintf(char *s_src1, float f, char *s_src2);
  extern void	spA_xprintf(char *s_src1, char *s_in, char *s_src2);
  extern int	IMXPutv_OpenLotusFile(char *filename);
  extern void	IMXPutv_CloseLotusFile(void);
  extern int	IMXExport_PutHeader(void);
  extern void	IMXExport_HeaderSet(void);
  
#ifdef __cplusplus 
}
#endif

#endif  // #ifndef __IMX_EXPORT_FILE_H
