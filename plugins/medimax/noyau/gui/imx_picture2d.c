/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*!
//  \file    : imx_picture2d.c
//  Project : GTK-Imagix
//
//  Description:
//      ...
//
//  Copyright (c) 1993, ULP-IPB Strasbourg
//  All rights are reserved
//
//  >>  Created on ... ..th, .... by F. BOULEAU
//  >>  Last user action on Oct 29th, 2001 by R. CHEVRIER
//  >> Last user action on May 2007 by Alain BONNET
//
*/
#include <config.h>
#include <errno.h> 
#include <string.h>

#include "noyau/imagix.h" 
#include "noyau/gtkimagix.h"
#include "noyau/imx_2d.h"
#include "noyau/imx_lang.h"
#include "outils/imx_misc.h"
#ifdef __GTK_GUI
# include "noyau/gui/imx_logical_cmap.h"
# include "noyau/gui/imx_picture_common.h"
#endif /* __GTK_GUI */
#include "noyau/gui/imx_picture2d.h"
#include "noyau/gui/imx_picture3d.h"
#ifdef __GTK_GUI
# include "noyau/gui/imx_text.h"
# include "noyau/gui/interface.h"
# include "noyau/imx_zone.h"
# include "noyau/gui/imx_cmap.h"

void (*IMXPict2d_ShowBar)(grphic *pImage, UCHAR *pDest);
#endif /* __GTK_GUI */

#ifndef __DOXYGEN_SHOULD_SKIP_THIS__
static void IMXPict2d_8Zoom1Mask_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pBuff);
static void IMXPict2d_8Zoom1_Dynamic_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pBuff);
static void IMXPict2d_8Zoom1_Dynamic_Sym(ptr_grphic pImage, TMaskType mask_type, UCHAR *pBuff);
static void IMXPict2d_8Zoom1_Fixed_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pBuff);
static void IMXPict2d_8Zoom1_Fixed_Sym(ptr_grphic pImage, TMaskType mask_type, UCHAR *pBuff);

static void IMXPict2d_8Zoom2Mask_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pBuff);
static void IMXPict2d_8Zoom2_Dynamic_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pBuff);
static void IMXPict2d_8Zoom2_Dynamic_Sym(ptr_grphic pImage, TMaskType mask_type, UCHAR *pBuff);
static void IMXPict2d_8Zoom2_Fixed_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pBuff);
static void IMXPict2d_8Zoom2_Fixed_Sym(ptr_grphic pImage, TMaskType mask_type, UCHAR *pBuff);

static void IMXPict2d_8Zoom4Mask_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pBuff);
static void IMXPict2d_8Zoom4_Dynamic_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pBuff);
static void IMXPict2d_8Zoom4_Dynamic_Sym(ptr_grphic pImage, TMaskType mask_type, UCHAR *pBuff);
static void IMXPict2d_8Zoom4_Fixed_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pBuff);
static void IMXPict2d_8Zoom4_Fixed_Sym(ptr_grphic pImage, TMaskType mask_type, UCHAR *pBuff);

static void IMXPict2d_GenericMask_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pBuff);
static void IMXPict2d_GenericTrans_Fixed_Sym(ptr_grphic pImage, TMaskType mask_type, UCHAR *pBuff);
static void IMXPict2d_GenericTrans_Fixed_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pBuff);
static void IMXPict2d_GenericTrans_Dynamic_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pBuff);
static void IMXPict2d_GenericTrans_Dynamic_Sym(ptr_grphic pImage, TMaskType mask_type, UCHAR *pBuff);
static void IMXPict2d_Generic_Dynamic_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pBuff);
static void IMXPict2d_Generic_Dynamic_Sym(ptr_grphic pImage, TMaskType mask_type, UCHAR *pBuff);
static void IMXPict2d_Generic_Fixed_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pBuff);
static void IMXPict2d_Generic_Fixed_Sym(ptr_grphic pImage, TMaskType mask_type, UCHAR *pBuff);

#endif /*__DOXYGEN_SHOULD_SKIP_THIS__*/
#ifdef __GTK_GUI
static void IMXPict2d_ShowBar_Gen(grphic *pImage, UCHAR *pDest);
static void IMXPict2d_ShowBar_8(grphic *pImage, UCHAR *pDest);
    
static TPict2dDisplayFuncItem *tblDisplayFuncs = NULL;

static TPict2dDisplayFuncItem tblDisplayFuncs8bits[] =
  {
    { MASK_OPAQUE,      ZOOM_1,  MAP_DYNAMIC, DSP_SYMMETRICAL, IMXPict2d_8Zoom1_Dynamic_Sym     },
    { MASK_OPAQUE,      ZOOM_1,  MAP_DYNAMIC, DSP_ALL,         IMXPict2d_8Zoom1_Dynamic_Gen     },
    { MASK_OPAQUE,      ZOOM_1,  MAP_FIXED,   DSP_SYMMETRICAL, IMXPict2d_8Zoom1_Fixed_Sym       },
    { MASK_OPAQUE,      ZOOM_1,  MAP_FIXED,   DSP_ALL,         IMXPict2d_8Zoom1_Fixed_Gen       },
    { MASK_SUPERPOSED,  ZOOM_1,  MAP_DYNAMIC, DSP_ALL,         IMXPict2d_8Zoom1_Dynamic_Gen     },
    { MASK_SUPERPOSED,  ZOOM_1,  MAP_FIXED,   DSP_ALL,         IMXPict2d_8Zoom1_Fixed_Gen       },
    { MASK_BINARY,      ZOOM_1,  MAP_ALL,     DSP_ALL,         IMXPict2d_8Zoom1Mask_Gen         },
    { MASK_OPAQUE,      ZOOM_2,  MAP_DYNAMIC, DSP_SYMMETRICAL, IMXPict2d_8Zoom2_Dynamic_Sym     },
    { MASK_OPAQUE,      ZOOM_2,  MAP_DYNAMIC, DSP_ALL,         IMXPict2d_8Zoom2_Dynamic_Gen     },
    { MASK_OPAQUE,      ZOOM_2,  MAP_FIXED,   DSP_SYMMETRICAL, IMXPict2d_8Zoom2_Fixed_Sym       },
    { MASK_OPAQUE,      ZOOM_2,  MAP_FIXED,   DSP_ALL,         IMXPict2d_8Zoom2_Fixed_Gen       },
    { MASK_SUPERPOSED,  ZOOM_2,  MAP_DYNAMIC, DSP_ALL,         IMXPict2d_8Zoom2_Dynamic_Gen     },
    { MASK_SUPERPOSED,  ZOOM_2,  MAP_FIXED,   DSP_ALL,         IMXPict2d_8Zoom2_Fixed_Gen       },
    { MASK_BINARY,      ZOOM_2,  MAP_ALL,     DSP_ALL,         IMXPict2d_8Zoom2Mask_Gen         },
    { MASK_OPAQUE,      ZOOM_4,  MAP_DYNAMIC, DSP_SYMMETRICAL, IMXPict2d_8Zoom4_Dynamic_Sym     },
    { MASK_OPAQUE,      ZOOM_4,  MAP_DYNAMIC, DSP_ALL,         IMXPict2d_8Zoom4_Dynamic_Gen     },
    { MASK_OPAQUE,      ZOOM_4,  MAP_FIXED,   DSP_SYMMETRICAL, IMXPict2d_8Zoom4_Fixed_Sym       },
    { MASK_OPAQUE,      ZOOM_4,  MAP_FIXED,   DSP_ALL,         IMXPict2d_8Zoom4_Fixed_Gen       },
    { MASK_SUPERPOSED,  ZOOM_4,  MAP_DYNAMIC, DSP_ALL,         IMXPict2d_8Zoom4_Dynamic_Gen     },
    { MASK_SUPERPOSED,  ZOOM_4,  MAP_FIXED,   DSP_ALL,         IMXPict2d_8Zoom4_Fixed_Gen       },
    { MASK_BINARY,      ZOOM_4,  MAP_ALL,     DSP_ALL,         IMXPict2d_8Zoom4Mask_Gen         },
    { MASK_OPAQUE,      GENERIC, MAP_DYNAMIC, DSP_SYMMETRICAL, IMXPict2d_Generic_Dynamic_Sym    },
    { MASK_OPAQUE,      GENERIC, MAP_DYNAMIC, DSP_ALL,         IMXPict2d_Generic_Dynamic_Gen    },
    { MASK_OPAQUE,      GENERIC, MAP_FIXED,   DSP_SYMMETRICAL, IMXPict2d_Generic_Fixed_Gen      },
    { MASK_OPAQUE,      GENERIC, MAP_FIXED,   DSP_ALL,         IMXPict2d_Generic_Fixed_Gen      },
    { MASK_SUPERPOSED,  GENERIC, MAP_DYNAMIC, DSP_ALL,         IMXPict2d_Generic_Dynamic_Gen    },
    { MASK_SUPERPOSED,  GENERIC, MAP_FIXED,   DSP_ALL,         IMXPict2d_Generic_Fixed_Gen      },
    { MASK_BINARY,      GENERIC, MAP_ALL,     DSP_ALL,         IMXPict2d_GenericMask_Gen        }
  };

static TPict2dDisplayFuncItem tblDisplayFuncs32bits[] =
  {
    { MASK_OPAQUE,      GENERIC, MAP_DYNAMIC, DSP_SYMMETRICAL, IMXPict2d_Generic_Dynamic_Sym    },
    { MASK_OPAQUE,      GENERIC, MAP_DYNAMIC, DSP_ALL,         IMXPict2d_Generic_Dynamic_Gen    },
    { MASK_OPAQUE,      GENERIC, MAP_FIXED,   DSP_SYMMETRICAL, IMXPict2d_Generic_Fixed_Gen      },
    { MASK_OPAQUE,      GENERIC, MAP_FIXED,   DSP_ALL,         IMXPict2d_Generic_Fixed_Gen      },
    { MASK_SUPERPOSED,  GENERIC, MAP_DYNAMIC, DSP_ALL,         IMXPict2d_Generic_Dynamic_Gen    },
    { MASK_SUPERPOSED,  GENERIC, MAP_FIXED,   DSP_ALL,         IMXPict2d_Generic_Fixed_Gen      },
    { MASK_BINARY,      GENERIC, MAP_ALL,     DSP_ALL,         IMXPict2d_GenericMask_Gen        }
  };

static TPict2dDisplayFuncItem tblDisplayFuncsGen[] =
  {
    { MASK_OPAQUE,      GENERIC, MAP_DYNAMIC, DSP_SYMMETRICAL, IMXPict2d_Generic_Dynamic_Sym    },
    { MASK_OPAQUE,      GENERIC, MAP_DYNAMIC, DSP_ALL,         IMXPict2d_Generic_Dynamic_Gen    },
    { MASK_OPAQUE,      GENERIC, MAP_FIXED,   DSP_SYMMETRICAL, IMXPict2d_Generic_Fixed_Sym      },
    { MASK_OPAQUE,      GENERIC, MAP_FIXED,   DSP_ALL,         IMXPict2d_Generic_Fixed_Gen      },
    { MASK_SUPERPOSED,  GENERIC, MAP_DYNAMIC, DSP_SYMMETRICAL, IMXPict2d_Generic_Dynamic_Sym    },
    { MASK_SUPERPOSED,  GENERIC, MAP_DYNAMIC, DSP_ALL,         IMXPict2d_Generic_Dynamic_Gen    },
    { MASK_SUPERPOSED,  GENERIC, MAP_FIXED,   DSP_SYMMETRICAL, IMXPict2d_Generic_Fixed_Sym      },
    { MASK_SUPERPOSED,  GENERIC, MAP_FIXED,   DSP_ALL,         IMXPict2d_Generic_Fixed_Gen      },
    { MASK_TRANSPARENCY,GENERIC, MAP_DYNAMIC, DSP_SYMMETRICAL, IMXPict2d_GenericTrans_Dynamic_Sym },
    { MASK_TRANSPARENCY,GENERIC, MAP_DYNAMIC, DSP_ALL,         IMXPict2d_GenericTrans_Dynamic_Gen },
    { MASK_TRANSPARENCY,GENERIC, MAP_FIXED,   DSP_SYMMETRICAL, IMXPict2d_GenericTrans_Fixed_Sym },
    { MASK_TRANSPARENCY,GENERIC, MAP_FIXED,   DSP_ALL,         IMXPict2d_GenericTrans_Fixed_Gen },
    { MASK_BINARY,      GENERIC, MAP_ALL,     DSP_ALL,         IMXPict2d_GenericMask_Gen        }
  };

//
//  ---- IMXPict2d_Initialize() -----------------------------------------------
//
//  Description:
//      Initialization of functions vectors array
//
//  >>  Created on Sep 18th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//

void IMXPict2d_Initialize()
{
  IMXPict2d_ShowBar = IMXPict2d_ShowBar_Gen;

  switch(_visual->depth)
    {
    case 8:
      tblDisplayFuncs = tblDisplayFuncs8bits;
      IMXPict2d_ShowBar = IMXPict2d_ShowBar_8;
      break;
    case 16:
    case 24:
      tblDisplayFuncs = tblDisplayFuncsGen;
      break;
    case 32:
      tblDisplayFuncs = tblDisplayFuncs32bits;
      break;
    default:
      PUT_ERRFAT( "IMXPict2d_Initialize(): wrong display depth\n");
    };
}

//
//  ---- Common code parts ----------------------------------------------------
//
//  Description:
//      ...
//
//  >>  Created on Sep 24th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//

#define IMXPICT2D_STOREPIXEL_GENERIC                                        \
{                                                                           \
    if(mask_type == MASK_OPAQUE || !bBlack )                                 \
    {                                                                       \
        for(k = 0 ; k < nZoom ; k++)                                        \
        {                                                                   \
            for(l = 0 ; l < nZoom ; l++)                                    \
            {                                                               \
                for(m = 0 ; m < _nBytesPerRGB ; m++)                        \
                {                                                           \
                    n = (j * _pic_size_width * _nBytesPerPixel			    \
                            + i * _nBytesPerPixel) * nZoom         			\
                        + k * _pic_size_width * _nBytesPerPixel             \
                        + l * _nBytesPerPixel + m;                          \
																			\
					pDest[n] = (pix >> (sizeof(UCHAR) * 8 * m)) & 255;      \
																			\
              }                                                             \
            }                                                               \
        }                                                                   \
    }                                                                       \
}

#define IMXPICT2D_STOREPIXEL_TRANSPARENCY                                        \
{                                                                           \
    if(mask_type == MASK_TRANSPARENCY || !bBlack )                                 \
    {                                                                       \
        for(k = 0 ; k < nZoom ; k++)                                        \
        {                                                                   \
            for(l = 0 ; l < nZoom ; l++)                                    \
            {                                                               \
                for(m = 0 ; m < _nBytesPerRGB ; m++)                        \
                {                                                           \
                    n = (j * _pic_size_width * _nBytesPerPixel              \
                            + i * _nBytesPerPixel) * nZoom                  \
                        + k * _pic_size_width * _nBytesPerPixel             \
                        + l * _nBytesPerPixel + m;                          \
                                                                            \
                    if (!pDest[n])                                          \
                    {                                                       \
                      pDest[n] = (((pix >> (sizeof(UCHAR) * 8 * m))& 255) *pImage->alpha_channel)/100;;    \
                    }                                                       \
                    else                                                    \
                    {                                                       \
                      unsigned char  color_ancienne=pDest[n];               \
                      unsigned char  color_nouvelle=0;                      \
                      if (pix)                                              \
                      {                                                     \
                        color_ancienne = (pDest[n]*(100-pImage->alpha_channel))/100;\
                        color_nouvelle = (((pix >> (sizeof(UCHAR) * 8 * m))& 255) *pImage->alpha_channel)/100;\
                      }                                                     \
                                                                            \
                      pDest[n] = ( color_nouvelle + color_ancienne);        \
                    }                                                       \
                                                                            \
              }                                                             \
            }                                                               \
        }                                                                   \
    } 		                                                                \
}

#define IMXPICT2D_INITIALIZEDATA                                            \
{                                                                           \
    pix     = 0l;                                                           \
    bBlack  = FALSE;                                                        \
    cmap    = pImage->cmapinfo.cmap;                                        \
    nZoom   = pImage->zoom;                                                 \
    start   = 0.;                                                           \
    end     = cmap->nbColors - 1;                                           \
                                                                            \
    if(pImage->cmapinfo.dspStyle == DSP_POSITIVE)                           \
    {                                                                       \
      if (pImage->cutoff_min <=0.)                                          \
        pImage->cutoff_min = 0.;                                            \
    }                                                                       \
		 											\
    if(pImage->cmapinfo.dspStyle == DSP_NEGATIVE)                           \
    {                                                                       \
      if (pImage->cutoff_max >=0.)                                          \
        pImage->cutoff_max = 0.;                                            \
    }                                                                       \
																			\
    imin    = pImage->cutoff_min;                                           \
    imax    = pImage->cutoff_max;                                           \
                                                                            \
    /*if(imin == imax)                                                      \
      return;*/                                                           \
                                                                            \
    iImgStartX = pImage->zoom_x;                                            \
    iImgStartY = pImage->zoom_y;                                            \
																			\
    iImgEndX = MINI(_pic_size_width / nZoom, iImgStartX + pImage->width);   \
    iImgEndY = MINI(_pic_size_height / nZoom, iImgStartY + pImage->height); \
                                                                            \
    OverBlancking = MAXMRI;                                      \
    UnderBlancking= -MAXMRI;                                      \
}

#define IMXPICT2D_DECLAREVARS                                               \
    static int i, j, k, l, m, n;                                            \
    static LONG start, end;                                                 \
    static LONG imin, imax;                                                 \
    static ptrColormap cmap;                                                \
    static int nZoom;                                                       \
    static int offset_x,offset_y;                                           \
    static LONG pix;                                                        \
    static double pt;                                                       \
    static BOOL trouve;                                                     \
    static UCHAR *pDest2, *pDest3, *pDest4;                                 \
    static UINT iImgStartX, iImgStartY;                                     \
    static UINT iImgEndX, iImgEndY;                                         \
    static double valrank;                                                  \
    static BOOL bBlack;                                                     \
    static int OverBlancking;                                               \
    static int UnderBlancking;

#define IMXPICT2D_CALCPIXEL_DYNAMIC                                         \
  if (!bBlack)                                                              \
  {                                                                         \
    if( (imax - imin) + start != 0 )                                        \
      pix = (double) ( pt - imin) * (end - start)                           \
          / (imax - imin) + start;                                          \
    else                                                                    \
      pix = 0l;                                                             \
    bBlack = (cmap->tblColors[pix].red == 0                                 \
            && cmap->tblColors[pix].green == 0                              \
            && cmap->tblColors[pix].blue == 0);                             \
   	pix = cmap->tblColors[pix].pixel;                                       \
  }

#define IMXPICT2D_CALC_OVERBLANCKING                                        \
{                                                                           \
    switch (pImage->cmapinfo.dspStyle)                                      \
    {                                                                       \
        case DSP_POSITIVE : if (pt<0) { UnderBlancking = 0;}                \
        break;                                                              \
        case DSP_NEGATIVE : if (pt>0) { OverBlancking = 0;}                 \
        break;                                                              \
        default : ;/*do nothing*/                                           \
    }                                                                       \
                                                                            \
    switch (pImage->cmapinfo.dspOverBlancking)                              \
    {                                                                       \
        case OVERBLANCKING_ALL :                                            \
            OverBlancking = imax;                                           \
            UnderBlancking = imin;                                          \
        break;                                                              \
        case OVERBLANCKING_NONE : /* do nothing*/                           \
        break;                                                              \
        default : /* do nothing*/                                           \
        break;                                                              \
    }                                                                       \
    if ( (  pt > OverBlancking ) || (pt < UnderBlancking ) )                \
    {                                                                       \
      bBlack = TRUE; pix = 0l;                                              \
    }                                                                       \
    else bBlack = FALSE;                                                    \
}                                                                           



    
#define IMXPICT2D_TEST_THRESHOLD                                            \
    if(_threshold_draw)                                                     \
    {                                                                       \
        if(pt >= _maxmin_min && pt <= _maxmin_max)                          \
		{				                            \
			if (gdk_visual_get_system()->type == GDK_VISUAL_TRUE_COLOR) 	\
            {										\
				pix = _threshold_color24b;                              \
			}								\
			else															\
			{ 																\
            	pix = _threshold_color;										\
			}\
		}                                     							\
    }

#define IMXPICT2D_CALCPIXEL_FIXED                                       \
  if (!bBlack)                                                          \
  {                                                                     \
    for(trouve = k = 0 ; !trouve && k < cmap->nbStage ; k++)            \
    {                                                                   \
        if((pt >= cmap->tblStageFrom[k]                                 \
                && pt < cmap->tblStageTo[k])                            \
            ||  (cmap->tblStageFrom[k] == cmap->tblStageTo[k]           \
                && pt == cmap->tblStageFrom[k]))                        \
        {                                                               \
            trouve = TRUE;                                              \
            pix    = k;                                                 \
            break ;		                                                  \
        }                                                               \
    }                                                                   \
    if (trouve)                                                         \
    {                                                                   \
    bBlack = (cmap->tblColors[pix].red == 0                             \
            && cmap->tblColors[pix].green == 0                          \
            && cmap->tblColors[pix].blue == 0);                         \
                                                                        \
    pix = cmap->tblColors[pix].pixel;                                   \
    }                                                                   \
    else pix = 0;                                                       \
  }
//
//  ---- IMXPict2d_ShowIndexes() ----------------------------------------------------
//
//  Description:
//      ...
//
//  Parameters:
//      - 
//
//  Return value: always 0
//
//  >>  Created on Oct 05th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//
void IMXPict2d_ShowIndexes(int pos, ptrColormap cmap, TDisplayIndexPos nIndexPos)
{
  int i, h;
  IMXPoint ptfrom;
  double valfrom;
  char str[255];

  if(cmap->noType == MAP_FIXED)
    {
      h = _pic_size_height;

      for(i = 0 ; i < cmap->nbColors ; i++)
        {
	  switch(nIndexPos)
            {
	    case TEXT_ALIGN_RIGHT:
	      ptfrom.x = _pic_size_width - 30;
	      break;
	    case TEXT_ALIGN_LEFT:
	      ptfrom.x = BAR_LENGTH + 6;
	      break;
	    default : // do nothin
	      break;
            }
    
	  for(i = 0 ; i < cmap->nbColors ; i++)
            {
	      valfrom     = cmap->tblStageFrom[i];
	      ptfrom.y    = (cmap->nbColors - i) * ((8 + h) / cmap->nbColors);
                
	      if(ptfrom.y > h - 2)
		ptfrom.y = h - 2;
                
	      sprintf(str,"%s","");
                
	      switch((int) valfrom)
                {
		case -LONGMAX:
		case +LONGMAX:
		  break;
		default:
		  sprintf(str, "%3.1f", valfrom);
		  break;
                }
                
	      IMXText_DrawText(IMXWidget_PosToWidget(pos)
			       , _imx_img[pos - 1].gc,
			       ptfrom, str, strlen(str));
            }
        }
    }
}

//
//  ---- IMXPict2d_ShowBar() --------------------------------------------------------
//
//  Description:
//      ...
//
//  Parameters:
//      - 
//
//  Return value: always 0
//
//  >>  Created on Oct 05th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//

static void IMXPict2d_ShowBar_8(grphic *pImage, UCHAR *pDest)
{
  int i, j;
  static GdkColor tbltmp[IMXCOLOR__MAX_LOGICAL_COLORS];
  static GdkColor tbltmp2[IMXCOLOR__MAX_LOGICAL_COLORS];
  static int         noCMap;
  static ptrColormap cmap; 

  noCMap = pImage->cmapinfo.noColormap;
  cmap  = pImage->cmapinfo.cmap;

  IMXColor_spGetColormap(tbltmp, tbltmp2,cmap->nbColors, cmap->iOrigin, cmap->iCoeff, cmap->bOverblk);

  for(i = 0 ; i < cmap->nbColors ; i++)
    {
      tbltmp2[i].pixel = IMXCMap_RGBToPixel(
					    i * (cmap->iEnd - cmap->iStart + 1) / cmap->nbColors + cmap->iStart, 
					    tbltmp2[i].red, 
					    tbltmp2[i].green, 
					    tbltmp2[i].blue);
    }

  for(i = 0 ; i < BAR_LENGTH ; i++)
    {
      for(j = 0 ; j < _pic_size_height ; j++)
        {
	  pDest[(_pic_size_height - j - 1) * _pic_size_width + i] 
	    = tbltmp2[cmap->nbColors * j / _pic_size_height].pixel;
        }
    }
}

//
//  ---- IMXPict2d_ShowBar_Gen() --------------------------------------
//
//  Description:
//      ...      Affichage du la barre de couleur
//              Il faut penser a decouper cette fonction en plusieurs sous fonctions
//              suivant les differents cas (SYMETRIC... POSITIV... etc...)
//
//
static void IMXPict2d_ShowBar_Gen(grphic *pImage, UCHAR *pDest)
{
  int i, j, k, n;
  static int         noCMap;
  static ptrColormap cmap;
  static int index,max,min;
  static int index_min;
  static int index_max;
  static int OverBlancking;
  static int UnderBlancking;
  static BOOL bBlack;
  static LONG pix;
  static TYPEMRI pt;
  static LONG imin, imax;
  static TYPEMRI minpixel,maxpixel;
  static TYPEMRI interval;

  noCMap = pImage->cmapinfo.noColormap;
  cmap   = pImage->cmapinfo.cmap;
    
  //    IMXColor_spGetColormap(tbltmp, tbltmp2, 
  //    cmap->nbColors, cmap->iOrigin, cmap->iCoeff, cmap->bOverblk);

  minpixel =pImage->min_pixel;
  maxpixel =pImage->max_pixel;
  interval = maxpixel - minpixel;

  if (pImage->max_pixel -  pImage->min_pixel)
    { // max et min de la barre (en coord pixel)
      max = (pImage->cutoff_max - minpixel) * _pic_size_height / interval;
      min = (pImage->cutoff_min - minpixel) * _pic_size_height / interval;
    }
  else
    {
      max = _pic_size_height;
      min = 0;
    }

  // pour palette par pallier
  if (cmap->noType == MAP_FIXED)
    {
      index_min = IMXColor_GetStageFromValue(pImage->cutoff_min*pImage->rcoeff,cmap);
      index_max = IMXColor_GetStageFromValue(pImage->cutoff_max*pImage->rcoeff,cmap);
    }

  // pour overblancking
  bBlack = 0;
  OverBlancking = MAXMRI;
  UnderBlancking= -MAXMRI;

  for(i = 0 ; i < BAR_LENGTH ; i++)
    {
      for(j = 0 ; j < _pic_size_height ; j++)
        {
	  n = (((_pic_size_height - j - 1) * _pic_size_width + i) * _nBytesPerPixel);

	  switch (cmap->noType)
            {
	    case  MAP_DYNAMIC :
	      imin    = minpixel + min*interval/_pic_size_height;
	      imax    = minpixel + max*interval/_pic_size_height;
	      pt =minpixel +  j*interval/_pic_size_height;
	      IMXPICT2D_CALC_OVERBLANCKING
		if (!bBlack)
		  {
		    if (j>max)
		      index = cmap->nbColors;
		    else if (j<=min)
		      index = 0;
		    else if ((j<max) && (j>min))
		      {
			index = cmap->nbColors * (j-min)  /(max-min);
		      }      
		  }
		else index = 0;
              break;
	    case MAP_FIXED :
	      index = cmap->nbColors * j /_pic_size_height;
	      imin    = pImage->cutoff_min;
	      imax    = pImage->cutoff_max;
	      pt = (TYPEMRI)(cmap->tblStageFrom[index]/pImage->rcoeff);
	      switch (pImage->cmapinfo.dspStyle)
		{
		case DSP_SYMMETRICAL :
		  if ((TYPEMRI)(cmap->tblStageTo[index]/pImage->rcoeff)<imin
		      || (TYPEMRI)(cmap->tblStageTo[index]/pImage->rcoeff)>imax)
		    bBlack = FALSE;
		  else
		    bBlack = TRUE;
		  break;
		default :
		  IMXPICT2D_CALC_OVERBLANCKING;
		  if (!bBlack)
		    {
		      if (index>index_max) index = index_max;
		      if (index<index_min) index = index_min;
		    }
		  break;
		}
	      ; break;
	    default : break ; //do nothing
            }
 
	  if (index >= cmap->nbColors)
	    index = cmap->nbColors - 1;
	  if (index<0)
	    index = 0; 
	  if (!bBlack) //computed in IMXPICT2D_CALC_OVERBLANCKING
	    pix = cmap->tblColors[index].pixel;
	  else
	    pix = 0l;

	  for(k = 0 ; k < _nBytesPerRGB ; k++)
	    pDest[n+k] = (pix >> (sizeof(UCHAR) *8  * k));
        }
    }
}

//
// >>>> GENERIC PART <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//

//
//  ---- IMXPict2d_GenericMask_Gen() --------------------------------------
//
//  Description:
//      ...      Affichage du mask en mode binaire
//
//  >>  Created on Sep 18th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//

static void IMXPict2d_GenericMask_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pDest)
{
  static int i, j, k, l, m, n;                                            
  static LONG start, end;                                                 
  static LONG imin, imax;                                                 
  static ptrColormap cmap;                                                
  static int nZoom;							    
  static LONG pix;                                                        
  static double pt;                                                       
  static UINT iImgStartX, iImgStartY;                                     
  static UINT iImgEndX, iImgEndY;                                         
  static BOOL bBlack;                                                     
  static int OverBlancking;
  static int UnderBlancking;

  IMXPICT2D_INITIALIZEDATA
    for(i = 0/* iImgStartX */ ; i < iImgEndX ; i++)
      {
        for(j = 0/*iImgStartY */ ; j < iImgEndY ; j++)
	  {
            pt = pImage->mri[i+iImgStartX/nZoom][j+iImgStartY/nZoom];

            if (gdk_visual_get_system()->type == GDK_VISUAL_TRUE_COLOR)
	      pix = _mask_color24b;
            else
	      pix = _mask_color;
            if (pt)
	      bBlack = FALSE;
            else
	      bBlack = TRUE;
            IMXPICT2D_STOREPIXEL_GENERIC
	      }
      }
}

//
//  ---- IMXPict2d_Generic_Dynamic_Gen() --------------------------------------
//
//  Description:
//      ...      Affichage image palette dynamique en mode POSITIVE NEGATIVE et DISSYMETRIQUE
//
//  >>  Created on Sep 18th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//
static void IMXPict2d_Generic_Dynamic_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pDest)
{
  static int i, j, k, l, m, n;
  static LONG start, end;
  static LONG imin, imax;
  static ptrColormap cmap;
  static int nZoom;
  static LONG pix;
  static double pt;
  static UINT iImgStartX, iImgStartY;
  static UINT iImgEndX, iImgEndY;
  static BOOL bBlack;
  static int OverBlancking;
  static int UnderBlancking;
  IMXPICT2D_INITIALIZEDATA

    for(i = 0/* iImgStartX */ ; i < iImgEndX ; i++)
    {
      for(j = 0/*iImgStartY */ ; j < iImgEndY ; j++)
      {
            pt = pImage->mri[i+iImgStartX/nZoom][j+iImgStartY/nZoom];

            IMXPICT2D_CALC_OVERBLANCKING

            if(pt < imin) pt = imin;
            if(pt > imax) pt = imax;

            IMXPICT2D_CALCPIXEL_DYNAMIC
            IMXPICT2D_TEST_THRESHOLD
            IMXPICT2D_STOREPIXEL_GENERIC
      }
    }
}

//
//  ---- IMXPict2d_Generic_Dynamic_Sym() --------------------------------------
//
//  Description:
//      ...        palette dynamique symetrique
//
//  >>  Created on Sep 18th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//

static void IMXPict2d_Generic_Dynamic_Sym(ptr_grphic pImage, TMaskType mask_type, UCHAR *pDest)
{
  static int i, j, k, l, m, n;
  static LONG start, end;
  static LONG imin, imax;
  static ptrColormap cmap;
  static int nZoom;
  static LONG pix;
  static double pt;
  static UINT iImgStartX, iImgStartY;
  static UINT iImgEndX, iImgEndY;
  static BOOL bBlack;
  static int OverBlancking;
  static int UnderBlancking;
  IMXPICT2D_INITIALIZEDATA

    for(i =0/* iImgStartX */ ; i < iImgEndX ; i++)
      {
        for(j = 0/*iImgStartY */ ; j < iImgEndY ; j++)
	  {
            pt = pImage->mri[i+iImgStartX/nZoom][j+iImgStartY/nZoom];


            IMXPICT2D_CALC_OVERBLANCKING
	      if(pt < imin) pt = imin;
            if(pt > imax) pt = imax;

            if(pt <= 0.)
	      {
                end     = (end - start) / 2.;
                imax    = 0.;
	      }
            else if(pt > 0.)
	      {
                start   = (end - start) / 2.;
                imin    = 0.;
	      }
            IMXPICT2D_CALCPIXEL_DYNAMIC

	      IMXPICT2D_TEST_THRESHOLD
	      IMXPICT2D_STOREPIXEL_GENERIC

	      //			il faut remettre a jour les valeurs de start, end, imin et imax qui ont ete changees

	      start   = 0.;
	    end     = cmap->nbColors - 1;
	    imin    = pImage->cutoff_min;
	    imax    = pImage->cutoff_max;

	  }
      }
}

//
//  ---- IMXPict2d_Generic_Fixed_Gen() ----------------------------------------
//
//  Description:
//                   Palette par pallier POSITIVE NEGATIVE DISSYMETRIQUE SYMETRIQUE
//      ...
//
//  >>  Created on Sep 24th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//
static void IMXPict2d_Generic_Fixed_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pDest)
{
  static int i, j, k, l, m, n;
  static LONG start, end;
  static LONG imin, imax;
  static ptrColormap cmap;
  static int nZoom;
  static LONG pix;
  static double pt;
  static BOOL trouve;
  static UINT iImgStartX, iImgStartY;
  static UINT iImgEndX, iImgEndY;
  static BOOL bBlack;
  static int OverBlancking;
  static int UnderBlancking;
  IMXPICT2D_INITIALIZEDATA

    for(i =0; i < iImgEndX ; i++)
      {
        for(j = 0; j <iImgEndY ; j++)
	  {
            pt = pImage->mri[i+iImgStartX/nZoom][j+iImgStartY/nZoom];

            IMXPICT2D_CALC_OVERBLANCKING

	      if(pt < imin) pt = imin;
            if(pt > imax) pt = imax;
            pt   *= pImage->rcoeff;

            IMXPICT2D_CALCPIXEL_FIXED
            
	      IMXPICT2D_TEST_THRESHOLD
	      IMXPICT2D_STOREPIXEL_GENERIC
	      }
      }
}

//
//  ---- IMXPict2d_Generic_Fixed_Sym() ----------------------------------------
//
//  Description:
//      ...         palettes par pallier en mode symetrique ???
//                  en fait cette fonction a ete recyclee, maintenant elle sert
//                  a n'afficher que les valeurs qui sont comprises entre [minpixel ; cutoff_min]
//                  et [cutoff_max ; max_pixel]
//
//  >>  Created on Sep 24th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//

static void IMXPict2d_Generic_Fixed_Sym(ptr_grphic pImage, TMaskType mask_type, UCHAR *pDest)
{
  static int i, j, k, l, m, n;
  static LONG start, end;
  static LONG imin, imax;
  static ptrColormap cmap;
  static int nZoom;
  static LONG pix;
  static double pt;
  static BOOL trouve;
  static UINT iImgStartX, iImgStartY;
  static UINT iImgEndX, iImgEndY;
  static BOOL bBlack;
  static int OverBlancking;
  static int UnderBlancking;
  IMXPICT2D_INITIALIZEDATA

    for(i = 0 ; i < iImgEndX ; i++)
      {
        for(j = 0 ; j < iImgEndY ; j++)
	  {
            pt = pImage->mri[i+iImgStartX/nZoom][j+iImgStartY/nZoom];
                                                                            
            OverBlancking = imax;                                           
            UnderBlancking = imin;                                          
            if ( (  pt > imin ) && (pt < imax ) )                
	      {                                                                       
		bBlack = TRUE; pix = 0l;                                              
	      }                                                                       
            else bBlack = FALSE;                                                    

            pt   *= pImage->rcoeff;

            IMXPICT2D_CALCPIXEL_FIXED

	      IMXPICT2D_TEST_THRESHOLD
	      IMXPICT2D_STOREPIXEL_GENERIC
	      }
      }
}



//
//  ---- IMXPict2d_GenericTrans_Fixed_Gen() ----------------------------------------
//
//  Description:
//      ...      Transparence des palettes par pallier
//
//  >>  Created on Sep 24th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//

static void IMXPict2d_GenericTrans_Fixed_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pDest)
{
  static int i, j, k, l, m, n;
  static LONG start, end;
  static LONG imin, imax;
  static ptrColormap cmap;
  static int nZoom;
  static LONG pix;
  static double pt;
  static BOOL trouve;
  static UINT iImgStartX, iImgStartY;
  static UINT iImgEndX, iImgEndY;
  static BOOL bBlack;
  static int OverBlancking;
  static int UnderBlancking;
  IMXPICT2D_INITIALIZEDATA

    for(i =0/* iImgStartX */; i < iImgEndX ; i++)
      {
        for(j = 0/*iImgStartY */; j <iImgEndY ; j++)
	  {
            pt = pImage->mri[i+iImgStartX/nZoom][j+iImgStartY/nZoom];

            IMXPICT2D_CALC_OVERBLANCKING

	      if(pt < imin) pt = imin;
            if(pt > imax) pt = imax;

            pt = pt * pImage->rcoeff;

            IMXPICT2D_CALCPIXEL_FIXED

	      IMXPICT2D_TEST_THRESHOLD
	      IMXPICT2D_STOREPIXEL_TRANSPARENCY
	      }
      }
}

//
//  ---- IMXPict2d_GenericTrans_Fixed_Sym() ----------------------------------------
//
//  Description:
//      ...      Transparence des palettes par pallier
//
//  >>  Created on Sep 24th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//

static void IMXPict2d_GenericTrans_Fixed_Sym(ptr_grphic pImage, TMaskType mask_type, UCHAR *pDest)
{
  static int i, j, k, l, m, n;
  static LONG start, end;
  static LONG imin, imax;
  static ptrColormap cmap;
  static int nZoom;
  static LONG pix;
  static double pt;
  static BOOL trouve;
  static UINT iImgStartX, iImgStartY;
  static UINT iImgEndX, iImgEndY;
  static BOOL bBlack;
  static int OverBlancking;
  static int UnderBlancking;
  IMXPICT2D_INITIALIZEDATA

    for(i =0/* iImgStartX */; i < iImgEndX ; i++)
      {
        for(j = 0/*iImgStartY */; j <iImgEndY ; j++)
	  {
            pt = pImage->mri[i+iImgStartX/nZoom][j+iImgStartY/nZoom];

            OverBlancking = imax;
            UnderBlancking = imin;
            if ( (  pt > imin ) && (pt < imax ) )
	      {
		bBlack = TRUE; pix = 0l;
	      }
            else bBlack = FALSE;

            pt   *= pImage->rcoeff;

            IMXPICT2D_CALCPIXEL_FIXED

	      IMXPICT2D_TEST_THRESHOLD
	      IMXPICT2D_STOREPIXEL_TRANSPARENCY
	      }
      }
}

//
//  ---- IMXPict2d_GenericTrans_Dynamic_Sym() ----------------------------------------
//
//  Description:
//      ...      Transparence des palettes dynamiques SYMETRIQUE
//
//  >>  Created on Sep 24th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//
static void IMXPict2d_GenericTrans_Dynamic_Sym(ptr_grphic pImage, TMaskType mask_type, UCHAR *pDest)
{
  static int i, j, k, l, m, n;
  static LONG start, end;
  static LONG imin, imax;
  static ptrColormap cmap;
  static int nZoom;
  static LONG pix;
  static double pt;
  static UINT iImgStartX, iImgStartY;
  static UINT iImgEndX, iImgEndY;
  static BOOL bBlack;
  static int OverBlancking;
  static int UnderBlancking;
  IMXPICT2D_INITIALIZEDATA

    for(i =0/* iImgStartX */ ; i < iImgEndX ; i++)
      {
        for(j = 0/*iImgStartY */ ; j < iImgEndY ; j++)
	  {
            pt = pImage->mri[i+iImgStartX/nZoom][j+iImgStartY/nZoom];


            IMXPICT2D_CALC_OVERBLANCKING

	      if(pt < imin) pt = imin;
            if(pt > imax) pt = imax;

            if(pt <= 0.)
	      {
                end     = (end - start) / 2.;
                imax    = 0.;
	      }
            else if(pt > 0.)
	      {
                start   = (end - start) / 2.;
                imin    = 0.;
	      }
            IMXPICT2D_CALCPIXEL_DYNAMIC

	      IMXPICT2D_TEST_THRESHOLD
	      IMXPICT2D_STOREPIXEL_TRANSPARENCY

	      //			il faut remettre a jour les valeurs de start, end, imin et imax qui ont ete changees

	      start   = 0.;
	    end     = cmap->nbColors - 1;
	    imin    = pImage->cutoff_min;
	    imax    = pImage->cutoff_max;

	  }
      }
}

//
//  ---- IMXPict2d_GenericTrans_Dynamic_Gen() ----------------------------------------
//
//  Description:
//      ...      Transparence des palettes dynamiques POSITIVE NEGATIVE DISSYMETRIQUE
//
//  >>  Created on Sep 24th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//
static void IMXPict2d_GenericTrans_Dynamic_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pDest)
{
  static int i, j, k, l, m, n;
  static LONG start, end;
  static LONG imin, imax;
  static ptrColormap cmap;
  static int nZoom;
  static LONG pix;
  static double pt;
  static UINT iImgStartX, iImgStartY;
  static UINT iImgEndX, iImgEndY;
  static BOOL bBlack;
  static int OverBlancking;
  static int UnderBlancking;
  IMXPICT2D_INITIALIZEDATA

    for(i = 0/* iImgStartX */ ; i < iImgEndX ; i++)
      {
        for(j = 0/*iImgStartY */ ; j < iImgEndY ; j++)
	  {
            pt = pImage->mri[i+iImgStartX/nZoom][j+iImgStartY/nZoom];

            IMXPICT2D_CALC_OVERBLANCKING

	      if(pt < imin) pt = imin;
            if(pt > imax) pt = imax;
            IMXPICT2D_CALCPIXEL_DYNAMIC

	      IMXPICT2D_TEST_THRESHOLD
	      IMXPICT2D_STOREPIXEL_TRANSPARENCY
	      }
      }
}




//
// >>>> 8 BITS ZOOM 1 PART <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//

//
//  ---- IMXPict2d_8Zoom1_Dynamic_Gen() ---------------------------------------
//
//  Description:
//      ...
//
//  >>  Created on Sep 18th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//

static void IMXPict2d_8Zoom1Mask_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pDest)
{
  static int i, j, k, l, m, n;
  static LONG start, end;
  static LONG imin, imax;
  static ptrColormap cmap;
  static int nZoom;
  static LONG pix;
  static double pt;
  static UCHAR *pDest2;
  static UINT iImgStartX, iImgStartY;
  static UINT iImgEndX, iImgEndY;
  static BOOL bBlack;
  static int OverBlancking;
  static int UnderBlancking;
  IMXPICT2D_INITIALIZEDATA

    pDest2  = pDest + _pic_size_width * _nBytesPerPixel;

  for(j = 0/*iImgStartY */ ; j < iImgEndY ; j++)
    {
      for(i = 0/* iImgStartX */ ; i < iImgEndX ; i++)
        {
	  pt = pImage->mri[i+iImgStartX/nZoom][j+iImgStartY/nZoom];

	  if(pt < imin) pt = imin;
	  if(pt > imax) pt = imax;

	  IMXPICT2D_CALCPIXEL_DYNAMIC
            pix = _mask_color;
	  IMXPICT2D_STOREPIXEL_GENERIC
	    }
    }
}

static void IMXPict2d_8Zoom1_Dynamic_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pDest)
{
  static int i, j, k, l, m, n;
  static LONG start, end;
  static LONG imin, imax;
  static ptrColormap cmap;
  static int nZoom;
  static LONG pix;
  static double pt;
  static UCHAR *pDest2;
  static UINT iImgStartX, iImgStartY;
  static UINT iImgEndX, iImgEndY;
  static BOOL bBlack;
  static BOOL OverBlancking;
  static BOOL UnderBlancking;
  IMXPICT2D_INITIALIZEDATA

    pDest2  = pDest + _pic_size_width * _nBytesPerPixel;

  for(j = 0/*iImgStartY */ ; j < iImgEndY ; j++)
    {
      for(i = 0/* iImgStartX */ ; i < iImgEndX ; i++)
        {
	  pt = pImage->mri[i+iImgStartX/nZoom][j+iImgStartY/nZoom];

	  if(pt < imin) pt = imin;
	  if(pt > imax) pt = imax;

	  IMXPICT2D_CALCPIXEL_DYNAMIC
	    IMXPICT2D_TEST_THRESHOLD
            IMXPICT2D_STOREPIXEL_GENERIC
	    }
    }
}

//
//  ---- IMXPict2d_8Zoom1_Dynamic_Sym() ---------------------------------------
//
//  Description:
//      ...
//
//  >>  Created on Sep 18th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//

static void IMXPict2d_8Zoom1_Dynamic_Sym(ptr_grphic pImage, TMaskType mask_type, UCHAR *pDest)
{
  static int i, j, k, l, m, n;
  static LONG start, end;
  static LONG imin, imax;
  static ptrColormap cmap;
  static int nZoom;
  static LONG pix;
  static double pt;
  static UCHAR *pDest2;
  static UINT iImgStartX, iImgStartY;
  static UINT iImgEndX, iImgEndY;
  static BOOL bBlack;
  static int OverBlancking;
  static int UnderBlancking;
  IMXPICT2D_INITIALIZEDATA

    pDest2  = pDest + _pic_size_width * _nBytesPerPixel;

  for(j = 0/*iImgStartY */ ; j < iImgEndY ; j++)
    {
      for(i = 0/* iImgStartX */ ; i < iImgEndX ; i++)
        {
	  pt = pImage->mri[i+iImgStartX/nZoom][j+iImgStartY/nZoom];

	  if(pt < imin) pt = imin;
	  if(pt > imax) pt = imax;

	  if(pt <= 0.)
            {
	      end     = (end - start) / 2.;
	      imax    = 0.;
            }
	  else if(pt > 0.)
            {
	      start   = (end - start) / 2.;
	      imin    = 0.;
            }

	  IMXPICT2D_CALCPIXEL_DYNAMIC
            IMXPICT2D_TEST_THRESHOLD
            IMXPICT2D_STOREPIXEL_GENERIC
	    }
    }
}

//
//  ---- IMXPict2d_8Zoom1_Fixed_Gen() -----------------------------------------
//
//  Description:
//      ...
//
//  >>  Created on Sep 24th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//

static void IMXPict2d_8Zoom1_Fixed_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pDest)
{
  static int i, j, k, l, m, n;
  static LONG start, end;
  static LONG imin, imax;
  static ptrColormap cmap;
  static int nZoom;
  static LONG pix;
  static double pt;
  static BOOL trouve;
  static UCHAR *pDest2;
  static UINT iImgStartX, iImgStartY;
  static UINT iImgEndX, iImgEndY;
  static BOOL bBlack;
  static int OverBlancking;
  static int UnderBlancking;
  IMXPICT2D_INITIALIZEDATA

    pDest2  = pDest + _pic_size_width * _nBytesPerPixel;

  for(j = 0/*iImgStartY */ ; j < iImgEndY ; j++)
    {
      for(i = 0/* iImgStartX */ ; i < iImgEndX ; i++)
        {
	  pt = pImage->mri[i+iImgStartX/nZoom][j+iImgStartY/nZoom];

	  if(pt < imin) pt = imin;
	  if(pt > imax) pt = imax;

	  pt = pt * pImage->rcoeff;

	  /* management of the different types of stages */

	  IMXPICT2D_CALCPIXEL_FIXED
            IMXPICT2D_TEST_THRESHOLD
            IMXPICT2D_STOREPIXEL_GENERIC
	    }
    }
}

//
//  ---- IMXPict2d_8Zoom1_Fixed_Sym() -----------------------------------------
//
//  Description:
//      ...
//
//  >>  Created on Sep 24th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//

static void IMXPict2d_8Zoom1_Fixed_Sym(ptr_grphic pImage, TMaskType mask_type, UCHAR *pDest)
{
  static int i, j, k, l, m, n;
  static LONG start, end;
  static LONG imin, imax;
  static ptrColormap cmap;
  static int nZoom;
  static LONG pix;
  static double pt;
  static BOOL trouve;
  static UCHAR *pDest2;
  static UINT iImgStartX, iImgStartY;
  static UINT iImgEndX, iImgEndY;
  static BOOL bBlack;
  static int OverBlancking;
  static int UnderBlancking;
  IMXPICT2D_INITIALIZEDATA

    pDest2  = pDest + _pic_size_width * _nBytesPerPixel;

  for(j = 0/*iImgStartY */ ; j < iImgEndY ; j++)
    {
      for(i = 0/* iImgStartX */ ; i < iImgEndX ; i++)
        {
	  pt = pImage->mri[i+iImgStartX/nZoom][j+iImgStartY/nZoom];

	  if(pt < imin) pt = imin;
	  if(pt > imax) pt = imax;

	  if(pt <= 0.)
            {
	      end     = (end - start) / 2.;
	      imax    = 0.;
            }
	  else if(pt > 0.)
            {
	      start   = (end - start) / 2.;
	      imin    = 0.;
            }

	  IMXPICT2D_CALCPIXEL_FIXED
            IMXPICT2D_TEST_THRESHOLD
            IMXPICT2D_STOREPIXEL_GENERIC
	    }
    }
}

//
// >>>> 8 BITS ZOOM 2 PART <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//

//
//  ---- IMXPict2d_8Zoom2_Dynamic_Gen() ---------------------------------------
//
//  Description:
//      ...
//
//  >>  Created on Sep 18th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//

static void IMXPict2d_8Zoom2Mask_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pDest)
{
  static int i, j, k, l, m, n;
  static LONG start, end;
  static LONG imin, imax;
  static ptrColormap cmap;
  static int nZoom;
  static LONG pix;
  static double pt;
  static UCHAR *pDest2;
  static UINT iImgStartX, iImgStartY;
  static UINT iImgEndX, iImgEndY;
  static BOOL bBlack;
  static int OverBlancking;
  static int UnderBlancking;
  IMXPICT2D_INITIALIZEDATA

    pDest2  = pDest + _pic_size_width * _nBytesPerPixel;

  for(j = 0/*iImgStartY */ ; j < iImgEndY ; j++)
    {
      for(i = 0/* iImgStartX */ ; i < iImgEndX ; i++)
        {
	  pt = pImage->mri[i+iImgStartX/nZoom][j+iImgStartY/nZoom];

	  if(pt < imin) pt = imin;
	  if(pt > imax) pt = imax;

	  IMXPICT2D_CALCPIXEL_DYNAMIC
            pix = _mask_color;
	  IMXPICT2D_STOREPIXEL_GENERIC
	    }
    }
}

static void IMXPict2d_8Zoom2_Dynamic_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pDest)
{
  static int i, j, k, l, m, n;
  static LONG start, end;
  static LONG imin, imax;
  static ptrColormap cmap;
  static int nZoom;
  static LONG pix;
  static double pt;
  static UCHAR *pDest2;
  static UINT iImgStartX, iImgStartY;
  static UINT iImgEndX, iImgEndY;
  static BOOL bBlack;
  static int OverBlancking;
  static int UnderBlancking;
  IMXPICT2D_INITIALIZEDATA

    pDest2  = pDest + _pic_size_width * _nBytesPerPixel;

  for(j = 0/*iImgStartY */ ; j < iImgEndY ; j++)
    {
      for(i = 0/* iImgStartX */ ; i < iImgEndX ; i++)
        {
	  pt = pImage->mri[i+iImgStartX/nZoom][j+iImgStartY/nZoom];

	  if(pt < imin) pt = imin;
	  if(pt > imax) pt = imax;

	  IMXPICT2D_CALCPIXEL_DYNAMIC
            IMXPICT2D_TEST_THRESHOLD
	    if(_threshold_draw)
	      {
		if(pt >= _maxmin_min && pt <= _maxmin_max)
		  pix = _threshold_color;
	      }
	  IMXPICT2D_STOREPIXEL_GENERIC
	    }
    }
}

//
//  ---- IMXPict2d_8Zoom2_Dynamic_Sym() ---------------------------------------
//
//  Description:
//      ...
//
//  >>  Created on Sep 18th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//

static void IMXPict2d_8Zoom2_Dynamic_Sym(ptr_grphic pImage, TMaskType mask_type, UCHAR *pDest)
{
  static int i, j, k, l, m, n;
  static LONG start, end;
  static LONG imin, imax;
  static ptrColormap cmap;
  static int nZoom;
  static LONG pix;
  static double pt;
  static UCHAR *pDest2;
  static UINT iImgStartX, iImgStartY;
  static UINT iImgEndX, iImgEndY;
  static BOOL bBlack;
  static int OverBlancking;
  static int UnderBlancking;
  IMXPICT2D_INITIALIZEDATA

    pDest2  = pDest + _pic_size_width * _nBytesPerPixel;

  for(j = 0/*iImgStartY */ ; j < iImgEndY ; j++)
    {
      for(i = 0/* iImgStartX */ ; i < iImgEndX ; i++)
        {
	  pt = pImage->mri[i+iImgStartX/nZoom][j+iImgStartY/nZoom];

	  if(pt < imin) pt = imin;
	  if(pt > imax) pt = imax;

	  if(pt <= 0.)
            {
	      end     = (end - start) / 2.;
	      imax    = 0.;
            }
	  else if(pt > 0.)
            {
	      start   = (end - start) / 2.;
	      imin    = 0.;
            }

	  /***********************************************************/
	  if( (imax - imin) + start != 0 )
	    pix = (LONG) ((double) pt - imin) * (end - start)
	      / (imax - imin) + start;
	  else
	    pix = 1l;
	  bBlack = (cmap->tblColors[pix].red == 0
		    && cmap->tblColors[pix].green == 0
		    && cmap->tblColors[pix].blue == 0);
	  pix = cmap->tblColors[pix].pixel;
	  /***********************************************************/

	  //IMXPICT2D_CALCPIXEL_DYNAMIC
	  IMXPICT2D_TEST_THRESHOLD
            IMXPICT2D_STOREPIXEL_GENERIC
	    }
    }
}

//
//  ---- IMXPict2d_8Zoom2_Fixed_Gen() -----------------------------------------
//
//  Description:
//      ...
//
//  >>  Created on Sep 24th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//

static void IMXPict2d_8Zoom2_Fixed_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pDest)
{
  static int i, j, k, l, m, n;
  static LONG start, end;
  static LONG imin, imax;
  static ptrColormap cmap;
  static int nZoom;
  static LONG pix;
  static double pt;
  static BOOL trouve;
  static UCHAR *pDest2;
  static UINT iImgStartX, iImgStartY;
  static UINT iImgEndX, iImgEndY;
  static BOOL bBlack;
  static int OverBlancking;
  static int UnderBlancking;
  IMXPICT2D_INITIALIZEDATA

    pDest2  = pDest + _pic_size_width * _nBytesPerPixel;

  for(j = 0/*iImgStartY */ ; j < iImgEndY ; j++)
    {
      for(i = 0/* iImgStartX */ ; i < iImgEndX ; i++)
        {
	  pt = pImage->mri[i+iImgStartX/nZoom][j+iImgStartY/nZoom];

	  if(pt < imin) pt = imin;
	  if(pt > imax) pt = imax;

	  pt = pt * pImage->rcoeff;

	  /* management of the different types of stages */

	  IMXPICT2D_CALCPIXEL_FIXED
            IMXPICT2D_TEST_THRESHOLD
            IMXPICT2D_STOREPIXEL_GENERIC
	    }
    }
}

//
//  ---- IMXPict2d_8Zoom2_Fixed_Sym() -----------------------------------------
//
//  Description:
//      ...
//
//  >>  Created on Sep 24th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//

static void IMXPict2d_8Zoom2_Fixed_Sym(ptr_grphic pImage, TMaskType mask_type, UCHAR *pDest)
{
  static int i, j, k;
  static LONG start, end;
  static LONG imin, imax;
  static ptrColormap cmap;
  static int nZoom;
  static LONG pix;
  static double pt;
  static BOOL trouve;
  static UCHAR *pDest2;
  static UINT iImgStartX, iImgStartY;
  static UINT iImgEndX, iImgEndY;
  static BOOL bBlack;
  static int OverBlancking;
  static int UnderBlancking;
  IMXPICT2D_INITIALIZEDATA

    pDest2  = pDest + _pic_size_width * _nBytesPerPixel;

  for(j = 0/*iImgStartY */ ; j < iImgEndY ; j++)
    {
      for(i = 0/* iImgStartX */ ; i < iImgEndX ; i++)
        {
	  pt = pImage->mri[i+iImgStartX/nZoom][j+iImgStartY/nZoom];

	  if(pt < imin) pt = imin;
	  if(pt > imax) pt = imax;

	  if(pt <= 0.)
            {
	      end     = (end - start) / 2.;
	      imax    = 0.;
            }
	  else if(pt > 0.)
            {
	      start   = (end - start) / 2.;
	      imin    = 0.;
            }

	  IMXPICT2D_CALCPIXEL_FIXED
            IMXPICT2D_TEST_THRESHOLD
	    }
    }
}

//
// >>>> 8 BITS ZOOM 4 PART <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//

//
//
//  ---- IMXPict2d_8Zoom4_Dynamic_Gen() ---------------------------------------
//
//  Description:
//      ...
//
//  >>  Created on Sep 18th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//

static void IMXPict2d_8Zoom4Mask_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pDest)
{
  static int i, j, k, l, m, n;
  static LONG start, end;
  static LONG imin, imax;
  static ptrColormap cmap;
  static int nZoom;
  static LONG pix;
  static double pt;
  static UCHAR *pDest2, *pDest3, *pDest4;
  static UINT iImgStartX, iImgStartY;
  static UINT iImgEndX, iImgEndY;
  static BOOL bBlack;
  static int OverBlancking;
  static int UnderBlancking;
  IMXPICT2D_INITIALIZEDATA

    pDest2  = pDest + _pic_size_width * _nBytesPerPixel;
  pDest3  = pDest2 + _pic_size_width * _nBytesPerPixel;
  pDest4  = pDest3 + _pic_size_width * _nBytesPerPixel;

  for(j = 0/*iImgStartY */ ; j < iImgEndY ; j++)
    {
      for(i = 0/* iImgStartX */ ; i < iImgEndX ; i++)
        {
	  pt = pImage->mri[i+iImgStartX/nZoom][j+iImgStartY/nZoom];

	  if(pt < imin) pt = imin;
	  if(pt > imax) pt = imax;

	  IMXPICT2D_CALCPIXEL_DYNAMIC
            pix = _mask_color;
	  IMXPICT2D_STOREPIXEL_GENERIC
	    }
    }
}

static void IMXPict2d_8Zoom4_Dynamic_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pDest)
{
  static int i, j, k, l, m, n;
  static LONG start, end;
  static LONG imin, imax;
  static ptrColormap cmap;
  static int nZoom;
  static LONG pix;
  static double pt;
  static UCHAR *pDest2, *pDest3, *pDest4;
  static UINT iImgStartX, iImgStartY;
  static UINT iImgEndX, iImgEndY;
  static BOOL bBlack;
  static int OverBlancking;
  static int UnderBlancking;
  IMXPICT2D_INITIALIZEDATA

    pDest2  = pDest + _pic_size_width * _nBytesPerPixel;
  pDest3  = pDest2 + _pic_size_width * _nBytesPerPixel;
  pDest4  = pDest3 + _pic_size_width * _nBytesPerPixel;

  for(j = 0/*iImgStartY */ ; j < iImgEndY ; j++)
    {
      for(i = 0/* iImgStartX */ ; i < iImgEndX ; i++)
        {
	  pt = pImage->mri[i+iImgStartX/nZoom][j+iImgStartY/nZoom];

	  if(pt < imin) pt = imin;
	  if(pt > imax) pt = imax;

	  IMXPICT2D_CALCPIXEL_DYNAMIC
            IMXPICT2D_TEST_THRESHOLD
            IMXPICT2D_STOREPIXEL_GENERIC
	    }
    }
}

//
//  ---- IMXPict2d_8Zoom4_Dynamic_Sym() ---------------------------------------
//
//  Description:
//      ...
//
//  >>  Created on Sep 18th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//

static void IMXPict2d_8Zoom4_Dynamic_Sym(ptr_grphic pImage, TMaskType mask_type, UCHAR *pDest)
{
  static int i, j, k, l, m, n;
  static LONG start, end;
  static LONG imin, imax;
  static ptrColormap cmap;
  static int nZoom;
  static LONG pix;
  static double pt;
  static UCHAR *pDest2, *pDest3, *pDest4;
  static UINT iImgStartX, iImgStartY;
  static UINT iImgEndX, iImgEndY;
  static BOOL bBlack;
  static int OverBlancking;
  static int UnderBlancking;
  IMXPICT2D_INITIALIZEDATA

    pDest2  = pDest + _pic_size_width * _nBytesPerPixel;
  pDest3  = pDest2 + _pic_size_width * _nBytesPerPixel;
  pDest4  = pDest3 + _pic_size_width * _nBytesPerPixel;

  for(j = 0/*iImgStartY */ ; j < iImgEndY ; j++)
    {
      for(i = 0/* iImgStartX */ ; i < iImgEndX ; i++)
        {
	  pt = pImage->mri[i+iImgStartX/nZoom][j+iImgStartY/nZoom];

	  if(pt < imin) pt = imin;
	  if(pt > imax) pt = imax;

	  if(pt <= 0.)
            {
	      end     = (end - start) / 2.;
	      imax    = 0.;
            }
	  else if(pt > 0.)
            {
	      start   = (end - start) / 2.;
	      imin    = 0.;
            }

	  IMXPICT2D_CALCPIXEL_DYNAMIC
            IMXPICT2D_TEST_THRESHOLD
            IMXPICT2D_STOREPIXEL_GENERIC
	    }
    }
}

//
//  ---- IMXPict2d_8Zoom4_Fixed_Gen() -----------------------------------------
//
//  Description:
//      ...
//
//  >>  Created on Sep 24th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//

static void IMXPict2d_8Zoom4_Fixed_Gen(ptr_grphic pImage, TMaskType mask_type, UCHAR *pDest)
{
  static int i, j, k, l, m, n;
  static LONG start, end;
  static LONG imin, imax;
  static ptrColormap cmap;
  static int nZoom;
  static LONG pix;
  static double pt;
  static BOOL trouve;
  static UCHAR *pDest2, *pDest3, *pDest4;
  static UINT iImgStartX, iImgStartY;
  static UINT iImgEndX, iImgEndY;
  static BOOL bBlack;
  static int OverBlancking;
  static int UnderBlancking;
  IMXPICT2D_INITIALIZEDATA

    pDest2  = pDest + _pic_size_width * _nBytesPerPixel;
  pDest3  = pDest2 + _pic_size_width * _nBytesPerPixel;
  pDest4  = pDest3 + _pic_size_width * _nBytesPerPixel;

  for(j = 0/*iImgStartY */ ; j < iImgEndY ; j++)
    {
      for(i = 0/* iImgStartX */ ; i < iImgEndX ; i++)
        {
	  pt = pImage->mri[i+iImgStartX/nZoom][j+iImgStartY/nZoom];

	  if(pt < imin) pt = imin;
	  if(pt > imax) pt = imax;

	  pt = pt * pImage->rcoeff;

	  /* management of the different types of stages */

	  IMXPICT2D_CALCPIXEL_FIXED
            IMXPICT2D_TEST_THRESHOLD
            IMXPICT2D_STOREPIXEL_GENERIC
	    }
    }
}

//
//  ---- IMXPict2d_8Zoom4_Fixed_Sym() -----------------------------------------
//
//  Description:
//      ...
//
//  >>  Created on Sep 24th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//

static void IMXPict2d_8Zoom4_Fixed_Sym(ptr_grphic pImage, TMaskType mask_type, UCHAR *pDest)
{
  static int i, j, k, l, m, n;
  static LONG start, end;
  static LONG imin, imax;
  static ptrColormap cmap;
  static int nZoom;
  static LONG pix;
  static double pt;
  static BOOL trouve;
  static UCHAR *pDest2, *pDest3, *pDest4;
  static UINT iImgStartX, iImgStartY;
  static UINT iImgEndX, iImgEndY;
  static BOOL bBlack;
  static int OverBlancking;
  static int UnderBlancking;
  IMXPICT2D_INITIALIZEDATA
    
    pDest2  = pDest + _pic_size_width * _nBytesPerPixel;
  pDest3  = pDest2 + _pic_size_width * _nBytesPerPixel;
  pDest4  = pDest3 + _pic_size_width * _nBytesPerPixel;
    
  for(j = 0/*iImgStartY */ ; j < iImgEndY ; j++)
    {
      for(i = 0/* iImgStartX */ ; i < iImgEndX ; i++)
        {
	  pt = pImage->mri[i+iImgStartX/nZoom][j+iImgStartY/nZoom];
            
	  if(pt < imin) pt = imin;
	  if(pt > imax) pt = imax;
            
	  if(pt <= 0.)
            {
	      end     = (end - start) / 2.; 
	      imax    = 0.;
            }
	  else if(pt > 0.)
            {
	      start   = (end - start) / 2.; 
	      imin    = 0.;
            }
          
	  IMXPICT2D_CALCPIXEL_FIXED
            IMXPICT2D_TEST_THRESHOLD
            IMXPICT2D_STOREPIXEL_GENERIC
	    }
    }
}

//
//
//  ---- IMXPict2d_BuildDisplayableImageData() --------------------------------
//
//  Description:
//      Builds displayable image data.
//
//  >>  Created on Sep 18th, 2001 by F. BOULEAU
//  >>  Last user action on ... ..th, 2001 by ...
//

void IMXPict2d_BuildDisplayableImageData(ptr_grphic pImage, UCHAR *pDest, TMaskType mask_type)
{
  int i;
  BOOL found;

  if(mask_type == MASK_OPAQUE )
    {
      memset(pDest,
	     0, 
	     _pic_size_width * _pic_size_height * _nBytesPerPixel);
    }

  // Correction R. Chevrier, le 29 octobre 2001
  /*    if(pImage->cutoff_min == 0 && pImage->cutoff_max == 0 
	|| pImage->cutoff_min > pImage->cutoff_max)
        return;
  */

  if( pImage->cutoff_min > pImage->cutoff_max )
    return;

  for(found = FALSE, i = 0 ; !found ; i++)
    {
      if((tblDisplayFuncs[i].id == pImage->zoom 
	  || tblDisplayFuncs[i].id == GENERIC)
	 && (tblDisplayFuncs[i].mode == pImage->cmapinfo.dspStyle
	     || tblDisplayFuncs[i].mode == DSP_ALL)
	 && (tblDisplayFuncs[i].type == pImage->cmapinfo.cmap->noType
	     || tblDisplayFuncs[i].type == MAP_ALL)
	 && tblDisplayFuncs[i].masktype == mask_type)
        {
	  found = TRUE;
	  tblDisplayFuncs[i].func(pImage, mask_type, pDest);
        }
    }

  if(pImage->mask_active && pImage->mask)
    IMXPict2d_BuildDisplayableImageData(pImage->mask, pDest, pImage->mask_type);
}
#endif /* __GTK_GUI */


/******getmri_convert()************************
 * Cette fonction convertit le buff precedement lu dans getmri()
 * suivant le type de donnees et l'architecture machine/fichier
 *(little ou big Endian)
 * Laurent Burger, decembre 1999
 ************************************************/
long *getmri_convert(long int *buff, void *buffer, int data_size, int format, int num_items)
{
  unsigned char *buffchar=(unsigned char)NULL;
  unsigned short *buffshort=(unsigned short *)NULL;
  int i=0;


  switch(data_size) {
  case sizeof(char):
    buffchar = (unsigned char *)buffer;
    for(i=0; i<num_items; i++)
      buff[i] = (unsigned long)buffchar[i];
    FREE(buffer);
    break;
  case sizeof(short):
    buffshort = (unsigned short *)buffer;
    for(i=0; i<num_items; i++)
      switch(format)
	{
	case ACRNEMA2B_FORMAT :
	  if (!_is_bigEndian) buff[i] = (long) shEndianConversion(buffshort[i]);
	  else buff[i] = (long)buffshort[i];
	  break;

	case ACRNEMA1B_FORMAT :
	  if (!_is_bigEndian) buff[i] = (long) shEndianConversion(buffshort[i]);
	  else buff[i] = (long)buffshort[i];
	  break;

	case ACRNEMA2L_FORMAT :
	  if (_is_bigEndian) buff[i] = (long) shEndianConversion(buffshort[i]);
	  else buff[i] = (long)buffshort[i];
	  break;

	case ACRNEMA1L_FORMAT :
	  if (_is_bigEndian) buff[i] = (long) shEndianConversion(buffshort[i]);
	  else buff[i] = (long)buffshort[i];
	  break;

	case SMIS_FORMAT :
	  if (_is_bigEndian) buff[i] = (long) shEndianConversion(buffshort[i]);
	  else buff[i] = (long)buffshort[i];
	  break;

	case DICMBIMP_FORMAT :
	case DICMLIMP_FORMAT :
	case DICMBEXP_FORMAT :
	case DICMLEXP_FORMAT :
	  if (_is_bigEndian) buff[i] = (long) shEndianConversion(buffshort[i]);
	  else buff[i] = (long)buffshort[i];
	  break;

	default :
	  if (!_is_bigEndian) buff[i] = (long) shEndianConversion(buffshort[i]);
	  else buff[i] = (long)buffshort[i];
	  break;
	}
    FREE(buffer);
    break;

  case sizeof(long):
    /* picture's data is already in the good format */
    buff = (long *)buffer;
    for(i=0; i<num_items; i++)
      switch(format)
	{
	case ACRNEMA2B_FORMAT :
	  if (!_is_bigEndian) buff[i] = (long) longEndianConversion(buff[i]);
	  else buff[i] = (long)buff[i];
	  break;

	case ACRNEMA1B_FORMAT :
	  if (!_is_bigEndian) buff[i] = (long) longEndianConversion(buff[i]);
	  else buff[i] = (long)buff[i];
	  break;

	case ACRNEMA2L_FORMAT :
	  if (_is_bigEndian) buff[i] = (long) longEndianConversion(buff[i]);
	  else buff[i] = (long)buff[i];
	  break;

	case ACRNEMA1L_FORMAT :
	  if (!_is_bigEndian) buff[i] = (long) longEndianConversion(buff[i]);
	  else buff[i] = (long)buff[i];
	  break;

	case SMIS_FORMAT :
	  if (_is_bigEndian) buff[i] = (long) longEndianConversion(buff[i]);
	  else buff[i] = (long)buff[i];
	  break;

	case DICMBIMP_FORMAT :
	case DICMLIMP_FORMAT :
	case DICMBEXP_FORMAT :
	case DICMLEXP_FORMAT :
	  if (_is_bigEndian) buff[i] = (long) longEndianConversion(buff[i]);
	  else buff[i] = (long)buff[i];
	  break;

	default :
	  if (!_is_bigEndian) buff[i] = (long) longEndianConversion(buff[i]);
	  else buff[i] = (long)buff[i];
	  break;
	}	
    /*        FREE(buffer);
	      if(data_size != sizeof(long))
	      FREE(buff);
    */      
  }
  return buff;
}


/*  -- show_picture() ---------------------------------------------------------
**
**  Regenerate picture display number pos
**
**  Arguments:
**  - pos:  index of picture to regenerate
**
**  Return value: -1 if an error occurred, pos otherwise
**
*/
int show_picture(int pos)
{
#ifdef __GTK_GUI
  if(pos < -1 || pos > MAX_PICTURES) 
    {
      _imagix_error   = IMX_BAD_POSITION;
      _imagix_err_ext = pos;
      return(-1);
    }

  if(MESG_DEBUG)
    printf("Show_picture> %d\n", pos);
    
  if(pos > 0 && pos <= MAX_PICTURES) 
    {
      if(MESG_DEBUG)
	printf("Show picture at %d\n", pos);
      IMXWidget_ShowMRI(pos, (grphic *) &_imx_img[pos - 1]);
    }
    
#endif /* __GTK_GUI */
  return(pos);
}

/******************************************************************************
 ** -- Refresh_Order() ---------------------------------------------------------
 **
 ** >> Last user action by: F. BOULEAU on July 22nd, 1999
 **
 ******************************************************************************/
void Refresh_Order(void)
{ 
  /* Redraw all pictures... */
  int i;

  for(i=1;i<=MAX_PICTURES;i++)
    show_picture(i);

  /* TODO:GTK
     if(_imx_roi.pos>0&& _imx_roi.pos<=MAX_PICTURES)
     IMXWidget_ShowROI(_imx_roi.pos,&_imx_roi,0);
  */
}

#ifdef __GTK_GUI

/* --- unzoom_img() ------------------------
***
***
***
***	pour enlever le zoom sur une position
***
***
*** ------------------------------------------- */
int	unzoom_img(int pos)
{

  /* Pas de deplacement de l'image */
  ptr_img(pos)->zoom_x=0;
  ptr_img(pos)->zoom_y=0;
  ptr_img(pos)->zoom=1;
  show_picture(pos);
  return(1);
}

#endif /* __GTK_GUI */
