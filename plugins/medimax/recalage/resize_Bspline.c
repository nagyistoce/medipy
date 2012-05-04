/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*******************************************************************************/
/*!	   \file:		resize_Bspline.c
***
***		project:	Imagix 2.01
***
***
***		\brief description:
***   methode de resize en BSpline avec moindres carres
***   decrite par unser et al.
***   adapte du plugin java:
***   http://bigwww.epfl.ch/algorithms/ijplugins/resize/index.html    
***
***
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***
********************************************************************************/

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_lang.h"
#include "noyau/imx_3d.h"
#include "noyau/io/imx_head.h"

#include "noyau/imx_2d.h"
#include "noyau/extern.h"
#include "noyau/mani_3d.h"
#include "noyau/io/imx_log.h"
#include "math/imx_matrix.h"
#include "recalage/chps_3d.h"
#include "noyau/gui/imx_picture3d.h"

#include "recalage/resize_Bspline.h"

int query_interpolation_zoom_BSpline()
{
	char *query[3];
	query[0]="lineaire";
	query[1]="cubique";
	query[2]=NULL;
	return GETV_QCM("interpolation",(char **)query);
}

int query_methode_zoom_BSpline()
{
	char *query[3];
	query[0]="interpolation";
	query[1]="moindres carres";
	query[2]=NULL;
	return GETV_QCM("methode",(char **)query);
}
 
void compute_zoom_BSpline_3d()
{
 int im_deb, im_res;
 int analyDegree, syntheDegree, interpDegree;
 double zoomX, zoomY, zoomZ;
 int interpolation, methode;
 int err=0;

 im_deb = GET_WND3D("image depart");
 im_res = GET_WND3D("image resultat");

 interpolation = query_interpolation_zoom_BSpline();
 methode = query_methode_zoom_BSpline();

 interpDegree = 3;
 syntheDegree = 3;
 analyDegree = 3;
 if(interpolation == 0)
 {
     interpDegree = 1;
     syntheDegree = 1;
     analyDegree = 1;
 }
 if(methode == 0) analyDegree = -1;

 zoomX=GET_FLOAT("zoomX (reduction 2 : 0.5)", 1.0, &err); if (err) return;
 zoomY=GET_FLOAT("zoomY (reduction 2 : 0.5)", 1.0, &err); if (err) return;
 zoomZ=GET_FLOAT("zoomZ (reduction 2 : 0.5)", 1.0, &err); if (err) return;
    
 imx_compute_zoom_BSpline_3d(im_deb, im_res, analyDegree, syntheDegree, interpDegree, zoomX, zoomY, zoomZ);

 show_picture_3d(im_res);
}

void imx_compute_zoom_BSpline_3d(int im_deb, int im_res, int analyDegree, int syntheDegree, int interpDegree, double zoomX, double zoomY, double zoomZ)
{
	grphic3d *imdeb,*imres;

	imdeb = ptr_img_3d(im_deb);
	imres = ptr_img_3d(im_res);

  imx_compute_zoom_BSpline_3d_p(imdeb, imres, analyDegree, syntheDegree, interpDegree, zoomX, zoomY, zoomZ); 
}

void imx_compute_zoom_BSpline_3d_p(grphic3d *imdeb, grphic3d *imres, int analyDegree, int syntheDegree, int interpDegree, double zoomX,
        double zoomY, double zoomZ)
{
  grphic3d *imtemp=NULL;
  int imdebWidth = imdeb->width;
  int imdebHeight = imdeb->height;
  int imdebDepth = imdeb->depth;
  int l2 = interpDegree + analyDegree + 1;
  int * ai = calculatefinalsize_3d(imdebWidth, imdebHeight, imdebDepth, zoomX, zoomY, zoomZ);
  int reducX = ai[3];
  int reducY = ai[4];
  int reducZ = ai[5];
  int i, j, k, dim_filtre=0;
  double d_index=0.0, d_pow=0.0;
  float ***tmp_res;
  double *bufferBigBordX=NULL, *bufferReducBordX=NULL, *bufferBigBordY=NULL, *bufferReducBordY=NULL, *bufferBigBordZ=NULL, *bufferReducBordZ=NULL;
  double *OutBufferX=NULL, *InBufferX=NULL, *OutBufferY=NULL, *InBufferY=NULL, *OutBufferZ=NULL, *InBufferZ=NULL;
  int *indexMinWidth=NULL, *indexMaxWidth=NULL;
  int *indexMinHeight=NULL, *indexMaxHeight=NULL;
  int *indexMinDepth=NULL, *indexMaxDepth=NULL;
  double *splineArrayWidth=NULL, *splineArrayHeight=NULL, *splineArrayDepth=NULL;
      
  int corrDegree = analyDegree + syntheDegree + 1;
  double halfSupport = ((double)l2 + 1.0) / 2.0;

  int borderX, reducBordX, bigBordX, borderY, reducBordY, bigBordY, borderZ, reducBordZ, bigBordZ;
  int lExtBordY, lExtBordX, sExtBordY, sExtBordX, lExtBordZ, sExtBordZ;
  
  //calcul des coefficients de filtrage en X
  borderX = border(reducX, corrDegree);
  if (borderX < l2) borderX += l2;
  reducBordX = reducX + borderX;
  bigBordX = imdebWidth + (int)ceil((double)borderX / zoomX);
  indexMinWidth = CALLOC(reducBordX, int);
  indexMaxWidth = CALLOC(reducBordX, int);
  dim_filtre = reducBordX * (2 + l2);
  k = 0;
  d_pow = pow(zoomX, analyDegree + 1);
  splineArrayWidth = CALLOC(dim_filtre, double);
  for(i = 0; i < reducBordX; i++)
  {
      d_index = (double)i / zoomX + (((double)analyDegree + 1.0) / 2.0 - floor(((double)analyDegree + 1.0) / 2.0)) * (1.0 / zoomX - 1.0);
      indexMinWidth[i] = (int)ceil(d_index - halfSupport);
      indexMaxWidth[i] = (int)floor(d_index + halfSupport);
      for(j = indexMinWidth[i]; j <= indexMaxWidth[i]; j++)
      {
          splineArrayWidth[k] = d_pow * beta(d_index - (double)j, l2);
          k++;
      }

  }

  //calcul des coefficients de filtrage en Y
  borderY = border(reducY, corrDegree);
  if (borderY < l2) borderY += l2;
  reducBordY = reducY + borderY;
  bigBordY = imdebHeight + (int)ceil((double)borderY / zoomY);
  indexMinHeight = CALLOC(reducBordY, int);
  indexMaxHeight = CALLOC(reducBordY, int);
  dim_filtre = reducBordY * (2 + l2);
  k = 0;
  d_pow = pow(zoomY, analyDegree + 1);
  splineArrayHeight = CALLOC(dim_filtre, double);
  for(i = 0; i < reducBordY; i++)
  {
      d_index = (double)i / zoomY + (((double)analyDegree + 1.0) / 2.0 - floor(((double)analyDegree + 1.0) / 2.0)) * (1.0 / zoomY - 1.0);
      indexMinHeight[i] = (int)ceil(d_index - halfSupport);
      indexMaxHeight[i] = (int)floor(d_index + halfSupport);
      for(j = indexMinHeight[i]; j <= indexMaxHeight[i]; j++)
      {
          splineArrayHeight[k] = d_pow * beta(d_index - (double)j, l2);
          k++;
      }

  }

  //calcul des coefficients de filtrage en Z
  borderZ = border(reducZ, corrDegree);
  if (borderZ < l2) borderZ += l2;
  reducBordZ = reducZ + borderZ;
  bigBordZ = imdebDepth + (int)ceil((double)borderZ / zoomZ);
  indexMinDepth = CALLOC(reducBordZ, int);
  indexMaxDepth = CALLOC(reducBordZ, int);
  dim_filtre = reducBordZ * (2 + l2);
  k = 0;
  d_pow = pow(zoomZ, analyDegree + 1);
  splineArrayDepth = CALLOC(dim_filtre, double);
  for(i = 0; i < reducBordZ; i++)
  {
      d_index = (double)i / zoomZ + (((double)analyDegree + 1.0) / 2.0 - floor(((double)analyDegree + 1.0) / 2.0)) * (1.0 / zoomZ - 1.0);
      indexMinDepth[i] = (int)ceil(d_index - halfSupport);
      indexMaxDepth[i] = (int)floor(d_index + halfSupport);
      for(j = indexMinDepth[i]; j <= indexMaxDepth[i]; j++)
      {
          splineArrayDepth[k] = d_pow * beta(d_index - (double)j, l2);
          k++;
      }

  }

  //calcul des bords de filtrage  
  lExtBordY = 2 * imdebHeight - 2;
  lExtBordX = 2 * imdebWidth - 2;
  sExtBordY = 2 * imdebHeight - 3;
  sExtBordX = 2 * imdebWidth - 3;
  lExtBordZ = 2 * imdebDepth - 2;
  sExtBordZ = 2 * imdebDepth - 3;  

  //---Filtrage---//

  //creation de l'image temporaire en float
  tmp_res=alloc_matrix_3d(reducX, MAXI(reducY, imdebHeight), MAXI(reducZ, imdebDepth));

  //filtrage en X
  bufferBigBordX = CALLOC(bigBordX, double);
  bufferReducBordX = CALLOC(reducBordX, double);
  OutBufferX = CALLOC(reducX, double);
  InBufferX = CALLOC(imdebWidth, double);
  for(j = 0; j < imdebHeight; j++)
   for (k=0; k<imdebDepth; k++)
   {
      for (i=0; i<imdebWidth; i++) InBufferX[i] = (double)imdeb->mri[i][j][k];
      getInterpolationCoefficients(InBufferX, interpDegree, imdebWidth);
      resamplingLine(InBufferX, OutBufferX, bufferBigBordX, bufferReducBordX, lExtBordX, sExtBordX, indexMinWidth, indexMaxWidth, splineArrayWidth, analyDegree, corrDegree, syntheDegree, imdebWidth, reducX, bigBordX, reducBordX);
      for (i=0; i<reducX; i++) tmp_res[i][j][k]=(float)OutBufferX[i];
   }
  FREE(bufferBigBordX); FREE(bufferReducBordX); FREE(OutBufferX); FREE(InBufferX);  

  //filtrage en Y
  bufferBigBordY = CALLOC(bigBordY, double);
  bufferReducBordY = CALLOC(reducBordY, double);
  OutBufferY = CALLOC(reducY, double);
  InBufferY = CALLOC(imdebHeight, double);
  for(i = 0; i < reducX; i++)
   for (k=0; k<imdebDepth; k++)
   {
      for (j=0; j<imdebHeight; j++) InBufferY[j]=tmp_res[i][j][k];
      getInterpolationCoefficients(InBufferY, interpDegree, imdebHeight);
      resamplingLine(InBufferY, OutBufferY, bufferBigBordY, bufferReducBordY, lExtBordY, sExtBordY, indexMinHeight, indexMaxHeight, splineArrayHeight, analyDegree, corrDegree, syntheDegree, imdebHeight, reducY, bigBordY, reducBordY);
      for (j=0; j<reducY; j++) tmp_res[i][j][k]=(float)OutBufferY[j];
   }
  FREE(bufferBigBordY); FREE(bufferReducBordY); FREE(OutBufferY); FREE(InBufferY);
  
   //filtrage en Z
  bufferBigBordZ = CALLOC(bigBordZ, double);
  bufferReducBordZ = CALLOC(reducBordZ, double);
  OutBufferZ = CALLOC(reducZ, double);
  InBufferZ = CALLOC(imdebDepth, double);
  for (i=0; i<reducX; i++)
   for (j=0; j<reducY; j++)
   {
     for (k=0; k<imdebDepth; k++) InBufferZ[k]=tmp_res[i][j][k];
     getInterpolationCoefficients(InBufferZ, interpDegree, imdebDepth);
     resamplingLine(InBufferZ, OutBufferZ, bufferBigBordZ, bufferReducBordZ, lExtBordZ, sExtBordZ, indexMinDepth, indexMaxDepth, splineArrayDepth, analyDegree, corrDegree, syntheDegree, imdebDepth, reducZ, bigBordZ, reducBordZ);
     for (k=0; k<reducZ; k++) tmp_res[i][j][k]=(float)OutBufferZ[k];
   }      
  FREE(bufferBigBordZ); FREE(bufferReducBordZ); FREE(OutBufferZ); FREE(InBufferZ);

 //creation de l'image finale 
 if (imres==imdeb)
 {
  imtemp=cr_grphic3d_modif(reducX, reducY, reducZ, 1.0, 1.0, 0);
  if (!imtemp)
  { fprintf (stderr, "erreur d'allocation memoire dans imx_compute_zoom_BSpline_3d_p\n"); goto end_func; }
  imtemp->dx=((float)imdebWidth/(float)reducX)*imdeb->dx;
  imtemp->dy=((float)imdebHeight/(float)reducY)*imdeb->dy;
  imtemp->dz=((float)imdebDepth/(float)reducZ)*imdeb->dz;
  normaliser_float_3d(imdeb, imtemp, tmp_res, reducX, reducY, reducZ);
  imx_copie_3d_p(imtemp, imres);
 } 
 else
 {
  imx_copie_param_3d_p(imdeb, imres);
  imres->width=reducX; imres->height=reducY; imres->depth=reducZ;
  imres->dx=((float)imdebWidth/(float)reducX)*imdeb->dx;
  imres->dy=((float)imdebHeight/(float)reducY)*imdeb->dy;
  imres->dz=((float)imdebDepth/(float)reducZ)*imdeb->dz;
  //normalisation
  normaliser_float_3d(imdeb, imres, tmp_res, reducX, reducY, reducZ);
 }
 
 imx_inimaxminpixel_3d_p(imres);   

end_func:
 
 // --- Free the temporary image ---
 free_matrix_3d(tmp_res);
 FREE(indexMinWidth); FREE(indexMaxWidth); FREE(indexMinHeight); FREE(indexMaxHeight); FREE(indexMinDepth); FREE(indexMaxDepth);
 FREE(splineArrayWidth); FREE(splineArrayHeight); FREE(splineArrayDepth);
 FREE(ai);
 if (imtemp) free_grphic3d(imtemp);   
}

int * calculatefinalsize_3d(int width, int height, int depth, double zoomX, double zoomY, double zoomZ)
{
  int *finalSizes=CALLOC(6, int);
  finalSizes[0]=width; finalSizes[1]=height; finalSizes[2]=depth;

  finalSizes[3] = (int)floor((double)finalSizes[0] * zoomX + 0.5);
  finalSizes[4] = (int)floor((double)finalSizes[1] * zoomY + 0.5);
  finalSizes[5] = (int)floor((double)finalSizes[2] * zoomZ + 0.5);
  
  return finalSizes;
}

double beta(double d, int i)
{
  double d1 = 0.0;
  switch(i)
  {
  default:
      break;

  case 0: // '\0'
      if(fabs(d) < 0.5)
      {
          d1 = 1.0;
          break;
      }

      if(d == -0.5)
          d1 = 1.0;
      break;

  case 1: // '\001'
      d = fabs(d);
      if(d < 1.0)
          d1 = 1.0 - d;
      break;

  case 2: // '\002'
      d = fabs(d);
      if(d < 0.5)
      {
          d1 = 0.75 - d * d;
          break;
      }
      if(d < 1.5)
      {
          d -= 1.5;
          d1 = d * d * 0.5;
      }
      break;

  case 3: // '\003'
      d = fabs(d);
      if(d < 1.0)
      {
          d1 = d * d * (d - 2.0) * 0.5 + 0.66666666666666663;
          break;
      }
      if(d < 2.0)
      {
          d -= 2.0;
          d1 = d * d * d * -0.16666666666666666;
      }
      break;

  case 4: // '\004'
      d = fabs(d);
      if(d < 0.5)
      {
          d *= d;
          d1 = d * (d * 0.25 - 0.625) + 0.59895833333333337;
          break;
      }
      if(d < 1.5)
      {
          d1 = d * (d * (d * (0.83333333333333337 - d * 0.16666666666666666) - 1.25) + 0.20833333333333334) + 0.57291666666666663;
          break;
      }
      if(d < 2.5)
      {
          d -= 2.5;
          d *= d;
          d1 = d * d * 0.041666666666666664;
      }
      break;

  case 5: // '\005'
      d = fabs(d);
      if(d < 1.0)
      {
          double d2 = d * d;
          d1 = d2 * (d2 * (0.25 - d * 0.083333333333333329) - 0.5) + 0.55000000000000004;
          break;
      }
      if(d < 2.0)
      {
          d1 = d * (d * (d * (d * (d * 0.041666666666666664 - 0.375) + 1.25) - 1.75) + 0.625) + 0.42499999999999999;
          break;
      }
      if(d < 3.0)
      {
          double d3 = 3.0 - d;
          d = d3 * d3;
          d1 = d3 * d * d * 0.0083333333333333332;
      }
      break;

  case 6: // '\006'
      d = fabs(d);
      if(d < 0.5)
      {
          d *= d;
          d1 = d * (d * (0.14583333333333334 - d * 0.027777777777777776) - 0.40104166666666669) + 0.51102430555555556;
          break;
      }
      if(d < 1.5)
      {
          d1 = d * (d * (d * (d * (d * (d * 0.020833333333333332 - 0.14583333333333334) + 0.328125) - 0.12152777777777778) - 0.35546875) - 0.0091145833333333339) + 0.51178385416666672;
          break;
      }
      if(d < 2.5)
      {
          d1 = d * (d * (d * (d * (d * (0.11666666666666667 - d * 0.0083333333333333332) - 0.65625) + 1.8472222222222223) - 2.5703125) + 1.3197916666666667) + 0.17955729166666667;
          break;
      }
      if(d < 3.5)
      {
          d -= 3.5;
          d *= d * d;
          d1 = d * d * 0.0013888888888888889;
      }
      break;

  case 7: // '\007'
      d = fabs(d);
      if(d < 1.0)
      {
          double d4 = d * d;
          d1 = d4 * (d4 * (d4 * (d * 0.0069444444444444441 - 0.027777777777777776) + 0.1111111111111111) - 0.33333333333333331) + 0.47936507936507938;
          break;
      }
      if(d < 2.0)
      {
          d1 = d * (d * (d * (d * (d * (d * (0.050000000000000003 - d * 0.0041666666666666666) - 0.23333333333333334) + 0.5) - 0.3888888888888889) - 0.10000000000000001) - 0.077777777777777779) + 0.49047619047619045;
          break;
      }
      if(d < 3.0)
      {
          d1 = d * (d * (d * (d * (d * (d * (d * 0.0013888888888888889 - 0.027777777777777776) + 0.23333333333333334) - 1.0555555555555556) + 2.7222222222222223) - 3.8333333333333335) + 2.411111111111111) - 0.22063492063492063;
          break;
      }
      if(d < 4.0)
      {
          double d5 = 4.0 - d;
          d = d5 * d5 * d5;
          d1 = d * d * d5 * 0.00019841269841269841;
      }
      break;
  }
  return d1;
}

void resamplingLine(double ad[], double ad1[], double ad2[], double ad3[], int i, int j, int *indexMins, int *indexMaxs, double *splineArray, int analyDegree, int corrDegree, int syntheDegree, int ad_length, int ad1_length, int ad2_length, int ad3_length)
{
  int inc;
  int k = ad_length;
  int l = ad1_length;
  int i1 = ad2_length;
  int j1 = ad3_length;
  double d1 = 0.0;
  int analyEven=0;
  int l1, l2, i3, j3;
  int k1;
 
  if(analyDegree != -1)
      d1 = doInteg(ad, analyDegree + 1, ad_length);
      
  if(((analyDegree + 1) / 2) * 2 == analyDegree + 1) analyEven = 1;
      
//  System.arraycopy(ad, 0, ad2, 0, k);
  for (inc=0;inc<k;inc++) { ad2[inc]=ad[inc]; }
  
  for(l2 = k; l2 < i1; l2++)
      if(analyEven == 1)
      {
          int i2 = l2;
          if(l2 >= i)
              i2 = (int)fabs(fmod(l2, i));
          if(i2 >= k)
              i2 = i - i2;
          ad2[l2] = ad[i2];
      } else
      {
          int j2 = l2;
          if(l2 >= j)
              j2 = (int)fabs(fmod(l2, j));
          if(j2 >= k)
              j2 = j - j2;
          ad2[l2] = -ad[j2];
      }

  k1 = 0;
  for(i3 = 0; i3 < j1; i3++)
  {
      ad3[i3] = 0.0;
      for(j3 = indexMins[i3]; j3 <= indexMaxs[i3]; j3++)
      {
          int k2 = j3;
          double d = 1.0;
          if(j3 < 0)
          {
              k2 = -j3;
              if(analyEven == 0)
              {
                  k2--;
                  d = -1.0;
              }
          }
          if(j3 >= i1)
              k2 = i1 - 1;
          ad3[i3] += d * ad2[k2] * splineArray[k1];
          k1++;
      }

  }

  if(analyDegree != -1)
  {
      doDiff(ad3, analyDegree + 1, ad3_length);
      for(l1 = 0; l1 < j1; l1++)
          ad3[l1] += d1;

      getInterpolationCoefficients(ad3, corrDegree, ad3_length);
      getSamples(ad3, syntheDegree, ad3_length);
  }
  
//  System.arraycopy(ad3, 0, ad1, 0, l);
  for (inc=0;inc<l;inc++) { ad1[inc]=ad3[inc]; }
}


double doInteg(double ad[], int i, int ad_length)
{
  int j = ad_length;
  double d = 0.0;
  double d1 = 0.0;
  int k, l, i1, k1, j1, l1;
  
  switch(i)
  {
  default:
      break;

  case 1: // '\001'
      for(k = 0; k < j; k++)
          d1 += ad[k];

      d1 = (2.0 * d1 - ad[j - 1] - ad[0]) / (double)(2 * j - 2);
      integSA(ad, d1, ad_length);
      break;

  case 2: // '\002'
      for(l = 0; l < j; l++)
          d1 += ad[l];

      d1 = (2.0 * d1 - ad[j - 1] - ad[0]) / (double)(2 * j - 2);
      integSA(ad, d1, ad_length);
      integAS(ad, ad, ad_length);
      break;

  case 3: // '\003'
      for(i1 = 0; i1 < j; i1++)
          d1 += ad[i1];

      d1 = (2.0 * d1 - ad[j - 1] - ad[0]) / (double)(2 * j - 2);
      integSA(ad, d1, ad_length);
      integAS(ad, ad, ad_length);
      for(j1 = 0; j1 < j; j1++)
          d += ad[j1];

      d = (2.0 * d - ad[j - 1] - ad[0]) / (double)(2 * j - 2);
      integSA(ad, d, ad_length);
      break;

  case 4: // '\004'
      for(k1 = 0; k1 < j; k1++)
          d1 += ad[k1];

      d1 = (2.0 * d1 - ad[j - 1] - ad[0]) / (double)(2 * j - 2);
      integSA(ad, d1, ad_length);
      integAS(ad, ad, ad_length);
      for(l1 = 0; l1 < j; l1++)
          d += ad[l1];

      d = (2.0 * d - ad[j - 1] - ad[0]) / (double)(2 * j - 2);
      integSA(ad, d, ad_length);
      integAS(ad, ad, ad_length);
      break;
  }
  return d1;
}

void integSA(double ad[], double d, int ad_length)
{
  int j;
  int i = ad_length;
  ad[0] = (ad[0] - d) * 0.5;
  for(j = 1; j < i; j++)
      ad[j] = (ad[j] - d) + ad[j - 1];
}

void integAS(double ad[], double ad1[], int ad_length)
{
  int j;
  int i = ad_length;
  double *ad2=CALLOC(i, double);
  
//  System.arraycopy(ad, 0, ad2, 0, i);
  for (j=0; j<i; j++) ad2[j]=ad[j];
  
  ad1[0] = ad2[0];
  ad1[1] = 0.0;
  for(j = 2; j < i; j++)
      ad1[j] = ad1[j - 1] - ad2[j - 1];

  FREE(ad2);
}

void doDiff(double ad[], int i, int ad_length)
{
  switch(i)
  {
  case 1: // '\001'
      diffAS(ad, ad_length);
      break;

  case 2: // '\002'
      diffSA(ad, ad_length);
      diffAS(ad, ad_length);
      break;

  case 3: // '\003'
      diffAS(ad, ad_length);
      diffSA(ad, ad_length);
      diffAS(ad, ad_length);
      break;

  case 4: // '\004'
      diffSA(ad, ad_length);
      diffAS(ad, ad_length);
      diffSA(ad, ad_length);
      diffAS(ad, ad_length);
      break;
  }
}

void diffSA(double *ad, int ad_length)
{
  int j;
  int i = ad_length;
  double d = ad[i - 2];
  for(j = 0; j <= i - 2; j++)
      ad[j] = ad[j] - ad[j + 1];

  ad[i - 1] = ad[i - 1] - d;
}

void diffAS(double *ad, int ad_length)
{
  int j;
  int i = ad_length;
  for(j = i - 1; j > 0; j--)
      ad[j] = ad[j] - ad[j - 1];

  ad[0] = 2.0 * ad[0];
}

int border(int i, int j)
{
  int k = i;
  double d=0;
  switch(j)
  {
  case 0: // '\0'
  case 1: // '\001'
      return 0;

  case 2: // '\002'
      d = sqrt(8.0) - 3.0;
      break;

  case 3: // '\003'
      d = sqrt(3.0) - 2.0;
      break;

  case 4: // '\004'
      d = (sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0)) - 19.0;
      break;

  case 5: // '\005'
      d = (sqrt(67.5 - sqrt(4436.25)) + sqrt(26.25)) - 6.5;
      break;

  case 6: // '\006'
      d = -0.48829458930304476;
      break;

  case 7: // '\007'
      d = -0.53528043079643817;
      break;

  default:
      fprintf(stderr, "Invalid interpDegree degree (should be [0..7])");
  }
  k = 2 + (int)floor(log(1.0000000000000001E-009) / log(fabs(d)) + 0.5);
  k = k >= i ? i : k;
  return k;
}

void getInterpolationCoefficients(double *ad, int i, int ad_length)
{
  double *ad1=NULL;
  double d = 1.0;
  int ad1_length=0;
  int j, k, l, i1, j1;

  switch(i)
  {
  case 0: // '\0'
  case 1: // '\001'
      return;

  case 2: // '\002'
      ad1 = CALLOC(1, double); ad1_length=1;
      ad1[0] = sqrt(8.0) - 3.0;
      break;

  case 3: // '\003'
      ad1 = CALLOC(1, double); ad1_length=1;
      ad1[0] = sqrt(3.0) - 2.0;
      break;

  case 4: // '\004'
      ad1 = CALLOC(2, double); ad1_length=2;
      ad1[0] = (sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0)) - 19.0;
      ad1[1] = sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0;
      break;

  case 5: // '\005'
      ad1 = CALLOC(2, double); ad1_length=2;
      ad1[0] = (sqrt(67.5 - sqrt(4436.25)) + sqrt(26.25)) - 6.5;
      ad1[1] = sqrt(67.5 + sqrt(4436.25)) - sqrt(26.25) - 6.5;
      break;

  case 6: // '\006'
      ad1 = CALLOC(3, double); ad1_length=3;
      ad1[0] = -0.48829458930304476;
      ad1[1] = -0.081679271076237514;
      ad1[2] = -0.0014141518083258177;
      break;

  case 7: // '\007'
      ad1 = CALLOC(3, double); ad1_length=3;
      ad1[0] = -0.53528043079643817;
      ad1[1] = -0.12255461519232669;
      ad1[2] = -0.0091486948096082769;
      break;

  default:
      fprintf(stderr, "Invalid spline degree (should be [0..7])");
  }

  if(ad_length == 1)
      { FREE(ad1); return; }
  for(j = 0; j < ad1_length; j++)
      d = d * (1.0 - ad1[j]) * (1.0 - 1.0 / ad1[j]);

  for(k = 0; k < ad_length; k++)
      ad[k] = ad[k] * d;

  for(l = 0; l < ad1_length; l++)
  {
      ad[0] = getInitialCausalCoefficient(ad, ad1[l], 1.0000000000000001E-009, ad_length);
      for(i1 = 1; i1 < ad_length; i1++)
          ad[i1] = ad[i1] + ad1[l] * ad[i1 - 1];

      ad[ad_length - 1] = getInitialAntiCausalCoefficient(ad, ad1[l], 1.0000000000000001E-009, ad_length);
      for(j1 = ad_length - 2; j1 >= 0; j1--)
          ad[j1] = ad1[l] * (ad[j1 + 1] - ad[j1]);

  }

  if (ad1) FREE(ad1);
}

void getSamples(double ad[], int i, int ad_length)
{
  double *ad1 = NULL;
  double *ad2 = CALLOC(ad_length, double);
  int ad1_length=0;
  int j;
  
  switch(i)
  {
  case 0: // '\0'
  case 1: // '\001'
      { FREE(ad2); return; }

  case 2: // '\002'
      ad1 = CALLOC(2, double); ad1_length=2;
      ad1[0] = 0.75;
      ad1[1] = 0.125;
      break;

  case 3: // '\003'
      ad1 = CALLOC(2, double); ad1_length=2;
      ad1[0] = 0.66666666666666663;
      ad1[1] = 0.16666666666666666;
      break;

  case 4: // '\004'
      ad1 = CALLOC(3, double); ad1_length=3;
      ad1[0] = 0.59895833333333337;
      ad1[1] = 0.19791666666666666;
      ad1[2] = 0.0026041666666666665;
      break;

  case 5: // '\005'
      ad1 = CALLOC(3, double); ad1_length=3;
      ad1[0] = 0.55000000000000004;
      ad1[1] = 0.21666666666666667;
      ad1[2] = 0.0083333333333333332;
      break;

  case 6: // '\006'
      ad1 = CALLOC(4, double); ad1_length=4;
      ad1[0] = 0.51102430555555556;
      ad1[1] = 0.22879774305555556;
      ad1[2] = 0.015668402777777778;
      ad1[3] = 2.170138888888889E-005;
      break;

  case 7: // '\007'
      ad1 = CALLOC(4, double); ad1_length=4;
      ad1[0] = 0.47936507936507938;
      ad1[1] = 0.2363095238095238;
      ad1[2] = 0.023809523809523808;
      ad1[3] = 0.00019841269841269841;
      break;

  default:
      fprintf(stderr, "Invalid spline degree (should be [0..7])");
  }
  symmetricFir(ad1, ad, ad2, ad1_length, ad_length, ad_length);
//  System.arraycopy(ad2, 0, ad, 0, ad_length);
  for (j=0; j<ad_length; j++) ad[j]=ad2[j];

  if (ad1) FREE(ad1);
  FREE(ad2);
}

double getInitialAntiCausalCoefficient(double ad[], double d, double d1, int ad_length)
{
  return ((d * ad[ad_length - 2] + ad[ad_length - 1]) * d) / (d * d - 1.0);
}

double getInitialCausalCoefficient(double ad[], double d, double d1, int ad_length)
{
  double d2 = d;
  double d3 = pow(d, ad_length - 1);
  double d4 = ad[0] + d3 * ad[ad_length - 1];
  int i = ad_length;
  int j;
  
  if(d1 > 0.0)
  {
      i = 2 + (int)(log(d1) / log(fabs(d)));
      i = i >= ad_length ? ad_length : i;
  }
  d3 *= d3;
  for(j = 1; j < i - 1; j++)
  {
      d3 /= d;
      d4 += (d2 + d3) * ad[j];
      d2 *= d;
  }

  return d4 / (1.0 - pow(d, 2 * ad_length - 2));
}

void symmetricFir(double ad[], double ad1[], double ad2[], int ad_length, int ad1_length, int ad2_length)
{
 int i, j, k;
 
 if(ad1_length != ad2_length)
     fprintf(stderr, "Incompatible size");
 switch(ad_length)
 {
 case 2: // '\002'
     if(ad1_length >= 2)
     {
         ad2[0] = ad[0] * ad1[0] + 2.0 * ad[1] * ad1[1];
         for(i = 1; i < ad1_length - 1; i++)
             ad2[i] = ad[0] * ad1[i] + ad[1] * (ad1[i - 1] + ad1[i + 1]);

         ad2[ad2_length - 1] = ad[0] * ad1[ad1_length - 1] + 2.0 * ad[1] * ad1[ad1_length - 2];
         break;
     }
     switch(ad1_length)
     {
     case 1: // '\001'
         ad2[0] = (ad[0] + 2.0 * ad[1]) * ad1[0];
         break;

     default:
         fprintf(stderr, "Invalid length of data");
     }
     break;

 case 3: // '\003'
     if(ad1_length >= 4)
     {
         ad2[0] = ad[0] * ad1[0] + 2.0 * ad[1] * ad1[1] + 2.0 * ad[2] * ad1[2];
         ad2[1] = ad[0] * ad1[1] + ad[1] * (ad1[0] + ad1[2]) + ad[2] * (ad1[1] + ad1[3]);
         for(j = 2; j < ad1_length - 2; j++)
             ad2[j] = ad[0] * ad1[j] + ad[1] * (ad1[j - 1] + ad1[j + 1]) + ad[2] * (ad1[j - 2] + ad1[j + 2]);

         ad2[ad2_length - 2] = ad[0] * ad1[ad1_length - 2] + ad[1] * (ad1[ad1_length - 3] + ad1[ad1_length - 1]) + ad[2] * (ad1[ad1_length - 4] + ad1[ad1_length - 2]);
         ad2[ad2_length - 1] = ad[0] * ad1[ad1_length - 1] + 2.0 * ad[1] * ad1[ad1_length - 2] + 2.0 * ad[2] * ad1[ad1_length - 3];
         break;
     }
     switch(ad1_length)
     {
     case 3: // '\003'
         ad2[0] = ad[0] * ad1[0] + 2.0 * ad[1] * ad1[1] + 2.0 * ad[2] * ad1[2];
         ad2[1] = ad[0] * ad1[1] + ad[1] * (ad1[0] + ad1[2]) + 2.0 * ad[2] * ad1[1];
         ad2[2] = ad[0] * ad1[2] + 2.0 * ad[1] * ad1[1] + 2.0 * ad[2] * ad1[0];
         break;

     case 2: // '\002'
         ad2[0] = (ad[0] + 2.0 * ad[2]) * ad1[0] + 2.0 * ad[1] * ad1[1];
         ad2[1] = (ad[0] + 2.0 * ad[2]) * ad1[1] + 2.0 * ad[1] * ad1[0];
         break;

     case 1: // '\001'
         ad2[0] = (ad[0] + 2.0 * (ad[1] + ad[2])) * ad1[0];
         break;

     default:
         fprintf(stderr, "Invalid length of data");
     }
     break;

 case 4: // '\004'
     if(ad1_length >= 6)
     {
         ad2[0] = ad[0] * ad1[0] + 2.0 * ad[1] * ad1[1] + 2.0 * ad[2] * ad1[2] + 2.0 * ad[3] * ad1[3];
         ad2[1] = ad[0] * ad1[1] + ad[1] * (ad1[0] + ad1[2]) + ad[2] * (ad1[1] + ad1[3]) + ad[3] * (ad1[2] + ad1[4]);
         ad2[2] = ad[0] * ad1[2] + ad[1] * (ad1[1] + ad1[3]) + ad[2] * (ad1[0] + ad1[4]) + ad[3] * (ad1[1] + ad1[5]);
         for(k = 3; k < ad1_length - 3; k++)
             ad2[k] = ad[0] * ad1[k] + ad[1] * (ad1[k - 1] + ad1[k + 1]) + ad[2] * (ad1[k - 2] + ad1[k + 2]) + ad[3] * (ad1[k - 3] + ad1[k + 3]);

         ad2[ad2_length - 3] = ad[0] * ad1[ad1_length - 3] + ad[1] * (ad1[ad1_length - 4] + ad1[ad1_length - 2]) + ad[2] * (ad1[ad1_length - 5] + ad1[ad1_length - 1]) + ad[3] * (ad1[ad1_length - 6] + ad1[ad1_length - 2]);
         ad2[ad2_length - 2] = ad[0] * ad1[ad1_length - 2] + ad[1] * (ad1[ad1_length - 3] + ad1[ad1_length - 1]) + ad[2] * (ad1[ad1_length - 4] + ad1[ad1_length - 2]) + ad[3] * (ad1[ad1_length - 5] + ad1[ad1_length - 3]);
         ad2[ad2_length - 1] = ad[0] * ad1[ad1_length - 1] + 2.0 * ad[1] * ad1[ad1_length - 2] + 2.0 * ad[2] * ad1[ad1_length - 3] + 2.0 * ad[3] * ad1[ad1_length - 4];
         break;
     }
     switch(ad1_length)
     {
     case 5: // '\005'
         ad2[0] = ad[0] * ad1[0] + 2.0 * ad[1] * ad1[1] + 2.0 * ad[2] * ad1[2] + 2.0 * ad[3] * ad1[3];
         ad2[1] = ad[0] * ad1[1] + ad[1] * (ad1[0] + ad1[2]) + ad[2] * (ad1[1] + ad1[3]) + ad[3] * (ad1[2] + ad1[4]);
         ad2[2] = ad[0] * ad1[2] + (ad[1] + ad[3]) * (ad1[1] + ad1[3]) + ad[2] * (ad1[0] + ad1[4]);
         ad2[3] = ad[0] * ad1[3] + ad[1] * (ad1[2] + ad1[4]) + ad[2] * (ad1[1] + ad1[3]) + ad[3] * (ad1[0] + ad1[2]);
         ad2[4] = ad[0] * ad1[4] + 2.0 * ad[1] * ad1[3] + 2.0 * ad[2] * ad1[2] + 2.0 * ad[3] * ad1[1];
         break;

     case 4: // '\004'
         ad2[0] = ad[0] * ad1[0] + 2.0 * ad[1] * ad1[1] + 2.0 * ad[2] * ad1[2] + 2.0 * ad[3] * ad1[3];
         ad2[1] = ad[0] * ad1[1] + ad[1] * (ad1[0] + ad1[2]) + ad[2] * (ad1[1] + ad1[3]) + 2.0 * ad[3] * ad1[2];
         ad2[2] = ad[0] * ad1[2] + ad[1] * (ad1[1] + ad1[3]) + ad[2] * (ad1[0] + ad1[2]) + 2.0 * ad[3] * ad1[1];
         ad2[3] = ad[0] * ad1[3] + 2.0 * ad[1] * ad1[2] + 2.0 * ad[2] * ad1[1] + 2.0 * ad[3] * ad1[0];
         break;

     case 3: // '\003'
         ad2[0] = ad[0] * ad1[0] + 2.0 * (ad[1] + ad[3]) * ad1[1] + 2.0 * ad[2] * ad1[2];
         ad2[1] = ad[0] * ad1[1] + (ad[1] + ad[3]) * (ad1[0] + ad1[2]) + 2.0 * ad[2] * ad1[1];
         ad2[2] = ad[0] * ad1[2] + 2.0 * (ad[1] + ad[3]) * ad1[1] + 2.0 * ad[2] * ad1[0];
         break;

     case 2: // '\002'
         ad2[0] = (ad[0] + 2.0 * ad[2]) * ad1[0] + 2.0 * (ad[1] + ad[3]) * ad1[1];
         ad2[1] = (ad[0] + 2.0 * ad[2]) * ad1[1] + 2.0 * (ad[1] + ad[3]) * ad1[0];
         break;

     case 1: // '\001'
         ad2[0] = (ad[0] + 2.0 * (ad[1] + ad[2] + ad[3])) * ad1[0];
         break;

     default:
         fprintf(stderr, "Invalid length of data");
     }
     break;

 default:
     fprintf(stderr, "Invalid filter half-length (should be [2..4])");
 }
}



