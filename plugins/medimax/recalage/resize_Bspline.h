/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*******************************************************************************/
/*!	   \file:		resize_Bspline.h
***
***		project:	Imagix 2.01
***
***
***		\brief description:    
***
***
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***
********************************************************************************/

#ifndef __RESIZE_BSPLINE_H__
#define __RESIZE_BSPLINE_H__

#include "recalage/definitions_types_recalage.h"

int query_interpolation_zoom_BSpline();
int query_methode_zoom_BSpline();

void compute_zoom_BSpline_3d();
void imx_compute_zoom_BSpline_3d(int im_deb, int im_res, int analyDegree, int syntheDegree, int interpDegree, double zoomX, double zoomY, double zoomZ);
void imx_compute_zoom_BSpline_3d_p(grphic3d *imdeb, grphic3d *imres, int analyDegree, int syntheDegree, int interpDegree, double zoomX, double zoomY, double zoomZ);

int * calculatefinalsize_3d(int width, int height, int depth, double zoomX, double zoomY, double zoomZ);

double beta(double d, int i);
void resamplingLine(double ad[], double ad1[], double ad2[], double ad3[], int i, int j, int *indexMins, int *indexMaxs, double *splineArray, int analyDegree, int corrDegree, int syntheDegree, int ad_length, int ad1_length, int ad2_length, int ad3_length);
double doInteg(double ad[], int i, int ad_length);
void integSA(double ad[], double d, int ad_length);
void integAS(double ad[], double ad1[], int ad_length);
void doDiff(double ad[], int i, int ad_length);
void diffSA(double ad[], int ad_length);
void diffAS(double ad[], int ad_length);
int border(int i, int j);
void getInterpolationCoefficients(double *ad, int i, int ad_length);
void getSamples(double ad[], int i, int ad_length);
double getInitialAntiCausalCoefficient(double ad[], double d, double d1, int ad_length);
double getInitialCausalCoefficient(double ad[], double d, double d1, int ad_length);
void symmetricFir(double ad[], double ad1[], double ad2[], int ad_length, int ad1_length, int ad2_length);



#endif /*__RESIZE_BSPLINE_H__*/
