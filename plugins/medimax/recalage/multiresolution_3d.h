/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/******************************************************************************
**/	
/*!	\file:		multiresolution_3d.h
***
***			
***
***	\brief description:    Recalage rigide avec multiresolution
***	
***	
***	Copyright (c) 2003, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***
********************************************************************************/

#ifndef __MULTIRESOLUTION_3D_H__
#define __MULTIRESOLUTION_3D_H__

extern int imx_ReSampleDimension_3d_p(grphic3d *imsrc, grphic3d *imdst, float dx, float dy, float dz, int inter_type, int bFiltrage);
extern int imx_ReSample_3d_p(grphic3d *imsrc, grphic3d *imdst, int inter_type, float nZoomOut);

//-----------reduction d'une image avec un filtrage de Butterworth----------//
/*
extern int imx_ReSample_Butterworth_3d_p(grphic3d *imsrc, grphic3d *imdst, double nZoomOut);
extern int imx_ReSample_Butterworth_3d(int im_src, int im_dst, double nZoomOut);
extern void ReSample_Butterworth_3d();*/

//---------------reduction d'une image avec un filtrage gaussien-------------//

extern int imx_ReSample_gaussien_3d_p(grphic3d *imsrc, grphic3d *imdst, float nZoomOut);
extern int imx_ReSample_gaussien_3d(int im_src, int im_dst, float nZoomOut);
extern void ReSample_gaussien_3d();

extern int imx_filtrage_gaussien_relatif_3d_p(grphic3d *imref, grphic3d *imsrc, grphic3d *imdst);
extern int imx_filtrage_gaussien_relatif_3d(int im_ref, int im_src, int im_dst);
extern void filtrage_gaussien_relatif_3d();

extern int pyramide_gaussienne_recalage_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imrefres, grphic3d *imrecares, double subsamplingfactor);

//reduit une image trop grosse a la taille maximale admise pour le recalage
extern int imx_reduire_image_recalage(grphic3d *im);

//----reduit une image par un facteur puissance de 2 par la methode de reduction BSpline+moindres carres d'Unser-------//

//--toolbox--//
void PyramidFilterSplinel2(double g[],long *ng,double *h,long *nh,long Order);
void PyramidFilterSplineL2(double g[],long *ng,double h[],long *nh,long Order);
void PyramidFilterCentered(double g[],long *ng,double h[],long *nh,long Order);
void PyramidFilterCenteredL2(double g[],long *ng,double h[],long *nh,long Order);
void PyramidFilterCenteredL2Derivate(double g[],long *ng,double h[],long *nh,long Order);
int GetPyramidFilter(char *Filter, long Order, double g[], long *ng, double h[], long *nh, short *IsCentered);
void ReduceCentered_1D(	double In[], long NxIn, double Out[], double g[], long ng);
void ReduceStandard_1D(double In[], long NxIn, double Out[], double g[], long ng);
void Reduce_1D(double In[], long NxIn, double Out[], double g[], long ng, short IsCentered);

//--resampling--//
extern int imx_ReSample_BSpline_3d_p(grphic3d *imsrc, grphic3d *imdst, int nbResolX, int nbResolY, int nbResolZ);
extern int imx_ReSample_BSpline_3d(int im_src, int im_dst,  int nbResolX, int nbResolY, int nbResolZ);
extern void ReSample_BSpline_3d();

extern void ReSample_3d();
extern int imx_ReSample_3d(int im_src, int im_dst, int inter_type, float nZoomOut);
extern int imx_ReSample_3d_p(grphic3d *imsrc, grphic3d *imdst, int inter_type, float nZoomOut);

#endif /* __MULTIRESOLUTION_3D_H__*/



