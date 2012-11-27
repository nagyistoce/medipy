/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef FONCTIONS_WRAPPEE_H
#define FONCTIONS_WRAPPEE_H

#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/


#include "recalage/definitions_types_recalage.h"

DllExport int LinearRegistration(grphic3d *imref,grphic3d *imreca,grphic3d *imres, int registration_type,int dist_type,int inter_type,int min_type,int multistart,int save_type,char *nomfichres,int precision,int start_resol, int end_resol);

DllExport int ApplyTransfo3d(grphic3d *imdeb, char *nomfichier, grphic3d *imres, int inter_type);

DllExport double SimilarityMeasure3d(grphic3d *im1, grphic3d *im2, int err_type);


DllExport int BsplineRegistration3d(grphic3d *imref, grphic3d *imreca, grphic3d *imres,int dist_type,int reg_type,double reg_factor,int min_type,int save_type,char *nomfichres,int resolf, double Jmin,double Jmax,int normalisation_type, int nb_classe,int biais, int symetrique);


DllExport int CombineTransfo3d(char *nomfichier1, char *nomfichier2, char *nomfichierres, int interpolation);

DllExport int InvertTransfo3d(char *nomfichier, char *nomfichres, int wdthref, int hghtref, int dpthref, float dxref, float dyref, float dzref, float prec);

DllExport int LoadTrfFile(char *nomfichier, grphic3d *ux, grphic3d *uy,grphic3d *uz, grphic3d *imSource);
DllExport int SaveTrfFile(char *nomfichier, grphic3d *ux, grphic3d *uy,grphic3d *uz, grphic3d *imSource);

DllExport void MriInfo3D(grphic3d *p);
DllExport void VisuTrfFile(char *nomfichier, grphic3d *output, int type);

DllExport void simulationAtrophie(grphic3d *imref,grphic3d *mask, char *nomfichres, int resolf, float lambda);

#ifdef __cplusplus
}
#endif /*__cplusplus*/

#endif
