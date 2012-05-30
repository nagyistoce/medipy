/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*******************************************************************************/
/*!	   \file:		detection_anomalies_3d.h
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

#ifndef __DETECTION_ANOMALIES_3D__
#define __DETECTION_ANOMALIES_3D__

enum NORMLINTYPE {REGLIN1,REGLIN2,REGTLS};

typedef double (*lin_energie_func_t)(double, void *);

//-------------Methodes de choix d'un region de reference pour la normalisation-------------------//

//-----Methode par z Map-----//
extern void selection_Z_Map_3d();
extern void imx_selection_Z_Map_3d(int im_1, int im_2);
extern int imx_selection_Z_Map_3d_p(grphic3d *im1, grphic3d *im2);

//-------------Methodes de detection de foyers epileptiques-------------------//

//calcul d'une image soustraction d'un spect ictal et d'un spect interictal recale
//methode de Koo (http://www.blackwell-synergy.com/links/doi/10.1046/j.1528-1157.2003.29402.x/full/)
extern void soustraction_spect_Koo_3d();
extern void imx_soustraction_spect_Koo_3d(int im_ictal, int im_interictal, int im_res);
extern void imx_soustraction_spect_Koo_3d_p(grphic3d *imictal, grphic3d *iminterictal, grphic3d *imres);

//-------------Methodes de normalisation entre images-------------------//

//extern void normalisation_histo_joint_3d(void);
//extern int imx_normalisation_histo_joint_3d(int im_ref, int im_src, int im_res);
//extern int imx_normalisation_histo_joint_3d_p(grphic3d *imref, grphic3d *imsrc, grphic3d *imres);
//
//extern void normalisation_lineaire_masque_3d();
//extern int imx_normalisation_lineaire_masque_3d(int im_ref, int im_src, int im_res, int norm_type);
//extern int imx_normalisation_lineaire_masque_3d_p(grphic3d *imref, grphic3d *imsrc, grphic3d *imres, enum NORMLINTYPE normLinType);
//
////toolbox pour la normalisation lineaire
////regression lineaire a 1 parametre
//extern void regression_lineaire_masque_1_param_3d();
//extern int imx_regression_lineaire_masque_1_param_3d(int im_ref, int im_src);
//extern int imx_regression_lineaire_masque_1_param_3d_p(grphic3d *imref, grphic3d *imsrc, double *alpha, double *beta);
//
////regression lineaire a 2 parametres
//extern void regression_lineaire_masque_2_param_3d();
//extern int imx_regression_lineaire_masque_2_param_3d(int im_ref, int im_src);
//extern int imx_regression_lineaire_masque_2_param_3d_p(grphic3d *imref, grphic3d *imsrc, double *alpha, double *beta);
//
//regression lineaire par total least squares
//extern void regression_TLS_masque_3d();
//extern int imx_regression_TLS_masque_3d(int im_ref, int im_src);
//extern int imx_regression_TLS_masque_3d_p(grphic3d *imref, grphic3d *imsrc, double *alpha, double *beta);
extern void imx_regression_TLS_roi_3d_p(grphic3d *imsrc, grphic3d *imref, grphic3d *imres);
//
//extern int imx_query_normalisation_lineaire_type_3d();
//extern enum NORMLINTYPE imx_choose_normalisation_lineaire_type_3d(int norm_type);

#endif /*__DETECTION_ANOMALIES_3D__*/
