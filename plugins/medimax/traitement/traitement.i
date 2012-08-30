/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

%module traitement
%{
#include "correction.h"
#include "fonctions_wrappees_traitement.h"
#include "trai_3d.h"

%}

%include ../wrapping/grphic3d.i



//-------------------------------------------------------------
//    Correction coupes suivant l'axe Y alternance sombre et claire 
//-------------------------------------------------------------
int imx_corr_sagittal_3d_p(grphic3d *imori,grphic3d *imres);
MEDIMAX_FUNCTION_MACRO(imx_corr_sagittal_3d_p,
"""
Correct Dark Bright Alternation of sections along Y axis

    imori   : Image to correct
    imres   : resulting image   
""", , imori, imres)


//-------------------------------------------------------------
//     Normalize intensity between two registered images 
//-------------------------------------------------------------
void IntensityNormalisationForRegisteredImages(grphic3d *imsrc, grphic3d *imref, grphic3d *imres, grphic3d *mask_imsrc,grphic3d * mask_imref, int method);
MEDIMAX_FUNCTION_MACRO(IntensityNormalisationForRegisteredImages,
"""
Normalize intensity between two registered images

    imsrc : Image to correct
    imref : Reference image
    imres : Corrected image
    mask_imsrc : mask of Image to correct
    mask_imref : mask of Reference image
    method : 0->oneParameterLinearRegression, 1->twoParameterLinearRegression, 2->totalLeastSquare
    
""",imres.copy_information(imsrc), imsrc, imref, imres, mask_imsrc, mask_imref, method )


//-------------------------------------------------------------
//     Computer chamfer distance transform
//-------------------------------------------------------------
void imx_chamfer_distance_3d_p(grphic3d *imdeb, grphic3d *imres);
MEDIMAX_FUNCTION_MACRO(imx_chamfer_distance_3d_p,
"""
Computer chamfer distance transform
imdeb : Input image
imres : Output image

""",, imdeb, imres)

