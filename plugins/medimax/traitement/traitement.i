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
%}

%include ../wrapping/grphic3d.i



//-------------------------------------------------------------
//    Correction coupes suivant l'axe Y alternance sombre et claire 
//-------------------------------------------------------------
int imx_corr_sagittal_3d_p(grphic3d *imori,grphic3d *imres);
MEDIMAX_FUNCTION_MACRO(imx_corr_sagittal_3d_p,
"""
Correct Dark Bright Alternation of sections along Y axis

	imori  	: Image to correct
	imres	: resulting image 	
""", imres.data = numpy.ndarray(imori.shape, imori.data.dtype), imori, imres)



