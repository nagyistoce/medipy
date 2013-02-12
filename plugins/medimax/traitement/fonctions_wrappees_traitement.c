/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "fonctions_wrappees_traitement.h"
#include "norm_3d.h"

/*******************************************************************************
**  IntensityNormalisationForRegisteredImages
*******************************************************************************/
void IntensityNormalisationForRegisteredImages(grphic3d *imsrc, grphic3d *imref, grphic3d *imres, grphic3d *mask_imsrc,grphic3d * mask_imref, int method)
{
    norm_func_t normalisation=NULL;
    imsrc->mask=mask_imsrc;
    imref->mask=mask_imref;
  
    normalisation=imx_choose_normalisation_roi_reca_fct(method);
    normalisation(imsrc,imref,imres);  
    imsrc->mask=NULL;
    imref->mask=NULL;
    
}
