##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import itk
import medipy.itk
import numpy as np

def fractional_anisotropy(image):
    """ Compute the fractional anisotropy from a second order symmetric tensor field.
    
        <gui>
            <item name="image" type="Image" label="Input"/>
            <item name="output" type="Image" initializer="output=True" 
                  role="return" label="Output"/>
        </gui>
    """
    
    itk_image = medipy.itk.medipy_image_to_itk_image(image, False)
    ScalarImage = itk.Image[itk.template(itk_image)[1]]
    
    fa_filter = itk.FractionalAnisotropyImageFilter[itk_image, ScalarImage].New()
    fa_filter.SetInput(0,itk_image)
    fa_filter.Update()
    itk_output = fa_filter.GetOutput(0)
    output = medipy.itk.itk_image_to_medipy_image(itk_output,None,True)
    output.data[np.isnan(output.data)] = 0
    output.data[np.isinf(output.data)] = 0

    return output

def mean_diffusivity(image):
    """ Compute the mean diffusivity from a second order symmetric tensor field.
    
        <gui>
            <item name="image" type="Image" label="Input"/>
            <item name="output" type="Image" initializer="output=True" 
                  role="return" label="Output"/>
        </gui>
    """
    
    itk_image = medipy.itk.medipy_image_to_itk_image(image, False)
    ScalarImage = itk.Image[itk.template(itk_image)[1]]
    
    md_filter = itk.MeanDiffusivityImageFilter[itk_image, ScalarImage].New()
    md_filter.SetInput(0,itk_image)
    md_filter.Update()
    itk_output = md_filter.GetOutput(0)
    output = medipy.itk.itk_image_to_medipy_image(itk_output,None,True)

    return output
