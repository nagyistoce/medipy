##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import itk
import medipy.itk

def fractional_anisotropy(image):
    """ Compute the fractional anisotropy from a second order symmetric tensor field.
    
        <gui>
            <item name="images" type="Image" label="Input"/>
            <item name="output" type="Image" initializer="output=True" 
                  role="return" label="Output"/>
        </gui>
    """
    
    fa_filter = itk.FractionalAnisotropyImageFilter[itk.VectorImage[itk.F,3], itk.Image[itk.F,3]].New()
    
    itk_tensor = medipy.itk.medipy_image_to_itk_image(tensor, False)
    fa_filter.SetInput(0,itk_tensor)
    fa_filter.Update()
    itk_output = fa_filter.GetOutput(0)
    output = medipy.itk.itk_image_to_medipy_image(itk_output,None,True)

    return output

def mean_diffusivity(image):
    """ Compute the mean diffusivity from a second order symmetric tensor field.
    
        <gui>
            <item name="images" type="Image" label="Input"/>
            <item name="output" type="Image" initializer="output=True" 
                  role="return" label="Output"/>
        </gui>
    """
    
    md_filter = itk.MeanDiffusivityImageFilter[itk.VectorImage[itk.F,3], itk.Image[itk.F,3]].New()
    itk_tensor = medipy.itk.medipy_image_to_itk_image(tensor, False)
    md_filter.SetInput(0,itk_tensor)
    md_filter.Update()
    itk_output = md_filter.GetOutput(0)
    output = medipy.itk.itk_image_to_medipy_image(itk_output,None,True)

    return output
