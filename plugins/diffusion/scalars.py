##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import itk
import medipy.itk
import numpy

def axial_diffusivity(image):
    """ Compute the axial diffusivity from a second order symmetric tensor field.
    
        <gui>
            <item name="image" type="Image" label="Input"/>
            <item name="output" type="Image" initializer="output=True" 
                  role="return" label="Output"/>
        </gui>
    """
    
    return _compute_scalar(image, "AxialDiffusivity")

def fractional_anisotropy(image):
    """ Compute the fractional anisotropy from a second order symmetric tensor field.
    
        <gui>
            <item name="image" type="Image" label="Input"/>
            <item name="output" type="Image" initializer="output=True" 
                  role="return" label="Output"/>
        </gui>
    """
    
    return _compute_scalar(image, "FractionalAnisotropy")

def mean_diffusivity(image):
    """ Compute the mean diffusivity from a second order symmetric tensor field.
    
        <gui>
            <item name="image" type="Image" label="Input"/>
            <item name="output" type="Image" initializer="output=True" 
                  role="return" label="Output"/>
        </gui>
    """
    
    return _compute_scalar(image, "MeanDiffusivity")

def radial_diffusivity(image):
    """ Compute the radial diffusivity from a second order symmetric tensor field.
    
        <gui>
            <item name="image" type="Image" label="Input"/>
            <item name="output" type="Image" initializer="output=True" 
                  role="return" label="Output"/>
        </gui>
    """
    
    return _compute_scalar(image, "RadialDiffusivity")

def _compute_scalar(image, scalar_name):
    """ Compute a scalar from a tensor image.
    """
    
    Filter = getattr(itk, "{0}ImageFilter".format(scalar_name))
    
    itk_image = medipy.itk.medipy_image_to_itk_image(image, False)
    ScalarImage = itk.Image[itk.template(itk_image)[1]]
    
    filter_ = Filter[itk_image, ScalarImage].New(Input=itk_image)
    filter_()
    itk_output = filter_[0]
    
    output = medipy.itk.itk_image_to_medipy_image(itk_output, None, True)
    output.data[numpy.isnan(output.data)] = 0
    output.data[numpy.isinf(output.data)] = 0
    
    return output
