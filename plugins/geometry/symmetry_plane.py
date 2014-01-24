##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import itk
import medipy.itk

def flip_around_symmetry_plane(input, reference=None) :
    """ Compute the symmetry plane of input image and flip it around this plane.
        If reference is given, compute the symmetry plane of reference, and flip
        input around it.
    
        <gui>
            <item name="input" type="Image" label="Input"/>
            <item name="reference" type="Image" 
                  initializer="may_be_empty=True, may_be_empty_checked=True"
                  label="Reference" />
            <item name="output" type="Image" initializer="output=True" 
                  role="return" label="Output"/>
        </gui>
    """
    
    reference = reference or input
    itk_reference = medipy.itk.medipy_image_to_itk_image(reference, False)
    calculator = itk.SymmetryPlaneCalculator[itk_reference].New(Image=itk_reference)
    calculator.Compute()
    
    transform = calculator.GetPlane().GetReflectionTransform()

    itk_input = medipy.itk.medipy_image_to_itk_image(input, False)
    resample = itk.ResampleImageFilter[itk_input, itk_input].New(
        Input=itk_input, Transform=transform,
        Interpolator=itk.LinearInterpolateImageFunction[itk_input, itk.D].New(),
        DefaultPixelValue=0, ReferenceImage=itk_input, UseReferenceImage=True)
    flipped_itk = resample()[0]

    flipped = medipy.itk.itk_image_to_medipy_image(flipped_itk, None, True)
    return flipped
