##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import itk
import medipy.itk

def flip_around_symmetry_plane(input) :
    """ Flip the input image around its symmetry plane.
    
        <gui>
            <item name="input" type="Image" label="Input"/>
            <item name="output" type="Image" initializer="output=True" 
                  role="return" label="Output"/>
        </gui>
    """
    itk_input = medipy.itk.medipy_image_to_itk_image(input, False)

    symmetry_plane = itk.SymmetryPlaneCalculator[itk_input].New(Image=itk_input)
    symmetry_plane.Compute()
    transform = symmetry_plane.GetPlane().GetReflectionTransform()

    resample = itk.ResampleImageFilter[itk_input, itk_input].New(
        Input=itk_input, Transform=transform,
        Interpolator=itk.LinearInterpolateImageFunction[itk_input, itk.D].New(),
        DefaultPixelValue=0, ReferenceImage=itk_input, UseReferenceImage=True)
    flipped_itk = resample()[0]

    flipped = medipy.itk.itk_image_to_medipy_image(flipped_itk, None, True)
    return flipped
