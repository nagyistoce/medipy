/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

%module mask
%{
#include "Base.includes"
#include "mask.h"
%}

%include "../../lib/itk/function_wrapper.i"

void apply_mask(itkImageF3* input,
        itkImageF3* mask, 
        float background, float outside,
         itkImageF3* output);
ITK_FUNCTION_MACRO(apply_mask, """ Apply mask
        
        <gui>
            <item name="input" type="Image" label="Input"/>
            <item name="mask" type="Image" label="Mask"/>
            <item name="background" type="Float" initializer="0" label="Background"/>
            <item name="outside" type="Float" initializer="0" label="Outside"/>
            <item name="output" type="Image" initializer="output = True"
                role="output" label="Output"/>
        </gui>
    """, input, mask, background, outside, output)
