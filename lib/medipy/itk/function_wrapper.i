/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

// Release the GIL while in C-code
%exception {
    Py_BEGIN_ALLOW_THREADS
    $action
    Py_END_ALLOW_THREADS
}

%define ITK_FUNCTION_MACRO(function, docstring, ...)

%pythoncode %{
import medipy.base
from medipy.itk import itk_image_to_medipy_image, medipy_image_to_itk_image

unwrapped_ ## function = function
def function ## ( ##__VA_ARGS__) :
    docstring
    
    args = [ ##__VA_ARGS__]
    itk_args = []
    
    # Transform all arguments of type base.Image to itk.Image, keep the other
    for arg in args :
        if isinstance(arg, medipy.base.Image) :
            itk_image = medipy_image_to_itk_image(arg, False)
            itk_args.append(itk_image)
        else : 
        	itk_args.append(arg)
    
    # Call the C/C++ function
    result = unwrapped_ ## function(*itk_args)
    
    # Transform back all arguments of type itk.Image to base.Image
    for i in range(len(args)) :
        if isinstance(args[i], medipy.base.Image) :
            itk_image_to_medipy_image(itk_args[i], args[i], True)
    
    return result
%}
%enddef