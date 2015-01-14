/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifdef SWIGPYTHON

%{
#include <cmath>
#include <numpy/arrayobject.h>
#include <noyau/imx_3d.h>
#include <wrapping/conversion.hpp>
%}

%include exception.i

%typemap(in) ptr_grphic3d (PyObject* pythonImage)
{
    pythonImage = $input;
    $1 = MedimaxWrapper::imageToGrphic3D($input);
}

%typemap(freearg) ptr_grphic3d
{
    MedimaxWrapper::grphic3DToImage($1, pythonImage$argnum);
	free_grphic3d($1);
}

// Release the GIL while in C-code
%exception {
    Py_BEGIN_ALLOW_THREADS
    $action
    Py_END_ALLOW_THREADS
}

%apply ptr_grphic3d { grphic3d* };

%define MEDIMAX_FUNCTION_MACRO(function, docstring, init_code, ...)
// function is the name of the function
// init_code is a piece of python code to initialize the return images
%pythoncode %{
unwrapped_ ## function = function
def function ## ( ##__VA_ARGS__) :
    docstring
    import numpy
    from medipy.base import Image
    init_code
    result = unwrapped_ ## function(##__VA_ARGS__)
    return result
%}
%enddef

#endif /* SWIGPYTHON */
