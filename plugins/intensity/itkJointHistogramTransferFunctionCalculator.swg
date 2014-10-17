/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

%include exception.i

%{
#include "vtkPiecewiseFunction.h"
#include "vtkPythonUtil.h"
%}

%{

#include "vtkPythonUtil.h"

#if (VTK_MAJOR_VERSION > 5 ||((VTK_MAJOR_VERSION == 5) && (VTK_MINOR_VERSION > 6)))

#define vtkPythonGetObjectFromPointer vtkPythonUtil::GetObjectFromPointer
#define vtkPythonGetPointerFromObject vtkPythonUtil::GetPointerFromObject

#endif

%}

%typemap(out) vtkPiecewiseFunction* {
    PyImport_ImportModule("vtk");
    $result = vtkPythonGetObjectFromPointer(
        reinterpret_cast<vtkPiecewiseFunction*>($1));
}

%typemap(in) vtkPiecewiseFunction* {
    $1 = NULL;
    $1 = reinterpret_cast<vtkPiecewiseFunction*>(
        vtkPythonGetPointerFromObject($input, "vtkPiecewiseFunction"));
    if($1 == NULL) 
    { 
        SWIG_fail; 
    }
}
