/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

%module generate_dicomdir
%{
#include "generate_dicomdir.h"
%}

%include "generate_dicomdir.h"

%pythoncode
%{
import os
import medipy.base

def generate_dicomdir(files, root, dicomdir, patient_extra_attributes=None,
    study_extra_attributes=None, series_extra_attributes=None):
    """ Create a DICOMDIR.
    """
    
    # Perform argument checking
    
    if not hasattr(files, "__getitem__") or not all(isinstance(x, str) for x in files):
        raise medipy.base.Exception("files must be a list of strings")
    
    if not isinstance(root, str):
        raise medipy.base.Exception("root must be a string")
    
    files = [os.path.abspath(x) for x in files]
    root = os.path.join(os.path.abspath(root), "")
    if not all(x.startswith(root) for x in files):
        raise medipy.base.Exception("all files must be under root")
    files = [x[len(root):] for x in files]
    
    if not isinstance(dicomdir, str):
        raise medipy.base.Exception("dicomdir must be a string")
    
    patient_extra_attributes = patient_extra_attributes or []
    if not hasattr(patient_extra_attributes, "__getitem__") or not all(isinstance(x, int) for x in patient_extra_attributes):
        raise medipy.base.Exception("patient_extra_attributes must be a list of ints")
    
    study_extra_attributes = study_extra_attributes or []
    if not hasattr(study_extra_attributes, "__getitem__") or not all(isinstance(x, int) for x in study_extra_attributes):
        raise medipy.base.Exception("study_extra_attributes must be a list of ints")
        
    series_extra_attributes = series_extra_attributes or []
    if not hasattr(series_extra_attributes, "__getitem__") or not all(isinstance(x, int) for x in series_extra_attributes):
        raise medipy.base.Exception("series_extra_attributes must be a list of ints")
    
    # Call C-function
    oldcwd = os.getcwd()
    os.chdir(root)
    generate_dicomdir_cpp(files, root, dicomdir, 
        patient_extra_attributes, study_extra_attributes, 
        series_extra_attributes)
    os.chdir(oldcwd)
%}
