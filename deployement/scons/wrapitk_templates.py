##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

""" Templates for files used in the WrapITK modules.
    
    Unless specified otherwise, all template parameters are strings
"""

# Template for the wrap_class.cxx file, (the Cable wrapper)
# cf. http://public.kitware.com/Cable/HTML/Running.html for more informations
# WrapITK seems to use one such file per wrapped class, although the 
# documentation on Cable says it's possible to do otherwise
# Parameters :
#   * includes
#   * swig_class_name (e.g. itkMedianImageFilter)
#   * typedefs (e.g. 
#     typedef itk::MedianImageFilter<
#         itk::Image<float,3>, itk::Image<float,3> 
#     >::MedianImageFilter itkMedianImageFilterIF3IF3;
#   * forced_instantiations (according to above webpage, "This is necessary to 
#     get GCC-XML to dump proper descriptions of class template instantiations".
#     seldom used in WrapITK (as of WrapITK 0.3.0, it is used in 
#     vnl_unary_function, vcl_complex and itkSpatialObjectBase)
#   * extra : according to the above webpage, "Additional code to ensure 
#     implicit instantiation of class templates that are to be wrapped".

wrap_class_cxx = """%(includes)s

#ifdef CABLE_CONFIGURATION
namespace _cable_
{
   const char* const group = "%(swig_class_name)s";
   namespace wrappers
   {
%(typedefs)s
   }
void force_instantiate()
{
  using namespace _cable_::wrappers;
%(forced_instantiations)s
}
%(extra)s
}
#endif
"""

# Template for the wrap_class_ext.i file.
# Parameters :
#   * swig_class_name
#   * ref_count_declarations (a set of DECLARE_REF_COUNT_CLASS(itkMedianImageFilterIF3IF3)
#     one for each instantiation)
wrap_class_ext_i = """%%import wrap_pyBase.i

%%include wrap_%(swig_class_name)s_doc.i

%(ref_count_declarations)s
"""

# Template for the Module.i file
# Parameters :
#   * module_name
#   * includes (a set of #include "OtherModule.includes", where the current 
#     module depends on OtherModule)
#   * interface_content (seems unused, set to 
#     #"${deps_imports}${SWIG_INTERFACE_MODULE_CONTENT}")
#     in cmake
module_i = """
%%module %(module_name)sPython

%%{
%(includes)s
#include "%(module_name)s.includes"

%%}

%%include itk.i

%%include _%(module_name)s_ext.i

%(interface_content)s
"""

# Template for the Module_ext.i file
# Parameters :
#   * declarations (a list of extern "C" int init_itkCurvatureFlowImageFilterPython();)
#   * calls (a list of init_itkCurvatureFlowImageFilterPython();)
#   * imports (import the modules on which this one depends (import BasePython)
#     and a list of from itkCurvatureFlowImageFilterPython import *)
module_ext_i = """%%{
%(declarations)s
%%}

%%init %%{
%(calls)s
%%}

%%pythoncode %%{
import ITKPyBasePython
%(imports)s
%%}
"""

# Template for the Module.py file. 
# Parameters : 
#   * module_name
module_py = """import itkBase
itkBase.LoadModule("%(module_name)s")
del itkBase
"""