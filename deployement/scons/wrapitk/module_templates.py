# Template for the Module.py file. 
# Parameters : 
#   * module_name
loader = """import itkBase
itkBase.LoadModule("{module_name}")
del itkBase
"""

# Template for the Module.i file
# Parameters :
#   * module_name
#   * includes (a set of #include "OtherModule.includes", where the current 
#     module depends on OtherModule)
#   * interface_content (seems unused, set to 
#     #"${deps_imports}${SWIG_INTERFACE_MODULE_CONTENT}")
#     in cmake
swig = """
%module {module_name}Python

%{{
{includes}
#include "{module_name}.includes"

%}}

%include itk.i

%include _{module_name}_ext.i

{interface_content}
"""

# Template for the Module_ext.i file
# Parameters :
#   * declarations (a list of extern "C" int init_itkCurvatureFlowImageFilterPython();)
#   * calls (a list of init_itkCurvatureFlowImageFilterPython();)
#   * imports (import the modules on which this one depends (import BasePython)
#     and a list of from itkCurvatureFlowImageFilterPython import *)
swig_extra = """%{{
{declarations}
%}}

%init %{{
{calls}
%}}

%pythoncode %{{
import ITKPyBasePython
{imports}
%}}
"""
