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

cable_swig = """{includes}

#ifdef CABLE_CONFIGURATION
namespace _cable_
{{
   const char* const group = "{swig_class_name}";
   namespace wrappers
   {{
{typedefs}
   }}
void force_instantiate()
{{
  using namespace _cable_::wrappers;
{forced_instantiations}
}}
{extra}
}}
#endif
"""

# Template for the wrap_class_ext.i file.
# Parameters :
#   * swig_class_name
#   * ref_count_declarations (a set of DECLARE_REF_COUNT_CLASS(itkMedianImageFilterIF3IF3)
#     one for each instantiation)
swig_extra = """%import wrap_pyBase.i

// Don't use Doxygen stuff
//%include wrap_{swig_class_name}_doc.i

{ref_count_declarations}
"""
