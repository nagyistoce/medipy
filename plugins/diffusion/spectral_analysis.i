%module spectral_analysis
%{
#include "spectral_analysis.h"
%}

%include "../../lib/medipy/itk/function_wrapper.i"

void spectral_analysis(itk::VectorImage<float, 3>* dt6, 
                       itk::VectorImage<float, 3>* eigval, 
                       itk::VectorImage<float, 3>* eigvec);

ITK_FUNCTION_MACRO(spectral_analysis, 
    """ Spectral decomposition of a DTI image. ``eigval`` and ``eigvec`` are
        output parameters containing respectively the eigenvalues and 
        eigenvectors of the tensor field.""", 
    dt6, eigval, eigvec)
