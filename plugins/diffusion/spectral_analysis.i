%module spectral_analysis
%{
#include "Base.includes"
#include "spectral_analysis.h"
%}

%include "../../lib/medipy/itk/function_wrapper.i"

void spectral_analysis(itkVectorImageF3* dt6, 
                       itkVectorImageF3* eigval, 
                       itkVectorImageF3* eigvec);

ITK_FUNCTION_MACRO(spectral_analysis, 
    """ Spectral decomposition of an image of dt6 tensors """, 
    dt6, eigval, eigvec)

