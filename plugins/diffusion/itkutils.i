%module itkutils
%{
#include "itkutils.h"
%}

%include "../../lib/medipy/itk/function_wrapper.i"

void dtiInvMatrix(itk::VectorImage<float, 3>* im );
void gradient(itk::Image<float, 3>* im, itk::VectorImage<float, 3>* grad);

ITK_FUNCTION_MACRO(dtiInvMatrix, """ Invert matrix """, im)
ITK_FUNCTION_MACRO(gradient, """ Compute gradient of an image """, im,grad)


