%module itkutils
%{
#include "Base.includes"
#include "itkutils.h"
%}

%include "../../lib/medipy/itk/function_wrapper.i"

void dtiInvMatrix( itkVectorImageF3* im );
void gradient(itkImageF3* im, itkVectorImageF3* grad);

ITK_FUNCTION_MACRO(dtiInvMatrix, """ Invert matrix """, im)
ITK_FUNCTION_MACRO(gradient, """ Compute gradient of an image """, im,grad)


