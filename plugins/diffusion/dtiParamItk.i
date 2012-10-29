%module dtiParamItk
%{
#include "Base.includes"
#include "dtiParamItk.h"
%}

%include "../../lib/medipy/itk/function_wrapper.i"

void dtiParamItk( itkVectorImageF3* dt6, itkVectorImageF3* mean, itkVectorImageF3* var, itkImageF3* mask, unsigned int w_size_plane, unsigned int w_size_depth, bool masked );

ITK_FUNCTION_MACRO(dtiParamItk, """ Compute mean and variance from a tensor 2 image """, dt6,mean,var,mask,w_size_plane,w_size_depth,masked)

