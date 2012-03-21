%module segmentation
%{
#include "Base.includes"
#include "threshold.hpp"
%}

%include "../../lib/medipy/itk/function_wrapper.i"

void threshold(itkImageF2* input,float low, float high,float inside, float outside,itkImageF2* output);
ITK_FUNCTION_MACRO(threshold,"""Binarize an input image by thresholding""",input,low,high,inside,outside,output)
