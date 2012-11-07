
#ifndef UTILS_h
#define UTILS_h

#include <itkVectorImage.h>
#include <itkImage.h>
  
void dtiInvMatrix( itk::VectorImage<float, 3>::Pointer im);
void gradient(itk::Image<float, 3>::Pointer im, itk::VectorImage<float, 3>::Pointer grad);

#endif
