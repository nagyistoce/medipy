
#ifndef dtiParam_itk_h
#define dtiParam_itk_h

#include <itkImage.h>
#include <itkVectorImage.h>

void parameter_estimation(itk::VectorImage<float, 3>::Pointer dt6, 
    itk::VectorImage<float, 3>::Pointer mean, itk::VectorImage<float, 3>::Pointer var, 
    itk::Image<float, 3>::Pointer mask, 
    unsigned int w_size_plane, unsigned int w_size_depth, bool masked);

#endif
