#ifndef segmentation_threshold_hpp
#define segmentation_threshold_hpp
 
#include <itkImage.h>
 
void threshold(itk::Image<float, 2>::Pointer input,
               float low, float high,
               float inside, float outside,
               itk::Image<float, 2>::Pointer output);
int max(itk::Image<float, 2>::Pointer input, int i,int j,float low); 
#endif // segmentation_threshold_hpp