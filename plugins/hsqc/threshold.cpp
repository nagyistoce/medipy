#include "threshold.hpp"
#include "itkImage.h" 
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <cstdio>
#include <iostream>
#include <iomanip>

void threshold(itk::Image<float, 2>::Pointer input,float low, float high,float inside, float outside,itk::Image<float, 2>::Pointer output)
{
    // Resize output if needed
 /*
    itk::ImageRegionIterator<itk::Image<float, 2> >
        inputIterator(input, input->GetRequestedRegion());
    itk::ImageRegionIterator<itk::Image<float, 2> >
        outputIterator(output, output->GetRequestedRegion());*/
    //int k;
    //ImageType::IndexType pixelIndexF3;
    //k=input->GetLargestPossibleRegion().GetSize(0);
    itk::Image<float, 2>::IndexType pixelIndexF2;
    for (unsigned int i=30; i<input->GetLargestPossibleRegion().GetSize(0)-30; i++) {
	for (unsigned int j=30; j<input->GetLargestPossibleRegion().GetSize(1)-30; j++) {
		
            pixelIndexF2[0] = i;       
            pixelIndexF2[1] = j;
           //outputIterator.Set(i,j,33);
           output->SetPixel(pixelIndexF2,max(input,i,j,low));
           // printf( "the max is %d \n", output->GetMaximum ());
		
        }
    }

 
    /*while(!inputIterator.IsAtEnd())
    {
        if(inputIterator.Get() >= low && inputIterator.Get() <= high)
        {
            outputIterator.Set(k);
        }
        else
        {
            outputIterator.Set(k);
        }
 
        ++inputIterator;
        ++outputIterator;
    }*/
}
int max(itk::Image<float, 2>::Pointer input, int i,int j,float low)
{   int cond;
    float max;
    itk::Image<float, 2>::IndexType pixelIndexF2;
    pixelIndexF2[0] = i;       
    pixelIndexF2[1] = j;
    cond=1;
    max=input->GetPixel(pixelIndexF2);
    if (max>low){
    //printf( "the max is %d \n", max);
     for (int ii=i-2; ii<i+3; ii++) {
	for (int jj=j-2; jj<j+3; jj++) {
		
            pixelIndexF2[0] = ii;       
            pixelIndexF2[1] = jj;
           
           if(max <= input->GetPixel(pixelIndexF2))
           {
            if(ii==i && jj==j){
                
                
            }
            else{
            cond=0;
            break;
            }
           }
          
		
        }
    }
    }
    else
    {cond=0;}
    return cond;
}