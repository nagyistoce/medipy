/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "itkutils.h"

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_math.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/vnl_inverse.h>
#include <vnl/vnl_vector_ref.h>

#include "itkImageRegionIterator.h"

void dtiInvMatrix( itk::VectorImage<float, 3>::Pointer im )
{ 
    vnl_matrix<double> matrix(3,3);
    itk::ImageRegionIterator< itk::VectorImage<float, 3> > it(im, im->GetLargestPossibleRegion());
    it.GoToBegin();

    while( !it.IsAtEnd() ) {  
		
        itk::VectorImage<float, 3>::PixelType vec = it.Get();

	    matrix[0][0] =  vec[0];
	    matrix[0][1] =  vec[1];
	    matrix[0][2] =  vec[2];
	    matrix[1][0] =  vec[3];
	    matrix[1][1] =  vec[4];
	    matrix[1][2] =  vec[5];
	    matrix[2][0] =  vec[6];
	    matrix[2][1] =  vec[7];
	    matrix[2][2] =  vec[8];			  			
	
	    matrix = vnl_inverse<double>( matrix );

    	vec[0] = matrix(0,0);
    	vec[1] = matrix(0,1);
    	vec[2] = matrix(0,2);
    	vec[3] = matrix(1,0);
    	vec[4] = matrix(1,1);
    	vec[5] = matrix(1,2);
   		vec[6] = matrix(2,0);
    	vec[7] = matrix(2,1);
    	vec[8] = matrix(2,2);	

        ++it;
    }
}


void gradient(itk::Image<float, 3>::Pointer im, itk::VectorImage<float, 3>::Pointer grad)
{
    int fx,fy,fz,bx,by,bz;
    float dx,dy,dz,v1,v2;
    itk::Image<float, 3>::SizeType size = im->GetLargestPossibleRegion().GetSize();
    itk::Image<float, 3>::IndexType pixelIndex;
    itk::VectorImage<float, 3>::IndexType voxelIndex;
    assert( grad->GetNumberOfComponentsPerPixel()==3 );

    for ( int z=0; z<(int)size[2]; z++) {
        fz = z == ((int)size[2]-1) ? 0 :  1;
        bz = z == 0 ? 0 : -1;
        dz = (fz==0 || bz==0) ? 1.0 :  2.0;

        for ( int y=0; y<(int)size[1]; y++) {
            fy = y ==((int)size[1]-1) ? 0 :  1;
            by = y == 0 ? 0 : -1;
            dy = (fy==0 || by==0) ? 1.0 :  2.0;

            for ( int x=0; x<(int)size[0]; x++) {
                fx = x ==((int)size[0]-1) ? 0 :  1;
                bx = x == 0 ? 0 : -1;
                dx = (fx==0 || bx==0) ? 1.0 :  2.0;

                voxelIndex[0] = x; voxelIndex[1] = y; voxelIndex[2] = z;
                itk::VariableLengthVector<float> g = grad->GetPixel(voxelIndex);

                pixelIndex[0] = x+fx; pixelIndex[1] = y; pixelIndex[2] = z;
                v1 = im->GetPixel(pixelIndex);
                pixelIndex[0] = x+bx; pixelIndex[1] = y; pixelIndex[2] = z;
                v2 = im->GetPixel(pixelIndex);		
                g[0] = (v1 - v2)/dx;

                pixelIndex[0] = x; pixelIndex[1] = y+fy; pixelIndex[2] = z;
                v1 = im->GetPixel(pixelIndex);
                pixelIndex[0] = x; pixelIndex[1] = y+by; pixelIndex[2] = z;
                v2 = im->GetPixel(pixelIndex);
                g[1] = (v1 - v2)/dy;

                pixelIndex[0] = x; pixelIndex[1] = y; pixelIndex[2] = z+fz;
                v1 = im->GetPixel(pixelIndex);
                pixelIndex[0] = x; pixelIndex[1] = y; pixelIndex[2] = z+bz;
                v2 = im->GetPixel(pixelIndex);
                g[2] = (v1 - v2)/dz;
            }
        }
    }
}

