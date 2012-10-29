import itk
import medipy.itk
import numpy as np
from medipy.base import Image
from medipy.diffusion.utils import spectral_decomposition,voxel_parameters

def voxel_test(tensor1,tensor2,*args,**kwargs):
    """Multivariate Statistical Tests at a voxel level.

    <gui>
        <item name="tensor1" type="Image" label="Input"/>
        <item name="tensor2" type="Image" label="Input"/>
        <item name="size_plane" type="Int" initializer="3" label="Neighborhood plane size"/>
        <item name="size_depth" type="Int" initializer="3" label="Neighborhood depth size"/>
        <item name="test_flag" type="Enum" initializer="('unrestricted', 'eigenvalues')" label="Choose test"/>
        <item name="output" type="Image" initializer="output=True" role="return" label="Output"/>
    </gui> 
    """
    
    spacing = tensor1.spacing
    origin = tensor1.origin
    direction = tensor1.direction
    M1,S1,N1 = voxel_parameters(tensor1,kwargs['size_plane'],kwargs['size_depth'])
    M2,S2,N2 = voxel_parameters(tensor2,kwargs['size_plane'],kwargs['size_depth'])
    T,s,df = dtiLogTensorTestAS(kwargs['test_flag'],M1.data,M2.data,S1.data.squeeze(),S2.data.squeeze(),N1,N2,spacing,origin,direction)
    output = T
    return output


def dtiLogTensorTestAS(test_flag,M1,M2,S1,S2,N1,N2,spacing,origin,direction):
    """ Computes voxel-wise test statistics for two groups
    Source: Armin Schwatzman, "RANDOM ELLIPSOIDS AND FALSE DISCOVERY RATES: STATISTICS FOR DIFFUSION TENSOR IMAGING" June 2006

    Input:
    test_flag:	'unrestricted': 	H0: M1=M2
                'eiginvalues': 		H0: both groups have the same eigenvalues, with possible different unknown eigenvectors
    M1,M2:		(Z,Y,X,6):	        arrays of mean log tensors for each group
    S1,S2:		(Z,Y,X):		    arrays of the standard deviations for each group
    N1,N2:				            number of subjects in each groups

    Output:
    T:		    (Z,Y,X):		    array of test statistics
    s		    (Z,Y,X):		    standard deviation
	df:				                degrees of freedom of the distribution """

    # Define some constants
    N = N1+N2
    q = 6
    p = 3
    df = []

    # Test statistics
    if test_flag=='unrestricted':
        T = N1*N2/N * np.sum( ((M1 - M2)**2),3 ) # chi2(q)
        df.insert(0,q)
        df.insert(1,q*(N-2))
    elif test_flag=='eigenvalues':
        L1,V1 = spectral_decomposition(Image(data=M1,data_type="vector"))
        L2,V2 = spectral_decomposition(Image(data=M2,data_type="vector"))
        L1 = L1.data
        L2 = L2.data
        T = N1*N2/N * np.sum( (L1 - L2)**2,3 ) # chi2(p)
        df.insert(0,p)
        df.insert(1,q*(N-2))
    else:
        raise medipy.base.Exception("Unknown test flag : %s"%(test_flag,))

    # Variance
    s = ( (N1-1)*S1**2 + (N2-1)*S2**2 )/(N-2)

    # Statistic
    T = df[1]/df[0] * T/(q*(N-2)*s)
    s = np.sqrt(s)

    return Image(data=T,spacing=spacing,origin=origin,direction=direction),Image(data=s,spacing=spacing,origin=origin,direction=direction),df


