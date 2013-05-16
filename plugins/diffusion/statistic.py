##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import itk
import numpy

import medipy.base
import medipy.itk
from utils import (spectral_decomposition, exp_transformation, 
                   log_transformation, get_diffusion_information)

def spatial_voxel_parameter_estimation(tensor, size_plane=3, size_depth=3, mask=None):
    """ Neighborhood-based estimation of the mean and standard deviation of 
        second-order tensors.

        <gui>
            <item name="tensor" type="Image" label="DWI data"/>
            <item name="size_plane" type="Int" initializer="3" 
                  label="Neighborhood plane size"/>
            <item name="size_depth" type="Int" initializer="3" 
                  label="Neighborhood depth size"/>
            <item name="mask" type="Image" 
                  initializer="may_be_empty=True, may_be_empty_checked=True" 
                  label="Mask"/>
            <item name="mean" type="Image" initializer="output=True" 
                  role="return" label="Mean image"/>
            <item name="stdev" type="Image" initializer="output=True" 
                  role="return" label="Standard deviation image"/>
        </gui>
    """
    
    log_tensor = log_transformation(tensor)
    log_tensor_itk = medipy.itk.medipy_image_to_itk_image(log_tensor, False)
    ScalarImage = itk.Image[itk.template(log_tensor_itk)[1]]
    
    mask_itk = None
    MaskImage = ScalarImage
    if mask :
        mask_itk = medipy.itk.medipy_image_to_itk_image(mask, False)
        MaskImage = mask_itk.__class__
    
    estimation_filter = itk.SpatialDWIStatisticsImageFilter[
        log_tensor_itk, log_tensor_itk, ScalarImage, MaskImage].New(
        Input=log_tensor_itk, SizePlane=size_plane, SizeDepth=size_depth)
    if mask :
        estimation_filter.SetMaskImage(mask_itk)
    estimation_filter()
    
    mean_itk = estimation_filter.GetMeanImage()
    mean = medipy.itk.itk_image_to_medipy_image(mean_itk, None, True)
    mean.image_type = "tensor_2"
    mean = exp_transformation(mean)
    
    standard_deviation_itk = estimation_filter.GetStandardDeviationImage()
    standard_deviation = medipy.itk.itk_image_to_medipy_image(standard_deviation_itk, None, True)
    
    return mean, standard_deviation

def bootstrap_parameter_estimation(images, size_plane=3, size_depth=3, 
                                   nb_bootstrap=10, bootstrap_type="spatial") :
    """ Local or Spatial bootstrap parameter estimation.

        <gui>
            <item name="images" type="ImageSerie" label="Input diffusion data"/>
            <item name="size_plane" type="Int" initializer="3" 
                  label="Neighborhood plane size (for spatial bootstrap only)"/>
            <item name="size_depth" type="Int" initializer="3" 
                  label="Neighborhood depth size (for spatial bootstrap only)"/>
            <item name="nb_bootstrap" type="Int" initializer="10" 
                  label="Number of bootstrap samples"/>
            <item name="bootstrap_type" type="Enum" 
                  initializer="('spatial', 'local')" label="Choose bootstrap strategy"/>
            <item name="mean" type="Image" initializer="output=True" 
                  role="return" label="Mean tensor image"/>
            <item name="var" type="Image" initializer="output=True" 
                  role="return" label="Variance image"/>
        </gui>
    """

    ScalarImage = itk.Image[itk.F, images[0].ndim]
    VectorImage = itk.VectorImage[itk.F, images[0].ndim]
    EstimationFilter = itk.BootstrapParameterEstimationImageFilter[ScalarImage, VectorImage]
    
    estimation_filter = EstimationFilter.New(
        BVal=float(get_diffusion_information(images[1])["diffusion_bvalue"]),
        NumberOfBootstrap=nb_bootstrap, SizePlane=size_plane, 
        SizeDepth=size_depth, UseSpatialBootstrap=(bootstrap_type=='spatial'))

    for cnt,image in enumerate(images) :
        image.data = numpy.cast[numpy.single](image.data)
        itk_image = medipy.itk.medipy_image_to_itk_image(image, False)
        grad = get_diffusion_information(image)["diffusion_gradient_orientation"]
        estimation_filter.SetInput(cnt,itk_image)
        estimation_filter.SetGradientDirection(cnt,grad.tolist())

    estimation_filter.Update()

    itk_mean = estimation_filter.GetOutput(0)
    mean = medipy.itk.itk_image_to_medipy_image(itk_mean,None,True)
    mean.image_type = "tensor_2"
    mean = exp_transformation(mean)
    
    itk_var = estimation_filter.GetOutput(1)
    var = medipy.itk.itk_image_to_medipy_image(itk_var,None,True)
    # The variance image must be a VectorImage to be fed to the ITK filter, but
    # is really a scalar image
    var.data = var.data.reshape(images[0].shape)
    var.data_type = "scalar"

    return mean,var

def spatial_voxel_test_gui(tensor1,tensor2, size_plane, size_depth, test_flag) :
    """Multivariate Statistical Tests at a voxel level.

    <gui>
        <item name="tensor1" type="Image" label="First time point tensor image"/>
        <item name="tensor2" type="Image" label="Second time point tensor image"/>
        <item name="size_plane" type="Int" initializer="3" label="Neighborhood plane size"/>
        <item name="size_depth" type="Int" initializer="3" label="Neighborhood depth size"/>
        <item name="test_flag" type="Enum" initializer="('unrestricted', 'eigenvalues')" label="Choose test"/>
        <item name="output" type="Image" initializer="output=True" role="return" label="Output"/>
    </gui> 
    """
    
    spacing = tensor1.spacing
    origin = tensor1.origin
    direction = tensor1.direction
    M1,S1 = spatial_voxel_parameter_estimation(tensor1,size_plane,size_depth)
    M2,S2 = spatial_voxel_parameter_estimation(tensor2,size_plane,size_depth)
    
    N1 = size_plane*size_plane*size_depth
    N2 = N1
    
    T,s,df = _dtiLogTensorTestAS(test_flag, M1, M2, S1, S2, N1, N2)
    
    return T


def _dtiLogTensorTestAS(test_flag, M1, M2, S1, S2, N1, N2):
    """ Computes voxel-wise test statistics for two groups
    Source: Armin Schwatzman, "RANDOM ELLIPSOIDS AND FALSE DISCOVERY RATES: STATISTICS FOR DIFFUSION TENSOR IMAGING" June 2006

    Input:
    test_flag:	'unrestricted': 	H0: M1=M2
                'eiginvalues': 		H0: both groups have the same eigenvalues, with possible different unknown eigenvectors
    M1,M2:		(Z,Y,X,6):	        images of mean tensors for each group
    S1,S2:		(Z,Y,X):		    images of the standard deviations for each group
    N1,N2:				            number of subjects in each groups

    Output:
    T:		    (Z,Y,X):		    image of test statistics
    s		    (Z,Y,X):		    image of standard deviation
	df:				                degrees of freedom of the distribution """

    # Define some constants
    N = N1+N2
    q = 6
    p = 3
    df = []
    
    # Transform the mean images to the Log domain
    log_M1 = log_transformation(M1)
    log_M2 = log_transformation(M2)

    # Test statistics
    if test_flag=='unrestricted':
        T = N1*N2/N * numpy.sum( ((log_M1.data - log_M2.data)**2),3 ) # chi2(q)
        sign = numpy.sign( numpy.sum( (log_M1.data - log_M2.data),3 ) )
        df.insert(0,q)
        df.insert(1,q*(N-2))
    elif test_flag=='eigenvalues':
        L1,V1 = spectral_decomposition(log_M1)
        L2,V2 = spectral_decomposition(log_M2)
        T = N1*N2/N * numpy.sum( (L1.data - L2.data)**2,3 ) # chi2(p)
        sign = numpy.sign( numpy.sum( (L1.data - L2.data),3 ) )
        df.insert(0,p)
        df.insert(1,q*(N-2))
    else:
        raise medipy.base.Exception("Unknown test flag : %s"%(test_flag,))

    # Variance
    s = ( (N1-1)*S1.data**2 + (N2-1)*S2.data**2 )/(N-2)

    # Statistic
    #T = df[1]/df[0] * T/(q*(N-2)*s)
    index = numpy.where(s!=0)
    temp = numpy.zeros(T.shape,dtype=numpy.single)
    temp[index] = df[1]/df[0] * T[index]/(q*(N-2)*s[index])
    T = temp

    s = numpy.sqrt(s)

    T[numpy.isnan(T)] = -1
    T[numpy.isinf(T)] = -1

    T = medipy.base.Image(data=sign*T, 
        spacing=M1.spacing, origin=M1.origin, direction=M1.direction)
    s = medipy.base.Image(data=s, 
        spacing=M1.spacing, origin=M1.origin, direction=M1.direction)

    return T, s, df
