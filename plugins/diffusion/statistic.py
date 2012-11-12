##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################


import itk
import medipy.itk
import numpy as np
from medipy.base import Image
from medipy.diffusion.utils import spectral_decomposition, exp_transformation, log_transformation
from medipy.diffusion.tensor import ls_SecondOrderSymmetricTensorEstimation_


def spatial_voxel_parameter_estimation_gui(images,*args,**kwargs) :
    """ Spatial parameter estimation.

    <gui>
        <item name="images" type="ImageSerie" label="Input diffusion data"/>
        <item name="size_plane" type="Int" initializer="3" label="Neighborhood plane size"/>
        <item name="size_depth" type="Int" initializer="3" label="Neighborhood depth size"/>
        <item name="mean" type="Image" initializer="output=True" role="return" label="Mean tensor image"/>
        <item name="var" type="Image" initializer="output=True" role="return" label="Variance image"/>
    </gui>
    """    
    size_plane = kwargs['size_plane']
    size_depth = kwargs['size_depth']

    model = ls_SecondOrderSymmetricTensorEstimation_(images)
    mean,var,N = spatial_voxel_parameter_estimation(model,size_plane,size_depth)
    # zero at the image edges == ellispoid display break
    mean = exp_transformation(mean)
    var.data[np.isnan(var.data)] = 0.0
    var.data[np.where(var.data>1.0)] = 0.0

    return mean,var

def spatial_voxel_parameter_estimation(tensor,w_size_plane=3,w_size_depth=3,mask=None):

    mean = medipy.base.Image(data=np.zeros(tensor.shape+(6,),dtype=np.single),spacing=tensor.spacing,origin=tensor.origin,direction=tensor.direction,data_type="vector")
    var = medipy.base.Image(data=np.zeros(tensor.shape+(1,),dtype=np.single),spacing=tensor.spacing,origin=tensor.origin,direction=tensor.direction,data_type="vector")
    log_tensor = log_transformation(tensor)

    if mask==None :
        mask = medipy.base.Image(data=np.zeros((1,1,1),dtype=np.single))
        medipy.diffusion.dtiParamItk(log_tensor,mean,var,mask,w_size_plane,w_size_depth,False)
    else :
        medipy.diffusion.dtiParamItk(log_tensor,mean,var,mask,w_size_plane,w_size_depth,True)
    var.data = np.sqrt( var.data/6.0 ) # compute standard deviation
    mean.image_type = "tensor_2"

    return mean,var,w_size_plane*w_size_plane*w_size_depth


def bootstrap_parameter_estimation_gui(images,*args,**kwargs) :
    """ Local or Spatial bootstrap parameter estimation.

    <gui>
        <item name="images" type="ImageSerie" label="Input diffusion data"/>
        <item name="size_plane" type="Int" initializer="3" label="Neighborhood plane size (for spatial bootstrap only)"/>
        <item name="size_depth" type="Int" initializer="3" label="Neighborhood depth size (for spatial bootstrap only)"/>
        <item name="nb_bootstrap" type="Int" initializer="10" label="Number of bootstrap samples"/>
        <item name="bootstrap_type" type="Enum" initializer="('spatial', 'local')" label="Choose bootstrap strategy"/>
        <item name="mean" type="Image" initializer="output=True" role="return" label="Mean tensor image"/>
        <item name="var" type="Image" initializer="output=True" role="return" label="Variance image"/>
    </gui>
    """    
    bootstrap_type = kwargs['bootstrap_type']
    size_plane = kwargs['size_plane']
    size_depth = kwargs['size_depth']
    nb_bootstrap = kwargs['nb_bootstrap']

    mean,var = bootstrap_parameter_estimation(images,bootstrap_type,size_plane,size_depth,nb_bootstrap)
    # zero at the image edges == ellispoid display break
    mean = exp_transformation(mean)
    var.data[np.isnan(var.data)] = 0.0
    var.data[np.where(var.data>1.0)] = 0.0

    return mean,var

def bootstrap_parameter_estimation(images,bootstrap_type,size_plane,size_depth,nb_bootstrap) :

    ndim = len(images[0].shape)
    estimation_filter = itk.BootstrapParameterEstimationImageFilter[itk.Image[itk.F,ndim], itk.VectorImage[itk.F,ndim]].New()
    estimation_filter.SetBVal(images[1].metadata["mr_diffusion_sequence"][0].diffusion_bvalue)
    estimation_filter.SetNumberOfBootstrap(nb_bootstrap)
    estimation_filter.SetSizePlane(size_plane)
    estimation_filter.SetSizeDepth(size_depth)
    if bootstrap_type=='spatial' :
        estimation_filter.SetUseSpatialBootstrap(True)
    else :
        estimation_filter.SetUseSpatialBootstrap(False)

    for cnt,image in enumerate(images) :
        image.data = np.cast[np.single](image.data)
        itk_image = medipy.itk.medipy_image_to_itk_image(image, False)
        grad = image.metadata["mr_diffusion_sequence"][0].diffusion_gradient_direction_sequence[0].diffusion_gradient_orientation
        itk_grad = itk.Point[itk.F,3]()
        itk_grad[0] = float(grad[0])
        itk_grad[1] = float(grad[1])
        itk_grad[2] = float(grad[2])
        estimation_filter.SetInput(cnt,itk_image)
        estimation_filter.SetGradientDirection(cnt,itk_grad)

    estimation_filter.Update()

    itk_mean = estimation_filter.GetOutput(0)
    itk_var = estimation_filter.GetOutput(1)
    mean = medipy.itk.itk_image_to_medipy_image(itk_mean,None,True)
    var = medipy.itk.itk_image_to_medipy_image(itk_var,None,True)
    mean.image_type = "tensor_2"

    return mean,var






def spatial_voxel_test_gui(tensor1,tensor2,*args,**kwargs):
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
    M1,S1,N1 = spatial_voxel_parameter_estimation(tensor1,kwargs['size_plane'],kwargs['size_depth'])
    M2,S2,N2 = spatial_voxel_parameter_estimation(tensor2,kwargs['size_plane'],kwargs['size_depth'])
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
    #T = df[1]/df[0] * T/(q*(N-2)*s)
    index = np.where(s!=0)
    temp = np.zeros(T.shape,dtype=np.single)
    temp[index] = df[1]/df[0] * T[index]/(q*(N-2)*s[index])
    T = temp

    s = np.sqrt(s)

    T[np.isnan(T)] = -1
    T[np.isinf(T)] = -1

    return Image(data=T,spacing=spacing,origin=origin,direction=direction),Image(data=s,spacing=spacing,origin=origin,direction=direction),df


