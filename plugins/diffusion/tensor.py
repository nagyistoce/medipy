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

def ls_SecondOrderSymmetricTensorEstimation(images):
    """ Least Square Second Order Symmetric Tensor Estimation.
    A diffusion serie is composed of a float reference image (first element in the list) and a set of float diffusion weighted images (on shell, i.e. one bval).
    """
    
    # We're in the same package as itkSecondOrderSymmetricTensorReconstructionFilter, so it has already been included in itk by __init__
    
    estimation_filter = itk.SecondOrderSymmetricTensorReconstructionFilter[itk.Image[itk.F,3], itk.VectorImage[itk.F,3]].New()
    estimation_filter.SetBVal(images[1].metadata["mr_diffusion_sequence"][0].diffusion_bvalue)
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
    itk_output = estimation_filter()[0]
    tensors = medipy.itk.itk_image_to_medipy_image(itk_output,None,True)

    return tensors


def save_SecondOrderSymmetricTensor(tensors,fout):
    """ Save .vtk tensor image  
    """
    import vtk

    data = dti6to33(tensors.data)
    spacing = tensors.spacing

    # Convert numpy -> VTK table
    vtk_array = vtk.vtkFloatArray()
    vtk_array.SetNumberOfComponents(9)
    vtk_array.SetVoidArray(data, len(data.ravel()), 1)

    # Create VTK dataset
    structured_points = vtk.vtkStructuredPoints()
    structured_points.SetDimensions(data.shape[2],data.shape[1],data.shape[0])
    structured_points.SetSpacing(spacing[2],spacing[1],spacing[0])
    structured_points.GetPointData().SetTensors(vtk_array)

    # Write VTK file
    writer = vtk.vtkStructuredPointsWriter()
    writer.SetFileName(fout)
    writer.SetFileTypeToBinary()
    writer.SetInput(structured_points)
    writer.Update()

    return 0

def dti6to33(dt6):
    """ Full second order symmetric tensor from the six idependent components
    """
    dt33 = np.zeros(dt6.shape[:3]+(3,3),dtype=np.single)
    dt33[:,:,:,0,0] = dt6[:,:,:,0]
    dt33[:,:,:,0,1] = dt6[:,:,:,1]
    dt33[:,:,:,0,2] = dt6[:,:,:,2]
    dt33[:,:,:,1,0] = dt6[:,:,:,1]
    dt33[:,:,:,1,1] = dt6[:,:,:,3]
    dt33[:,:,:,1,2] = dt6[:,:,:,4]
    dt33[:,:,:,2,0] = dt6[:,:,:,2]
    dt33[:,:,:,2,1] = dt6[:,:,:,4]
    dt33[:,:,:,2,2] = dt6[:,:,:,5]

    return dt33
