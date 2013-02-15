##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import itk
import medipy.itk
import medipy.base
import numpy as np

def least_squares(images,accu):
    """ Least Square Second Order Symmetric Tensor Estimation.
        A diffusion serie is composed of a float reference image (first element 
        in the list) and a set of float diffusion weighted images (on shell, 
        i.e. one bval).
        
        All images must have the same shape and the same dtype, and must 
        contain diffusion metadata.
        
        <gui>
            <item name="images" type="ImageSerie" label="Input"/>
            <item name="accu" type="Enum" initializer="('First', 'Mean')" label="Accumulation Type" />
            <item name="output" type="Image" initializer="output=True" 
                  role="return" label="Output"/>
        </gui>
    """
    
    # We're in the same package as itkSecondOrderSymmetricTensorReconstructionFilter, 
    # so it has already been included in itk by __init__

    # reference must be at the begining of the diffusion images list
    bval = [image.metadata["mr_diffusion_sequence"][0].diffusion_bvalue for image in images]
    if 0 in bval :
        nb_ref = bval.count(0)
        argb = np.argsort(bval)
        if accu=="First" :
            reference = images[argb[0]]
            images.remove(reference)
            images.insert(0,reference)
        elif accu=="Mean" :
            references = []
            for index in argb[:nb_ref] :
                references.append( images[index] )
            reference = references[0]
            images.remove(reference)
            for r in references[1:] :
                images.remove(r)
                reference.data += r
            reference.data /= nb_ref        
            images.insert(0,reference)
        elif accu=="Overall Mean" :
            if len(bval)%nb_ref==0 :
                # ToDo: track gradient directions if possible (not working yet)
                f = len(bval)/nb_ref
                sequences = []
                for i in range(nb_ref) :
                    sequences.append( images[i*nb_ref:(i+1)*nb_ref] )
                images = sequences[0]
                for accu_images in sequences[1:] :
                    for cnt,image in enumerate(accu_images) :
                        images[cnt].data += image.data
                for cnt in range(len(images)) :
                    images[cnt].data /= nb_ref
            else :
                raise medipy.base.Exception("Impossible to average the whole diffusion serie")

        else :
            raise medipy.base.Exception("Uknown accumlation option")
    else :
        raise medipy.base.Exception("No reference image in diffusion sequence")
    
    PixelType = medipy.itk.dtype_to_itk[images[0].dtype.type]
    Dimension = images[0].ndim
    InputImage = itk.Image[PixelType, Dimension]
    OutputImage = itk.VectorImage[PixelType, Dimension]
    EstimationFilter = itk.SecondOrderSymmetricTensorReconstructionFilter[
        InputImage, OutputImage]
    
    estimation_filter = EstimationFilter.New()
    estimation_filter.SetBVal(images[1].metadata["mr_diffusion_sequence"][0].diffusion_bvalue.value)
    for cnt,image in enumerate(images) :
        itk_image = medipy.itk.medipy_image_to_itk_image(image, False)
        estimation_filter.SetInput(cnt,itk_image)
        
        gradient = image.metadata["mr_diffusion_sequence"][0].\
            diffusion_gradient_direction_sequence.value[0].diffusion_gradient_orientation.value
        estimation_filter.SetGradientDirection(cnt, [float(x) for x in gradient])
    
    itk_output = estimation_filter()[0]
    tensors = medipy.itk.itk_image_to_medipy_image(itk_output,None,True)
    tensors.image_type = "tensor_2"
    return tensors
