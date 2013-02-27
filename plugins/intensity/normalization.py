##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import numpy
import medipy.base
import medipy.base.exception
import medipy.logic

def mean_stdev_normalization(reference, image, mask_ref=None, mask_image=None):
    """ Return a normalized version of image, so that the mean and standard 
        deviation match those of reference.
        
        <gui>
            <item name="reference" type="Image" label="Reference"/>
            <item name="image" type="Image" label="Image"/>
            <item name="mask_ref" type="Image" initializer="may_be_empty=True" label="Mask of reference image" />
            <item name="mask_image" type="Image" initializer="may_be_empty=True" label="Mask of the image to normalize" />       
            <item name="output" type="Image" initializer="output=True"
                role="return" label="Output"/>
        </gui>
    """
    
    if mask_ref :
        meanref=reference.data[numpy.nonzero(reference.data)].mean()
        stdref=reference.data[numpy.nonzero(reference.data)].std()
    else :
        meanref=reference.data.mean()
        stdref=reference.data.std()
        
    if mask_image :
        meanimage=image.data[numpy.nonzero(image.data)].mean()
        stdimage=image.data[numpy.nonzero(image.data)].std()
    else :
        meanimage=image.data.mean()
        stdimage=image.data.std()
        
    means = (meanref, meanimage)
    stdevs = ( stdref,stdimage )
    
    alpha = stdevs[0]/stdevs[1]
    beta = means[0]-means[1]*alpha
    
    data = alpha*image.data+beta
    output = medipy.base.Image(data=data)
    output.copy_information(image)
    
    if mask_image :
        output=medipy.logic.apply_mask( output, mask_image.astype(numpy.float32), 0, 0)
    
    
    return output

def one_parameter_linear_regression_normalization(src,ref):
    """ Return a normalized version of image, so that the mean and standard 
        deviation match those of reference.
        
        <gui>
            <item name="src" type="Image" label="Image to normalize"/>
            <item name="ref" type="Image" label="Reference"/>
            <item name="output" type="Image" initializer="output=True"
                role="return" label="Output"/>
        </gui>
    """
    
    if src.shape != ref.shape :
        raise medipy.base.exception.Exception("Images must have the same size") 
         
    alpha=numpy.sum(numpy.multiply(src,ref),dtype=float)/numpy.sum(numpy.multiply(src,src),dtype=float)
    data=alpha*src.data
    output = medipy.base.Image(data=data)
    output.copy_information(src)
    
    return output
