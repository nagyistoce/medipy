##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import numpy
import medipy.base

def mean_stdev_normalization(reference, image):
    """ Return a normalized version of image, so that the mean and standard 
        deviation match those of reference.
        
        <gui>
            <item name="reference" type="Image" label="Reference"/>
            <item name="image" type="Image" label="Image"/>
            <item name="output" type="Image" initializer="output=True"
                role="return" label="Output"/>
        </gui>
    """
    
    means = (reference.data.mean(), image.data.mean())
    stdevs = (reference.data.std(), image.data.std())
    
    alpha = stdevs[0]/stdevs[1]
    beta = means[0]-means[1]*alpha
    
    data = alpha*image.data+beta
    output = medipy.base.Image(data=data)
    output.copy_information(image)
    
    return output
