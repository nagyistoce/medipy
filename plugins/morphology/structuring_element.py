##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import itk

def name_to_structuring_element(shape, dimension, radius) :
    """ Create an ITK structuring element given a name and a shape
    """
    
    StructuringElement = itk.FlatStructuringElement[dimension]
    
    if shape == "ball" :
        return StructuringElement.Ball(dimension*(radius,))
    elif shape == "box" :
        return StructuringElement.Box(dimension*(radius,))
    elif shape == "cross" :
        return StructuringElement.Cross(dimension*(radius,))
    else :
        raise Exception("Unknown structuring element shape : \"{0}\"".format(shape))