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

def replace_labels(image, change_map, get_modified_indices=False):
    """ Replace the labels in ``image`` as specified by ``change_map``. In the
        following code sample, all voxels with value ``1`` will be replaced by
        ``2``, and all voxels with value ``42`` will be replaced by 0: ::
        
        >>> image = medipy.base.Image((2,2), dtype=numpy.int8, value=0)
        >>> image[0,0] = 1
        >>> image[1,1] = 42
        >>> modified_indices = replace_labels(image, {1: 2, 42: 3}, True)
        >>> dict([[key, numpy.transpose(value).tolist()] for key, value in modified_indices.items()])
        {1: [[0, 0]], 42: [[1, 1]]}
        >>> image[0,0]
        2
        >>> image[1,1]
        3
        >>> (image[0,1], image[1,0])
        (0, 0)
        
        If ``get_modified_indices`` is set to ``True``, then a dictionary of 
        modified indices is returned. Otherwise, no value is returned.
        
        .. warning:: 
        
            This function only works for images of integer types.
    """
    
    if not numpy.issubdtype(image.dtype, int) :
        raise medipy.base.Exception("Image must be of integer type")
    
    indices = []
    
    if get_modified_indices :
        for value in change_map :
            indices.append((value, numpy.where(image == image.dtype.type(value))))
    
    itk_image = medipy.itk.medipy_image_to_itk_image(image, False)
    
    filter = itk.ChangeLabelImageFilter[itk_image, itk_image].New(
        Input = itk_image, InPlace = True, ChangeMap = change_map)
    filter()
    
    if get_modified_indices :
        return dict(indices)
