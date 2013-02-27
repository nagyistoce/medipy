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

import types

def array_to_itk_matrix(array, flip):
    """ Create an ``itk.Matrix`` matching the contents and type of given array. If
        ``flip`` is ``True``, then the rows and columns of the ``itk.Matrix`` will be 
        flipped with respect to the numpy array to reflect the different
        coordinate order between ITK and numpy.
    """
    
    itk_type = types.dtype_to_itk[array.dtype.type]
    
    matrix_type = itk.Matrix[itk_type, array.shape[-1], array.shape[-2]]
    matrix_bridge = itk.MatrixBridge[matrix_type]
    
    itk_matrix = matrix_bridge.GetMatrixFromArray(
        numpy.flipud(numpy.fliplr(array)) if flip else array)
    
    return itk_matrix

def itk_matrix_to_array(matrix, flip):
    """ Create an ``numpy.ndarray`` matching the contents of given ``itk.Matrix``. If
        ``flip`` is ``True``, then the rows and columns of the array will be 
        flipped with respect to the ITK matrix to reflect the different
        coordinate order between ITK and numpy.
    """
    
    array = itk.MatrixBridge[matrix].GetArrayFromMatrix(matrix)
    if flip :
        array = numpy.fliplr(numpy.flipud(array))
    return array

########################################################
# Convert itk.Image, itk.VectorImage and numpy.ndarray #
########################################################

def array_to_itk_image(array, transferOwnership):
    """ Create an ``itk.Image`` matching the contents and type of given ``array``. If
        ``transferOwnership`` is ``True``, then the ``itk.Image`` will own the data, and the
        array will not. Otherwise, the ``itk.Image`` does not own the image, and the 
        array is unchanged.
    """
    
    itk_type = types.dtype_to_itk[array.dtype.type]
    image_type = itk.Image[itk_type, array.ndim]
    return itk.NumpyBridge[image_type].GetImageFromArray(array, transferOwnership)

def array_to_itk_vector_image(array, transferOwnership):
    """ Create an ``itk.VectorImage`` matching the contents and type of given array.
        If ``transferOwnership`` is ``True``, then the ``itk.Image`` will own the data, and 
        the array will not. Otherwise, the ``itk.Image`` does not own the image, and
        the array is unchanged.
    """
    
    itk_type = types.dtype_to_itk[array.dtype.type]
    image_type = itk.VectorImage[itk_type, array.ndim-1]
    return itk.NumpyBridge[image_type].GetImageFromArray(array, transferOwnership)

def itk_image_to_array(image, transferOwnership):
    """ Create an ``numpy.ndarray`` matching the contents and type of given ``image``. If
        ``transferOwnership`` is ``True``, then the array will own the data, and the
        image will not. Otherwise, the array does not own the image, and the 
        image is unchanged.
    """
    
    return itk.NumpyBridge[image].GetArrayFromImage(image, transferOwnership)

def itk_vector_image_to_array(vector_image, transferOwnership):
    """ Create an ``numpy.ndarray`` matching the contents and type of given 
        ``vector_image``. If ``transferOwnership`` is ``True``, then the array will own the 
        data, and the image will not. Otherwise, the array does not own the 
        image, and the image is unchanged.
    """
    
    return itk.NumpyBridge[vector_image].GetArrayFromImage(vector_image, transferOwnership)

############################################################
# Convert itk.Image, itk.VectorImage and medipy.base.Image #
############################################################

def medipy_image_to_itk_image(image, transferOwnership):
    """ Create an ``itk.Image`` or ``itk.VectorImage`` matching the contents and type 
        of given ``medipy.base.Image``. If ``transferOwnership`` is ``True``, then the 
        ``itk.Image`` will own the data, and the image will not. Otherwise, the 
        ``itk.Image`` does not own the image, and the image is unchanged.
    """
    
    if image.data_type == "scalar" :
        itk_image = array_to_itk_image(image.data, transferOwnership)
    elif image.data_type == "vector": 
        itk_image = array_to_itk_vector_image(image.data, transferOwnership)
    else :
        raise medipy.base.Exception("Unknown image data_type : %s"%image.data_type)

    matrix_type = itk.Matrix[itk.D, len(image.shape), len(image.shape)]    
    matrix_bridge = itk.MatrixBridge[matrix_type]
    itk_direction = matrix_bridge.GetMatrixFromArray(
        numpy.flipud(numpy.fliplr(image.direction)).astype(float))
    itk_image.SetDirection(itk_direction)
    
    origin = [float(x) for x in reversed(image.origin)]
    itk_image.SetOrigin(origin)
    
    spacing = [float(x) for x in reversed(image.spacing)]
    itk_image.SetSpacing(spacing)
#    itk_metadata = dict_to_itk_metadatadictionary(image.metadata)
#    itk_image.SetMetaDataDictionary(itk_metadata)
    return itk_image

def itk_image_to_medipy_image(itk_image, medipy_image, transferOwnership):
    """ Modify a ``medipy.base.Image`` to match the contents and type of given ITK image. If
        ``transferOwnership`` is ``True``, then the image will own the data, and the
        ``itk.Image`` will not. Otherwise, the image does not own the data, and the 
        ``itk.Image`` is unchanged. If ``medipy_image`` is ``None``, then a new image is 
        created. In any case, the medipy Image is returned
    """
    
    if medipy_image is None :
        itk_type = itk.template(itk_image)[1][0]
        dimension = itk.template(itk_image)[1][1]
        medipy_image = medipy.base.Image(dimension*(0,), types.itk_to_dtype[itk_type])
    
    if itk_image.GetNameOfClass() == "Image" :
        if not itk.NumpyBridge[itk_image].IsBufferShared(medipy_image.data, itk_image) :
            medipy_image.data = itk_image_to_array(itk_image, transferOwnership)
        medipy_image.data_type = "scalar"
    elif itk_image.GetNameOfClass() == "VectorImage" :
        if not itk.NumpyBridge[itk_image].IsBufferShared(medipy_image.data, itk_image) :
            medipy_image.data = itk_vector_image_to_array(itk_image, transferOwnership)
        medipy_image.data_type = "vector"
    
    matrix_type = itk.Matrix[itk.D, medipy_image.ndim, medipy_image.ndim]    
    matrix_bridge = itk.MatrixBridge[matrix_type]
    itk_direction = matrix_bridge.GetArrayFromMatrix(itk_image.GetDirection())
    medipy_image.direction = numpy.fliplr(numpy.flipud(itk_direction))
    
    medipy_image.origin =  [x for x in reversed(itk_image.GetOrigin())]
    medipy_image.spacing = [x for x in reversed(itk_image.GetSpacing())]

    return medipy_image

def itk_image_type(medipy_image):
    """ Return the ITK image type corresponding to the given ``medipy.base.Image``
    """
    
    if medipy_image.data_type == "scalar" :
        itk_type = types.dtype_to_itk[medipy_image.dtype.type]
        image_type = itk.Image[itk_type, medipy_image.ndim]
    elif medipy_image.data_type == "vector" :
        itk_type = types.dtype_to_itk[medipy_image.dtype.type]
        image_type = itk.VectorImage[itk_type, medipy_image.ndim-1]
    else :
        raise medipy.base.Exception("Unknown image data_type : %s"%medipy_image.data_type)
    
    return image_type
