##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import itk
import numpy
import medipy.base

itk_to_dtype_table = [
    # Unsigned int types (char, short, long)
    (itk.UC, numpy.ubyte),
    (itk.US, numpy.ushort),
    (itk.UL, numpy.uint32),
    # Signed int types (char, short, long)
    (itk.SC, numpy.byte),
    (itk.SS, numpy.short),
    (itk.SL, numpy.int32),
    # Float types (float, double, long double)
    (itk.F, numpy.single),
    (itk.D, numpy.double),
    (itk.LD, numpy.longdouble),
    # Complex types
    #(itk.COMPLEX_REALS, ),
    # Vector types
    #(itk.VECTOR_REALS, ),
    # RGB types
    #(itk.RGBS, ),
    # RGBA types
    #(itk.RGBAS, ),
    # Covariant vector types
    #(itk.COV_VECTOR_REALS, ),
]

itk_to_dtype = {}
dtype_to_itk = {}

for itk_type, dtype in itk_to_dtype_table :
    itk_to_dtype[itk_type] = dtype
    dtype_to_itk[dtype] = itk_type

########################################################
# Convert itk.Image, itk.VectorImage and numpy.ndarray #
########################################################

def array_to_itk_image(array, transferOwnership):
    """ Create an itk.Image matching the contents and type of given array. If
        transferOwnership is True, then the itk.Image will own the data, and the
        array will not. Otherwise, the itk.Image does not own the image, and the 
        array is unchanged.
    
        >>> import medipy.itk
        >>> array = numpy.ndarray((40, 128, 128), dtype=numpy.float32)
        >>> image = array_to_itk_image(array, False)
        >>> print image.GetNameOfClass(), image.GetRequestedRegion().GetSize()
        Image itkSize3([128, 128, 40])
    """
    
    itk_type = dtype_to_itk[array.dtype.type]
    image_type = itk.Image[itk_type, array.ndim]
    return itk.NumpyBridge[image_type].GetImageFromArray(array, transferOwnership)

def array_to_itk_vector_image(array, transferOwnership):
    """ Create an itk.VectorImage matching the contents and type of given array.
        If transferOwnership is True, then the itk.Image will own the data, and 
        the array will not. Otherwise, the itk.Image does not own the image, and
        the array is unchanged.
    
        >>> import medipy.itk
        >>> array = numpy.ndarray((40, 128, 128, 6), dtype=numpy.float32)
        >>> image = array_to_itk_vector_image(array, False)
        >>> print image.GetNameOfClass(), image.GetRequestedRegion().GetSize(), image.GetNumberOfComponentsPerPixel()
        VectorImage itkSize3([128, 128, 40]) 6
    """
    
    itk_type = dtype_to_itk[array.dtype.type]
    image_type = itk.VectorImage[itk_type, array.ndim-1]
    return itk.NumpyBridge[image_type].GetImageFromArray(array, transferOwnership)

def itk_image_to_array(image, transferOwnership):
    """ Create an numpy.ndarray matching the contents and type of given image. If
        transferOwnership is True, then the array will own the data, and the
        image will not. Otherwise, the array does not own the image, and the 
        image is unchanged.
    
        >>> import medipy.itk
        >>> image = itk.Image[itk.F, 3].New()
        >>> region = itk.ImageRegion[3]()
        >>> region.SetIndex((0,0,0))
        >>> region.SetSize((128,128,40))
        >>> image.SetRegions(region)
        >>> image.Allocate()
        >>> array = itk_image_to_array(image, False)
        >>> print array.shape, array.dtype
        (40, 128, 128) float32
    """
    
    return itk.NumpyBridge[image].GetArrayFromImage(image, transferOwnership)

def itk_vector_image_to_array(image, transferOwnership):
    """ Create an numpy.ndarray matching the contents and type of given 
        VectorImage. If transferOwnership is True, then the array will own the 
        data, and the image will not. Otherwise, the array does not own the 
        image, and the image is unchanged.
    
        >>> import medipy.itk
        >>> image = itk.VectorImage[itk.F, 3].New()
        >>> image.SetNumberOfComponentsPerPixel(6)
        >>> region = itk.ImageRegion[3]()
        >>> region.SetIndex((0,0,0))
        >>> region.SetSize((128,128,40))
        >>> image.SetRegions(region)
        >>> image.Allocate()
        >>> array = itk_vector_image_to_array(image, False)
        >>> print array.shape, array.dtype
        (40, 128, 128, 6) float32
    """
    
    return itk.NumpyBridge[image].GetArrayFromImage(image, transferOwnership)

############################################################
# Convert itk.Image, itk.VectorImage and medipy.base.Image #
############################################################

def medipy_image_to_itk_image(image, transferOwnership):
    """ Create an itk.Image or itk.VectorImage matching the contents and type 
        of given medipy.base.Image. If transferOwnership is True, then the 
        itk.Image will own the data, and the image will not. Otherwise, the 
        itk.Image does not own the image, and the image is unchanged.
    
        >>> import medipy.base
        >>> import medipy.itk
        >>> medipy_image = medipy.base.Image((40, 128, 128), dtype=numpy.float32, spacing=[1,2,3])
        >>> itk_image = medipy_image_to_itk_image(medipy_image, False)
        >>> print itk_image.GetNameOfClass(), itk_image.GetRequestedRegion().GetSize(), itk_image.GetSpacing()
        Image itkSize3([128, 128, 40]) itkVectorD3([3, 2, 1])
        
        >>> medipy_image = medipy.base.Image((40, 128, 128), dti="tensor_2", dtype=numpy.float32, spacing=[1,2,3])
        >>> itk_image = medipy_image_to_itk_image(medipy_image, False)
        >>> print itk_image.GetNameOfClass(), itk_image.GetRequestedRegion().GetSize(), itk_image.GetNumberOfComponentsPerPixel(), itk_image.GetSpacing()
        VectorImage itkSize3([128, 128, 40]) 6 itkVectorD3([3, 2, 1])
    """
    
    if image.data_type == "scalar" :
        itk_image = array_to_itk_image(image.data, transferOwnership)
        matrix_type = itk.Matrix[itk.D, len(image.shape), len(image.shape)]
    elif image.data_type == "vector": 
        itk_image = array_to_itk_vector_image(image.data, transferOwnership)
        matrix_type = itk.Matrix[itk.D, len(image.shape)-1, len(image.shape)-1]
    else :
        raise Exception("Unknown image data_type : %s"%image.data_type)
    
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
    """ Modify a medipy.base.Image to match the contents and type of given ITK image. If
        transferOwnership is True, then the image will own the data, and the
        itk.Image will not. Otherwise, the image does not own the data, and the 
        itk.Image is unchanged. If medipy_image is None, then a new image is 
        created. In any case, the medipy Image is returned
    
        >>> import medipy.base
        >>> import medipy.itk
        >>> itk_image = itk.Image[itk.F, 3].New()
        >>> region = itk.ImageRegion[3]()
        >>> region.SetIndex((0,0,0))
        >>> region.SetSize((128,128,40))
        >>> itk_image.SetRegions(region)
        >>> itk_image.Allocate()
        >>> medipy_image = medipy.base.Image()
        >>> itk_image_to_medipy_image(itk_image, medipy_image, False)
        >>> print medipy_image.shape, medipy_image.dtype
        (40, 128, 128) float32
        
        >>> itk_vector_image = itk.VectorImage[itk.F, 3].New()
        >>> itk_vector_image.SetNumberOfComponentsPerPixel(6)
        >>> itk_vector_image.SetRegions(region)
        >>> itk_vector_image.Allocate()
        >>> medipy_vector_image = medipy.base.Image()
        >>> itk_image_to_medipy_image(itk_vector_image, medipy_vector_image, False)
        >>> print medipy_vector_image.shape, medipy_vector_image.dtype, medipy_vector_image.data_type
        (40, 128, 128, 6) float32 vector
    """
    
    if medipy_image is None :
        itk_type = itk.template(itk_image)[1][0]
        dimension = itk.template(itk_image)[1][1]
        medipy_image = medipy.base.Image(dimension*(0,), itk_to_dtype[itk_type])
    
    if itk_image.GetNameOfClass() == "Image" :
        if not itk.NumpyBridge[itk_image].IsBufferShared(medipy_image.data, itk_image) :
            medipy_image.data = itk_image_to_array(itk_image, transferOwnership)
        medipy_image.data_type = "scalar"
        matrix_type = itk.Matrix[itk.D, len(medipy_image.shape), len(medipy_image.shape)]
    elif itk_image.GetNameOfClass() == "VectorImage" :
        if not itk.NumpyBridge[itk_image].IsBufferShared(medipy_image.data, itk_image) :
            medipy_image.data = itk_vector_image_to_array(itk_image, transferOwnership)
        medipy_image.data_type = "vector"
        matrix_type = itk.Matrix[itk.D, len(medipy_image.shape)-1, len(medipy_image.shape)-1]    

    matrix_bridge = itk.MatrixBridge[matrix_type]
    itk_direction = matrix_bridge.GetArrayFromMatrix(itk_image.GetDirection())
    medipy_image.direction = numpy.fliplr(numpy.flipud(itk_direction))
    
    medipy_image.origin =  [x for x in reversed(itk_image.GetOrigin())]
    medipy_image.spacing = [x for x in reversed(itk_image.GetSpacing())]

    return medipy_image