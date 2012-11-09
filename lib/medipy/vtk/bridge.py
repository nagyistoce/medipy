##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import logging
import operator

import numpy

from vtk import vtkDataArray, vtkImageData, vtkImageExport, vtkImageImport, vtkLookupTable
import vtk.util.numpy_support
from vtk.util import vtkConstants

import medipy.base
from medipy.vtk import vtkColorTransferFunctionWithAlpha, vtkEnhancedLookupTable 

def get_vtk_array_type(numeric_array_type) :
    """Return a VTK typecode from a numpy array."""
    
    # Mapping from numpy array types to VTK array types
    numpy_to_vtk = {
        numpy.dtype(numpy.character):vtkConstants.VTK_UNSIGNED_CHAR,
        numpy.dtype(numpy.uint8):vtkConstants.VTK_UNSIGNED_CHAR,
        numpy.dtype(numpy.uint16):vtkConstants.VTK_UNSIGNED_SHORT,
#       numpy.dtype(ULONG_TYPE_CODE):vtkConstants.VTK_UNSIGNED_LONG,
        numpy.dtype(numpy.int8):vtkConstants.VTK_CHAR,
        numpy.dtype(numpy.int16):vtkConstants.VTK_SHORT,
        numpy.dtype(numpy.int32):vtkConstants.VTK_INT,
#       numpy.dtype(LONG_TYPE_CODE):vtkConstants.VTK_LONG,
        numpy.dtype(numpy.float32):vtkConstants.VTK_FLOAT,
        numpy.dtype(numpy.float64):vtkConstants.VTK_DOUBLE,
        numpy.dtype(numpy.complex64):vtkConstants.VTK_FLOAT,
        numpy.dtype(numpy.complex128):vtkConstants.VTK_DOUBLE
    }
    
    try :
        return numpy_to_vtk[numeric_array_type]
    except KeyError :
        for key in numpy_to_vtk :
            if numpy.issubdtype(numeric_array_type, key):
                return numpy_to_vtk[key]
    raise TypeError, "Couldn't translate array's type %s to VTK"%(numeric_array_type)

def build_vtk_image(array, vtk_image=None, save=1) :
    """ Build a vtkImageData from the given array.
        If vtk_image is None, a new vtkImageData is created, otherwise vtk_image
        is used as a container. 
    """
    
    numpy_array = numpy.asarray(array)
    
    number_of_components = 1
    shape = numpy_array.shape
    if len(numpy_array.shape) == 2 :
        shape = [1, shape[0], shape[1]]
    if len(numpy_array.shape)==4 :
        number_of_components = numpy_array.shape[3]
        shape = shape[:3]
    
    vtk_array_type = get_vtk_array_type(numpy_array.dtype)
    
    # Adjust the number of components for complex types
    if numpy_array.dtype in [numpy.complex64, numpy.complex128] :
        number_of_components = number_of_components*2

    if vtk_image is None :
        vtk_image = vtkImageData()
    
    size = reduce(lambda x,y:x*y, shape)
    extent = (0, shape[2]-1,
              0, shape[1]-1,
              0, shape[0]-1)
    
    vtk_image.SetScalarType(vtk_array_type)
    vtk_image.SetNumberOfScalarComponents(number_of_components)
    vtk_image.SetDimensions(shape[2], shape[1], shape[0])
    
    scalars = vtkDataArray.CreateDataArray(vtk_array_type)
    scalars.SetNumberOfComponents(1)
    scalars.SetVoidArray(numpy_array, size, save)
    vtk_image.GetPointData().SetScalars(scalars) 
    
    vtk_image.SetWholeExtent(extent)
    vtk_image.SetUpdateExtentToWholeExtent()
    vtk_image.SetExtent(extent)
    
    # If we have an image, use its origin and spacing
    if hasattr(array, "spacing") :
        spacing = list(reversed(array.spacing))
        if len(spacing) == 2 :
            spacing.append(1)
        vtk_image.SetSpacing(spacing)
    # 
    if hasattr(array, "origin") :
        origin = list(reversed(array.origin))
        if len(origin) == 2 :
            origin.append(0)
        vtk_image.SetOrigin(origin)
    
    return vtk_image

def array_to_vtk_image(array, copy_data, data_type="scalar"):
    """ Create an vtkImage matching the contents and type of given array. If
        copy_data is True, then the data of the will be copied. Otherwise the
        data will be shared, and the array MUST NOT be destroyed before the 
        vtkImage. data_type specifies how the array should be interpreted : 
        either as a n-array of scalars (data_type="scalar") or as an n-1 array
        of vectors (data_type="vector").
    """
    
    if data_type not in ["scalar", "vector"] :
        raise medipy.base.Exception("Unknown data_type: {0}".format(repr(data_type)))
    
    if data_type == "scalar" :
        ndim = array.ndim
    elif data_type == "vector" :
        ndim = array.ndim-1
    
    if ndim > 3 :
        raise medipy.base.Exception(
            "Cannot convert a {0} array of dimension {1}".format(data_type, 
                                                                 array.ndim))
    
    importer = vtkImageImport()
    
    if numpy.iscomplexobj(array) :
        # Get the first element of the array
        element = array.flat.next()
        scalar_type = vtk.util.numpy_support.get_vtk_array_type(element.real.dtype)
    else :
        scalar_type = vtk.util.numpy_support.get_vtk_array_type(array.dtype)
    importer.SetDataScalarType(scalar_type)
    
    if data_type == "scalar" :
        number_of_components = 1
    elif data_type == "vector" :
        number_of_components = array.shape[ndim]
    if numpy.iscomplexobj(array) :
        number_of_components *= 2
    importer.SetNumberOfScalarComponents(number_of_components)
    
    extent = 6*[0]
    extent[1] = array.shape[ndim-1]-1
    if ndim >= 2 :
        extent[3] = array.shape[ndim-2]-1
    if ndim >= 3 :
        extent[5] = array.shape[ndim-3]-1
    importer.SetDataExtent(extent)
    importer.SetWholeExtent(extent)
    
    size = array.itemsize*reduce(operator.mul, array.shape, 1)
    if copy_data :
        importer.CopyImportVoidPointer(array, size)
    else :
        importer.SetImportVoidPointer(array, size)
    
    importer.Update()

    return importer.GetOutput()

def vtk_image_to_array(vtk_image) :
    """ Create an numpy.ndarray matching the contents and type of given image. 
        If the number of scalars components in the image is greater than 1, then
        the ndarray will be 4D, otherwise it will be 3D. 
    """
    
    exporter = vtkImageExport()
    exporter.SetInput(vtk_image)
    
    # Create the destination array
    extent = vtk_image.GetWholeExtent()
    shape = [extent[5]-extent[4]+1,
             extent[3]-extent[2]+1,
             extent[1]-extent[0]+1]
    if vtk_image.GetNumberOfScalarComponents() > 1:
        shape += [vtk_image.GetNumberOfScalarComponents()]
    dtype = vtk.util.numpy_support.get_numpy_array_type(vtk_image.GetScalarType())
    array = numpy.zeros(shape, dtype=dtype)
    
    exporter.Export(array)
    
    return array

def vtk_image_to_medipy_image(vtk_image, medipy_image):
    if medipy_image is None :
        medipy_image = medipy.base.Image(
            (0,0,0), 
            vtk.util.numpy_support.get_numpy_array_type(vtk_image.GetScalarType()))
    
    medipy_image.data = vtk_image_to_array(vtk_image)
    if vtk_image.GetNumberOfScalarComponents() > 1 :
        medipy_image.data_type = "vector"
    else :
        medipy_image.data_type = "scalar"
    medipy_image.origin = [x for x in reversed(vtk_image.GetOrigin())]
    medipy_image.spacing = [x for x in reversed(vtk_image.GetSpacing())]
    # VTK images are not oriented, assume identity
    medipy_image.direction = numpy.identity(3)
    
    return medipy_image

def build_vtk_colormap(colormap, vtk_colormap=None) :
    """ Build either a vtkLookupTable or a vtkColorTransferFunctionWithAlpha 
        from the given colormap. The colormap is specified as a  
        custom table -- which must respect the formats of the dictionaries 
        mentionned above. If vtk_colormap is None, a new vtk object is created, 
        otherwise vtk_colormap is used as a container.
    """
    
    if type(colormap) not in [list, tuple] :
        colormap_type = str(type(colormap))
        raise medipy.base.Exception("Cannot process colormap of type %s"%(colormap_type))
    
    if type(colormap[0][0]) in [list, tuple] :
        # Stage colormap
        vtk_colormap = build_color_transfer_function_with_alpha(colormap, vtk_colormap)
    else :
        # "Regular" colormap 
        vtk_colormap = build_lookup_table(colormap, vtk_colormap)
    
    return vtk_colormap

def build_lookup_table(colormap, vtk_colormap=None):
    """ Build a vtkLookupTable from a colormap, given as an array of colors.
    """
    
    if vtk_colormap is None :
        vtk_colormap = vtkEnhancedLookupTable()
        vtk_colormap.SetRampToLinear()
    # Allocate a new colormap to avoid numerous ModifiedEvents to be fired by
    # the original colormap
    new_colormap = vtkEnhancedLookupTable()
    new_colormap.DeepCopy(vtk_colormap)
    
    new_colormap.Allocate(len(colormap), len(colormap))
    for i in range(len(colormap)) :
        new_colormap.SetTableValue(i, colormap[i])
    
    vtk_colormap.DeepCopy(new_colormap)
    vtk_colormap.Modified()
    return vtk_colormap

def build_color_transfer_function_with_alpha(colormap, vtk_colormap=None):
    """ Build a vtkColorTransferFunctionWithAlpha from an array of (range, color)
    """
    
    if vtk_colormap is None :
        vtk_colormap = vtkColorTransferFunctionWithAlpha()
    vtk_colormap.RemoveAllPoints()
    
    special_points = []
    for extent, (r,g,b,a) in colormap:
        if extent[0] == extent[1]:
            special_points.append((extent, (r,g,b,a)))
        else:
            try :
                vtk_colormap.AddRGBASegment(extent[0], r, g, b, a, 
                                            extent[1], r, g, b, a)
            except :
                logging.debug("%s %s"%(extent, (r,g,b,a)))
    
    for extent, color in special_points: 
        vtk_colormap.AddRGBAPoint(extent[0], r, g, b, a)

    return vtk_colormap
