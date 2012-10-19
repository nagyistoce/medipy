##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import logging
import numpy
from vtk import vtkDataArray, vtkImageData, vtkLookupTable
from vtk.util import vtkConstants

from wx import GetTranslation as _

import medipy.base
from medipy.vtk import vtkColorTransferFunctionWithAlpha, vtkEnhancedLookupTable 

import vtk
from vtk import vtkImageExport
from vtk import vtkImageImport

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



# Useful constants for VTK arrays.
VTK_ID_TYPE_SIZE = vtk.vtkIdTypeArray().GetDataTypeSize()
if VTK_ID_TYPE_SIZE == 4:
    ID_TYPE_CODE = numpy.int32
elif VTK_ID_TYPE_SIZE == 8:
    ID_TYPE_CODE = numpy.int64

VTK_LONG_TYPE_SIZE = vtk.vtkLongArray().GetDataTypeSize()
if VTK_LONG_TYPE_SIZE == 4:
    LONG_TYPE_CODE = numpy.int32
    ULONG_TYPE_CODE = numpy.uint32
elif VTK_LONG_TYPE_SIZE == 8:
    LONG_TYPE_CODE = numpy.int64
    ULONG_TYPE_CODE = numpy.uint64


def get_vtk_to_numpy_typemap():
    """Returns the VTK array type to numpy array type mapping."""
    _vtk_np = {vtkConstants.VTK_BIT:numpy.bool,
                vtkConstants.VTK_CHAR:numpy.int8,
                vtkConstants.VTK_UNSIGNED_CHAR:numpy.uint8,
                vtkConstants.VTK_SHORT:numpy.int16,
                vtkConstants.VTK_UNSIGNED_SHORT:numpy.uint16,
                vtkConstants.VTK_INT:numpy.int32,
                vtkConstants.VTK_UNSIGNED_INT:numpy.uint32,
                vtkConstants.VTK_LONG:LONG_TYPE_CODE,
                vtkConstants.VTK_UNSIGNED_LONG:ULONG_TYPE_CODE,
                vtkConstants.VTK_ID_TYPE:ID_TYPE_CODE,
                vtkConstants.VTK_FLOAT:numpy.float32,
                vtkConstants.VTK_DOUBLE:numpy.float64}
    return _vtk_np

def get_numpy_array_type(vtk_array_type):
    """Returns a numpy array typecode given a VTK array type."""
    return get_vtk_to_numpy_typemap()[vtk_array_type]

def vtk_to_numpy_array(vtk_image) :
    """ Build a numpy array from the given vtk image data. 
    """
    vtk_image_shape = vtk_image.GetDimensions()
    vtk_type = vtk_image.GetScalarType()
    extent = vtk_image.GetWholeExtent()
    numComponents = vtk_image.GetNumberOfScalarComponents()
    dim = (extent[5]-extent[4]+1,
           extent[3]-extent[2]+1,
           extent[1]-extent[0]+1)
    assert vtk_type in get_vtk_to_numpy_typemap().keys(), "Unsupported array type %s"%vtk_type
    assert vtk_type != vtkConstants.VTK_BIT, 'Bit arrays are not supported.'

    if (numComponents > 1):
        dim = dim + (numComponents,)

    dtype = get_numpy_array_type(vtk_type)
    array = numpy.zeros(dim,dtype=dtype)
    export = vtkImageExport()
    export.SetInput(vtk_image)
    export.Export(array)

    return array

def vtk_type_size():
    _vtk_size = { vtkConstants.VTK_SIGNED_CHAR:1,
                vtkConstants.VTK_UNSIGNED_CHAR:1,
                vtkConstants.VTK_SHORT:2,
                vtkConstants.VTK_UNSIGNED_SHORT:2,
                vtkConstants.VTK_INT:4,
                vtkConstants.VTK_UNSIGNED_INT:4,
                vtkConstants.VTK_FLOAT:4,
                vtkConstants.VTK_DOUBLE:8 }
    return _vtk_size

def get_vtk_type_size(vtk_array_type):
    """Returns the size of a given VTK array type."""
    return vtk_type_size()[vtk_array_type]

def numpy_array_to_vtk(array):
    numComponents = 1
    dim = array.shape
    if len(dim) == 0:
        dim = (1,1,1)
    elif len(dim) == 1:
        dim = (1, 1, dim[0])
    elif len(dim) == 2:
        dim = (1, dim[0], dim[1])
    elif len(dim) == 4:
        numComponents = dim[3]
        dim = (dim[0],dim[1],dim[2])

    vtk_extent = (0, dim[2]-1,
              0, dim[1]-1,
              0, dim[0]-1)

    vtk_array_type = get_vtk_array_type(array.dtype)

    # Adjust the number of components for complex types
    complexComponents = 1
    if array.dtype in [numpy.complex64, numpy.complex128] :
        number_of_components = number_of_components*2
        complexComponents = 2

    _import = vtkImageImport()
    size = len(array.flat)*get_vtk_type_size(vtk_array_type)*complexComponents
    _import.SetDataScalarType(vtk_array_type)
    _import.SetNumberOfScalarComponents(numComponents)
    _import.SetDataExtent(vtk_extent)
    _import.SetWholeExtent(vtk_extent)
    _import.CopyImportVoidPointer(array, size)
    _import.Update()

    return _import.GetOutput()


def build_vtk_colormap(colormap, vtk_colormap=None) :
    """ Build either a vtkLookupTable or a vtkColorTransferFunctionWithAlpha 
        from the given colormap. The colormap is specified as a  
        custom table -- which must respect the formats of the dictionaries 
        mentionned above. If vtk_colormap is None, a new vtk object is created, 
        otherwise vtk_colormap is used as a container.
    """
    
    if type(colormap) not in [list, tuple] :
        colormap_type = str(type(colormap))
        raise medipy.base.Exception(_("Cannot process colormap of type %s")%(colormap_type))
    
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
