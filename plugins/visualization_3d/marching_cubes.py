##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import re

import numpy
from vtk import (vtkDiscreteMarchingCubes, vtkMarchingCubes, vtkMatrix4x4, 
                 vtkTransform, vtkTransformPolyDataFilter)

from medipy.base import Object3D
import medipy.vtk

def marching_cubes(image, isovalue, compute_normals):
    """ Generate a triangle mesh of the surface using the marching cubes 
        algorithm.
        
        <gui>
            <item name="image" type="Image" label="Image"/>
            <item name="isovalue" type="Float" initializer="0.5" label="Isovalue"/>
            <item name="compute_normals" type="Bool" initializer="False" label="Compute normals"/>
            <item name="output" type="Object3D" role="return" label="Output"/>
        </gui>
    """
    
    return _marching_cubes(image, isovalue, compute_normals, vtkMarchingCubes)
    
def discrete_marching_cubes(image, isovalue, compute_normals):
    """ Generate a triangle mesh of the surface using the discrete marching cubes 
        algorithm.
        
        <gui>
            <item name="image" type="Image" label="Image"/>
            <item name="isovalue" type="Float" initializer="0.5" label="Isovalue"/>
            <item name="compute_normals" type="Bool" initializer="False" label="Compute normals"/>
            <item name="output" type="Object3D" role="return" label="Output"/>
        </gui>
    """

    return _marching_cubes(image, isovalue, compute_normals, vtkDiscreteMarchingCubes)

def _marching_cubes(image, isovalue, compute_normals, Filter) :
    """ Create an Object3D using vtkMarchingCubes or one of its derived class.
    """
    
    # Run the VTK filter
    vtk_image = medipy.vtk.array_to_vtk_image(image.data, False, image.data_type)
    vtk_filter = Filter()
    vtk_filter.SetInput(vtk_image)
    vtk_filter.SetValue(0, isovalue)
    vtk_filter.SetComputeNormals(compute_normals)
    vtk_filter.Update()

    # Marching cubes result is in index space, convert to physical space
    polydata = _to_physical_space(image, vtk_filter.GetOutput())
    
    # Use class name (without VTK prefix) for Object3D name
    name = " ".join(re.split(r"([A-Z][a-z]+)", vtk_filter.GetClassName())[1::2])

    return Object3D(polydata, name, image)

def _to_physical_space(image, polydata) :
    """ Transform a ``polydata`` from index space to physical space using the 
        transformation in ``image``.
    """
    
    matrix = numpy.identity(4)
    matrix[:3,:3] = numpy.dot(numpy.diag(image.spacing), image.direction)[::-1,::-1]
    matrix[:3,3] = image.origin[::-1]
    
    transform = vtkTransform()
    transform.SetMatrix(matrix.ravel())

    transform_filter = vtkTransformPolyDataFilter()
    transform_filter.SetTransform(transform)
    transform_filter.SetInput(polydata)
    
    transform_filter.Update()
    
    return transform_filter.GetOutput()
