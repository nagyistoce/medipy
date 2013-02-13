##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

from vtk import vtkDiscreteMarchingCubes, vtkMarchingCubes, vtkMatrix4x4, vtkTransform, vtkTransformPolyDataFilter

from medipy.base import Object3D
import medipy.vtk
from medipy.gui import get_colormap_from_name

def marching_cubes(input, isovalue, compute_normals, output):
    """ Generate a triangle mesh of the surface using the marching cubes 
        algorithm.
        
        <gui>
            <item name="input" type="Image" label="Image"/>
            <item name="isovalue" type="Float" initializer="0.5" label="Isovalue"/>
            <item name="compute_normals" type="Bool" initializer="False" label="Compute normals"/>
            <item name="output" type="Object3D" role="output" label="Output"/>
        </gui>
    """
    
    vtk_image = medipy.vtk.array_to_vtk_image(input.data, False, input.data_type)
    marching_cubes = vtkMarchingCubes()
    marching_cubes.SetInput(vtk_image)
    marching_cubes.SetValue(0, isovalue)
    marching_cubes.SetComputeNormals(compute_normals)
    marching_cubes.Update()
    
    output.image = input
    output.dataset = marching_cubes.GetOutput()
    
def discrete_marching_cubes(input,*args,**kwargs):
    """ Generate a triangle mesh of the surface using the discrete marching cubes 
        algorithm.
        
        <gui>
            <item name="input" type="Image" label="Image"/>
            <item name="isovalue" type="Float" initializer="0.5" label="Isovalue"/>
            <item name="output" type="Object3D" role="return" label="Output"/>
        </gui>
    """

    isovalue = kwargs['isovalue']

    vtk_image = medipy.vtk.array_to_vtk_image(input.data, False, input.data_type)
    marching_cubes = vtkDiscreteMarchingCubes()
    marching_cubes.SetInput(vtk_image)
    marching_cubes.SetValue(0, isovalue)
    marching_cubes.Update()

    vtk_polydata = marching_cubes.GetOutput()

    matrix = vtkMatrix4x4()
    matrix.Identity()
    direction = input.direction[::-1,::-1]
    for i in range(3):
        for j in range(3):
            matrix.SetElement(i,j,direction[i,j])

    transform = vtkTransform()
    transform.Translate(input.origin[::-1])
    transform.Concatenate(matrix)
    transform.Scale(input.spacing[::-1])

    transformPD = vtkTransformPolyDataFilter()
    transformPD.SetTransform(transform)
    transformPD.SetInput(vtk_polydata)

    output = Object3D(transformPD.GetOutput(),"Discrete Marching Cubes")
    output._set_image(input)

    return output

    
