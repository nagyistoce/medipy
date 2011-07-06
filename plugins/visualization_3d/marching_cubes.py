from vtk import vtkDiscreteMarchingCubes, vtkMarchingCubes

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
    
    vtk_image = medipy.vtk.build_vtk_image(input)
    vtk_image.SetSpacing(*reversed(input.spacing))
    marching_cubes = vtkMarchingCubes()
    marching_cubes.SetInput(vtk_image)
    marching_cubes.SetValue(0, isovalue)
    marching_cubes.SetComputeNormals(compute_normals)
    marching_cubes.Update()
    
    output.image = input
    output.dataset = marching_cubes.GetOutput()
    
def discrete_marching_cubes(input, isovalue, output):
    """ Generate a triangle mesh of the surface using the discrete marching cubes 
        algorithm.
        
        <gui>
            <item name="input" type="Image" label="Image"/>
            <item name="isovalue" type="Float" initializer="0.5" label="Isovalue"/>
            <item name="output" type="Object3D" role="output" label="Output"/>
        </gui>
    """
    
    vtk_image = medipy.vtk.build_vtk_image(input)
    marching_cubes = vtkDiscreteMarchingCubes()
    marching_cubes.SetInput(vtk_image)
    marching_cubes.SetValue(0, isovalue)
    marching_cubes.Update()
    
    output.dataset = marching_cubes.GetOutput()
    output.image = input
    