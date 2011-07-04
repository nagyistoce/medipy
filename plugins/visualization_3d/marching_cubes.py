from vtk import vtkDiscreteMarchingCubes, vtkMarchingCubes

from medipy.base import Object3D
import medipy.vtk
from medipy.gui import get_colormap_from_name

def marching_cubes(input, isovalue, compute_normals, output):
    """ Generate a triangle mesh of the surface using the marching cubes 
        algorithm.
        
        :gui:
            input : Image
                Image
            isovalue : Float : 0.5
                Isovalue
            compute_normals : Bool : False
                Compute normals
            output : Object3D
                Output
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
        
        :gui:
            input : Image
                Image
            isovalue : Float
                Isovalue
            output : Object3D
                Output
    """
    
    vtk_image = medipy.vtk.build_vtk_image(input)
    marching_cubes = vtkDiscreteMarchingCubes()
    marching_cubes.SetInput(vtk_image)
    marching_cubes.SetValue(0, isovalue)
    marching_cubes.Update()
    
    output.dataset = marching_cubes.GetOutput()
    output.image = input
    
    # Try to color the object according to colormap
    if "colormap" in input.metadata :
        colormap_name = input.metadata["colormap"]
        colormap = None
        try :
            colormap = get_colormap_from_name(colormap_name)
        except Exception, e :
            pass
        else :
            if len(colormap[0]) == 2 :
                # Stage colormap only
                i = 0
                while not (colormap[i][0][0] <= isovalue < colormap[i][0][1]) and i<len(colormap) :
                    i = i+1
                if i<len(colormap) :
                    color = colormap[i][1]
                    output.diffuse_color = color[:3]
                    output.opacity = color[3]