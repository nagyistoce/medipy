import numpy
import scipy.ndimage
from vtk import vtkFloatArray

def texture_from_depth(object_3d, image, depth, center_of_mass = None):
    """ Color the vertices of the object using color from  image with
        a depth factor related to the image center of mass.
    """
    
    # TODO : move this out of the function once the migration to VTK 5.4 is done
    import vtk.util.numpy_support
    
    # Compute the center of mass if necessary
    if center_of_mass is None :
        center_of_mass = scipy.ndimage.center_of_mass(image)
    
    # Vertices coordinates, as VTK and numpy arrays
    vtk_array = object_3d.dataset.GetPoints().GetData()
    points = vtk.util.numpy_support.vtk_to_numpy(vtk_array)
    
    # VTK->Numpy coordinates
    points = numpy.fliplr(points)
    
    # World coordinates->voxel coordinates
    points = numpy.subtract(points, image.origin)
    points = numpy.divide(points, image.spacing)
    
    # Scale each point w.r.t depth factor
    points -= center_of_mass
    points *= depth
    points += center_of_mass
        
    # Probe image : must supply list of shape (3,n) instead of (n,3)
    values = image.data[points.T.tolist()]
        
    # Store values in scalars
    values_vtk = vtk.util.numpy_support.numpy_to_vtk(values)
    if object_3d.dataset.GetPointData().GetScalars() is None :
        object_3d.dataset.GetPointData().SetScalars(vtkFloatArray())
    object_3d.dataset.GetPointData().GetScalars().DeepCopy(values_vtk)