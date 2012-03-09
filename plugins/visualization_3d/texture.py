##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import logging
import sys
import numpy
import scipy.ndimage

from vtk import vtkFloatArray
from vtk import vtkVersion
if vtkVersion.GetVTKMajorVersion() == 5 and vtkVersion.GetVTKMinorVersion() < 2 :
    logging.warning("This plugin may not work with {0}".format(vtkVersion.GetVTKVersion()))
else :
    import vtk.util.numpy_support

def texture_from_depth(object_3d, image, depth, center_of_mass = None):
    """ Color the vertices of the object using color from  image with
        a depth factor related to the image center of mass.
    """
    
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
    
    # Clip to image boundary
    numpy.clip(points, image.ndim*[0], numpy.subtract(image.shape, 1), points)
    
    # Probe image : must supply list of shape (3,n) instead of (n,3)
    values = image.data[points.T.tolist()]
        
    # Store values in scalars
    values_vtk = vtk.util.numpy_support.numpy_to_vtk(values)
    if object_3d.dataset.GetPointData().GetScalars() is None :
        object_3d.dataset.GetPointData().SetScalars(vtkFloatArray())
    object_3d.dataset.GetPointData().GetScalars().DeepCopy(values_vtk)