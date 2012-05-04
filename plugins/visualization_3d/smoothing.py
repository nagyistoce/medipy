##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

from vtk import *

def windowed_sinc_smoothing(input, iterations, output):
    """ Use a vtkWindowedSincPolyDataFilter to smooth the input image
    """
    
    smoothing_filter = vtkWindowedSincPolyDataFilter()
    smoothing_filter.SetInput(input.dataset)
    smoothing_filter.SetNumberOfIterations(iterations)
    smoothing_filter.NormalizeCoordinatesOn()
    smoothing_filter.SetFeatureAngle(180)
    smoothing_filter.Update()
    
    output.dataset = smoothing_filter.GetOutput()

def generate_normals(input, splitting, feature_angle, output):
    """ Generate normals
    """
    
    normals_filter = vtkPolyDataNormals()
    normals_filter.SetInput(input.dataset)
    normals_filter.SetSplitting(splitting)
    normals_filter.SetFeatureAngle(feature_angle)
    
    normals_filter.Update()
    output.dataset = normals_filter.GetOutput()

if __name__ == "__main__" :
    from vtk import *
    from medipy.base import Object3D
    
    sphere_source = vtkSphereSource()
    sphere_source.Update()
    
    object_3d = Object3D(sphere_source.GetOutput())
    
    windowed_sinc_smoothing(object_3d, 15, object_3d)