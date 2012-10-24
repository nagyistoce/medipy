##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import fnmatch
import logging
import os

import vtk
import numpy

from medipy.diffusion.utils import dti6to33,rotation33todt6

import medipy.base
from medipy.base import ObservableList
from io_base import IOBase

class Vtk(IOBase) :
    """I/O class for Vtk format.
    
    Uses vtk library."""
    
    filenames = ["*.vtk"]
    
    def __init__(self, *args, **kwargs):
        IOBase.__init__(self, *args, **kwargs)
        self._filename = str(self.filename)
    
    def can_load(self) :
        loadable = False
        #for pattern in self.filenames :
        #    if fnmatch.fnmatch(self.filename, pattern):
        #        loadable = True
        #        break
        return loadable
    
    def number_of_images(self) :
        return 1
    
    def load_data(self, index=0) :
        return None
    
    def can_save(self, image) :
        found = False
        for pattern in self.filenames :
            if fnmatch.fnmatch(self.filename, pattern):
                found = True
                break
        return found
    
    def save(self, image) :
        if image.image_type == "tensor_2" :
            data = dti6to33(self.to_axis_aligned_lps_space(image))
            spacing = image.spacing
            origin = image.origin

            # Convert numpy -> VTK table
            vtk_array = vtk.vtkFloatArray()
            vtk_array.SetNumberOfComponents(9)
            vtk_array.SetVoidArray(data, len(data.ravel()), 1)

            # Create VTK dataset
            structured_points = vtk.vtkStructuredPoints()
            structured_points.SetDimensions(data.shape[2],data.shape[1],data.shape[0])
            structured_points.SetSpacing(spacing[2],spacing[1],spacing[0])
            structured_points.SetOrigin(origin[2],origin[1],origin[0])
            structured_points.GetPointData().SetTensors(vtk_array)

            # Write VTK file
            writer = vtk.vtkStructuredPointsWriter()
            writer.SetFileName(self._filename)
            writer.SetFileTypeToBinary()
            writer.SetInput(structured_points)
            writer.Update()
        else :
            raise medipy.base.Exception("Can only save tensor 2 image in .vtk format")    

    def to_axis_aligned_lps_space(self, image):
        """ Transform a tensor image to the closest axis-aligned approximation of LPS (i.e. DICOM) space
        """
        
        # Transform a point from the (R,A,S) NIfTI frame to the (L,P,S) DICOM frame
        if image.ndim == 3 and image._get_number_of_components() == 6 :
            original_direction = image.direction #image.direction[1:,1:]
            direction = numpy.linalg.inv(original_direction)
        else :
            raise medipy.base.Exception("Expect a 4D diffusion tensor image")     
        
        return rotation33todt6(image.data,direction)

