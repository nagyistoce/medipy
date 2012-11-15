##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import numpy
from vtk import (vtkActor, vtkAssembly, vtkGlyph3D, vtkImageActor, 
                 vtkLineSource, vtkPolyDataMapper, vtkSphereSource, 
                 vtkStripper, vtkTensorGlyph)
 
import medipy.itk
import medipy.gui.image
import medipy.vtk
from medipy.diffusion.utils import rotation33todt6, spectral_decomposition

class PrincipalDirectionVoxelPipeline(object):
    def __init__(self):
        self.actor = vtkImageActor()
        self.actor.InterpolateOff()
    
    def update(self, image):
        _, eigenvectors = spectral_decomposition(image)
        numpy_principal_direction = numpy.ascontiguousarray(
            (numpy.abs(eigenvectors[...,6:])*255.0).astype(numpy.uint8))
        vtk_principal_direction = medipy.vtk.bridge.array_to_vtk_image(
            numpy_principal_direction, True, "vector")
        vtk_principal_direction.Modified()
        
        self.actor.SetInput(vtk_principal_direction)
        self.actor.SetPosition(image.origin[::-1])
        self.actor.SetScale(image.spacing[::-1])
        
class PrincipalDirectionLinePipeline(object):
    def __init__(self):
        self._glyph_line = vtkGlyph3D()
        self._mapper_line = vtkPolyDataMapper()
        self.actor = vtkActor()
        
        self._line = vtkLineSource()
        
        self._glyph_line.ScalingOn()
        self._glyph_line.SetVectorModeToUseVector()
        self._glyph_line.SetScaleModeToScaleByVector()
        self._glyph_line.SetColorModeToColorByVector()
        self._glyph_line.SetScaleFactor(1.0)
        self._glyph_line.SetSource(self._line.GetOutput())
        self._glyph_line.ClampingOff()
        
        self._mapper_line.ScalarVisibilityOn()
        self._mapper_line.SetScalarModeToUsePointFieldData()
        
        self.actor.SetMapper(self._mapper_line)
    
    def update(self, image):
        eigenvalues, eigenvectors = spectral_decomposition(image)
        val = numpy.log(numpy.maximum(eigenvalues[...,2],1e-4)*1e4)
        scale = 1.0/val.max()
        
        numpy_principal_direction = numpy.ascontiguousarray(
            eigenvectors[...,6:]*val.repeat(3).reshape(
                eigenvalues.shape+(3,))*scale)
        
        vtk_principal_direction = medipy.vtk.bridge.array_to_vtk_image(
            numpy_principal_direction, True, "vector")

        name = vtk_principal_direction.GetPointData().GetScalars().GetName()

        self._glyph_line.SetInput(vtk_principal_direction)
        self._glyph_line.SetInputArrayToProcess(1,0,0,0,name) # vectors
          
        # Generate a lookup table for coloring by vector components or magnitude
        lut = vtk.vtkLookupTable()
        lut.SetValueRange(val.min()*scale, 1.0)
        lut.SetVectorModeToMagnitude()
        lut.Build()

        # now set up the mapper
        self._mapper_line.SetInput(self._glyph_line.GetOutput())
        self._mapper_line.SetLookupTable(lut)
        self._mapper_line.SelectColorArray(name)
        
        self._glyph_line.Modified() 

class EllipsoidPipeline(object):
    def __init__(self):
        self._glyph_tensor = vtkTensorGlyph()
        self._stripper = vtkStripper()
        self._mapper_tensor = vtkPolyDataMapper()
        self.actor = vtkActor()
        
        self._sphere = vtkSphereSource()
        
        self._sphere.SetThetaResolution(6)
        self._sphere.SetPhiResolution(4)

        self._glyph_tensor.SetScaleFactor(700)
        self._glyph_tensor.ColorGlyphsOn()
        self._glyph_tensor.ExtractEigenvaluesOn()
        self._glyph_tensor.SetColorModeToEigenvalues()
        self._glyph_tensor.SetSource(self._sphere.GetOutput())
        self._glyph_tensor.ClampScalingOff()
        
        self._stripper.SetInput(self._glyph_tensor.GetOutput())
        self._mapper_tensor.SetInput(self._stripper.GetOutput())
        self.actor.SetMapper(self._mapper_tensor)
    
    def update(self, image):
        numpy_slice_tensors_ = numpy.zeros(image.shape+(9,), dtype=numpy.single)
        # TODO : dti6to33 ?
        numpy_slice_tensors_[...,:3] = image[...,:3]
        numpy_slice_tensors_[...,4:6] = image[...,3:5]
        numpy_slice_tensors_[...,8] = image[...,5]
        numpy_slice_tensors_[...,3] = image[...,1]
        numpy_slice_tensors_[...,6] = image[...,2]
        numpy_slice_tensors_[...,7] = image[...,4]
        vtk_tensor = medipy.vtk.bridge.array_to_vtk_image(numpy_slice_tensors_, 
                                                          True, "vector")
        vtk_tensor.GetPointData().SetActiveTensors(
            vtk_tensor.GetPointData().GetScalars().GetName())

        self._glyph_tensor.SetInput(vtk_tensor)
        
        self._glyph_tensor.Modified()

class Tensor2Layer(medipy.gui.image.Layer) :
    """ Layer showing its diffusion data. The image will be positionned in the
        viewport such that its origin (i.e. lower left corner) matches the 
        transformed origin of the layer's image (with an altitude of 0).
    """
    
    @staticmethod
    def can_create(image):
        return (image.data_type == "vector" and
                image.image_type == "tensor_2" and
                image.ndim <= 3)
      
    def __init__(self, world_to_slice, tensor_image, display_coordinates="physical",
                 colormap=None, opacity = 1.0, display_mode="principal_direction_voxel") :
        
        self._display_mode = "principal_direction_voxel"
        self.display_coordinates_ = display_coordinates

        medipy.gui.image.Layer.__init__(self, world_to_slice, tensor_image, 
                                        display_coordinates, colormap, opacity)
        self.add_allowed_event("display_mode")

        ######################
        # Initialize members #
        ######################
        
        self._pipelines = {
            "principal_direction_voxel" : PrincipalDirectionVoxelPipeline(),
            "principal_direction_line" : PrincipalDirectionLinePipeline(),
            "ellipsoid" : EllipsoidPipeline()
        } 

        ############################
        # Property-related members #
        ############################
        
        self._actor = vtkAssembly()

        ###################
        # Private members #
        ###################

        self._set_display_mode(display_mode)
        self.add_observer("any", self._update)
        
    def _update(self, event) :

        # We need the actual output of _change_information
        self._change_information.Update()

        numpy_slice_tensors = medipy.vtk.bridge.vtk_image_to_medipy_image(
            self._change_information.GetOutput(), None)
        if self.display_coordinates_=="physical" :          
            numpy_slice_tensors.data = rotation33todt6(numpy_slice_tensors.data,self.world_to_slice[::-1,::-1])
        numpy_slice_tensors.data_type = "vector"
        numpy_slice_tensors.image_type = "tensor_2"

        self._pipelines[self._display_mode].update(numpy_slice_tensors)

    def _get_actor(self):
        "VTK ImageActor."
        return self._actor

    ##############
    # Properties #
    ##############

    def _get_display_mode(self) :
        return self._display_mode
    
    def _set_display_mode(self, display_mode) :
        if display_mode not in["principal_direction_voxel", "principal_direction_line", "ellipsoid"] :
            raise medipy.base.Exception("Unknown display mode : %s"%(display_mode,))

        if self._display_mode!=display_mode :
            self._display_mode = display_mode
            for pipeline in self._pipelines.values() :
                self._actor.RemovePart(pipeline.actor)
            self._actor.AddPart(self._pipelines[self._display_mode].actor)
        
    display_mode = property(_get_display_mode, _set_display_mode)

medipy.gui.image.Layer.derived_classes.append(Tensor2Layer)
