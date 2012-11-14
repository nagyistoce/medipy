##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import numpy

from vtk import vtkGlyph3D, vtkPolyDataMapper, vtkActor, vtkLineSource, vtkImageActor, vtkTensorGlyph, vtkSphereSource

from layer import Layer
 
import medipy.itk
import medipy.vtk
from medipy.diffusion.utils import spectral_decomposition

import vtk

from medipy.diffusion.utils import rotation33todt6

class Tensor2Layer(Layer) :
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

        super(Tensor2Layer, self).__init__(world_to_slice, tensor_image, display_coordinates, colormap, opacity)
        self.add_allowed_event("display_mode")

        self._change_information.Update()
        self.vtk_slice_tensors = self._change_information.GetOutput()

        ######################
        # Initialize members #
        ###################### 

        self._actor = vtk.vtkAssembly() #vtkActor()     
        self._actor_1 = vtkImageActor() 
        self._actor_2 = vtkActor()

        self._glyph_line = vtkGlyph3D()
        self._glyph_tensor = vtkTensorGlyph()

        self._line = vtkLineSource()
        self._sphere = vtkSphereSource()

        self._mapper_line = vtkPolyDataMapper()
        self._mapper_tensor = vtkPolyDataMapper()

        ############################
        # Property-related members #
        ############################

        self._actor_1.InterpolateOff() 
        self._origin = self.vtk_slice_tensors.GetOrigin()
        self._spacing = self.vtk_slice_tensors.GetSpacing()                    

        if display_mode=="principal_direction_voxel" :
            self._actor.AddPart(self._actor_1)
        else :
            self._actor.AddPart(self._actor_2)

        self._glyph_line.ScalingOn()
        self._glyph_line.SetVectorModeToUseVector()
        self._glyph_line.SetScaleModeToScaleByVector()
        self._glyph_line.SetColorModeToColorByVector()
        self._glyph_line.SetScaleFactor(1.0)
        self._glyph_line.SetSource(self._line.GetOutput())
        self._glyph_line.ClampingOff()

        self._sphere.SetThetaResolution(6)
        self._sphere.SetPhiResolution(4)

        self._glyph_tensor.SetScaleFactor(700)
        self._glyph_tensor.ColorGlyphsOn()
        self._glyph_tensor.ExtractEigenvaluesOn()
        self._glyph_tensor.SetColorModeToEigenvalues()
        self._glyph_tensor.ClampScalingOff() 

        self._mapper_line.ScalarVisibilityOn()
        self._mapper_line.SetScalarModeToUsePointFieldData()     
        
        ###################
        # Private members #
        ###################

        self._set_display_mode(display_mode)
        self.add_observer("position", self.on_position)           


    def on_position(self, event) :

        self.numpy_slice_tensors = medipy.vtk.bridge.vtk_image_to_medipy_image(
            self.vtk_slice_tensors, None)
        if self.display_coordinates_=="physical" :          
            #self._reslicer_axes_inverse
            self.numpy_slice_tensors.data = rotation33todt6(self.numpy_slice_tensors.data,self.world_to_slice[::-1,::-1])
        self.numpy_slice_tensors.data_type = "vector"
        self.numpy_slice_tensors.image_type = "tensor_2"

        if self._display_mode == "principal_direction_voxel" :
            numpy_eigenvalues,numpy_eigenvectors = spectral_decomposition(self.numpy_slice_tensors)
            numpy_principal_direction = medipy.base.Image(
                    data=numpy.ascontiguousarray(numpy.cast[numpy.uint8](
                        numpy.abs(numpy_eigenvectors[...,6:])*255.0)),
                    data_type="vector")
            numpy_principal_direction.copy_information(self.numpy_slice_tensors) 
            vtk_principal_direction = medipy.vtk.bridge.array_to_vtk_image(
                numpy_principal_direction.data, True, 
                numpy_principal_direction.data_type)
            self._actor_1.SetInput(vtk_principal_direction)
            if self._display_coordinates=="physical" :
                self._actor_1.SetPosition(self._origin)
                self._actor_1.SetScale(self._spacing)
            else :
                self._actor_1.SetPosition((0,0,0))
                self._actor_1.SetScale((1,1,1))

        elif self._display_mode=="principal_direction_line" :
            numpy_eigenvalues,numpy_eigenvectors = spectral_decomposition(self.numpy_slice_tensors)
            val = numpy.log(numpy.maximum(numpy_eigenvalues[...,2],1e-4)*1e4)
            scale = 1.0/val.max()
            
            numpy_principal_direction = medipy.base.Image(
                data=numpy.ascontiguousarray(
                    numpy_eigenvectors[...,6:]*val.repeat(3).reshape(
                        numpy_eigenvalues.shape+(3,))*scale),
                    data_type="vector")
            numpy_principal_direction.copy_information(self.numpy_slice_tensors)
            
            numpy_principal_diffusion = medipy.base.Image(
                data=numpy.ascontiguousarray(numpy_eigenvalues[...,0]),
                data_type="vector")
            numpy_principal_diffusion.copy_information(self.numpy_slice_tensors)
            
            vtk_principal_direction = medipy.vtk.bridge.array_to_vtk_image(
                numpy_principal_direction.data, True, 
                numpy_principal_direction.data_type)
            vtk_principal_diffusion = medipy.vtk.bridge.array_to_vtk_image(
                numpy_principal_diffusion.data, True,
                numpy_principal_direction.data_type)

            #vtk_principal_direction.GetPointData().SetActiveVectors(vtk_principal_direction.GetPointData().GetScalars().GetName())
            name = vtk_principal_direction.GetPointData().GetScalars().GetName()

            self._glyph_line.SetInput(vtk_principal_direction)
            #self._glyph.SetInputArrayToProcess(0,0,0,0,name) # scalars
            self._glyph_line.SetInputArrayToProcess(1,0,0,0,name) # vectors

            # set up a stripper for faster rendering
            stripper = vtk.vtkStripper()
            stripper.SetInput(self._glyph_line.GetOutput())
              
            # Generate a lookup table for coloring by vector components or magnitude
            lut = vtk.vtkLookupTable()
            lut.SetValueRange(val.min()*scale, 1.0)
            #lut.SetSaturationRange(0.1, 1.0)
            #lut.SetHueRange(0.4,0.6)
            #lut.SetRampToLinear()
            # When using a vector component for coloring
            #lut.SetVectorModeToComponent()
            #lut.SetVectorComponent(0)
            # When using vector magnitude for coloring
            lut.SetVectorModeToMagnitude()
            lut.Build()

            # now set up the mapper
            self._mapper_line.SetInput(stripper.GetOutput())
            self._mapper_line.SetLookupTable(lut)
            self._mapper_line.SelectColorArray(name) 

            self._actor_2.SetMapper(self._mapper_line)

        elif self._display_mode=="ellipsoid" :
            numpy_slice_tensors_ = medipy.base.Image(
                data=numpy.zeros(self.numpy_slice_tensors.shape[:3]+(9,),
                                 dtype=numpy.single),
                data_type="vector")
            numpy_slice_tensors_.copy_information(self.numpy_slice_tensors)
            numpy_slice_tensors_[...,:3] = self.numpy_slice_tensors[...,:3]
            numpy_slice_tensors_[...,4:6] = self.numpy_slice_tensors[...,3:5]
            numpy_slice_tensors_[...,8] = self.numpy_slice_tensors[...,5]
            numpy_slice_tensors_[...,3] = self.numpy_slice_tensors[...,1]
            numpy_slice_tensors_[...,6] = self.numpy_slice_tensors[...,2]
            numpy_slice_tensors_[...,7] = self.numpy_slice_tensors[...,4]
            vtk_tensor = medipy.vtk.bridge.array_to_vtk_image(
                numpy_slice_tensors_.data, True, numpy_slice_tensors_.data_type)
            vtk_tensor.GetPointData().SetActiveTensors(vtk_tensor.GetPointData().GetScalars().GetName())

            self._glyph_tensor.SetInput(vtk_tensor)
            self._glyph_tensor.SetSource(self._sphere.GetOutput())
            
            # set up a stripper for faster rendering
            stripper = vtk.vtkStripper()
            stripper.SetInput(self._glyph_tensor.GetOutput())
              
            # now set up the mapper
            self._mapper_tensor.SetInput(stripper.GetOutput())

            self._actor_2.SetMapper(self._mapper_tensor)
            
    
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
        
            old_display_mode = self._display_mode
            self._display_mode = display_mode

            if self._display_mode=="principal_direction_voxel" :
                self._actor.AddPart(self._actor_1)
                self._actor.RemovePart(self._actor_2)
            elif old_display_mode=="principal_direction_voxel" : 
                self._actor.AddPart(self._actor_2)
                self._actor.RemovePart(self._actor_1)


    display_mode = property(_get_display_mode, _set_display_mode)
    
    #####################
    # Private interface #
    #####################
    
    def _on_vtk_colormap(self, dummy) :
        """ Event handler called when the ID of colormap.vtk_colormap changes.
        """
        
        self._set_colormap(self._colormap)

Layer.derived_classes.append(Tensor2Layer)
