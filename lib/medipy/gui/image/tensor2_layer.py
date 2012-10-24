##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
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
      
    def __init__(self, world_to_slice, tensor_image, display_coordinates="physical",
                 colormap=None, opacity = 1.0, display_mode="principal_direction_voxel") :
        

        self.display_mode = display_mode
        self.display_coordinates_ = display_coordinates

        super(Tensor2Layer, self).__init__(world_to_slice, tensor_image, display_coordinates, colormap, opacity)

        self._change_information.Update()
        self.vtk_slice_tensors = self._change_information.GetOutput()

        if self.display_mode=="principal_direction_voxel" :
            ############################
            # Property-related members #
            ############################
            
            self._actor = vtkImageActor() 

            self._actor.InterpolateOff()         
            
            ######################
            # Initialize members #
            ######################            
            
            self.add_observer("position", self.on_position)
            

        elif self.display_mode=="principal_direction_line" :

            ############################
            # Property-related members #
            ############################
            
            self._actor = vtkActor()

            ###################
            # Private members #
            ###################

            self._glyph = vtkGlyph3D()
            self._line = vtkLineSource()
            self._mapper = vtkPolyDataMapper()

            self._glyph.ScalingOn()
            self._glyph.SetVectorModeToUseVector()
            self._glyph.SetScaleModeToScaleByVector()
            self._glyph.SetColorModeToColorByVector()
            self._glyph.SetScaleFactor(1.0)
            self._glyph.SetSource(self._line.GetOutput())
            self._glyph.ClampingOff()

            ######################
            # Initialize members #
            ######################
                       
            self.add_observer("position", self.on_position)
            # set up the actor
            self._actor.SetMapper(self._mapper)
            

        elif self.display_mode=="ellipsoid" :

            ############################
            # Property-related members #
            ############################
            
            self._actor = vtkActor()

            ###################
            # Private members #
            ###################

            self._glyph = vtkTensorGlyph()
            self._sphere = vtkSphereSource()
            self._mapper = vtkPolyDataMapper()

            self._sphere.SetThetaResolution(8)
            self._sphere.SetPhiResolution(8)

            self. _glyph.SetScaleFactor(700)
            self._glyph.ColorGlyphsOn()
            self._glyph.ExtractEigenvaluesOn()
            self._glyph.SetColorModeToEigenvalues()
            self._glyph.ClampScalingOn()

            ######################
            # Initialize members #
            ######################

            self.add_observer("position", self.on_position)
            # set up the actor
            self._actor.SetMapper(self._mapper)

        else :

            raise medipy.base.Exception("Unknown display mode : %s"%(self.display_mode,))



    def on_position(self, event) :

        self.numpy_slice_tensors = medipy.vtk.bridge.vtk_to_numpy_array(self.vtk_slice_tensors)
        if self.display_coordinates_=="physical" :          
            #self._reslicer_axes_inverse
            self.numpy_slice_tensors.data = rotation33todt6(self.numpy_slice_tensors.data,self.world_to_slice[::-1,::-1])
        self.numpy_slice_tensors.data_type = "vector"
        self.numpy_slice_tensors.image_type = "tensor_2"

        if self.display_mode == "principal_direction_voxel" :
            numpy_eigenvalues,numpy_eigenvectors = spectral_decomposition(self.numpy_slice_tensors)
            numpy_principal_direction = medipy.base.Image(data=numpy.ascontiguousarray(numpy.cast[numpy.uint8](numpy_eigenvectors[...,::3]*255.0)))
            numpy_principal_direction.copy_information(self.numpy_slice_tensors) 
            #print numpy_principal_direction.spacing, numpy_principal_direction.origin, numpy_principal_direction.shape
            vtk_principal_direction = medipy.vtk.bridge.numpy_array_to_vtk(numpy_principal_direction)
            #print vtk_principal_direction.GetDimensions(), vtk_principal_direction.GetSpacing(), vtk_principal_direction.GetOrigin()
            self._actor.SetInput(vtk_principal_direction)

        elif self.display_mode=="principal_direction_line" :
            numpy_eigenvalues,numpy_eigenvectors = spectral_decomposition(self.numpy_slice_tensors)
            val = numpy.abs(numpy.log(numpy.maximum(numpy_eigenvalues[...,0],1e-6)))
            scale = 1.0/val.max()
            print scale,val.max()
            numpy_principal_direction = medipy.base.Image(data=numpy.ascontiguousarray(numpy_eigenvectors[...,::3]*\
                                                          val.repeat(3).reshape(numpy_eigenvalues.shape+(3,))*scale))
            numpy_principal_direction.copy_information(self.numpy_slice_tensors)
            numpy_principal_diffusion = medipy.base.Image(data=numpy.ascontiguousarray(numpy_eigenvalues[...,0]))
            numpy_principal_diffusion.copy_information(self.numpy_slice_tensors)
            vtk_principal_direction = medipy.vtk.bridge.numpy_array_to_vtk(numpy_principal_direction)
            vtk_principal_diffusion = medipy.vtk.bridge.numpy_array_to_vtk(numpy_principal_diffusion)

            #vtk_principal_direction.GetPointData().SetActiveVectors(vtk_principal_direction.GetPointData().GetScalars().GetName())
            name = vtk_principal_direction.GetPointData().GetScalars().GetName()

            self._glyph.SetInput(vtk_principal_direction)
            #self._glyph.SetInputArrayToProcess(0,0,0,0,name) # scalars
            self._glyph.SetInputArrayToProcess(1,0,0,0,name) # vectors

            # set up a stripper for faster rendering
            stripper = vtk.vtkStripper()
            stripper.SetInput(self._glyph.GetOutput())
              
            # Generate a lookup table for coloring by vector components or magnitude
            lut = vtk.vtkLookupTable()
            lut.SetValueRange(0.0, 1.0)
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
            self._mapper.SetInput(stripper.GetOutput())
            self._mapper.SetLookupTable(lut)
            self._mapper.ScalarVisibilityOn()
            self._mapper.SetScalarModeToUsePointFieldData()
            self._mapper.SelectColorArray(name) 

        elif self.display_mode=="ellipsoid" :
            numpy_slice_tensors_ = medipy.base.Image(data=numpy.zeros(self.numpy_slice_tensors.shape[:3]+(9,),dtype=numpy.single))
            numpy_slice_tensors_.copy_information(self.numpy_slice_tensors)
            numpy_slice_tensors_[...,:3] = self.numpy_slice_tensors[...,:3]
            numpy_slice_tensors_[...,4:6] = self.numpy_slice_tensors[...,3:5]
            numpy_slice_tensors_[...,8] = self.numpy_slice_tensors[...,5]
            numpy_slice_tensors_[...,3] = self.numpy_slice_tensors[...,1]
            numpy_slice_tensors_[...,6] = self.numpy_slice_tensors[...,2]
            numpy_slice_tensors_[...,7] = self.numpy_slice_tensors[...,4]
            vtk_tensor = medipy.vtk.bridge.numpy_array_to_vtk(numpy_slice_tensors_)
            vtk_tensor.GetPointData().SetActiveTensors(vtk_tensor.GetPointData().GetScalars().GetName())

            self._glyph.SetInput(vtk_tensor)
            self._glyph.SetSource(self._sphere.GetOutput())
            
            # set up a stripper for faster rendering
            stripper = vtk.vtkStripper()
            stripper.SetInput(self._glyph.GetOutput())
              
            # now set up the mapper
            self._mapper.SetInput(stripper.GetOutput())
            

      

   
    ##############
    # Properties #
    ##############
    
    #def _set_colormap(self, colormap):
    #    super(ImageLayer, self)._set_colormap(colormap)
    #    self._image_map_to_colors.SetLookupTable(self._colormap.vtk_colormap)
    #    self.colormap.add_observer("vtk_colormap", self._on_vtk_colormap)
    
    #def _set_opacity(self, opacity) :
    #    super(Tensor2Layer, self)._set_opacity(opacity)
   #     self._actor.GetProperty().SetOpacity(opacity)
    
    def _get_actor(self):
        "VTK ImageActor."
        return self._actor
    
    #####################
    # Private interface #
    #####################
    
    def _on_vtk_colormap(self, dummy) :
        """ Event handler called when the ID of colormap.vtk_colormap changes.
        """
        
        self._set_colormap(self._colormap)
