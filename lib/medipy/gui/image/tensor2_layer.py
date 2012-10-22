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

class Tensor2Layer(Layer) :
    """ Layer showing its diffusion data. The image will be positionned in the
        viewport such that its origin (i.e. lower left corner) matches the 
        transformed origin of the layer's image (with an altitude of 0).
    """
    
    def __init__(self, world_to_slice, tensor_image, display_coordinates="physical",
                 colormap=None, opacity = 1.0, display_mode="principal_direction_voxel") :
        

        display_mode = "ellipsoid"

        if display_mode=="principal_direction_voxel" :
            ############################
            # Property-related members #
            ############################
            
            self._actor = vtkImageActor()          
            
            ######################
            # Initialize members #
            ######################
            
            super(Tensor2Layer, self).__init__(world_to_slice, tensor_image, display_coordinates, colormap, opacity)
            
            self._change_information.Update()
            vtk_slice_tensors = self._change_information.GetOutput()
            print ":", vtk_slice_tensors.GetDimensions(),vtk_slice_tensors.GetSpacing(),vtk_slice_tensors.GetOrigin()
            numpy_slice_tensors = medipy.vtk.bridge.vtk_to_numpy_array(vtk_slice_tensors)
            numpy_eigenvalues,numpy_eigenvectors = spectral_decomposition(numpy_slice_tensors)
            numpy_principal_direction = numpy.ascontiguousarray(numpy_eigenvectors[...,::3])
            numpy_principal_direction[:,:,:] = (0,1,0)

            vtk_principal_direction = medipy.vtk.bridge.numpy_array_to_vtk(numpy.cast[numpy.uint8](numpy_principal_direction*255))
            vtk_principal_direction.SetSpacing(vtk_slice_tensors.GetSpacing())
            vtk_principal_direction.SetOrigin(vtk_slice_tensors.GetOrigin())
            print "::", vtk_principal_direction.GetDimensions(),vtk_principal_direction.GetSpacing(),vtk_principal_direction.GetOrigin()

            self._actor.SetInput(vtk_principal_direction)
            self._actor.InterpolateOff()

        elif display_mode=="principal_direction_line" :

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

            ######################
            # Initialize members #
            ######################
            
            super(Tensor2Layer, self).__init__(world_to_slice, tensor_image, display_coordinates, colormap, opacity)
            
            self._change_information.Update()
            vtk_slice_tensors = self._change_information.GetOutput()
            print ":", vtk_slice_tensors.GetDimensions(),vtk_slice_tensors.GetSpacing(),vtk_slice_tensors.GetOrigin()
            numpy_slice_tensors = medipy.vtk.bridge.vtk_to_numpy_array(vtk_slice_tensors)
            numpy_eigenvalues,numpy_eigenvectors = spectral_decomposition(numpy_slice_tensors)
            numpy_principal_direction = numpy.ascontiguousarray(numpy_eigenvectors[...,::3])
            numpy_principal_direction[:,:,:] = (0,1,0)

            vtk_principal_direction = medipy.vtk.bridge.numpy_array_to_vtk(numpy.cast[numpy.uint8](numpy_principal_direction*255))
            vtk_principal_direction.SetSpacing(vtk_slice_tensors.GetSpacing())
            vtk_principal_direction.SetOrigin(vtk_slice_tensors.GetOrigin())
            print "::", vtk_principal_direction.GetDimensions(),vtk_principal_direction.GetSpacing(),vtk_principal_direction.GetOrigin()

            self._glyph.ScalingOn()
            self._glyph.SetScaleModeToScaleByVector()
            self._glyph.SetColorModeToColorByVector()
            self._glyph.SetScaleFactor(0.5)
            self._glyph.SetSource(self._line.GetOutput())
            self._glyph.SetInput(vtk_principal_direction)
            self._glyph.ClampingOff()


            # set up a stripper for faster rendering
            stripper = vtk.vtkStripper()
            stripper.SetInput(self._glyph.GetOutput())
              
            # get the maximum norm of the data
            maxNorm = 1.0
              
            # now set up the mapper
            self._mapper.SetInput(stripper.GetOutput())
            self._mapper.SetScalarRange(0, maxNorm) 

            # set up the actor
            self._actor.SetMapper(self._mapper)

        elif display_mode=="ellipsoid" :

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

            ######################
            # Initialize members #
            ######################
            
            super(Tensor2Layer, self).__init__(world_to_slice, tensor_image, display_coordinates, colormap, opacity)

            self._change_information.Update()
            vtk_slice_tensors = self._change_information.GetOutput()
            print ":", vtk_slice_tensors.GetDimensions(),vtk_slice_tensors.GetSpacing(),vtk_slice_tensors.GetOrigin()
            numpy_slice_tensors = medipy.vtk.bridge.vtk_to_numpy_array(vtk_slice_tensors)
            numpy_slice_tensors_ = numpy.zeros(numpy_slice_tensors.shape[:3]+(9,),dtype=numpy.single)
            numpy_slice_tensors_[...,:3] = numpy_slice_tensors[...,:3]
            numpy_slice_tensors_[...,4:6] = numpy_slice_tensors[...,3:5]
            numpy_slice_tensors_[...,8] = numpy_slice_tensors[...,5]
            numpy_slice_tensors_[...,3] = numpy_slice_tensors[...,1]
            numpy_slice_tensors_[...,6] = numpy_slice_tensors[...,2]
            numpy_slice_tensors_[...,7] = numpy_slice_tensors[...,4]
            vtk_tensor = medipy.vtk.bridge.numpy_array_to_vtk(numpy_slice_tensors_)
            vtk_tensor.GetPointData().SetActiveTensors(vtk_tensor.GetPointData().GetScalars().GetName())

            self._sphere.SetThetaResolution (8)
            self._sphere.SetPhiResolution (8)

            self._glyph.SetInput(vtk_tensor)
            self._glyph.SetSource(self._sphere.GetOutput())
            self. _glyph.SetScaleFactor (700)
            self._glyph.SetColorGlyphs(1)

            # set up a stripper for faster rendering
            stripper = vtk.vtkStripper()
            stripper.SetInput(self._glyph.GetOutput())
              
            # get the maximum norm of the data
            maxNorm = 1.0
              
            # now set up the mapper
            self._mapper.SetInput(stripper.GetOutput())
            self._mapper.SetScalarRange(0, maxNorm) 

            # set up the actor
            self._actor.SetMapper(self._mapper)

            

      

   
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
