##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import numpy
from vtk import vtkImageChangeInformation, vtkImageReslice, vtkMatrix4x4

import medipy.base
import medipy.base.array
from medipy.base import LateBindingProperty
import medipy.gui
import medipy.vtk

class Layer(object) :
    """ Representation of a plane slice through a 2D or 3D image.
    
        This abstract class encapsulates an image and a colormap, and leaves the
        concrete representations to derived classes, through the "actor" 
        property.
        
        The slice plane is given by the world_to_slice matrix, and the position
        of the plane can be given either in physical or index coordinates.
        
        The layer can be displayed either using the physical coordinates, or the
        index coordinates. In the former case, anisotropic voxels will appear
        anisotropically, and the origin of the sliced image will be the origin
        of the image. In the latter case, anisotropic voxels will appear
        isotropically, and the origin of the sliced image will be 0. A global
        layer opacity (between 0 and 1, inclusive) can also be specified.
        
        Note that no change of position is performed when changing the display
        coordinates : it is up to the client to do this if needed.
    """
    
    def __init__(self, world_to_slice, image, display_coordinates="physical",
                 colormap=None, opacity = 1.0) :
        
        ############################
        # Property-related members #
        ############################
        
        self._world_to_slice = None
        self._slice_to_world = None
        
        self._image = None
        
        self._index_position = None
        self._physical_position = None
        
        self._display_coordinates = None
        
        self._colormap = None
        
        self._opacity = None
        
        ###################
        # Private members #
        ###################
        
        # VTK image with the same content as self._image
        self._vtk_image = None
        self._reslicer = vtkImageReslice()
        self._reslicer_axes_inverse = vtkMatrix4x4()
        # The reslice matrix will be from voxel space to slice space. Since we
        # want to position the sliced image at its transformed origin, we must 
        # set the origin of the sliced image after the slicing.
        self._change_information = vtkImageChangeInformation()
        
        ##################
        # Initialization #
        ##################
        
        self._reslicer.SetInterpolationModeToNearestNeighbor()
        self._reslicer.SetOutputDimensionality(2)
        self._change_information.SetInputConnection(self._reslicer.GetOutputPort())
        
        self._set_world_to_slice(world_to_slice)
        self._set_image(image)
        self._set_display_coordinates(display_coordinates)
        if colormap is None :
            colormap = medipy.gui.Colormap(medipy.gui.colormaps["gray"], None)
        self._set_colormap(colormap)
        self._set_opacity(opacity)
        
        if colormap is not None and colormap.display_range is None and image is not None :
            # Default to window center and width if present in the metadata, 
            # or to the full image range
            window_center = self._image.metadata.get("window_center", None)
            window_width = self._image.metadata.get("window_width", None)
            if None not in [window_center, window_width] :
                display_range = (
                    window_center-window_width/2,window_center+window_width/2)
            else :
                display_range = (image.data.min(), image.data.max())
            self._colormap.display_range = display_range
    
    def world_to_index(self, world) :
        """ Convert a VTK world coordinate (VTK order) to the corresponding 
            image index (numpy order).
        """
        
        # Make sure the pipeline is up-to-date
        self._change_information.Update()
        
        # Convert to index coordinate in resliced image (VTK order)
        index = numpy.divide(
            numpy.subtract(world, self._change_information.GetOutputOrigin()),
            self._change_information.GetOutputSpacing())
        # Set height to 0, since the picked value will depend on the position
        # of the actor
        index[2] = 0
        
        # Apply the reslicer transform (homogeneous coordinates, VTK order),
        # converting to the non-sliced image
        physical = numpy.add(
            numpy.multiply(index, self._reslicer.GetOutput().GetSpacing()),
            self._reslicer.GetOutput().GetOrigin())
        physical = numpy.hstack((physical, 1.))
        physical = self._reslicer.GetResliceAxes().MultiplyPoint(physical)
        physical = [physical[i]/physical[3] for i in range(3)]
        
        # Convert to index coordinate in the non-sliced image (VTK order)
        index = numpy.divide(
            numpy.subtract(physical, self._vtk_image.GetOrigin()),
            self._vtk_image.GetSpacing())
        
        # VTK order -> NumPy order
        index = index[::-1]
        
        return index

    def world_to_physical(self, world) :
        """ Convert a VTK world coordinate (VTK order) to the corresponding 
            image physical coordinate (numpy order).
        """
        
        index = self.world_to_index(world)
        return self._image.index_to_physical(index) 
    
    def index_to_world(self, index) :
        """ Convert an image index (numpy order) to the corresponding VTK world
            coordinate (VTK order).
        """
        
        index = medipy.base.array.reshape(numpy.asarray(index), (3,), 
            "constant", False, value=0)
        
        # Make sure the pipeline is up-to-date
        self._change_information.Update()
        
        # NumPy Order -> VTK order
        index = index[::-1]
        
        # Convert from the non-sliced image point coordinates (VKT order)
        physical = numpy.add(
            numpy.multiply(index, self._vtk_image.GetSpacing()),
            self._vtk_image.GetOrigin())
        
        # Apply the inverse reslicer transform (homogeneous coordinates, VTK 
        # order), converting to the sliced image
        physical = numpy.hstack((physical, 1.))
        physical = self._reslicer_axes_inverse.MultiplyPoint(physical)
        physical = [physical[i]/physical[3] for i in range(3)]
        
        # Convert to index coordinate in resliced image (VTK order)
        index = numpy.divide(
            numpy.subtract(physical, self._reslicer.GetOutput().GetOrigin()),
            self._reslicer.GetOutput().GetSpacing())
        
        # Convert to world coordinates
        world = numpy.add(
            numpy.multiply(index, self._change_information.GetOutputSpacing()),
            self._change_information.GetOutputOrigin())
        
        return world
    
    def physical_to_world(self, physical) :
        """ Convert an image physical coordinate (numpy order) to the 
            corresponding VTK world coordinate (VTK order).
        """
        
        index = self._image.physical_to_index(physical)
        return self.index_to_world(index)
    
    ##############
    # Properties #
    ##############
    
    def _get_world_to_slice(self) :
        "Projection matrix from world frame to 2D slice frame."
        return self._world_to_slice
    
    def _set_world_to_slice(self, world_to_slice) :
        self._world_to_slice = world_to_slice
        self._slice_to_world = numpy.linalg.inv(world_to_slice)
        
        self._update_reslicer_matrix()
        self._update_change_information()
        
    def _get_slice_to_world(self) :
        "Inverse of projection matrix from world frame to 2D slice frame."
        return self._slice_to_world
    
    def _get_image(self) :
        "Image to be sliced."
        return self._image
    
    def _set_image(self, image) :
        if self._image :
            self._image.remove_observer("modified", self._on_image_modified)
            
        self._image = image
        
        self._vtk_image = medipy.vtk.build_vtk_image(self._image, self._vtk_image)
        
        self._image.add_observer("modified", self._on_image_modified)
        
        # Reset the origin and spacing of the VTK image since the reslice matrix
        # will be from voxel space to slice space
        self._vtk_image.SetOrigin(0,0,0)
        self._vtk_image.SetSpacing(1,1,1)
        # Update the pipeline
        self._update_reslicer_matrix()
        self._update_change_information()
        
        self._reslicer.SetInput(self._vtk_image)

    def _get_physical_position(self) :
        "Physical position through which the slicing plane passes."
        return self._physical_position
    
    def _set_physical_position(self, physical_position) :
        
        physical_position = numpy.asarray(physical_position)
        
        # Normalize dimension of physical_position w.r.t. to the image
        physical_position_image = medipy.base.array.reshape(
            physical_position, (self._image.ndim,), "constant", False, value=0)
        
        self._physical_position = physical_position_image
        self._index_position = self._image.physical_to_index(physical_position_image)
        
        self._reslicer.SetResliceAxesOrigin(self._index_position[::-1])
        vtkMatrix4x4.Invert(
            self._reslicer.GetResliceAxes(), self._reslicer_axes_inverse)
    
    def _get_index_position(self) :
        "Index position through which the slicing plane passes."
        return self._index_position
    
    def _set_index_position(self, index_position) :
        index_position = numpy.asarray(index_position)
        index_position_image = medipy.base.array.reshape(
            index_position, (self._image.ndim,), "constant", False, value=0)
        
        physical_position_image = self._image.index_to_physical(index_position_image)
        self._set_physical_position(physical_position_image)
    
    def _get_display_coordinates(self) :
        "Display image using physical or index coordinates."
        return self._display_coordinates
    
    def _set_display_coordinates(self, display_coordinates) :
        if display_coordinates not in ["physical", "nearest_axis_aligned", "index"] :
            raise medipy.base.Exception("Unknown display coordinates : %s"%(display_coordinates,))
        
        self._display_coordinates = display_coordinates
        
        self._update_reslicer_matrix()
        self._update_change_information()
        
        if self._physical_position is not None :
            self._set_physical_position(self._get_physical_position())
    
    def _get_colormap(self) :
        "Colormap to be applied to the image."
        return self._colormap
    
    def _set_colormap(self, colormap) :
        self._colormap = colormap
    
    def _get_opacity(self) :
        "Global opacity of the layer."
        return self._opacity
    
    def _set_opacity(self, opacity) :
        self._opacity = opacity
        
    def _get_actor(self) :
        "VTK actor of the layer. Must be defined in concrete derived classes."
        
        raise NotImplementedError()
    
    world_to_slice = LateBindingProperty(_get_world_to_slice, 
                                         _set_world_to_slice)
    slice_to_world = LateBindingProperty(_get_slice_to_world)
    image = LateBindingProperty(_get_image, _set_image)
    physical_position = LateBindingProperty(_get_physical_position, 
                                            _set_physical_position)
    index_position = LateBindingProperty(_get_index_position, 
                                         _set_index_position)
    display_coordinates = LateBindingProperty(_get_display_coordinates, 
                                              _set_display_coordinates)
    colormap = LateBindingProperty(_get_colormap, _set_colormap)
    opacity = LateBindingProperty(_get_opacity, _set_opacity)
    actor = LateBindingProperty(_get_actor)
    
    #####################
    # Private interface #
    #####################
    
    def _on_image_modified(self, event):
        self._vtk_image.Modified()
    
    def _update_reslicer_matrix(self):
        """ Update the reslicer matrix with respect to the world_to_slice matrix
            and the display coordinates
        """
        
        if None in [self._world_to_slice, self._image, self._display_coordinates] :
            return
        
        # Reshape to 3x3 matrix
        world_to_slice_3d = medipy.base.array.reshape(self.world_to_slice, (3,3),
            "constant", False, value=0)
        # Add ones on the diagonal when necessary
        for rank in range(3) :
            if numpy.less_equal(self.world_to_slice.shape, rank).all() : 
                world_to_slice_3d[3-rank-1, 3-rank-1] = 1.
        
        # Same for direction
        direction_3d = medipy.base.array.reshape(self.image.direction, (3,3),
            "constant", False, value=0)
        # Add ones on the diagonal when necessary
        for rank in range(3) :
            if numpy.less_equal(self.image.direction.shape, rank).all() : 
                direction_3d[3-rank-1, 3-rank-1] = 1.
        
        # The reslice matrix from orignal voxel space to slice space
        # (i.e. input->output) 
        if self._display_coordinates == "index" :
            matrix = world_to_slice_3d
        elif self._display_coordinates == "nearest_axis_aligned" :
            nearest = medipy.base.coordinate_system.best_fitting_axes_aligned_matrix(direction_3d)
            matrix = numpy.dot(world_to_slice_3d, nearest)
        else :
            matrix = numpy.dot(world_to_slice_3d, direction_3d)
        
        self._input_index_to_output_index = matrix.copy()
        self._output_index_to_input_index = numpy.linalg.inv(matrix)
        
        # vtkImageReslice.SetResliceAxesDirectionCosines fills the rotation
        # part of the 3x3 matrix column wise, and this rotation part must be
        # the transform from output to input. We then must pass the transpose
        # of the inverse of matrix, i.e. (M^{-1})^T. For orthogonal matrices, 
        # this is equal to M itself, but we're never too careful.
        matrix = numpy.linalg.inv(matrix).T
        
        # Update reslicer, numpy axes -> VTK axes
        self._reslicer.SetResliceAxesDirectionCosines(matrix[::-1,::-1].ravel())
        vtkMatrix4x4.Invert(
            self._reslicer.GetResliceAxes(), self._reslicer_axes_inverse)
    
    def _update_change_information(self) :
        """ Update the origin of the ImageChangeInformation filter to the
            transformed origin of the image. The origin of the image is set to
            the coordinate-wise minimum of the transformed bounding box, 
            the spacing is set to the absolute value of the transformed spacing.
        """
        
        if None in [self._world_to_slice, self._image, self._display_coordinates] :
            return
        
        if self.display_coordinates in ["physical", "nearest_axis_aligned"] :
            # Transform the bounding box of the image to find the coordinates
            # of the first voxel
            # TODO : is this correct for nearest_axis_aligned or should we 
            # compute nearest*spacing*index+origin ?
            
            begin = numpy.dot(self._world_to_slice, 
                              self._image.index_to_physical((0,0,0)))
            end = numpy.dot(self._world_to_slice, 
                            self._image.index_to_physical(self._image.shape))
            
            changed_origin = numpy.minimum(begin, end)
            # Set altitude to 0
            changed_origin[0] = 0

            changed_spacing = numpy.abs(numpy.dot(self._world_to_slice, 
                                                  self._image.spacing))
        else :
            changed_origin = (0,0,0)
            changed_spacing = (1,1,1)
        
        self._change_information.SetOutputOrigin(changed_origin[::-1])
        self._change_information.SetOutputSpacing(changed_spacing[::-1])
