##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import numpy
from vtk import vtkImageChangeInformation, vtkImageReslice

import medipy.base
from medipy.base import LateBindingProperty
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
        
        FIXME : image direction is not (yet) used. Since VTK images are not
        oriented, as opposed to MediPy's or ITK's, the world_to_slice should
        include direction information.
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
        # The reslicer works with physical coordinates, and thus its output
        # origin /must/ be calculated automatically. However, since we want to
        # position the sliced image at its transformed origin, the origin of the
        # sliced image fed to the actor /must/ be the transformed layer's image
        # origin (here with an altitude of 0).
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
    
    ##############
    # Properties #
    ##############
    
    def _get_world_to_slice(self) :
        "Projection matrix from world frame to 2D slice frame."
        return self._world_to_slice
    
    def _set_world_to_slice(self, world_to_slice) :
        self._world_to_slice = world_to_slice
        self._slice_to_world = numpy.linalg.inv(world_to_slice)
        
        # Normalize to 3D if needed
        world_to_slice_3d = None
        if world_to_slice.shape == (2,2) :
            temp = numpy.vstack(([0,0], world_to_slice))
            column = numpy.asarray([1,0,0]).reshape(3,1)
            world_to_slice_3d = numpy.hstack((column, temp))
        else :
            world_to_slice_3d = world_to_slice.copy()
        
        # Update reslicer, numpy axes -> VTK axes
        vtk_matrix = world_to_slice_3d
        
        temp = vtk_matrix[0,:].copy()
        vtk_matrix[0,:] = vtk_matrix[2,:]
        vtk_matrix[2,:] = temp
        temp = vtk_matrix[:,0].copy()
        vtk_matrix[:,0] = vtk_matrix[:,2]
        vtk_matrix[:,2] = temp
        self._reslicer.SetResliceAxesDirectionCosines(vtk_matrix.ravel())
        
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
        
        self._update_vtk_image_position()
        self._update_change_information()
        
        self._reslicer.SetInput(self._vtk_image)

    def _get_physical_position(self) :
        "Physical position through which the slicing plane passes."
        return self._physical_position
    
    def _set_physical_position(self, physical_position) :
        
        # Normalize dimension of physical_position w.r.t. to the image
        if len(physical_position) > self._image.ndim :
            physical_position = physical_position[-self._image.ndim:]
        elif len(physical_position) < self._image.ndim :
            prefix = (self._image.ndim-len(physical_position))*(0,)
            physical_position = prefix+tuple(physical_position)
        # else : nothing to do, dimension matches
        
        if len(physical_position) < 3 :
            prefix = (3-len(physical_position))*(0,)
            physical_position_3d = prefix+tuple(physical_position)
        else :
            physical_position_3d = physical_position
        index_position_3d = numpy.subtract(physical_position, self._image.origin)/self._image.spacing
        
        self._physical_position = physical_position
        self._index_position = numpy.subtract(physical_position, self._image.origin)/self._image.spacing
        
        if self._display_coordinates == "physical" :
            self._reslicer.SetResliceAxesOrigin(*reversed(physical_position_3d))
        else :
            self._reslicer.SetResliceAxesOrigin(*reversed(index_position_3d))
    
    def _get_index_position(self) :
        "Index position through which the slicing plane passes."
        return self._index_position
    
    def _set_index_position(self, index_position) :
        physical_position = self._image.origin + self._image.spacing*index_position
        self._set_physical_position(physical_position)
    
    def _get_display_coordinates(self) :
        "Display image using physical or index coordinates."
        return self._display_coordinates
    
    def _set_display_coordinates(self, display_coordinates) :
        if display_coordinates not in ["physical", "index"] :
            raise medipy.base.Exception("Unknown display coordinates : %s"%(display_coordinates,))
        
        self._display_coordinates = display_coordinates
        
        self._update_vtk_image_position()
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
    
    def _update_vtk_image_position(self) :
        """ Update the origin and spacing of the vtk image with respect to image
            and display coordinates.
        """
        
        if None in [self._vtk_image, self._display_coordinates] :
            return
            
        if self._display_coordinates == "physical" :
            origin = list(reversed(self.image.origin))
            if len(origin) == 2 :
                origin.append(0)
            self._vtk_image.SetOrigin(origin)
            
            spacing = list(reversed(self.image.spacing))
            if len(spacing) == 2 :
                spacing.append(1)
            self._vtk_image.SetSpacing(spacing)
        else :
            self._vtk_image.SetOrigin(0,0,0)
            self._vtk_image.SetSpacing(1,1,1)
    
    def _update_change_information(self) :
        "Set the origin of the actor input for correct placement."
        
        if None in [self._world_to_slice, self._image] :
            return
        
        # Projected origin, altitude 0
        matrix = numpy.reshape(self._reslicer.GetResliceAxesDirectionCosines(),
                               (3,3))
        origin = numpy.dot(matrix, self._vtk_image.GetOrigin())
        origin[2] = 0
        self._change_information.SetOutputOrigin(origin)
