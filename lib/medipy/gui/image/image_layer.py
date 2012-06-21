##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import numpy

from vtk import vtkImageActor, vtkImageMapToColors

from layer import Layer 

class ImageLayer(Layer) :
    """ Layer showing its data as an image. The image will be positionned in the
        viewport such that its origin (i.e. lower left corner) matches the 
        transformed origin of the layer's image (with an altitude of 0).
    """
    
    def __init__(self, world_to_slice, image, display_coordinates="physical",
                 colormap=None, opacity = 1.0) :
        
        ############################
        # Property-related members #
        ############################
        
        self._actor = vtkImageActor()
        
        ###################
        # Private members #
        ###################
        
        self._image_map_to_colors = vtkImageMapToColors()
        
        ######################
        # Initialize members #
        ######################
        
        super(ImageLayer, self).__init__(
            world_to_slice, image, display_coordinates, colormap, opacity)
        
        self._image_map_to_colors.SetInputConnection(self._change_information.GetOutputPort())
        self._actor.SetInput(self._image_map_to_colors.GetOutput())
        self._actor.InterpolateOff()
    
    def world_to_index(self, world) :
        # Convert to index coordinate in resliced image (VTK order)
        index = numpy.divide(
            numpy.subtract(world, self._change_information.GetOutputOrigin()),
            self._change_information.GetOutputSpacing())
        # Set height to 0, since the picked value will depend on the position
        # of the actor
        index[2] = 0
        
        # We're off by 1/2 pixel ! TODO : adjust actor position
        # cf. http://www.vtk.org/pipermail/vtkusers/2005-May/079848.html
        index[0] += 0.5
        index[1] += 0.5
        
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
        index = self.world_to_index(world)
        return self._image.index_to_physical(index)    
    
    ##############
    # Properties #
    ##############
    
    def _set_colormap(self, colormap):
        super(ImageLayer, self)._set_colormap(colormap)
        self._image_map_to_colors.SetLookupTable(self._colormap.vtk_colormap)
        self.colormap.add_observer("vtk_colormap", self._on_vtk_colormap)
    
    def _set_opacity(self, opacity) :
        super(ImageLayer, self)._set_opacity(opacity)
        self._actor.SetOpacity(opacity)
    
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
