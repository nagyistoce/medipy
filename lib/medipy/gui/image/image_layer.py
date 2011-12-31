##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

from vtk import vtkImageActor, vtkImageMapToColors

from layer import Layer #from medipy.gui.image.layer import Layer 

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
