##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import math

import wx

import medipy.base
import medipy.gui.image

class ImageGrid(wx.Panel, medipy.base.Observable):
    """ wx container for a list of layers. It synchronizes the view mode (i.e.
        slice_mode, interpolation, display_coordinates, scalar_bar_visibility,
        orientation_visibility, corner_annotations_visibility and convention)
        of the images it contains. Additionally, the cursor_position, center,
        zoom and display_range of the images might be synchronized.
    """
    
    def __init__(self, parent, slice_mode="multiplanar", layers=None,
                 interpolation=False,
                 display_coordinates="physical", scalar_bar_visibility = False,
                 orientation_visibility=True, corner_annotations_visibility=False,
                 convention="radiological", *args, **kwargs):
        
        # Display settings
        self._slice_mode = None
        self._interpolation = None
        self._display_coordinates = None
        self._scalar_bar_visibility = None
        self._orientation_visibility = None
        self._corner_annotations_visibility = None
        self._convention = None
        
        # Images on the grid
        self._images = []
        # Index of the active image
        self._active = None
        # Synchronization state
        self._synchronization = {
            "cursor_position" : False,
            "center" : False,
            "zoom" : False,
            "display_range" : False
        }
        
        self._sizer = wx.GridSizer(cols=0, rows=0, vgap=0, hgap=0)
        
        ##################
        # Initialization #
        ##################
        wx.Panel.__init__(self, parent, *args, **kwargs)
        medipy.base.Observable.__init__(self, ["active"])
        self.SetSizer(self._sizer)
        
        self._set_slice_mode(slice_mode)
        self._set_interpolation(interpolation)
        self._set_display_coordinates(display_coordinates)
        self._set_scalar_bar_visibility(scalar_bar_visibility)
        self._set_orientation_visibility(orientation_visibility)
        self._set_corner_annotations_visibility(corner_annotations_visibility)
        self._set_convention(convention)
        
    def append(self, layers):
        """ Append an image (described by its layers' attributes) to the grid.
        """
        
        self.insert(len(self), layers)
        
    def insert(self, index, layers):
        """ Insert an image (described by its layers' attributes) in the grid,
            before the specified index
        """
        
        # Create the gui.Image
        image = medipy.gui.image.Image(self, 
            self.slice_mode, 
            layers, 
            layers[0]["image"].annotations, 
            self.interpolation,
            self.display_coordinates, self.scalar_bar_visibility, 
            self.orientation_visibility, self.corner_annotations_visibility,
            self.convention
        )
        image.Bind(wx.EVT_LEFT_DOWN, self.OnImageClicked)
        image.Bind(wx.EVT_MIDDLE_DOWN, self.OnImageClicked)
        image.Bind(wx.EVT_RIGHT_DOWN, self.OnImageClicked)
        image.Bind(wx.EVT_MOUSEWHEEL, self.OnImageClicked)
        
        # Add it to the list
        self._images.insert(index, image)
        
        # Set the synchronization
        if self._images :
            reference_image = self._images[0]
            if self._synchronization["cursor_position"] :
                image.cursor_physical_position = reference_image.cursor_physical_position
            if self._synchronization["center"] :
                image.center_physical_position = reference_image.center_physical_position
            if self._synchronization["zoom"] :
                image.zoom = reference_image.zoom
            if self._synchronization["display_range"] :
                display_range = reference_image.get_layer_colormap(0).display_range
                image.get_layer_colormap(0).display_range = display_range
        image.synchronize_on_event["zoom"] = self._synchronization["zoom"]
        image.add_observer("cursor_position", self._on_cursor_position)
        image.add_observer("center", self._on_center)
        image.add_observer("colormap_display_range", self._on_colormap_display_range)
            
        # Set children for synchronization :
        for image_index, other_image in enumerate(self._images) :
            next_image = self._images[(image_index+1)%len(self._images)]
            other_image.delete_all_children()
            other_image.append_child(next_image)
        
        # Synchronize if more than one image
        if len(self) > 1 :
            reference = self[(index+1)%len(self)]
            if self._synchronization["cursor_position"] :
                reference.notify_observers("cursor_position")
            if self._synchronization["center"] :
                reference.notify_observers("center")
            if self._synchronization["zoom"] :
                reference.notify_observers("zoom")
            if self._synchronization["display_range"] :
                reference.notify_observers("display_range")
        
        # Refresh the GUI
        self._sizer.Insert(index, image, 1, wx.EXPAND)
        self._adjust_sizer()
        
        self.active = index
        
    def delete(self, index):
        """ Delete an image from the grid.
        """
        
        image = self._images[index]
        self._sizer.Detach(image)
        self._adjust_sizer()
        
        self._images[index].Close()
        del self._images[index]
        
        # Adjust the active image
        if len(self._images) > 0 and self._active >= len(self._images) : 
            self.active = len(self._images)-1
        elif len(self._images) > 0 :
            self.active = self.active
    
    def __len__(self):
        """ Return the number of images in the grid.
        """
        
        return len(self._images)
    
    def __iter__(self):
        """ Iterate through the grid images.
        """
        
        return self._images.__iter__()
    
    def __getitem__(self, index):
        """ Return an image from the grid
        """
        
        return self._images[index]
    
    def get_synchronization(self, attribute):
        """ Return the synchronization state of an image attribute (may be
            cursor_position, center, zoom or display_range).
        """
        
        return self._synchronization[attribute]
    
    def set_synchronization(self, attribute, value):
        """ Set the synchronization state of an image attribute (may be
            cursor_position, center, zoom or display_range).
        """
        
        if attribute not in self._synchronization :
            raise medipy.base.Exception("Invalid synchronization attribute: {0}".format(attribute))
        self._synchronization[attribute] = value
        
        # Other attributes are handled in observers since they are not properties
        # of the image
        if attribute == "zoom" :
            for image in self._images :
                if hasattr(image, attribute) :
                    image.synchronize_on_event[attribute] = value
    
    ##############
    # Properties #
    ##############
    
    def _get_slice_mode(self) :
        return self._slice_mode
    
    def _set_slice_mode(self, value) :
        self._slice_mode = value
        for image in self._images :
            image.slice_mode = value
            image.render()
    
    def _get_interpolation(self) :
        return self._interpolation
    
    def _set_interpolation(self, value) :
        self._interpolation = value
        for image in self._images :
            image.interpolation = value
            image.render()
    
    def _get_display_coordinates(self) :
        return self._display_coordinates
    
    def _set_display_coordinates(self, value) :
        self._display_coordinates = value
        for image in self._images :
            image.display_coordinates = value
            image.render()
    
    def _get_scalar_bar_visibility(self) :
        return self._scalar_bar_visibility
    
    def _set_scalar_bar_visibility(self, value) :
        self._scalar_bar_visibility = value
        for image in self._images :
            image.scalar_bar_visibility = value
            image.render()
    
    def _get_orientation_visibility(self) :
        return self._orientation_visibility
    
    def _set_orientation_visibility(self, value) :
        self._orientation_visibility = value
        for image in self._images :
            image.orientation_visibility = value
            image.render()
    
    def _get_corner_annotations_visibility(self) :
        return self._corner_annotations_visibility
    
    def _set_corner_annotations_visibility(self, value) :
        self._corner_annotations_visibility = value
        for image in self._images :
            image.corner_annotations_visibility = value
            image.render()
    
    def _get_convention(self) :
        return self._convention
    
    def _set_convention(self, value) :
        self._convention = value
        for image in self._images :
            image.convention = value
            image.render()
    
    def _get_active(self):
        return self._active
    
    def _set_active(self, index):
        self._active = index
        
        for item in self._sizer.GetChildren() :
            other_image = item.GetWindow()
            if other_image == self._images[index] :
                other_image.SetBackgroundColour(wx.GREEN)
            else :
                other_image.SetBackgroundColour(wx.BLACK)
        self.notify_observers("active")
    
    active = property(_get_active, _set_active)
    slice_mode = property(_get_slice_mode, _set_slice_mode)
    interpolation = property(_get_interpolation, _set_interpolation)
    display_coordinates = property(_get_display_coordinates, 
                                   _set_display_coordinates)
    scalar_bar_visibility = property(_get_scalar_bar_visibility, 
                                     _set_scalar_bar_visibility)
    orientation_visibility = property(_get_orientation_visibility, 
                                      _set_orientation_visibility)
    corner_annotations_visibility = property(_get_corner_annotations_visibility, 
                                             _set_corner_annotations_visibility)
    convention = property(_get_convention, _set_convention)

    ##################
    # Event handlers #
    ##################
    
    def OnImageClicked(self, event):
        if self.active is None or event.GetEventObject() != self._images[self.active] :
            self.active = self._images.index(event.GetEventObject())
    
    def _on_cursor_position(self, event) :
        
        if not self._synchronization["cursor_position"] :
            return
        
        for image in self._images :
            if image != event.object :
                image.remove_observer("cursor_position", self._on_cursor_position)
                if self.display_coordinates == "index" :
                    image.cursor_index_position = event.object.cursor_index_position
                else :
                    image.cursor_physical_position = event.object.cursor_physical_position
                image.render()
                image.add_observer("cursor_position", self._on_cursor_position)
    
    def _on_center(self, event) :
        
        if not self._synchronization["center"] :
            return
        
        for image in self._images :
            if image != event.object :
                image.remove_observer("center", self._on_center)
                if self.display_coordinates == "index" :
                    image.center_index_position = event.object.center_index_position
                else :
                    image.center_physical_position = event.object.center_physical_position
                image.render()
                image.add_observer("center", self._on_center)

    def _on_colormap_display_range(self, event):
        
        if not self._synchronization["display_range"] :
            return
        
        for image in self._images :
            if image != event.object :
                image.remove_observer("colormap_display_range", self._on_colormap_display_range)
                source_colormap = event.object.get_layer_colormap(event.layer_index)
                destination_colormap = image.get_layer_colormap(event.layer_index)
                destination_colormap.display_range = source_colormap.display_range
                image.render()
                image.add_observer("colormap_display_range", self._on_colormap_display_range)

    #####################
    # Private interface #
    #####################
    
    def _adjust_sizer(self):
        """ Re-arrange the sizer so that images are layed-out on a grid as 
            square as possible.
        """
        
        nb_objects = len(self._sizer.GetChildren())
        rows = max(int(math.ceil(math.sqrt(nb_objects))),1)
        self._sizer.SetRows(rows)
        self._sizer.SetCols(rows)
        self._sizer.Layout()
    