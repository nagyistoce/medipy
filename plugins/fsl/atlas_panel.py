##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import os

import wx

import medipy.base
import medipy.io

from atlas import Atlas
import utils

class AtlasPanel(wx.Panel):
    def __init__(self, parent, *args, **kwargs):
        self._atlas = None
        self._label_image = None
        self._image = None
        
        wx.Panel.__init__(self, parent, *args, **kwargs)
        
        self._atlases = wx.Choice(self)
        self._atlases.Bind(wx.EVT_CHOICE, self.OnAtlases)
        
        self._label = wx.TextCtrl(self, style=wx.TE_READONLY)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self._atlases, flag=wx.EXPAND)
        sizer.Add(self._label, flag=wx.EXPAND)
        self.SetSizer(sizer)
        
        self._setup_gui()
        
    ##############
    # Properties #
    ##############
    
    def _get_image(self):
        return self._image
    
    def _set_image(self, image):
        if self._image is not None :
            self._image.remove_observer("cursor_position", self.on_cursor_position)
        self._image = image
        if self._image is not None :
            self._image.add_observer("cursor_position", self.on_cursor_position)
        self._update_label()
    
    def _get_atlas(self):
        return self._atlas
    
    def _set_atlas(self, atlas):
        self._atlas = atlas
        self._label_image = medipy.io.load(self._atlas.images[0][1])
        
        self._atlases.SetSelection(-1)
        for index in range(self._atlases.GetCount()) :
            if self._atlases.GetString(index) == self._atlas.name :
                self._atlases.SetSelection(index)
        
        self._update_label()
        
    image = property(_get_image, _set_image)
    atlas = property(_get_atlas, _set_atlas)
    
    ##################
    # Event handlers #
    ##################
    
    def OnAtlases(self, event):
        self.atlas = event.GetClientData()
    
    def on_cursor_position(self, event):
        self._update_label()
    
    #####################
    # Private interface #
    #####################
    
    def _find_atlases(self):
        """ Return a list of FSL atlases
        """
        
        fsldir = utils.environment()["FSLDIR"]
        if fsldir is None :
            return []
        
        atlases_dir = os.path.join(fsldir, "data", "atlases")
        if not os.path.isdir(atlases_dir) :
            return []
        
        atlases = []
        for entry in os.listdir(atlases_dir) :
            if not os.path.isfile(os.path.join(atlases_dir, entry)) :
                continue
            
            try :
                atlas = Atlas.read(os.path.join(atlases_dir, entry))
            except medipy.base.Exception, e :
                # Cannot read the atlas
                print e
                continue
            else :
                atlases.append(atlas)
        
        return atlases
        
    def _setup_gui(self):
        """ Setup the GUI
        """
        
        self._atlases.Clear()
        
        atlases = sorted(self._find_atlases(), key = lambda x:x.name)
        
        if len(atlases) == 0 :
            self._atlases.SetSelection(-1)
            return
        
        for atlas in atlases :
            self._atlases.Append(atlas.name, atlas)
        
        self.atlas = self._atlases.GetClientData(0)
    
    def _update_label(self):
        """ Update the label informations given the currently selected atlas
            and the current image.
        """
        
        if None in [self._atlas, self._image] :
            return
        else :
            position_physical = self._image.cursor_physical_position
            position_index = self._label_image.physical_to_index(position_physical)
            position_index = tuple([int(x) for x in position_index])
            label = self._label_image[position_index]
            
            if self.atlas.type == Atlas.Type.probabilistic and label != 0 :
                name = self._atlas.labels[label-1]
            elif self.atlas.type == Atlas.Type.label and label != 0 :
                name = self._atlas.labels[label]
            else :
                name = ""
            
            self._label.ChangeValue(name)