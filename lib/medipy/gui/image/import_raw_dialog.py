##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

#import os
#from struct import pack
#import xml.dom.minidom as minidom
#
#import numpy
#import wx
#import wx.xrc
#
#from medipy.base import find_resource
#from medipy.gui.control import Int
#from medipy.base import Image

import operator
import os
import struct

import numpy
import wx

import medipy.io

import medipy.gui.base
import medipy.gui.control
import medipy.gui.xrc_wrapper

class ImportRawDialog(medipy.gui.base.Dialog):
    
    class UI(medipy.gui.base.UI):
        def __init__(self):
            self.filename = None
            self.filesize = None
            
            self.ndim = None
            self.sizes = None
            self.sizes_sizer = wx.BoxSizer()
            
            self.dtype = None
            self.endianness = None
            
            self.offset = None
            self.skip_end = None
            
            self.message = None
            
            self.import_ = None
            self.cancel = None
            
            self.controls = ["filename", "filesize",
                             "ndim", "sizes",
                             "dtype", "endianness",
                             "offset", "skip_end",
                             "message",
                             "import_", "cancel"]
        
        def from_window(self, window, names):
            medipy.gui.base.UI.from_window(self, window, names)
            
            self.sizes.SetSizer(self.sizes_sizer)
            self.import_.SetId(wx.ID_OK)
            self.cancel.SetId(wx.ID_CANCEL)
    
    def __init__(self, parent=None, *args, **kwargs):
        self.ui = ImportRawDialog.UI()
        
        xrc_file = medipy.base.find_resource(
            os.path.join("resources", "gui", "import_raw_dialog.xrc"))
        wrappers = [medipy.gui.xrc_wrapper.BoolXMLHandler(),
                    medipy.gui.xrc_wrapper.EnumXMLHandler(),
                    medipy.gui.xrc_wrapper.FileXMLHandler(),
                    medipy.gui.xrc_wrapper.IntXMLHandler(),
                   ]
        medipy.gui.base.Dialog.__init__(self, xrc_file, "import_raw_dialog", 
            wrappers, self.ui, self.ui.controls, parent, *args, **kwargs)
        
        self.ui.filename.add_observer("value", self._on_filename)
        self.ui.dtype.add_observer("value", self._on_dtype)
        self.ui.ndim.add_observer("value", self._on_ndim)
        self.ui.endianness.add_observer("value", self._on_endianness)
        self.ui.offset.add_observer("value", self._on_offset)
        self.ui.skip_end.add_observer("value", self._on_skip_end)
        
        self._filesize = None
        self._dtype = None
        self._endianness = None
        
        self._image = None
        
        self.ui.ndim.value = 3
        
        # Update the dtype control using numpy types
        self._types = {}
        for type_ in sorted(numpy.nbytes.keys()):
            if (type_.__name__.startswith("int") or 
                type_.__name__.startswith("uint") or
                type_.__name__.startswith("float")) :
                self._types[type_.__name__] = type_
        self.ui.dtype.choices = self._types.keys()
        self.ui.dtype.value = self.ui.dtype.choices[0]
        
        # Update the endianness control
        self._endiannesses = {
            "native" : "=",
            "little" : "<",
            "big" : ">",
        }
        self.ui.endianness.choices = self._endiannesses.keys()
        self.ui.endianness.value = "native"
    
    def ShowModal(self, *args, **kwargs):
        return_code = medipy.gui.base.Dialog.ShowModal(self, *args, **kwargs)
    
        if return_code == wx.ID_OK :
            shape = [x.GetWindow().value 
                     for x in self.ui.sizes_sizer.GetChildren()]
            shape = list(reversed(shape))
            
            self._image = medipy.io.raw.load(self.ui.filename.value,
                shape, self._dtype, self._endianness, self.ui.offset.value,
                self.ui.skip_end.value)
        else :
            self._image = None

        return return_code
    
    ##############
    # Properties #
    ##############
    
    def _get_image(self):
        return self._image    
    
    image = property(_get_image)
    
    ##################
    # Event handlers #
    ##################
    
    def _on_filename(self, event):
        
        self._filesize = os.path.getsize(event.object.value)
        self.ui.filesize.SetLabel("{0} byte{1}".format(
            self._filesize, "s" if self._filesize>1 else ""))
        
        self._update_status()
    
    def _on_dtype(self, event):
        self._dtype = self._types[event.object.value]
        self._update_status()
    
    def _on_ndim(self, event):
        if not event.object.validate() :
            return
        
        if event.object.value < len(self.ui.sizes_sizer.GetChildren()) :
            while event.object.value < len(self.ui.sizes_sizer.GetChildren()) :
                control = self.ui.sizes_sizer.GetItem(len(self.ui.sizes_sizer.GetChildren())-1)
                control.GetWindow().remove_observer("value", self._on_size)
                self.ui.sizes_sizer.Remove(control)
        elif event.object.value > len(self.ui.sizes_sizer.GetChildren()) :
            while event.object.value > len(self.ui.sizes_sizer.GetChildren()) :
                control = medipy.gui.control.Int(self.ui.sizes)
                control.value = 0
                control.add_observer("value", self._on_size)
                self.ui.sizes_sizer.Add(control, 1)
        
        self.Layout()
    
    def _on_size(self, event):
        self._update_status()
    
    def _on_endianness(self, event):
        self._endianness = self._endiannesses[event.object.value]
        self._update_status()
    
    def _on_offset(self, event):
        self._update_status()
    
    def _on_skip_end(self, event):
        self._update_status()
    
    #####################
    # Private interface #
    #####################
    
    def _update_status(self):
        controls = ["filename", 
                    "dtype", "endianness",
                    "offset", "skip_end"]
        controls = [getattr(self.ui, control) for control in controls]
        
        controls.extend([x.GetWindow() for x in self.ui.sizes_sizer.GetChildren()])
        
        controls_valid = all([control.validate() for control in controls])
        
        if controls_valid :
            
            shape = [x.GetWindow().value 
                     for x in self.ui.sizes_sizer.GetChildren()]
            shape = list(reversed(shape))
            
            can_load = medipy.io.raw.can_load(self.ui.filename.value,
                shape, self._dtype, self._endianness, self.ui.offset.value,
                self.ui.skip_end.value)
            
            if can_load :
                self.ui.message.SetLabel("Parameters match the file size")
                self.ui.message.SetForegroundColour(wx.GREEN)
                self.ui.import_.Enable()
            else :
                format = medipy.io.raw.get_struct_format(
                     shape, self._dtype, self._endianness, self.ui.offset.value)
                size = struct.calcsize(format)
                
                message = "Parameters match the file size ({0} vs. {1})".format(
                    size, self._filesize)
                self.ui.message.SetLabel(message)
                self.ui.message.SetForegroundColour(wx.RED)
                self.ui.import_.Disable()
        else :
            self.ui.message.SetLabel("")
            self.ui.import_.Disable()
    
if __name__ == "__main__" :
    import wx
    
    app = wx.PySimpleApp()
    
    raw_dialog = ImportRawDialog(None)
    if raw_dialog.ShowModal() == wx.ID_OK :
        print raw_dialog.image
    else :
        print "canceled"
