##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import os
import xml.dom.minidom as md

import numpy
import wx
import wx.xrc

from medipy.base import find_resource, ImageAnnotation, ObservableList
import medipy.io
import medipy.io.rbnmr as rbnmr
import medipy.gui.xrc_wrapper

class SpectroDialog(medipy.gui.xrc_wrapper.Dialog):
    """ Dialog allowing the user to choose a spectroscopy image
        within a directory.
        Il also gives the choice either to open a 1D spectrum or compute
        the projection histogram.
    """
    def __init__(self, parent=None, *args, **kwargs):
    
        resource = wx.xrc.EmptyXmlResource()
        resource.InsertHandler(medipy.gui.xrc_wrapper.DirectoryXMLHandler())
        resource.InsertHandler(medipy.gui.xrc_wrapper.FileXMLHandler())
        
        file = open(find_resource("resources/gui/spectro_dialog.xrc"))
        resource.LoadFromString(file.read())
        
        dialog = resource.LoadDialog(parent, "main_dialog")
        medipy.gui.xrc_wrapper.Dialog.__init__(self, dialog, *args, **kwargs)
        
        controls = ["dir_dialog", "dir_listbox", "image_listbox",
                    "reference_listbox", "annotations_listbox", "annotations_checkbox",
                     "file_dialog", "open_button", "cancel_button"]
        
        for control in controls : 
            setattr(self, "_"+control, wx.xrc.XRCCTRL(self, control))
            
        self.SetTitle("Load spectroscopy image")    
        self._open_button.Disable()
            
        # Attributes initialization    
        self._patient_dirname = None
        self._image_dict = {}
        self._rect_dict = {}
        self._image_path = None
        self._reference_path = None
        self._annotations_path = None
        self._file_dialog._wildcard = "annotation file|peaklist.xml"
        self._file_dialog._button.Disable()
        
        # Events
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        self._dir_listbox.Bind(wx.EVT_LISTBOX, self.OnDirChosen)
        self._image_listbox.Bind(wx.EVT_LISTBOX, self.OnImageChosen)
        self._reference_listbox.Bind(wx.EVT_LISTBOX, self.OnReferenceChosen)
        self._annotations_listbox.Bind(wx.EVT_LISTBOX, self.OnAnnotationChosen)
        self._annotations_checkbox.Bind(wx.EVT_CHECKBOX, self.OnCustomChecked)
        self._open_button.Bind(wx.EVT_BUTTON, self.OnOpenClicked)
        self._cancel_button.Bind(wx.EVT_BUTTON, self.OnCancelClicked)

        self._dir_dialog.add_observer("value", self.OnPathChanged)
        self._file_dialog.add_observer("value", self.OnCustomPathChanged)
      
    def update_information(self):
        """ Enable the open button if all the necessary pieces of information
            have been gathered. Disable it if not.
        """
        if self._image_path is not None:
            self._open_button.Enable() 
        else:
            self._open_button.Disable() 
        
        self.Fit()
        self.GetSizer().SetSizeHints(self)
        
    ##########
    # Events #
    ##########  
      
    def OnClose(self, event):
        """ Shut the window
        """
        self.Destroy()
        
    def OnPathChanged(self, event):
        """ Set up the directories and reference spectra listboxes
        """
        # Clean all the listboxes
        self._dir_listbox.Clear()
        self._image_listbox.Clear()
        self._reference_listbox.Clear()
        self._annotations_listbox.Clear()
        
        # Set up the directories listbox  
        self._patient_dirname = self._dir_dialog._text.GetValue()      
        dir_list = []
        for dirpath, dirnames, filenames in os.walk(self._patient_dirname):
            dir_list.append((dirpath, dirnames, filenames))
        dir_list[0][1].sort()
        self._dir_listbox.InsertItems(dir_list[0][1], 0)
        
        # Define the path splitter
        if '/' in self._patient_dirname:
            splitter = '/'
        else:
            splitter = '\\'
        
        # Set up the reference spectra and annotations listboxes
        self._ref_dict ={}
        self._annotations_dict ={}

        for i in dir_list:
            for filename in i[2]:
                # A reference spectrum has been found
                if filename in ['1r', '1i']:
                    self._ref_dict[i[0].split(splitter)[-3] +'->'+ filename] = i[0]
                    
                # An annotation file has been found
                if filename == "peaklist.xml":
                    self._annotations_dict[i[0].split(splitter)[-3] +'->'+ filename] = os.path.join(i[0], "peaklist.xml")
        
        sorted_spectra = self._ref_dict.keys()
        sorted_spectra.sort()
        self._reference_listbox.InsertItems(sorted_spectra, 0)
        self._reference_listbox.Insert("None", 0)
        
        if self._annotations_dict.keys() != []:
            sorted_annotations = self._annotations_dict.keys()
            sorted_annotations.sort()
            self._annotations_listbox.InsertItems(sorted_annotations, 0)
        self._annotations_listbox.Insert("None", 0)
            
        self.update_information()
        
    def OnCustomPathChanged(self, event):
        """ Set up the path to the custom annotation file
        """
        self._annotations_path = self._file_dialog.value
        
    def OnDirChosen(self, event):
        """ Display the available images within the selected directory
        """
        self._image_listbox.Clear()
        image_list = []
        for dirpath, dirnames, filenames in os.walk(os.path.join(self._patient_dirname, self._dir_listbox.GetStringSelection())):
            image_list.append((dirpath, dirnames, filenames))
        
        # Dictionary associating a filename with its parent directory    
        self._image_dict = {}
        for i in image_list:
            for filename in i[2]:
                if filename in ['2rr','2ri','2ir','2ii']:
                    self._image_dict[filename]=i[0]
                    
        sorted_images = self._image_dict.keys()
        sorted_images.sort()
        sorted_images.reverse()
        self._image_listbox.InsertItems(sorted_images, 0)
        
        self.update_information()
    
    def OnImageChosen(self, event):
        """ Set the full path to the selected image
        """
        if self._image_listbox.GetStringSelection() != '':
            self._image_path = os.path.join(self._image_dict[self._image_listbox.GetStringSelection()],self._image_listbox.GetStringSelection())
        
        self.update_information()
    
    def OnReferenceChosen(self, event):
        """Set the full path to the selected reference spectrum

        """
        if (self._reference_listbox.GetStringSelection() != "None") and (self._reference_listbox.GetStringSelection() != ''):
            self._reference_path = os.path.join(self._ref_dict[self._reference_listbox.GetStringSelection()], self._reference_listbox.GetStringSelection()[-2:])
        else:
            self._reference_path = None
        
        self.update_information()
        
    def OnAnnotationChosen(self, event):
        """ Set the full path to the selected annotation file
        """
        if (self._annotations_listbox.GetStringSelection() != "None") and (self._annotations_listbox.GetStringSelection() != ''):
            self._annotations_path = self._annotations_dict[self._annotations_listbox.GetStringSelection()]
            
    def OnCustomChecked(self, event):   
        """ Allow the user to load an annotation file from another patient directory
        """
        if self._annotations_checkbox.IsChecked():
            self._annotations_listbox.Disable()
            self._file_dialog._button.Enable()
            if self._file_dialog.validate():
                self._annotations_path = self._file_dialog.value
            else:
                self._annotations_path = None
        else:
            self._annotations_listbox.Enable()
            self._file_dialog._button.Disable()
            if self._annotations_listbox.IsEmpty():
                self._annotations_path = None
            elif self._annotations_listbox.GetStringSelection() != '':
                self._annotations_path = self._annotations_dict[self._annotations_listbox.GetStringSelection()]
        
    def OnOpenClicked(self, event):
        """ Load the spectrum with either a reference spectrum or a computed histogram
        """
        
        # Create the image
        image = medipy.io.load(self._image_path) #, 0, loader_class= nmr2D.Nmr2D)
        
        # Insert a reference spectrum into the image if one has been specified
        if self._reference_path  is not None:
            spectrum = numpy.fromfile(self._reference_path, numpy.int32)
            image.metadata["header"]["proton_spectrum"] = spectrum
        
        # Load a list of annotations if an annotation file has been specified
        if self._annotations_path is not None:
            image.metadata["Data"] = image.data
            image.metadata["annotations"] = []
            dom = md.parse(self._annotations_path)
            peaks = dom.getElementsByTagName("Peak2D")
            for peak in peaks:
                annotation = ImageAnnotation()
                ppm = (float(peak.getAttribute("F1")),float(peak.getAttribute("F2")))
                point = rbnmr.ppm_to_point(ppm, 
                    image.metadata["header"]["Procs"],
                    image.metadata["header"]["Proc2s"])
                annotation.position = [0, point[-2], point[-1]]
                annotation.label = peak.getAttribute("annotation")
                annotation.shape = ImageAnnotation.Shape.cross
                annotation.size = 10
                annotation.color = [0, 1., 0.]
                annotation.filled = False
                annotation.depth = 10    
                image.metadata["annotations"].append(annotation)  
                
            image.annotations = ObservableList(image.metadata["annotations"])
        
        wx.GetApp().append_image(image)
        
        # Close the window
        self.Destroy()
        
    def OnCancelClicked(self, event):
        """ Abort
        """
        self.OnClose(event)
        
if __name__ == "__main__" :
    app = wx.App()
    
    dlg = SpectroDialog()
    dlg.ShowModal()
    dlg.GetSizer().SetSizeHints(dlg)
    
    app.MainLoop()