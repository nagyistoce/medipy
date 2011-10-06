##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import wx

import medipy.gui.image
from medipy.gui.dicom.hierarchy_tree import HierarchyTree
from medipy.gui.dicom.dataset_list import DataSetList

from medipy.io.dicom.misc import load_dicomdir_records
from medipy.io.dicom import split, sort

class Explorer(wx.Panel):
    
    def __init__(self, parent, ID=wx.ID_ANY, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=wx.TAB_TRAVERSAL, name="explorer"):
        
        self._series_position = {}
        self._number_of_frames = {}
        self._loaded_datasets = {}
        self._is_multiframe_image = False
        
        wx.Panel.__init__(self, parent, ID, pos, size, style, name)
        
        # Create widgets
        main_splitter = wx.SplitterWindow(self)
        left_splitter = wx.SplitterWindow(main_splitter)
        self._hierarchy_tree = HierarchyTree(left_splitter)
        bottom_panel = wx.Panel(left_splitter)
        self._slider = wx.Slider(bottom_panel)
        self._image = medipy.gui.image.Image(bottom_panel, "axial")
        self._dataset_list = DataSetList(main_splitter)
        
        # Setup widgets
        self._slider.Disable()
        self._image.use_origin_and_spacing = False
        self._image.orientation_visibility = False
        left_splitter.SplitHorizontally(self._hierarchy_tree, bottom_panel)
        main_splitter.SplitVertically(left_splitter, self._dataset_list) 

        # Layout
        bottom_panel.SetSizer(wx.BoxSizer(wx.VERTICAL))
        bottom_panel.GetSizer().Add(self._slider, flag=wx.EXPAND)
        bottom_panel.GetSizer().Add(self._image, 1, wx.EXPAND)
        self.SetSizer(wx.BoxSizer())
        self.GetSizer().Add(main_splitter, 1, wx.EXPAND)
        
        # Bind events
        self._hierarchy_tree.Bind(wx.EVT_TREE_SEL_CHANGED, self.OnHierarchy)
        self._slider.Bind(wx.EVT_SLIDER, self.OnSlider)
    
    def set_datasets(self, datasets):
        self._hierarchy_tree.set_datasets(datasets)
    
    ##############
    # Properties #
    ##############
    def _get_hierarchy_tree(self):
        return self._hierarchy_tree
    
    hierarchy_tree = property(_get_hierarchy_tree)
    
    ##################
    # Event handlers #
    ##################
    
    def OnHierarchy(self, dummy):
        """ Handler for hierarchy tree selection """

        if self._hierarchy_tree.selected_datasets :
            datasets = self._hierarchy_tree.selected_datasets
            
            self._setup_slider(datasets)
            self._slider.Enable()
            
            slider_value = self._series_position.get(datasets[0].series_instance_uid, 0)
            self._slider.SetValue(slider_value)
            self.OnSlider(None)
        else :
            self._slider.Disable()
    
    def OnSlider(self, dummy):
        """ Handler for slider events """
        
        self.Freeze()
        
        datasets = self._hierarchy_tree.selected_datasets
        series_uid = datasets[0].series_instance_uid
        
        loaded_datasets = self._loaded_datasets[series_uid]
        if "number_of_frames" in loaded_datasets[0] :
            dataset = loaded_datasets[0]
            frame_number = self._slider.GetValue()
            slice = slice = dataset.pixel_array[frame_number]
        else :
            dataset = loaded_datasets[self._slider.GetValue()]
            frame_number = None
            slice = dataset.pixel_array
        
        slice = medipy.base.Image(data=slice.reshape((1,)+slice.shape))
        
        if self._image.number_of_layers == 0 :
            self._image.append_layer(slice)
            self._image.reset_view()
        else :
            self._image.delete_layer(0)
            self._image.insert_layer(0, slice)
            self._image.render()
            #self._image.set_layer_image(0, slice)
            # Contrast data
        
        scroll = (self._dataset_list.GetScrollPos(wx.HORIZONTAL), 
                  self._dataset_list.GetScrollPos(wx.VERTICAL))
        self._dataset_list.set_dataset(dataset, frame_number)
        self._dataset_list.SetScrollPos(wx.HORIZONTAL, scroll[0])
        self._dataset_list.SetScrollPos(wx.VERTICAL, scroll[1])
        
        self._series_position[series_uid] = self._slider.GetValue()
        
        self.Thaw()
    
    #####################
    # Private interface #
    #####################
    
    def _setup_slider(self, datasets):
        """ Setup the image slider.
        """
        
        series_uid = datasets[0].series_instance_uid
        
        if series_uid not in self._loaded_datasets :
            loaded_datasets = [x for x in datasets if "directory_record_type" not in x]
            
            records = [x for x in datasets if "directory_record_type" in x]
            dialog_max = len(records)
            progress_dialog = wx.ProgressDialog(
                "Loading files ...", "Loading files ...", dialog_max)
            for index, record in enumerate(records) :
                loaded_datasets.extend(load_dicomdir_records([record]))
                progress_dialog.Update(index)
            progress_dialog.Destroy()
            
            images = split.images(loaded_datasets)
            sort.sort(images)
            
            self._loaded_datasets[series_uid] = images 
            
        
        loaded_datasets = self._loaded_datasets[series_uid]
        
        slider_max = 0
        for dataset in loaded_datasets :
            if "number_of_frames" in dataset :
                slider_max += dataset.number_of_frames
            else :
                slider_max += 1
        
        self._slider.SetRange(0, slider_max-1)

def main():
    import os
    
    class DicomExplorerFrame(wx.Frame):
        def __init__(self, *args, **kwargs):
            self._config = wx.Config("MediPy")
            
            wx.Frame.__init__(self, *args, **kwargs)
        
            # Widgets and layout
            self._explorer = Explorer(self)
            self.SetSizer(wx.BoxSizer())
            self.GetSizer().Add(self._explorer, 1, wx.EXPAND)
            
            # Menu Bar
            self.SetMenuBar(wx.MenuBar())
            
            # File Menu
            file_menu = wx.Menu()
            self.GetMenuBar().Append(file_menu, "&File")
            
            open_dicomdir_item = file_menu.Append(wx.ID_FILE, "O&pen DICOMDIR ...\tCtrl+O")
            self.Bind(wx.EVT_MENU, self.OnOpenDicomDir, open_dicomdir_item)
            
            import_directory_item = file_menu.Append(wx.ID_OPEN, "Import directory ...\tCtrl+I")
            self.Bind(wx.EVT_MENU, self.OnImportDirectory, import_directory_item)
            
            file_menu.AppendSeparator()
            
            quit_item = file_menu.Append(wx.ID_EXIT, "E&xit")
            self.Bind(wx.EVT_MENU, self.OnQuit, quit_item) 
        
        ##################
        # Event handlers #
        ##################
        
        def OnOpenDicomDir(self, dummy):
            file_dialog = wx.FileDialog(None, style=wx.FD_OPEN)
            file_dialog.SetDirectory(
                self._config.Read("ImageFileDialog/DefaultPath"))
            
            if file_dialog.ShowModal() == wx.ID_OK :
                datasets = [medipy.io.dicom.parse(file_dialog.GetPath())]
                self._explorer.set_datasets(datasets)
                
                self._config.Write("ImageFileDialog/DefaultPath",
                                   file_dialog.GetDirectory())
        
        def OnImportDirectory(self, dummy):
            directory_dialog = wx.DirDialog(None, style=wx.FD_OPEN)
            directory_dialog.SetPath(
                self._config.Read("ImageFileDialog/DefaultPath"))
            
            if directory_dialog.ShowModal() == wx.ID_OK :
                
                directory = directory_dialog.GetPath()
                files = [os.path.join(directory, x) for x in os.listdir(directory)]
                datasets = [medipy.io.dicom.parse(x) for x in files]
                self._explorer.set_datasets(datasets)
                
                self._config.Write("ImageFileDialog/DefaultPath", 
                                   directory_dialog.GetPath())
        
        def OnQuit(self, dummy):
            self.Close()
    
    app = wx.PySimpleApp()
    app.SetAppName("DICOM Explorer")
    
    frame = DicomExplorerFrame(None, size=(1000,900))
    frame.Show()
    
    app.MainLoop()

if __name__ == "__main__" :
    main()