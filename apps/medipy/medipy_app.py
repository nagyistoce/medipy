# encoding: utf-8
##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import os
import sys

from vtk import (vtkDataReader, vtkPolyDataReader, vtkPolyDataWriter,
    vtkRectilinearGridReader, vtkStructuredGridReader, 
    vtkUnstructuredGridReader, vtkVRMLImporter,)
import wx

import medipy.gui.base
from medipy.gui.image.cine_dialog import CineDialog
from medipy.gui.function_gui_builder import FunctionGUIBuilder

from medipy.base import ObservableList, Object3D

from main_frame import MainFrame
import menu_builder

class MediPyApp(medipy.gui.base.Application) :
    
    _application_name = "MediPy"
    _vendor_name = "Université de Strasbourg"
    
    def __init__(self, *args, **kwargs):
        
        # Public interface
        self.viewer_3ds = ObservableList()
        
        # Private interface
        self._frame = None
        self._full_screen = False
        
        super(MediPyApp, self).__init__(*args, **kwargs)
    
    ####################
    # Public interface #
    ####################
        
    def insert_viewer_3d(self, index, viewer_3d):
        self.viewer_3ds.insert(index, viewer_3d)
        viewer_3d.add_observer("close", self._on_viewer_3d_close)
    
    def append_viewer_3d(self, viewer_3d):
        return self.insert_viewer_3d(len(self.viewer_3ds), viewer_3d)
    
    def toggle_full_screen(self):
        self._full_screen = not self._full_screen
        
        if self._full_screen : 
            self._frame.full_screen(self.active_image)
        else :
            self._frame.full_screen(None)      
    
    def execute_script(self, filename):
        execfile(filename, globals(), locals())
    
    def quit(self):
        for viewer in self.viewer_3ds :
            viewer.Close()
    
    def load_object_3d(self, path, viewer_3d):
        generic_reader = vtkDataReader()
        generic_reader.SetFileName(path)
        
        if generic_reader.OpenVTKFile() and generic_reader.ReadHeader() :
            if generic_reader.IsFileStructuredPoints() :
                raise Exception("Cannot read VTK structured points")
            elif generic_reader.IsFilePolyData() :
                reader = vtkPolyDataReader()
                reader.SetFileName(path)
                object_3d = Object3D(reader.GetOutput(), path)
            elif generic_reader.IsFileStructuredGrid() :
                reader = vtkStructuredGridReader()
                reader.SetFileName(path)
                object_3d = Object3D(reader.GetOutput(), path)
            elif generic_reader.IsFileUnstructuredGrid() :
                reader = vtkUnstructuredGridReader()
                reader.SetFileName(path)
                object_3d = Object3D(reader.GetOutput(), path)
            elif generic_reader.IsFileRectilinearGrid() :
                reader = vtkRectilinearGridReader()
                reader.SetFileName(path)
                object_3d = Object3D(reader.GetOutput(), path)
            else : 
                raise Exception("Cannot read VTK file containing type %i"%generic_reader.GetFileType())
            viewer_3d.objects_3d.append(object_3d)
        else :
            importer = vtkVRMLImporter()
            importer.SetFileName(path)
            importer.Update()
            actors = importer.GetRenderer().GetActors()
            number_of_actors = actors.GetNumberOfItems()
            actors.InitTraversal()
            for i in range(number_of_actors) :
                actor = actors.GetNextItem()
                object_3d = Object3D(path + ", %i"%i)
                object_3d.actor.SetProperty(actor.GetProperty())
                object_3d.actor.SetMapper(actor.GetMapper())
                viewer_3d.objects_3d.append(object_3d)
#        else :
#            raise Exception("Cannot load file %s of unknown type"%path)
    
    def save_object_3d(self, object_3d, path) :
        writer = vtkPolyDataWriter()
        writer.SetFileName(path)
        writer.SetFileTypeToBinary()
        writer.SetHeader(object_3d.name)
        writer.SetInput(object_3d.dataset)
        writer.Update()
    

    
    ##################
    # Event handlers #
    ##################
    def OnInit(self) :
        self.SetAppName(self._application_name)
        self.SetVendorName(self._vendor_name)
        
        if hasattr(sys, "frozen") :
            directory = os.path.dirname(medipy.__file__)
            os.environ["MEDIPY_PLUGINS_PATH"] = "{0}{1}{2}".format(
                directory, os.pathsep, os.environ.get("MEDIPY_PLUGINS_PATH", ""))
        
        if self.options.ensure_value("menu_file", None) is not None :
            menu = menu_builder.from_file.build_menu(self.options.menu_file)
        else :
            menu = []
            for directory in medipy.Importer().plugins_path :
                menu.extend(menu_builder.from_api.build_menu(directory))
        
        self._frame = MainFrame(menu, None, title=self._application_name, size=(1000,800))
        self._frame.Show()
        self.SetTopWindow(self._frame)
        
        return True

    def _on_viewer_3d_close(self, event):
        viewer = event.object
        index = self.viewer_3ds.index(viewer)
        del self.viewer_3ds[index]
