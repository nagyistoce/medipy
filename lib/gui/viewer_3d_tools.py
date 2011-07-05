##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import medipy

import numpy
from vtk import *
import wx

class Tool(object):
    def __init__(self, viewer_3d):
        self._viewer_3d = viewer_3d
    
    def select(self) :
        """ Action to perform when the user selects this tool """
        pass
    
    def deselect(self) :
        """ Action to perform when the user deselects this tool """ 
        pass
    
    def start_interaction(self) :
        """ Action to perform when the user clicks on the mouse button """
        pass
    
    def dispatch_interaction(self) :
        """ Action to perform when the user moves the mouse with the button 
            pressed 
        """
        pass
    
    def stop_interaction(self) :
        """ Action to perform when the user releases the mouse button """
        pass
    
    def _pick_point(self, position, collection=None):
        picker = self._viewer_3d.render_window_interactor.GetPicker()
        
        if collection is not None :
            on_something = picker.PickProp(position[0], position[1], 
                                           self._viewer_3d.renderer, collection)
        else : 
            on_something = picker.PickProp(position[0], position[1], 
                                           self._viewer_3d.renderer)
        
        if on_something :
            path = []
            picker.GetPath().InitTraversal()
            for i in range(picker.GetPath().GetNumberOfItems()) :
                path.append(picker.GetPath().GetNextNode().GetViewProp())
            return (picker.GetPickPosition(), path)
        else : 
            return (None, [])

class Draw(object):
    """ Add points to the selected drawings, based on other objects
    """
    
    def __init__(self, viewer_3d):
        self._viewer_3d = viewer_3d
    
    def select(self):
        pass
    
    def deselect(self):
        pass
    
    def start_interaction(self):
        self.draw_point()
    
    def stop_interaction(self):
        pass
    
    def dispatch_interaction(self):
        self.draw_point()
    
    def draw_point(self):
        selected_object = viewer_3d.selected_object
        if selected_object is None :
            return
        
        rwi = self._viewer_3d.render_window_interactor
        
        picker = vtkPropPicker()
        on_something = picker.PickProp(rwi.GetEventPosition()[0], rwi.GetEventPosition()[1], self._renderer)
        
        if on_something :
            picked_point = picker.GetPickPosition()
            
            point_id = drawing.dataset.GetPoints().InsertNextPoint(picked_point)

            drawing.dataset.Modified()
            
            rwi.Render()

class Rotate(Tool):
    def __init__(self, viewer_3d):
        Tool.__init__(self, viewer_3d)
    
    def dispatch_interaction(self):
        rwi = self._viewer_3d.render_window_interactor
        # Code taken from VTK/Examples/GUI/Python/CustomInteraction.py
        # Similar code available in vtkInteractorStyleTrackballCamera::Rotate
        camera = self._viewer_3d.renderer.GetActiveCamera()
        camera.Azimuth(rwi.GetLastEventPosition()[0] - rwi.GetEventPosition()[0])
        camera.Elevation(rwi.GetLastEventPosition()[1] - rwi.GetEventPosition()[1])
        camera.OrthogonalizeViewUp()
    
        rwi.Render()

class Pan(Tool):
    def __init__(self, viewer_3d):
        Tool.__init__(self, viewer_3d)
    
    def dispatch_interaction(self):
        # Code taken from VTK/Examples/GUI/Python/CustomInteraction.py
        # Similar code available in vtkInteractorStyleTrackballCamera::Pan
        
        camera = self._viewer_3d.renderer.GetActiveCamera()
        focal_point = camera.GetFocalPoint()
        position = camera.GetPosition()
        
        self._viewer_3d.renderer.SetWorldPoint(focal_point[0], focal_point[1], focal_point[2], 1.0)
        self._viewer_3d.renderer.WorldToDisplay()
        display_point = self._viewer_3d.renderer.GetDisplayPoint()
        focal_depth = display_point[2]
        
        rwi = self._viewer_3d.render_window_interactor
        center = numpy.divide(rwi.GetRenderWindow().GetSize(), 2.)
        a_point = numpy.add(center, 
                            numpy.subtract(rwi.GetEventPosition(), 
                                           rwi.GetLastEventPosition()))
        
        self._viewer_3d.renderer.SetDisplayPoint(a_point[0], a_point[1], focal_depth)
        self._viewer_3d.renderer.DisplayToWorld()
        r_point = list(self._viewer_3d.renderer.GetWorldPoint())

        if r_point[3] != 0.0:
            r_point = [x/r_point[3] for x in r_point[:3]]
    
        focal_point = numpy.add(numpy.subtract(focal_point, r_point[:3])/2., focal_point)
        camera.SetFocalPoint(*focal_point)
        
        position = numpy.add(numpy.subtract(focal_point, r_point[:3])/2., position)
        camera.SetPosition(*position)
    
        rwi.Render()

class Dolly(Tool) :
    def __init__(self, viewer_3d):
        Tool.__init__(self, viewer_3d)
        
        self._motion_factor = 0.5
        
    def dispatch_interaction(self, factor=None):
        # Code taken from VTK/Examples/GUI/Python/CustomInteraction.py
        # Similar code available in vtkInteractorStyleTrackballCamera::Dolly
        rwi = self._viewer_3d.render_window_interactor
        
        if factor is None :
            dy = rwi.GetEventPosition()[1] - rwi.GetLastEventPosition()[1]
            dyf = self._motion_factor * float(dy)
            
            factor = 1.02**dyf
        
        camera = self._viewer_3d.renderer.GetActiveCamera()
        if camera.GetParallelProjection() :
            camera.SetParallelScale(camera.GetParallelScale() / factor)
        else :
            camera.Dolly(factor)
            self._viewer_3d.renderer.ResetCameraClippingRange()
      
        rwi.Render()

class WheelDolly(Tool):
    def __init__(self, viewer_3d, direction):
        Tool.__init__(self, viewer_3d) 
        self._direction = 1.3 if direction=="forward" else 1/1.3
        self._dolly = Dolly(viewer_3d)
    
    def start_interaction(self):
        self._dolly.dispatch_interaction(self._direction)
    
class Pick(Tool):
    def __init__(self, viewer_3d):
        Tool.__init__(self, viewer_3d)
    
    def start_interaction(self):
        rwi = self._viewer_3d.render_window_interactor
        
        picker = vtkPropPicker()
        on_something = picker.PickProp(rwi.GetEventPosition()[0], rwi.GetEventPosition()[1], self._viewer_3d.renderer)
        
        if on_something :
            
            actor = picker.GetActor()
            for object in self._viewer_3d.objects_3d :
                if actor is object.actor :
                    # Pick in the associated image, if it exists
                    if object.gui_image is not None : 
                        picked_point = numpy.asarray(picker.GetPickPosition()).round()
                        picked_point = [ x for x in reversed(picked_point)]
                        object.gui_image.position = picked_point
                        object.gui_image.render()
                    # Update the GUI
                    self._viewer_3d.active_object = object

# tool for distance calculation on the 3d objects on the scene (geodesic, euclidean ...)
class Distance(Tool):
    def __init__(self, viewer_3d):
        Tool.__init__(self, viewer_3d)
        
        self._start_annotation_sphere = self._create_annotation_sphere(1.5, (0, 1, 0))
        self._start_annotation_sphere.VisibilityOff()
        self._end_annotation_sphere = self._create_annotation_sphere(1.5, (1, 0, 0))
        self._end_annotation_sphere.VisibilityOff()
        
        self._selected_sphere = None
        
        self._picker = medipy.vtk_addons.vtkModifiedPropPicker()
        self._props = None
        self._locators = {}
        
        self._start_object = None
        self._end_object = None
        
        self._leader = vtkLeaderActor2D()
        self._leader.GetPositionCoordinate().SetCoordinateSystemToWorld()
        self._leader.GetPosition2Coordinate().SetCoordinateSystemToWorld()
        self._leader.SetArrowStyleToFilled()
        self._leader.GetProperty().SetColor(0, 0, 0) #black
        self._leader.VisibilityOff()
        
        self._viewer_3d.renderer.AddActor(self._start_annotation_sphere)
        self._viewer_3d.renderer.AddActor(self._end_annotation_sphere)
        self._viewer_3d.renderer.AddActor(self._leader)
        self._viewer_3d.render_window_interactor.Render()
    
    
    def select(self):
        self._dialog = wx.Dialog(self._viewer_3d, -1, "Distance type", style=wx.CAPTION)
        label = wx.StaticText(self._dialog, label="Distance type : ")
        self._choice = wx.Choice(self._dialog, -1, choices=("Euclidean", "Geodesic", "Automatic"))
        self._choice.Bind(wx.EVT_CHOICE, self._on_choice_changed)
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(label, flag=wx.ALIGN_CENTER)
        sizer.Add(self._choice, 0, wx.ALIGN_CENTER|wx.ALL, 5)
        self._dialog.SetSizer(sizer)
        sizer.SetSizeHints(self._dialog)
        self._dialog.Show()
        
    
    def deselect(self):
        self._viewer_3d.renderer.RemoveActor(self._start_annotation_sphere)
        self._viewer_3d.renderer.RemoveActor(self._end_annotation_sphere)
        self._viewer_3d.renderer.RemoveActor(self._leader)
        self._dialog.Hide()
        self._dialog.Destroy()
        self._viewer_3d.render_window_interactor.Render()
    
    
    def start_interaction(self):
        rwi = self._viewer_3d.render_window_interactor
        rend = self._viewer_3d.renderer
        
        self._props = vtkPropCollection()
        actors = rend.GetActors()
        actors.InitTraversal()
        for i in range(actors.GetNumberOfItems()) :
            actor = actors.GetNextItem()
            for _object in self._viewer_3d.objects_3d :
                if actor is _object.actor :
                    self._props.AddItem(actor)
                    
        screen_position = rwi.GetEventPosition()
        picked_point, path = self._pick_point(screen_position)
        
        if picked_point is not None : 
            if prop in [self._start_annotation_sphere, self._end_annotation_sphere] :
                if prop == self._start_annotation_sphere :
                    self._selected_sphere = 1
                elif prop == self._end_annotation_sphere :
                    self._selected_sphere = 2
            else :
                self._selected_sphere = None
                if not self._start_annotation_sphere.GetVisibility() :
                    self._start_annotation_sphere.VisibilityOn()
                    self._start_annotation_sphere.SetPosition(picked_point)
                    self._start_object = prop
                else :
                    if not self._end_annotation_sphere.GetVisibility() :
                        self._end_annotation_sphere.VisibilityOn()
                        self._end_annotation_sphere.SetPosition(picked_point)
                        self._end_object = prop
            
            if self._start_annotation_sphere.GetVisibility() and self._end_annotation_sphere.GetVisibility() :
                self._display_distance()
                        
        self._viewer_3d.render_window_interactor.Render()
                        
    def dispatch_interaction(self):
        
        screen_position = self._viewer_3d.render_window_interactor.GetEventPosition()
        picked_point, path = self._pick_point(screen_position, self._props)
        
        if picked_point is not None :
            if self._selected_sphere == 1 :
                self._start_annotation_sphere.SetPosition(picked_point)
                self._start_object = prop
            elif self._selected_sphere == 2 :
                self._end_annotation_sphere.SetPosition(picked_point)
                self._end_object = prop
            
            if self._start_annotation_sphere.GetVisibility() and self._end_annotation_sphere.GetVisibility() :
                self._display_distance()
                    
        self._viewer_3d.render_window_interactor.Render()
        

    def _on_choice_changed(self, event):
        if self._start_annotation_sphere.GetVisibility() and self._end_annotation_sphere.GetVisibility() :
            self._display_distance()
        self._viewer_3d.render_window_interactor.Render()
        
        
    def _display_distance(self):
        a = self._start_annotation_sphere.GetCenter()
        b = self._end_annotation_sphere.GetCenter()
        
        if self._choice.GetStringSelection() == "Euclidean" :
            distance = self._compute_euclidean_distance(a, b)
            color = (1, 1, 0) # yellow
        elif self._choice.GetStringSelection() == "Geodesic" :
            if self._start_object == self._end_object :
                distance = self._compute_geodesic_distance(a, b)
            else :
                distance = "Invalid"
            color = (0, 1, 1) # cyan
        elif self._choice.GetStringSelection() == "Automatic" :
            if self._start_object == self._end_object :
                # compute geodesic distance
                distance = self._compute_geodesic_distance(a, b)
                color = (0, 1, 1)
            else :
                # compute euclidean distance
                distance = self._compute_euclidean_distance(a, b)
                color = (1, 1, 0)
                
        self._leader.VisibilityOn()
        self._leader.GetProperty().SetColor(color)
        self._leader.GetPositionCoordinate().SetValue(a)
        self._leader.GetPosition2Coordinate().SetValue(b)
        if distance == "Invalid" :
            self._leader.SetLabel("Invalid")
        else :
            self._leader.SetLabel("%.2f" % distance)
    
    
    def _compute_euclidean_distance(self, a, b):
        v = numpy.subtract(b, a)
        return numpy.linalg.norm(v)
    
    
    def _compute_geodesic_distance(self, a, b):
        for _object in self._viewer_3d.objects_3d :
                if self._start_object is _object.actor :
                    _dataset = _object.dataset
        geo_path = vtkDijkstraGraphGeodesicPath()
        geo_path.SetInput(_dataset)
            
        if not self._locators.has_key(self._start_object) :
            _locator = vtkPointLocator()
            _locator.SetDataSet(_dataset)
            _locator.BuildLocator()
            self._locators[self._start_object] = _locator
                                
        locator = self._locators[self._start_object]
        start_point = locator.FindClosestPoint(a)
        end_point = locator.FindClosestPoint(b)
            
        geo_path.SetStartVertex(start_point)
        geo_path.SetEndVertex(end_point)
            
        #compute geodesic path boundaries
        start_point_coords = _dataset.GetPoint(start_point)
        end_point_coords = _dataset.GetPoint(end_point)
           
        v = numpy.subtract(start_point_coords, a)
        start_boundary = numpy.linalg.norm(v)
        v = numpy.subtract(b, end_point_coords)
        end_boundary = numpy.linalg.norm(v)
            
        return self._geodesic_distance(geo_path) + start_boundary + end_boundary
    
    
    def _create_annotation_sphere(self, radius, color):
        sphere = vtkSphereSource()
        sphere.SetPhiResolution(16)
        sphere.SetThetaResolution(16)
        sphere.SetRadius(radius)

        sphere_mapper = vtkPolyDataMapper()
        sphere_mapper.SetInputConnection(sphere.GetOutputPort())

        actor = vtkActor()
        actor.SetMapper(sphere_mapper)
        actor.GetProperty().SetColor(color)
    
        return actor
    
    
    def _geodesic_distance(self, geodesic_path):
        geodesic_path.Update()
        data = geodesic_path.GetOutput()
        points = data.GetPoints()
        n = points.GetNumberOfPoints()
        geodesic_path = 0
        for i in range(n-1):
            a = data.GetPoint(i)
            b = data.GetPoint(i+1)
            v = numpy.subtract(b, a)
            path = numpy.linalg.norm(v)
            geodesic_path = geodesic_path + path        
        
        return geodesic_path