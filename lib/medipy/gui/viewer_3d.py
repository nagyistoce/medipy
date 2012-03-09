##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import numpy
from vtk import *
from math import *
import wx

import medipy.base
from medipy.base import Observable, ObservableList
from medipy.gui import wxVTKRenderWindowInteractor

import viewer_3d_tools

class Viewer3D(wx.Panel, Observable):
    
    def __init__(self, *args, **kwargs):
        
        wx.Panel.__init__(self, *args, **kwargs)
        Observable.__init__(self, ["active_object"])
        
        self._objects_3d = ObservableList()
        self._active_object = None
        
        # Currently pressed button, None if no button is pressed
        self._button = None
        
        # Bindings from mouse button to tool
        self.tools = {
            "LeftButton" : viewer_3d_tools.Rotate(self), 
            "MiddleButton" : viewer_3d_tools.Pan(self),
            "RightButton" : viewer_3d_tools.Dolly(self),
            "MouseWheelForward" : viewer_3d_tools.WheelDolly(self, "forward"),
            "MouseWheelBackward" : viewer_3d_tools.WheelDolly(self, "backward")
        }
        
        # Bindings from key to tool
        self.key_bindings = {
            # Nothing by default
        }
        
        # VTK objects for rendering
        self._renderer = vtkRenderer()
        self._renderer.SetBackground(82./255.,87./255.,110./255.)
        self._renderer.GetActiveCamera().ParallelProjectionOff()
        
        self._rwi = wxVTKRenderWindowInteractor(self, wx.ID_ANY)
        self._rwi.GetRenderWindow().AddRenderer(self._renderer)
        self._sizer = wx.BoxSizer()
        self._sizer.Add(self._rwi, 1, wx.EXPAND)
        self.SetSizer(self._sizer)
        
        # Interaction styles : customize events
        self._rwi.SetInteractorStyle(None)
        
        self._rwi.AddObserver("LeftButtonPressEvent", self._start_interaction)
        self._rwi.AddObserver("MiddleButtonPressEvent", self._start_interaction)
        self._rwi.AddObserver("RightButtonPressEvent", self._start_interaction)
        
        self._rwi.AddObserver("LeftButtonReleaseEvent", self._stop_interaction)
        self._rwi.AddObserver("MiddleButtonReleaseEvent", self._stop_interaction)
        self._rwi.AddObserver("RightButtonReleaseEvent", self._stop_interaction)
        
        self._rwi.AddObserver("MouseWheelForwardEvent", self._start_interaction)
        self._rwi.AddObserver("MouseWheelBackwardEvent", self._start_interaction)
        
        self._rwi.AddObserver("MouseMoveEvent", self._dispatch_interaction)
        
        self._rwi.AddObserver("KeyPressEvent", self._keypress)
        self._rwi.AddObserver("KeyReleaseEvent", self._keyrelease)
        
        # Clean up when the window is closed
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        
        self.objects_3d.add_observer("any", self._on_objects_3d_modified)
    
    
    def view_all(self):
        self.reset_camera()
        self._rwi.Render()
        
        
    def get_bounds(self):    
        # get the bounding box of all of the objects of the current scene (visible or not)
        bounds = []
        
        for _object in self._objects_3d :
            actor = _object.actor
            bounds.append(actor.GetBounds())
       
        Xmin = 0
        Xmax = 0
        Ymin = 0
        Ymax = 0
        Zmin = 0
        Zmax = 0
        
        if len(bounds) > 0 :
            Xmin = min(map(lambda x:x[0], bounds))
            Xmax = max(map(lambda x:x[1], bounds))
            Ymin = min(map(lambda x:x[2], bounds))
            Ymax = max(map(lambda x:x[3], bounds))
            Zmin = min(map(lambda x:x[4], bounds))
            Zmax = max(map(lambda x:x[5], bounds))
            
        return ((Xmin, Ymin, Zmin), (Xmax, Ymax, Zmax))
        

    def _get_distance(self):
        """get the greatest dimension of the enclosing bounding box of all the element in the scene\
           (visible or not)
        """
        
        Pmin, Pmax = self.get_bounds()
        
        w1 = Pmax[0] - Pmin[0]
        w2 = Pmax[1] - Pmin[1]
        w3 = Pmax[2] - Pmin[2]
        
        return max(w1, w2, w3)
    
    
    def reset_camera(self):
        # view angle
        view_angle = 53.13
        self._renderer.GetActiveCamera().SetViewAngle(view_angle)
        
        Pmin, Pmax = self.get_bounds()
        
        # target : focal point
        focal_point = numpy.add(Pmin, Pmax)/2.
        
        distance = self._get_distance()
        
        position = self._renderer.GetActiveCamera().GetPosition()
        v = numpy.subtract(position, focal_point)
        
        norme = numpy.linalg.norm(v)
        
        remoteness = 1.5
        factor = (distance*remoteness)/norme
        new_pos = numpy.add(numpy.multiply(v, factor), focal_point)
        
        # update the camera ...
        self._renderer.GetActiveCamera().SetFocalPoint(*focal_point)
        
        self._renderer.GetActiveCamera().SetPosition(new_pos)
        
        self._renderer.ResetCameraClippingRange(Pmin[0], Pmax[0], Pmin[1], Pmax[1], Pmin[2], Pmax[2])
    
    ##############
    # Properties #
    ##############
    
    def _set_objects_3d(self, objects_3d):
        collection = self._renderer.GetViewProps()
        collection.InitTraversal()
        for i in range(collection.GetNumberOfItems()) :
            actor = collection.GetNextProp()
            self._renderer.RemoveActor(actor)
        
        self._objects_3d = objects_3d
        for objet in self.objects_3d :
            self._renderer.AddActor(objet.actor)
        
        self._objects_3d.remove_observer("any", self._on_objects_3d_modified)
        self._objects_3d = objects_3d
        self.objects_3d.add_observer("any", self._on_objects_3d_modified)
        
    
    def _set_active_object(self, object):
        if object not in self.objects_3d :
            raise medipy.base.Exception("Object is not in viewer's objects")
        
        self._active_object = object
        self.notify_observers("active_object")
    
    objects_3d = property(lambda x:x._objects_3d, _set_objects_3d)
    active_object = property(lambda x:x._active_object, _set_active_object)
    renderer = property(lambda x:x._renderer)
    render_window_interactor = property(lambda x:x._rwi)
    
    ##################
    # Event handlers #
    ##################
    
    def OnClose(self, event) :
        self._renderer.RemoveAllViewProps()
        del self._renderer
        self._rwi.GetRenderWindow().Finalize()
        self._rwi.SetRenderWindow(None)
        self._rwi.Close()
        del self._rwi
        self.Destroy()
    
    
    def _on_objects_3d_modified(self, event):
        if event.event == "append" :
            for objet in self.objects_3d :
                if not self._renderer.HasViewProp(objet.actor) :
                    self._renderer.AddActor(objet.actor)
        elif event.event == "delete_item" :
            self._renderer.RemoveActor(event.old_value.actor)
        else :
            raise medipy.base.Exception("Unmanaged event : %s" % event.event)
        
        self._rwi.Render()
    
    
    def _start_interaction(self, object, event):
        if event.startswith("MouseWheel") :
            wheel_event = event[:-len("Event")]
            if wheel_event in self.tools :
                self.tools[wheel_event].start_interaction()
                self.tools[wheel_event].stop_interaction()
        elif event.endswith("PressEvent") :
            button = event[:-len("PressEvent")]
            
            event = button
            if object.GetControlKey() :
                event += "Control"
            if object.GetShiftKey() :
                event += "Shift"
            
            # In case no modifier is configured, fall back to the basic behavior
            if event not in self.tools :
                event = button
            
            if event in self.tools :
                self._button = button
                self.tools[event].start_interaction()

    def _stop_interaction(self, object, event):
        button = event[:-len("ReleaseEvent")]
        
        event = self._button
        if object.GetControlKey() :
            event += "Control"
        if object.GetShiftKey() :
            event += "Shift"
        
        for key, tool in self.tools.items() :
            if key.startswith(button) : 
                tool.stop_interaction()
        
        self._button = None
    
    def _dispatch_interaction(self, object, event):
        event = self._button
        
        if event is None :
            return
        
        if object.GetControlKey() :
            event += "Control"
        if object.GetShiftKey() :
            event += "Shift"
        
        # In case no modifier is configured, fall back to the basic behavior
        if event not in self.tools :
            event = self._button
        if event in self.tools :
            self.tools[event].dispatch_interaction()
    
    def _keypress(self, object, event):
        ordinal = ord(object.GetKeyCode()) 
        if ordinal in self.key_bindings :
            self.key_bindings[ordinal].press()
    
    def _keyrelease(self, object, event):
        ordinal = ord(object.GetKeyCode()) 
        if ordinal in self.key_bindings :
            self.key_bindings[ordinal].release()
    
    #####################
    # Private interface #
    #####################
    
    def _update_actors(self):
        actors_from_renderer = []
        collection = self._renderer.GetViewProps()
        collection.InitTraversal()
        for i in range(collection.GetNumberOfItems()) :
            actor = collection.GetNextProp()
            actors_from_renderer.append(actor)
        
        actors_from_objects = [object.actor for object in self.objects_3d]
        
        for actor in actors_from_renderer :
            if actor not in actors_from_objects :
                self._renderer.RemoveActor(actor)
        
        for object in self.objects_3d :
            if not self._renderer.HasViewProp(object.actor) :
                self._renderer.AddActor(object.actor)
                
                
        

if __name__ == "__main__" :
    import sys
    
    if len(sys.argv) != 2 :
        print "Syntax : %s <filename>"%sys.argv[0]
        sys.exit(1)
    
    from vtk import *
    import wx
    
    from medipy.base import Object3D
    
    # Build basic application
    app = wx.App()
    
    frame = wx.Frame(None, size=(800,800))
    frame.Show()
    app.SetTopWindow(frame)
    
    sizer = wx.BoxSizer()
    frame.SetSizer(sizer)
    sizer.SetSizeHints(frame)
    
    # Read the object
    reader = vtkPolyDataReader()
    reader.SetFileName(sys.argv[1])
    object = Object3D(reader.GetOutput(), sys.argv[1])
    
    viewer_3d = Viewer3D(frame)
    viewer_3d.objects_3d.append(object)
    sizer.Add(viewer_3d, 1, wx.EXPAND)
    sizer.Layout()
    
    app.MainLoop()