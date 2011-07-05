import re
import numpy
import wx

import medipy.gui.control
from medipy.base import Image, Object3D, ObservableList
from medipy.gui import PeriodicProgressDialog, WorkerThread
from medipy.gui.viewer_3d_frame import Viewer3DFrame

from builder import parse_docstring

class FunctionGUIBuilder(object):
    """ Create the GUI of a function based on its docstring.
    
        When clicking the "Run" button, the function will be called with arguments
        extracted from the controls on this GUI.
    """
    
    def __init__(self, parent, function, images, viewer_3ds):
        
        self.panel = wx.Panel(parent)
        
        self._controls = {}
        self._function = function
        self._images = images
        self._viewer_3ds = viewer_3ds
        
        self._parameters = parse_docstring(self._function.__doc__)

        # Run & reset buttons
        self._run_button = wx.Button(self.panel, label="Run")
        reset_button = wx.Button(self.panel, label="Reset")
        buttons_sizer = wx.BoxSizer(wx.HORIZONTAL)
        buttons_sizer.AddStretchSpacer(1)
        buttons_sizer.Add(self._run_button) 
        buttons_sizer.Add(reset_button)        
        # Controls
        controls_sizer = wx.BoxSizer(wx.VERTICAL)
        self._create_gui(controls_sizer)
        # Main layout
        panel_sizer = wx.BoxSizer(wx.VERTICAL)
        panel_sizer.Add(controls_sizer, flag=wx.EXPAND)
        panel_sizer.Add(buttons_sizer, flag=wx.EXPAND)
        self.panel.SetSizer(panel_sizer)
        
        # Events
        self._run_button.Bind(wx.EVT_BUTTON, self.OnRunClicked)
        reset_button.Bind(wx.EVT_BUTTON, self.OnResetClicked)
        self._images.add_observer("any", self._on_images_modified)
        self._viewer_3ds.add_observer("any", self._on_viewer_3ds_modified)
    
    def __call__(self):
        
        # Build the expression namespace
        namespace = {"function" : self._function}
        for parameter in self._parameters :
            if parameter.get("role", "") == "return" :
                namespace[parameter["name"]] = None
            elif parameter.get("role", "") == "output" :
                # Special case for images and 3D viewers
                if parameter["type"] == "Image" and self._controls[parameter["name"]].output_checked :
                    dummy = Image(shape=(1,1,1), value=0, dtype=numpy.single)
                    namespace[parameter["name"]] = dummy
                elif parameter["type"] == "Object3D" :
                    if self._controls[parameter["name"]].output_checked :
                        viewer = Viewer3DFrame(parent=None, objects_3d = ObservableList())
                        viewer.Show()
                        wx.GetApp().append_viewer_3d(viewer)
                        getattr(self, parameter["name"]).value = self._viewer_3ds[-1]
                    dummy = Object3D(name="New object")
                    namespace[parameter["name"]] = dummy
                else :
                    namespace[parameter["name"]] = self._controls[parameter["name"]].value
            else :
                # Parameter will not be modified
                namespace[parameter["name"]] = self._controls[parameter["name"]].value
        
        # Build the expression
        args = []
        return_values = []
        for parameter in self._parameters :
            if parameter.get("role", "") == "return" :
                return_values.append(parameter["name"])
            else :
                if parameter["type"] == "Object3D" :
                    args.append("{0} = objects_3d[\"{0}\"]".format(parameter["name"]))
                else :
                    args.append("{0} = {0}".format(parameter["name"]))
        expression = "function({0})".format(", ".join(args))
        if return_values : 
            return_expression = ", ".join(return_values)
            expression = "{0} = {1}".format(return_expression, expression)
        
        # Execute it
        def f():
            exec expression in namespace
        
        periodic_progress_dialog = PeriodicProgressDialog(0.2, 
            "Running %s ..."%self._function.func_name,
            "Running %s ..."%self._function.func_name)
        worker_thread = WorkerThread(periodic_progress_dialog, target=f)
        worker_thread.start()
        periodic_progress_dialog.start()
        worker_thread.join()
        periodic_progress_dialog.Destroy()
        
        if worker_thread.exception :
            wx.MessageBox("Could not run function : {0}".format(worker_thread.exception),
                          "Could not run function")
        
        # Update controls and application
        for parameter in self._parameters :
            if parameter.get("role", "") not in ["return", "output"] :
                # Parameter has not been modified
                continue 
            else :
                # Special case for Image and Object3D
                if parameter["type"] == "Image" :
                    if self._controls[parameter["name"]].output_checked :
                        wx.GetApp().append_image(namespace[parameter["name"]])
                    else :
                        index = wx.GetApp().images.index(self._controls[parameter["name"]].value)
                        wx.GetApp().close_image(wx.GetApp().gui_images[index])
                        wx.GetApp().insert_image(index, namespace[parameter["name"]])
                elif parameter["type"] == "Object3D" :
                    # TODO
                    # If object has an associated image, set the GUI image
                    object_3d = self._controls[parameter["name"]].value.objects_3d[-1] 
                    image = object_3d.image
                    if image is not None :
                        index = wx.GetApp().images.index(image)
                        object_3d.gui_image = wx.GetApp().gui_images[index]
                    if len(self._controls[parameter["name"]].value.objects_3d) == 1 :
                        getattr(self, parameter["name"]).value.view_all()
                        getattr(self, parameter["name"]).value.update_object_editor()
                # In any case, update the control
                self._controls[parameter["name"]].value = namespace[parameter["name"]]
    
    def validate_form(self):
        all_controls_valid = True
        for parameter in self._parameters :
            control = self._controls[parameter["name"]]
            all_controls_valid = (all_controls_valid and control.validate())
        
        self._run_button.Enable(all_controls_valid)
    
    ##########
    # Events #
    ##########
    def OnRunClicked(self, dummy):
        """ Display a progress dialog and call the function with the correct
            parameters
        """
        
        self()
        
    def OnResetClicked(self,dummy):
        for parameter in self._parameters:
            self._controls[parameter["name"]].reset()
        self.validate_form()
    
    def _on_images_modified(self, dummy):
        self.panel.Fit()
    
    def _on_viewer_3ds_modified(self, dummy):
        self.panel.Fit()
    
    def _on_control_value_changed(self, dummy) :
        self.validate_form()
    
    #####################
    # Private interface #
    #####################
    
    def _create_gui(self, sizer):
        """Create the graphical interface to get back the data 
        """
        
        for parameter in self._parameters :
            initializer = parameter.get("initializer", "")
            # Check if the initializer depends on another parameter
            pattern = re.compile(r"\$\{([a-zA-Z_][a-zA-Z0-9_]*)\}")
            elements = pattern.split(initializer)
            if len(elements) > 1 : 
                initializer = ""
                for i in range(len(elements)) :
                    if i%2 == 0 :
                        initializer += elements[i]
                    else : # i%2 == 1
                        initializer += "self._controls[{0}].value".format(elements[i])
            # Create the control
            expression = "medipy.gui.control.%s(self.panel"%parameter["type"]
            
            if parameter["type"] == "Image" : 
                expression += ", self._images"
            elif parameter["type"] == "Object3D" :
                expression += ", self._viewer_3ds"
                
            if initializer != "" :
                expression += ", " + initializer
            expression += ")"
            
            try :
                control = eval(expression)
                self._controls[parameter["name"]] = control
            except : 
                return
            
            # Add the control to the sizer, and create a corresponding attribute to self
            label = wx.StaticText(self.panel, label=parameter["label"])
            sizer.Add(label)
            sizer.Add(control, flag=wx.EXPAND)
            control.SetToolTipString(parameter.get("tooltip", ""))
            label.SetToolTipString(parameter.get("tooltip", ""))
            setattr(self, parameter["name"], control)
            # Create the event handler in the case of a parameter depending on another
            if len(elements) > 1 :
                observed = set()
                for name in elements[1:len(elements):2] :
                    observed.add(name)
                for name in observed :
                    statement = "self._controls[%s].%s"%(parameter["name"], initializer) 
                    def f(self, event):
                        exec statement in globals(), locals()
                    setattr(FunctionGUIBuilder, parameter["name"] + "_observer", f)
                    getattr(self, name).add_observer("any", getattr(self, parameter["name"] + "_observer"))
            control.add_observer("value", self._on_control_value_changed)
        
        self.validate_form()
