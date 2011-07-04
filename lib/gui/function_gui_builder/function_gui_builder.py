import logging
import re
import sys
import traceback

import numpy

import wx

import medipy.gui.control
from medipy.base import Image, Object3D, ObservableList
from medipy.gui.viewer_3d_frame import Viewer3DFrame

from builder import parseRSTAsParameters

class FunctionGUIBuilder(object):
    """ Create the GUI of a function based on its docstring.
    
        When clicking the "Run" button, the function will be called with arguments
        extracted from the controls on this GUI.
    """
    
    def __init__(self, parent, function, images, viewer_3ds): 
        ##############
        # Properties #
        ##############
        self._panel = wx.Panel(parent)
        
        #####################
        # Private interface #
        #####################
        self._function = function
        self._images = images
        self._viewer_3ds = viewer_3ds
        self._parameters = None
        self._controls_sizer = None
        self._timer = wx.Timer(self._panel)
        self._progress_dialog = None
        
        ##################
        # Initialization #
        ##################
        self._parameters = parseRSTAsParameters(self._function.__doc__)
        # Widgets
        self._run_button = wx.Button(self._panel, label="Run")
        reset_button = wx.Button(self._panel, label="Reset")
        # Layout for the buttons
        buttons_sizer = wx.BoxSizer(wx.HORIZONTAL)
        buttons_sizer.AddStretchSpacer(1)
        buttons_sizer.Add(self._run_button) 
        buttons_sizer.Add(reset_button)
        # Layout for the controls
        self._controls_sizer = wx.BoxSizer(wx.VERTICAL)
        # Layout for the panel
        panel_sizer = wx.BoxSizer(wx.VERTICAL)
        panel_sizer.Add(self._controls_sizer, flag=wx.EXPAND)
        panel_sizer.Add(buttons_sizer, flag=wx.EXPAND)
        self._panel.SetSizer(panel_sizer)
        
        # Populate the controls sizer
        self._create_gui_elements()
        
        # Events
        self._run_button.Bind(wx.EVT_BUTTON, self.OnRunClicked)
        reset_button.Bind(wx.EVT_BUTTON, self.OnResetClicked)
        self._images.add_observer("any", self._on_images_modified)
        self._viewer_3ds.add_observer("any", self._on_viewer_3ds_modified)
        self._panel.Bind(wx.EVT_TIMER, self.OnProgressTimer)
    
    ####################
    # Public interface #
    ####################
    def __call__(self):
        
        objects_3d = {}
        
        for parameter in self._parameters :
            if parameter.type == "Image" and getattr(self, parameter.name).output_checked :
                dummy = Image(shape=(1,1,1), value=0, dtype=numpy.single)
                wx.GetApp().append_image(dummy)
                getattr(self, parameter.name).value = dummy
                del dummy
            elif parameter.type == "Object3D" :
                if getattr(self, parameter.name).output_checked :
                    viewer = Viewer3DFrame(parent=None, objects_3d = ObservableList())
                    viewer.Show()
                    wx.GetApp().append_viewer_3d(viewer)
                    getattr(self, parameter.name).value = self._viewer_3ds[-1]
                dummy = Object3D(name="New object")
                getattr(self, parameter.name).value.objects_3d.append(dummy)
                objects_3d[parameter.name] = dummy
                del dummy
        
        args = []
        for parameter in self._parameters :
            name = parameter.name
            if parameter.type == "Object3D" :
                args.append(name + " = objects_3d[\"" + name + "\"]")
            else :
                args.append(name + " = self." + name + ".value")
        expression = 'self._function(' + ', '.join(args) + ')'
        try :
            eval(expression)
        except Exception, e :
            exc_info = sys.exc_info()
            logging.error("".join(traceback.format_exception(*exc_info)))
            wx.MessageBox(str(e), "Error")
        else : 
            # Update output images which are not new
            for parameter in self._parameters :
                if parameter.type == "Image" :
                    #if getattr(self, parameter.name).output_checked :
                    #    wx.GetApp().append_image(getattr(self, parameter.name).value)
                    #elif getattr(self, parameter.name).output : 
                    if getattr(self, parameter.name).output :
                        image = getattr(self, parameter.name).value
                        index = wx.GetApp().images.index(image)
                        wx.GetApp().close_image(wx.GetApp().gui_images[index])
                        wx.GetApp().insert_image(index, image)
                elif parameter.type == "Object3D" :
                    # If object has an associated image, set the GUI image
                    object_3d = getattr(self, parameter.name).value.objects_3d[-1] 
                    image = object_3d.image
                    if image is not None :
                        index = wx.GetApp().images.index(image)
                        object_3d.gui_image = wx.GetApp().gui_images[index]
                    if len(getattr(self, parameter.name).value.objects_3d) == 1 :
                        getattr(self, parameter.name).value.view_all()
                        getattr(self, parameter.name).value.update_object_editor()
    
    def validate_form(self):
        all_controls_valid = True
        for parameter in self._parameters :
            control = getattr(self, parameter.name)
            all_controls_valid = (all_controls_valid and control.validate())
        
        self._run_button.Enable(all_controls_valid)
    
    ##############
    # Properties #
    ##############
    panel = property(lambda x:x._panel)
    
    ##########
    # Events #
    ##########
    def OnRunClicked(self, event):
        """ Display a progress dialog and call the function with the correct
            parameters
        """
        
        message = "Executing %s"%self._function.func_name
        self._progress_dialog = wx.ProgressDialog(message, message)
        self._timer.Start(100)
        self._timer.Notify()
        
        self()
        
        self._timer.Stop()
        self._progress_dialog.Update(100)
        self._progress_dialog.Destroy()
        self._progress_dialog = None
        
    def OnResetClicked(self,event):
        for parameter in self._parameters:
            getattr(self, parameter.name).reset()
        self.validate_form()
    
    def OnProgressTimer(self, event):
        if self._progress_dialog is not None and self._progress_dialog.IsShown() :
            self._progress_dialog.Pulse()
    
    def _on_images_modified(self, event):
        self._panel.Fit()
    
    def _on_viewer_3ds_modified(self, event):
        self._panel.Fit()
    
    def _on_control_value_changed(self, event) :
        self.validate_form()

    #####################
    # Private interface #
    #####################
    def _create_gui_elements(self):
        """Create the graphical interface to get back the data 
        """
        
        for parameter in self._parameters :
            initializer = parameter.initializer if parameter.initializer is not None else ""
            # Check if the initializer depends on another parameter
            pattern = re.compile(r"\$\{([a-zA-Z_][a-zA-Z0-9_]*)\}")
            elements = pattern.split(initializer)
            if len(elements) > 1 : 
                initializer = ""
                for i in range(len(elements)) :
                    if i%2 == 0 :
                        initializer += elements[i]
                    else : # i%2 == 1
                        initializer += "self." + elements[i] + ".value"
            # Create the control
            expression = "medipy.gui.control.%s(self._panel"%parameter.type
            
            if parameter.type == "Image" : 
                expression += ", self._images"
            elif parameter.type == "Object3D" :
                expression += ", self._viewer_3ds"
                
            if initializer != "" :
                expression += ", " + initializer
            expression += ")"
            
            try :
                control = eval(expression)
            except : 
                return
            
            # Add the control to the sizer, and create a corresponding attribute to self
            label = wx.StaticText(self._panel, label=parameter.label)
            self._controls_sizer.Add(label)
            self._controls_sizer.Add(control, flag=wx.EXPAND)
            control.SetToolTipString(parameter.tooltip or "")
            label.SetToolTipString(parameter.tooltip or "")
            setattr(self, parameter.name, control)
            # Create the event handler in the case of a parameter depending on another
            if len(elements) > 1 :
                observed = set()
                for name in elements[1:len(elements):2] :
                    observed.add(name)
                for name in observed :
                    statement = "self.%s.%s"%(parameter.name, initializer) 
                    def f(self, event):
                        exec statement in globals(), locals()
                    setattr(FunctionGUIBuilder, parameter.name + "_observer", f)
                    getattr(self, name).add_observer("any", getattr(self, parameter.name + "_observer"))
            control.add_observer("value", self._on_control_value_changed)
        
        self.validate_form()
   
if __name__=='__main__':
    
    import wx
    
    def append_image(self, image):
        images.append(image)
    wx.App.append_image = append_image
    
    def close_image(self, gui_image):
        index = self.gui_images.index(gui_image)
        del self.images[index]
    wx.App.close_image = close_image
    
    def insert_image(self, index, image):
        self.images.insert(index, image)
    wx.App.insert_image = insert_image

    app = wx.App()
    
    from medipy.base import ObservableList
    from medipy.components.io import load
    
    images = ObservableList()
    app.images = images
    app.gui_images = images
    
#    for file in ["avg152T1_LR_nifti.nii", "avg152T1_LR_nifti.nii", "kayassseh_irm_cer_nifti.nii", "zstat1.nii", "lobo_cra_nifti.nii" ] :
#        image = load(os.path.join("/home/lamy/Images/nifti", file))
#        images.append(image)

    frame = wx.Frame(None)
    
    print images
    
    def something(input, value, output) :
        """ Do something
            :gui:
                input : Image
                    Input
                    
                    Input tooltip
                value : Int : range = (${input}.data.min(), ${input}.data.max())
                    Value
                    
                    Value tooltip
                output : Image : output = True
                    Output
                    
                    Output tooltip
        """
    
        print "input=%s, value=%s, output=%s"%(input, value, output)
    
    from medipy.components.visualization_3d import mesh_surface
    from medipy.gui.viewer_3d_frame import Viewer3DFrame
    
    objects = ObservableList()
    
    viewer_3ds = ObservableList()
    #viewer_3ds.append(Viewer3DFrame(None, objects)) 
    
    gui_builder = FunctionGUIBuilder(frame, mesh_surface, images, viewer_3ds)
    sizer = wx.BoxSizer()
    sizer.Add(gui_builder.panel, 1, wx.EXPAND)
    frame.SetSizer(sizer)
    sizer.SetSizeHints(frame)
    frame.Show(True)
    app.MainLoop()
    