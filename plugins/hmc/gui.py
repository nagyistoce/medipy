import sys

import wx

import medipy.gui.control
import medipy.base
import medipy.io

from HMC import hmc


class HMCPanel(wx.Panel):

    def __init__(self, parent, *args, **kwargs):
        wx.Panel.__init__(self, parent, *args, **kwargs)


        self._number_input_images = medipy.gui.control.Int(self)
        self._atlas_bool = wx.CheckBox(self, label="Atlas")
        self._atlas_WM = medipy.gui.control.File(self)
        self._atlas_WM_text = wx.StaticText(self, label="Atlas WM")
        self._atlas_DM = medipy.gui.control.File(self)
        self._atlas_DM_text = wx.StaticText(self, label="Atlas DM")
        self._atlas_CSF = medipy.gui.control.File(self)
        self._atlas_CSF_text = wx.StaticText(self, label="Atlas CSF")
        self._number_iter = medipy.gui.control.Int(self)
        self._number_classes = medipy.gui.control.Int(self)
        self._criterion_outliers = medipy.gui.control.Enum(self, ["percentage","threshold"])
        self._criterion_outliers_value = medipy.gui.control.Float(self)
        self._output_seg_image = medipy.gui.control.Image(self, wx.GetApp().images, output=True)
        self._output_outliers_image = medipy.gui.control.Image(self, wx.GetApp().images, output=True)


        self._run = wx.Button(self, label="Run")
        reset = wx.Button(self, label="Reset")

        self._bool_atlas = True

        self._atlas_images = []


        # Layout
        number_input_images_sizer = wx.BoxSizer(wx.VERTICAL)
        number_input_images_sizer.Add(wx.StaticText(self, label="Number of input images"))
        number_input_images_sizer.Add(self._number_input_images, flag=wx.EXPAND)

        self._input_images_sizer = wx.BoxSizer(wx.VERTICAL)

        atlas_sizer = wx.BoxSizer(wx.VERTICAL)
        atlas_sizer.Add(self._atlas_bool, flag=wx.EXPAND)

        self._atlas_images_sizer = wx.BoxSizer(wx.VERTICAL)
        self._atlas_images_sizer.Add(self._atlas_WM_text)
        self._atlas_images_sizer.Add(self._atlas_WM, flag=wx.EXPAND)
        self._atlas_images_sizer.Add(self._atlas_DM_text)
        self._atlas_images_sizer.Add(self._atlas_DM, flag=wx.EXPAND)
        self._atlas_images_sizer.Add(self._atlas_CSF_text)
        self._atlas_images_sizer.Add(self._atlas_CSF, flag=wx.EXPAND)
        self._atlas_WM.Hide()
        self._atlas_WM_text.Hide()
        self._atlas_DM.Hide()
        self._atlas_DM_text.Hide()
        self._atlas_CSF.Hide()
        self._atlas_CSF_text.Hide()
        

        number_iter_classes_sizer = wx.BoxSizer(wx.VERTICAL)
        number_iter_classes_sizer.Add(wx.StaticText(self, label="Number of iterations"))
        number_iter_classes_sizer.Add(self._number_iter, flag=wx.EXPAND)
        number_iter_classes_sizer.Add(wx.StaticText(self, label="Number of classes"))
        number_iter_classes_sizer.Add(self._number_classes, flag=wx.EXPAND)
        number_iter_classes_sizer.Add(wx.StaticText(self, label="Criterion for outliers"))

        outliers_sizer = wx.BoxSizer(wx.HORIZONTAL)
        outliers_sizer.Add(self._criterion_outliers, flag=wx.EXPAND)
        outliers_sizer.Add(self._criterion_outliers_value, flag=wx.EXPAND)

        output_images_sizer = wx.BoxSizer(wx.VERTICAL)
        output_images_sizer.Add(wx.StaticText(self, label="Output Segmentation Image"))
        output_images_sizer.Add(self._output_seg_image, flag=wx.EXPAND)
        output_images_sizer.Add(wx.StaticText(self, label="Output Outliers Image"))
        output_images_sizer.Add(self._output_outliers_image, flag=wx.EXPAND)

        buttons_sizer = wx.BoxSizer(wx.HORIZONTAL)
        buttons_sizer.AddStretchSpacer(1)
        buttons_sizer.Add(self._run) 
        buttons_sizer.Add(reset)


        sizer = wx.BoxSizer(wx.VERTICAL)    
        sizer.Add(number_input_images_sizer, flag=wx.EXPAND)
        sizer.Add(self._input_images_sizer, flag=wx.EXPAND)
        sizer.Add(atlas_sizer, flag=wx.EXPAND)
        sizer.Add(self._atlas_images_sizer, flag=wx.EXPAND)
        sizer.Add(number_iter_classes_sizer, flag=wx.EXPAND)
        sizer.Add(outliers_sizer, flag=wx.EXPAND)
        sizer.Add(output_images_sizer, flag=wx.EXPAND)
        sizer.Add(buttons_sizer, flag=wx.EXPAND)
        self.SetSizer(sizer) 


        # Events
        self._number_input_images.add_observer("value", self.OnNumberInputImages)
        self._atlas_bool.Bind(wx.EVT_CHECKBOX, self.OnAtlas)
        self._atlas_WM.add_observer("value", self.OnFile)
        self._atlas_DM.add_observer("value", self.OnFile)
        self._atlas_CSF.add_observer("value", self.OnFile)
        self._number_iter.add_observer("value", self._on_control)
        self._number_classes.add_observer("value", self._on_control)
        self._criterion_outliers_value.add_observer("value", self._on_control)
        wx.GetApp().images.add_observer("any", self._on_image)
        self._run.Bind(wx.EVT_BUTTON, self.OnRun)
        reset.Bind(wx.EVT_BUTTON, self.OnReset)


        self._update_gui()


    def run(self):

        self._input_images = []
        for image in self._input_images_object:
            self._input_images.append(image.value)

        if self._atlas_bool.GetValue() : 
            self._atlas_images = [medipy.io.load(self._atlas_WM.value), medipy.io.load(self._atlas_DM.value), medipy.io.load(self._atlas_CSF.value)]


        self._output_images = hmc(self._input_images, self._atlas_bool.GetValue(), self._atlas_images, self._number_iter.value, self._number_classes.value, self._criterion_outliers.value, self._criterion_outliers_value.value)

        if self._output_seg_image.value is None :
            wx.GetApp().append_image(self._output_images[0])
        else :
            index = wx.GetApp().images.index(self._output_seg_image.value)
            wx.GetApp().close_image(wx.GetApp().gui_images[index])
            wx.GetApp().insert_image(index, self._output_images[0])

        if self._output_outliers_image.value is None :
            wx.GetApp().append_image(self._output_images[1])
        else :
            index = wx.GetApp().images.index(self._output_outliers_image.value)
            wx.GetApp().close_image(wx.GetApp().gui_images[index])
            wx.GetApp().insert_image(index, self._output_images[1])



    def reset(self):
        self._number_input_images.reset()
        self._atlas_bool.SetValue(False)
        self._atlas_WM.reset()
        self._atlas_DM.reset()
        self._atlas_CSF.reset()
        self._update_input_atlas()
        self._number_iter.reset()
        self._number_classes.reset()
        self._criterion_outliers_value.reset()
        self._output_seg_image.reset()
        self._output_outliers_image.reset()
        


    ##################
    # Event handlers #
    ##################

    def OnNumberInputImages(self, event):
        self._update_input_images()

    def OnAtlas(self, event):
        self._update_input_atlas()

    def _on_control(self, event):
        self._update_gui()

    def OnRun(self, event):
        self.run()
    
    def OnReset(self, event):
        self.reset()

    def _on_image(self, event):
        self._update_gui()

    def OnFile(self, event):
        self._bool_atlas = self._atlas_WM.validate() and self._atlas_DM.validate() and self._atlas_CSF.validate()
        self._update_gui()

    #####################
    # Private interface #
    #####################

    def _update_input_images(self):
        self._input_images_sizer.Clear(True)

        self._input_images_object = []
        if self._number_input_images.value is not None :
            for i in range(0, int(self._number_input_images.value)) :
                image = medipy.gui.control.Image(self, wx.GetApp().images)
                self._input_images_sizer.Add(wx.StaticText(self, label="Input"))
                self._input_images_sizer.Add(image, flag=wx.EXPAND)
                self._input_images_object.append(image)

        self._update_gui()

    def _update_input_atlas(self):
        
        if self._atlas_bool.GetValue() :
            self._atlas_WM.Show()
            self._atlas_WM_text.Show()
            self._atlas_DM.Show()
            self._atlas_DM_text.Show()
            self._atlas_CSF.Show()
            self._atlas_CSF_text.Show()
            self._bool_atlas = False
        else :
            self._atlas_WM.Hide()
            self._atlas_WM_text.Hide()
            self._atlas_DM.Hide()
            self._atlas_DM_text.Hide()
            self._atlas_CSF.Hide()
            self._atlas_CSF_text.Hide()
            self._bool_atlas = True

        self._update_gui()
        
    
    def _update_gui(self):
        self.Fit() 

        self._run.Enable(
            self._number_input_images.validate() and self._number_iter.validate() and self._number_classes.validate() and self._criterion_outliers_value.validate() and
                self._number_input_images.value != 0 and self._number_iter.value != 0 and self._number_classes.value != 0 and self._criterion_outliers_value.value != 0 and self._bool_atlas and
                    len(wx.GetApp().images) > 0)




