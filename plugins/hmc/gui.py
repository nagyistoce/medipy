##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import sys

import wx

import medipy.gui.control
import medipy.base
import medipy.io

from HMC import segmentation


class HMCPanel(wx.Panel):

    def __init__(self, parent, *args, **kwargs):
        wx.Panel.__init__(self, parent, *args, **kwargs)


        self._number_input_images = medipy.gui.control.Int(self)
        self._mask = medipy.gui.control.Image(self, wx.GetApp().images)
        self._number_classes = medipy.gui.control.Int(self)
        self._atlas_bool = wx.CheckBox(self, label="Atlas")
        """
        self._atlas_WM = medipy.gui.control.File(self)
        self._atlas_WM_text = wx.StaticText(self, label="Atlas WM")
        self._atlas_DM = medipy.gui.control.File(self)
        self._atlas_DM_text = wx.StaticText(self, label="Atlas DM")
        self._atlas_CSF = medipy.gui.control.File(self)
        self._atlas_CSF_text = wx.StaticText(self, label="Atlas CSF")
        """
        self._number_iter = medipy.gui.control.Int(self)
        self._criterion_outliers = medipy.gui.control.Enum(self, ["percentage","threshold"])
        self._criterion_outliers_value = medipy.gui.control.Float(self)
        self._flair_bool = wx.CheckBox(self)
        self._flair_position = medipy.gui.control.Int(self)
        self._output_seg_image = medipy.gui.control.Image(self, wx.GetApp().images, output=True)
        self._output_outliers_image = medipy.gui.control.Image(self, wx.GetApp().images, output=True)


        self._run = wx.Button(self, label="Run")
        reset = wx.Button(self, label="Reset")

        
        self._bool_atlas = not(self._atlas_bool.GetValue())
        self._bool_flair = not(self._flair_bool.GetValue())


        # Layout
        number_input_images_sizer = wx.BoxSizer(wx.VERTICAL)
        number_input_images_sizer.Add(wx.StaticText(self, label="Number of input images"))
        number_input_images_sizer.Add(self._number_input_images, flag=wx.EXPAND)

        self._input_images_sizer = wx.BoxSizer(wx.VERTICAL)

        mask_sizer = wx.BoxSizer(wx.VERTICAL)
        mask_sizer.Add(wx.StaticText(self, label="Mask Image"))
        mask_sizer.Add(self._mask, flag=wx.EXPAND)

        number_classes_sizer = wx.BoxSizer(wx.VERTICAL)
        number_classes_sizer.Add(wx.StaticText(self, label="Number of classes"))
        number_classes_sizer.Add(self._number_classes, flag=wx.EXPAND)        

        atlas_sizer = wx.BoxSizer(wx.VERTICAL)
        atlas_sizer.Add(self._atlas_bool, flag=wx.EXPAND)

        self._atlas_images_sizer = wx.BoxSizer(wx.VERTICAL)

        """
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
        """
        

        number_iter_sizer = wx.BoxSizer(wx.VERTICAL)
        number_iter_sizer.Add(wx.StaticText(self, label="Number of iterations"))
        number_iter_sizer.Add(self._number_iter, flag=wx.EXPAND)
        number_iter_sizer.Add(wx.StaticText(self, label="Criterion for outliers"))

        outliers_sizer = wx.BoxSizer(wx.HORIZONTAL)
        outliers_sizer.Add(self._criterion_outliers, flag=wx.EXPAND)
        outliers_sizer.Add(self._criterion_outliers_value, flag=wx.EXPAND)

        flair_sizer = wx.BoxSizer(wx.HORIZONTAL)
        flair_sizer.Add(self._flair_bool, flag=wx.EXPAND)
        flair_sizer.Add(self._flair_position, flag=wx.EXPAND)
        self._flair_position.Hide()
        

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
        sizer.Add(mask_sizer, flag=wx.EXPAND)
        sizer.Add(number_classes_sizer, flag=wx.EXPAND)
        sizer.Add(atlas_sizer, flag=wx.EXPAND)
        sizer.Add(self._atlas_images_sizer, flag=wx.EXPAND)
        sizer.Add(number_iter_sizer, flag=wx.EXPAND)
        sizer.Add(outliers_sizer, flag=wx.EXPAND)
        sizer.Add(wx.StaticText(self, label="Outliers on one image only"))
        sizer.Add(flair_sizer)
        sizer.Add(output_images_sizer, flag=wx.EXPAND)
        sizer.Add(buttons_sizer, flag=wx.EXPAND)
        self.SetSizer(sizer) 


        # Events
        self._number_input_images.add_observer("value", self.OnNumberInputImages)
        self._atlas_bool.Bind(wx.EVT_CHECKBOX, self.OnAtlas)
        self._number_classes.add_observer("value", self.OnNumberAtlasImages)
        """
        self._atlas_WM.add_observer("value", self.OnFile)
        self._atlas_DM.add_observer("value", self.OnFile)
        self._atlas_CSF.add_observer("value", self.OnFile)
        """
        self._number_iter.add_observer("value", self._on_control)
        self._criterion_outliers_value.add_observer("value", self._on_control)
        self._flair_bool.Bind(wx.EVT_CHECKBOX, self.OnFlair)
        self._flair_position.add_observer("value", self.OnFlair)
        wx.GetApp().images.add_observer("any", self._on_image)
        self._run.Bind(wx.EVT_BUTTON, self.OnRun)
        reset.Bind(wx.EVT_BUTTON, self.OnReset)


        self._update_gui()


    def run(self):

        self._input_images = []
        for image in self._input_images_object:
            self._input_images.append(image.value)

        self._atlas_images = []
        if self._atlas_bool.GetValue() :
            for atlas in self._atlas_images_object:
                self._atlas_images.append(medipy.io.load(atlas.value))
            """ 
            self._atlas_images = [medipy.io.load(self._atlas_WM.value), medipy.io.load(self._atlas_DM.value), medipy.io.load(self._atlas_CSF.value)]
            """


        self._output_images = hmc(self._input_images, self._mask.value, self._atlas_bool.GetValue(), self._atlas_images, self._number_iter.value, self._number_classes.value, self._criterion_outliers.value, self._criterion_outliers_value.value, self._flair_bool.GetValue(), self._flair_position.value)

        if self._output_seg_image.value is None :
            wx.GetApp().frame.append_image([{"image":self._output_images[0]}])
        else :
            index = wx.GetApp().images.index(self._output_seg_image.value)
            wx.GetApp().frame.delete_image(index)
            wx.GetApp().frame.insert_image(index, [{"image":self._output_images[0]}])
            #wx.GetApp().frame.close_image(wx.GetApp().gui_images[index])
            #wx.GetApp().frame.insert_image(index, self._output_images[0])

        if self._output_outliers_image.value is None :
            wx.GetApp().frame.append_image([{"image":self._output_images[1]}])
        else :
            index = wx.GetApp().images.index(self._output_outliers_image.value)
            wx.GetApp().frame.delete_image(index)
            wx.GetApp().frame.insert_image(index, [{"image":self._output_images[1]}])
            #wx.GetApp().frame.close_image(wx.GetApp().gui_images[index])
            #wx.GetApp().frame.insert_image(index, self._output_images[1])



    def reset(self):
        self._number_input_images.reset()
        self._number_classes.reset()
        self._atlas_bool.SetValue(False)
        """
        self._atlas_WM.reset()
        self._atlas_DM.reset()
        self._atlas_CSF.reset()
        """
        self._update_input_atlas()
        self._number_iter.reset()
        self._criterion_outliers_value.reset()
        self._output_seg_image.reset()
        self._output_outliers_image.reset()
        


    ##################
    # Event handlers #
    ##################

    def OnNumberInputImages(self, event):
        self._update_input_images()

    def OnNumberAtlasImages(self, event):
        self._update_input_atlas()

    def OnAtlas(self, event):
        self._update_atlas()
    
    def OnFlair(self, event):
        self._update_flair()

    def _on_control(self, event):
        self._update_gui()

    def OnRun(self, event):
        self.run()
    
    def OnReset(self, event):
        self.reset()

    def _on_image(self, event):
        self._update_gui()

    def OnFile(self, event):
        self._bool_atlas = True
        for atlas in self._atlas_images_object: 
          self._bool_atlas = self._bool_atlas and atlas.validate()
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
                self._input_images_sizer.Add(wx.StaticText(self, label="Input "+str(i)))
                self._input_images_sizer.Add(image, flag=wx.EXPAND)
                self._input_images_object.append(image)

        self._update_gui()


    def _update_input_atlas(self):
        self._atlas_images_sizer.Clear(True)

        self._atlas_images_object = []
        self._atlas_text_object = []
        if self._number_classes.value is not None:
            for i in range(0, int(self._number_classes.value)):
                atlas = medipy.gui.control.File(self)
                atlas_text = wx.StaticText(self, label="Atlas")
                self._atlas_images_sizer.Add(atlas_text)
                self._atlas_images_sizer.Add(atlas, flag=wx.EXPAND)
                self._atlas_images_object.append(atlas)
                self._atlas_text_object.append(atlas_text)
                atlas.add_observer("value", self.OnFile)
                if self._atlas_bool.GetValue() :
                    atlas.Show()
                    atlas_text.Show()
                else :
                    atlas.Hide()
                    atlas_text.Hide()

        self._update_gui()


    def _update_atlas(self):
        
        """
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
        """
        if self._atlas_bool.GetValue() :
            self._bool_atlas = True
            for atlas in self._atlas_images_object:
                atlas.Show()
                self._bool_atlas = self._bool_atlas and atlas.validate()
            for atlas_text in self._atlas_text_object:
                atlas_text.Show()
        else :
            self._bool_atlas = True
            for atlas in self._atlas_images_object:
                atlas.Hide()
            for atlas_text in self._atlas_text_object:
                atlas_text.Hide()
            

        self._update_gui()

    def _update_flair(self):
        if self._flair_bool.GetValue():
            self._bool_flair = self._flair_position.validate() and self._flair_position.value < self._number_input_images.value
            self._flair_position.Show()
        else:
            self._bool_flair = True
            self._flair_position.Hide()

        self._update_gui()
        
    
    def _update_gui(self):
        self.Fit()

        self._run.Enable(
            self._number_input_images.validate() and self._number_iter.validate() and self._number_classes.validate() and self._criterion_outliers_value.validate() and
                self._number_input_images.value != 0 and self._number_iter.value != 0 and self._number_classes.value != 0 and self._criterion_outliers_value.value != 0 and self._bool_atlas and len(wx.GetApp().images) > 0 and self._bool_flair)




