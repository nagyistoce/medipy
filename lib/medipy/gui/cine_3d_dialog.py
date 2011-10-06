import os
import xml.dom.minidom

import wx
import wx.xrc

from medipy.base import find_resource

class Cine3dDialog(wx.Dialog):
    """ Display the dialog to control the cine 3d mode of 3d
    """
    
    def __init__(self, parent, *args, **kwargs):
        
        self._timer = wx.Timer()
        self._delay = None
        self._parent = parent
        self._animated_objects_3d = []
        self._visible_objects_3d = []
        self._stop_flag = True
        
        wx.Dialog.__init__(self, parent, style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER, *args, **kwargs)
        
        dialog_file = find_resource(os.path.join("resources", "gui", "cine_3d_dialog.xrc"))
        self._dialog_xml = xml.dom.minidom.parse(dialog_file)
        
        self._dialog_xrc = wx.xrc.XmlResource("")
        self._dialog_xrc.LoadFromString(str(self._dialog_xml.toxml()))
        
        # Widgets
        self._panel = self._dialog_xrc.LoadPanel(self, "main_panel")
        
        for control in ["objects_check_list_box", "type_combo_box", "cursor_slider", 
                        "cursor_text", "speed_slider", "speed_text", "play_pause_button",
                        "stop_button", "close_button"
                       ]:
            setattr(self, "_"+control, wx.xrc.XRCCTRL(self, control))
        
        self._speed_slider.SetValue(5)
        self._speed_text.SetValue(str(self._speed_slider.GetValue()))
        
        self._play_bitmap = wx.Bitmap(find_resource(os.path.join("resources", "gui", "media-playback-start.png")))
        self._pause_bitmap = wx.Bitmap(find_resource(os.path.join("resources", "gui", "media-playback-pause.png")))
        stop_bitmap = wx.Bitmap(find_resource(os.path.join("resources", "gui", "media-playback-stop.png")))
        
        self._play_pause_button.SetBitmapLabel(self._play_bitmap)
        
        self._stop_button.SetBitmapLabel(stop_bitmap)
        
        # Layout
        global_sizer = wx.BoxSizer(wx.VERTICAL)
        global_sizer.Add(self._panel, 1, wx.EXPAND)

        self.SetTitle("Cine 3D")
        self.SetSizer(global_sizer)
        global_sizer.Layout()
        self.Fit()
        
        # Events
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        self._timer.Bind(wx.EVT_TIMER, self.OnTimer)
        self._type_combo_box.Bind(wx.EVT_COMBOBOX, self.OnComboBox)
        self._cursor_slider.Bind(wx.EVT_SLIDER, self.OnCursorSlider)
        self._cursor_text.Bind(wx.EVT_TEXT, self.OnCursorText)
        self._speed_slider.Bind(wx.EVT_SLIDER, self.OnSpeedSlider)
        self._speed_text.Bind(wx.EVT_TEXT, self.OnSpeedText)
        self._play_pause_button.Bind(wx.EVT_BUTTON, self.OnPlayPauseButtonClick)
        self._stop_button.Bind(wx.EVT_BUTTON, self.OnStopButtonClick)
        self._close_button.Bind(wx.EVT_BUTTON, lambda event : self.Close())
        
        self._parent._objects_check_list_box.Bind(wx.EVT_CHECKLISTBOX, self._on_cine_objects_3d_modified)
        self._parent.objects_3d.add_observer("any", self._on_cine_objects_3d_modified)
        
        self._objects_check_list_box.Bind(wx.EVT_CHECKLISTBOX, self._on_fixed_objects_3d_modified)
        
        self._set_animated_objects()
        self._set_cursor_slider()
        self._update_objects_check_list_box()
        self._reset_animation()
        
        self._set_fps(self._speed_slider.GetValue())
        
        self.stop()
    
    ####################
    # Public interface #
    ####################
    def play(self) :
        self._set_fps(self._speed_slider.GetValue())
        self._timer.Start(self._delay)
        self._stop_flag = False
        
        self._objects_animation()
        
        self._play_pause_button.SetBitmapLabel(self._pause_bitmap)
    
    def pause(self):
        self._timer.Stop()
        
        self._play_pause_button.SetBitmapLabel(self._play_bitmap)
    
    def stop(self):
        self._timer.Stop()
        self._stop_flag = True
        self._set_cursor(0)
        
        self._clear_all_objects()
        self._show_visible_objects()
        self._show_fixed_objects()
        
        self._parent._viewer_3d.render_window_interactor.Render()
        
        self._play_pause_button.SetBitmapLabel(self._play_bitmap)
    
    #####################
    # Private interface #
    #####################
    
    def _set_fps(self, fps):
        self._delay = 1000. / float(fps)
        if self._timer.IsRunning() :
            self._timer.Start(self._delay)
    
    
    def _set_animated_objects(self):
        self._visible_objects_3d = []
        for i in range(len(self._parent.objects_3d)) :
            if self._parent._objects_check_list_box.IsChecked(i) :
                self._visible_objects_3d.append(self._parent.objects_3d[i])
        
        self._animated_objects_3d = []
        for i in range(len(self._visible_objects_3d)) :
            if not (self._objects_check_list_box.IsChecked(i)) : # _object fixe
                self._animated_objects_3d.append(i)
        
        if self._type_combo_box.GetValue() == "Step by step" :
            if self._animated_objects_3d[-1] is None :
                self._animated_objects_3d.pop()
        elif self._type_combo_box.GetValue() == "Incremental" :
            if self._animated_objects_3d[-1] is not None :
                self._animated_objects_3d.append(None)
    
    
    def _set_cursor_slider(self):
        if len(self._animated_objects_3d) > 1 :
            self._cursor_slider.SetMax(len(self._animated_objects_3d)-1)
            self._cursor_slider.Enable()
        else :
            self._cursor_slider.Disable() 
    
    def _set_cursor(self, cursor):
        self._cursor = cursor
        self._cursor_slider.SetValue(self._cursor)
        self._cursor_text.ChangeValue(str(self._cursor+1))
        self._cursor_text.SetInsertionPointEnd()
        
    def _clear_all_objects(self):
        for i in range(len(self._parent.objects_3d)) :
            self._parent.objects_3d[i].visibility = False
    
    def _show_visible_objects(self) :
        for _object in self._visible_objects_3d :
            _object.visibility = True
            
    def _show_fixed_objects(self):
        for i in range(len(self._visible_objects_3d)) :
            if self._objects_check_list_box.IsChecked(i) :
                self._visible_objects_3d[i].visibility = True
            
    def _reset_animation(self):
        val = self._cursor_slider.GetValue()
        self._set_cursor(val)
        self._clear_all_objects()
        self._show_fixed_objects()
        
        if self._type_combo_box.GetValue() == "Step by step" :
            self._visible_objects_3d[self._animated_objects_3d[self._cursor]].visibility = True
        elif self._type_combo_box.GetValue() == "Incremental" :
            self._incremental_animation()
                
        self._parent._viewer_3d.render_window_interactor.Render()      
    
    
    def _objects_animation(self):
        if self._cursor == 0 :
            self._clear_all_objects()
            self._show_fixed_objects()
        
        if self._type_combo_box.GetValue() == "Step by step" :
            self._visible_objects_3d[self._animated_objects_3d[self._cursor]].visibility = True
            if self._cursor > 0 :
                self._visible_objects_3d[self._animated_objects_3d[self._cursor-1]].visibility = False
        elif self._type_combo_box.GetValue() == "Incremental" :
            self._incremental_animation()
        
        self._parent._viewer_3d.render_window_interactor.Render()
        
        
    def _incremental_animation(self):
        index = self._animated_objects_3d[self._cursor]
        if index is None :
            self._clear_all_objects()
            self._show_fixed_objects()
        else :
            self._visible_objects_3d[self._animated_objects_3d[self._cursor]].visibility = True
            for i in range(self._cursor) :
                self._visible_objects_3d[self._animated_objects_3d[i]].visibility = True
        
        
    def _on_cine_objects_3d_modified(self, event):
        self._set_animated_objects()
        self._set_cursor_slider()
        self._update_objects_check_list_box()
        self._reset_animation()
        if self._stop_flag :
            self.stop()
            
    
    def _update_objects_check_list_box(self): 
        self._objects_check_list_box.Clear()
            
        for _object in self._visible_objects_3d :
            self._objects_check_list_box.Append(_object.name)
            #if _object in self._fixed_objects :
            #    self._objects_check_list_box.Checked()
            
            
    def _on_fixed_objects_3d_modified(self, event):
        self._set_animated_objects()
        self._set_cursor_slider()
        
        
            
    ##################
    # Event handlers #
    ##################
    
    def OnPlayPauseButtonClick(self, event):
        if self._play_pause_button.GetBitmapLabel() == self._play_bitmap :
            self.play()
        else :
            self.pause()
    
    def OnStopButtonClick(self, event):
        self.stop()
    
    def OnSpeedSlider(self, event):
        self._speed_text.SetValue(str(self._speed_slider.GetValue()))
        self._set_fps(self._speed_slider.GetValue())
    
    def OnSpeedText(self, event):
        try :
            value = int(self._speed_text.GetValue())
        except :
            self._speed_text.SetBackgroundColour(wx.RED)
        else :
            if self._speed_slider.GetMin() <= value <= self._speed_slider.GetMax() : 
                self._speed_slider.SetValue(value)
                self._speed_text.SetBackgroundColour(None)
                self._set_fps(self._speed_slider.GetValue())
            else :  
                self._speed_text.SetBackgroundColour(wx.RED)
    
    def OnClose(self, event):
        self._parent._objects_check_list_box.Unbind(wx.EVT_CHECKLISTBOX)
        self._parent.objects_3d.remove_observer("any", self._on_cine_objects_3d_modified)
        self._parent._objects_check_list_box.Bind(wx.EVT_CHECKLISTBOX, self._parent.OnObjectChecked)
        
        self.stop()
        event.Skip()
    
    def OnTimer(self, event):
        cursor = (self._cursor+1)%len(self._animated_objects_3d)
        self._set_cursor(cursor)
        self._objects_animation()
    
    def OnCursorSlider(self, event):
        self._stop_flag = False
        self._reset_animation()
        
    def OnComboBox(self, event):
        self._set_animated_objects()
        self._set_cursor_slider()
        if (not self._stop_flag) :
            self._reset_animation()
    
    def OnCursorText(self, event):
        try :
            value = int(self._cursor_text.GetValue())-1
        except :
            self._cursor_text.SetBackgroundColour(wx.RED)
        else :
            if self._cursor_slider.GetMin() <= value <= self._cursor_slider.GetMax() : 
                self._cursor_text.SetBackgroundColour(None)
                self._set_cursor(value)
                self._stop_flag = False
                self._reset_animation()
            else :  
                self._cursor_text.SetBackgroundColour(wx.RED)
