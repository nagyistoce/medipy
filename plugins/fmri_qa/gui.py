import datetime

import wx

import medipy.io.dicom.misc
import medipy.gui
import medipy.gui.control

import io
import cross_sectional

class CrossSectionalPanel(wx.Panel):
    
    def __init__(self, parent, *args, **kwargs):
        wx.Panel.__init__(self, parent, *args, **kwargs)
        
        self._input_directory = medipy.gui.control.Directory(self)
        self._output_directory = medipy.gui.control.Directory(self)
        
        self._run = wx.Button(self, label="Run")
        reset = wx.Button(self, label="Reset")
        
        # Layout
        controls_sizer = wx.BoxSizer(wx.VERTICAL)
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        controls_sizer.Add(wx.StaticText(self, label="Input"))
        controls_sizer.Add(self._input_directory, flag=wx.EXPAND)
        controls_sizer.Add(wx.StaticText(self, label="Output"))
        controls_sizer.Add(self._output_directory, flag=wx.EXPAND)
        
        buttons_sizer = wx.BoxSizer(wx.HORIZONTAL)
        buttons_sizer.AddStretchSpacer(1)
        buttons_sizer.Add(self._run) 
        buttons_sizer.Add(reset)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(controls_sizer, flag=wx.EXPAND)
        sizer.Add(buttons_sizer, flag=wx.EXPAND)
        self.SetSizer(sizer)
        
        # Events
        self._input_directory.add_observer("value", self._on_control)
        self._output_directory.add_observer("value", self._on_control)
        self._run.Bind(wx.EVT_BUTTON, self.OnRun)
        reset.Bind(wx.EVT_BUTTON, self.OnReset)
        
        self._update_gui()
    
    def run(self):
        
        steps = [
            "Loading data...", 
            "1/8 : Computing signal image...",
            "2/8 : Computing temporal fluctuation noise image...",
            "3/8 : Computing signal-to-fluctuation-noise ratio image and summary...",
            "4/8 : Computing static spatial noise...",
            "5/8 : Computing signal to noise ratio...",
            "6/8 : Computing fluctuation and drift...",
            "7/8 : Computing spectrum of residuals...",
            "8/8 : Weisskoff analysis ...",
            "Saving summary..."
        ]
        
        progress_dialog = medipy.gui.PeriodicProgressDialog(0.1, 
            "Running cross-sectional QA", "")
        
        progress_dialog.Pulse("\n".join(steps[:1]))
        progress_dialog.Layout()
        
        worker_thread = medipy.gui.WorkerThread(progress_dialog, 
            target=io.load_data, args=(self._input_directory.value,))
        worker_thread.start()
        progress_dialog.start()
        worker_thread.join()
        if worker_thread.exception :
            pass
        image = worker_thread.result
        
        # Discard the first two volumes (Scanning Protocol, p. 828)
        image.data = image.data[2:,...]
        
        progress_dialog.Pulse("\n".join(steps[:2]))
        progress_dialog.Layout()
        
        def f():
            signal = cross_sectional.get_signal_image(image)
            io.save_signal(signal, self._output_directory.value)
            return signal
        worker_thread = medipy.gui.WorkerThread(progress_dialog, target=f)
        worker_thread.start()
        progress_dialog.start()
        worker_thread.join()
        if worker_thread.exception :
            wx.MessageBox("Could not run cross-sectional QA : {0}".format(worker_thread.exception),
                          "Could not run cross-sectional QA")
            progress_dialog.Destroy()
            return
        signal = worker_thread.result
        
        progress_dialog.Pulse("\n".join(steps[:3]))
        progress_dialog.Layout()
        
        def f():
            tfn = cross_sectional.get_temporal_fluctuation_noise_image(image)
            io.save_temporal_fluctuation_noise(tfn, self._output_directory.value)
            return tfn
        worker_thread = medipy.gui.WorkerThread(progress_dialog, target=f)
        worker_thread.start()
        progress_dialog.start()
        worker_thread.join()
        if worker_thread.exception :
            wx.MessageBox("Could not run cross-sectional QA : {0}".format(worker_thread.exception),
                          "Could not run cross-sectional QA")
            progress_dialog.Destroy()
            return
        tfn = worker_thread.result

        progress_dialog.Pulse("\n".join(steps[:4]))
        progress_dialog.Layout()

        def f():
            sfnr_image, sfnr_summary = cross_sectional.get_sfnr_image(signal, tfn)
            io.save_sfnr(sfnr_image, self._output_directory.value)
            return sfnr_image, sfnr_summary
        worker_thread = medipy.gui.WorkerThread(progress_dialog, target=f)
        worker_thread.start()
        progress_dialog.start()
        worker_thread.join()
        if worker_thread.exception :
            wx.MessageBox("Could not run cross-sectional QA : {0}".format(worker_thread.exception),
                          "Could not run cross-sectional QA")
            progress_dialog.Destroy()
            return
        sfnr_image, sfnr_summary = worker_thread.result
        
        progress_dialog.Pulse("\n".join(steps[:5]))
        progress_dialog.Layout()
        
        def f():
            ssn = cross_sectional.get_static_spatial_noise_image(image)
            io.save_static_spatial_noise(ssn, self._output_directory.value)
            return ssn
        worker_thread = medipy.gui.WorkerThread(progress_dialog, target=f)
        worker_thread.start()
        progress_dialog.start()
        worker_thread.join()
        if worker_thread.exception :
            wx.MessageBox("Could not run cross-sectional QA : {0}".format(worker_thread.exception),
                          "Could not run cross-sectional QA")
            progress_dialog.Destroy()
            return
        ssn = worker_thread.result
        
        progress_dialog.Pulse("\n".join(steps[:6]))
        progress_dialog.Layout()
        
        def f():
            snr = cross_sectional.get_snr(signal, ssn, image.shape[0])
            return snr
        worker_thread = medipy.gui.WorkerThread(progress_dialog, target=f)
        worker_thread.start()
        progress_dialog.start()
        worker_thread.join()
        if worker_thread.exception :
            wx.MessageBox("Could not run cross-sectional QA : {0}".format(worker_thread.exception),
                          "Could not run cross-sectional QA")
            progress_dialog.Destroy()
            return
        snr = worker_thread.result
        
        progress_dialog.Pulse("\n".join(steps[:7]))
        progress_dialog.Layout()
        
        def f():
            time_series, polynomial, residuals, fluctuation, drift = cross_sectional.get_fluctuation_and_drift(image)
            io.save_fluctuation_and_drift(time_series, polynomial, residuals, self._output_directory.value)
            return time_series, polynomial, residuals, fluctuation, drift
        worker_thread = medipy.gui.WorkerThread(progress_dialog, target=f)
        worker_thread.start()
        progress_dialog.start()
        worker_thread.join()
        if worker_thread.exception :
            wx.MessageBox("Could not run cross-sectional QA : {0}".format(worker_thread.exception),
                          "Could not run cross-sectional QA")
            progress_dialog.Destroy()
            return
        time_series, polynomial, residuals, fluctuation, drift = worker_thread.result

        progress_dialog.Pulse("\n".join(steps[:8]))
        progress_dialog.Layout()
        
        def f():
            spectrum = cross_sectional.get_residuals_spectrum(residuals, image.metadata["repetition_time"]/1000.)
            io.save_residuals_spectrum(spectrum, self._output_directory.value)
            return spectrum
        worker_thread = medipy.gui.WorkerThread(progress_dialog, target=f)
        worker_thread.start()
        progress_dialog.start()
        worker_thread.join()
        if worker_thread.exception :
            wx.MessageBox("Could not run cross-sectional QA : {0}".format(worker_thread.exception),
                          "Could not run cross-sectional QA")
            progress_dialog.Destroy()
            return
        spectrum = worker_thread.result

        progress_dialog.Pulse("\n".join(steps[:9]))
        progress_dialog.Layout()
        
        def f():
            fluctuations, theoretical_fluctuations, rdc = cross_sectional.get_weisskoff_analysis(image)
            io.save_weisskoff_analysis(fluctuations, theoretical_fluctuations, rdc, self._output_directory.value)
            return fluctuations, theoretical_fluctuations, rdc
        worker_thread = medipy.gui.WorkerThread(progress_dialog, target=f)
        worker_thread.start()
        progress_dialog.start()
        worker_thread.join()
        if worker_thread.exception :
            wx.MessageBox("Could not run cross-sectional QA : {0}".format(worker_thread.exception),
                          "Could not run cross-sectional QA")
            progress_dialog.Destroy()
            return
        fluctuations, theoretical_fluctuations, rdc = worker_thread.result
        
        progress_dialog.Pulse("\n".join(steps[:10]))
        progress_dialog.Layout()
        
        def f():
            if image.metadata.get("series_date", "") :
                date = medipy.io.dicom.misc.parse_da(image.metadata["study_date"]) 
                if image.metadata.get("series_time", "") :
                    time = medipy.io.dicom.misc.parse_tm(image.metadata["series_time"])
                    date = datetime.datetime.combine(date, time)
            else :
                logging.warning("No Series Date present in image metadata, "
                                "using current date and time")
                date = datetime.datetime.now()
        
            io.save_summary(date, sfnr_summary, snr, fluctuation, drift, rdc, self._output_directory.value)
            io.save_report(date, sfnr_summary, snr, fluctuation, drift, rdc, self._output_directory.value)
        worker_thread = medipy.gui.WorkerThread(progress_dialog, target=f)
        worker_thread.start()
        progress_dialog.start()
        worker_thread.join()
        if worker_thread.exception :
            wx.MessageBox("Could not run cross-sectional QA : {0}".format(worker_thread.exception),
                          "Could not run cross-sectional QA")
            progress_dialog.Destroy()
            return

        progress_dialog.Destroy()
    
    def reset(self):
        self._input_directory.reset()
        self._output_directory.reset()
    
    ##################
    # Event handlers #
    ##################
    
    def _on_control(self, event):
        self._update_gui()
    
    def OnRun(self, event):
        self.run()
    
    def OnReset(self, event):
        self.reset()
    
    #####################
    # Private interface #
    #####################
    
    def _update_gui(self):
        self._run.Enable(
            self._input_directory.validate() and self._output_directory.validate())
