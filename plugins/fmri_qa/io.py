# coding: utf-8
##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import csv
import datetime
import inspect
import os

import matplotlib.pyplot
import numpy

from medipy.io import save

def save_signal(image, directory):
    """ Save the signal image to directory/signal.nii.gz
    """
    
    save(image, os.path.join(directory, "signal.nii.gz"))

def save_temporal_fluctuation_noise(image, directory):
    """ Save the temporal fluctuation noise image to directory/tfn.nii.gz
    """
    
    save(image, os.path.join(directory, "tfn.nii.gz"))

def save_sfnr(image, directory):
    """ Save the signal-to-fluctuation-noise ratio image to directory/sfnr.nii.gz
    """
    
    save(image, os.path.join(directory, "sfnr.nii.gz"))
    
def save_static_spatial_noise(image, directory):
    """ Save the static spatial noise image to directory/ssn.nii.gz
    """
    
    save(image, os.path.join(directory, "ssn.nii.gz"))

def save_fluctuation_and_drift(time_series, polynomial, residuals, directory) :
    """ Save fluctuation and drift data to given directory : 
          * time_series.csv : a CSV file containing the average of the signal on
            a given ROI, for each volume.
          * time_series_model : a 2nd order polynomial modeling the time series,
            stored as its coefficients, highest power first.
          * residuals.csv : a CSV file containing the difference between the
            time series and its model.
          * time_series.png : a plot of the fluctuation and its model.
    """
    
    time_series_file = open(os.path.join(directory, "time_series.csv"), "w")
    for x in time_series :
        time_series_file.write("%f\n"%x)
    
    polynomial_file = open(os.path.join(directory, "time_series_model"), "w")
    polynomial_file.write("%f %f %f"%tuple(polynomial))
    
    residuals_file = open(os.path.join(directory, "residuals.csv"), "w")
    for x in residuals :
        residuals_file.write("%f\n"%x)
    
    time_series_figure(time_series, polynomial).savefig(os.path.join(directory, "time_series.png"))

def save_residuals_spectrum(spectrum, directory):
    """ Save the spectrum of the residuals to directory :
          * spectrum.csv : CSV file containing the spectrum data (frequency, magnitude).
          * spectrum_no_dc.png : a plot of the spectrum without the DC component.
    """
    
    spectrum_writer = csv.writer(open(os.path.join(directory, "spectrum.csv"), "w"))
    for x in spectrum.T :
        spectrum_writer.writerow(x)
    spectrum_figure(spectrum).savefig(os.path.join(directory, "spectrum_no_dc.png"))

def save_weisskoff_analysis(fluctuations, theoretical_fluctuations, rdc, directory):
    """ Save the Weisskoff analysis to directory :
          * fluctuations.csv : CSV file containing the fluctuations w.r.t. ROI size.
          * theoretical_fluctuations.csv : CSV file containing the theoretical 
            fluctuations w.r.t. ROI size.
          * weisskoff.png : a plot of the fluctuations and theoretical fluctuations.
    """
    fluctuations_writer = csv.writer(open(os.path.join(directory, "fluctuations.csv"), "w"))
    for x in fluctuations.T :
        fluctuations_writer.writerow(x)
    
    theoretical_fluctuations_writer = csv.writer(
        open(os.path.join(directory, "theoretical_fluctuations.csv"), "w"))
    for x in theoretical_fluctuations.T :
        theoretical_fluctuations_writer.writerow(x)
    
    weisskoff_figure(fluctuations, theoretical_fluctuations, rdc).savefig(
        os.path.join(directory, "weisskoff.png"))

def load_summary(filename):
    """ Return the SFNR, SNR, fluctuation, drift and radius of decorrelation
        from the given file.
    """
    
    fd = open(filename)
    
    result = {}
    for line in fd.readlines() :
        key, value = line.strip().split("=")
        result[key] = float(value)
    
    return result

def save_summary(sfnr, snr, fluctuation, drift, rdc, directory):
    """ Save summary values to directory/summary :
          * Signal-to-fluctuation-noise ratio
          * Signal-to-noise ratio
          * Fluctuation percentage
          * Drift percentage
          * Radius of decorrelation (in pixels)
    """
    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    
    summary_file = open(os.path.join(directory, "summary"), "w")
    
    for name in args :
        if name == "directory" :
            continue
        value = values[name]
        summary_file.write("%s=%s\n"%(name, value))
    
    summary_file.close()

def save_report(date, sfnr, snr, drift, fluctuation, rdc, directory):
    """ Save an HTML report of the QA analysis in directory/report.html
    """
    
    template = """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" dir="ltr" lang="fr-FR">
<head>
    <title>Contrôle qualité du %(date)s</title>
    <meta http-equiv="content-type" content="text/html; charset=UTF-8" />
<body>
    <h1>Contrôle qualité du %(date)s</h1>
    
    <h2>Mesures</h2>
    
    <table border="1">
        <tr><td>Rapport signal-bruit de fluctuation</td><td>%(sfnr).1f</td></tr>
        <tr><td>Rapport signal-bruit</td><td>%(snr).1f</td></tr>
        <tr><td>Fluctuation</td><td>%(fluctuation).2f %%</td></tr>
        <tr><td>Dérive</td><td>%(drift).2f %%</td></tr>
        <tr><td>Rayon de décorrélation</td><td>%(rdc).1f pixels</td></tr>
    </table>
    
    <h2>Graphes</h2>
    
    <img src="%(time_series)s"/><br/>
    <img src="%(spectrum)s"/><br/>
    <img src="%(weisskoff)s"/><br/>
</body>
</html>"""

    data = template%{
        "date" : date.strftime("%d %B %Y"),
        "time_series" : "time_series.png",
        "spectrum" : "spectrum_no_dc.png",
        "weisskoff" : "weisskoff.png",
        "sfnr" : sfnr,
        "snr" : snr, 
        "drift" : drift,
        "fluctuation" : fluctuation,
        "rdc" : rdc
    }
    
    report_file = open(os.path.join(directory, "report.html"), "w")
    report_file.write(data)

def save_longitudinal(snr, sfnr, fluctuation, drift, directory) :
    for name in ["snr", "sfnr", "fluctuation", "drift"] :
        data = locals()[name]
        writer = csv.writer(open(os.path.join(directory, "%s.csv"%name), "w"))
        for date, value in data :
            writer.writerow((date.year, date.month, date.day, date.hour, date.minute, date.second, value))

def save_longitudinal_figures(snr, sfnr, fluctuation, drift, directory) :
    
    long_names = {
        "snr" : "Signal-to-noise ratio",
        "sfnr" : "Signal-to-fluctuation-noise ratio",
        "fluctuation" : "Fluctuation (%)",
        "drift" : "Drift (%)",
    }
    
    for name in ["snr", "sfnr", "fluctuation", "drift"] :
        data = locals()[name]
        figure = matplotlib.pyplot.figure()
        plot = figure.add_subplot(111)
        
        plot.plot([d[0] for d in data], [d[1] for d in data], "k-")
        
        date_range = data[-1][0]-data[0][0]
        
        if date_range.days > 365 :
            major_locator = matplotlib.dates.YearLocator()
            minor_locator = matplotlib.dates.MonthLocator()
            format = "%Y"
        elif date_range.days > 30 :
            major_locator = matplotlib.dates.MonthLocator()
            minor_locator = matplotlib.dates.DayLocator()
            format = "%Y-%m"
        else :
            major_locator = matplotlib.dates.DayLocator()
            minor_locator = matplotlib.dates.HourLocator()
            format = "%Y-%m %H"

        plot.xaxis.set_major_locator(major_locator)
        plot.xaxis.set_major_formatter(matplotlib.dates.DateFormatter(format))
        plot.xaxis.set_minor_locator(minor_locator)
        
        plot.axes.set_xlabel("Date")
        plot.axes.set_ylabel(long_names[name])
        
        figure.savefig(os.path.join(directory, "%s.png"%name))

def spectrum_figure(spectrum):
    """ Return a matplotlib figure containing the Fourier spectrum, without its
        DC coefficient.
    """
    
    figure = matplotlib.pyplot.figure()
    plot = figure.add_subplot(111)
    plot.plot(spectrum[0,1:], spectrum[1,1:], "k-")
    plot.axes.set_xlabel("Frequency (Hz)")
    plot.axes.set_ylabel("Magnitude")
    
    return figure

def weisskoff_figure(fluctuations, theoretical_fluctuations, rdc):
    """ Return a matplotlib figure containing the Weisskoff analysis.
    """
    
    figure = matplotlib.pyplot.figure()
    plot = figure.add_subplot(111)
    
    plot.plot(fluctuations[0], fluctuations[1], "ko-", fillstyle="full")
    plot.plot(theoretical_fluctuations[0], theoretical_fluctuations[1], "ko-", markerfacecolor="w")
    
    plot.axes.loglog()
    plot.axes.set_xlabel("ROI width (pixels)")
    plot.axes.set_ylabel("Fluctuation (%)")
    plot.xaxis.set_major_formatter(matplotlib.pyplot.FormatStrFormatter("%.2f"))
    plot.yaxis.set_major_formatter(matplotlib.pyplot.FormatStrFormatter("%.2f"))
    plot.legend(("Measured", "Theoretical"), "upper right")
    
    return figure

def time_series_figure(time_series, polynomial) :
    """ Return a matplotlib figure containing the time series and its polynomial
        model.
    """
    
    figure = matplotlib.pyplot.figure()
    plot = figure.add_subplot(111)
    
    x = numpy.arange(2, 2+len(time_series))
    model = numpy.polyval(polynomial, x)
    
    plot.plot(x, time_series, "k-")
    plot.plot(x, model, "k-")
    
    plot.axes.set_xlabel("Volume number")
    plot.axes.set_ylabel("Intensity")
    
    return figure