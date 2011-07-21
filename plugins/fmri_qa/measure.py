##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import os
import logging
import sys

import numpy

from medipy.base import Image
import medipy.io.dicom
import medipy.io.dicom.misc
import medipy.io.dicom.split

from medipy.fmri_qa import io

def measure(image, output_directory):
    """ Measure the quality of a set of images.
    """
    
    # Discard the first two volumes (Scanning Protocol, p. 828)
    logging.info("Discarding the first two volumes ... ")
    image.data = image.data[2:,...]
    
    logging.info("1/8 : Computing signal image ... ")
    signal = get_signal_image(image)
    io.save_signal(signal, output_directory)
    
    logging.info("2/8 : Computing temporal fluctuation noise image ... ")
    tfn = get_temporal_fluctuation_noise_image(image)
    io.save_temporal_fluctuation_noise(tfn, output_directory)
    
    logging.info("3/8 : Computing signal-to-fluctuation-noise ratio image and summary ... ")
    sfnr_image, sfnr_summary = get_sfnr_image(signal, tfn)
    io.save_sfnr(sfnr_image, output_directory)
    
    logging.info("4/8 : Computing static spatial noise ... ")
    ssn = get_static_spatial_noise_image(image)
    io.save_static_spatial_noise(ssn, output_directory)
    
    logging.info("5/8 : Computing signal to noise ratio ... ")
    snr = get_snr(signal, ssn, image.shape[0])
    
    logging.info("6/8 : Computing fluctuation and drift ... ")
    time_series, polynomial, residuals, fluctuation, drift = get_fluctuation_and_drift(image)
    io.save_fluctuation_and_drift(time_series, polynomial, residuals, output_directory)
    
    logging.info("7/8 : Computing spectrum of residuals ... ")
    spectrum = get_residuals_spectrum(residuals, image.metadata["repetition_time"]/1000.)
    io.save_residuals_spectrum(spectrum, output_directory)
    
    logging.info("8/8 : Weisskoff analysis ... ")
    fluctuations, theoretical_fluctuations, rdc = get_weisskoff_analysis(image)
    io.save_weisskoff_analysis(fluctuations, theoretical_fluctuations, rdc, output_directory)
    
    io.save_summary(sfnr_summary, snr, fluctuation, drift, rdc, output_directory)
    
    date = medipy.io.dicom.misc.parse_da(image.metadata["series_date"])
    io.save_report(date, sfnr_summary, snr, fluctuation, drift, rdc, output_directory)

def get_signal_image(image) :
    """ The signal image is the simple average, voxel by voxel, across [all] 
        volumes (p. 828).
    """
    
    return Image(data=numpy.average(image, 0).astype(numpy.single), 
                 origin=image.origin[1:], spacing=image.spacing[1:],
                 direction=image.direction[1:,1:])

def get_temporal_fluctuation_noise_image(image, verbose=True) :
    """ To calculate the fluctuation noise image, the time-series across [all]
        images for each voxel is detrended with a second-order polynomial.
        The fluctuation noise image is an image of the standard deviation (SD) 
        of the residuals, voxel by voxel, after this detrending step.
    """
    
    flat = image.data.ravel()
    
    voxels_per_volume = reduce(lambda x,y:x*y, image.shape[1:], 1)
    tfn = numpy.ndarray((voxels_per_volume,))
    
    x = numpy.arange(image.shape[0])
    
    for i in range(voxels_per_volume) :
        y = flat[i::voxels_per_volume]
        polynomial = numpy.polyfit(x, y, 2)
        model = numpy.polyval(polynomial, x)
        residuals = y-model
        tfn[i] = numpy.std(residuals)
    
    return Image(data=tfn.reshape(image.shape[1:]).astype(numpy.single), 
                 origin=image.origin[1:], spacing=image.spacing[1:],
                 direction=image.direction[1:,1:])

def get_sfnr_image(signal_image, tfn_image, roi_radius=10) :
    """ The signal image and the temporal fluctuation image are divided voxel
        by voxel to create the SFNR image. A 21 x 21 voxel ROI, placed in the
        center of the image, is created. The average SFNR across these 441 
        voxels is the SFNR summary value. (p. 828)
    """
    
    sfnr_image = Image(
        data = numpy.divide(signal_image, tfn_image+numpy.finfo(float).eps).astype(numpy.single),
        origin=signal_image.origin, spacing=signal_image.spacing,
        direction=signal_image.direction)
    
    center = numpy.divide(sfnr_image.shape, 2.).round()
    roi = sfnr_image[center[0],
                     center[1]-roi_radius:center[1]+roi_radius+1,
                     center[2]-roi_radius:center[2]+roi_radius+1]
    
    summary = numpy.average(roi)
    
    return (sfnr_image, summary)

def get_static_spatial_noise_image(image) :
    """ The first step is to sum all of the odd-numbered images (sumODD image)
        and separately sum all of the even-numbered images (sumEVEN image). The
        difference between the sum of the odd images and the sum of the even
        images (DIFF = sumODD - sumEVEN) is taken as a raw measure of static 
        spatial noise. (p. 828-829)
    """
    
    image_odd = image[range(1, image.shape[0],2)].astype(numpy.single)
    sum_odd = numpy.sum(image_odd, 0)
    
    image_even = image[range(0, image.shape[0],2)].astype(numpy.single)
    sum_even = numpy.sum(image_even, 0)
    
    diff = sum_odd-sum_even
    
    return Image(data=diff,
                 origin=image.origin[1:], spacing=image.spacing[1:],
                 direction=image.direction[1:,1:])

def get_snr(signal_image, ssn_image, nb_time_points, roi_radius=10) :
    """ SNR = (signal summary value)/sqrt((variance summary value)/[all] time points).
        p. 829
    """

    c = numpy.divide(signal_image.shape, 2.).round()
    
    roi = (c[0], 
           slice(c[1]-roi_radius, c[1]+roi_radius+1), 
           slice(c[2]-roi_radius, c[2]+roi_radius+1))
    
    signal_summary = numpy.average(signal_image[roi])
    variance_summary = numpy.var(ssn_image[roi])
    
    return signal_summary/numpy.sqrt(variance_summary/nb_time_points)

def get_fluctuation_and_drift(image, roi_radius=10) :
    """ A time-series of the average intensity within a 21 x 21 voxel ROI
        centered in the image is obtained. A second-order polynomial trend is
        fit to these data. The mean signal intensity of the time-series (prior
        to detrending) and SD of the residuals after subtracting the fit line
        from the data, are computed. Percent fluctuation equals 
        100*(SD of the residuals)/(mean signal intensity). Drift is computed by
        subtracting the minimum fit value from the maximum fit and dividing by
        the mean signal intensity. This drift value is also multiplied by 100
        to obtain a percentage.
    """

    c = numpy.divide(image.shape[1:], 2.).round()
    
    roi = image[:,c[0], 
                slice(c[1]-roi_radius, c[1]+roi_radius+1), 
                slice(c[2]-roi_radius, c[2]+roi_radius+1)].astype(numpy.single)
    
    intensity = numpy.average(roi)
    
    time_series = numpy.sum(numpy.sum(roi, 1), 1)
    time_series /= roi.shape[-2]*roi.shape[-1]
    
    x = numpy.arange(image.shape[0])
    polynomial = numpy.polyfit(x, time_series, 2)
    time_series_model = numpy.polyval(polynomial, x)
    
    residuals = time_series-time_series_model
    
    fluctuation = 100.*numpy.std(residuals)/intensity
    drift = 100.*(time_series_model.max()-time_series_model.min())/intensity
    
    return time_series, polynomial, residuals, fluctuation, drift

def get_residuals_spectrum(residuals, repetition_time) :
    """ Fourier analysis of mean signal intensity in the ROI over time (volume
        number). After the data are detrended with a second-order polynomial,
        the residuals are submitted to a fast Fourier transform (FFT). 
        p. 829.
    """

    fft = numpy.fft.rfft(residuals)
    if residuals.size%2 == 0 :
        # Discard real term for frequency n/2
        fft = fft[:-1]
    fftfreq = numpy.fft.fftfreq(len(residuals), repetition_time)
    return numpy.vstack((fftfreq[:len(fft)], numpy.abs(fft)))

def get_weisskoff_analysis(image, max_size=21):
    """ The Weisskoff analysis (8) provides another measure of scanner 
        stability included in the GSQAP.
        p. 829.
    """
    
    c = numpy.divide(image.shape[1:], 2.).round()
    
    x = numpy.arange(1, max_size+1)
    fluctuation = [get_fluctuation_and_drift(image, r)[-2] for r in x]
    theoretical_fluctuation = fluctuation[0]/x
    rdc = fluctuation[0]/fluctuation[-1]
    
    return numpy.vstack((x, fluctuation)), numpy.vstack((x, theoretical_fluctuation)), rdc