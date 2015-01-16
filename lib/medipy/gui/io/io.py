##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import glob
import logging
import os
import sys
import traceback

import numpy
import wx

import medipy.gui.base
import medipy.io
import medipy.io.dicom
from medipy.gui import PeriodicProgressDialog, WorkerThread
from medipy.gui.dicom import reconstruction
from medipy.gui.io import ImageFileDialog

def load(parent=None, dtype=numpy.single, multiple=False, wildcard="*"):
    """ Load images with appropriate dialogs (file selector, and DICOM selector
        if needed).
        
        parent : parent wx Window for dialogs
        dtype : type to cast images to, or None to preserve original type
        multiple : allow selection of several images
        wildcard : dialog wildcard
        
        Return a list of loaded images
    """

    style = (wx.FD_OPEN | wx.FD_MULTIPLE) if multiple else wx.FD_OPEN
    dialog = ImageFileDialog(parent, style=style)
    dialog.SetWildcard(wildcard)
    
    if dialog.ShowModal() != wx.ID_OK :
        return []
    
    paths = dialog.GetPaths() if multiple else [dialog.GetPath()]
    
    images = []
    
    dicomdirs, non_dicomdirs = _separate_dicomdirs(paths)
    
    # Load all dicomdirs
    for path in dicomdirs :
        
        periodic_progress_dialog = PeriodicProgressDialog(
            0.2, "Loading DICOMDIR", "Loading ...")
        worker_thread = WorkerThread(periodic_progress_dialog,
                                     target=medipy.io.dicom.read, 
                                     args=(path,))
        worker_thread.start()
        periodic_progress_dialog.start()
        worker_thread.join()
        periodic_progress_dialog.Destroy()
            
        if worker_thread.exception is not None :
            wx.MessageBox("Could not load file {0} : {1}".format(
                path, worker_thread.exception), "Could not load DICOMDIR")
        else :
            new_images = reconstruction.images(
                [worker_thread.result], parent, dtype)
            if new_images:
                images.append(new_images)
    
    # Load other images
    periodic_progress_dialog = PeriodicProgressDialog(0.2, "Loading image", 
                                                           "Loading ...")
    for path in non_dicomdirs :
        periodic_progress_dialog.Pulse("Loading {0}...".format(path))
        periodic_progress_dialog.Fit()
        
        worker_thread = WorkerThread(periodic_progress_dialog,
                                    target=medipy.io.load, args=(path, dtype))
        worker_thread.start()
        periodic_progress_dialog.start()
        
        worker_thread.join()
        
        if worker_thread.exception is not None :
            wx.MessageBox("Could not load file {0} : {1}".format(
                path, worker_thread.exception), "Could not load image")
        else :
            if worker_thread.result:
                images.append(worker_thread.result)
            
    periodic_progress_dialog.Destroy()
    
    return images

def import_dicom_directory(parent = None, dtype = numpy.single, recursive=False,
                           load_all_images=False):
    """ Import an image from a directory containing DICOM files
        
        parent : parent wx Window for dialogs
        dtype : type to cast images to, or None to preserve original type
        recursive : set to True to recursively search for DICOM files
        load_all_images : load all images contained in a file or just the first one
    """

    preferences = medipy.gui.base.Preferences(
        wx.GetApp().GetAppName(), wx.GetApp().GetVendorName())
    path = preferences.get("IO/load_path", "")
    
    dialog = wx.DirDialog(parent, defaultPath=path, style=wx.FD_OPEN)
    if dialog.ShowModal() != wx.ID_OK :
        return []
       
    files = []
    if recursive : 
        for dirpath, dirnames, filenames in os.walk(dialog.GetPath()) :
            files.extend(os.path.join(dirpath, x) for x in filenames)
    else :
        files = [x for x in glob.glob(os.path.join(dialog.GetPath(), "*"))
                 if os.path.isfile(x)]
    
    def loader(files) :
        datasets = []
        for file in files :
            if medipy.io.dicom.can_read(str(file)) :
                datasets.append(medipy.io.dicom.read(str(file)))
        return datasets
    
    periodic_progress_dialog = PeriodicProgressDialog(
        0.2, "Loading files", "Loading ...")
    worker_thread = WorkerThread(periodic_progress_dialog,
                                 target=loader, args=(files,))
    worker_thread.start()
    periodic_progress_dialog.start()
    worker_thread.join()
    periodic_progress_dialog.Destroy()
    
    images = reconstruction.images(worker_thread.result, parent, dtype)
    
    if images :
        preferences.set("IO/load_path", dialog.GetPath())
        return images
    else :
        return []

def save(image, parent=None):
    """ Save an image with appropriate dialogs (file selector)
        
        Return the chosen save path or None
        
        parent : parent wx Window for dialogs
    """

    dialog = ImageFileDialog(parent, style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
    
    if dialog.ShowModal() != wx.ID_OK :
        return None
    
    path = dialog.GetPath()

    if isinstance(image, list) :
        periodic_progress_dialog = PeriodicProgressDialog(0.2, "Saving files", "Saving ...")
        worker_thread = WorkerThread(periodic_progress_dialog,target=medipy.io.save_serie, args=(image,path,))
        worker_thread.start()
        periodic_progress_dialog.start()
        worker_thread.join()
        periodic_progress_dialog.Destroy()
    else :
        try :
            medipy.io.save(image, path)
            return path
        except :
            # Try by adding a ".nii" suffix
            try :
                path += ".nii"
                medipy.io.save(image, path)
                return path
            except :
                exc_info = sys.exc_info()
                logging.error("".join(traceback.format_exception(*exc_info)))
                wx.MessageBox("Could not save image to %s: %s"%(path, exc_info[1]),
                              "Could not save image", wx.OK|wx.ICON_ERROR)

def _separate_dicomdirs(paths):
    """ Separate dicomdirs from "regular" images
        
        Return a pair (dicomdirs, non_dicomdirs)
    """
    
    dicomdirs = []
    non_dicomdirs = []
    
    for path in paths :
        path = str(path)
        is_dicomdir = False
        if medipy.io.dicom.can_read(str(path)) :
            dataset = medipy.io.dicom.read(str(path))
            is_dicomdir = (dataset.header.get("media_storage_sop_class_uid", None) == 
                           "1.2.840.10008.1.3.10") # Media Storage Directory Storage
        
        if is_dicomdir :
            dicomdirs.append(path)
        else :
            non_dicomdirs.append(path)
    
    return dicomdirs, non_dicomdirs
