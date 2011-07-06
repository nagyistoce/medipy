##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import logging
import os
import sys
import traceback

import numpy
import wx

import medipy.io
import medipy.io.dicom
from medipy.gui import PeriodicProgressDialog, WorkerThread
from medipy.gui.dicom import reconstruction
from medipy.gui.io import ImageFileDialog

def load(parent=None, dtype=numpy.single, multiple=False, load_all_images=False,
         wildcard="*"):
    """ Load images with appropriate dialogs (file selector, and DICOM selector
        if needed).
        
        parent : parent wx Window for dialogs
        dtype : type to cast images to, or None to preserve original type
        multiple : allow selection of several images
        load_all_images : load all images contained in a file or just the first one
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
                                     target=medipy.io.dicom.parse, 
                                     args=(path,))
        worker_thread.start()
        periodic_progress_dialog.start()
        worker_thread.join()
        periodic_progress_dialog.Destroy()
            
        if worker_thread.exception is not None :
            wx.MessageBox(
                "Could not load file %s : %s"%(path, worker_thread.exception), 
                "Could not load DICOMDIR")
        else :
            images.extend(
                reconstruction.images([worker_thread.result], parent, dtype, True))
    
    # Load other images
    periodic_progress_dialog = PeriodicProgressDialog(0.2, "Loading image", 
                                                           "Loading ...")
    for path in non_dicomdirs :
        periodic_progress_dialog.Pulse("Loading %s ..."%path)
        periodic_progress_dialog.Fit()
        
        if load_all_images :
            indices = range(medipy.io.number_of_images(path))
        else :
            indices = [0] 
        
        for index in indices :
            worker_thread = WorkerThread(periodic_progress_dialog,
                                        target=medipy.io.load, 
                                        args=(path, index),
                                        kwargs={"type":dtype})
            worker_thread.start()
            periodic_progress_dialog.start()
            
            worker_thread.join()
            
            if worker_thread.exception is not None :
                wx.MessageBox(
                    "Could not load file %s : %s"%(
                        path, worker_thread.exception), 
                    "Could not load image")
            else :
                images.append(worker_thread.result)
            
    periodic_progress_dialog.Destroy()
    
    return images

def import_dicom_directory(parent = None, dtype = numpy.single):
    """ Import an image from a directory containing DICOM files
        
        parent : parent wx Window for dialogs
        dtype : type to cast images to, or None to preserve original type
    """
    
    config = wx.Config("MediPy")
    path = config.Read("ImageFileDialog/DefaultPath")
    
    dialog = wx.DirDialog(parent, style=wx.FD_OPEN)
    dialog.SetPath(path)
    
    if dialog.ShowModal() != wx.ID_OK :
        return []
        
    directory = dialog.GetPath()
    files = os.listdir(directory)
    files = [str(os.path.join(directory, x)) for x in files]
    files.sort()
    
    def loader(files):
        return [medipy.io.dicom.parse(x) for x in files]
    
    periodic_progress_dialog = PeriodicProgressDialog(
        0.2, "Loading files", "Loading ...")
    worker_thread = WorkerThread(periodic_progress_dialog,
                                 target=loader, args=(files,))
    worker_thread.start()
    periodic_progress_dialog.start()
    worker_thread.join()
    periodic_progress_dialog.Destroy()
        
    if worker_thread.exception is not None :
        wx.MessageBox(
            "Could not load files : %s"%(worker_thread.exception, ), 
            "Could not load files")
        return []
    
    images = reconstruction.images(worker_thread.result, parent, dtype, True)
    
    if images :
        config.Write("ImageFileDialog/DefaultPath", dialog.GetPath())
        config.Flush() 
        return images[0]

def save(image, parent=None):
    """ Save an image with appropriate dialogs (file selector)
        
        Return the chosen save path or None
        
        parent : parent wx Window for dialogs
    """
    
    dialog = ImageFileDialog(parent, style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
    
    if dialog.ShowModal() != wx.ID_OK :
        return None
    
    path = dialog.GetPath()
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
        
        dataset = medipy.io.dicom.parse(path)
        if dataset is None :
            is_dicomdir = False
        else :
            is_dicomdir = (dataset.get("media_storage_sop_class_uid", None) == 
                           "1.2.840.10008.1.3.10") # Media Storage Directory Storage
        
        if is_dicomdir :
            dicomdirs.append(path)
        else :
            non_dicomdirs.append(path)
    
    return dicomdirs, non_dicomdirs