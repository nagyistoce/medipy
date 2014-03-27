##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import numpy
import wx

from medipy.gui import PeriodicProgressDialog, WorkerThread
from explorer_dialog import ExplorerDialog
from stacks_dialog import StacksDialog

import medipy.io.dicom

def images(datasets, parent, dtype=numpy.single, 
           size=(700,700), style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER) :
    """ Return a list of medipy.base.Image objects from a list of DataSets.
    """
    
    dialog = ExplorerDialog(parent, size=size,style=style)
    dialog.set_datasets(datasets)
    
    # Size must be specified explicitely on Windows
    dialog.SetSize((700,700))
    
    if dialog.ShowModal() != wx.ID_OK :
        dialog.Destroy()
        return []
    dialog.Destroy()
    
    # Get selected series from dialog, load it.
    series = dialog.get_selected_datasets()
    
    periodic_progress_dialog = PeriodicProgressDialog(
            0.2, "Loading files", "Loading files ...")
    worker_thread = WorkerThread(
        periodic_progress_dialog, 
        target=medipy.io.dicom.load_dicomdir_records, args=(series,))
    worker_thread.start()
    periodic_progress_dialog.start()
    worker_thread.join()
    periodic_progress_dialog.Destroy()
    if worker_thread.exception is not None :
        wx.MessageBox(
            "Could not load series : %s"%(worker_thread.exception,), 
            "Could not load series")
        return []
    
    series = worker_thread.result
     
    # Reconstruct one image per stack
    series = medipy.io.dicom.split.images(series)
    series = medipy.io.dicom.normalize.normalize(series)
    stacks = medipy.io.dicom.stacks(series)
    
    result = []
    periodic_progress_dialog = PeriodicProgressDialog(
            0.2, "Reconstructing images", 
            "Reconstructing images (%i/%i) ..."%(0, len(stacks)))
    for index, stack in enumerate(stacks) :
        periodic_progress_dialog.Pulse(
            "Reconstructing images (%i/%i) ..."%(index+1, len(stacks)))
        
        worker_thread = WorkerThread(periodic_progress_dialog,
            target=medipy.io.dicom.image, args=(stack,))
        worker_thread.start()
        periodic_progress_dialog.start()
        worker_thread.join()
        if worker_thread.exception is not None :
            wx.MessageBox(
            "Could not reconstruct image : %s"%(worker_thread.exception,), 
            "Could not reconstruct image")
        else :
            image = worker_thread.result
            if dtype is not None :
                image.data = image.data.astype(dtype)
            result.append(image)
    periodic_progress_dialog.Destroy()
        
    return result
