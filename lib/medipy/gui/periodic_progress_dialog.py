##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import logging
import sys
import threading
import time
import traceback

import wx

class PeriodicProgressDialog(wx.ProgressDialog) : 
    """ Pulsing progress dialog driven by a long-running thread. Typical use : 
    
        >>> import time
        >>> app = wx.PySimpleApp()
        >>> periodic_progress_dialog = PeriodicProgressDialog(0.2, "title", "message")
        >>> worker_thread = WorkerThread(
        ...    periodic_progress_dialog,
        ...    target=time.sleep,
        ...    args=(2,))
        >>> worker_thread.start() # Returns immediately
        >>> periodic_progress_dialog.start() # Returns only when the progress is stopped by the working thread
        >>> worker_thread.join()
        >>> periodic_progress_dialog.Destroy()
        >>> print "all done"
        all done
    """
    
    def __init__(self, period, *args, **kwargs):
        wx.ProgressDialog.__init__(self, *args, **kwargs)
        self._period = period
        self._active = True
    
    def start(self):
        while self._active:
            self.Pulse()
            time.sleep(self._period)

    def stop(self):
        self._active = False

class WorkerThread(threading.Thread) :
    """ Worker thread associated with a periodic progress dialog
    """
    
    def __init__(self, periodic_progress_dialog, group=None, 
                 target=None, name=None, args=(), kwargs={}) :
        self._periodic_progress_dialog = periodic_progress_dialog
        self._exception = None
        self._result = None
        self._target = target
        self._args = args
        self._kwargs = kwargs
        threading.Thread.__init__(self, group, target, name, args, kwargs)
    
    def run(self) :
        try:
            if self._target:
                self._result = self._target(*self._args, **self._kwargs)
        except Exception, e :
            exc_info = sys.exc_info()
            logging.error("".join(traceback.format_exception(*exc_info)))
            self._exception = e
        finally:
            # Avoid a refcycle if the thread is running a function with
            # an argument that has a member that points to the thread.
            del self._target, self._args, self._kwargs

        self._periodic_progress_dialog.stop()
    
    exception = property(lambda x:x._exception)
    result = property(lambda x:x._result)

if __name__ == "__main__":
    import doctest
    doctest.testmod()