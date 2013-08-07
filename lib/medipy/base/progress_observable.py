##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import medipy.base

def progress_observable(f):
    """ Decorator allowing observers to report the progress of a function. The 
        observer must be a function taking a single argument: the progress, 
        a number between 0 (processing has not yet begun) and 1 (processing is 
        done). The observer *should be removed* after use to avoid multiple
        calls. ::
        
            import itk

            import medipy.base
            import medipy.io
            import medipy.itk

            @medipy.base.progress_observable
            def median(image, size):
                itk_image = medipy.itk.medipy_image_to_itk_image(image, False)
                filter = itk.MedianImageFilter[itk_image, itk_image].New(
                    Input=itk_image, Radius=3*(int(size),))
                
                filter.AddObserver(itk.ProgressEvent(), lambda: median.progress(filter.GetProgress()))
                filter()
                
                return medipy.itk.itk_image_to_medipy_image(filter[0], None, True)

            def reporter(progress):
                print progress
            median.add_progress_observer(reporter)

            image = medipy.io.load("image.nii.gz")
            filtered = median(image, 1)
            
            median.remove_progress_observer(reporter)
    """
    
    def progress(value):
        f._observable.notify_observers("progress", progress=float(value))
    
    def add_progress_observer(observer):
        # Define a wrapper to hide the `event` object to the user
        def wrapper(event):
            observer(event.progress)
        f._observable.add_observer("progress", wrapper)
        f._observers[observer] = wrapper
    
    def remove_progress_observer(observer):
        f._observable.remove_observer("progress", f._observers[observer])
        del f._observers[observer]
    
    f._observable = medipy.base.Observable(["progress"])
    f._observers = {}
    f.progress = progress
    f.add_progress_observer = add_progress_observer
    f.remove_progress_observer = remove_progress_observer
    
    return f
