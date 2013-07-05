##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import weakref

import exception

class EventObject(object):
    """ Dummy class for events.
    """
    pass

class Observable(object) :
    """Base class for observable classes.
    This class defines a simple interface to add or remove observers 
    on an object"""
    
    def __init__(self, allowed_events=None) :
        self._allowed_events = []
        self._observers = {}
        
        # A locked object will not notify observers. This will avoid cycles in
        # observable lists
        self._locked = False
        
        if allowed_events is not None :
            for event in allowed_events :
                self.add_allowed_event(event)
        else :
            self._allowed_events = []
        
        if "any" not in self._allowed_events : 
            self.add_allowed_event("any")
    
    def add_allowed_event(self, event):
        """ Add an allowed event.
            Return False if the event was already allowed, else return True.
        """
        
        self._remove_dead_observers()
        
        if event not in self._allowed_events :
            self._allowed_events.append(event)
            self._observers[event] = []
            return True
        else :
            return False
    
    def add_observer(self, event, observer) :
        """Add an observer to the object.
        Return False if the observer was already in the list, 
        otherwise return True"""
        
        self._remove_dead_observers()
        
        if not self.is_allowed_event(event) :
            raise exception.Exception("Event " + event + " is not allowed for type " + str(type(self)))
        
        if "im_self" in dir(observer) and "im_func" in dir(observer) :
            return self._add_bound_method_observer(event, observer)
        else :
            return self._add_function_observer(event, observer)

    def has_observer(self, event, observer):
        """ Test if the object has a specific observer
        """
        
        self._remove_dead_observers()
        
        if "im_self" in dir(observer) and "im_func" in dir(observer) :
            for im_self, im_func in self._observers[event] :
                if im_self() == observer.im_self and im_func == observer.im_func :
                    return True
            return False
        else : 
            return observer in self._observers[event] 
        
        self._observers[event].append((weakref.ref(observer.im_self), observer.im_func))

    def remove_observer(self, event, observer) :
        """Remove an observer from the object.
        Return False if the observer was not in the list. Otherwise return True
        """
        
        self._remove_dead_observers()
        
        if not self.is_allowed_event(event) :
            raise exception.Exception("Event " + event + " is not allowed for type " + str(type(self)))
        
        if "im_self" in dir(observer) and "im_func" in dir(observer) :
            return self._remove_bound_method_observer(event, observer)
        else :
            return self._remove_function_observer(event, observer)
    
    def notify_observers(self, event, *args, **kwargs) :
        """Notify observers of a given event.
        The event must be a member function of the observer. Arguments will
        depend on the event. Events are documented in the concrete observable
        classes""" 
        
        self._remove_dead_observers()
        
        # We are already processing an event
        if self._locked :
            return
        
        self._locked = True
        
        event_info = EventObject()
        setattr(event_info, "object", self)
        setattr(event_info, "event", event)
        for name, value in kwargs.items() :
            setattr(event_info, name, value)
        
        for observer in self._observers.get(event, []) :
            if type(observer) is tuple : 
                self._notify_bound_method_observer(observer, event_info)
            else : 
                self._notify_function_observer(observer, event_info)
        if event != "any" :
            for observer in self._observers["any"] :
                if type(observer) is tuple : 
                    self._notify_bound_method_observer(observer, event_info)
                else : 
                    self._notify_function_observer(observer, event_info)
        
        self._locked = False
    
    def is_allowed_event(self, event):
        """ Test if ``event`` is allowed for the current object.
        """
        
        self._remove_dead_observers()
        return event in self._allowed_events
    
    ##############
    # Properties #
    ##############
    
    def _get_allowed_events(self):
        """ Events allowed for the current object.
        """
        
        return self._allowed_events
    
    allowed_events = property(_get_allowed_events)
    
    #####################
    # Private interface #
    #####################
    
    def _add_bound_method_observer(self, event, observer):
        observers = [(o[0](), o[1]) for o in self._observers[event] if isinstance(o, tuple)]
        if (observer.im_self, observer.im_func) in observers :
            return False
        else :
            self._observers[event].append((weakref.ref(observer.im_self), observer.im_func))
            return True
    
    def _add_function_observer(self, event, observer):
        if observer in self._observers[event] : 
            return False
        else : 
            self._observers[event].append(observer)

    def _remove_bound_method_observer(self, event, observer):
        observers = [(o[0](), o[1]) for o in self._observers[event]]
        if (observer.im_self, observer.im_func) in observers :
            index = observers.index((observer.im_self, observer.im_func))
            del self._observers[event][index]
            return True
        else :
            return False
    
    def _remove_function_observer(self, event, observer):
        if observer in self._observers[event] :
            index = self._observers[event].index(observer)
            del self._observers[event][index]
            return True
        else :
            return False

    def _notify_bound_method_observer(self, observer, event_info):
        im_self, im_func = observer[0], observer[1]
        if im_self() is not None : 
            im_func(im_self(), event_info)
    
    def _notify_function_observer(self, observer, event_info):
        observer(event_info)
    
    def _remove_dead_observers(self):
        for event, observers in self._observers.items() :
            living_observers = []
            for observer in observers :
                if isinstance(observer, tuple) :
                    if observer[0]() is not None :
                        # Weakref is alive, keep observer
                        living_observers.append(observer)
                    # Otherwise discard it
                else :
                    # No weakref involved, keep observer
                    living_observers.append(observer)
            self._observers[event] = living_observers
