##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import weakref

from observable import Observable

class PropertySynchronized(Observable) :
    """ A PropertySynchronized object will replicate any event to a child.
        This class expects the following formalism :
          * each value to be synchronized is a property
          * each property generate an event bearing its name
        
        >>> class Image(PropertySynchronized) :
        ...     def __init__(self, position=None, zoom=None) :
        ...         self._position = position
        ...         self._zoom = zoom
        ...         
        ...         PropertySynchronized.__init__(self, ["position", "zoom"])
        ...     
        ...     ##############
        ...     # Properties #
        ...     ##############
        ...     
        ...     def _set_position(self, position) :
        ...         self._position = position
        ...         self.notify_observers("position")
        ...     
        ...     def _set_zoom(self, zoom) :
        ...         self._zoom = zoom
        ...         self.notify_observers("zoom")
        ...     
        ...     position = property(lambda x:x._position, _set_position)
        ...     zoom = property(lambda x:x._zoom, _set_zoom)
        >>> image1 = Image("somewhere", 2.0)
        >>> image2 = Image("nowhere", 0.5)
        
        Upon creation, the images are not synchronized
        >>> image1.position = "north_pole"
        >>> image2.position
        'nowhere'
        
        Define each image to be the child of the other : image1 will propagate
        events to image2, and vice versa. No infinite loop will occur
        >>> image1.append_child(image2)
        >>> image2.append_child(image1)
        
        The position and zoom properties are now synchronized.
        >>> image1.position = "here"
        >>> image2.position
        'here'
        >>> image2.zoom = 1.0
        >>> image1.zoom
        1.0
        
        An event can be de-synchronized and the re-synchronized
        >>> image1.synchronize_on_event["position"] = False
        >>> image1.position = "away"
        >>> image2.position
        'here'
        >>> image1.synchronize_on_event["position"] = True
        >>> image1.position = "back"
        >>> image2.position
        'back'
        
        It is possible to deactivate all synchronization events and 
        reactivate them all :
        >>> image1.synchronize_on_no_event()
        >>> image1.zoom = 0.5
        >>> image2.zoom
        1.0
        >>> image1.synchronize_on_all_events()
        >>> image1.zoom = 1.0
        >>> image2.zoom
        1.0
    """
    
    def __init__(self, *args, **kwargs) :
        self._children = []
        self.synchronize_on_event = {}
        
        Observable.__init__(self, *args, **kwargs)
        
        self.add_observer("any", self._replicate_event)
        self.synchronize_on_all_events()
    
    def synchronize_on_all_events(self) :
        """ Enable the synchronization for all events 
        """
        for event in self.allowed_events :
            self.synchronize_on_event[event] = True
    
    def synchronize_on_no_event(self) :
        """ Disable the synchronization for all events
        """
        for event in self.allowed_events :
            self.synchronize_on_event[event] = False
    
    def append_child(self, child):
        self.insert_child(len(self._children), child)
    
    def insert_child(self, index, child):
        self._children.insert(index, weakref.ref(child))
    
    def delete_child(self, index) :
        del self._children[index]
    
    def delete_all_children(self):
        while self._children :
            self.delete_child(0)
    
    def has_child(self, object) :
        have_it = False
        for child in self._children :
            if child() is object :
                have_it = True
                break
        return have_it
    
    def child_index(self, object):
        have_it = False
        for index, child in enumerate(self._children) :
            if child() is object :
                have_it = True
                break
        if not have_it :
            raise ValueError("PropertySynchronized.child_index(x) : x not in children")
        else :
            return index
    
    #####################
    # Private interface #
    #####################
    
    def _replicate_event(self, event) :
        if self.synchronize_on_event.get(event.event, False) :
            value = getattr(self, event.event)
            for child in self._children :
                if child() is not None and not child()._locked :
                    setattr(child(), event.event, value)
