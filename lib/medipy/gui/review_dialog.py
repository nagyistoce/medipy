##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import wx
import medipy.base

class ReviewDialog(wx.Dialog, medipy.base.Observable):
    """ A dialog displaying a list of items to be reviewed. Each item is
        specified by a label (displayed in the dialog) and user data. Each item
        has a status, which can be either ``not_reviewed``, ``in_progress``, or
        ``reviewed``. 
        
        This class fires the "item" event when the current item is changed, and 
        the "status" event when the status of an event is changed. The following
        example illustrates a simple use case. ::
        
            import wx
            import medipy.gui
            
            def item_changed(event) :
                print "{0!r} selected".format(event.item[0])
                print "Load {0!r}".format(event.item[1])
            
            def status_changed(event) :
                print "status of {0!r} changed to {1!r}".format(
                    event.item, medipy.gui.ReviewDialog.Status[event.status])
            
            app = wx.PySimpleApp()
            dialog = medipy.gui.ReviewDialog(None, 
                style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
            
            dialog.add_observer("item", item_changed)
            dialog.add_observer("status", status_changed)
            
            dialog.items = [("Patient 1", "image1.nii.gz"),
                            ("Patient 2", "image2.nii.gz"),
                            ("Patient 3", "image3.nii.gz")]
            
            
            dialog.ShowModal()
    """
    
    #: Review status of an item
    Status = medipy.base.enum("Status", "not_reviewed", "in_progress", "reviewed")
    
    def __init__(self, parent, *args, **kwargs):
        wx.Dialog.__init__(self, parent, *args, **kwargs)
        medipy.base.Observable.__init__(self, ["item", "status"])
        
        # Widgets
        self._items = wx.ListBox(self)
        self._previous = wx.Button(self, label="Previous")
        self._next = wx.Button(self, label="Next") 
        
        # Layout
        buttons_sizer = wx.BoxSizer(wx.HORIZONTAL)
        buttons_sizer.Add(self._previous, 1, wx.EXPAND)
        buttons_sizer.Add(self._next, 1, wx.EXPAND)
        
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(self._items, 1, wx.EXPAND)
        main_sizer.Add(buttons_sizer, flag=wx.EXPAND)
        self.SetSizer(main_sizer)
        
        # Events
        self._items.Bind(wx.EVT_LISTBOX, self._set_item)
        self._previous.Bind(wx.EVT_BUTTON, self.previous)
        self._next.Bind(wx.EVT_BUTTON, self.next)
        
        # Other members
        self._labels_to_index = {}
    
    def get_status(self, index_or_label):
        """ Return the review status of an item specified either by its index
            or by its label.
        """
        
        index = self._get_index(index_or_label)
        return self._items.GetClientData(index)["status"]
    
    def set_status(self, index_or_label, value):
        """ Set the review status of an item specified either by its index or 
            by its label.
        """
        
        index = self._get_index(index_or_label)
        
        data = self._items.GetClientData(index)
        data["status"] = value
        self._items.SetClientData(index, data)
        
        self._items.SetString(index, "{0} ({1})".format(
            data["label"], ReviewDialog.Status[value]))
        
        item = (self._items.GetClientData(index)["label"],
                self._items.GetClientData(index)["data"])
        self.notify_observers("status", item=item, status=data["status"])
    
    def previous(self, *args):
        """ Select the previous item in the list.
        """
        
        if args and isinstance(args[0], wx.Event) :
            self.previous()
        else :
            index = self._items.GetSelection()
            if index < 0 :
                raise medipy.base.Exception("Cannot select previous item: no item selected")
            elif index == 0 :
                raise medipy.base.Exception("Cannot select previous item: first item selected")
            else :
                index = self._items.GetSelection()
                self.item = index-1
    
    def next(self, *args):
        """ Mark the current item as reviewed and go to next item.
        """
        
        if args and isinstance(args[0], wx.Event) :
            self.next()
        else :
            index = self._items.GetSelection()
            if index < 0 :
                raise medipy.base.Exception("Cannot select next item: no item selected")
            else :
                self.set_status(index, ReviewDialog.Status.reviewed)
                if index != self._items.GetCount()-1 :
                    self.item = index+1
                    self.set_status(index+1, ReviewDialog.Status.in_progress)
                else :
                    self.Close()
    
    ##############
    # Properties #
    ##############
    
    def _get_items(self):
        """ List of items (as a pair (label, data)) to be reviewed.
        """
        
        items = []
        for i in range(self._items.GetCount()) :
            items.append((self._items.GetClientData(i)["label"], 
                          self._items.GetClientData(i)["data"]))
        return items
    
    def _set_items(self, items):
        self._labels_to_index.clear()
        self._items.Clear()
        
        for index, (label, data) in enumerate(items) :
            data = {"label" : label, "data" : data, 
                    "status" : ReviewDialog.Status.not_reviewed}
            
            label = "{0} ({1})".format(
                data["label"], ReviewDialog.Status[data["status"]])
            self._items.Append(label, data)
            self._labels_to_index[data["label"]] = index
        
        self.item = 0
        self.set_status(0, ReviewDialog.Status.in_progress)
    
    def _get_item(self):
        """ Item being currently reviewed.
        """
        
        index = self._items.GetSelection()
        item = (self._items.GetClientData(index)["label"],
                self._items.GetClientData(index)["data"])
        return item
    
    def _set_item(self, value):
        if isinstance(value, wx.Event) :
            index = value.GetInt()
            if index != wx.NOT_FOUND :
                self._set_item(index)
        else :
            index = self._get_index(value)
            self._items.SetSelection(index)
            
            if index == self._items.GetCount()-1 :
                self._next.SetLabel("Done")
            else :
                self._next.SetLabel("Next")
            
            if index == 0 :
                self._previous.Disable()
            else :
                self._previous.Enable()
            
            item = (self._items.GetClientData(index)["label"],
                    self._items.GetClientData(index)["data"])
            self.notify_observers("item", item=item)
    
    items = property(_get_items, _set_items)
    item = property(_get_item, _set_item)
    
    #####################
    # Private interface #
    #####################
    
    def _get_index(self, value):
        """ Return a index in the wx.Listbox from either an int (with the same
            same semantics as the indices in Python lists), or a string. In the
            latter case, the string must be a label present in the wx.Listbox.
        """
        
        if isinstance(value, basestring) :
            index = self._labels_to_index[value]
        elif isinstance(value, int) :
            index = value
            if index < 0 :
                index = self._items.GetCount()+index
        else :
            raise medipy.base.Exception("Cannot get index from {0!r}".format(type(value)))
        return index
    
