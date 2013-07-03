##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import cPickle
import wx

class Preferences(object) :
    """ Store the preferences of an application. If not provided, the 
        ``application_name`` and ``company_name`` parameters default to the 
        corresponding ``wx.GetApp()`` members.
    
        This class is designed to store anything that can be pickled. Note 
        however that the stored version is slightly larger than the pickled
        version.
        
        Sample usage : ::
        
            app = wx.PySimpleApp()
            app.SetAppName("MediPy")
            app.SetVendorName("Universite de Strasbourg")
            
            # Storing a preference
            preferences = medipy.gui.base.Preferences()
            preferences.set("/io/default_path", "/tmp")
            
            # Retrieving a preference
            preferences = medipy.gui.base.Preferences()
            path = preferences.get("/io/default_path", "/home/user")
            
            # Testing if a preference is stored
            print "/io/default_path" in preferences
    """
    
    def __init__(self, application_name=None, company_name=None) :
        self._config = wx.Config(
            application_name or wx.GetApp().GetAppName(), 
            company_name or wx.GetApp().GetVendorName())
    
    def get(self, key, *args) :
        """ Return a previously-stored preference value. This function can be 
            called either with the name of the preference to retrieve, or with
            the name and a default value. If no default value is provided, and
            the preference is missing, an exception is raised : ::
            
                preferences = medipy.gui.base.Preferences()
                # Return the value of the preference, raise an exception if this
                # preference does not exist
                path = preferences.get("/io/default_path")
                # Return the value of the preference, return "/home/user" if
                # this preference does not exist
                path = preferences.get("/io/default_path", "/home/user")
        """
        if args :
            default_provided = True
            default = args[0]
        else :
            default_provided = False
        
        if key not in self :
            if default_provided :
                return default
            else :
                raise medipy.base.Exception("No such preference: {0}".format(repr(key)))
        else :
            string = self._config.Read(key)
            return cPickle.loads(cPickle.loads(string.encode("utf-8")))
    
    def set(self, key, value) :
        """ Store a preference value.
        """
        
        # wxPython converts all strings to unicode ; however, the result of
        # cPickle.dumps cannot be converted to unicode if its input is already
        # unicode. Using a double pickling solves this problem at the expense of
        # a slightly longer value
        string = cPickle.dumps(cPickle.dumps(value))
        self._config.Write(key, string)
        self._config.Flush()
    
    def delete(self, key) :
        """ Delete an existing preference. If no such preference is stored, 
            raise an exception.
        """
        
        if key not in self :
            raise medipy.base.Exception("No such preference: {0}".format(repr(key)))
        self._config.DeleteEntry(key)
    
    def __contains__(self, key) :
        return self._config.HasEntry(key)
