##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import cPickle
import wx

class Preferences(object) :
    """ Store the preferences of an application.
    
        This class is designed to store anything that can be pickled. Note 
        however that the stored version is slightly higher than the pickled
        version.
    """
    
    def __init__(self, application_name, company_name) :
        self._config = wx.Config(application_name, company_name)
    
    def get(self, key, *args) :
        """ Three forms
              * p.get("foo") : throw if "foo" is not in p
              * p.get("foo", "bar") : return "bar" if "foo" is not in p
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
        # wxPython converts all strings to unicode ; however, the result of
        # cPickle.dumps cannot be converted to unicode if its input is already
        # unicode. Using a double pickling solves this problem at the expense of
        # a slightly longer value
        string = cPickle.dumps(cPickle.dumps(value))
        self._config.Write(key, string)
        self._config.Flush()
    
    def delete(self, key) :
        if key not in self :
            raise medipy.base.Exception("No such preference: {0}".format(repr(key)))
        self._config.DeleteEntry(key)
    
    def __contains__(self, key) :
        return self._config.HasEntry(key)
