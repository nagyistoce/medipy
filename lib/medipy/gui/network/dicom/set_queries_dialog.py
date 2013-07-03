##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import wx
import wx.gizmos

import os

import medipy.base
import medipy.gui.base

class SetQueriesDialog(medipy.gui.base.Panel):

    _queries_fields = "network/dicom/queries"
    _queries_list = []

    class UI(medipy.gui.base.UI):
        def __init__(self):
            self.queries = None
            self.save = None
                        
            self.controls = ["queries","save"]
            
    def __init__(self, parent=None, *args, **kwargs):
        
        # User interface
        self.ui = SetQueriesDialog.UI()
        
        xrc_file = medipy.base.find_resource(os.path.join("resources", "gui","set_queries_dialog.xrc"))
        wrappers = []
        medipy.gui.base.Panel.__init__(self, xrc_file, "set_queries", 
            wrappers, self.ui, self.ui.controls, parent, *args, **kwargs)

            
        self.queries = wx.gizmos.EditableListBox(self,label="Queries List")
        sizer = self.ui.queries.GetSizer()
        sizer.Add(self.queries,1,wx.EXPAND)

        self.ui.save.Bind(wx.EVT_BUTTON,self.save_queries)
        self.queries.GetNewButton().Bind(wx.EVT_BUTTON,self.OnNew)
        
        preferences = medipy.gui.base.Preferences(
            wx.GetApp().GetAppName(), wx.GetApp().GetVendorName())
        self._queries_list = preferences.get(self._queries_fields,[])

        for row,query in enumerate(self._queries_list) :        
            self.queries.GetListCtrl().InsertStringItem(row,query)
    
    def OnNew(self,_):
        self.queries.GetListCtrl().InsertStringItem(0,"new query")
          
    def save_queries(self,*args):
        
        listctrl = self.queries.GetListCtrl()
        self._queries_list=[]
        
        for row in range(listctrl.GetItemCount()):
            item = listctrl.GetItem(row)
            if item.GetText()!='':
                self._queries_list.append(item.GetText())

        preferences = medipy.gui.base.Preferences(
            wx.GetApp().GetAppName(), wx.GetApp().GetVendorName())
        preferences.set(self._queries_fields, self._queries_list)
