##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import wx
import wx.grid

import os

import medipy.network.dicom
import medipy.gui.base
import medipy.base

class PreferencesDialog(medipy.gui.base.Panel):

    _connections = "network/dicom/connections"
    
    class UI(medipy.gui.base.UI):
        def __init__(self):
            self.add = None
            self.delete = None
            self.validate = None
            self.connections = None
                        
            self.controls = ["add","validate","delete","connections"]
            
    def __init__(self, parent=None, *args, **kwargs):
        
        self.list_connections = medipy.base.ObservableList()
        
        #Set TextCtrl and labels associated
        self.headers=['Hostname','Port','Calling AE','Called AE','Description',
            'Retrieve','Retrieve Data']
        self.shortnames={
            "Hostname" : 'host', "Port" : 'port', "Description" : 'description',
            "Calling AE" : 'calling_ae_title', "Called AE" : 'called_ae_title',
            "Retrieve" : "retrieve", "Retrieve Data" : "retrieve_data"}
            
        self.retrieve = ['wado', 'get', 'move']
          
        # User interface
        self.ui = PreferencesDialog.UI()
        
        xrc_file = medipy.base.find_resource(os.path.join("resources","gui","preferences_dialog.xrc"))
        wrappers = []
        medipy.gui.base.Panel.__init__(self, xrc_file, "preferences", 
            wrappers, self.ui, self.ui.controls, parent, *args, **kwargs)
                         
        #Show known connections
        self.ui.connections.SetRowLabelSize(0)
        self.ui.connections.SetDefaultColSize(180)
        self.ui.connections.CreateGrid(0,7)
        self.ui.connections.SetColFormatNumber(1)
        
        attr = wx.grid.GridCellAttr()
        attr.SetRenderer(wx.grid.GridCellEnumRenderer(str(self.retrieve)))
        attr.SetEditor(wx.grid.GridCellChoiceEditor(self.retrieve))
        self.ui.connections.SetColAttr(5, attr)
        
        for column, header in enumerate(self.headers) :
            self.ui.connections.SetColLabelValue(column,header)
 
        #Set Events
        self.ui.add.Bind(wx.EVT_BUTTON,self.OnAdd)
        self.ui.delete.Bind(wx.EVT_BUTTON,self.OnDelete)
        self.ui.validate.Bind(wx.EVT_BUTTON,self.OnValidate)
        self.ui.connections.Bind(wx.grid.EVT_GRID_CELL_CHANGE,
                self._update_connections)
        self.list_connections.add_observer("any", self.save_connections)
        self.list_connections.add_observer("any", self._update_listview)

        preferences = medipy.gui.base.Preferences(
                wx.GetApp().GetAppName(), wx.GetApp().GetVendorName())
        # An observable list cannot be pickled, so a list is stored in the
        # preferences : fill the content of self.list_connections instead
        # of re-assigning it.
        self.list_connections[:] = preferences.get(self._connections, [])
        self._update_listview()
        
        self.Show(True)
    
    def _update_connections(self,event):
        """ Update and save connections observable list from ui entry
        """
        row = event.GetRow()
        col = event.GetCol()
        
        shortname = self.shortnames[self.headers[col]]

        if shortname == "port" :
            try :
                value = int(self.ui.connections.GetCellValue(row,col))
            except :
                value = self.ui.connections.GetCellValue(row,col)
        else : 
            value = self.ui.connections.GetCellValue(row,col)
            
        if shortname == "description":
            self.list_connections[row][0]=value
        elif shortname == "retrieve":
            self.list_connections[row][2]=value
            # On 'get' choice, we don't need any value in retrieve_data
            # Updating listview is necessary to refresh the grid
            if value == 'get' and self.list_connections[row][3]!='':
                self.list_connections[row][3]=''
                self._update_listview()
        elif shortname == "retrieve_data":
            if self.list_connections[row][2]=='get':
                pass
            else:
                self.list_connections[row][3]=value
        else :
            setattr(self.list_connections[row][1],shortname, value)
        self.save_connections()
           
    def OnAdd(self,_):
        """ Add a new connection set with default parameters
        """
        connection = medipy.network.dicom.Connection("----", 0, "----", "----")
        self.list_connections.append(["----", connection, "wado", "----"])

    def OnDelete(self,_):
        """ Delete selected connection from the connections ObservableList
        """
        row = self.ui.connections.GetGridCursorRow()
        del self.list_connections[row]

    def OnValidate(self,_):
        """ Test the selected connection
        """
        row = self.ui.connections.GetGridCursorRow()
        echo = medipy.network.dicom.scu.Echo(self.list_connections[row][1])
        try:
            echo()
        except medipy.base.Exception, e:
            dlg = wx.MessageDialog(self, "Cannot contact entity.\nMake sure you entered the right parameters.",'Connection failure',wx.OK|wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
        else:
            dlg = wx.MessageDialog(self, "Parameters seem to be correct.\nYou may be able to work on this connection",'Connection succeeded',wx.OK|wx.ICON_INFORMATION)
            dlg.ShowModal()
            dlg.Destroy()

    def save_connections(self, *args) :
        """ Save the connections to the preferences. This function can also be
            used as an event handler.
        """
        preferences = medipy.gui.base.Preferences(
                wx.GetApp().GetAppName(), wx.GetApp().GetVendorName())
        preferences.set(self._connections, self.list_connections[:])

    def _update_listview(self, *args, **kwargs) :
        """ Refresh gui on any connections ObservableList modification
        """
        self.ui.connections.ClearGrid()
        delta = self.ui.connections.GetNumberRows()-len(self.list_connections)
        if delta<0 :
            self.ui.connections.AppendRows(abs(delta))
        elif delta>0:
            self.ui.connections.DeleteRows(0,delta)
            
        for row,(description,connection,retrieve,retrieve_data) in enumerate(self.list_connections) :
            for column, header in enumerate(self.headers) :
                shortname = self.shortnames[header]
                if shortname == "description" :
                    value = description
                elif shortname == "retrieve" :
                    value = retrieve
                elif shortname == "retrieve_data":
                    value = retrieve_data
                else :
                    value = getattr(connection, shortname)
                
                self.ui.connections.SetCellValue(row,column, str(value))
