##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import os
import socket

import wx
import wx.grid

import medipy.network.dicom
import medipy.gui.base
import medipy.base

import logging
try :
    import paramiko
except ImportError :
    logging.warning("medipy.io.dicom.SSHTunnelConnect not usable: cannot find paramiko package")

class PreferencesDialog(medipy.gui.base.Panel):

    _connections = "network/dicom/connections"
    
    class UI(medipy.gui.base.UI):
        def __init__(self):
            self.add = None
            self.close = None
            self.connections = None
                        
            self.controls = ["add","close","connections"]
            
    def __init__(self, parent=None, *args, **kwargs):
        
        self.list_connections = medipy.base.ObservableList()
        
        #Set TextCtrl and labels associated
        self.headers=['SSH','Username','Hostname','Port','Calling AE',
            'Called AE','Description','Retrieve','Retrieve Data']
            
        self.shortnames={"SSH" : 'ssh', "Username" : 'username',
            "Hostname" : 'host', "Port" : 'port', "Description" : 'description',
            "Calling AE" : 'calling_ae_title', "Called AE" : 'called_ae_title',
            "Retrieve" : "retrieve", "Retrieve Data" : "retrieve_data"}
                  
        # User interface
        self.ui = PreferencesDialog.UI()
        
        xrc_file = medipy.base.find_resource(os.path.join("resources","gui","connection_preferences.xrc"))
        wrappers = []
        medipy.gui.base.Panel.__init__(self, xrc_file, "ConnectionPreferences", 
            wrappers, self.ui, self.ui.controls, parent, *args, **kwargs)

        self.sizer = self.ui.connections.GetSizer()
        
        #Set Events
        self.ui.add.Bind(wx.EVT_BUTTON,self.OnAdd)
        self.list_connections.add_observer("any", self._save_connections)
        self.list_connections.add_observer("any", self._update_list)

        preferences = medipy.gui.base.Preferences(
                wx.GetApp().GetAppName(), wx.GetApp().GetVendorName())

        # An observable list cannot be pickled, so a list is stored in the
        # preferences : fill the content of self.list_connections instead
        # of re-assigning it.
        self.list_connections[:] = preferences.get(self._connections, [])
        self._update_list()
        self.Show(True)
    
    def _update_list(self,*args,**kwargs):
        """ Initialize gui from connection in list_connections
            Add observer on any connection/retrieve panel
        """
        self.sizer.Clear(True)
        self.panels=[]
        
        for row,(description,connection,method,option) in enumerate(self.list_connections):
            
            connection_panel = medipy.gui.network.dicom.Connection(self,connection)
            retrieve_panel = medipy.gui.network.dicom.Retrieve(self,method,option)
            description_text = wx.TextCtrl(self.ui.connections,value=description,name=str(row))
            delete = wx.Button(self.ui.connections,label='Delete',name=str(row))
            validate = wx.Button(self.ui.connections,label='Validate',name=str(row))
            
            box = wx.StaticBox(self.ui.connections, label=description)
            sizer = wx.StaticBoxSizer(box,wx.VERTICAL)
                        
            self.panels.append([connection_panel,
                                retrieve_panel,
                                description_text,
                                delete,
                                validate,
                                box])
            
            connection_panel.add_observer("modify",self._update_connections)
            retrieve_panel.add_observer("modify",self._update_connections)
            description_text.Bind(wx.EVT_TEXT,self._update_connections)
            delete.Bind(wx.EVT_BUTTON,self.OnDelete)
            validate.Bind(wx.EVT_BUTTON,self.OnValidate)
            
            hrz_sizer = wx.BoxSizer(wx.HORIZONTAL)
            hrz_sizer.AddMany([delete,validate])
            
            sizer.AddMany([(description_text,0,wx.EXPAND), (connection_panel,1,wx.EXPAND),
                retrieve_panel, hrz_sizer])
            
            self.sizer.Add(sizer,0,wx.EXPAND)
        
        self.sizer.Layout()
    
    def _update_label(self,event):
        source = event.GetEventObject()
        row = int(source.GetName())
        description = self.panels[row][2].GetValue()
        self.panels[row][5].SetLabel(description)
    
    def _update_connections(self,event=None,*args,**kwargs):
        """ Update and save connections observable list from ui entry
        """
        
        for row,item in enumerate(self.panels):
            self.list_connections[row][0] = item[2].GetValue()
            self.list_connections[row][1] = item[0].connection
            self.list_connections[row][2] = item[1].retrieve[0]
            self.list_connections[row][3] = item[1].retrieve[1]

        self._save_connections()
        if isinstance(event,wx.Event):
            self._update_label(event)
    
    def _save_connections(self, *args) :
        """ Save the connections to the preferences. This function can also be
            used as an event handler.
        """
        preferences = medipy.gui.base.Preferences(
                wx.GetApp().GetAppName(), wx.GetApp().GetVendorName())
        preferences.set(self._connections, self.list_connections[:])

    ##################
    # Event handlers #
    ##################
           
    def OnAdd(self,_):
        """ Add a new connection set with default parameters
        """
        connection =medipy.network.dicom.Connection('----',0,socket.gethostname(),'----')
        self.list_connections.append(["----", connection, "get", ' '])

    def OnDelete(self,event):
        """ Delete selected connection from the connections ObservableList
        """
        source = event.GetEventObject()
        del self.list_connections[int(source.GetName())]

    def OnValidate(self,event):
        """ Test the selected connection
        """
        source = event.GetEventObject()
        row = int(source.GetName())
        
        connection = self.list_connections[row][1]
        
        if isinstance(connection,medipy.network.dicom.SSHTunnelConnection):
            #Ask Password to user
            dlg = wx.PasswordEntryDialog(self,'Enter Your Password','SSH Connection, {0}'.format(connection.user))
            dlg.ShowModal()
            connection.password = dlg.GetValue()
            dlg.Destroy()
        
        connection.connect()
        echo = medipy.network.dicom.scu.Echo(connection)
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
        connection.disconnect()
