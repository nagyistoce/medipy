##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import wx
import os

import medipy.network.dicom
import medipy.io.dicom
import medipy.gui.base
import medipy.base

    ########################
    # Retrieve Panel Class #
    ########################

class Retrieve(wx.Panel,medipy.base.Observable):
    def __init__(self, parent=None, retrieve_by="get", retrieve_option='',
                    *args,**kwargs):
        wx.Panel.__init__(self,parent,*args,**kwargs)
        medipy.base.Observable.__init__(self,["modify"])
        
        #retrieve = [retrieve_by,retrieve_option] such as [wado,url_wado]
        self._retrieve=[]
        self.choice = ["wado","move","get"]
        self.sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.SetSizer(self.sizer)
        self._set_retrieve([retrieve_by,retrieve_option])
    
    def modify(self,event):
        """ Event handler
            Modify retrieve object as specified in gui
            Notify modification to observer (in preferences_dialog)
        """
        if self.choice[self.choicebox.GetSelection()]=='get':
            self.sizer.Clear(True)
            self._set_retrieve(["get",' '])
            self.sizer.Layout()
        elif self._retrieve[1]==' ':
            index = self.choicebox.GetSelection()
            self.sizer.Clear(True)
            self._set_retrieve([self.choice[index],'option'])
            self.sizer.Layout()
        else:
            index = self.choicebox.GetSelection()
            self._retrieve = [self.choice[index],self.option.GetValue()]


        self.notify_observers("modify")
    
    def _get_retrieve(self):
        """ Retrieve list, retrieve_by may be either [wado,move,get],
            retrieve_option store usefull information for specified retrieve
            such as wado url or move destination
        """
        return self._retrieve

    def _set_retrieve(self,retrieve):
        self.choicebox = wx.Choice(self,choices=self.choice)
        self.choicebox.SetSelection(self.choice.index(retrieve[0]))
        self.choicebox.Bind(wx.EVT_CHOICE,self.modify)
        self.sizer.Add(self.choicebox,0,wx.EXPAND)
        
        if retrieve[0]!='get':
            self.option = wx.TextCtrl(self,value=retrieve[1])
            self.option.Bind(wx.EVT_TEXT,self.modify)
            self.sizer.Add(self.option,1,wx.EXPAND)

        self._retrieve = retrieve

    retrieve = property(_get_retrieve,_set_retrieve)

    ##########################
    # Connection Panel Class #
    ##########################

class Connection(wx.Panel,medipy.base.Observable):
   
    def __init__(self,parent=None,connection=None,*args,**kwargs):

        wx.Panel.__init__(self,parent,*args,**kwargs)
        medipy.base.Observable.__init__(self,["modify"])

        #connection = medipy.network.dicom.Connection or SSHTunnelConnection Object
        self._connection = None       
        
        self.headers=['SSH','Username','Hostname','Port','Calling AE','Called AE']
            
        self.shortnames={"SSH" : 'ssh', "Username" : 'username',
        "Hostname" : 'host',"Port" : 'port',"Calling AE" : 'calling_ae_title',
        "Called AE" : 'called_ae_title'}

        self.sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.SetSizer(self.sizer)
        self._set_connection(connection) 
        

    def modify(self,event):
        """ Event handler
            Modify connection object as specified in gui
            Notify modification to observer (in preferences_dialog)
        """
        if self.checkbox.GetValue()==True and isinstance(self.connection,medipy.network.dicom.Connection):
            connection = medipy.network.dicom.SSHTunnelConnection(
                remote_host = self.text['Hostname'].GetValue(),
                remote_port = int(self.text['Port'].GetValue()),
                calling_ae_title = self.text['Calling AE'].GetValue(),
                called_ae_title = self.text['Called AE'].GetValue(),
                username = '',password='')
            
            self.sizer.Clear(True)
            self._set_connection(connection)
            self.sizer.Layout()
        
        elif self.checkbox.GetValue()==False and isinstance(self.connection,medipy.network.dicom.SSHTunnelConnection):
            connection = medipy.network.dicom.Connection(
                host = self.text['Hostname'].GetValue(),
                port = int(self.text['Port'].GetValue()),
                calling_ae_title = self.text['Calling AE'].GetValue(),
                called_ae_title = self.text['Called AE'].GetValue())
            
            self.sizer.Clear(True)
            self._set_connection(connection)
            self.sizer.Layout()
            
        else:
            if self.checkbox.GetValue()==True:
                self._connection.user = self.user.GetValue()
            for header in self.headers[2:]:
                name = self.shortnames[header]
                value = self.text[header].GetValue()
                if name == "port" and value!='':
                    self._connection.__setattr__(name,int(value))
                else:
                    self._connection.__setattr__(name,value)

        self.notify_observers("modify")
    
    ##############
    # Properties #
    ##############
    
    def _get_connection(self):
        """ Connection object, may be either a simple medipy.network.dicom.Connection
            or a medipy.network.dicom.SSHTunnelConnection
        """
        return self._connection
    
    def _set_connection(self,connection):

        self.checkbox = wx.CheckBox(self)
        self.checkbox.Bind(wx.EVT_CHECKBOX,self.modify)
        self.sizer.Add(wx.StaticText(self,label='SSH :'),0,wx.ALIGN_CENTER)
        self.sizer.Add(self.checkbox,0,wx.EXPAND)

        if isinstance(connection,medipy.network.dicom.SSHTunnelConnection):
            self.checkbox.SetValue(True)
            self.user = wx.TextCtrl(self,value=connection.user)
            self.user.Bind(wx.EVT_TEXT,self.modify)
            self.sizer.Add(wx.StaticText(self,label='Username :'),0,wx.ALIGN_CENTER)
            self.sizer.Add(self.user,1,wx.EXPAND)
        else:
            self.checkbox.SetValue(False)

        self.text={}
        for header in self.headers[2:]:
            name=self.shortnames[header]
            self.text[header] = wx.TextCtrl(self, value=str(connection.__getattribute__(name)))
            self.text[header].Bind(wx.EVT_TEXT,self.modify)
            if header!="Port":
                self.sizer.Add(wx.StaticText(self,label=header+' :'),
                    0,wx.ALIGN_CENTER)
                self.sizer.Add(self.text[header],1,wx.EXPAND)
            else :
                self.sizer.Add(wx.StaticText(self,label=':'),0,wx.ALIGN_CENTER)
                self.sizer.Add(self.text[header],0,wx.EXPAND)

        self._connection = connection

    connection = property(_get_connection, _set_connection)
