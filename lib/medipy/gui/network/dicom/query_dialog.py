##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################
import wx
import os

import medipy.gui.base
import medipy.base

import medipy.gui.control
import medipy.io.dicom
import medipy.network.dicom
import medipy.gui.base
import medipy.base
import medipy.gui.network.dicom

class QueryDialog(medipy.gui.base.Panel):

    _connections = "network/dicom/connections"
    _current_connection = "network/dicom/current_connection"
    _ssh_connections = "network/dicom/ssh"
    _queries_fields = "network/dicom/queries"
    
    tree_headers = ['acquisition_date','acquisition_time', 'modality']
    tree_levels = ['patient_id', 'study_description', 'series_description']

    class UI(medipy.gui.base.UI):
        def __init__(self):
            self.qu_panel = None
            self.search = None
            self.download = None
            self.directory = None
            self.preferences = None
            self.selected_connection = None
            self.results = None
            self.radio = None
            self.set_queries = None
            self.ssh = None
                        
            self.controls = ["qu_panel","search","download","directory",
                "preferences","selected_connection","results","radio","set_queries",
                "ssh"]
            
    def __init__(self, parent=None, *args, **kwargs):
        
        self.connection = None
        self.query_ctrl={}      #Dictionary for wx.controls
        
        # User interface
        self.ui = QueryDialog.UI()
        
        xrc_file = medipy.base.find_resource(os.path.join("resources", "gui", "query_dialog.xrc"))
        wrappers = []
        medipy.gui.base.Panel.__init__(self, xrc_file, "Search and Download", 
            wrappers, self.ui, self.ui.controls, parent, *args, **kwargs)

        #Set Controls
        self.destination = medipy.gui.control.Directory(self.ui.directory)
        sizer=self.ui.directory.GetSizer()
        sizer.Add(self.destination,1,wx.EXPAND|wx.BOTTOM)
        self.queries = self.ui.qu_panel.GetSizer()
        
        self.tree = wx.gizmos.TreeListCtrl(self,name='Patient',
            style=wx.TR_DEFAULT_STYLE|wx.TR_HIDE_ROOT|wx.TR_FULL_ROW_HIGHLIGHT|wx.TR_EXTENDED|wx.TR_MULTIPLE)
        treesizer=self.ui.results.GetSizer()
        treesizer.Add(self.tree,1,wx.EXPAND)
        
        #Set Events
        self.ui.search.Bind(wx.EVT_BUTTON,self.OnSearch)
        self.ui.download.Bind(wx.EVT_BUTTON,self.OnDownLoad)
        self.destination.add_observer("value", self.OnDestination)
        self.ui.results.Bind(wx.EVT_LIST_ITEM_SELECTED,self.OnResult)
        self.ui.results.Bind(wx.EVT_LIST_ITEM_DESELECTED,self.OnResult)
        self.ui.preferences.Bind(wx.EVT_BUTTON,self.OnPreferences)
        self.ui.selected_connection.Bind(wx.EVT_CHOICE,self.OnChoice)
        self.ui.radio.Bind(wx.EVT_RADIOBOX,self._update_download)
        self.ui.set_queries.Bind(wx.EVT_BUTTON,self.OnSetQueries)
        self.ui.ssh.Bind(wx.EVT_CHECKBOX,self._update_choice)

        self._update_choice()
        self._update_queries()
        self.update_tree_column()
        
        self.Show(True)
    
    #------------------
    #   Update GUI
    #------------------
        
    def _update_queries(self):

        self.queries.Clear(True)
        self.query_ctrl = {}
        self.update_tree()
        
        #Load query preferences
        preferences = medipy.gui.base.Preferences(
                wx.GetApp().GetAppName(), wx.GetApp().GetVendorName())
        self.fields = preferences.get(self._queries_fields,[])
        
        if self.fields == []:
            self.fields = ['patients_name','series_description','study_description']
        self.queries.SetRows(len(self.fields))
                    
        for field in self.fields:
            tag = medipy.io.dicom.dictionary.name_dictionary[field]
            label = medipy.io.dicom.dictionary.data_dictionary[tag][2]
            self.query_ctrl[field] = wx.TextCtrl(self.ui.qu_panel)
            self.queries.Add(wx.StaticText(self.ui.qu_panel,label=label),
                    flag=wx.ALIGN_CENTER_VERTICAL)
            self.queries.Add(self.query_ctrl[field],proportion=1,
                    flag=wx.EXPAND|wx.ALIGN_CENTER_VERTICAL)
            
            self.queries.Layout()
        
    def _update_choice(self, *args):
        
        preferences = medipy.gui.base.Preferences(
                wx.GetApp().GetAppName(), wx.GetApp().GetVendorName())
        self.ui.selected_connection.Clear()
        
        (choice,self.connection,_,__) = preferences.get(self._current_connection,[])
        
        if not self.ui.ssh.GetValue():
            list_connections = preferences.get(self._connections, [])        
            for connection in list_connections:
                self.ui.selected_connection.Append(connection[1].host+' --- '+
                        str(connection[1].port)+' --- '+connection[0])

        else:
            list_connections = preferences.get(self._ssh_connections, [])
            for connection in list_connections:
                self.ui.selected_connection.Append(connection[1].host+' --- '+
                    str(connection[1].port)+' --- '+connection[1].username+' --- '
                    +connection[0])
        
        self.ui.selected_connection.SetSelection(choice)
            
        
    def _update_download(self, *args, **kwargs):
        self.ui.download.Enable(
            (self.destination.validate() or self.ui.radio.GetSelection()==1)
            and self.tree.GetSelections()!=-1)

    def update_tree_column(self):
        self.tree.AddColumn("Patient Tree",width=250)
        for header in self.tree_headers:
            tag = medipy.io.dicom.dictionary.name_dictionary[header]
            label = medipy.io.dicom.dictionary.data_dictionary[tag][2]
            self.tree.AddColumn(label)

    def update_tree(self,datasets=[]):
        self.tree.DeleteAllItems()
        self.root = self.tree.AddRoot(text='Patient')
        for dataset in datasets:
            patient = dataset['patient_id'].value
            study = dataset['study_description'].value
            serie = dataset['series_description'].value
        
            if self.tree.ItemHasChildren(self.root) :
                patient_found,patient_item = self.IsChild(self.root,patient)
                    
                if not patient_found :
                    patient_item = self.tree.AppendItem(self.root,text=patient)
                    study_item = self.tree.AppendItem(patient_item,text=study)
                    serie_item = self.tree.AppendItem(study_item,text=serie)
                                            
                else :
                    study_found,study_item = self.IsChild(patient_item,study)
                    if not study_found :
                        study_item = self.tree.AppendItem(patient_item,text=study)
                        serie_item = self.tree.AppendItem(study_item,text=serie)
                                                    
                    else:
                        serie_found,serie_item = self.IsChild(study_item,serie)
                        if not serie_found :
                            serie_item = self.tree.AppendItem(study_item,text=serie)
            else:
                patient_item = self.tree.AppendItem(self.root,text=patient)
                study_item = self.tree.AppendItem(patient_item,text=study)
                serie_item = self.tree.AppendItem(study_item,text=serie)
            
            for row,header in enumerate(self.tree_headers):       
                self.tree.SetItemText(serie_item,dataset[header].value,row+1)
            
            self.tree.SortChildren(self.root)
            self.tree.SortChildren(patient_item)
            self.tree.SortChildren(study_item)
            self.tree.SortChildren(serie_item)

    def GetItemTreeLevel(self,item):
        """ Count every parent that exists over item
            Return count.
        """
        count = 0
        while self.tree.GetItemText(item)!="Patient": #Text associated to root
            item = self.tree.GetItemParent(item)
            count = count+1
            
        return count

    def IsChild(self,itemid,text):
        """ Search text in item.text of any child in itemid
            Return boolean and focused item if found
        """
        item,cookie = self.tree.GetFirstChild(itemid)
        count = self.tree.GetChildrenCount(itemid,False)
        found = False
        for index in range(count):
            if text == self.tree.GetItemText(item):
                found = True
                break
            else:
                item,cookie = self.tree.GetNextChild(itemid,cookie)
        if found == False:
            item = None
            
        return found,item
            
    #------------------
    #   Event handlers
    #------------------
    def OnResult(self, _):
        self._update_download()
    
    def OnDestination(self, _):
        self._update_download()

    def OnSetQueries(self,_):
        self.quer_dlg = wx.Dialog(self,size=(250,300))
        self.quer_panel = medipy.gui.network.dicom.SetQueriesDialog(self.quer_dlg)
        sizer = wx.BoxSizer()
        sizer.Add(self.quer_panel,1,wx.EXPAND)
        
        self.quer_dlg.ShowModal()
        self.quer_dlg.Destroy()
        
        self._update_queries()

    def OnChoice(self,_):
        
        preferences = medipy.gui.base.Preferences(
                wx.GetApp().GetAppName(), wx.GetApp().GetVendorName())
                
        if not self.ui.ssh.GetValue():
            list_connections = preferences.get(self._connections, [])        
        else :
            list_connections = preferences.get(self._ssh_connections, [])
                   
        choice = self.ui.selected_connection.GetCurrentSelection()
        self.connection = list_connections[choice][1]
        retrieve = list_connections[choice][2]
        retrieve_data = list_connections[choice][3]
        preferences.set(self._current_connection,[choice, self.connection, retrieve,
                retrieve_data])

    def OnPreferences(self,_):
       
        self.pref_dlg = wx.Dialog(self,size=(900,400),
                    style=wx.DEFAULT_DIALOG_STYLE|wx.THICK_FRAME)
        if not self.ui.ssh.GetValue():
            self.pref_panel = medipy.gui.network.dicom.PreferencesDialog(self.pref_dlg)
        else :
            self.pref_panel = medipy.gui.network.dicom.SSHDialog(self.pref_dlg)

        sizer = wx.BoxSizer()
        sizer.Add(self.pref_panel, 1, wx.EXPAND)
        self.pref_dlg.SetSizer(sizer)
        
        self.pref_dlg.ShowModal()
        self.pref_dlg.Destroy()
        
        self._update_choice()
        
    def OnSearch(self,_):
        """ Send specified query on dicom.Dataset and show results in ListCtrl
        """
        list_queries={}
        for key, control in self.query_ctrl.items():
            list_queries[key]=control.GetValue()

        echo = medipy.network.dicom.scu.Echo(self.connection)
        try:
            echo()
        except medipy.base.Exception, e:
            dlg = wx.MessageDialog(self, "Cannot contact entity.\nMake sure you entered the right parameters.",'Connection failure',wx.OK|wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
        else:
            query = medipy.io.dicom.DataSet(**list_queries)
            
            for key in ["patient_id", "study_description","series_description"] :
                query.setdefault(key,"")
            
            for key in self.tree_headers:
                query.setdefault(key,"")
                
            datasets =  medipy.network.dicom.query.relational(
                        self.connection,"patient","patient",query)

            if datasets==[]:
                dlg = wx.MessageDialog(self, "No Results",'Try again',wx.OK|wx.ICON_INFORMATION)
                dlg.ShowModal()
                dlg.Destroy()
            else:
                self.update_tree(datasets)

    
    def OnDownLoad(self,_):
        """ DownLoad selected object in ListCtrl
            A path should be specified with control.Directory
        """

        preferences = medipy.gui.base.Preferences(
                wx.GetApp().GetAppName(), wx.GetApp().GetVendorName())
        (_,self.connection,retrieve,retrieve_data) = preferences.get(self._current_connection, []) 

        query = self.build_retrieve_query()
        if retrieve == 'wado':
            datasets = self.wado_dl(retrieve_data,query)
        elif retrieve == 'move':
            datasets = self.move_dl(retrieve_data,query)
        elif retrieve == 'get':
            pass

        if self.ui.radio.GetSelection()==0:
            save = medipy.io.dicom.routing.SaveDataSet(
                                str(self.destination.value),mode="hierarchical")
            for dataset in datasets:
                save(dataset)
        else:
            datasets_dict={}
            
            for item in query:
                key = item.patient_id.value + '_' +                         item.study_instance_uid.value + '_' + item.series_instance_uid.value
                for dataset in datasets:
                    if item.sop_instance_uid.value == dataset.sop_instance_uid.value:
                        if not key in datasets_dict:
                            datasets_dict[key]=[dataset]
                        else :
                            datasets_dict[key].append(dataset)

            for key in datasets_dict:
                images = []
                layers = []
                images.append(medipy.io.dicom.image(datasets_dict[key]))
                for image in images:
                    layers.append({"image":image})
                wx.GetApp().frame.append_image(layers)
                            
        dlg = wx.MessageDialog(self, "Successful DownLoad",'Success',
                    wx.OK|wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()
            
            
    #------------------
    #   Retrieve
    #------------------

    def build_retrieve_query(self):
        """ Build a list of queries based on selected area in ListCtrl
        """
        retrieve_query = []
        
        for item in self.tree.GetSelections():
            list_queries={}
            level = self.GetItemTreeLevel(item)
            for each in reversed(self.tree_levels[:level]):
                list_queries[each]=self.tree.GetItemText(item)
                item = self.tree.GetItemParent(item)
            
            #Query for any useful uid
            query = medipy.io.dicom.DataSet(**list_queries)
            for key in ["patient_id", "study_instance_uid", 
                            "series_instance_uid", "sop_instance_uid"] :
                query.setdefault(key, None)
            datasets =  medipy.network.dicom.query.relational(
                    self.connection,"patient","patient",query)
            for dataset in datasets:
                retrieve_query.append(medipy.io.dicom.DataSet(
                    patient_id = dataset.patient_id.value,
                    study_instance_uid = dataset.study_instance_uid.value,
                    series_instance_uid = dataset.series_instance_uid.value,
                    sop_instance_uid = dataset.sop_instance_uid.value))
                        
        return retrieve_query

    def wado_dl(self,wado_url,retrieve_query):
        datasets_wado = []
        for query in retrieve_query:
            datasets_wado.append(medipy.network.dicom.wado.get(wado_url,query))
            
        return datasets_wado
    
    def move_dl(self,destination,retrieve_query):
        move_query = medipy.io.dicom.DataSet(sop_instance_uid='')
        for query in retrieve_query:
            sop_uid = str(query.sop_instance_uid.value)
            mv_sop = str(move_query.sop_instance_uid.value) + '\\' + sop_uid
            move_query.__setattr__('sop_instance_uid',mv_sop)
            
        move = medipy.network.dicom.scu.Move(self.connection, "patient", "image",
                destination, move_query)
        results = move()
        return results
