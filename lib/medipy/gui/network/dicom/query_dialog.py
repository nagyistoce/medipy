##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################
import os
import wx

import medipy.gui.base
import medipy.gui.network.dicom
import medipy.gui.control
import medipy.gui.xrc_wrapper
import medipy.base
import medipy.io.dicom
import medipy.network.dicom


class QueryDialog(medipy.gui.base.Dialog):

    _connections = "network/dicom/connections"
    _current_connection = "network/dicom/current_connection"
    _queries_fields = "network/dicom/queries"
    _hierarchy = "network/dicom/hierarchy"
  
    tree_headers = ['patients_birth_date','patients_sex',
    'modalities_in_study','study_date',
    'study_time','modality','number_of_series_related_instances']

    class UI(medipy.gui.base.UI):
        def __init__(self):
            self.node = None
            self.patient_based = None
            self.trial_based = None
            self.edit_nodes = None
            self.edit_elements = None
            self.elements_container = None
            self.search = None
            self.results = None
            self.cancel = None
            self.save = None
            self.open = None
                        
            self.controls = ["node", "patient_based", "trial_based",
                "edit_nodes", "edit_elements", "elements_container",
                "search", "results", "cancel", "save", "open"]
            
    def __init__(self, parent=None, *args, **kwargs):        
        self.query_ctrl={}      #Dictionary for wx.controls
        
        # User interface
        self.ui = QueryDialog.UI()
        
        xrc_file = medipy.base.find_resource(os.path.join("resources", "gui", "dicom_qr.xrc"))
        wrappers = [medipy.gui.xrc_wrapper.TreeListCtrlXmlHandler()]
        medipy.gui.base.Dialog.__init__(self, xrc_file, "dicom_qr", 
            wrappers, self.ui, self.ui.controls, parent, *args, **kwargs)

        self.ui.cancel.SetId(wx.ID_CANCEL)
        
        #Set Events
        self.ui.search.Bind(wx.EVT_BUTTON,self.OnSearch)
        self.ui.open.Bind(wx.EVT_BUTTON,self.OnOpen)
        self.ui.save.Bind(wx.EVT_BUTTON,self.OnSave)
        self.ui.results.Bind(wx.EVT_TREE_SEL_CHANGED,self._update_download)
        self.ui.edit_nodes.Bind(wx.EVT_BUTTON,self.OnEditNodes)
        self.ui.node.Bind(wx.EVT_CHOICE,self.OnNode)
        self.ui.edit_elements.Bind(wx.EVT_BUTTON,self.OnEditElements)
        self.ui.patient_based.Bind(wx.EVT_RADIOBUTTON,self.OnRadio)
        self.ui.trial_based.Bind(wx.EVT_RADIOBUTTON,self.OnRadio)

        self._update_hierarchy()
        self._update_choice()
        self._update_queries()
        self.update_tree_column()
    
    ##############
    # GUI Update #
    ##############
    
    def  _update_hierarchy(self):
        preferences = medipy.gui.base.Preferences(
                wx.GetApp().GetAppName(), wx.GetApp().GetVendorName())
        
        hierarchy = preferences.get(self._hierarchy,None)
        
        if hierarchy!=None:
            getattr(self.ui, "{0}_based".format(hierarchy)).SetValue(1)
        else :
            self.ui.patient_based.SetValue(1)
        
    def _update_queries(self):
        """ Query fields update based on stored queries preferences
        """
        #TODO : Catch exception when an unknown label is entered for name_dictionary
        
        queries = self.ui.elements_container.GetSizer()
        queries.Clear(True)
        self.query_ctrl = {}
        self.update_tree()
        
        #Load query preferences
        preferences = medipy.gui.base.Preferences(
                wx.GetApp().GetAppName(), wx.GetApp().GetVendorName())
        self.fields = preferences.get(self._queries_fields,[])
        
        if self.fields == []:
            self.fields = ['patients_name','series_description','study_description']
        queries.SetRows(len(self.fields))
                    
        for field in self.fields:
            tag = medipy.io.dicom.dictionary.name_dictionary[field]
            label = medipy.io.dicom.dictionary.data_dictionary[tag][2]
            self.query_ctrl[field] = wx.TextCtrl(self.ui.elements_container)
            queries.Add(wx.StaticText(self.ui.elements_container,label=label+" :"),
                    flag=wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_RIGHT)
            queries.Add(self.query_ctrl[field],proportion=1,
                    flag=wx.EXPAND|wx.ALIGN_CENTER_VERTICAL)
            
        queries.Layout()
        self.Fit()
        
    def _update_choice(self, *args):
        """ Connection list update based on stored preferences
        """
        preferences = medipy.gui.base.Preferences(
                wx.GetApp().GetAppName(), wx.GetApp().GetVendorName())
        self.ui.node.Clear()

        choice,_ = preferences.get(self._current_connection, [None, None])

        list_connections = preferences.get(self._connections, [])
        if list_connections != []:
            for connection in list_connections:
                if isinstance(connection[1],medipy.network.dicom.SSHTunnelConnection):
                    self.ui.node.Append(connection[1].remote_host+' --- '+
                        str(connection[1].remote_port)+' --- '+connection[0])
                else :
                    self.ui.node.Append(connection[1].host+' --- '+
                        str(connection[1].port)+' --- '+connection[0])

        if choice!=None :
            self.ui.node.SetSelection(int(choice))
            self._update_download()
            self.OnNode()

    def _update_download(self, *args, **kwargs):
        preferences = medipy.gui.base.Preferences(
                wx.GetApp().GetAppName(), wx.GetApp().GetVendorName())
        _,connection = preferences.get(self._current_connection, [None, None])

        self.ui.open.Enable(self.ui.results.GetSelections()!=[]
            and connection!=None)
        self.ui.save.Enable(self.ui.results.GetSelections()!=[]
            and connection!=None)

    def update_tree_column(self):    
        self.ui.results.AddColumn("Patient Tree",width=150)
        self.ui.results.AddColumn("",width=65)
        self.ui.results.AddColumn("",width=65)
        self.ui.results.AddColumn("",width=20)
    
    def update_tree(self,datasets=[]):
        """ TreeCtrl update based on selected view (study or patient):
        """
        self.ui.results.DeleteAllItems()
        self.root = self.ui.results.AddRoot(text='Root')
        for dataset in datasets:
            # Build the hierarchy as a list of (level, label, key)
            if self.ui.patient_based.GetValue():
                hierarchy = [
                    ("patient", dataset.patients_name.value, dataset.patient_id.value),
                    ("study", dataset.study_description.value, dataset.study_instance_uid.value),
                    ("series", dataset.series_description.value, dataset.series_instance_uid.value),
                ]
            else :
                trial, time_point = dataset.study_description.value.rsplit("^", 1)
                
                hierarchy = [(None, item, item) for item in trial.split("^")]
                hierarchy.extend([
                    ("patient", dataset.patients_name.value, dataset.patient_id.value),
                    ("study", time_point, dataset.study_instance_uid.value),
                    ("series", dataset.series_description.value, dataset.series_instance_uid.value),
                ])
            
            item = self.root
            for level, label, key in hierarchy : 
                found, child = self.IsChild(item, (level,key))
                if not found :
                    # Add the node in the tree
                    child = self.ui.results.AppendItem(item, label)
                    self.ui.results.SetItemPyData(child, (level,key))
                    self.SetInformations(child, level, dataset)
                    self.ui.results.SortChildren(item)
                item = child
            
        self.ui.results.SortChildren(self.root)
    
    ##########################
    # Tree related functions #
    ##########################
  
    def SetInformations(self, item, level, dataset):
        """ Set informations related to item
            Format into a more readable piece of information (date, hour...)
        """
        if level == "patient":
            date = str(dataset['patients_birth_date'].value)
            date = date[6:]+'/'+date[4:6]+'/'+date[:4]
            self.ui.results.SetItemText(item,str(dataset['patients_sex'].value),1)
            self.ui.results.SetItemText(item,date,2)
        
        elif level == "study":
            date = str(dataset['study_date'].value)
            date = date[6:]+'/'+date[4:6]+'/'+date[:4]
            hour = str(dataset['study_time'].value)
            hour = hour[:2]+':'+hour[2:4]+':'+hour[4:6]
            self.ui.results.SetItemText(item,date,1)
            self.ui.results.SetItemText(item,hour,2)
            
            modalities = dataset['modalities_in_study'].value
            if isinstance(modalities, list) :
                modalities = ", ".join(sorted(modalities))
            self.ui.results.SetItemText(item, modalities,3)
        
        elif level == "series":
            self.ui.results.SetItemText(item,
               str(dataset['number_of_series_related_instances'].value)+' instances',1)
            self.ui.results.SetItemText(item,str(dataset['modality'].value),2)

    def ItemQuery(self,item):
        """ Build query (dictionary) based on treectrl item
            If item has child, recursive call until not
            Return a list of queries
        """
        queries=[]
        if self.ui.results.ItemHasChildren(item):
            count = self.ui.results.GetChildrenCount(item,False)
            child,cookie = self.ui.results.GetFirstChild(item)
            for index in range(count):
                results = self.ItemQuery(child)
                queries = queries + results
                child,cookie = self.ui.results.GetNextChild(item,cookie)
        else:
            query={}
            while self.ui.results.GetItemText(item)!="Root":
                level,key = self.ui.results.GetItemPyData(item)
                if level!=None:
                    if level == "patient":
                        query["patient_id"]=key
                    else :
                        query["{0}_instance_uid".format(level)]=key
                item = self.ui.results.GetItemParent(item)
            queries.append(query)
        return queries

    def IsChild(self,itemid,key):
        """ Search text in item.text of any child in itemid
            Return boolean and focused item if found, None if not
        """
        item,cookie = self.ui.results.GetFirstChild(itemid)
        count = self.ui.results.GetChildrenCount(itemid,False)
        found = False
        for index in range(count):
            if key == self.ui.results.GetItemPyData(item):
                found = True
                break
            else:
                item,cookie = self.ui.results.GetNextChild(itemid,cookie)
        if found == False:
            item = None
            
        return found,item
            
    ##################
    # Event handlers #
    ##################  
    def OnEditElements(self,_):
        self.quer_dlg = wx.Dialog(self,size=(250,300))
        self.quer_panel = medipy.gui.network.dicom.SetQueriesDialog(self.quer_dlg)
        sizer = wx.BoxSizer()
        sizer.Add(self.quer_panel,1,wx.EXPAND)
        
        self.quer_dlg.ShowModal()
        self.quer_dlg.Destroy()
        
        self._update_queries()

    def OnNode(self,*args,**kwargs):
        """ Store current connection in preferences (index)
        """
        preferences = medipy.gui.base.Preferences(
                wx.GetApp().GetAppName(), wx.GetApp().GetVendorName())
        
        list_connections = preferences.get(self._connections,[])
        choice = self.ui.node.GetCurrentSelection()
        
        if choice == -1 :
            return
        
        preferences.set(self._current_connection,(choice,list_connections[choice]))

    def OnEditNodes(self,_):       
        self.pref_dlg = wx.Dialog(self,style=wx.DEFAULT_DIALOG_STYLE|wx.THICK_FRAME)
        self.pref_panel = medipy.gui.network.dicom.PreferencesDialog(self.pref_dlg)

        sizer = wx.BoxSizer()
        sizer.Add(self.pref_panel, 1, wx.EXPAND)
        self.pref_dlg.SetSizer(sizer)
        
        self.pref_dlg.ShowModal()
        self.pref_dlg.Destroy()
        
        self._update_choice()
        
    def OnSearch(self,_):
        """ Use relational to retrieve specified query
            Call update_tree to show results
        """
        
        preferences = medipy.gui.base.Preferences(
                wx.GetApp().GetAppName(), wx.GetApp().GetVendorName())
        _,current = preferences.get(self._current_connection,[None, None])
        connection = current[1]

        if isinstance(connection,medipy.network.dicom.SSHTunnelConnection):
            #Ask Password to user
            dlg = wx.PasswordEntryDialog(self,'Enter Your Password','SSH Connection : {0}'.format(connection.user))
            dlg.ShowModal()
            connection.password = dlg.GetValue()
            dlg.Destroy()
        
        connection.connect()
        
        list_queries={}
        for key, control in self.query_ctrl.items():
            list_queries[key]=control.GetValue()

        echo = medipy.network.dicom.scu.Echo(connection)
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
                        connection,"patient","patient",query)

            if datasets==[]:
                dlg = wx.MessageDialog(self, "No Results",'Try again',wx.OK|wx.ICON_INFORMATION)
                dlg.ShowModal()
                dlg.Destroy()
            else:
                self.update_tree(datasets)

        connection.disconnect()
        self._update_download()
    
    def OnSave(self,_):
        """ Ask user for a path directory
            Then save selected objects in query results in a DICOMDIR
            at this location
        """
        #Dialog Directory
        dlg = wx.DirDialog(self, style=wx.DD_DIR_MUST_EXIST)
        dlg.ShowModal()
        path = dlg.GetPath()
        dlg.Destroy()
        
        preferences = medipy.gui.base.Preferences(
                wx.GetApp().GetAppName(), wx.GetApp().GetVendorName())
        _,current = preferences.get(self._current_connection,[None, None])
        connection = current[1]
        retrieve = current[2]
        retrieve_data = current[3]
        
        if isinstance(connection,medipy.network.dicom.SSHTunnelConnection):
            #Ask Password to user
            dlg = wx.PasswordEntryDialog(self,'Enter Your Password','SSH Connection, {0}'.format(connection.user))
            dlg.ShowModal()
            connection.password = dlg.GetValue()
            dlg.Destroy()

        connection.connect()
        
        query = self.build_retrieve_query(connection)
        retrieve_function = getattr(self, "{0}_dl".format(retrieve))
        datasets = retrieve_function(connection,retrieve_data, query)
        
        save = medipy.io.dicom.routing.SaveDataSet(str(path),mode="hierarchical")
        for dataset in datasets:
            save(dataset)
        
        dlg = wx.MessageDialog(self, "Successful DownLoad",'Success',
                    wx.OK|wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()
        
        connection.disconnect()
        
    def OnOpen(self,_):
        """ Open with MediPy selected objects in query results
        """
        preferences = medipy.gui.base.Preferences(
                wx.GetApp().GetAppName(), wx.GetApp().GetVendorName())
        _,current = preferences.get(self._current_connection,[None, None])
        connection = current[1]
        retrieve = current[2]
        retrieve_data = current[3]
        
        if isinstance(connection,medipy.network.dicom.SSHTunnelConnection):
            #Ask Password to user
            dlg = wx.PasswordEntryDialog(self,'Enter Your Password','SSH Connection, {0}'.format(connection.user))
            dlg.ShowModal()
            connection.password = dlg.GetValue()
            dlg.Destroy()

        connection.connect()
        
        query = self.build_retrieve_query(connection)
        retrieve_function = getattr(self, "{0}_dl".format(retrieve))
        datasets = retrieve_function(connection,retrieve_data, query)
        
        series = medipy.io.dicom.series(datasets)
        for serie in series:
            stacks = medipy.io.dicom.split.stacks(serie)
            # Display dialog
            if len(stacks)>1:
                dialog = medipy.gui.dicom.StacksDialog(self,False)
                dialog.set_stacks(stacks)
                if dialog.ShowModal() != wx.ID_OK :
                    stacks=[]
                stacks = dialog.get_selected_stacks()
            images = [medipy.io.dicom.image(stack) for stack in stacks]
            wx.GetApp().frame.append_image([{"image":image} for image in images])
                        
        dlg = wx.MessageDialog(self, "Successful DownLoad",'Success',
                    wx.OK|wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()
        
        connection.disconnect()
    
    def OnRadio(self,_):
        preferences = medipy.gui.base.Preferences(
            wx.GetApp().GetAppName(), wx.GetApp().GetVendorName())
            
        if self.ui.patient_based.GetValue():
            hierarchy = "patient"
        else :
            hierarchy = "trial"
        preferences.set(self._hierarchy, hierarchy)
    
    ############
    # Retrieve #
    ############

    def build_retrieve_query(self,connection):
        """ Build a list of queries based on selected area in ListCtrl
            Return a list of DataSet
        """
        retrieve_query = []
        
        for item in self.ui.results.GetSelections():
            list_queries = self.ItemQuery(item)
            for queries in list_queries:
                #Query for any useful uid
                query = medipy.io.dicom.DataSet(**queries)
                for key in ["patient_id", "study_instance_uid", 
                                "series_instance_uid", "sop_instance_uid"] :
                    query.setdefault(key, None)
                datasets =  medipy.network.dicom.query.relational(
                        connection,"patient","patient",query)
                for dataset in datasets:
                    retrieve_query.append(medipy.io.dicom.DataSet(
                        patient_id = dataset.patient_id.value,
                        study_instance_uid = dataset.study_instance_uid.value,
                        series_instance_uid = dataset.series_instance_uid.value,
                        sop_instance_uid = dataset.sop_instance_uid.value))
                        
        return retrieve_query

    def wado_dl(self,connection,wado_url,retrieve_query):
        """ Download data specified in query from wado_url
            Return a list of DataSets
        """
        datasets_wado = []
        for query in retrieve_query:
            datasets_wado.append(medipy.network.dicom.wado.get(wado_url,query))
            
        return datasets_wado
    
    def move_dl(self,connection,destination,retrieve_query):
        """ Move SCU call to download specified query to desination
            Return a list of DataSets
        """
        move_query = medipy.io.dicom.DataSet(
            sop_instance_uid="\\".join(
                x.sop_instance_uid.value for x in retrieve_query))
            
        move = medipy.network.dicom.scu.Move(connection, "patient", "image",
                destination, move_query)
        results = move()
        return results
