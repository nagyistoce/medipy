##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import datetime
import math
import os

import wx

import medipy.gui.base
import medipy.gui.network.dicom
import medipy.gui.control
import medipy.gui.xrc_wrapper
from medipy.gui import PeriodicProgressDialog, WorkerThread
import medipy.base
import medipy.io.dicom
import medipy.io.dicom.misc
import medipy.network.dicom


class QueryDialog(medipy.gui.base.Dialog):

    _connections = "network/dicom/connections"
    _current_connection = "network/dicom/current_connection"
    _queries_fields = "network/dicom/queries"
    _hierarchy = "network/dicom/hierarchy"
  
    tree_headers = [
        'patients_birth_date','patients_sex',
        'modalities_in_study','study_date', 'study_time', 
        "number_of_study_related_series",
        "series_date", "series_time", 
        'modality','number_of_series_related_instances']

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
        self.ui.open.Bind(wx.EVT_BUTTON,self.OnDownLoad)
        self.ui.save.Bind(wx.EVT_BUTTON,self.OnDownLoad)
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
                item = child
            
            # Adjust the columns width
            self.ui.results.ExpandAll(self.root)
            for i in range(self.ui.results.GetColumnCount()):
                self.ui.results.SetColumnWidth(i, wx.LIST_AUTOSIZE)
                width = self.ui.results.GetColumnWidth(i)
                self.ui.results.SetColumnWidth(i, wx.LIST_AUTOSIZE_USEHEADER)
                width = max(width, self.ui.results.GetColumnWidth(i))
                self.ui.results.SetColumnWidth(i, width)
    
    ##########################
    # Tree related functions #
    ##########################
  
    def SetInformations(self, item, level, dataset):
        """ Set informations related to item
            Format into a more readable piece of information (date, hour...)
        """
        
        # Gather the informations according to the level
        informations = []
        
        if level == "patient":
            birth_date = medipy.io.dicom.misc.parse_da(dataset["patients_birth_date"])
            
            informations.append(dataset["patients_sex"].value)
            informations.append(birth_date.strftime("%d/%m/%Y"))
        elif level == "study":
            modalities = dataset['modalities_in_study'].value
            if isinstance(modalities, list) :
                modalities = ", ".join(sorted(modalities))
            
            study_date = medipy.io.dicom.misc.parse_da(dataset["study_date"])
            study_time = medipy.io.dicom.misc.parse_tm(dataset["study_time"])
            study_date_time = datetime.datetime.combine(study_date, study_time)
            
            informations.append(modalities)
            informations.append(study_date_time.strftime("%d/%m/%Y %H:%M"))
            informations.append("{0} series".format(
                dataset["number_of_study_related_series"].value))
        elif level == "series":
            series_date = medipy.io.dicom.misc.parse_da(dataset["series_date"])
            series_time = medipy.io.dicom.misc.parse_tm(dataset["series_time"])
            series_date_time = datetime.datetime.combine(series_date, series_time)
            
            informations.append(dataset["modality"].value)
            informations.append(series_date_time.strftime("%d/%m/%Y %H:%M"))
            informations.append("{0} instances".format(
                dataset["number_of_series_related_instances"].value))
        
        # Set the informations in the tree
        for index, value in enumerate(informations):
            self.ui.results.SetItemText(item, str(value), 1+index)

    def ItemQuery(self,item):
        """ Build query (dictionary) based on treectrl item
            If item has child, recursive call until not
            Return a dataset
        """
        query={}
        while self.ui.results.GetItemText(item)!="Root":
            level,key = self.ui.results.GetItemPyData(item)
            if level!=None:
                if level == "patient":
                    query["patient_id"]=key
                else :
                    query["{0}_instance_uid".format(level)]=key
            item = self.ui.results.GetItemParent(item)

        return medipy.io.dicom.DataSet(**query)

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
        self.pref_dlg = wx.Dialog(self, size = (800,600),
            style=wx.DEFAULT_DIALOG_STYLE|wx.THICK_FRAME)
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
            if not self._get_SSHPasswd(connection):
                return
        
        connection.connect()
        
        list_queries={}
        for key, control in self.query_ctrl.items():
            list_queries[key]=control.GetValue()

        echo = medipy.network.dicom.scu.Echo(connection)
        try:
            echo()
        except medipy.base.Exception, e:
            dlg = wx.MessageDialog(self,
            "Cannot contact entity.\nMake sure you entered the right parameters.",
            'Connection failure',wx.OK|wx.ICON_ERROR)
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
                # Sort the datasets here, since the Sort function of the tree is
                # not customizable
                def datasets_key(dataset):
                    patient = dataset.get("patients_name", 
                        dataset.get("patient_id", medipy.io.dicom.CS(""))).value
                    study = (
                        dataset.get("study_date", medipy.io.dicom.DA("")).value +
                        dataset.get("study_time", medipy.io.dicom.DA("")).value)
                    series = (
                        dataset.get("series_date", medipy.io.dicom.DA("")).value +
                        dataset.get("series_time", medipy.io.dicom.DA("")).value)
                    
                    return (patient, study, series)
                datasets.sort(key=datasets_key)
                self.update_tree(datasets)

        connection.disconnect()
        self._update_download()
    
    def OnDownLoad(self,event):
        """ Ask user for a path directory
            Then save selected objects in query results in a DICOMDIR
            at this location
        """
        source = event.GetEventObject()
        
        if source.GetName()== "save":        
            #Dialog Directory
            dlg = wx.DirDialog(self, style=wx.DD_DIR_MUST_EXIST)
            if dlg.ShowModal() != wx.ID_OK:
                dlg.Destroy()
                return
            path = dlg.GetPath()
            dlg.Destroy()
            
        preferences = medipy.gui.base.Preferences(
                wx.GetApp().GetAppName(), wx.GetApp().GetVendorName())
        _,current = preferences.get(self._current_connection,[None, None])
        connection = current[1]
        retrieve = current[2]
        retrieve_data = current[3]
              
        if isinstance(connection,medipy.network.dicom.SSHTunnelConnection):
            if not self._get_SSHPasswd(connection):
                return
        
        connection.connect()
        
        retrieve_function = getattr(self, "{0}_dl".format(retrieve))
        datasets = retrieve_function(connection,retrieve_data)
        
        if source.GetName()== "save": 
            save = medipy.io.dicom.routing.SaveDataSet(str(path),mode="hierarchical")
            for dataset in datasets:
                save(dataset)
        else :
            self._open(datasets)
        
        dlg = wx.MessageDialog(self, "Successful DownLoad",'Success',
                    wx.OK|wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()
        
        connection.disconnect()
    
    def _get_SSHPasswd(self,connection):
        #Ask Password to user
        dlg = wx.PasswordEntryDialog(self,'Enter Your Password',
                    'SSH Connection, {0}'.format(connection.user))
        if dlg.ShowModal() != wx.ID_OK:
            dlg.Destroy()
            return None
        connection.password = dlg.GetValue()
        dlg.Destroy()
    
    def _open(self,datasets):
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

    def build_wado_query(self,connection):
        """ Build a list of queries based on selected area in ListCtrl
            Return a list of DataSet
        """
        retrieve_query = []
        for item in self.ui.results.GetSelections():
            query = self.ItemQuery(item)
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

    def wado_dl(self,connection,wado_url):
        """ Download data specified in query from wado_url
            Return a list of DataSets
        """
        retrieve_query = self.build_wado_query(connection)
        
        datasets_wado = []
        progress = wx.ProgressDialog(
                    title="Retrieving data from server",
                    message="Downloading data...",
                    maximum=len(retrieve_query),
                    parent=self,
                    style=wx.PD_AUTO_HIDE)

        for index,query in enumerate(retrieve_query):
            value = math.ceil(float(index)/float(len(retrieve_query))*100)
            progress.Update(index,"Retrieving data : {0}%".format(value))
            datasets_wado.append(medipy.network.dicom.wado.get(wado_url,query))
        progress.Destroy()
        return datasets_wado
    
    def GetSelectedUids(self):
        """ Build query based on selected objects in TreeCtrl
            Returns query dataset
        """
        query=[]
        for item in self.ui.results.GetSelections():
            query.append(self.ItemQuery(item))
        
        return query
        
    def GetQueryLevel(self,dataset):
        """ Returns query_level associated to specified dataset
        """
        if "sop_instance_uid" in dataset:
            return "image"
        elif "series_instance_uid" in dataset:
            return "series"
        elif "study_instance_uid" in dataset:
            return "study"
        elif "patient_id" in dataset:
            return "patient"
        
    def move_dl(self,connection,destination):
        """ Move SCU call to send selected objects to specified desination
            Return a list of DataSets
        """
        periodic_progress_dialog = PeriodicProgressDialog(0.2, "Move Retrieve",
                                    "Retrieving Data From Server ...")
        move_query = self.GetSelectedUids()
        
        results = []
        for query in move_query:
            query_level = self.GetQueryLevel(query)
            move = medipy.network.dicom.scu.Move(connection, "patient", query_level,
                destination, query)
            worker_thread = WorkerThread(periodic_progress_dialog,target=move)
            worker_thread.start()
            periodic_progress_dialog.start()
            worker_thread.join()
            
            results.extend(worker_thread.result)

        periodic_progress_dialog.Destroy()

        return results
    
    def get_dl(self,connection,dummy):
        """ Get SCU call, download selected objects
            Return a list of DataSets
        """
        periodic_progress_dialog = PeriodicProgressDialog(0.2, "Get Retrieve", 
                                    "Retrieving Data From Server ...")
        get_query = self.GetSelectedUids()

        results=[]
        for query in get_query:
            query_level = self.GetQueryLevel(query)
            get = medipy.network.dicom.scu.Get(connection, "patient", query_level, query)
            worker_thread = WorkerThread(periodic_progress_dialog,target=get)
            worker_thread.start()
            periodic_progress_dialog.start()
            worker_thread.join()
            
            results.extend(worker_thread.result)
        
        periodic_progress_dialog.Destroy()
        
        return results      
