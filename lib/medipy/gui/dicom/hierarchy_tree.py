##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import datetime

import wx

import medipy.io.dicom
import medipy.io.dicom.misc

class HierarchyTree(wx.TreeCtrl):
    """ Hierarchical tree of Patients, Studies and Series
    """
    
    def __init__(self, parent, ID=wx.ID_ANY, pos=wx.DefaultPosition, 
                 size=wx.DefaultSize, style=wx.TR_HIDE_ROOT|wx.TR_HAS_BUTTONS, 
                 validator=wx.DefaultValidator, name="hierarchytree"):
        
        wx.TreeCtrl.__init__(self, parent, ID, pos, size, style, validator, name)
    
    def set_datasets(self, datasets):
        """ Set the list of datasets to be displayed in the tree
        """
        
        patients = self._get_patients(datasets)
        
        self.DeleteAllItems()
        root = self.AddRoot("")
        
        # Fill the tree
        for patient in sorted(patients.keys()) :
            patient_item = self.AppendItem(root, patient)
            
            studies = sorted(patients[patient].items(), self._study_comparator)
            for study_uid, study in studies :
                study_item = self.AppendItem(patient_item, study["description"])
                
                series = sorted(patients[patient][study_uid]["series"].items(), 
                                self._series_comparator)
                for _, serie in series :
                    series_item = self.AppendItem(study_item, serie["description"])
                    self.SetItemPyData(series_item, serie["datasets"])
        
        self.ExpandAll()
    
    ##############
    # Properties #
    ##############
    
    def _get_selected_datasets(self):
        """ Return the datasets corresponding to the selected series.
        """
        
        item = self.GetSelection()
        if item.IsOk() :
            return self.GetItemPyData(item)
        else :
            return None
    
    selected_datasets = property(_get_selected_datasets)
    
    #####################
    # Private interface #
    #####################
    
    @staticmethod
    def _get_patients(datasets):
        """ Return the patient/study/series hierarchy contained in the datasets.
            The following structure is returned :
            
            {
                patient_name_or_id :
                {
                    study_instance_uid :
                    {
                        "date" : study date and time or None,
                        "description" : study description or None,
                        "series" :
                        {
                            series_instance_uid :
                            {
                                "date" : series date and time or None,
                                "description" : series description or None,
                                "datasets" : [...]
                            }
                            ... other series
                        }
                    }
                    ... other studies
                }
                ... other patients
            }
        """
        
        patients = {}
        
        # Split the datasets in series
        series = medipy.io.dicom.series(datasets)
        for serie in series :
            
            # Get a full dataset from file
            dataset = serie[0]
            
            if "directory_record_type" in dataset :
                dataset = medipy.io.dicom.load_dicomdir_records([dataset.children[0]])[0]
            
            # Fill patient information : name or ID
            patient = dataset.get("patients_name", dataset.get("patient_id", "(no name)")) 
            if patient not in patients :
                patients[patient] = {}
            
            # Fill study information : UID, description and date
            if dataset.study_instance_uid not in patients[patient] :
                patients[patient][dataset.study_instance_uid] = {
                    "description" : "(no description)",
                    "date" : None,
                    "series" : {}
                }
            
            study = patients[patient][dataset.study_instance_uid] 
            
            if "study_description" in dataset : 
                study["description"] = dataset.study_description
            if "study_date" in dataset and dataset.study_date :
                study["date"] = medipy.io.dicom.misc.parse_da(dataset.study_date) 
                if "study_time" in dataset and dataset.study_time :
                    time = medipy.io.dicom.misc.parse_tm(dataset.study_time)
                    study["date"] = datetime.datetime.combine(study["date"], time)
            
            # Fill series information : UID, description and date
            study["series"][dataset.series_instance_uid] = {
                "description" : "(no description)",
                "date" : None,
                "datasets" : serie
            }
            
            series = study["series"][dataset.series_instance_uid] 
            
            if "series_description" in dataset : 
                series["description"] = dataset.series_description
            if "series_date" in dataset and dataset.series_date :
                series["date"] = medipy.io.dicom.misc.parse_da(dataset.series_date)
                if "series_time" in dataset and dataset.series_time :
                    time = medipy.io.dicom.misc.parse_tm(dataset.series_time)
                    series["date"] = datetime.datetime.combine(series["date"], time)
        
        return patients
    
    @staticmethod
    def _study_comparator(s1, s2):
        """ Compare two studies by date, as computed in set_datasets.
        """
        
        if (s1[1]["date"], s2[1]["date"]) == [None, None] :
            return 0
        elif s1[1]["date"] is None : 
            return -1
        elif s2[1]["date"] is None :
            return +1
        else :
            td = (s1[1]["date"] - s2[1]["date"])
            return td.days*24*3600+td.seconds
    
    @staticmethod
    def _series_comparator(s1, s2):
        """ Compare two series by date, as computed in set_datasets.
        """
        
        if (s1[1]["date"], s2[1]["date"]) == [None, None] :
            return 0
        elif s1[1]["date"] is None : 
            return -1
        elif s2[1]["date"] is None :
            return +1
        else :
            td = (s1[1]["date"] - s2[1]["date"])
            return td.days*24*3600+td.seconds
