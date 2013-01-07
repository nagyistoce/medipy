##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import os
import threading

import medipy.base
import medipy.io.dicom

class Action(object):
    """ Base class for all routing actions. All concrete derived classes must
        implement the __call__ method.
    """
    
    def __call__(self, dataset):
        raise NotImplementedError()

class SetElement(Action):
    """ Set an element in the dataset, inserting it if it does not exist.
    """
    
    def __init__(self, tag, value):
        self.tag = tag
        self.value = value
    
    def __call__(self, dataset):
        dataset[self.tag] = self.value

class DeleteElement(Action) :
    """ Remove an element from a dataset. If the element is not in the dataset,
        nothing happens.
    """
    
    def __init__(self,tag):
        self.tag = tag
        
    def __call__(self, dataset):
        if self.tag in dataset :
            del dataset[self.tag]

class EmptyElement(Action) :
    """ Set an element in the dataset to an empty value. If the element is not
        present in the dataset, it is inserted.
        
        The VR of the tag MUST NOT be one of AT, FL, FD, SL, SS, UL, US
    """
    
    def __init__(self, tag, vr=None):
        
        if isinstance(tag, basestring):
            tag = medipy.io.dicom.dictionary.name_dictionary.get(tag, None)
            if tag is not None :
                tag = medipy.io.dicom.Tag(tag)
        else :
            tag = Tag(tag)
        
        self.tag = tag
        if vr is None :
            if tag not in medipy.io.dicom.dictionary.data_dictionary :
                raise medipy.base.Exception("Cannot guess VR of {0!r}: "
                                            "not in dictionary".format(tag))
            vr = medipy.io.dicom.dictionary.data_dictionary[tag][0]
        self.vr = vr
    
    def __call__(self, dataset):
        if self.vr == "SQ" :
            dataset[self.tag] = []
        else :
            dataset[self.tag] = ""

class SaveDataSet(Action):
    """ Save the dataset to a file. The root directory must be empty the first
        time this rule is called.
        
        This action is thread-safe when called on different data-sets.
    
        mode can be either "flat" or "hierarchical"
    """
    
    _lock = threading.Lock()
    
    def __init__(self, root, mode="flat"):
        self.root = root
        self.mode = mode
        
        self._next_file = 0
        self._hierarchy = {}
        self._next_series_file = {}
    
    def __call__(self, dataset):
        self._lock.acquire(True)
        
        if self.mode == "flat" :
            filename = os.path.join(self.root, "{0:08}".format(self._next_file))
            self._next_file += 1
        elif self.mode == "hierarchical" :
            patient_key = dataset.get("patient_id", None).value
            study_key = dataset.get("study_instance_uid", None).value
            series_key = dataset.get("series_instance_uid", None).value
            
            patient = self._hierarchy.setdefault(patient_key, (len(self._hierarchy), {}))
            study = patient[1].setdefault(study_key, (len(patient[1]), {}))
            series = study[1].setdefault(series_key, len(study[1]))
            
            dirname = os.path.join(self.root, 
                "{0:08}".format(patient[0]), "{0:08}".format(study[0]), 
                "{0:08}".format(series))
            
            if not os.path.isdir(dirname) :
                os.makedirs(dirname)
            
            filename = os.path.join(dirname, 
                "{0:08}".format(self._next_series_file.setdefault(series, 0)))
            
            self._next_series_file[series] += 1
        
        medipy.io.dicom.write(dataset, filename) 
        
        self._lock.release()

# TODO StoreDataset(connection)
# TODO stop_rule_processing
