##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import os
import threading
import zlib

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
    """ Save the dataset to a file. 
        
        The destination filename is base on Patient ID, Study Instance UID,
        Series Instance UID and SOP Instance UID. The organization of the files
        in ``root`` is according to ``mode`` :
        
        * ``"flat"``: datasets are saved directly in ``root``
        * ``"hierarchical"``: datasets are saved in sub-directories according to
          their Patient ID, Study Instance UID and Series Instance UID.
        
        This action is thread-safe.
    """
    
    _lock = threading.Lock()
    
    def __init__(self, root, mode="flat"):
        self.root = root
        self.mode = mode
        
    def __call__(self, dataset):
        
        patient = dataset.get("patient_id", medipy.io.dicom.CS("")).value
        study = dataset.get("study_instance_uid", medipy.io.dicom.UI("")).value
        series = dataset.get("series_instance_uid", medipy.io.dicom.UI("")).value
        instance = dataset.get("sop_instance_uid", medipy.io.dicom.UI("")).value
        
        if self.mode == "flat" :
            filename = "{0:X}".format(zlib.crc32(patient+study+series+instance))
        elif self.mode == "hierarchical" :
            filename = os.path.join("{0:X}".format(zlib.crc32(patient)&0xffffffff),
                                    "{0:X}".format(zlib.crc32(study)&0xffffffff),
                                    "{0:X}".format(zlib.crc32(series)&0xffffffff),
                                    "{0:X}".format(zlib.crc32(instance)&0xffffffff))
        
        destination = os.path.join(self.root, filename)
        
        destination_dir = os.path.dirname(destination)
        # Make sure we have a lock before creating the dir, otherwise we might
        # have a race condition and an exception thrown in makedirs
        self._lock.acquire(True)
        if not os.path.isdir(destination_dir) :
            os.makedirs(destination_dir)
        self._lock.release()
        
        medipy.io.dicom.write(dataset, destination)

# TODO StoreDataset(connection)
# TODO stop_rule_processing
