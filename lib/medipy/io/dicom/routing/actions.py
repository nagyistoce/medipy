##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import datetime
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
    """ Set an element in the dataset. If the element is not present, it will 
        be inserted : ::
    
            >>> dataset = medipy.io.dicom.DataSet()
            >>> dataset.patients_name = "Brouchard^Georges"
            >>> set_patients_name = medipy.io.dicom.routing.SetElement("patients_name", "Doe^John")
            >>> set_patients_name(dataset)
            >>> dataset.patients_name.value
            'Doe^John'
        
        This overwrites the element. If you want to keep track of modifications,
        use :class:`ModifyDataSet`.
    """
    
    def __init__(self, tag, value):
        self.tag = tag
        self.value = value
    
    def __call__(self, dataset):
        dataset[self.tag] = self.value

class DeleteElement(Action) :
    """ Remove an element from a dataset. If the element is not in the dataset, 
        nothing happens.::
        
            >>> dataset = medipy.io.dicom.DataSet()
            >>> dataset.patients_name = "Brouchard^Georges"
            >>> delete_patients_name = medipy.io.dicom.routing.DeleteElement("patients_name")
            >>> delete_patients_name(dataset)
            >>> "patients_name" in dataset
            False
    """
    
    def __init__(self,tag):
        self.tag = tag
        
    def __call__(self, dataset):
        if self.tag in dataset :
            del dataset[self.tag]

class EmptyElement(Action) :
    """ Set an element in the dataset to an empty value. If the element is not 
        present in the dataset, it is inserted. The VR of the tag MUST NOT be 
        one of AT, FL, FD, SL, SS, UL, US. ::
        
            >>> dataset = medipy.io.dicom.DataSet()
            >>> dataset.patients_name = "Brouchard^Georges"
            >>> clear_patients_name = medipy.io.dicom.routing.EmptyElement("patients_name")
            >>> clear_patients_name(dataset)
            >>> "patients_name" in dataset
            True
            >>> dataset.patients_name.value
            ''
    """
    
    def __init__(self, tag, vr=None):
        
        if isinstance(tag, basestring):
            tag = medipy.io.dicom.dictionary.name_dictionary.get(tag, None)
            if tag is not None :
                tag = medipy.io.dicom.Tag(tag)
        else :
            tag = medipy.io.dicom.Tag(tag)
        
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

class ModifyDataSet(Action):
    """ Modify a dataset, storing the modifications in Original Attributes 
        Sequence (0400, 0561).
    
        * ``modifications`` can be either a dictionary mapping tags (numerical
          or named) to their modified values, or a callable object. In the 
          latter case, the object is called with the dataset as its only 
          argument, and must return a dictionary mapping tags to their modified
          values. Original Attributes Sequence (0400, 0561) MUST NOT be present
          in the dictionary
        * ``source`` : the source that provided the SOP Instance prior to the 
          removal or replacement of the values. For example, this might be the 
          Institution from which imported SOP Instances were received.
        * ``system`` : identification of the system which removed and/or 
          replaced the attributes.
        * ``reason`` : reason for the attribute modification, must be one of the 
          following :
        
          * ``"COERCE"`` : replace values of attributes such as Patient Name, 
            ID,  Accession Number, for example, during import of media from an 
            external institution, or reconciliation against a master patient 
            index.
          * ``"CORRECT"`` : replace incorrect values, such as Patient Name or 
            ID, for example, when incorrect worklist item was chosen or 
            operator input error.
        * ``modification_datetime`` : date and time the attributes were removed 
          and/or replaced. If not specified, ``datetime.datetime.now()`` is used.
        
        Example with a dictionary : ::
        
            >>> dataset = medipy.io.dicom.DataSet()
            >>> dataset.patients_name = "Doe^John"
            >>> dataset.series_description = "3D T1 SPGR 180 slices 1mm"
            >>> modify = medipy.io.dicom.routing.ModifyDataSet({"series_description": "T1_3D"}, "source", "system", "COERCE")
            >>> modify(dataset)
            >>> dataset.series_description.value
            'T1_3D'
        
        Example with a callable : ::
        
            >>> def modify_series_description(dataset) :
            ...     series_description = dataset.get(
            ...         "series_description", medipy.io.dicom.LO(None)).value
            ...     if series_description == "3D T1 SPGR 180 slices 1mm" :
            ...         return { "series_description" : "T1_3D" }
            ...     else :
            ...         return {}
            >>> dataset = medipy.io.dicom.DataSet()
            >>> dataset.patients_name = "Doe^John"
            >>> dataset.series_description = "3D T1 SPGR 180 slices 1mm"
            >>> modify = medipy.io.dicom.routing.ModifyDataSet(modify_series_description, "source", "system", "COERCE")
            >>> modify(dataset)
            >>> dataset.series_description.value
            'T1_3D'
    """
    
    def __init__(self, modifications, source, system, reason, modification_datetime=None):
        self.modifications = modifications
        self.source = source
        self.system = system
        self.reason = reason
        self.modification_datetime = modification_datetime
    
    def __call__(self, dataset):
        # Check that reason is a defined value
        if self.reason not in ["COERCE", "CORRECT"] :
            raise medipy.base.Exception(
                "Reason must be one of \"COERCE\" or \"CORRECT\", "
                "{0!r} found".format(self.reason))
        
        # Get the modifications dictionary if we got a function
        if callable(self.modifications) :
            modifications = self.modifications(dataset)
        else :
            modifications = self.modifications
        
        if not modifications :
            return
        
        # Check that Original Attributes Sequence (0400, 0561) is not in the
        # modifications dictionary
        original_attributes_sequence = [
            medipy.io.dicom.Tag(0x0400, 0x0561), (0x0400,0x0561), 0x04000561,
            "original_attributes_sequence"]
        if set(modifications.keys()).intersection(original_attributes_sequence) :
            raise medipy.base.Exception("Original Attributes Sequence (0400, 0561) "
                                        "may not be modified")
        
        # Set the modification datetime if not specified
        if self.modification_datetime is None :
            modification_datetime = datetime.datetime.now()
        else :
            modification_datetime = self.modification_datetime
        # Format it to the DICOM DT VR
        modification_datetime = modification_datetime.strftime("%Y%m%d" "%H%M%S")
        
        # Modify the attributes, storing the old values in modified_attributes
        # if necessary.
        modified_attributes = medipy.io.dicom.DataSet()
        for tag, value in modifications.items() :
            if tag in dataset :
                modified_attributes[tag] = dataset[tag]
            dataset[tag] = value
        
        # Preserve Original Attributes Sequence (0400, 0561) if needed
        if "original_attributes_sequence" in dataset :
            modified_attributes.original_attributes_sequence = dataset.original_attributes_sequence
        
        # Create the Original Attributes Sequence (0400,0561) and store it in the
        # dataset
        original_attributes = medipy.io.dicom.DataSet()
        original_attributes.source_of_previous_values = self.source
        original_attributes.attribute_modification_datetime = modification_datetime
        original_attributes.modifying_system = self.system
        original_attributes.reason_for_the_attribute_modification = self.reason
        
        original_attributes.modified_attributes_sequence = [modified_attributes]
        dataset.original_attributes_sequence = [original_attributes]

class RestoreDataSet(Action):
    """ Restore the original attributes of a dataset that has been modified by
        :class:`ModifyDataSet`. The values of ``source``, ``system``, and 
        ``reason`` must match those passed to :class:`ModifyDataSet` for the
        restore operation to succeed. ``modification_datetime`` will be checked
        only if it is not ``None`` (its default value). :: 
        
            >>> dataset = medipy.io.dicom.DataSet()
            >>> dataset.series_description = "3D T1 SPGR 180 slices 1mm"
            >>> modify = medipy.io.dicom.routing.ModifyDataSet({"series_description": "T1_3D"}, "source", "system", "COERCE")
            >>> modify(dataset)
            >>> restore = medipy.io.dicom.routing.RestoreDataSet("source", "system", "COERCE")
            >>> restore(dataset)
            >>> dataset.series_description.value
            '3D T1 SPGR 180 slices 1mm'
            >>> "original_attributes_sequence" in dataset
            False
    """
    
    def __init__(self, source, system, reason, modification_datetime=None):
        self.source = source
        self.system = system
        self.reason = reason
        self.modification_datetime = modification_datetime
    
    def __call__(self, dataset):
        if "original_attributes_sequence" not in dataset :
            # Nothing to do
            return
        
        # Check that reason is a defined value
        if self.reason not in ["COERCE", "CORRECT"] :
            raise medipy.base.Exception(
                "Reason must be one of \"COERCE\" or \"CORRECT\", "
                "{0!r} found".format(self.reason))
        
        original_attributes = dataset.original_attributes_sequence.value[0]
        
        # Check that modifications parameters match
        source = original_attributes.source_of_previous_values.value
        if source != self.source :
            raise medipy.base.Exception(
                "Cannot restore attributes with source {0!r}: "
                "{1!r} expected".format(source, self.source))
        
        system = original_attributes.modifying_system.value
        if system != self.system :
            raise medipy.base.Exception(
                "Cannot restore attributes with system {0!r}: "
                "{1!r} expected".format(system, self.system)) 
        
        reason = original_attributes.reason_for_the_attribute_modification.value
        if reason != self.reason :
            raise medipy.base.Exception(
                "Cannot restore attributes with reason {0!r}: "
                "{1!r} expected".format(reason, self.reason))
        
        modification_datetime = original_attributes.attribute_modification_datetime.value
        if self.modification_datetime and modification_datetime != self.modification_datetime :
            raise medipy.base.Exception(
                "Cannot restore attributes with modification date time {0!r}: "
                "{1!r} expected".format(modification_datetime, self.modification_datetime))
        
        # Restore the attributes
        modified_attributes = original_attributes.modified_attributes_sequence.value[0]
        for key, value in modified_attributes.items() :
            dataset[key] = value
        
        if "original_attributes_sequence" in modified_attributes :
            dataset.original_attributes_sequence = modified_attributes.original_attributes_sequence
        else :
            del dataset["original_attributes_sequence"]

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
            filename = "{0:X}".format(zlib.crc32(patient+study+series+instance)&0xffffffff)
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
