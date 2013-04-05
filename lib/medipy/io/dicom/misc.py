##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import datetime
import os.path
import uuid

import medipy.base

import dataset_io
from vr import * 

def find_dicomdir_file(base_directory, referenced_file_id):
    """ Search the filesystem for upper-case and lower case version of a
        Referenced File ID in a DICOMDIR
    """
    
    if isinstance(referenced_file_id, basestring) :
        # Only one path component, normalize to a list with one element
        referenced_file_id = [referenced_file_id]
    
    filenames = [
        os.path.join(base_directory, *referenced_file_id),
        os.path.join(base_directory, *[x.lower() for x in referenced_file_id]),
        os.path.join(base_directory, *[x.upper() for x in referenced_file_id]),
    ]
    try :
        filename = [x for x in filenames if os.path.isfile(x)][0]
    except IndexError :
        raise medipy.base.Exception("Missing file : %s"%(os.path.join(base_directory, *referenced_file_id),))
    
    return filename

def uid_and_description(datasets, base_directory=None, 
                        default_description = "(no description)"):
    """ Return a pair (Series Instance UID, Series Description) for the Data 
        Sets.
        
        All Data Sets must belong to the same Series.
    """
    
    series_instance_uid = datasets[0].series_instance_uid
    series_description = default_description
    for dataset in datasets :
        if "series_description" in dataset :
            series_description = dataset.series_description
            break
        elif "directory_record_type" in dataset : # Directory Record Type
            found_in_child = False
            for child in dataset.children :
                filename = find_dicomdir_file(os.path.dirname(child.path),
                                              child.referenced_file_id)
                child_dataset = dataset_io.read(filename)
                if "series_description" in child_dataset :
                    series_description = child_dataset.series_description
                    found_in_child = True
                    break
            if found_in_child :
                break
    
    return series_instance_uid, series_description

def get_series_datetime(dataset):
    """ Return a datetime object containing the Series Date and Series Time, or
        None is they are not in the dataset.
    """
    
    series_datetime = None
    if "series_date" in dataset :
        series_datetime = parse_da(dataset.series_date)
        if "series_time" in dataset :
            series_time = parse_tm(dataset.series_time)
            series_datetime = datetime.datetime.combine(series_datetime, 
                                                        series_time)
    
    return series_datetime

def get_child_file_records(record):
    """ Return children of given DICOMDIR Record which contain 
        Referenced File ID (0004,1500).
    """
    
    result = []
    
    queue = [record]
    while queue :
        dataset = queue.pop()
        if "referenced_file_id" in dataset :
            if isinstance(dataset.referenced_file_id.value, basestring) :
                # Normalize to a list with one element
                dataset.referenced_file_id = [dataset.referenced_file_id.value]
            result.append(dataset)
        queue.extend(dataset.children)
    
    return result

def load_dicomdir_records(datasets):
    """ If a Data Set is a DICOMDIR Record, replace it by the file it 
        (or its children) references.
    """
    
    result = []
    
    file_ids = set()
    
    for dataset in datasets :
        if "directory_record_type" in dataset : # Directory Record Type
            children = get_child_file_records(dataset)
            file_ids.update([(child.path, tuple(child.referenced_file_id.value))
                              for child in children])
        else :
            result.append(dataset)

    for path, file_id in file_ids :
        filename = find_dicomdir_file(os.path.dirname(path), file_id)
        result.append(dataset_io.read(filename))
    
    return result

def parse_da(value):
    """ Return a datetime object from a DA value.
    """

    if isinstance(value, VR) :
        value = value.value
    
    if not value :
        return None
    
    return datetime.datetime.strptime(value, "%Y%m%d")
    
def parse_tm(value):
    """ Return a time object from a TM value.
    """
    
    if isinstance(value, VR) :
        value = value.value
    
    if not value :
        return None
    
    format = ""
    length = 0
    microseconds = 0
    
    if len(value) >= 2 :
        format += "%H"
        length += 2
    if len(value) >= 4 :
        format += "%M"
        length += 2
    if len(value) >= 6 :
        format += "%S"
        length += 2
    if len(value) >= 8 :
        microseconds = int(value[7:]+(13-len(value))*"0")
        
    time = datetime.datetime.strptime(value[:length], format)
    return datetime.time(time.hour, time.minute, time.second, microseconds)

def generate_uid(convert_to_vr=True) :
    """ Generate a DICOM Unique Identifier using the method 
        described on David Clunie's website
        http://www.dclunie.com/medical-image-faq/html/part2.html#UID
    """
    
    uid = "2.25.{0}".format(uuid.uuid4().int)
    # Make sure the generated UID is not larger than the 64 characters specified
    # by the DICOM standard
    uid = uid[:64]
    
    if convert_to_vr :
        uid = UI(uid)
    
    return uid
