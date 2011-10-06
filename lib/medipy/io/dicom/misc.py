##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import datetime
import os.path

from medipy.io.dicom import parse 

def find_dicomdir_file(base_directory, referenced_file_id):
    """ Search the filesystem for upper-case and lower case version of a
        Referenced File ID in a DICOMDIR
    """
    
    filenames = [
        os.path.join(base_directory, *referenced_file_id),
        os.path.join(base_directory, *[x.lower() for x in referenced_file_id]),
        os.path.join(base_directory, *[x.upper() for x in referenced_file_id]),
    ]
    try :
        filename = [x for x in filenames if os.path.isfile(x)][0]
    except IndexError :
        raise Exception("Missing file : %s"%(os.path.join(base_directory, *referenced_file_id),))
    
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
                child_dataset = parse(filename)
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

def load_dicomdir_records(datasets, base_directory=None):
    """ If a Data Set is a DICOMDIR Record, replace it by the file it 
        (or its children) references.
    """
    
    result = []
    
    file_ids = set()
    
    for dataset in datasets :
        if "directory_record_type" in dataset : # Directory Record Type
            # Find all children with Referenced File ID (0x0004,0x1500)
            queue = [dataset]
            while queue :
                d = queue.pop()
                if "referenced_file_id" in d :
                    file_ids.add((d.path, tuple(d.referenced_file_id)))
                queue.extend(d.children)
        else :
            result.append(dataset)

    for path, file_id in file_ids :
        filename = find_dicomdir_file(os.path.dirname(path), file_id)
        result.append(parse(filename))
    
    return result

def parse_da(value):
    """ Return a datetime object from a DA value.
    """
    
    if value == "" :
        return None
    
    return datetime.datetime.strptime(value, "%Y%m%d")
    
def parse_tm(value):
    """ Return a time object from a TM value.
    """
    
    if value == "" :
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