/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

%module dataset_io
%{
#include "dataset_io.h"
%}

%include "std_string.i"

%include "dataset_io.h"

%pythoncode
%{
import medipy.base
from dataset import DataSet

def _get_dicomdir_records_hierarchy(dataset):
    """ Return the hierarchy of records in a DICOMDIR.
    """
    
    offsets = set()
    for record in dataset.directory_record_sequence.value :
        offsets.add(record.offset_of_the_next_directory_record.value)
        offsets.add(record.offset_of_referenced_lowerlevel_directory_entity.value)
    offsets.add(0)
    if len(offsets) != len(dataset.directory_record_sequence.value) :
        raise medipy.base.Exception("Some records are not referenced")
    offsets = list(offsets)
    offsets.sort()
    
    for record in dataset.directory_record_sequence.value :
        children = []
        child_offset = record.offset_of_referenced_lowerlevel_directory_entity.value
        while child_offset != 0 :
            child = dataset.directory_record_sequence.value[offsets.index(child_offset)]
            children.append(child)
            child_offset = child.offset_of_the_next_directory_record.value
        record.children = children

unwrapped_read = read

def read(filename) :
    try :
        dataset = unwrapped_read(str(filename))
    except Exception, e :
        raise medipy.base.Exception(e)
    if dataset is not None :
        if "directory_record_sequence" in dataset :
            _get_dicomdir_records_hierarchy(dataset)
            for record in dataset.directory_record_sequence.value :
                record.path = filename
        return dataset
    else :
        return None
%}
