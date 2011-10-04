/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

%module parse
%{
#include "parse.h"
%}

%include "std_string.i"

%include "parse.h"

%pythoncode
%{
from dataset import DataSet

def _get_dicomdir_records_hierarchy(dataset):
    """ Return the hierarchy of records in a DICOMDIR.
    """
    
    offsets = set()
    for record in dataset.directory_record_sequence :
        offsets.add(record.offset_of_the_next_directory_record)
        offsets.add(record.offset_of_referenced_lowerlevel_directory_entity)
    offsets.add(0)
    if len(offsets) != len(dataset.directory_record_sequence) :
        raise Exception("Some records are not referenced")
    offsets = list(offsets)
    offsets.sort()
    
    for record in dataset.directory_record_sequence :
        children = []
        child_offset = record.offset_of_referenced_lowerlevel_directory_entity
        while child_offset != 0 :
            child = dataset.directory_record_sequence[offsets.index(child_offset)]
            children.append(child)
            child_offset = child.offset_of_the_next_directory_record
        record.children = children

def parse(filename) :
    dictionary = parse_file(str(filename))
    if dictionary is not None :
        dataset = DataSet.from_dict(dictionary)
        if "directory_record_sequence" in dataset :
            _get_dicomdir_records_hierarchy(dataset)
            for record in dataset.directory_record_sequence :
                record.path = filename
        return dataset
    else :
        return None
%}