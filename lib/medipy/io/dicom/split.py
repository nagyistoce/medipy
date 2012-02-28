##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import copy
import logging

import numpy

from medipy.io.dicom import DataSet, Tag

def images(datasets):
    """ Return all Data Sets that are images.
    """
    
    # TODO : use Media Storage SOP Class UID instead ?
    return [d for d in datasets if isinstance(d, DataSet) and "pixel_data" in d]

def non_images(datasets):
    """ Return all Data Sets that are not images.
    """
    
    # TODO : use Media Storage SOP Class UID instead ?
    return [d for d in datasets if isinstance(d, DataSet) and "pixel_data" not in d]

def series(datasets, keep_only_images = True):
    """ Split a sequence of datasets to Series, return a list of list of Data Sets.
        If keep_only_images, then only the datasets that contain images will be
        kept, otherwise, all datasets are kept.
    """
    
    result = {}
    
    if keep_only_images :
        filtered_datasets = images(datasets)
        for dataset in datasets :
            if isinstance(dataset, DataSet) and "directory_record_sequence" in dataset :
                filtered_datasets.append(dataset)
    else :
        filtered_datasets = datasets
    
    for dataset in filtered_datasets :
        if "directory_record_sequence" in dataset :
            for record in dataset.directory_record_sequence :
                if record.directory_record_type != "SERIES" :
                    continue
                if keep_only_images :
                    modified_record = copy.deepcopy(record)
                    modified_record.children = [x for x in record.children if x.directory_record_type == "IMAGE"]
                    
                    if modified_record.children :
                        if modified_record.series_instance_uid not in result :
                            result[modified_record.series_instance_uid] = []
                        
                        result[modified_record.series_instance_uid].append(modified_record)
                else :
                    if record.series_instance_uid not in result :
                        result[record.series_instance_uid] = []
                    result[record.series_instance_uid].append(record)
        else :
            if dataset.series_instance_uid not in result :
                result[dataset.series_instance_uid] = []
    
            result[dataset.series_instance_uid].append(dataset)
    
    return result.values()

def stacks(datasets, use_dimension_index_sequence=True):
    """ Return a list of stacks from the datasets. Stacks are defined as  
        "groups of frames that have a geometric relationship" (PS 3.3-2011,
        C.7.6.16.2.2.4). Stacks are formed using the following attributes : 
          * Multi-frame images : attributes in the Dimension Index Sequence.
            Private values are skipped since we don't know their semantics.
          * NM images : attributes in the Frame Increment Pointer.
          * Other images : Image Orientation (Patient), Echo Number, 
            Acquisition Number, Diffusion Gradient Orientation and 
            Diffusion b-Value.
        Set use_dimension_index_sequence to False to bypass the processing of
        multi-frame and NM images.
    """
    
    def orientation_comparator(o1, o2):
        """ Compare two values of Image Orientation Patient.
        """
        
        if (o1, o2) == (None, None) :
            return True
        elif None in (o1, o2) :
            return False
        else :
            return (numpy.linalg.norm(numpy.subtract(o1,o2), numpy.inf) <= 0.05)
    
    if use_dimension_index_sequence and "dimension_index_sequence" in datasets[0] :
        # A dimension is a set of attributes [...] which are especially 
        # intended for presentation. (PS 3.3-2011, 7.5.2).
        # Build a list of dimensions
        dimensions = []
        for item in datasets[0].dimension_index_sequence :
            value = {
                "dimension_index_pointer" : Tag(item.dimension_index_pointer)
            }
            if "functional_group_pointer" in item :
                value["functional_group_pointer"] = Tag(item.functional_group_pointer)
            dimensions.append(value)
        
        # Build a dictionary of stacks, where the key is the 
        # Dimension Index Values (0020,9157) minus its last value which is an
        # index in the stack
        stacks_dictionary = {}
        
        for dataset in datasets :
            index_values = dataset.frame_content_sequence[0].dimension_index_values
            index = []
            for dimension in dimensions :
                if "functional_group_pointer" in dimension :
                    if dimension["functional_group_pointer"].private :
                        # Skip private Functional Groups, since their semantics
                        # are unknown
                        continue
                    else :
                        # Look for the index value in a Functional Group
                        getter = dataset[dimension["functional_group_pointer"]][0].__getitem__
                else :
                    # Look for the index value in the dataset itself
                    getter = dataset.__getitem__
                
                if dimension["dimension_index_pointer"].private :
                    # Skip private indices, since their semantics are unknown
                    continue
                else :
                    index.append(getter(dimension["dimension_index_pointer"]))
            
            index = index[:-1]
            stacks_dictionary.setdefault(tuple(index), []).append(dataset)
        
        # Process stacks that are not related to Dimension Index sequence
        result = []
        for stack in stacks_dictionary.values() :
            result.extend(stacks(stack, False))
        return result
    elif use_dimension_index_sequence and "frame_increment_pointer" in datasets[0] :
        # Nuclear Medicine image. Almost like a multi-frame image, but
        # "dimensions" are stored in Frame Increment Pointer 
        if isinstance(datasets[0].frame_increment_pointer, list) :
            dimensions = datasets[0].frame_increment_pointer
        else :
            dimensions = [datasets[0].frame_increment_pointer]
        dimensions = [Tag(x) for x in dimensions]
        
        stacks_dictionary = {}
        for dataset in datasets :
            index = []
            for tag in dimensions :
                index.append(dataset[tag][dataset.instance_number])
            
            index = index[:-1]
            stacks_dictionary.setdefault(tuple(index), []).append(dataset)
        
        # Process stacks that are not related to Dimension Index sequence
        result = []
        for stack in stacks_dictionary.values() :
            result.extend(stacks(stack, False))
        return result
    else :
        # Regular, single-frame image.
        getters = [
            # Image Orientation (Patient) (0020,0037)
            lambda x:tuple(x.get("image_orientation_patient", None)),
            # Echo Number(s) (0018,0086)
            lambda x:x.get("echo_numbers", None),
            # Acquisition Number (0020,0012)
            lambda x:x.get("acquisition_number", None),
            # Diffusion Gradient Orientation (0018,9089)
            lambda x:x.get("mr_diffusion_sequence", [{}])[0].get("diffusion_gradient_orientation", None),
            # Diffusion b-value (0018,9087)
            lambda x:x.get("mr_diffusion_sequence", [{}])[0].get("diffusion_b_value", None)
        ]
        
        stacks_dictionary = {}
        for dataset in datasets :
            values = [x(dataset) for x in getters]
            
            # If a close value of Image Orientation Patient exists, use it
            o1 = values[0]
            for key in stacks_dictionary :
                o2 = key[0]
                if orientation_comparator(o1, o2) :
                    values[0] = o2
                    break
            
            stacks_dictionary.setdefault(tuple(values), []).append(dataset)
        
        return stacks_dictionary.values()
