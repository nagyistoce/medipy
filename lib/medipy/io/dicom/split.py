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

def stacks_dictionary(datasets):
    """ Return a dictionary of stacks from the normalized datasets. The 
        dictionary is keyed by the pairs (attribute_name, attribute_value)
        that differentiate the stacks.
        
        Refer to the stacks function for a list of the attributes that are
        used to form the stacks.
    """
    
    def orientation_comparator(o1, o2):
        """ Compare two values of Image Orientation Patient.
        """
        
        if (o1, o2) == ((), ()) :
            return True
        elif () in (o1, o2) :
            return False
        else :
            return (numpy.linalg.norm(numpy.subtract(o1,o2), numpy.inf) <= 0.05)
    
    getters = {
        # Image Orientation (Patient) (0020,0037)
        Tag(0x0020,0x0037) : lambda x:tuple(x.get("image_orientation_patient", ())),
        # Echo Number(s) (0018,0086)
        Tag(0x0018,0x0086) : lambda x:x.get("echo_numbers", None),
        # Acquisition Number (0020,0012)
        Tag(0x0020,0x0012) : lambda x:x.get("acquisition_number", None) if x.get("modality", None) != "CT" else None,
        # Frame Type (0008,9007)
        Tag(0x0008,0x9007) : lambda x:tuple(
            x.get("mr_image_frame_type_sequence", [{}])[0].
                get("frame_type", ())),
        # Temporal Position Index (0020,9128)
        Tag(0x0020,0x9128) : lambda x:x.get("frame_content_sequence", [{}])[0].
            get("temporal_position_index", None),
        # Diffusion Gradient Orientation (0018,9089)
        Tag(0x0018,0x9089) : lambda x:tuple(
            x.get("mr_diffusion_sequence", [{}])[0].
                get("diffusion_gradient_direction_sequence", [{}])[0].
                    get("diffusion_gradient_orientation", ())),
        # Diffusion b-value (0018,9087)
        Tag(0x0018,0x9087) : lambda x:x.
            get("mr_diffusion_sequence", [{}])[0].
                get("diffusion_bvalue", None)
    }
    
    dictionary = {}
    for dataset in datasets :
        dataset_key = dict([(key, getter(dataset)) 
                            for key, getter in getters.items()])
       
        # If a close value of Image Orientation Patient exists, use it
        o1 = dataset_key[Tag(0x0020,0x0037)] #values[0]
        for key in dictionary.keys() :
            key = dict(key)
            o2 = key[Tag(0x0020,0x0037)]
            if orientation_comparator(o1, o2) :
                dataset_key[Tag(0x0020,0x0037)] = o2
                break
        
        dictionary.setdefault(tuple(dataset_key.items()), []).append(dataset)

    
    return dictionary

def stacks(datasets):
    """ Return a list of stacks from the normalized datasets. Stacks are 
        defined as "groups of frames that have a geometric relationship" 
        (PS 3.3-2011, C.7.6.16.2.2.4). Stacks are formed using the following 
        attributes : 
          * Image Orientation (Patient)
          * Echo Number, 
          * Acquisition Number
          * Frame Type
          * Temporal Position Index
          * Diffusion Gradient Orientation 
          * Diffusion b-Value.
    """
    
    dictionary = stacks_dictionary(datasets)
    return dictionary.values()
