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
    
    def orientation_comparator(o1, o2):
        """ Compare two values of Image Orientation Patient.
        """
        
        if (o1, o2) == ((), ()) :
            return True
        elif () in (o1, o2) :
            return False
        else :
            return (numpy.linalg.norm(numpy.subtract(o1,o2), numpy.inf) <= 0.05)
    
    getters = [
        # Image Orientation (Patient) (0020,0037)
        lambda x:tuple(x.get("image_orientation_patient", ())),
        # Echo Number(s) (0018,0086)
        lambda x:x.get("echo_numbers", None),
        # Acquisition Number (0020,0012)
        lambda x:x.get("acquisition_number", None),
        # Frame Type (0008,9007)
        lambda x:tuple(x.get("mr_image_frame_type_sequence", [{}])[0].get("frame_type", ())),
        # Temporal Position Index (0020,9128)
        lambda x:x.get("frame_content_sequence", [{}])[0].get("temporal_position_index", None),
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
