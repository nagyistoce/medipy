##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import copy

import numpy

from medipy.io.dicom import DataSet

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
    """ Split a sequence of images to stacks. Return a list of datasets or
        pairs of (Data Set, frame number) if a dataset is a multi-frame
        
        All datasets must belong to the same Series.
    """
    
    def orientation_comparator(o1, o2):
        if (o1, o2) == (None, None) :
            return True
        elif None in (o1, o2) :
            return False
        else :
            return (numpy.linalg.norm(numpy.subtract(o1,o2), numpy.inf) <= 0.05)
    
    result = {}
    
    for dataset in datasets :
        if "number_of_frames" in dataset :
            # Number Of Frames : we have a Multi-frame image
            per_frame_functional_groups = dataset.perframe_functional_groups_sequence
            for frame_number, functional_group in enumerate(per_frame_functional_groups) :
                # Plane Orientation Sequence
                plane_orientation_sequence = functional_group.get("plane_orientation_sequence", None)
                if plane_orientation_sequence :
                    orientation = plane_orientation_sequence[0].get("image_orientation_patient", None)
                else :
                    orientation = None
                
                if orientation is not None :
                    orientation = tuple(orientation)
                    
                equal_values = [x for x in result.keys() if orientation_comparator(x, orientation)]
                if equal_values :
                    result[equal_values[0]].append((dataset, frame_number))
                else :
                    result[orientation] = [(dataset, frame_number)]
        else :
            # Single frame
            orientation = dataset.get("image_orientation_patient", None)
            
            if orientation is not None :
                orientation = tuple(orientation)
            
            equal_values = [x for x in result.keys() if orientation_comparator(x, orientation)]
            if equal_values :
                result[equal_values[0]].append(dataset)
            else :
                result[orientation] = [dataset]
    
    return result.values()

def acquisitions(datasets):
    """ Return a list of "acquisitions" present in the datasets. Acquisitions
        are defined as sets of datasets that have a similar :
          * Echo Number(s) (0018,0086)
          * Acquisition Number (0020,0012)
          * Diffusion Gradient Orientation (0018,9089)
          * Diffusion b-value (0018,9087)
          
        All datasets
        must belong to the same series.
        
    """
    
    result = {}
    for dataset in datasets :
        if "MOSAIC" not in dataset.image_type and "echo_numbers" in dataset :
            identifier = dataset.echo_numbers
        elif "number_of_frames" in dataset :
            # TODO
            # TODO : test for mosaic ?
            # Multi-frame : Diffusion Gradient Orientation (0018,9089) and 
            # Diffusion b-value (0018,9087) are in MR Diffusion Sequence (0018,9117)
            if "mr_diffusion_sequence" in dataset :
                pass
        elif "acquisition_number" in dataset :
            identifier = dataset.acquisition_number
        else : 
            identifier = None
        
        result.setdefault(identifier, []).append(dataset)
    
    return [result[x] for x in sorted(result.keys())]
