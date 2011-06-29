##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import logging
import operator

import numpy

def sort(datasets, sort_function = None):
    """ Sort the Data Sets (or pairs of (Data Set, frame number) in geometrical
        order.    
        
        All Data Sets must belong to the same stack.
    """
    
    if sort_function is not None :
        datasets.sort(sort_function)
        return True
    else :
        volumes = {}
        if _is_dwi(datasets) :
            # Find all values of Diffusion b-value (0018,9087) and
            # Diffusion Gradient Orientation (0018,9089)
            for dataset in datasets :
                # B-value and Gradient Orientation
                if isinstance(dataset, tuple) :
                    functional_group = dataset[0].perframe_functional_groups_sequence[dataset[1]]
                    mr_diffusion = functional_group.mr_diffusion_sequence[0]
                    key = (mr_diffusion.diffusion_bvalue, 
                           mr_diffusion.diffusion_gradient_direction_sequence[0].diffusion_gradient_orientation)
                else :
                    key = (dataset.diffusion_bvalue, tuple(dataset.diffusion_gradient_orientation))
                if key not in volumes :
                    volumes[key] = []
                volumes[key].append(dataset)
        else :
            # TODO : mosaic, fMRI
            volumes["all"] = datasets
        
        all_sorted = True
        for volume in volumes.values() :
            if not _sort_by_image_position_patient(volume) :
                logging.warning("Could not sort images using Image Position (Patient)")
                all_sorted = False
                if not _sort_by_image_number(volume) :
                    logging.warning("Could not sort images using Image Number")
                    # Sort by file name ?
                    all_sorted = False
        
        datasets[:] = reduce(operator.concat, volumes.values(), [])
        return all_sorted

def _is_dwi(datasets):
    """ Test if the Data Sets are a DWI series.
        All Data Sets must belong to the same stack.
    """
    
    for dataset in datasets :
        if (0x0018,0x9089) in dataset :
            # Diffusion Gradient Orientation
            return True
        if (0x0018,0x9117) in dataset :
            # MR Diffusion Sequence
            return True
        
    return False

def _sort_by_image_position_patient(datasets) :
    
    for dataset in datasets :
        if isinstance(dataset, tuple) :
            functional_group = dataset[0].perframe_functional_groups_sequence[dataset[1]]
            if "plane_position_sequence" not in functional_group :
                logging.warning("Plane Position Sequence is not in Per-Frame Functional Groups")
                return False
            if "image_position_patient" not in functional_group.plane_position_sequence[0] :
                logging.warning("Image Position (Patient) is not in Plane Position Sequence")
                return False
        else :
            if "image_position_patient" not in dataset :
                logging.warning("Image Position (Patient) is not in Data Set")
                return False
    
    # Compute normal from image_orientation_patient
    if isinstance(datasets[0], tuple) :
        functional_group = datasets[0][0].perframe_functional_groups_sequence[datasets[0][1]]
        v1 = functional_group.plane_orientation_sequence[0].image_orientation_patient[:3]
        v2 = functional_group.plane_orientation_sequence[0].image_orientation_patient[3:]
    else :
        v1 = datasets[0].image_orientation_patient[:3]
        v2 = datasets[0].image_orientation_patient[3:]
    normal = numpy.cross(v1, v2)
    
    def distance(dataset):
        if isinstance(dataset, tuple) :
            functional_group = dataset[0].perframe_functional_groups_sequence[dataset[1]]
            position = functional_group.plane_position_sequence[0].image_position_patient
        else :
            position = dataset.image_position_patient
        return numpy.dot(position, normal)
    
    distances = [distance(d) for d in datasets]
    distances.sort()
    if(distances[0] == distances[-1]) :
        logging.warning("All Data Sets have the same position. "
                        "Cannot sort using position")
        return False
    
    for index, d in enumerate(distances[:-1]) :
        if d == distances[index+1] :
            logging.warning("Two Data Sets have the same position. "
                            "Cannot sort using position")
            return False
    
    # Sort images using image_position_patient (0x0020,0x0032) or (0x0020, 0x0030)
    def sort_function(i1, i2):
        d1 = distance(i1)
        d2 = distance(i2)
        return int(numpy.sign(d1-d2))
    datasets.sort(sort_function)
    
    return True

def _sort_by_image_number(datasets) :
    for dataset in datasets :
        if "image_number" not in dataset :
            return False
    
    datasets.sort(lambda x,y : x.image_number-y.image_number)
    return True
