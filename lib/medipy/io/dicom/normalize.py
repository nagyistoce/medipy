##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

""" This module modifies datasets to facilitate their later reconstruction as
    images. Resulting datasets may not be DICOM-compliant.
"""

import numpy
import operator

import medipy.base
from dataset import DataSet

def normalize(dataset_or_datasets):
    """ Normalize a dataset or a sequence of datasets.
    """
    
    if isinstance(dataset_or_datasets, DataSet) :
        if dataset_or_datasets.normalized :
            return dataset_or_datasets
        if "perframe_functional_groups_sequence" in dataset_or_datasets :
            single_frames = multi_frame(dataset_or_datasets)
            return [single_frame(x) for x in single_frames]
        elif "frame_increment_pointer" in dataset_or_datasets :
            single_frames = nuclear_medicine(dataset_or_datasets)
            return [single_frame(x) for x in single_frames]
        else :
            return single_frame(dataset_or_datasets)
    else :
        single_frames = []
        for dataset in dataset_or_datasets :
            result = normalize(dataset)
            if isinstance(result, DataSet) :
                single_frames.append([result])
            else :
                single_frames.append(result)
        return reduce(operator.concat, single_frames, [])

def single_frame(dataset):
    """ Normalize private tags for whom a public counterpart exist. 
    """
    
    result = DataSet()
    for tag, value in dataset.items() :
        # cf. http://www.na-mic.org/Wiki/index.php/NAMIC_Wiki:DTI:DICOM_for_DWI_and_DTI
        result[tag] = value
    result.normalized = True

    return result

def multi_frame(dataset):
    """ Return a dataset for each frame (i.e. 2D pixel plane) of the input
        dataset.
    """
    
    result = []
    
    # Per-frame functional groups that will be merged in the frame dataset
    merge_groups = [
        (0x0028,0x9110), # Pixel Measures Sequence
        # (0020,9111) Frame Content Sequence : remove "Frame" from tag names and merge (?)
        (0x0020,0x9113), # Plane Position Sequence
        (0x0020,0x9116), # Plane Orientation Sequence
        #TODO (0+ items) (0008,1140) Referenced Image Sequence
        #TODO (0+ items) (0008,9124) Derivation Image Sequence
        # (0018,9118) Cardiac Synchronization Sequence : skip, attributes not present anywhere else 
        #(0020,9071) Frame Anatomy Sequence : remove "Frame" from tag names and merge (?)
        (0x0028,0x9145), # Pixel Value Transformation Sequence
        (0x0028,0x9132), # Frame VOI LUT Sequence 
        #TODO (1+ items) (0040,9096) Real World Value Mapping Sequence
        #TODO (1+ items) (0018,9341) Contrast/Bolus Usage Sequence
        #TODO (1+ items) (0028,9422) Pixel Intensity Relationship LUT Sequence
        #TODO (1+ items) (0028,9415) Frame Pixel Shift Sequence
        (0x0020,0x9450), # Patient Orientation in Frame Sequence
        #TODO (0018,9472) Frame Display Shutter Sequence
        #(0020,9253) Respiratory Synchronization Sequence : skip, attributes not present anywhere else
        #(0x0018,0x9477), # Irradiation Event Identification Sequence
        #TODO (1+ items) (0018,9737) Radiopharmaceutical Usage Sequence
        #TODO (0018,9771) Patient Physiological State Sequence
        #TODO (0020,930E) Plane Position (Volume) Sequence
        #TODO (0020,930F) Plane Orientation (Volume) Sequence
        #TODO (0020,9310) Temporal Position Sequence
        #(0018,9807) Image Data Type Sequence : skip, attributes not present anywhere else
        #
        #(0048,021A) Plane Position (Slide) Sequence
        #(0048,0207) Optical Path Identification Sequence
        #TODO (0048,0110) Specimen Reference Sequence
        #
        #(0018,9226) MR Image Frame Type Sequence : store Frame Type as Image Type (0008,0008), merge rest (?)
        #(0018,9112) MR Timing and Related Parameters Sequence
        #(0018,9125) MR FOV/Geometry Sequence
        #(0018,9114) MR Echo Sequence : store Effective Echo Time as Echo Time (0018,0081) (?)
        #(0018,9115) MR Modifier Sequence
        #(0018,9006) MR Imaging Modifier Sequence
        #(0018,9042) MR Receive Coil Sequence
        #(0018,9049) MR Transmit Coil Sequence
        #(0018,9117) MR Diffusion Sequence
        #(0018,9119) MR Averages Sequence
        #TODO (0018,9107) MR Spatial Saturation Sequence
        #(0018,9152) MR Metabolite Map Sequence
        #TODO (0018,9197) MR Velocity Encoding Sequence
        #TODO (0018,9251) MR Arterial Spin Labeling Sequence
        #
        #(0018,9103) MR Spectroscopy FOV/Geometry Sequence
        #
        #(0018,9227) MR Spectroscopy Frame Type Sequence
        #
        #(0018,9329) CT Image Frame Type Sequence
        #(0018,9301) CT Acquisition Type Sequence
        #(0018,9304) CT Acquisition Details Sequence
        #(0018,9308) CT Table Dynamics Sequence
        #(0018,9326) CT Position Sequence
        #(0018,9312) CT Geometry Sequence
        #(0018,9314) CT Reconstruction Sequence
        #(0018,9321) CT Exposure Sequence
        #(0018,9325) CT X-Ray Details Sequence
        #(0028,9145) Pixel Value Transformation Sequence
        #TODO (0018,9360) CT Additional X-Ray Source Sequence
        #
        #TODO (0022,0031) Ophthalmic Frame Location Sequence
        #
        #(0018,9412) XA/XRF Frame Characteristics Sequence
        #(0018,9432) Field of View Sequence
        #TODO (0018,9434) Exposure Control Sensing Regions Sequence
        #(0028,9443) Frame Pixel Data Properties Sequence
        #(0018,9451) Frame Detector Parameters Sequence
        #(0018,9455) Calibration Sequence
        #(0018,9456) Object Thickness Sequence
        #(0018,9417) Frame Acquisition Sequence
        #(0018,9401) Projection Pixel Calibration Sequence
        #(0018,9405) Positioner Position Sequence
        #(0018,9406) Table Position Sequence
        #(0018,9407) Collimator Shape Sequence
        #(0018,9462) Isocenter Reference System Sequence
        #(0018,9476) X-Ray Geometry Sequence
        #
        #(0062,000A) Segment Identification Sequence
        #
        #(0018,9504) X-Ray 3D Frame Type Sequence
        #
        #(0018,9751) PET Frame Type Sequence
        #(0018,9732) PET Frame Acquisition Sequence
        #(0018,9733) PET Detector Motion Details Sequence
        #(0018,9735) PET Position Sequence
        #(0018,9736) PET Frame Correction Factors Sequence
        #(0018,9749) PET Reconstruction Sequence
        #(0018,9734) PET Table Dynamics Sequence
        #
        #(0018,9806) US Image Description Sequence
        #
        #(0052,0029) Intravascular OCT Frame Content Sequence
    ]
    
    # Size in bytes of a frame
    frame_size = dataset.bits_allocated/8*dataset.rows*dataset.columns
    
    for frame_number in range(dataset.number_of_frames) :
        frame = DataSet()
        
        for key, value in dataset.items() :
            if key == (0x0028,0x0008) :
                # Number of Frames : do not store it in frame dataset
                pass
            elif key == (0x5200,0x9229) :
                # Shared Functional Groups Sequence : get the group and assign
                # elements in the group to the frame dataset
                for group_key, group_value in value[0].items() :
                    if group_key in merge_groups : 
                        # Copy the values of the only sequence item directly to
                        # the frame dataset
                        frame.update(group_value[0])
                    else :
                        # Otherwise, copy the group element to the frame dataset
                        frame[group_key] = group_value 
            elif key == (0x5200,0x9230) : 
                # Per-frame Functional Groups Sequence : get the group based
                # on the frame number, and assign elements in the group to
                # the frame dataset
                group = value[frame_number]
                for group_key, group_value in group.items() :
                    if group_key in merge_groups : 
                        # Copy the values of the only sequence item directly to
                        # the frame dataset
                        frame.update(group_value[0])
                    else :
                        # Otherwise, copy the group element to the frame dataset
                        frame[group_key] = group_value
            elif key == (0x7fe0, 0x0010) :
                # Pixel Data
                offset = frame_size*frame_number
                frame[key] = dataset.pixel_data[offset:offset+frame_size]
            else :
                # Otherwise add the element to the frame dataset as is
                frame[key] = value
        result.append(frame)
    
    return result

def nuclear_medicine(dataset):
    """ Return a dataset for each frame (i.e. 2D pixel plane) of the input
        dataset.
    """
    
    result = []
    
    # Size in bytes of a frame
    frame_size = dataset.bits_allocated/8*dataset.rows*dataset.columns
    
    # Orientation and normal of the frames, one per detector
    orientations = []
    normals = []
    for detector_information in dataset.detector_information_sequence :
        orientation = detector_information.image_orientation_patient
        orientations.append(orientation)
        
        normal = numpy.cross(orientation[:3], orientation[3:])
        normals.append(normal)
    
    # Z-spacing and vector between frames, one per detector
    slice_vectors = []
    z_spacing = dataset.get("slice_thickness", dataset.get("spacing_between_slices", 1.))
    for normal in normals :
        slice_vector = normal.copy()
        slice_vector[2] = z_spacing
        slice_vectors.append(slice_vector)
    
    for frame_number in range(dataset.number_of_frames) :
        detector = 0
        if (isinstance(dataset.frame_increment_pointer, list) and 
            (0x0054,0x0020) in dataset.frame_increment_pointer) :
            detector = dataset.detector_vector[frame_number]
        origin = dataset.detector_information_sequence[detector].image_position_patient
        
        frame = DataSet()
        frame.update(dataset)
        delattr(frame, "number_of_frames")
        
        # Instance Number : use the frame_number so we can extract information
        # from the elements referenced in Frame Increment Pointer
        frame.instance_number = frame_number
        
        # Image Orientation (Patient) : use the one from Detector Information Sequence (0054,0022)
        frame.image_orientation_patient = orientations[detector]
        
        # Image Position (Patient) : derive from dataset origin, frame normal and Z spacing
        frame.image_position_patient = origin+frame_number*slice_vectors[detector]
        
        # Pixel Data (7fe0,0010) : extract slice from 3D dataset
        offset = frame_size*frame_number
        frame.pixel_data = dataset.pixel_data[offset:offset+frame_size]
        
        result.append(frame)
    
    return result