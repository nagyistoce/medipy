##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import colorsys
import random

import itk

import medipy.base
import medipy.intensity
import medipy.itk

def mono_modal_statistical_change_detection(input1, input2, mask):
    """ Mono-modal statistical change detection.
    
        Implementation of the method described in "Automatic Chance Detection in
        Multi-Modal Serial MRI: Application to Multiple Sclerosis Lesion 
        Follow-up" by Bosc & al. (NeuroImage, 2(20), pp. 643-656, Oct. 2003).
    
        <gui>
            <item name="input1" type="Image" label="Input 1"/>
            <item name="input2" type="Image" label="Input 2"/>
            <item name="mask" type="Image" label="Mask"/>
            <item name="result" type="Image" initializer="output=True"
                  role="return" label="Result"/>
        </gui>
    """
    
    itk_input_1 = medipy.itk.medipy_image_to_itk_image(input1, False)
    itk_input_2 = medipy.itk.medipy_image_to_itk_image(input2, False)
    itk_mask = medipy.itk.medipy_image_to_itk_image(mask, False)
    change_detection_filter = itk.MonoModalStatisticalChangeDetectionImageFilter[
        itk_input_1, itk_input_2, itk_input_1].New(
            Input1=itk_input_1, Input2=itk_input_2, Mask=itk_mask)
    itk_output = change_detection_filter()[0]
    return medipy.itk.itk_image_to_medipy_image(itk_output, None, True)

def change_detection_clustering(input, mask, 
                                maximum_number_of_clusters, minimum_cluster_size):
    """ Create clusters on the change detection image.
        
        <gui>
            <item name="input" type="Image" label="Input"/>
            <item name="mask" type="Image" label="Mask"/>
            <item name="maximum_number_of_clusters" type="Int" label="Max # of clusters" 
                  initializer="30"/>
            <item name="minimum_cluster_size" type="Int" label="Min cluster size" 
                  initializer="5"/>
            <item name="result" type="Image" initializer="output=True"
                  role="return" label="Result"/>
        </gui>
    """
    
    itk_input = medipy.itk.medipy_image_to_itk_image(input, False)
    itk_mask = medipy.itk.medipy_image_to_itk_image(mask, False)
    clustering_filter = itk.ChangeDetectionClusteringImageFilter[
        itk_input, itk_mask, itk_input].New(
            Input=itk_input, Mask=itk_mask, 
            MaximumNumberOfClusters=maximum_number_of_clusters,
            MinimumClusterSize=minimum_cluster_size)
    itk_output = clustering_filter()[0]
    return medipy.itk.itk_image_to_medipy_image(itk_output, None, True)

def clusters_to_annotations(image):
    """
        <gui>
            <item name="image" type="Image" label="Image" role="output"/>
        </gui>
    """
    
    itk_image = medipy.itk.medipy_image_to_itk_image(image, False)
    annotations_calculator = itk.ClustersToAnnotationsCalculator[itk_image].New(
        Image=itk_image)
    annotations_calculator.Compute()
    
    image.annotations[:] = []
    
    for label in annotations_calculator.GetAnnotationsLabels() :
        position = [x for x in reversed(annotations_calculator.GetAnnotation(label))]
        color = colorsys.hsv_to_rgb(random.random(), 1, 1)
        
        annotation = medipy.base.ImageAnnotation(
            position, str(label), 
            medipy.base.ImageAnnotation.Shape.cross, 10, color)
        image.annotations.append(annotation)
    