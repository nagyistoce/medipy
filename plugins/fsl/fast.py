##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import os
import tempfile

import medipy.io

from fsl_tool import FSLTool

class FAST(FSLTool):
    """ Wrapper for FAST (FMRI's Automated Segmentation Tool) from FSL.
    
        Inputs : 
          * input
          * nb_tissues : number of tissue-type classes; default=3.
          * nb_iter_bias : number of main-loop iteraions during bias-field removal; default=4.
          * bias_smoothing : bias field smoothing extent (FWHM) in mm; default=20.
          * image_type : type of images (n=1 for T1, n=2 for T2, n=3 for PD).
          * spatial_smoothing : initial segmentation spatial smoothness (during bias field estimation); default=0.02.
          * create_separate_segment : outputs a separate binary image for each tissue type.
          * turn_off_pve : turn off PVE (partial volume estimation).
          * nb_channels : number of input images (channels); default 1.
          * nb_iter_init_segmentation : number of segmentation-initialisation iterations; default=15.
          * spatial_smooth_mixeltype : spatial smoothness for mixeltype; default=0.3.
          * nb_iter_after_bias : number of main-loop iterations after bias-field removal; default=4.
          * spatial_smooth_segmentation : segmentation spatial smoothness; default=0.1.
          * create_estimated_bias : flag to generate an estimated bias field.
          * create_restored_image : flag to generate the restored image (bias-corrected image).
    
        Outputs :
          * outputs.
    """
    
    def __init__(self, input=None, nb_tissues=3, nb_iter_bias=4, bias_smoothing=20, image_type=None, spatial_smoothing=0.02, create_separate_segment=False,
                turn_off_pve=False, nb_channels=1, nb_iter_init_segmentation=15, spatial_smooth_mixeltype=0.3,
                nb_iter_after_bias=4, spatial_smooth_segmentation=0.1, create_estimated_bias=False, create_restored_image=False, *args, **kwargs):
        
        super(FAST, self).__init__(*args, **kwargs)
        
        self.input = input
        self.nb_tissues = nb_tissues
        self.nb_iter_bias = nb_iter_bias
        self.bias_smoothing = bias_smoothing
        self.image_type = image_type
        self.spatial_smoothing = spatial_smoothing
        self.create_separate_segment = create_separate_segment
        self.turn_off_pve = turn_off_pve
        
        self.nb_channels = nb_channels
        self.nb_iter_init_segmentation = nb_iter_init_segmentation
        self.spatial_smooth_mixeltype = spatial_smooth_mixeltype
        self.nb_iter_after_bias = nb_iter_after_bias
        self.spatial_smooth_segmentation = spatial_smooth_segmentation
        
        self.create_estimated_bias = create_estimated_bias
        self.create_restored_image = create_restored_image

    def get_partial_volume_map(self, class_):
        pass

    ##############################
    # Property-related functions #
    ##############################
    
    def _get_separate_segment(self):
        if self.create_separate_segment :
            return os.path.extsep.join([os.path.splitext(self.input)[0]+"_separate_segment", self._suffixes[self.fsl_environment["FSLOUTPUTTYPE"]]])
        else :
            return None
            
    def _get_estimated_bias(self):
        if self.create_estimated_bias :
            return os.path.extsep.join([os.path.splitext(self.input)[0]+"_estimated_bias", 
                                        self._suffixes[self.fsl_environment["FSLOUTPUTTYPE"]]])
        else :
            return None
    
    def _get_restored_image(self):
        if self.create_restored_image :
            return os.path.extsep.join([os.path.splitext(self.input)[0]+"_restored_image", 
                                        self._suffixes[self.fsl_environment["FSLOUTPUTTYPE"]]])
        else :
            return None
            
    def _get_command(self):
        if None in [self.input, self.image_type] :
            raise Exception("Input and image type must be specified. Usage : fast [options] file(s)")
            
        command = ["fast"]
        
        if self.nb_tissues :
            command.extend(["-n", "{0}".format(self.nb_tissues)])
        if self.nb_iter_bias :
            command.extend(["-I", "{0}".format(self.nb_iter_bias)])
        if self.bias_smoothing :
            command.extend(["-l", "{0}".format(self.bias_smoothing)])
        if self.image_type :
            command.extend(["-t", "{0}".format(self.image_type)])
        if self.spatial_smoothing :
            command.extend(["-f", "{0}".format(self.spatial_smoothing)])
        if self.create_separate_segment :
            command.append("-g")
        if self.turn_off_pve :
            command.append("--nopve")
        if self.create_estimated_bias :
            command.append("-B")
        if self.create_restored_image :
            command.append("-b")
        
        if self.nb_channels :
            command.extend(["-S", "{0}".format(self.nb_channels)])
        if self.nb_iter_init_segmentation :
            command.extend(["-W", "{0}".format(self.nb_iter_init_segmentation)])
        if self.spatial_smooth_mixeltype :
            command.extend(["-R", "{0}".format(self.spatial_smooth_mixeltype)])
        if self.nb_iter_after_bias :
            command.extend(["-O", "{0}".format(self.nb_iter_after_bias)])
        if self.spatial_smooth_segmentation :
            command.extend(["-H", "{0}".format(self.spatial_smooth_segmentation)])
        
        if self.input :
            command.extend([self.input])
        
        return command

    def _get_binary_segmentation(self):
        pass
    
    separate_segment = property(_get_separate_segment)
    restored_image = property(_get_restored_image)
    estimated_bias = property(_get_estimated_bias)
    binary_segmentation = property(_get_binary_segmentation)


