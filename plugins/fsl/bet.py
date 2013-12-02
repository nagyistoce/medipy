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

class BET(FSLTool):
    """ Wrapper for BET (Brain Extraction Tool) from FSL.
    
        Inputs : 
          * input
          * create_outline : flag to generate an overlay of the brain surface
            outline on the input image.
          * create_brain_mask : flag to generate the binary brain mask.
          * create_skull : flag to generate the approximate skull image.
          * intensity_threshold : fractional intensity threshold (0->1). Default
            to 0.5.
          * gradient_threshold : vertical gradient in fractional intensity
            threshold (-1->1). Default to 0.
          * head_radius : head radius (mm not voxels); initial surface sphere 
            is set to half of this. Automatically computed if not specified.
          * center_of_gravity : center of gravity (voxels not mm) of initial
            mesh surface. Automatically computed if not specified.
          * robust: robust brain center estimation.
    
        Outputs :
          * output
          * outline : name of the file containing the brain surface outline 
            overlaid onto original image. Generated only if create_outline is 
            True. 
          * brain_mask : name of the file containing the binary brain mask.
            Generated only if create_brain_mask is True.  
          * skull : name of the file containing the approximate skull image.
            Generated only if create_skull is True.
    """
    
    def __init__(self, input=None, output=None, 
                 create_outline=False, create_brain_mask=False, create_skull=False,
                 intensity_threshold=None, gradient_threshold=None, 
                 head_radius=None, center_of_gravity=None, 
                 robust=False,
                 *args, **kwargs):
        
        super(BET, self).__init__(*args, **kwargs)
        
        self.input = input
        self.output = output
        
        self.create_outline = create_outline
        self.create_brain_mask = create_brain_mask
        self.create_skull = create_skull
        
        self.intensity_threshold = intensity_threshold
        self.gradient_threshold = gradient_threshold
        
        self.head_radius = head_radius
        self.center_of_gravity = center_of_gravity
        
        self.robust = robust
    
    def _get_outline(self):
        if self.create_outline :
            return os.path.extsep.join([os.path.splitext(self.input)[0]+"_overlay", 
                                        self._suffixes[self.fsl_environment["FSLOUTPUTTYPE"]]])
        else :
            return None
    
    def _get_brain_mask(self):
        if self.create_brain_mask :
            base_name = os.path.splitext(self.output)[0]
            if base_name.endswith(".nii") :
                base_name = os.path.splitext(base_name)[0]
            return os.path.extsep.join([base_name+"_mask", 
                                        self._suffixes[self.fsl_environment["FSLOUTPUTTYPE"]]])
        else :
            return None
    
    def _get_skull(self):
        if self.create_skull :
            return os.path.extsep.join([os.path.splitext(self.input)[0]+"_skull", 
                                        self._suffixes[self.fsl_environment["FSLOUTPUTTYPE"]]])
        else :
            return None
    
    def _get_command(self):
        if None in [self.input, self.output] :
            raise Exception("Both input and output must be specified") 
        command = ["bet", self.input, self.output]
        
        if self.create_outline :
            command.append("-o")
        if self.create_brain_mask :
            command.append("-m")
        if self.create_skull :
            command.append("-s")
        
        if self.intensity_threshold :
            command.extend(["-f", "{0}".format(self.intensity_threshold)])
        if self.gradient_threshold :
            command.extend(["-g", "{0}".format(self.gradient_threshold)])
        
        if self.head_radius :
            command.extend(["-r", "{0}".format(self.head_radius)])
        if self.center_of_gravity :
            command.extend(["-c", "{0}".format(self.center_of_gravity[0]), 
                            "{0}".format(self.center_of_gravity[1]),
                            "{0}".format(self.center_of_gravity[2])])
        
        if self.robust:
            command.append("-R")
        
        return command

    outline = property(_get_outline)
    brain_mask = property(_get_brain_mask)
    skull = property(_get_skull)

def bet(input, create_outline=False, create_brain_mask=False, create_skull=False,
        intensity_threshold=0.5, gradient_threshold=None, 
        head_radius=None, center_of_gravity=None, robust=False, *args, **kwargs) :
    """ BET (Brain Extraction Tool) from FSL.
    
        <gui>
            <item name="input" type="Image" label="Input" />
            <item name="create_outline" type="Bool" label="Create outline"
                  initializer="False" />
            <item name="create_brain_mask" type="Bool" label="Create brain mask"
                  initializer="False" />
            <item name="create_skull" type="Bool" label="Create skull"
                  initializer="False" />
            <item name="intensity_threshold" type="Float" label="Intensity threshold"
                  initializer="0.5" />
            <item name="gradient_threshold" type="Float" label="Gradient threshold"
                  initializer="0" />
            <item name="robust" type="Bool" label="Robust center estimation"
                  initializer="False" />
            <item name="output" type="Image" label="Output" 
                  initializer="output=True" role="return" />
        </gui>
    """
    
    fd, input_filename = tempfile.mkstemp(".nii.gz")
    os.close(fd)
    medipy.io.save(input, input_filename)
    
    fd, output_filename = tempfile.mkstemp(".nii.gz")
    os.close(fd)

    bet = BET(input_filename, output_filename, create_outline, create_brain_mask,
              create_skull, intensity_threshold, gradient_threshold, head_radius,
              center_of_gravity, robust, *args, **kwargs)
    bet()
    
    output = medipy.io.load(output_filename)
    
    os.unlink(input_filename)
    os.unlink(output_filename)
    
    return output
