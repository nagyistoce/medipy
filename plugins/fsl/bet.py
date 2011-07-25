import os

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
                 intensity_threshold=0.5, gradient_threshold=0, 
                 head_radius=None, center_of_gravity=None,
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
    
    def _get_outline(self):
        return os.path.extsep.join(os.path.splitext(self.input)+"_overlay", 
                                   self._suffixes[self.fsl_output_type])
    
    def _get_brain_mask(self):
        return os.path.extsep.join(os.path.splitext(self.input)+"_mask", 
                                   self._suffixes[self.fsl_output_type])
    
    def _get_skull(self):
        return os.path.extsep.join(os.path.splitext(self.input)+"_skull", 
                                   self._suffixes[self.fsl_output_type])
    
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
        
        command.append("-f {0:f}".format(self.intensity_threshold))
        command.append("-g {0:f}".format(self.gradient_threshold))
        
        if self.head_radius :
            command.append("-r {0:f}".format(self.head_radius))
        if self.center_of_gravity :
            command.append("-c {0[0]:f} {0[1]:f} {0[2]:f}".format(self.center_of_gravity))
        
        return command
