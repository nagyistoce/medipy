import copy
import os
import subprocess

import medipy.base

class FSLTool(object):
    """ Base class for wrappers around FSL tools.
        
        FSLOUTPUTTYPE can be specified in the environment or in the wrapper.
        If both are defined, the definition from the wrapper will override the
        one from the environment. FSLOUTPUTTYPE must be one of :
          * NIFTI_PAIR
          * NIFTI
          * NIFTI_GZ
          * NIFTI_PAIR_GZ
    """
    
    _suffixes = {
        "NIFTI_PAIR" : ".hdr",
        "NIFTI" : ".nii",
        "NIFTI_GZ" : ".nii.gz",
        "NIFTI_PAIR_GZ" : ".hdr.gz",
    }
    
    def __init__(self, fsl_dir=None, fsl_output_type="NIFTI"):
        self.fsl_dir = fsl_dir
        self.fsl_output_type = fsl_output_type
    
    def _get_command(self):
        """ FSL command that will be run. This function must be defined in
            derived classes.
        """
        
        raise NotImplementedError()
    
    command = medipy.base.LateBindingProperty(_get_command)
    
    def __call__(self):
        """ Run the specific FSL command.
        """
        
        env = copy.deepcopy(os.environ)

        if self.fsl_dir is None :
            if "FSLDIR" not in os.environ :
                raise Exception("No FSLDIR in environment, no fsl_dir defined in wrapper.")
        else :
            env["FSLDIR"] = self.fsl_dir
            
        env.setdefault("FSLOUTPUTTYPE", self.fsl_output_type)
        command = copy.deepcopy(self.command)
        command[0] = os.path.join(env["FSLDIR"], "bin", command[0])
        subprocess.call(command, env=env)
