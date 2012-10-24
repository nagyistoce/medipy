##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

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
        self._stdout = None
        self._stderr = None
        self._returncode = None
        
        self.fsl_dir = fsl_dir
        self.fsl_output_type = fsl_output_type
    
    def __call__(self):
        """ Run the specific FSL command.
        """
        
        env = copy.deepcopy(os.environ)

        if self.fsl_dir is None :
            if "FSLDIR" not in os.environ :
                raise Exception("No FSLDIR in environment, no fsl_dir defined in wrapper.")
        else :
            env["FSLDIR"] = self.fsl_dir
            
        path = env.get("PATH", "").split(os.pathsep)
        path.append(os.path.join(env["FSLDIR"], "bin"))
        env["PATH"] = os.pathsep.join(path)
        
        ld_library_path = env.get("LD_LIBRARY_PATH", "").split(os.pathsep)
        ld_library_path.insert(0, os.path.join(env["FSLDIR"], "bin"))
        env["LD_LIBRARY_PATH"] = os.pathsep.join(ld_library_path)

        env.setdefault("FSLOUTPUTTYPE", self.fsl_output_type)
        command = copy.deepcopy(self.command)
        command[0] = os.path.join(env["FSLDIR"], "bin", command[0])
        process = subprocess.Popen(
            command, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self._stdout, self._stderr = process.communicate()
        self._returncode = process.returncode

    ##############
    # Properties #
    ##############
    
    def _get_stdout(self):
        return self._stdout
    
    def _get_stderr(self):
        return self._stderr
    
    def _get_returncode(self):
        return self._returncode
    
    def _get_command(self):
        """ FSL command that will be run. This function must be defined in
            derived classes.
        """
        
        raise NotImplementedError()
    
    stdout = property(_get_stdout)
    stderr = property(_get_stderr)
    returncode = property(_get_returncode)
    command = medipy.base.LateBindingProperty(_get_command)
