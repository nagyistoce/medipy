##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import os
from fsl_tool import FSLTool

class EddyCorrect(FSLTool):
    """ Wrapper for Eddy current distortion correction from FSL.

        Inputs : 
          * 4D input
          * reference_no

        Outputs :
          * 4D output
    """

    def __init__(self, input=None, reference_no=None, *args, **kwargs):
        
        super(EDDY, self).__init__(*args, **kwargs)
        
        self.input = input
        base_name = os.path.splitext(self.input)[0]
        if base_name.endswith(".nii") :
            base_name = os.path.splitext(base_name)[0]
        self.output = os.path.extsep.join([
            base_name+"_fsl_eddy", 
            self._suffixes[self.fsl_environment["FSLOUTPUTTYPE"]]
        ]) 
        self.reference_no = reference_no
              
       
    def _get_command(self):
        if None in [self.input, self.output, self.reference_no] :
            raise Exception("Input, output, and reference must be specified") 
        command = ["eddy_correct", self.input, self.output, self.reference_no]
               
        return command
