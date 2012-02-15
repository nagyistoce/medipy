import os
import re

from fsl_tool import FSLTool

class Siena(FSLTool):
    """ Wrapper for Siena from FSL.
    
        Inputs
          * input1
          * input2
          * output_directory : default to <input1>_to_<input2>_siena
          * bet_options : options to pass to BET brain extraction (inside
            double-quotes), e.g. -B "-f 0.3"
          * two_class_segmentation : don't segment grey and white matter 
            separately
          * t2_weighted_input : T2-weighted input image (default T1-weighted)
          * masking : use standard-space masking as well as BET
          * ignore_upwards : ignore from t (mm) upwards in MNI152/Talairach
            space
          * ignore_downwards : ignore from b (mm) downwards in MNI152/Talairach
            space (b should probably be negative)
          * siena_diff_options : options to pass to siena_diff timepoint 
            differencing (inside double-quotes), e.g. -S "-s -i 20"
    """
    
    def __init__(self, input1=None, input2=None, output_directory=None,
                 bet_options=None, two_class_segmentation=False, 
                 t2_weighted_input=False, masking=False,
                 ignore_upwards=None, ignore_downwards=None,
                 siena_diff_options=None,
                 *args, **kwargs):
        
        super(Siena, self).__init__(*args, **kwargs)
        
        self.input1 = input
        self.input2 = input
        self.output_directory = output_directory
        
        self.bet_options = bet_options
        
        self.two_class_segmentation = two_class_segmentation
        self.t2_weighted_input = t2_weighted_input
        
        self.masking = masking
        
        self.ignore_upwards = ignore_upwards
        self.ignore_downwards = ignore_downwards
        
        self.siena_diff_options = siena_diff_options
    
    def parse_report(self, path=None):
        """ Return the volume informations contained in the SIENAX report. This
            is a dictionary with keys "version" (version number, as string), 
            "correction" (dictionary with keys "A_to_B", "B_to_A", "value", and
            float values), "A_to_B" (dictionary with keys "area" (in mm^2), 
            "volume_change" (in mm^3), "ratio" (in mm) and "pbvc" (in 
            percentage points) and float values), "B_to_A" (idem), "pbvc" 
            (percentage brain volume change, in percentage point)
            
            If path is None, use the file "report.sienax" in the default output
            directory. 
        """
        
        if path is None :
            path = os.path.join(self.output_directory, "report.siena")
        
        report = {}
        
        # Siena first call "siena_diff ${B} ${A}", then "siena_diff ${A} ${B}" 
        b_to_a_read = False
        
        fd = open(path)
        for line in fd.readlines() :
            version = re.match(
                r" running longitudinal atrophy measurement: siena version (\d+\.\d+)", 
                line)
            if version :
                report["version"] = version.group(1)
                continue
            
            correction = re.match(
                r"corr1=(.*) corr2=(.*) corr=(.*)",
                line)
            if correction :
                report["correction"] = {
                    "A_to_B" : float(correction.group(1)),
                    "B_to_A" : float(correction.group(2)),
                    "value" : float(correction.group(3))
                }
        
            area = re.match(
                r"AREA  (.*) mm\^2",
                line)
            if area :
                key = "A_to_B" if b_to_a_read else "B_to_A"
                report.setdefault(key, {})["area"] = float(area.group(1))
            
            volume_change = re.match(
                r"VOLC  (.*) mm\^3",
                line)
            if volume_change :
                key = "A_to_B" if b_to_a_read else "B_to_A"
                report.setdefault(key, {})["volume_change"] = float(volume_change.group(1))
            
            ratio = re.match(
                r"RATIO (.*) mm",
                line)
            if ratio :
                key = "A_to_B" if b_to_a_read else "B_to_A"
                report.setdefault(key, {})["ratio"] = float(ratio.group(1))
            
            pbvc_partial = re.match(
                r"PBVC  (.*) %",
                line)
            if pbvc_partial :
                key = "A_to_B" if b_to_a_read else "B_to_A"
                report.setdefault(key, {})["pbvc"] = float(pbvc_partial.group(1))
                b_to_a_read = True
        
            pbvc = re.match(
                r"finalPBVC (.*) %",
                line
            )
            if pbvc :
                report["pbvc"] = float(pbvc.group(1))
        
        return report
    
    ##############
    # Properties #
    ##############
    
    def _get_command(self):
        
        if None in [self.input1, self.input2] :
            raise Exception("Both inputs must be specified") 
        
        command = ["siena", self.input1, self.input2]
        
        if self.output_directory :
            command.extend(["-o", "{0}".format(self.output_directory)])
        
        if self.bet_options :
            command.extend(["-B", "\"{0}\"".format(self.bet_options)])
        
        if self.two_class_segmentation :
            command.append("-2")
        if self.t2_weighted_input :
            command.append("-t2")
        
        if self.masking :
            command.append("-m")
        
        if self.ignore_upwards :
            command.extend(["-t", "{0}".format(self.ignore_upwards)])
        if self.ignore_downwards :
            command.extend(["-b", "{0}".format(self.ignore_downwards)])
        
        if self.siena_diff_options :
            command.append(["-S", "\"{0}\"".format(self.siena_diff_options)])
        
        return command
    
    def _get_default_output_directory(self):
        root1 = os.path.splitext(os.path.basename(self.input1))[0]
        root2 = os.path.splitext(os.path.basename(self.input2))[0]
        return "{0}_to_{1}_siena".format(root1, root2)
    
    default_output_directory = property(_get_default_output_directory)