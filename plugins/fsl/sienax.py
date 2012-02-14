import os
import re

from fsl_tool import FSLTool

class Sienax(FSLTool):
    """ Wrapper for Sienax from FSL.
    
        Inputs
          * input
          * output_directory : default to <input>_sienax
          * two_class_segmentation : don't segment grey and white matter 
            separately
          * t2_weighted_input : T2-weighted input image (default T1-weighted)
          * ignore_upwards : ignore from t (mm) upwards in MNI152/Talairach
            space
          * ignore_downwards : ignore from b (mm) downwards in MNI152/Talairach
            space (b should probably be negative)
          * regional : regional - use standard-space masks to give peripheral
            cortex GM volume (3-class segmentation only) and ventricular CSF volume
          * lesion_mask : use lesion (or lesion+CSF) mask to remove incorrectly
            labelled "grey matter" voxels
          * bet_options : options to pass to BET brain extraction (inside
            double-quotes), e.g. -B "-f 0.3"
          * fast_options : options to pass to FAST segmentation (inside
            double-quotes), e.g. -S "-i 20"
    """
    
    def __init__(self, input=None, output_directory=None,
                 two_class_segmentation=False, t2_weighted_input=False,
                 ignore_upwards=None, ignore_downwards=None,
                 regional=False, lesion_mask=None,
                 bet_options=None, fast_options=None,
                 *args, **kwargs):
        
        super(Sienax, self).__init__(*args, **kwargs)
        
        self.input = input
        self.output_directory = output_directory
        
        self.two_class_segmentation = two_class_segmentation
        self.t2_weighted_input = t2_weighted_input
        
        self.ignore_upwards = ignore_upwards
        self.ignore_downwards = ignore_downwards
        
        self.regional = regional
        self.lesion_mask = lesion_mask
        
        self.bet_options = bet_options
        self.fast_options = fast_options
    
    def parse_report(self, path=None):
        """ Return the volume informations contained in the SIENAX report. This
            is a dictionary with keys "version" (version number, as string), 
            "vscale" (scale between image and MNI152), "pgrey" (peripheral gray
            matter), "vcsf" (ventricular CSF), "grey", "white", and "brain". 
            The informations for the different tissues is a dictionary with the
            normalized and raw values, in cubic millimeters.
            
            If path is None, use the file "report.sienax" in the default output
            directory. 
        """
        
        if path is None :
            path = os.path.join(self.output_directory, "report.sienax")
        
        report = {}
        
        fd = open(path)
        for line in fd.readlines() :
            version = re.match(
                r" running cross-sectional atrophy measurement: sienax version (\d+\.\d+)", 
                line)
            if version :
                report["version"] = version.group(1)
                continue
            
            for tissue in ["pgrey", "vcsf", "GREY", "WHITE", "BRAIN"] :
                pattern = tissue + r"\s+([\d+\.]+)\s+([\d+\.]+)"
                measure = re.match(pattern, line)
                if measure :
                    normalized = float(measure.group(1))
                    raw = float(measure.group(2))
                    report[tissue.lower()] = {"normalized" : normalized, "raw" : raw}
                    continue
            
            vscale = re.match("VSCALING ([\d\.]+)", line)
            if vscale :
                report["vscale"] = float(vscale.group(1))
        
        return report
    
    ##############
    # Properties #
    ##############
    
    def _get_command(self):
        
        if self.input is None :
            raise Exception("Input must be specified") 
        
        command = ["sienax", self.input]
        
        if self.output_directory :
            command.extend(["-o", "{0}".format(self.output_directory)])
        
        if self.two_class_segmentation :
            command.append("-2")
        if self.t2_weighted_input :
            command.append("-t2")
        
        if self.ignore_upwards :
            command.extend(["-t", "{0}".format(self.ignore_upwards)])
        if self.ignore_downwards :
            command.extend(["-b", "{0}".format(self.ignore_downwards)])
        
        if self.regional :
            command.append("-r")
        if self.lesion_mask :
            command.extend(["-lm", "{0}".format(self.lesion_mask)])
        
        if self.bet_options :
            command.extend(["-B", "\"{0}\"".format(self.bet_options)])
        if self.fast_options :
            command.append(["-S", "\"{0}\"".format(self.fast_options)])
        
        return command
    
    def _get_default_output_directory(self):
        root = os.path.splitext(os.path.basename(self.input))[0]
        return "{0}_sienax".format(root)
    
    default_output_directory = property(_get_default_output_directory)