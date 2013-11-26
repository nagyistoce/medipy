##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import os

import numpy

from base import Tool

class NewSegment(Tool):
    """ Extension of the default Unified segmentation
    """
    
    name = "tools.preproc8"
    
    # Default order of tissues in SPM Tissue probability map
    tissue_index = {
        "gray_matter" : 1,
        "white_matter" : 2,
        "csf" : 3,
        "bone" : 4,
        "soft_tissue" : 5,
        "air_background" : 6,
    }
    
    class Channel(Tool.Config):
        """ Processing channel (e.g. T1 or PD), containing
        
            * Volumes
            * Bias regularisation ; usual values are 
              
              * 0 (no regularisation),
              * 0.00001 (extremely light regularisation)
              * 0.0001 (very light regularisation)
              * 0.001 (light regularisation)
              * 0.01 (medium regularisation
              * 0.1 (heavy regularisation)
              * 1 (very heavy regularisation)
              * 10 (extremely heavy regularisation)

            * Bias FWHM in mm ; usual values are 30-150 in 10 mm steps, or 
              infinity
            * Save Bias Corrected ; options to save the bias field and the
              bias-corrected image
        """
        
        def __init__(self, vols=None, biasreg=0.0001, biasfwhm=60, 
                     save_bias_field=False, save_bias_corrected=False):
            self.vols = vols or []
            self.biasreg = biasreg
            self.biasfwhm = biasfwhm
            self.save_bias_field = save_bias_field
            self.save_bias_corrected = save_bias_corrected
        
        def _get_parameters(self):
            result = {
                "vols" : self.vols,
                "biasreg" : self.biasreg,
                "biasfwhm" : self.biasfwhm,
                "write" : numpy.asarray(
                    [self.save_bias_field, self.save_bias_corrected], 
                    dtype=int),
            }
            return result
    
    class Tissue(Tool.Config):
        """ Tissue class, containing :
        
            * Tissue probability map
            * Number of gaussians
            * Native space ; options to save the tissue class image in the
              original or DARTEL imported space
            * Warped Tissue ; options to the the tissue class image in the
              normalized space, with and without modulation
        """
        
        def __init__(self, tpm, ngaus=2, save_native=True, save_dartel=False,
                     save_warped_unmodulated=False, save_warped_modulated=False):
            self.tpm = tpm
            self.ngaus = ngaus
            
            self.save_native = save_native
            self.save_dartel = save_dartel
            
            self.save_warped_unmodulated = save_warped_unmodulated
            self.save_warped_modulated = save_warped_modulated
        
        def _get_parameters(self):
            result = {
                "tpm" : [self.tpm], # The tissue probability map must be a cell array
                "ngaus" : self.ngaus,
                "native" : numpy.asarray(
                    [self.save_native, self.save_dartel], 
                    dtype=int),
                "warped" : numpy.asarray(
                    [self.save_warped_unmodulated, self.save_warped_modulated],
                    dtype=int),
            }
            return result
    
    class Warp(Tool.Config):
        """ Warping options, containing :
        
            * Warping regularisation
            * Affine regularisation ; can be 

              * ICBM space template - European brains ("mni")
              * ICBM space template - East Asian brains ("eastern")
              * Average sized template ("subj")
              * No regularisation ("none")

            * Sampling distance
            * Deformation field writing, where the forward field normalises
              an image to the MNI space
        """
        
        def __init__(self, mrf=0, reg=4, affreg="mni", samp=3, 
                     save_inverse=False, save_forward=False):
            self.mrf = 0
            self.reg = 4
            self.affreg = "mni"
            self.samp = 3
            self.save_inverse = save_inverse
            self.save_forward = save_forward
        
        def _get_parameters(self):
            result = {
                "mrf": self.mrf,
                "reg" : self.reg,
                "affreg" : self.affreg,
                "samp" : self.samp,
                "write" : numpy.asarray(
                    [self.save_inverse, self.save_forward],
                    dtype=int),
            }
            return result
    
    def __init__(self, root, channels=None, tissues=None, warp=None):
        Tool.__init__(self, root)
        self.channels = channels or [NewSegment.Channel()]
        
        if tissues is None :
            # Generate default tissue classes using SPM tissue probability map
            self.tissues = []
            tpm = os.path.join(self.root, "toolbox", "Seg", "TPM.nii")
            for i in range(1,7) :
                tissue = NewSegment.Tissue("{0},{1}".format(tpm, i))
                self.tissues.append(tissue)
        else :
            self.tissues = tissues
        
        self.warp = warp or NewSegment.Warp()
    
    def _get_script(self):
        script = []
        
        for index, channel in enumerate(self.channels) :
            script.extend(Tool._generate_script(
                "{0}.channel({1})".format(self.name, 1+index), channel.parameters))
        for index, tissue in enumerate(self.tissues) :
            script.extend(Tool._generate_script(
                "{0}.tissue({1})".format(self.name, 1+index), tissue.parameters))
        script.extend(Tool._generate_script(
                "{0}.warp".format(self.name), self.warp.parameters))
        
        return script
