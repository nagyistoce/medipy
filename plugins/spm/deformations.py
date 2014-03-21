##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import medipy.base

from base import Tool

class Deformations(Tool):
    """ Utility for working with deformation fields.
    """
    
    name = "util.defs"
    
    Interpolation = medipy.base.enum("Deformations.Interpolation",
        "nearest_neighbor", "trilinear", "b_spline_2", "b_spline_3", 
        "b_spline_4", "b_spline_5", "b_spline_6", "b_spline_7"
    )
    
    class DeformationField(Tool.Config):
        def __init__(self, filename):
            Tool.Config.__init__(self)
            self._type = "def"
            self.filename = filename
        
        def _get_parameters(self):
            result = { self._type: [self.filename] }
            return result
    
    class Destination(Tool.Config):
        def __init__(self, mode="current", *args):
            Tool.Config.__init__(self)
            self.mode=mode
            if self.mode == "user_defined":
                self.output_directory = args[0]
        
        def _get_parameters(self):
            result = {}
            if self.mode == "current":
                result = {"savepwd": 1}
            elif self.mode == "image_source":
                result = {"savesrc": 1}
            elif self.mode == "deformation_source":
                result = {"savedef": 1}
            elif self.mode == "user_defined":
                result = {"saveusr": [self.output_directory]}
            
            return result
    
    def __init__(self, root, compositions, input_files, output_file, 
                 destination=None,
                 interpolation=Interpolation.trilinear):
        Tool.__init__(self, root)
        self.compositions = compositions
        self.input_files = input_files
        self.output_file = output_file
        self.destination = destination or Deformations.Destination()
        self.interpolation=interpolation
    
    def _get_script(self):
        script = []
        
        for index, composition in enumerate(self.compositions) :
            script.extend(Tool._generate_script(
                "{0}.comp{{{1}}}".format(self.name, 1+index), composition.parameters))
        
        script.extend(Tool._generate_script(self.name, {
            "ofname": self.output_file,
            "fnames": self.input_files,
        }))
        script.extend(Tool._generate_script("{0}.savedir".format(self.name), 
            self.destination.parameters))
        script.extend(Tool._generate_script(self.name, {
            "interp": self.interpolation,
        }))
        
        return script
