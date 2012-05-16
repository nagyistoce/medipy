##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import numpy

import medipy.base

class Tool(object):
    """ Base class for SPM tools.
    """
    
    class Config(object):
        """ Configuration element of an SPM tool.
        """
        
        def _get_parameters(self):
            """ Return a dictionary of configuration items using the names from SPM
                batches as keys.
            """
            
            raise NotImplementedError()
        
        parameters = medipy.base.LateBindingProperty(_get_parameters)
    
    # Name of the tool
    name = None
    
    def __init__(self, root):
        
        # SPM root directory, as provided by medipy.spm.utils.find
        self.root = root
    
    def _get_script(self):
        """ Return the Matlab script to run the tool.
        """
        
        raise NotImplementedError()
    
    script = medipy.base.LateBindingProperty(_get_script)
    
    @classmethod
    def _to_matlab(self, obj):
        """ Convert a Python object to a string which can be parsed by Matlab.
            The following types are supported :
            
              * numpy.ndarray of dimensions 1 and 2 : converted to a matrix
              * list : converted to a cell array
              * string and scalar types : converted to litterals
              
            For nested structures (arrays and lists), the elements are 
            recursively converted.
        """
        
        result = ""
        if isinstance(obj, numpy.ndarray) :
            
            normalized = obj
            if normalized.ndim == 1 :
                normalized = normalized.reshape((1, normalized.shape[0]))
            
            if normalized.dtype == numpy.object :
                result += "{ "
            else :
                result += "[ "
            
            if normalized.shape[0] > 1 :
                result += "\n"
            
            for row in normalized :
                result += " ".join(Tool._to_matlab(x) for x in row)
                if normalized.shape[0] > 1 :
                    result += "\n";
            
            if normalized.dtype == numpy.object :
                result += " }"
            else :
                result += " ]"
        elif isinstance(obj, list) :
            result = "{{ {0} }}".format(" ".join(Tool._to_matlab(x) for x in obj))
        elif isinstance(obj, basestring) :
            result = "'{0}'".format(obj)
        elif numpy.isscalar(obj) :
            result = "{0}".format(obj)
        else :
            raise medipy.base.Exception(
                "Cannot convert an object of type {0}".format(repr(type(obj).__name__)))
    
        return result
    
    @classmethod
    def _generate_script(self, prefix, obj) :
        """ Recursively generate a Matlab script from a nested dictionary of
            configuration parameters.
        """
        
        script = []
        
        for key, value in obj.items() :
            if isinstance(value, dict) :
                sub_script = Tool._generate_script(
                   "{0}.{1}".format(prefix, key), value)
                script.extend(sub_script)
            else :
                script.append("{0}.{1} = {2}".format(prefix, key, Tool._to_matlab(value)))
        
        return script

def script(tools, standalone=True, exit=True, modality="fmri"):
    """ Return a Matlab script running the provided tools. If standalone is 
        True, then the script will include SPM initialization commands and the
        resulting script can be directly fed to Matlab. If standalone is False,
        the resulting script can be loaded in SPM Batch Editor.
    """
    
    script = []
    
    if standalone :
        script.extend([
            "spm('defaults','{0}');".format(modality),
            "spm_jobman('initcfg');"
        ])
    
    for index, tool in enumerate(tools) :
        prefix = "matlabbatch{{{0}}}.spm".format(index+1)
        script.extend(["{0}.{1};".format(prefix, statement) 
                       for statement in tool.script])
    
    if standalone :
        script.append("spm_jobman('run',matlabbatch);")
        if exit :
            script.append("exit();")
    
    return "\n".join(script)
