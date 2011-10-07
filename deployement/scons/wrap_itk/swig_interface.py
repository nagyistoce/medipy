import os
import re

from SCons.Builder import Builder

import config
from utils import configure_file

def wrap_library(env, library):
    # Nothing to do
    pass

def wrap_module(env, module, library):
    includes = []
    includes.extend(["#include <{0}>".format(x) for x in library.includes])
    for wrapped_class in module.classes :
        includes.extend(["#include <{0}>".format(x) 
                         for x in wrapped_class.includes])
    includes = "\n".join(includes)
    
    typedefs = []
    for wrapped_class in module.classes :
        for typedef in wrapped_class.typedefs :
            typedefs.append("      typedef {0} {1};".format(
                typedef.full_name, typedef.mangled_name))
    typedefs = "\n".join(typedefs)
    def swig_interface_in_action(target, source, env):
        configure_file(source[0].path, target[0].path,
                       SWIG_INTERFACE_INCLUDES_CONTENT=includes,
                       SWIG_INTERFACE_TYPEDEFS=typedefs)
    source = os.path.join(config.wrapitk_root, "Configuration", "Languages", 
                          "SwigInterface", "module.includes.in")
    # TODO : path
    target = "{0}SwigInterface.h.in".format(module.name)
    Builder(action=swig_interface_in_action)(env, target, source)
    
    cableidx_action = "cableidx $SOURCE $TARGET"
    source = module.xml_node
    # TODO : path
    target = "{0}.idx".format(module.name)
    Builder(action=cableidx_action)(env, target, source)
    
    # TODO : igenerator
    
def wrap_class(env, wrapped_class, module, library):
    pass

class SwigInterface(object) :
    
    def __init__(self, env) :
        self.env = env
        
        self.library_depends = None
        
        self.IdxFiles = {}
        
#        self._mdx_content = None
#        self._module_content = None
#        self._idx_files = None
#        self._modules = None
        self.includes = None
        self.typedefs = None
    
    def wrap_library(self, library) :
        pass
        self._mdx_content = ""
#        self._module_content = ""
#        self._idx_files = []
#        self._modules = None
    
    def end_wrap_library(self) :
        # TODO ?
#        # Loop over the extra swig input files and copy them to the Typedefs directory
#        foreach(source ${WRAPPER_LIBRARY_SWIG_INPUTS})
#            get_filename_component(basename ${source} NAME)
#            set(dest "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/${basename}")
#            exec_program(${CMAKE_COMMAND} ARGS -E copy_if_different "\"${source}\"" "\"${dest}\"")
#            WRAP_ITK_INSTALL("/Configuration/Typedefs" "${dest}")
#            #    set(SWIG_INTERFACE_MODULE_CONTENT "${SWIG_INTERFACE_MODULE_CONTENT}%import ${basename}\n")
#        endforeach(source)
        files = []
        deps = []
        SwigFiles = {}
        for dep in self.library_depends :
            deps.extend(self.IdxFiles.get(dep, []))
            deps.extend(SwigFiles.get(dep, []))
            mdx_content = "{0}.mdx".format(dep)+mdx_content
        
    
    def wrap_module(self, module) :
        self.includes = []
        self.typedefs = []
    
    def end_wrap_module(self, module) :
        includes_content = ""
        
        if self.includes :
            for include_file in list(set(self.includes)) :
                if re.match("<.*>", include_file) :
                    includes_content += "#include {0}\n".format(include_file)
                else :
                    includes_content += "#include \"{0}\"\n".format(include_file)
        
        def swig_interface_h_in_action(target, source, env):
            configure_file(
                os.path.join(self.env["WRAP_ITK_SWIGINTERFACE_SOURCE_DIR"], 
                             "module.includes.in"),
                target[0].path,
                SWIG_INTERFACE_INCLUDES_CONTENT=includes_content,
                SWIG_INTERFACE_TYPEDEFS="\n".join(self.typedefs))
        includes_file = "{0}SwigInterface.h.in".format(module)
        Builder(action=swig_interface_h_in_action)(self.env, includes_file, "")
    
    def wrap_include(self, include_file) :
        self.includes.append(include_file)
    
    def add_one_typedef(self, wrap_class, swig_name, wrap_method="", template_params=None) :
        base_name = wrap_class.split("::")[-1]
        
        template_parameters = template_params
        if template_parameters :
            full_class_name = "{0}< {1} >".format(wrap_class, ", ".join(template_parameters))
        else :
            full_class_name = wrap_class
        
        if "2_SUPERCLASSES" in wrap_method :
            self.add_simple_typedef("{0}::Superclass::Superclass".format(full_class_name), 
                                    "{0}_Superclass_Superclass".format(swig_name))
            self.add_simple_typedef("{0}::Superclass::Superclass::Pointer".format(full_class_name), 
                                    "{0}_Superclass_Superclass_Pointer".format(swig_name))
        if "SUPERCLASS" in wrap_method :
            self.add_simple_typedef("{0}::Superclass".format(full_class_name), 
                                    "{0}_Superclass".format(swig_name))
            self.add_simple_typedef("{0}::Superclass::Pointer".format(full_class_name), 
                                    "{0}_Superclass_Pointer".format(swig_name))
        
        self.add_simple_typedef(full_class_name, swig_name)
        
        if "POINTER" in wrap_method :
            if wrap_method == "AUTOPOINTER" :
                self.add_simple_typedef("{0}::SelfAutoPointer".format(full_class_name),
                                        "{0}_AutoPointer".format(swig_name))
            else :
                self.add_simple_typedef("{0}::Pointer".format(full_class_name),
                                        "{0}_Pointer".format(swig_name))
    
    def add_simple_typedef(self, wrap_class, swig_name) :
        self.typedefs.append("typedef {0} {1};".format(wrap_class, swig_name))
