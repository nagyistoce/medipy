import os

from SCons.Builder import Builder

import config
from utils import configure_file

def wrap_library(env, library) :
    # Nothing to do
    pass

def wrap_module(env, module, library):
    def instantiations_action(target, source, env) :
        
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
        
        force_instantiate = []
        for wrapped_class in module.classes :
            for typedef in wrapped_class.typedefs :
                force_instantiate.append("  sizeof({0});".format(
                    typedef.mangled_name)) 
        force_instantiate = "\n".join(force_instantiate)
        
        configure_file(source[0].path, target[0].path,
            CONFIG_WRAPPER_INCLUDES=includes,
            CONFIG_WRAPPER_MODULE_NAME=module.name,
            CONFIG_WRAPPER_TYPEDEFS=typedefs,
            CONFIG_WRAPPER_FORCE_INSTANTIATE=force_instantiate,
            CONFIG_WRAPPER_EXTRAS="")

    # TODO : path
    source = os.path.join(config.wrapitk_root, "Configuration", "Languages", 
                          "GccXML", "wrap_.cxx.in")
    target = "{0}.cxx".format(module.name)
    
    instantiations_node = Builder(action=instantiations_action)(env, target, source)
    
    include_paths = " ".join(["-I{0}".format(dir) 
                              for dir in env["CPPPATH"]])
    gcc_xml_action = (
            "gccxml -fxml-start=_cable_ -fxml=$TARGET "
            "{0} -DCSWIG "
            "-DCABLE_CONFIGURATION -DITK_MANUAL_INSTANTIATION "
            "$SOURCE".format(include_paths)
    )
    
    #self.env.Depends(xml_file, gccxml_inc_file)
    source = instantiations_node
    target = "{0}.xml".format(module.name)
    Builder(action=gcc_xml_action)(env, target, source)

def wrap_class(env, wrapped_class, module, library):
    # Nothing to do
    pass
