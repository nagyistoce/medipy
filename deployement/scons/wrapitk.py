##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

""" This module contains utility functions to mimic in Python/scons the CMake
    scripts used for WrapITK.
    
    A WrapITK module called MyModule, wrapping myClass, myOtherClass and 
    yetAnotherClass has the following structure (files in the Configuration
    directory for myOtherClass and yetAnotherClass are omitted for readability)
    
    MyModule
     |- MyModule.py : itkBase.LoadModule("MyModule")
     |- MyModulePython.py : Swig
     |- _MyModulePython.so
     |- myClassPython.py         \
     |- myOtherClassPython.py    |- Swig
     |- yetAnotherClassPython.py /
     |- Configuration
         |- MyModuleConfig.py : depends, templates
         |- MyModule.i : include the .includes from the modules MyModule depends and MyModule.includes, then itk.i and MyModule_ext.i
         |- MyModule_ext.i : declarations and calls of init_<classes>Python(), import of <classes>Python
         |- MyModule.includes : dependancies, classes, typedefs
         |- MyModule.mdx : wrap_<classes>.idx
         |- wrap_myClass.i : includes as in MyModule.i, classes declarations
         |- wrap_myClass_ext.i : import wrap_pyBase.i, include wrap_myClass_doc.i, DECLARE_REF_COUNT_CLASS
         |- wrap_myClass_doc.i : %feature("docstring")
         |- wrap_myClass.idx : {full_cpp_name} {swig_name} {class_name}
         |- myClass.X : man
    
    For each templated class, we need to specify the types and dimensions of the
    instantiations. The template parameters are all the same : we have for now
    no mechanism to instantiate itk::MyWonderfulFilter<itk::Image<float, 3>, itk::Image<short, 2> >
    
    In the process, we also generate the following files :
      * wrap_myClass.cxx : cable wrapper for the typedefs to each instantiation
      * wrap_myClass.xml : gccxml-ized version of the above
      * the .cpp files for the above Swig files
      * the .os files for the above .cpp files (but not for wrap_myClass.cxx)
"""

import itertools
import os
import os.path
import re

import itkTypes

from SCons.Builder import Builder

import wrapitk_templates

def get_instantiated_types(wrapitk_root):
    """ Return the types WrapITK was compiled for.
    """
    
    result = []
    
    config_file = open(os.path.join(wrapitk_root, "WrapITKConfig.cmake"))
    for line in config_file.readlines() :
        match = re.match(r"SET\(WRAP_[a-z_]+ ON CACHE BOOL "
                         "\"Wrap ([a-z ]+) type\"\)", line)
        if match :
            result.append(match.group(1))
    return result
    

def get_unsigned_int_types(wrapitk_root):
    """ Return the unsigned integer types for which WrapITK was compiled.
    """
    
    result = []
    for type in get_instantiated_types(wrapitk_root) :
        if type.startswith("unsigned ") :
           result.append(type)
    return result

def get_signed_int_types(wrapitk_root):
    """ Return the signed integer types for which WrapITK was compiled.
    """
    
    result = []
    for type in get_instantiated_types(wrapitk_root) :
        if type.startswith("signed ") :
           result.append(type)
    return result

def get_real_types(wrapitk_root):
    """ Return the real types for which WrapITK was compiled.
    """
    
    result = []
    for type in get_instantiated_types(wrapitk_root) :
       if type in ["float", "double"]:
           result.append(type)
    return result

def get_int_types(wrapitk_root):
    """ Return the integer types for which WrapITK was compiled.
    """
    
    return get_unsigned_int_types(wrapitk_root)+get_signed_int_types(wrapitk_root)

def get_scalar_types(wrapitk_root):
    """ Return the scalar types for which WrapITK was compiled.
    """
    
    return get_int_types(wrapitk_root)+get_real_types(wrapitk_root)

def get_dimensions(wrapitk_root):
    """ Return the dimensions for which WrapITK was compiled.
    """
    
    result = []
    
    config_file = open(os.path.join(wrapitk_root, "WrapITKConfig.cmake"))
    for line in config_file.readlines() :
        match = re.match(r"SET\(WRAP_ITK_DIMS \"([0-9;]+)\" CACHE STRING "
                         "\"dimensions available separated by semicolons "
                         "\(;\)\"\)", line)
        if match :
            value = match.group(1)
            dimensions = [int(x) for x in value.split(";")]
            result.extend(dimensions)
    return result

############################################
# SCons build actions for the module class #
############################################

def module_includes(target, source, env):
    requirements = env["REQUIREMENTS"]
    classes_template_info = env["TEMPLATES"]
    pointer = env["POINTER"]
    content = get_module_includes(requirements, classes_template_info, pointer)
    file = open(str(target[0]), "w")
    file.write(content)
    file.close()

def module_swig_extra(target, source, env) :
    classes = env["CLASSES"]
    requirements = env["REQUIREMENTS"]
    
    declarations = ""
    calls = ""
    imports = ""
    for klass in classes :
        swig_class_name = get_swig_class_name(klass)
        declarations += "extern \"C\" int init_%sPython();\n"%swig_class_name
        calls += "init_%sPython();\n"%swig_class_name
    for requirement in requirements :
        imports += "import %sPython\n"%requirement
    for klass in classes :
        swig_class_name = get_swig_class_name(klass)
        imports += "from %sPython import *\n"%swig_class_name
        
    content = wrapitk_templates.module_ext_i%locals()
    
    file = open(str(target[0]), "w")
    file.write(content)
    file.close()

def module_swig(target, source, env) :
    requirements = env["REQUIREMENTS"]
    module_name = env["MODULE_NAME"]
    
    includes = "\n".join(["#include \"%s.includes\"\n"%x for x in requirements])
    content = wrapitk_templates.module_i%{"module_name":module_name, "includes":includes, "interface_content":""}
    
    file = open(str(target[0]), "w")
    file.write(content)
    file.close()
    
def module_loader(target, source, env):
    module_name = env["MODULE_NAME"]
    
    content = wrapitk_templates.module_py%{"module_name":module_name}
    file = open(str(target[0]), "w")
    file.write(content)
    file.close()

def module_config(target, source, env):
    requirements = env["REQUIREMENTS"]
    classes_template_info = env["TEMPLATES"]
    
    content = get_module_config_templates(requirements, classes_template_info)
    file = open(str(target[0]), "w")
    file.write(content)
    file.close()

###############################################
# SCons build actions for the CableSwig files #
###############################################

def class_wrapper(target, source, env) :
    """ cf. class_cable_swig_files_builder for doc about pointer
    """
    
    instantiations = env["INSTANTIATIONS"]
    class_name = env["CLASS_NAME"]
    pointer = env["POINTER"]
    
    swig_class_name = get_swig_class_name(class_name)
    unqualified_class_name = get_unqualified_class_name(class_name)
    
    includes = set()
    for instantiation in instantiations :
        for parameter in instantiation :
            includes.add("<%s.h>"%get_swig_class_name(parameter[0]))
    includes = list(includes)
    includes.append("%s.h"%swig_class_name)
    
    template_parameters = []
    for instantiation in instantiations :
        parameter_mangled_names = []
        for entry in instantiation :
            mangled_template_parameters = get_mangled_template_name(entry[0], *entry[1:])
            parameter_mangled_names.append(mangled_template_parameters)
        template_parameters.append(parameter_mangled_names)
    
    typedefs = []
    for parameter_set in template_parameters :
        cpp_name = [p[0] for p in parameter_set]
        cpp_name = ", ".join(cpp_name)
        cpp_name = "%s< %s >::%s"%(class_name, cpp_name, unqualified_class_name)
        
        mangled_name = [p[1] for p in parameter_set]
        swig_name = "%s%s"%(swig_class_name, "".join(mangled_name))
        
        typedefs.append((cpp_name, swig_name))

        if pointer == "pointer" :
            pointer_cpp_name = "%s::Pointer::SmartPointer"%cpp_name
            pointer_swig_name = "%s_Pointer"%swig_name
            typedefs.append((pointer_cpp_name, pointer_swig_name))
        elif pointer == "pointer_with_superclass" :
            superclass_cpp_name = "%s::Superclass::Self"%cpp_name
            superclass_swig_name = "%s_Superclass"%swig_name
            
            superclass_pointer_cpp_name = "%s::Superclass::Pointer::SmartPointer"%cpp_name
            superclass_pointer_swig_name = "%s_Superclass_Pointer"%swig_name
            
            typedefs.append((superclass_cpp_name, superclass_swig_name))
            typedefs.append((superclass_pointer_cpp_name, superclass_pointer_swig_name))
    
    wrapper = get_cable_wrapper(includes, swig_class_name, typedefs)
    
    file = open(str(target[0]), "w")
    file.write(wrapper)
    file.close()

##########################################
# SCons build actions for the Swig files #
##########################################

def master_index(target, source, env):
    content = "\n".join([x.name for x in source])
    file = open(str(target[0]), "w")
    file.write(content)
    file.close()
    
    return 0

def class_swig(source, target, env, for_signature) : 
    # This is a generator, hence the different prototype
    
    swig_class_name = get_swig_class_name(env["CLASS_NAME"])
    
    wrapitk_root = env["WRAPITK_ROOT"]
    
    xml_file = "wrap_%s.xml"%swig_class_name
    user_swig_file = "%s.swg"%swig_class_name
    swig_ext_file = "wrap_%s_ext.i"%swig_class_name
    module_includes_file = "%s.includes"%env["MODULE_NAME"]
    
    command = ["python",
        wrapitk_root+"/Configuration/Languages/SwigInterface/igenerator.py",
        "--mdx", wrapitk_root + "/Configuration/Typedefs/Base.mdx",
        "--include", "Base.includes",
        "--include", module_includes_file,
        "--swig-include", "itk.i",
        "--swig-include", user_swig_file,
        "--mdx", env["MASTER_INDEX"].path,
        "--swig-include", swig_ext_file, 
        "-w1", "-w3", "-w51", "-w52", "-w53", "-w54",
        str(source[0]), str(target[0])]
    
    return " ".join(command)
    
def class_swig_extra(target, source, env) :
    
    class_name = env["CLASS_NAME"]
    pointer = env["POINTER"]
    instantiations = env["INSTANTIATIONS"]
    
    swig_class_name = get_swig_class_name(class_name)
    if pointer == "pointer" :
        ref_count_declarations = ""
        for instantiation in instantiations :
            name = swig_class_name
            for parameter in instantiation :
                mangled_name = get_mangled_template_name(parameter[0], *parameter[1:])[1]
                name += mangled_name
            ref_count_declarations += "DECLARE_REF_COUNT_CLASS(%s)\n"%name
    else :
        ref_count_declarations = ""
    
    content = wrapitk_templates.wrap_class_ext_i%{
        "swig_class_name" : swig_class_name,
        "ref_count_declarations" : ref_count_declarations
    }
    file = open(str(target[0]), "w")
    file.write(content)
    file.close()

###################################################
# SCons build actions for the documentation files #
###################################################

def class_doxygen_config(target, source, env):
    
    dirname = os.path.dirname(__file__)
    template_file = open(os.path.join(dirname, "wrapitk_doxygen.config.in"))
    template = template_file.read()
    template_file.close()
    
    output_directory = os.path.dirname(source[0].path)
    
    headers = ""
    for s in source :
        headers += "\"%s\"\\\n"%str(s)
    
    generate_man = "NO"
    
    content = template%locals()
    
    file = open(str(target[0]), "w")
    file.write(content)
    file.close()

def class_doxy2swig_conf(target, source, env):
    class_name = env["CLASS_NAME"]
    instantiations = env["INSTANTIATIONS"]
    
    swig_class_name = get_swig_class_name(class_name)
    
    doxygen_class_name = re.sub(":", "_1", re.sub(r"([A-Z])", r"_\1", class_name)).lower()
    doxygen_file_name = "class%s.xml"%doxygen_class_name
    doxygen_file_name = os.path.join(os.path.dirname(target[0].abspath), "xml", doxygen_file_name)
    
    template_parameters = []
    for instantiation in instantiations :
        parameter_mangled_names = []
        for entry in instantiation :
            mangled_template_parameters = get_mangled_template_name(entry[0], *entry[1:])
            parameter_mangled_names.append(mangled_template_parameters)
        template_parameters.append(parameter_mangled_names)
    
    typedefs = []
    for parameter_set in template_parameters :
        
        mangled_name = [p[1] for p in parameter_set]
        swig_name = "%s%s"%(swig_class_name, "".join(mangled_name))
        
        typedefs.append(swig_name)
    
    content = "\n" + doxygen_file_name + "\t" + class_name + "\t"
    content += "\t".join(typedefs)
    
    file = open(str(target[0]), "w")
    file.write(content)
    file.close()

def class_swig_doc(source, target, env, for_signature) : 
    # This is a generator, hence the different prototype
    wrapitk_root = env["WRAPITK_ROOT"]
    
    itk_doxy2swig = os.path.join(wrapitk_root, "Configuration", "Languages", 
                                 "Doc", "itk_doxy2swig.py")
    
    doxy2swig_command = "python %s $SOURCE $TARGET"%itk_doxy2swig
    return doxy2swig_command

###################
# Pseudo builders #
###################

def module_files_builder(env, module_name, templates, requirements):
    """ Builder for module-specific files : Module.py, ModuleConfig.py, 
        Module.i, Module_ext.i, Module.includes
    """
    
    builders = {
        "module_loader" : Builder(action=module_loader),
        "module_config" : Builder(action=module_config),
        "module_swig" : Builder(action=module_swig),
        "module_swig_extra" : Builder(action=module_swig_extra),
        "module_includes" : Builder(action=module_includes),
    }
    
    targets = {
        "module_loader" : "%s.py"%module_name,
        "module_config" : "%sConfig.py"%module_name,
        "module_swig" : "_%s.i"%module_name,
        "module_swig_extra" : "_%s_ext.i"%module_name,
        "module_includes" : "%s.includes"%module_name, 
    }
    
    extra_arguments = {
        "module_loader" : { "MODULE_NAME" : module_name},
        "module_config" : { "REQUIREMENTS" : requirements, 
                            "TEMPLATES" : templates},
        "module_swig" : { "MODULE_NAME" : module_name, 
                          "REQUIREMENTS" : requirements},
        "module_swig_extra" : { "CLASSES" : [x[0] for x in templates], 
                                "REQUIREMENTS" : requirements},
        "module_includes" : { "POINTER" : False, "REQUIREMENTS" : requirements, 
                              "TEMPLATES" : templates},
    } 
    
    nodes = {}
    
    for name, builder in builders.items() :
        nodes[name] = builder(env, targets[name], [], **extra_arguments[name])
    
    env.Depends(nodes["module_swig"], nodes["module_swig_extra"])
    
    return nodes

def class_cable_swig_files_builder(env, class_name, header, instantiations, pointer=None):
    """ Builder for class-specific CableSwig files : wrap_MyClass.cxx, 
        wrap_MyClass.xml, wrap_MyClass.idx
        
        pointer can be None, "pointer" or "pointer_with_superclass", following
        WrapITK's semantics (Documentation/Guide.txt) :
            If POINTER is passed, then the class and the typedef'd class::Pointer
            type is wrapped. (Class::Pointer had better be a SmartPointer
            instantiation, or things won't work. This is always the case for
            ITK-style code.) If POINTER_WITH_SUPERCLASS is provided, then
            class::Pointer, class::Superclass and class::Superclass::Pointer are
            all wrapped. (Again, this only works for ITK-style code where the
            class has a typedef'd Superclass, and the superclass has Self and
            Pointer typedefs.)
    """
    
    swig_class_name = get_swig_class_name(class_name)
    
    includes = set(env["ITK"]["CPPPATH"] + env["Python"]["CPPPATH"] + env["CPPPATH"])
    includes = [x.replace("\\", "/") for x in includes]
    includes = [x.rstrip("/") for x in includes]
    include_directives = " ".join(["-I%s"%x for x in includes])
    
    class_xml = ("gccxml -fxml-start=_cable_ -fxml=$TARGET -DCSWIG " 
                 "-DCABLE_CONFIGURATION -DITK_MANUAL_INSTANTIATION "
                 "$INCLUDE_DIRECTIVES $SOURCE")
    
    builders = {
        "class_wrapper" : Builder(action=class_wrapper),
        "class_xml" : Builder(action=class_xml),
        "class_index" : Builder(action="cableidx $SOURCE $TARGET"),
    }
    
    targets = {
        "class_wrapper" : "wrap_%s.cxx"%swig_class_name,
        "class_xml" : "wrap_%s.xml"%swig_class_name,
        "class_index" : "wrap_%s.idx"%swig_class_name,
    }
    
    extra_arguments = {
        "class_wrapper" : { "INSTANTIATIONS" : instantiations,
                            "CLASS_NAME" : class_name,
                            "POINTER" : pointer }, 
        "class_xml" : { "INCLUDE_DIRECTIVES" : include_directives },
        "class_index" : { },
    } 
    
    nodes = {}
    
    source = header
    # Build order is important, cannot use dictionary iterator
    for name in ["class_wrapper", "class_xml", "class_index"] :
        builder = builders[name]
        nodes[name] = builder(env, targets[name], source, **extra_arguments[name])
        source = nodes[name]
    
    return nodes

def class_swig_files_builder(env, class_name, instantiations, pointer, class_nodes, 
                             module_name, module_nodes, wrapitk_root) :
    """ Builder for class-specific Swig files : wrap_MyClass.i and 
        wrap_MyClass_ext.i 
    """
    
    swig_class_name = get_swig_class_name(class_name)
    
    builders = {
        "class_swig" : Builder(generator=class_swig),
        "class_swig_extra" : Builder(action=class_swig_extra),
    }
    
    targets = {
        "class_swig" : "wrap_%s.i"%swig_class_name,
        "class_swig_extra" : "wrap_%s_ext.i"%swig_class_name,
    }
    
    extra_arguments = {
        "class_swig" : { "CLASS_NAME" : class_name,
                         "MODULE_NAME" : module_name, 
                         "MASTER_INDEX" : module_nodes["master_index"][0],
                         "WRAPITK_ROOT" : wrapitk_root }, 
        "class_swig_extra" : { "CLASS_NAME" : class_name,
                               "INSTANTIATIONS" : instantiations,
                               "POINTER" : pointer },
    } 
    
    nodes = {}
    
    for name, builder in builders.items() :
        nodes[name] = builder(env, targets[name], class_nodes["class_xml"], **extra_arguments[name])
    
    env.Depends(nodes["class_swig"], 
                [module_nodes["module_includes"], module_nodes["master_index"][0],
                 class_nodes["user_swig"], nodes["class_swig_extra"]])
    
    return nodes

def class_doc_files_builder(env, class_name, header, instantiations, wrapitk_root) :
    """ Builder for class-specific documentation files : MyClass_doxygen.config,
        xml/classmyclass.xml, MyClass.conf, wrap_MyClass_doc.i
    """
    
    swig_class_name = get_swig_class_name(class_name)
    doxygen_class_name = re.sub(":", "_1", re.sub(r"([A-Z])", r"_\1", class_name)).lower()
    
    builders = {
        "class_doxygen_config" : Builder(action=class_doxygen_config),
        "class_doxygen" : Builder(action="doxygen $SOURCE"),
        "class_doxy2swig_conf" : Builder(action=class_doxy2swig_conf),
        "class_swig_doc" : Builder(generator=class_swig_doc),
    }
    
    targets = {
        "class_doxygen_config" : "%s_doxygen.config"%swig_class_name,
        "class_doxygen" : "xml/class%s.xml"%doxygen_class_name,
        "class_doxy2swig_conf" : "%s.conf"%swig_class_name,
        "class_swig_doc" : "wrap_%s_doc.i"%swig_class_name
    }
    
    sources = {
        "class_doxygen_config" : header,
        "class_doxygen" : [targets["class_doxygen_config"], header],
        "class_doxy2swig_conf" : header,
        "class_swig_doc" : targets["class_doxy2swig_conf"]
    }
    
    extra_arguments = {
        "class_doxygen_config" : { "MODULE_NAME" : "module_name" }, # useful ?
        "class_doxygen" : { }, 
        "class_doxy2swig_conf" : { "CLASS_NAME" : class_name, 
                                   "INSTANTIATIONS" : instantiations },
        "class_swig_doc" : { "WRAPITK_ROOT" : wrapitk_root }
    } 
    
    nodes = {}
    
    # Build order is important, cannot use dictionary iterator
    for name in ["class_doxygen_config", "class_doxygen", "class_doxy2swig_conf", "class_swig_doc"] :
        builder = builders[name]
        nodes[name] = builder(env, targets[name], sources[name], **extra_arguments[name])
        source = nodes[name]
    
    return nodes

###################################################
# Utility functions for the SCons build functions #
###################################################

# Dictionary to map from C++ name to mangled WrapITK name, 
# e.g. "unsigned short" -> "US"
mangled_type_map = {}
    
# Build the mangled type map from the names in the itkTypes module. Objects
# which are type should have a name and a short_name. 
for name in dir(itkTypes) :
    object = getattr(itkTypes, name)
    if hasattr(object, "name") and hasattr(object, "short_name") :
        mangled_type_map[object.name] = object.short_name

# From WrapITKTypes.cmake
mangled_type_map["itk::Offset"] = "O"
mangled_type_map["itk::Vector"] = "V"
mangled_type_map["itk::CovariantVector"] = "CV"
mangled_type_map["itk::ContinuousIndex"] = "CI"
mangled_type_map["itk::Array"] = "A"
mangled_type_map["itk::FixedArray"] = "FA"
mangled_type_map["itk::RGBPixel"] = "RGB"
mangled_type_map["itk::RGBAPixel"] = "RGBA"
mangled_type_map["std::complex"] = "C"
mangled_type_map["itk::Image"] = "I"
mangled_type_map["itk::VectorImage"] = "VI"
mangled_type_map["itk::VariableLenghtVector"] = "VLV"
mangled_type_map["itk::Point"] = "P"
mangled_type_map["itk::LevelSetNode"] = "LSN"
mangled_type_map["itk::FlatStructuringElement"] = "SE"
mangled_type_map["itk::SpatialObject"] = "SO"
mangled_type_map["itk::Statistics::Histogram"] = "H"
mangled_type_map["itk::Matrix"] = "M"

def get_instantiations(class_name, *parameters) :
    """ Return all instantiations of template class with given set of parameters
        >>> get_instantiations("Image", ["float", "short"], [2, 3]) #doctest: +NORMALIZE_WHITESPACE
        [('Image', 'float', 2), ('Image', 'float', 3), 
         ('Image', 'short', 2), ('Image', 'short', 3)]
        
        >>> get_instantiations("float")
        ['float']
        
        This function can also be used recursively :
        >>> get_instantiations("vector", ["float", ["vector", ["int", "char"]]])  #doctest: +NORMALIZE_WHITESPACE
        [('vector', 'float'), 
         ('vector', ('vector', 'int')), 
         ('vector', ('vector', 'char'))]
        
    """

    if parameters :
        instantiated_parameters = []
        for parameter in parameters :
            instantiated_parameter = []
            for entry in parameter :
                if isinstance(entry, (list, tuple)) :
                    instantiated_parameter.extend(get_instantiations(*entry))
                else :
                    instantiated_parameter.extend(get_instantiations(entry))
            instantiated_parameters.append(instantiated_parameter)
        return [(class_name,)+x for x in itertools.product(*instantiated_parameters)]
    else :
        return [class_name]

def get_unqualified_class_name(class_name) :
    """ Return the non-qualified name of a C++ class name by stripping the 
        namespaces
        
        >>> get_unqualified_class_name("itk::Numerics::MyMarvelousOptimizer")
        'MyMarvelousOptimizer'
    """
    # Split the class name according to its namespaces
    class_name_elements = class_name.split("::")
    return class_name_elements[-1]

def get_swig_class_name(class_name) :
    """ Return the Swig name of a C++ class by deleting the scope separator ("::")
        operators from the name
        
        >>> get_swig_class_name("itk::Numerics::MyMarvelousOptimizer")
        'itkNumericsMyMarvelousOptimizer'
    """
    # Split the class name according to its namespaces
    class_name_elements = class_name.split("::")
    swig_class_name = "".join(class_name_elements)
    
    return swig_class_name

def get_mangled_template_name(class_name, *template_parameters):
    """ Return (C++ name, mangled name)
            
        >>> get_mangled_template_parameters("itk::Image", 
        ...    "unsigned short", 2) #doctest: +NORMALIZE_WHITESPACE
        ('itk::Image< unsigned short, 2 >', 'IUS2')
        
    """
    
    string_parameters = [str(x) for x in template_parameters]
    cpp_name = ("%s< %s >"%(class_name, ",".join(string_parameters)))
    
    mangled_parameters = [ mangled_type_map[x] if x in mangled_type_map else x 
                           for x in string_parameters]
    mangled_name = "%s%s"%(mangled_type_map[class_name], 
                           "".join(mangled_parameters))
    
    return (cpp_name, mangled_name)

#def get_image_filter_typedefs(class_name, types, dimensions, nb_template_parameters) :
#    typedefs = []
#    
#    swig_class_name = get_swig_class_name(class_name)
#    
#    for current_type in types :
#        for dimension in dimensions :
#            image_type = "itk::Image< %s, %i >"%(current_type, dimension)
#            template_parameter = ", ".join(nb_template_parameters*[image_type])
#            
#            cpp_name = "%s< %s >"%(class_name, template_parameter)
#            swig_name = "%sI%s%i"%(swig_class_name, mangled_type_map[current_type], dimension)
#            
#            typedefs.append((cpp_name, swig_name))
#    
#    return typedefs

def get_wrap_className_cxx_typedefs(class_name, template_parameters) :
    """ Generate a list of typedefs for the Cable wrapper of a class. 
    """
    
    swig_class_name = get_swig_class_name(class_name)
    unqualified_class_name = get_unqualified_class_name(class_name)
    
    cable_typedefs = []
    
    for i in range(len(template_parameters[0])) :
        parameters = [x[i] for x in template_parameters]
        
        cpp_name = "%s< %s >::%s"%(class_name, ", ".join([p[0] for p in parameters]), unqualified_class_name)
        swig_name = "%s%s"%(swig_class_name, "".join([p[1] for p in parameters]))
        cable_typedefs.append((cpp_name, swig_name))
    
    return cable_typedefs

def get_cable_wrapper(includes, swig_class_name, typedefs=[]) :
    """ Generate a Cable wrapper for a specific class.
    """ 
    
    data = {"includes" : "#include <itkCommand.h>\n",
            "swig_class_name" : swig_class_name, 
            "typedefs" : "", "forced_instantiations" : "", 
            "extra" : ""}

    for include in includes : 
        if re.match(r"<.*>", include) :
            data["includes"] += "#include %s\n"%include
        else :
            data["includes"] += "#include \"%s\"\n"%include
        
    for cpp_name, swig_name in typedefs : 
        data["typedefs"] += "typedef %s %s;\n"%(cpp_name, swig_name)
    
    return wrapitk_templates.wrap_class_cxx%data

def get_module_includes(requirements, classes_template_info, pointer):
    """ Generate the .includes file for the module. 
    """
    
    content = []
    
    for requirement in requirements :
        content.append("#include \"%s.includes\""%requirement)
    
    content.append("")
    
    includes = set()
    for class_name, instantiations, pointer in classes_template_info :
        for instantiation in instantiations :
            for parameter in instantiation :
                includes.add(parameter[0])
    includes = ["#include <%s.h>"%get_swig_class_name(x) for x in includes]
    content.extend(includes)
    
    for class_name, instantiations, pointer in classes_template_info :
        content.append("#include \"%s.h\""%get_swig_class_name(class_name))
    
    content.append("")
    
    for class_name, instantiations, pointer in classes_template_info :
        pointer = pointer or ""
        swig_class_name = get_swig_class_name(class_name)
        
        template_parameters = []
        for instantiation in instantiations :
            parameter_mangled_names = []
            for entry in instantiation :
                mangled_template_parameters = get_mangled_template_name(entry[0], *entry[1:])
                parameter_mangled_names.append(mangled_template_parameters)
            template_parameters.append(parameter_mangled_names)
        
        for parameter_set in template_parameters :
            cpp_name = [p[0] for p in parameter_set]
            cpp_name = ", ".join(cpp_name)
            cpp_name = "%s< %s >"%(class_name, cpp_name)
            
            mangled_name = [p[1] for p in parameter_set]
            swig_name = "%s%s"%(swig_class_name, "".join(mangled_name))
            
            content.append("typedef %s %s;"%(cpp_name, swig_name))
            if "pointer" in pointer :
                content.append("typedef %s::Pointer::SmartPointer %s_Pointer;"%(cpp_name, swig_name))
            if "superclass" in pointer :
                content.append("typedef %s::Superclass::Self %s_Superclass;"%(cpp_name, swig_name))
            if pointer == "pointer_with_superclass" :
                content.append("typedef %s::Superclass::Pointer::SmartPointer %s_Superclass_Pointer;"%(cpp_name, swig_name))
    
    return "\n".join(content)

def get_module_config_templates(requirements, classes_template_info):
    content = []
    content.append("depends = ('ITKPyBase', %s,)"%(",".join(["'%s'"%x for x in requirements])))
    content.append("templates = (")
    
    for class_name, instantiations, pointer in classes_template_info :
        
        template_parameters = []
        for template_parameter in instantiations :
            parameter_mangled_names = []
            for entry in template_parameter :
                mangled_template_parameters = get_mangled_template_name(entry[0], *entry[1:])
                parameter_mangled_names.append(mangled_template_parameters)
            template_parameters.append(parameter_mangled_names)
        
        unqualified_class_name = get_unqualified_class_name(class_name)
        swig_class_name = get_swig_class_name(class_name)
        
        for parameter_set in template_parameters :
            cpp_name = [p[0] for p in parameter_set]
            cpp_name = ", ".join(cpp_name)
            
            mangled_name = [p[1] for p in parameter_set]
            swig_name = "%s%s"%(swig_class_name, "".join(mangled_name))
            
            content.append("    ('%(unqualified_class_name)s', '%(class_name)s', "
                           "'%(swig_name)s', '%(cpp_name)s'),"%locals())

    content.append(")")
    
    return "\n".join(content)
