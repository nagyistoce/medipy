import module_templates
import utils

def get_loader(name):
    return module_templates.loader.format(module_name=name)

def get_config(requirements, classes_template_info):
    content = []
    
    depends = ["ITKPyBase"]+requirements
    content.append("depends = ({0})".format(
        ",".join(["{0}".format(repr(x)) for x in depends])))
    
    content.append("templates = (")
    for class_name, template_parameters_list, pointer in classes_template_info :
        unqualified_class_name = class_name.split("::")[-1]
        for template_parameters in template_parameters_list :
            full_type = utils.Instantiation(class_name, *template_parameters)
            mangled_name = full_type.mangled_name
            cpp_parameters = [x.cpp_name if hasattr(x, "cpp_name") else x
                              for x in template_parameters]
            template_info = "    ('{0}', '{1}', '{2}', '{3}'),".format(
                unqualified_class_name, class_name, 
                mangled_name, ", ".join(cpp_parameters))
            content.append(template_info)
    content.append(")")
    
    return "\n".join(content)

def get_swig(name, requirements):
    includes = "\n".join(["#include \"{0}.includes\"".format(x) 
                          for x in requirements])
    return module_templates.swig.format(
        module_name=name, includes=includes, interface_content="")

def get_swig_extra(requirements, classes):
    
    declarations = []
    calls = []
    for class_ in classes :
        function = "init_{0}Python()".format(utils.get_swig_class_name(class_))
        
        declarations.append("extern \"C\" int {0};".format(function))
        calls.append("{0};".format(function))
    
    imports = ["import {0}Python".format(x) for x in requirements]
    for class_ in classes :
        swig_class_name = utils.get_swig_class_name(class_)
        imports.append("from {0}Python import *".format(swig_class_name))
    imports = "\n".join(imports)
    
    return module_templates.swig_extra.format(
        declarations="\n".join(declarations), calls="\n".join(calls), 
        imports=imports)

def get_includes(requirements, classes_template_info):
    content = []
    
    # WrapITK includes
    for requirement in requirements :
        content.append("#include \"{0}.includes\"".format(requirement))
    content.append("")
    
    # Instantiations includes
    includes = set()
    for class_name, template_parameters_list, pointer in classes_template_info :
        for template_parameters in template_parameters_list :
            full_type = utils.Instantiation(class_name, *template_parameters)
            includes.update(utils.get_includes(full_type))
    includes = ["#include \"{0}.h\"".format(include) for include in includes]
    content.extend(includes)
    content.append("")
    
    # Typedefs
    for class_name, template_parameters_list, pointer in classes_template_info :
        for template_parameters in template_parameters_list :
            full_type = utils.Instantiation(class_name, *template_parameters)
            typedefs = utils.get_typedefs(full_type, pointer)
            typedefs = ["typedef {0} {1};".format(t1, t2) for t1, t2 in typedefs]
            content.extend(typedefs)
    
    return "\n".join(content)

def get_master_index(index_files) :
    return "\n".join(index_files)
    