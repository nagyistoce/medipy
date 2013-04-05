import class_templates
import utils

def get_cable_swig(class_name, template_parameters_list, pointer):
    # Instantiations includes
    includes = set()
    for template_parameters in template_parameters_list :
        full_type = utils.Instantiation(class_name, *template_parameters)
        includes.update(utils.get_includes(full_type))
    includes = ["#include \"{0}.h\"".format(include) for include in includes]
    
    # Typedefs
    typedefs = []
    for template_parameters in template_parameters_list :
        full_type = utils.Instantiation(class_name, *template_parameters)
        typedefs.extend(
            ["typedef {0} {1};".format(t1, t2)
             for t1, t2 in utils.get_typedefs(full_type, pointer, True)])
    
    includes.insert(0, "#include <itkCommand.h>")
    
    return class_templates.cable_swig.format(includes="\n".join(includes), 
        swig_class_name=utils.get_swig_class_name(class_name), 
        typedefs="\n".join(typedefs), forced_instantiations="", extra="")

def get_swig_extra(class_name, template_parameters_list, pointer):
    ref_count_declarations = []
    
    if pointer == "pointer" :
        for template_parameters in template_parameters_list :
            full_type = utils.Instantiation(class_name, *template_parameters)
            ref_count_declarations.append(
                "DECLARE_REF_COUNT_CLASS({0})".format(full_type.mangled_name))
    
    return class_templates.swig_extra.format(
        swig_class_name=utils.get_swig_class_name(class_name),
        ref_count_declarations="\n".join(ref_count_declarations))
