from SCons.Builder import Builder

import module_wrapper

def module_files(env, module_name, templates, requirements):
    """ Builder for module-specific files : Module.py, ModuleConfig.py, 
        Module.i, Module_ext.i, Module.includes
    """
    
    builders = {
        "loader" : Builder(action=loader),
        "config" : Builder(action=config),
        "swig" : Builder(action=swig),
        "swig_extra" : Builder(action=swig_extra),
        "includes" : Builder(action=includes),
    }
    
    targets = {
        "loader" : "{0}.py".format(module_name),
        "config" : "{0}Config.py".format(module_name),
        "swig" : "_{0}.i".format(module_name),
        "swig_extra" : "_{0}_ext.i".format(module_name),
        "includes" : "{0}.includes".format(module_name), 
    }
    
    extra_arguments = {
        "loader" : { "MODULE_NAME" : module_name},
        "config" : { "REQUIREMENTS" : requirements, "TEMPLATES" : templates},
        "swig" : { "MODULE_NAME" : module_name,  "REQUIREMENTS" : requirements},
        "swig_extra" : { "REQUIREMENTS" : requirements, "TEMPLATES" : templates},
        "includes" : { "REQUIREMENTS" : requirements, "TEMPLATES" : templates},
    } 
    
    nodes = {}
    
    for name, builder in builders.items() :
        node = builder(env, targets[name], [], **extra_arguments[name])
        nodes["module_{0}".format(name)] = node 
    
    env.Depends(nodes["module_swig"], nodes["module_swig_extra"])
    
    return nodes

def loader(target, source, env):
    content = module_wrapper.get_loader(env["MODULE_NAME"])
    file = open(str(target[0]), "w")
    file.write(content)
    file.close()

def config(target, source, env):
    content = module_wrapper.get_config(env["REQUIREMENTS"], env["TEMPLATES"])
    file = open(str(target[0]), "w")
    file.write(content)
    file.close()

def swig(target, source, env) :
    content = module_wrapper.get_swig(env["MODULE_NAME"], env["REQUIREMENTS"])
    file = open(str(target[0]), "w")
    file.write(content)
    file.close()

def swig_extra(target, source, env) :
    content = module_wrapper.get_swig_extra(
        env["REQUIREMENTS"], [x[0] for x in env["TEMPLATES"]])
    file = open(str(target[0]), "w")
    file.write(content)
    file.close()

def includes(target, source, env):
    content = module_wrapper.get_includes(env["REQUIREMENTS"], env["TEMPLATES"])
    file = open(str(target[0]), "w")
    file.write(content)
    file.close()

def master_index(target, source, env):
    content = module_wrapper.get_master_index([x.name for x in source])
    file = open(str(target[0]), "w")
    file.write(content)
    file.close()
    