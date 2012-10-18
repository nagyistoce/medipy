from SCons.Builder import Builder

import class_wrapper
import utils

def cable_swig_files(env, class_name, header, instantiations, pointer=None):
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
    
    swig_class_name = utils.get_swig_class_name(class_name)
    
    includes = set(env["ITK"]["CPPPATH"] + env["Python"]["CPPPATH"] + env["CPPPATH"])
    includes = [x.replace("\\", "/") for x in includes]
    includes = [x.rstrip("/") for x in includes]
    include_directives = " ".join(["-I%s"%x for x in includes])
    
    class_xml = ("gccxml -fxml-start=_cable_ -fxml=$TARGET -DCSWIG " 
                 "-DCABLE_CONFIGURATION -DITK_MANUAL_INSTANTIATION "
                 "$INCLUDE_DIRECTIVES $SOURCE")
    
    builders = {
        "wrapper" : Builder(action=cable_swig),
        "xml" : Builder(action=class_xml),
        "index" : Builder(action="cableidx $SOURCE $TARGET"),
    }
    
    targets = {
        "wrapper" : "wrap_%s.cxx"%swig_class_name,
        "xml" : "wrap_%s.xml"%swig_class_name,
        "index" : "wrap_%s.idx"%swig_class_name,
    }
    
    extra_arguments = {
        "wrapper" : { "INSTANTIATIONS" : instantiations,
                      "CLASS_NAME" : class_name,
                      "POINTER" : pointer }, 
        "xml" : { "INCLUDE_DIRECTIVES" : include_directives },
        "index" : { },
    } 
    
    nodes = {}
    
    source = header
    # Build order is important, cannot use dictionary iterator
    for name in ["wrapper", "xml", "index"] :
        builder = builders[name]
        node = builder(env, targets[name], source, **extra_arguments[name])
        nodes["class_{0}".format(name)] = node 
        source = node
    
    return nodes

def swig_files(env, class_name, instantiations, pointer, class_nodes, 
                             module_name, module_nodes, wrapitk_root) :
    """ Builder for class-specific Swig files : wrap_MyClass.i and 
        wrap_MyClass_ext.i 
    """
    
    swig_class_name = utils.get_swig_class_name(class_name)
    
    builders = {
        "swig" : Builder(generator=swig),
        "swig_extra" : Builder(action=swig_extra),
    }
    
    targets = {
        "swig" : "wrap_%s.i"%swig_class_name,
        "swig_extra" : "wrap_%s_ext.i"%swig_class_name,
    }
    
    extra_arguments = {
        "swig" : { "CLASS_NAME" : class_name,
                         "MODULE_NAME" : module_name, 
                         "MASTER_INDEX" : module_nodes["master_index"][0],
                         "WRAPITK_ROOT" : wrapitk_root }, 
        "swig_extra" : { "CLASS_NAME" : class_name,
                               "INSTANTIATIONS" : instantiations,
                               "POINTER" : pointer },
    } 
    
    nodes = {}
    
    for name, builder in builders.items() :
        node = builder(env, targets[name], class_nodes["class_xml"], **extra_arguments[name])
        nodes["class_{0}".format(name)] = node 
    
    env.Depends(nodes["class_swig"], 
                [module_nodes["module_includes"], module_nodes["master_index"][0],
                 class_nodes["user_swig"], nodes["class_swig_extra"]])
    
    return nodes

def cable_swig(target, source, env) :
    """ cf. class_cable_swig_files_builder for doc about pointer
    """
    
    content = class_wrapper.get_cable_swig(
        env["CLASS_NAME"], env["INSTANTIATIONS"], env["POINTER"])
    file = open(str(target[0]), "w")
    file.write(content)
    file.close()

def swig(source, target, env, for_signature) : 
    # This is a generator, hence the different prototype
    
    swig_class_name = utils.get_swig_class_name(env["CLASS_NAME"])
    
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

def swig_extra(target, source, env) :
    
    content = class_wrapper.get_swig_extra(
        env["CLASS_NAME"], env["INSTANTIATIONS"], env["POINTER"])
    file = open(str(target[0]), "w")
    file.write(content)
    file.close()
