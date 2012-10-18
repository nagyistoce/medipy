##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

""" Builders for ITK and WrapITK.
    
    The "ITK" dictionary will be added to the environment, with the following
    keys :
    * "CPPPATH"
    * "LIBPATH"
    * "LIBS"
    * "SWIGFLAGS"
    
    This module will add the following builders to the environment :
    * CableWrapper : build a CableSwig wrapper
    * ...
    * WrapITKPythonModule : 
      
"""

import os.path
import re
import string
import subprocess
import sys

import wrapitk
from utils import merge_construction_variables, split_environment_variable

import SCons
import SCons.Tool
import SCons.Util
from SCons.Builder import Builder

def itk_root():
    """ Return the root of the ITK installation
    """
    
    result = None

    for dir in _possible_itk_dirs() :
        if os.path.isfile(os.path.join(dir, "ITKConfig.cmake")) :
            result = os.path.dirname(os.path.dirname(dir))
            break
    
    return result

def wrapitk_root():
    """ Return the root of the WrapITK installation
    """
    
    result = None
    
    for dir in _possible_itk_dirs() :
        if os.path.isfile(os.path.join(dir, "WrapITK", "WrapITKConfig.cmake")) :
            result = os.path.join(dir, "WrapITK")
            break
    
    return result

def configuration_variables() :
    """ Return a dictionary containing the SCons configuration variables for
        ITK
    """
    
    possible_itk_dirs = _possible_itk_dirs()
    
    itk_variables = _itk(possible_itk_dirs)
    
    result = None
    if itk_variables is None : 
        print "Could not find ITK, some modules will not be compiled."
        print "Please add the directory containing ITKConfig.cmake to a standard path "
        print "or to the ITK_DIR environment variable"
    else : 
        wrap_itk_variables = _wrap_itk(possible_itk_dirs)
        if wrap_itk_variables is None :
            print "Could not find WrapITK, some modules will not be compiled."
            print "Please add the directory containing WrapITKConfig.cmake to a standard path "
            print "or to the ITK_DIR environment variable"
        else :
            result = merge_construction_variables(itk_variables, wrap_itk_variables)
    
    return result
    
def _possible_itk_dirs() :
    possible_itk_dirs = []
    # Look in some environment variables
    for element in split_environment_variable("PATH") : 
        possible_dir = os.path.join(os.path.dirname(element), "lib", "InsightToolkit")
        if possible_dir not in possible_itk_dirs :
            possible_itk_dirs.append(possible_dir)
    for element in split_environment_variable("ITK_DIR") :
        if element not in possible_itk_dirs : 
            possible_itk_dirs.append(element)
    # Look in predefined paths
    possible_itk_dirs.append(os.path.join("/usr/local/lib/InsightToolkit"))
    possible_itk_dirs.append(os.path.join("/usr/lib/InsightToolkit"))
    # Look in Windows' registry
    if sys.platform == "win32" :
        import _winreg
        key = "Software\\Kitware\\CMakeSetup\\Settings\\StartPath"
        for build_id in range(1, 11) :
            handle = _winreg.OpenKey(_winreg.HKEY_CURRENT_USER, key)
            try :
                value, key_type = _winreg.QueryValueEx(handle, "WhereBuild" + str(build_id))
            except WindowsError, e :
                pass
            else :
                if value not in possible_itk_dirs :
                    possible_itk_dirs.append(value)
    return possible_itk_dirs

def _itk(possible_itk_dirs) : 
    # Find ITKConfig.cmake

    result = None

    for dir in possible_itk_dirs :
        if os.path.isfile(os.path.join(dir, "ITKConfig.cmake")) :
            result = _parse_itk_config(os.path.join(dir, "ITKConfig.cmake"))
            break

    if result is not None :
        result["SWIGFLAGS"] = ["-I%s"%x for x in result["CPPPATH"]]
    
    return result

def _wrap_itk(possible_itk_dirs) :
    # Find WrapITKConfig.cmake
    wrap_itk_dir = None
    for dir in possible_itk_dirs :
        if os.path.isfile(os.path.join(dir, "WrapITK", "WrapITKConfig.cmake")) :
            wrap_itk_dir = os.path.join(dir, "WrapITK")
            break
    
    result = None
    
    if wrap_itk_dir is not None :
        result = {
            "CPPPATH" : [os.path.join(wrap_itk_dir, "Configuration", "Typedefs")],
            "SWIGFLAGS" : [],
        }
        # SWIG flags
        for dir in ["Languages", "Typedefs", os.path.join("Typedefs", "python")] :
            include_dir = os.path.join(wrap_itk_dir, "Configuration", dir)
            result["SWIGFLAGS"].append("-I%s"%include_dir)
    
    return result

def _walk_dependencies(library, dependencies) :
    """ Recursively walk the dependencies, given as a dictionnary.
    """
    
    result = []
    
    queue = [library]
    while queue :
        lib = queue.pop()
        for d in dependencies.get(lib, []) :
            queue.append(d)
        result.append(lib)
    
    # Make each name unique 
    result = list(set(result))
    
    return result
    

def _parse_itk_library_depends(filename, libs) :
    """ Parse an ITKLibraryDepends.cmake file and return a list of libraries.
    """
    
    fd = open(filename, "r")
    itk_library_depends = fd.readlines()
    fd.close()

    itk_library_depends = [x.strip() for x in itk_library_depends]

    # e.g. SET("ITKBasicFilters_LIB_DEPENDS" "ITKCommon;")
    dependency_pattern = r"SET\(\"(\w+)_LIB_DEPENDS\" \"([\w|;|-]+)\"\)"

    dependencies = {}
    for line in itk_library_depends :
        match = re.match(dependency_pattern, line)
        if match :
            name, value = match.groups()
            # Two syntaxes are possible : either a "general" token precedes all
            # libraries, or system libraries are given by "-lpthread" or "-lm". 
            value = [x.split("-l")[-1] for x in value.split(";") if x != "general" and x != ""]
            dependencies[name] = value
    
    result = []
    for lib in libs :
        result.extend(_walk_dependencies(lib, dependencies))
    
    # Make each name unique
    result = list(set(result))
    
    return result

def _parse_itk_config(filename):
    """ Parse an ITKConfig.cmake file and return a dictionary of CPPPATH,
        LIBPATH, CCFLAGS, CXXFLAGS, LIBS.
    """
    
    fd = open(filename, "r")
    itk_config = fd.readlines()
    fd.close()
    
    install_prefix = os.path.dirname(os.path.dirname(os.path.dirname(filename)))
    
    install_prefix_pattern = re.compile(r"\$\{ITK_INSTALL_PREFIX\}")
    
    patterns = [
        r"^SET\(ITK_(INCLUDE_DIRS) \"(.*)\"\)$", 
        r"^SET\(ITK_(LIBRARY_DIRS) \"(.*)\"\)$", 
        r"^SET\(ITK_(REQUIRED_C_FLAGS) \"(.*)\"\)$", 
        r"^SET\(ITK_(REQUIRED_CXX_FLAGS) \"(.*)\"\)$",
        r"^SET\(ITK_(LIBRARY_DEPENDS_FILE) \"(.*)\"\)$",
        r"^SET\(ITK_(LIBRARIES) (.*)\)$",
    ]
    
    patterns = [
        r"SET\(ITK_(\w+) \"(.*)\"\)",
        r"SET\(ITK_(\w+) ([^\"].*)\)",
    ]
    
    cmake_to_scons = {
        "INCLUDE_DIRS" : "CPPPATH", "LIBRARY_DIRS" : "LIBPATH", 
        "REQUIRED_C_FLAGS" : "CCFLAGS", "REQUIRED_CXX_FLAGS" : "CXXFLAGS", 
        "REQUIRED_LINK_FLAGS" : "LINKFLAGS", "LIBRARIES" : "LIBS"
    }
    
    itk_config = [x.strip() for x in itk_config]
    
    result = {}
    itk_library_depends = None
    for line in itk_config :
        if not line :
            continue
        if line.startswith("#") :
            continue
        
        for pattern in patterns :
            match = re.match(pattern, line)
            if match is None :
                continue
            
            name, value  = match.groups()
            
            value = re.sub(install_prefix_pattern, install_prefix, value)
            value = value.strip()
            
            if name == "LIBRARY_DEPENDS_FILE" :
                itk_library_depends = value
            elif name in cmake_to_scons :
                if name in ["INCLUDE_DIRS", "LIBRARY_DIRS"] :
                    value = [ x for x in value.split(";") if x]
                elif name in ["REQUIRED_C_FLAGS", "REQUIRED_CXX_FLAGS", "REQUIRED_LINK_FLAGS", "LIBRARIES"] :
                    value = [ x for x in value.split(" ") if x]
                result[cmake_to_scons[name]] = value
    
    if itk_library_depends is not None :
        result["LIBS"] = _parse_itk_library_depends(itk_library_depends, result["LIBS"])
    
    return result

from wrapitk.utils import get_swig_class_name

def module_builder(env, module_name, classes_template_info, **kwargs) :
    """ Build a WrapITK python module which can be loaded using 
        itkBase.LoadModule
        
        classes_template_info is a list of triplets :
            * class_name
            * instantiations, which is a list (one item per template 
                parameter) of (class_name, template_parameters)
            * use of pointer in typedefs (cf. class_cable_swig_files_builder for
              semantics)
        
        For example, we want to build the class itk::NumpyBridge. This class
        has one template parameter, and we want to instantiate it for 
        itk::Image and itk::VectorImage with the types "unsigned short", 
        "float" and "unsigned char", and the dimensions 2 and 3. We'll have :
        types = ["unsigned short", "float", "unsigned char"]
        dimensions = [2,3]
        instantiations = [[ ("Image", (types, dimensions)), ("VectorImage", (types, dimensions))]]
        
        Note the nested list : instantiations has one item per template 
        parameter, each of these items being the list of instantiations.
        
        For each C++ class called ns::ClassName, generate :
            * a C++ file with a set of typedefs to be instantiated (a Cable wrapper), wrap_nsClassName.cxx
            * a XML representation of the Cable wrapper, wrap_nsClassName.xml
            * an index of the names in the XML representation, wrap_nsClassName.idx
            * a Swig file to wrap that class, wrap_nsClassName.i
            * a compiled file, wrap_nsClassNamePython.o
            * a Python file, nsClassNamePython.py
        
        For the whole module called Module, containing several classes, generate :
            * a Swig file to wrap the module, ModuleName.i
            * a Swig file with extra stuff, ModuleName_ext.i
            * a compiled file, ModulePython.o
            * a Python file, ModulePython.py
            * a shared library _ModulePython.so (from all .o files)
    """
    
    # Use Instantiation if user still uses old framework
    normalized_classes_template_info = []
    for name, instantiations, pointer in classes_template_info :
        template_parameters_list = []
        for instantiation in instantiations :
            template_parameters = []
            for parameters in instantiation :
                if isinstance(parameters, (list, tuple)) and parameters[0] == "itk::Image" :
                    normalized_parameters = wrapitk.utils.Instantiation(
                        parameters[0], *parameters[1:])
                    template_parameters.append(normalized_parameters)
                else: 
                    template_parameters.append(parameters)
            template_parameters_list.append(template_parameters)
        normalized_classes_template_info.append((name, template_parameters_list, pointer))
    classes_template_info = normalized_classes_template_info
    
    suffixes = {
        "header" : [".h", ".H", ".hxx", ".hpp", ".hh"],
        "user_swig" : [".h", ".H", ".hxx", ".hpp", ".hh"],
        "template" : [".txx"]
    }

    module_nodes = wrapitk.module_builders.module_files(
        env, module_name, classes_template_info, ["Base"])
    
    class_nodes = {}
    for name, instantiations, pointer in classes_template_info :
        class_nodes[name] = {}
        swig_name = wrapitk.utils.get_swig_class_name(name)
        
        files = env.File(env.Glob("%s.*"%swig_name))
        
        for suffix_type in suffixes :
            for file in files :
                if os.path.splitext(file.path)[1] in suffixes[suffix_type] :
                    class_nodes[name][suffix_type] = file
        
        nodes = wrapitk.class_builders.cable_swig_files(
            env, name, class_nodes[name]["header"], instantiations, pointer)
        
        for n,v in nodes.items() :
            class_nodes[name][n] = v
        
        if "template" in class_nodes[name] :
            env.Depends(class_nodes[name]["class_xml"], class_nodes[name]["template"])
    
    module_nodes["master_index"] = Builder(action=wrapitk.module_builders.master_index)(
        env, "{0}.mdx".format(module_name), 
        [nodes["class_index"] for nodes in class_nodes.values()])
    
    for name, template_parameters_list, pointer in classes_template_info :
        swig_nodes = wrapitk.class_builders.swig_files(
            env, name, template_parameters_list, pointer, class_nodes[name], module_name, module_nodes, 
            wrapitk_root())
        
#        doc_nodes =wrapitk.class_builders.class_doc_files_builder(
#            env, name, class_nodes[name]["header"], template_parameters_list, 
#            wrapitk_root())
        
        for n,v in swig_nodes.items() :
            class_nodes[name][n] = v
#        for n,v in doc_nodes.items() :
#            class_nodes[name][n] = v

    sources = []
    sources.extend(module_nodes["module_swig"])
    for name, instantiations, pointer in classes_template_info :
        sources.extend(class_nodes[name]["class_swig"])
    
    swig_flags = ["-w%i"%x for x in [508,312,314,509,302,362,389,384,383,361,467]]
    swig_flags += ["-O", "-features", "autodoc=1", "-Werror"]
    
    shared_objects = env.SharedObject(
        source=sources,
        CPPDEFINES=env["CPPDEFINES"]+["SWIG_PYTHON_CAST_MODE"],
        CPPPATH=env["CPPPATH"]+env["Python"]["CPPPATH"]+env["ITK"]["CPPPATH"], 
        SWIGFLAGS=env["SWIGFLAGS"]+env["Swig"]["SWIGFLAGS"]+env["ITK"]["SWIGFLAGS"]+swig_flags,
    )
    
    if sys.platform != "win32" :
        command = ["gcc", "-dumpversion"]
        process = subprocess.Popen(command, stdout=subprocess.PIPE)
        stdout, _ = process.communicate()
        elements = stdout.split(".")
        major = int(elements[0])
        minor = int(elements[1])
        if major > 4 or (major == 4 and minor > 4) :
            env.Replace(CC="gcc-4.4", CXX="g++-4.4")
    
    env.PythonModule("_%sPython"%module_name, shared_objects,
                     LIBPATH=env["LIBPATH"]+env["ITK"]["LIBPATH"],
                     LIBS=env["LIBS"]+env["ITK"]["LIBS"])

# SCons.Tool.swig._find_module from scons 1.2.0.d20100117, used for instead of
# the buggy version of scons < 1.2.0

# Match '%module test', as well as '%module(directors="1") test'
# Also allow for test to be quoted (SWIG permits double quotes, but not single)
patched_reModule = re.compile(r'%module(\s*\(.*\))?\s+("?)(.+)\2')
def patched_find_modules(src):
    """Find all modules referenced by %module lines in `src`, a SWIG .i file.
       Returns a list of all modules, and a flag set if SWIG directors have
       been requested (SWIG will generate an additional header file in this
       case.)"""
    directors = 0
    mnames = []
    try:
        matches = patched_reModule.findall(open(src).read())
    except IOError:
        # If the file's not yet generated, guess the module name from the filename
        matches = []
        mnames.append(os.path.splitext(src)[0])

    for m in matches:
        mnames.append(m[2])
        directors = directors or string.find(m[0], 'directors') >= 0
    return mnames, directors    

# Same as patched_find_modules
def patched_add_director_header_targets(target, env):
    # Directors only work with C++ code, not C
    suffix = env.subst(env['SWIGCXXFILESUFFIX'])
    # For each file ending in SWIGCXXFILESUFFIX, add a new target director
    # header by replacing the ending with SWIGDIRECTORSUFFIX.
    for x in target[:]:
        n = x.name
        d = x.dir
        if n[-len(suffix):] == suffix:
            target.append(d.File(n[:-len(suffix)] + env['SWIGDIRECTORSUFFIX']))

# Same as patched_find_modules
def patched_swigEmitter(target, source, env):
    swigflags = env.subst("$SWIGFLAGS", target=target, source=source)
    flags = SCons.Util.CLVar(swigflags)
    for src in source:
        src = str(src.rfile())
        mnames = None
        if "-python" in flags and "-noproxy" not in flags:
            if mnames is None:
                mnames, directors = patched_find_modules(src)
            if directors:
                patched_add_director_header_targets(target, env)
            python_files = map(lambda m: m + ".py", mnames)
            outdir = env.subst('$SWIGOUTDIR', target=target, source=source)
            # .py files should be generated in SWIGOUTDIR if specified,
            # otherwise in the same directory as the target
            if outdir:
                python_files = map(lambda j, o=outdir, e=env:
                                   e.fs.File(os.path.join(o, j)),
                                   python_files)
            else:
                python_files = map(lambda m, d=target[0].dir:
                                   d.File(m), python_files)
            target.extend(python_files)
        if "-java" in flags:
            if mnames is None:
                mnames, directors = patched_find_modules(src)
            if directors:
                patched_add_director_header_targets(target, env)
            java_files = map(lambda m: [m + ".java", m + "JNI.java"], mnames)
            java_files = SCons.Util.flatten(java_files)
            outdir = env.subst('$SWIGOUTDIR', target=target, source=source)
            if outdir:
                java_files = map(lambda j, o=outdir: os.path.join(o, j), java_files)
            java_files = map(env.fs.File, java_files)
            for jf in java_files:
                t_from_s = lambda t, p, s, x: t.dir
                SCons.Util.AddMethod(jf, t_from_s, 'target_from_source')
            target.extend(java_files)
    return (target, source)

def exists(env):
    return (configuration_variables() is not None)

def generate(env):
    env["ITK"] = configuration_variables()
    
    # Patch the Swig emitter for old versions of SCons. cf. release notes
    # of 1.2.0.d20100117
    version = SCons.__version__.split(".")
    if int(version[0]) <= 1 and int(version[1]) <= 2 and version[-1] != "d20100117" :
        c_file, cxx_file = SCons.Tool.createCFileBuilders(env)
        c_file.emitter[".i"] = patched_swigEmitter
        cxx_file.emitter[".i"] = patched_swigEmitter 
    
    env.AddMethod(module_builder, "WrapITKPythonModule") 
