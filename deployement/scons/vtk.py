##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

""" Builders for VTK and VTK/Python.
    
    The "VTK" dictionary will be added to the environment, with the following
    keys :
    * "CPPPATH"
    * "LIBPATH"
    * "LIBS"
    * "PYTHON_LIBS"
    * "vtk_wrap_python", "vtk_wrap_python_init" : paths to the tools
    
    This module will add the following builders to the environment :
    * VTKClassWrapper : builder for vtkWrapPython
      $SOURCE : declaration of a VTK class (e.g. vtkMyClass.h)
      $TARGET : Python wrapper of the VTK class, defaults to vtkMyClass_wrap.cpp) 
      $HINT_FILE : optional, hints to be passed to vtkWrapPython 
    * VTKDataFile : build the datafile required by vtkWrapPythonInit
      $SOURCE : dummy variable
      $MODULE_NAME
      $CLASSES : classes to store in the data file
      $TARGET : set to $MODULE_NAME.data
    * VTKModuleInit : builder for vtkWrapPythonInit
      $SOURCE : data file
      $TARGET : module init file
    * vtkPythonModule : build a Python module from a set of classes and an 
      optional hint file. 
      
"""

import os.path
import re
import sys

from SCons.Builder import Builder

from utils import split_environment_variable

def configuration_variables() :
    possible_vtk_dirs = _possible_vtk_dirs()
    
    # Find VTKConfig.cmake
    vtk_dir = None

    for dir in possible_vtk_dirs :
        if os.path.isfile(os.path.join(dir, "VTKConfig.cmake")) :
            vtk_dir = dir
            break
    
    result = None
    
    if vtk_dir is None : 
        print "Could not find VTK, some modules will not be compiled."
        print "Please add the directory containing VTKConfig.cmake to a standard path "
        print "or to the VTK_DIR environment variable"
    else :     
        result = { 
                "CPPPATH" : [],
                "LIBPATH" : [],
                "LIBS" : [],
                "PYTHON_LIBS" : []
        }

        vtk_install_prefix = os.path.dirname(os.path.dirname(vtk_dir))
        
        headers_found = False
        for dir in ["vtk-5.0", "vtk-5.2", "vtk"] :
            path = os.path.join(vtk_install_prefix, "include", dir)
            if os.path.isdir(path) :
                result["CPPPATH"].append(path)
                headers_found = True
                break 
        
        if not headers_found : 
            raise Exception("VTK was found, but headers could not be found")
        
        if os.path.isdir(os.path.join(vtk_install_prefix, "lib")) :
            result["LIBPATH"].append(os.path.join(vtk_install_prefix, "lib"))
        for extension in ["", "-5.0", "-5.1", "-5.2", "-5.3", "-5.4"] :
            if os.path.isdir(os.path.join(vtk_install_prefix, "lib", "vtk%s"%extension)) :
                result["LIBPATH"].append(os.path.join(vtk_install_prefix, "lib", "vtk%s"%extension))
        
        for lib in ["Common", "Filtering", "Graphics", "Hybrid", "Imaging", "IO", 
                    "Parallel", "Rendering", "sys", "VolumeRendering", "Widgets"] :
            result["LIBS"].append("vtk"+lib)
            if lib != "sys" :
                result["PYTHON_LIBS"].append("vtk"+lib+"PythonD")

        vtk_wrap_python = os.path.join(vtk_install_prefix, "bin", "vtkWrapPython")
        vtk_wrap_python_init = os.path.join(vtk_install_prefix, "bin", "vtkWrapPythonInit")
        
        result["vtk_wrap_python"] = vtk_wrap_python
        result["vtk_wrap_python_init"] = vtk_wrap_python_init
        
        for command in ["vtk_wrap_python", "vtk_wrap_python_init"] :
            if sys.platform == "win32" and " " in result[command] :
                result[command] = "\"%s\""%result[command]
    
    return result

def _possible_vtk_dirs() :
    possible_vtk_dirs = []
    # Look in some environment variables
    for element in split_environment_variable("PATH") : 
        for extension in ["", "-5.0", "-5.1", "-5.2", "-5.3", "-5.4"] :
            posssible_dir = os.path.join(os.path.dirname(element), "lib", "vtk"+extension)
            if posssible_dir not in possible_vtk_dirs :
                possible_vtk_dirs.append(posssible_dir)
    for element in split_environment_variable("VTK_DIR") : 
        if element not in possible_vtk_dirs :
                possible_vtk_dirs.append(element)
    # Look in predefined paths
    possible_vtk_dirs.append(os.path.join("/usr/local/lib/vtk"))
    possible_vtk_dirs.append(os.path.join("/usr/lib/vtk"))
    # Look in Windows' registry
    if False and sys.platform == "win32" :
        import _winreg
        key = "Software\\Kitware\\CMakeSetup\\Settings\\StartPath"
        for build_id in range(1, 11) :
            handle = _winreg.OpenKey(_winreg.HKEY_CURRENT_USER, key)
            try :
                value, key_type = _winreg.QueryValueEx(handle, "WhereBuild" + str(build_id))
            except WindowsError, e :
                pass
            else :
                if value not in possible_vtk_dirs :
                    possible_vtk_dirs.append(value)
    return possible_vtk_dirs

def datafile_command(target, source, env):
    module_name = env["MODULE_NAME"].split(".")[-1]
    
    data_file = open(target[0].path, "w")
    data_file.write(module_name + os.linesep)
    for klass in env["CLASSES"] : 
        data_file.write(klass + os.linesep)
    data_file.close()
    
def datafile_emitter(target, source, env):
    source_directory = os.path.dirname(
        os.path.abspath(env["MODULE_NAME"].replace(".", os.path.sep)))
    target = [os.path.join(source_directory, 
                           env["MODULE_NAME"].split(".")[-1]+".data")]
    return target, source    

def module_init_command(target, source, env):
    import subprocess
    subprocess.call([env["VTK"]["vtk_wrap_python_init"],
                    source[0].abspath, target[0].abspath])
    
    module_name = env["MODULE_NAME"].split(".")[-1]
    if sys.platform != "win32" :
        path = target[0].abspath
        
        file = open(path, "r")
        lines = file.readlines()
        file.close()
        
        new_lines = [line.replace("lib%s"%module_name, "%s"%module_name)
                     for line in lines]
        
        file = open(path, "w")
        file.writelines(new_lines)
        file.close()

def module_init_emitter(target, source, env):
    source_directory = os.path.dirname(
        os.path.abspath(env["MODULE_NAME"].replace(".", os.path.sep)))
    target = [os.path.join(source_directory, 
                           env["MODULE_NAME"].split(".")[-1]+"_init.cpp")]
    return target, source    

def module_builder(env, module_name, classes, hint_file = ""):
    
    source_directory = os.path.dirname(
        os.path.abspath(module_name.replace(".", os.path.sep)))
    
    if hint_file != "" :
        hint_file = env.FindFile(hint_file, [source_directory]).path
    
    # 1. find files related to classes
    source_suffixes = [".c", ".C", ".cxx", ".cpp", ".c++", ".cc"]
    header_suffixes = [".h", ".H", ".hxx", ".hpp", ".hh"]
    
    files = {}
    for klass in classes :
        files[klass] = {}
        
        for suffix in source_suffixes :
            node = env.FindFile(klass+suffix, [source_directory])
            if node is not None :
                files[klass]["source"] = node.abspath
                
        for suffix in header_suffixes :
            node = env.FindFile(klass+suffix, [source_directory])
            if node is not None :
                files[klass]["header"] = node.abspath

    # 2. Create a wrapper for each class
    class_wrappers = []
    for klass, entries in files.items() :
        nodes = env.VTKClassWrapper(entries["header"], HINT_FILE=hint_file)
        class_wrappers.extend(nodes)
    
    # 3. Data file
    data_file = env.VTKDataFile([x["header"] for x in files.values()], 
                                CLASSES=classes, MODULE_NAME=module_name)

    # 4. Python module init function
    module_init = env.VTKModuleInit(data_file, MODULE_NAME=module_name)
    
    # 5. Python Module
    sources = class_wrappers + module_init 
    for klass in files :
        if "source" in files[klass] :
            sources.append(files[klass]["source"])

    module = env.PythonModule(module_name, sources, 
                              CPPPATH=env["VTK"]["CPPPATH"],
                              LIBPATH = env["VTK"]["LIBPATH"],
                              LIBS = env["VTK"]["LIBS"]+env["VTK"]["PYTHON_LIBS"])
    return module

def exists(env):
    return (configuration_variables() is not None)

def generate(env):
    env["VTK"] = configuration_variables()
    
    class_wrapper_command = "%s $SOURCE $HINT_FILE true $TARGET"%env["VTK"]["vtk_wrap_python"]
    env["BUILDERS"]["VTKClassWrapper"] = Builder(action=class_wrapper_command,
                                                 suffix="_wrap$CXXFILESUFFIX")
    
    env["BUILDERS"]["VTKDataFile"] = Builder(action=datafile_command,
                                             emitter=datafile_emitter)
    
    env["BUILDERS"]["VTKModuleInit"] = Builder(action=module_init_command,
                                               emitter=module_init_emitter) 
    
    env.AddMethod(module_builder, "vtkPythonModule")