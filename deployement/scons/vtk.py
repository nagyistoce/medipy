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
import subprocess
import sys
import tempfile

from SCons.Builder import Builder

from utils import split_environment_variable

def configuration_variables() :
    
    command = ["cmake", "-P", os.path.join(os.path.dirname(__file__), "vtk.cmake")]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode != 0 :
        print "Could not find VTK, some modules will not be compiled."
        print stderr
        return None
    else :
        result = {}
        # Result of cmake -P is on stderr
        for line in stderr.split("\n") :
            if not line :
                continue
            key, value = line.split("=", 1)
            if key in ["CPPPATH", "LIBPATH", "LIBS", "PYTHON_LIBS"] :
                value = list(set(value.split(";")))
            if key in ["LIBS", "PYTHON_LIBS"] :
                value = [x for x in value if x != "general"]
            if key in ["vtk_wrap_python", "vtk_wrap_python_init"] :
                if sys.platform == "win32" and " " in value :
                    value = "\"{0}\"".format(value)
            result[key] = value
        # Remove duplicate entries between LIBS and PYTHON_LIBS
        result["PYTHON_LIBS"] = [x for x in result["PYTHON_LIBS"] if x not in result["LIBS"]]
        # Remove absolute file paths from LIBS
        result["LIBS"] = [x for x in result["LIBS"] if os.path.sep not in x]
        return result

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
