##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

""" Builders for Python.
    
    The "Python" dictionary will be added to the environment, with the following
    keys :
    * "CPPPATH"
    * "LIBPATH"
    * "LIBS"
    
    This module will add the following builders to the environment :
    * Python : compile a Python source file to bytecode
    * PythonModule : create a Python module from Swig and source files
      
"""

import distutils.sysconfig
import os.path
import py_compile
import re
import sys

import numpy.distutils

from SCons.Builder import Builder

def configuration_variables() :
    result = {
        "CPPPATH" : [distutils.sysconfig.get_python_inc()] + numpy.distutils.misc_util.get_numpy_include_dirs(),
        "LIBPATH" : [],
        "LIBS" : [],
    }
    
    version_with_dots = ".".join([str(x) for x in sys.version_info[:2]])
    version_without_dots = "".join([str(x) for x in sys.version_info[:2]])
    if sys.platform == "linux2" :
        result["LIBS"] = ["python%s"%version_with_dots, "boost_python-mt-py%s"%version_without_dots]
    elif sys.platform == "win32" :
        # No dot in the library on windows
        result["LIBS"] = ["python%s"%version_without_dots]
        # Add Python libraries to the LIBPATH
        python_lib = distutils.sysconfig.get_python_lib()
        path = os.path.join(python_lib, "..", "..", "libs")
        path = os.path.normpath(path)
        result["LIBPATH"] = [path]
    
    return result

def python_command(target, source, env):
    """ Generate Python bytecode from source file
    """

    exit_code = 0

    for i, s in enumerate(source) :
        try :
            py_compile.compile(s.abspath, target[i].abspath)
        except Exception, e:
            print e
            exit_code = -1

    return exit_code

def python_emitter(target, source, env):
    target=[]
    for file in source:
        target.append(str(file) + "c")

    return target, source 

def module_builder(self, module_name, sources, **kwargs) :
    """ Build a Python module from the given sources
    """
    
    # Find the SWIG files in the sources, replace them by the wrappers
    swig_pattern = re.compile(r"(.*).i")
    swig_sources = []
    for source in sources :
        c_file = swig_pattern.sub("\1_wrap.c", str(source))
        cpp_file= swig_pattern.sub("\1_wrap.cc", str(source))
        if os.path.exists(c_file) :
            swig_sources.append(c_file)
        elif os.path.exists(cpp_file) :
            swig_sources.append(cpp_file)
        else : 
            swig_sources.append(source)
    env = self.Clone()

    # Add the necessary flags
    kwargs["SWIGFLAGS"] = kwargs.get("SWIGFLAGS", [])+env["Swig"]["SWIGFLAGS"]+env["SWIGFLAGS"]
    kwargs["SHLIBPREFIX"] = env["Swig"]["SHLIBPREFIX"]
    if "SHLIBSUFFIX" in env["Swig"] :
        kwargs["SHLIBSUFFIX"] = env["Swig"]["SHLIBSUFFIX"]
    
    for key in ["CPPPATH", "LIBPATH", "LIBS"] :
        kwargs[key] = kwargs.get(key, [])+env["Python"][key]+env[key]
        
    
    shared_library = env.SharedLibrary(module_name, swig_sources, **kwargs)
    if sys.platform == "win32" :
        env.AddPostAction(shared_library, 'mt.exe -nologo -manifest ${TARGET}.manifest -outputresource:$TARGET;2')

def exists(env):
    return (configuration_variables() is not None)

def generate(env):
    env["Python"] = configuration_variables()
    
    env["BUILDERS"]["Python"] = Builder(action=python_command,
                                        suffix=".pyc", src_suffix=".py",
                                        emitter=python_emitter)
    
    env.AddMethod(module_builder, "PythonModule") 