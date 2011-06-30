##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import glob
import imp
import os.path
import re
import subprocess
import sys
import traceback

from SCons.Script.SConscript import SConsEnvironment
from SCons.Builder import Builder

import openmp
from utils import merge_construction_variables

if sys.version_info[0] == 2 and sys.version_info[1] <= 5 :
    
    import itertools
    
    def itertools_product(*args, **kwds):
        # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
        # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
        pools = map(tuple, args) * kwds.get('repeat', 1)
        result = [[]]
        for pool in pools:
            result = [x+[y] for x in result for y in pool]
        for prod in result:
            yield tuple(prod)
    
    itertools.product = itertools_product

class Environment(SConsEnvironment) :
    def __init__(self, *args, **kwargs) :
        SConsEnvironment.__init__(self, *args, **kwargs)
        
        known_construction_variables = [
            "ASFLAGS", "CCFLAGS", "CFLAGS", "CPPDEFINES", "CPPPATH",
            "FRAMEWORKPATH", "FRAMEWORKS", "LIBPATH", "LIBS", 
            "LINKFLAGS", "RPATH", "SWIGFLAGS"
        ]

        for construction_variable in known_construction_variables :
            if construction_variable not in self :
                self[construction_variable] = []
            elif str(self[construction_variable]) :
                self[construction_variable] = [str(self[construction_variable])]
            else :
                self[construction_variable] = []
        
        if not self["verbose"] :
            self['CCCOMSTR'] = "Compiling $TARGET"
            self['CXXCOMSTR'] = "Compiling $TARGET"
            self['LINKCOMSTR'] = "Linking $TARGET"
            self["SHCCCOMSTR"] = "Compiling $TARGET"
            self["SHCXXCOMSTR"] = "Compiling $TARGET"
            self['SHLINKCOMSTR'] = "Linking $TARGET"
            self['SWIGCOMSTR'] = "Running SWIG for $TARGET"
        
        if sys.platform == "win32":
            self._msvc_win32()
        elif sys.platform == "linux2":
            self._gcc_linux2()
        
        self._configuration_variables = {}
        
        self._vtk_wrap_python = None
        self._vtk_wrap_python_init = None
        
        sys.path.append(os.path.dirname(__file__))
        
        self._addon_modules = {}
        for name in ["cython", "itk", "python", "swig", "vtk"] :
            try :
                f = None
                f, pathname, description = imp.find_module(
                      name, [os.path.dirname(__file__)]
                      )
                module = imp.load_module(name, f, pathname, description)
                
                if getattr(module, "exists")(self) :
                    getattr(module, "generate")(self)
                    self._addon_modules[name] = module
            except Exception, e:
                if f is not None :
                    f.close()
                print "Warning: could not enable module %s, exception caught."%name
                exc_info = sys.exc_info()
                print "".join(traceback.format_exception(*exc_info))
        
        self._find_itk()
        self._find_openmp()
        self._find_python()
        self._find_vtk()

    
    #############
    # Libraries #
    #############
    
    def HasITK(self):
        """ Test if ITK is usable
        """
        return len(self._configuration_variables["itk"])>0
    
    def UseITK(self):
        """ Add the configuration variables for ITK to the environment
        """
        for key, value in self._configuration_variables["itk"].items() :
            for entry in value :
                if entry not in self[key] :
                    self[key].append(entry)
    
    def UseOpenMP(self):
        """ Add the configuration variables for OpenMP to the environment
        """
        for key, value in self._configuration_variables["openmp"].items() :
            for entry in value :
                if entry not in self[key] :
                    self[key].append(entry)
    
    def UseVTK(self) : 
        """ Add the configuration variables for VTK to the environment
        """
        for key, value in self._configuration_variables["vtk"].items() :
            if key in self :
                for entry in value :
                    if entry not in self[key] :
                        self[key].append(entry)
    
    def HasVTK(self):
        """ Test if VTK is usable
        """
        return len(self._configuration_variables["vtk"])>0

    #####################
    # Private interface #
    #####################

    def _find_itk(self) :
        variables = self._addon_modules["itk"].configuration_variables()
        if variables is None : 
            variables = {}
        self._configuration_variables["itk"] = variables
    
    def _find_openmp(self):
        self._configuration_variables["openmp"] = openmp.configuration_variables()
    
    def _find_python(self) :
        self._configuration_variables["python"] = self._addon_modules["python"].configuration_variables()
    
    def _find_vtk(self) :
        variables = self._addon_modules["vtk"].configuration_variables()
        if variables is None : 
            variables = {}
        else :
            self._vtk_wrap_python = variables["vtk_wrap_python"]
            self._vtk_wrap_python_init = variables["vtk_wrap_python_init"]
            del variables["vtk_wrap_python"]
            del variables["vtk_wrap_python_init"]
        
        self._configuration_variables["vtk"] = variables
    
    def _gcc_linux2(self):
        """ Configure the environment with gcc on linux
        """
        
        if self["optimized"] :
            self["CCFLAGS"].extend(["-O3", "-Wall", "-ftree-vectorize", 
                                    "-DNDEBUG", 
                                    "-Werror=declaration-after-statement"])
        else :
            self["CCFLAGS"].extend(["-g", "-Wall",
                                    "-Werror=declaration-after-statement"])
    
    def _msvc_win32(self):
        """ Configure the environment with Visual Studio on Windows 32 bits
        """
        
        if self['optimized']:
            self["CCFLAGS"].extend(["/EHsc", "/Ox", "/MD", "-D_CRT_SECURE_NO_DEPRECATE", "-DWIN32"])
        else :
            self["CCFLAGS"].extend(["/EHsc", "/GS", "/Wall", "/MDd", "-DWIN32"])

        self["LIBPATH"].extend(os.environ.get("LIB", "").split(';'))
        self["CPPPATH"].extend(os.environ.get("INCLUDE", "").split(';'))

        self["SWIGFLAGS"].extend(['-DWIN32'])
        self["CCFLAGS"].extend(['-DWIN32'])
