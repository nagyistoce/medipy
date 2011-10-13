##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import distutils.core
import glob
import imp
import os
import os.path
import re
import shutil
import subprocess
import sys
import _winreg

import py2exe
import py2exe.mf

import medipy

def get_build_name():
    """ Return the build name, based on the platform and on the SVN revision.
    """
    
    command = ["hg", "tip", "--template", "{rev}"]
    revision = subprocess.Popen(command, stdout=subprocess.PIPE).communicate()[0]

    build_name = "%s-r%s"%(sys.platform, revision)
    
    return build_name

def patch_load_package(self, fqname, pathname):
    m = self.original_load_package(fqname, pathname)
    if fqname.startswith("medipy") :
        # See if there is a wrapitk module
        data = open(m.__file__).readlines()
        for line in data :
            match = re.match(r""".*load_wrapitk_module\(.*, "(.*)"\)""", line)
            if match :
                self.import_hook("{0}.{1}Config".format(fqname, match.group(1)))
                self.import_hook("{0}.{1}Python".format(fqname, match.group(1)))
    return m
py2exe.mf.ModuleFinder.original_load_package = py2exe.mf.ModuleFinder.load_package 
py2exe.mf.ModuleFinder.load_package = patch_load_package

def patch_find_module(self, name, path, parent=None) :
    try :
        result = self.original_find_module(name, path, parent)
    except :
        # Look for a MediPy plugin
        return imp.find_module(
            name, os.environ["MEDIPY_PLUGINS_PATH"].split(os.pathsep))
    else :
        return result
py2exe.mf.ModuleFinder.original_find_module = py2exe.mf.ModuleFinder.find_module
py2exe.mf.ModuleFinder.find_module = patch_find_module

def setup(project_name, main_script, includes=None, medipy_plugins=None):
    
    includes = includes or []
    medipy_plugins = medipy_plugins or []
    key = _winreg.OpenKey(
        _winreg.HKEY_LOCAL_MACHINE, 
        "SOFTWARE\\Microsoft\\VisualStudio\\9.0\\Setup\\VC")
    vs_root, key_type = _winreg.QueryValueEx(key, "ProductDir")
    key.Close()
    
    bin_directory = os.path.join("bin", "%s-%s"%(project_name, get_build_name()))
    
    # Remove destination dir if it exists
    if os.path.exists(bin_directory) :
        shutil.rmtree(bin_directory)
    
    # Process specific includes
    if "itk" in includes :
        import itk
        import itkTypes
        wrapitk_root = os.path.dirname(os.path.dirname(itkTypes.__file__))
        for file in os.listdir(os.path.join(wrapitk_root, "lib")) :
            # Only load "library" modules
            if not file.startswith("_") and not file.startswith("itk") :
                swig_module = os.path.splitext(file)[0]
                if swig_module not in includes :
                    includes.append(swig_module)
        # Copy configuration files from WrapITK
        configuration_dir = os.path.join(wrapitk_root, "Python", "Configuration")
        shutil.copytree(configuration_dir, os.path.join(bin_directory, "Configuration"))
    includes.extend(["medipy.{0}".format(plugin) for plugin in medipy_plugins])
    
    # Include Visual C runtime DLL
    data_files = [("Microsoft.VC90.CRT", 
        glob.glob(os.path.join(vs_root,  
                               "redist", "x86", "Microsoft.VC90.CRT", "*.*")))]
    print data_files

    # Main setup script
    sys.path.append(os.path.join(wrapitk_root, "lib"))
    distutils.core.setup(
        name = project_name,
        windows = [main_script],
        data_files=data_files, 
        options = { 
            "py2exe" : {
                "includes" : includes, 
                "dist_dir" : str(bin_directory), # py2exe does not like unicode strings 
                "verbose" : False,
                "skip_archive" : False,
                "excludes" : ["Tkconstants","Tkinter","tcl"],
                "packages" : ["gzip"],
                "skip_archive" : True
            }
        },
    )
    sys.path.pop()
    
    # Copy resources
    medipy_resources = list(os.walk(os.path.join(medipy.__path__[0], "resources")))
    project_resources = list(os.walk("resources"))
    for dirpath, dirnames, filenames in medipy_resources+project_resources :
        if ".svn" in dirnames :
            del dirnames[dirnames.index(".svn")]
        for filename in filenames :
            skip_file = (filename == "SConstruct" or
                         filename.endswith(".fbp"))
            if not skip_file :
                source = os.path.join(dirpath, filename)
            
                if dirpath.startswith(os.path.join(medipy.__path__[0])) :
                    prefix = dirpath[
                        len(os.path.join(medipy.__path__[0]))+len(os.path.sep):]
                    destination_dir = os.path.join(bin_directory, prefix)
                else :
                    destination_dir = os.path.join(bin_directory, dirpath)
                if not os.path.isdir(destination_dir) :
                    os.makedirs(destination_dir)
                
                destination = os.path.join(destination_dir, filename)
                print "copying {0} -> {1}".format(source, destination)
                shutil.copyfile(source, destination)
    