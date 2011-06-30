##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import distutils.core
import os
import os.path
import shutil
import subprocess
import sys
import xml.dom.minidom

import py2exe

# Duck-punch py2exe DLL inclusion to include msvcp71.dll
origIsSystemDLL = py2exe.build_exe.isSystemDLL
def isSystemDLL(pathname):
        if os.path.basename(pathname).lower() in ("msvcp71.dll", "gdiplus.dll"):
                return 0
        return origIsSystemDLL(pathname)
py2exe.build_exe.isSystemDLL = isSystemDLL

def get_build_name():
    """ Return the build name, based on the platform and on the SVN revision.
    """
    
    command = "svn info --xml --non-interactive".split(" ")
    data = subprocess.Popen(command, stdout=subprocess.PIPE).communicate()[0]
    
    document = xml.dom.minidom.parseString(data)
    entries = document.getElementsByTagName("entry")
    revision = entries[0].getAttribute("revision")

    build_name = "%s-r%s"%(sys.platform, revision)
    return build_name

def get_api_files(root) :
    """ Find all api.py files under given root
    """
    api_files = []
    for dirpath, dirnames, filenames in os.walk(root) :
        if "api.py" in filenames :
            api_files.append(os.path.join(dirpath, "api.py"))
    return api_files

def find_functions(api_files):
    """ Return the list of functions contained in api.py files
    """
    
    result = []
    
    for api_file in api_files :
        api_globals = {}
        api_locals = {}
        
        sys.path.append(dirpath)
        execfile(api_file, api_globals, api_locals)
        sys.path.pop()
        
        result.append([api_file, api_locals.keys()])
    return result

def setup(project_name, main_script, includes):
    bin_directory = os.path.join("bin", "%s-%s"%(project_name, get_build_name()))
    
    import itk
    import itkTypes
    wrapitk_root = os.path.dirname(os.path.dirname(itkTypes.__file__))
    
    # Remove destination dir if it exists
    if os.path.exists(bin_directory) :
        shutil.rmtree(bin_directory)
    
    # Process specific includes
    if "itk" in includes :
        for file in os.listdir(os.path.join(wrapitk_root, "lib")) :
            if not file.startswith("_") :
                swig_module = os.path.splitext(file)[0]
                if swig_module not in includes :
                    includes.append(swig_module)
    if "medipy.itk" in includes :
        modules = ["itkNumpyBridgePython", "MediPyBridge", "MediPyBridgeConfig",
                   "MediPyBridgePython", "numpy_bridge"]
        for module in modules :
            includes.append("medipy.itk.%s"%module)

    # Main setup script
    sys.path.append(os.path.join(wrapitk_root, "lib"))
    distribution = distutils.core.setup(
        name = project_name,
        windows = [main_script], 
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
    
    # Copy configuration files from WrapITK
    configuration_dir = os.path.join(wrapitk_root, "Python", "Configuration")
    shutil.copytree(configuration_dir, os.path.join(bin_directory, "Configuration"))
    
    # Copy resources
    for dirpath, dirnames, filenames in os.walk("resources") :
        if ".svn" in dirnames :
            del dirnames[dirnames.index(".svn")]
        for filename in filenames :
            skip_file = (filename == "SConstruct" or
                         filename.endswith(".fbp"))
            if not skip_file :
                destination = os.path.join(bin_directory, dirpath)
                if not os.path.isdir(destination) :
                    os.makedirs(destination)
                print os.path.join(dirpath, filename), "->", os.path.join(destination, filename)
                shutil.copyfile(os.path.join(dirpath, filename), 
                                os.path.join(destination, filename))
    