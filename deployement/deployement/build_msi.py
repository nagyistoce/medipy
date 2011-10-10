##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import os
import os.path
import shutil
import subprocess

from setup import get_build_name

def build_msi(project_name, wix_source = None):
    
    if wix_source is None :
        wix_source = "%s.wxs"%project_name
    
    build = get_build_name()

    source_dir = os.path.abspath(os.getcwd())
    working_dir = "bin\\%s-%s"%(project_name, build)

    # Harvest files
    heat_command = ["heat", "dir", ".", "-cg", "%s_files"%project_name, 
                    "-out", "files.wxs", "-ag", "-sfrag", "-dr", 
                    "APPLICATIONROOTDIRECTORY", "-srd"]
    subprocess.call(heat_command, cwd=working_dir)

    # Compile WiX source files
    source_files = [
        os.path.join(source_dir, "wix", wix_source),
        os.path.join(os.path.dirname(__file__), "wix", "WixUI_InstallDirNoLicense.wxs"),
        "files.wxs"
    ]
    
    for file in source_files :
        command = ["candle", "-nologo", "%s"%file]
        subprocess.call(command, cwd=working_dir)

    # Create the MSI
    msi_filename = "%s-%s.msi"%(project_name, build)
    light_command = ["light", "-nologo", "-ext", "WixUIExtension", "-cultures:fr-fr,en-us",
                     "-out", msi_filename,
                     "%s.wixobj"%project_name, "WixUI_InstallDirNoLicense.wixobj", 
                     "files.wixobj",
                     ]
    subprocess.call(light_command, cwd=working_dir)

    # Clean up
    to_delete = [
        "files.wxs",
        "%s.wixobj"%project_name,
        "WixUI_InstallDirNoLicense.wixobj",
        "files.wixobj",
        "%s-%s.wixpdb"%(project_name, build)
    ]
    for file in to_delete :
        os.remove(os.path.join(working_dir, file))
    shutil.move(os.path.join(working_dir, msi_filename), "bin")