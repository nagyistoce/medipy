#!/usr/bin/env python

import os
import shutil
import subprocess
import sys
import tarfile
import urllib
    
# Get the sources, unzip them, delete the archive
filename,_ = urllib.urlretrieve("http://sourceforge.net/projects/itk/files/itk/3.14/InsightToolkit-3.14.0.tar.gz/download")
archive = tarfile.open(filename)
archive.extractall()
archive.close()
os.remove(filename)

# Prepare the build directory
build_directory = os.path.join(os.getcwd(), "InsightToolkit-3.14.0-build")
os.makedirs(build_directory)

# Patch InsightToolkit-3.14.0/Utilities/MetaIO/CMakeLists.txt to use VTK and
# GDCM packages
source_directory = os.path.join(os.getcwd(), "InsightToolkit-3.14.0")
patch_file = os.path.join(os.path.dirname(__file__), "MetaIO_CMakeLists.patch")
subprocess.call(["patch", "-p1"], stdin=open(patch_file))

# Prepare the CMakeCache.txt
cmake_cache = open(os.path.join(build_directory, "CMakeCache.txt"), "w")
cmake_cache.write("""BUILD_TESTING:BOOL=OFF
CMAKE_BUILD_TYPE:STRING=RelWithDebInfo
BUILD_SHARED_LIBS:BOOL=ON
BUILD_DOXYGEN:BOOL=ON
BUILD_DOCUMENTATION:BOOL=ON
BUILD_EXAMPLES:BOOL=OFF
ITK_USE_REVIEW:BOOL=ON
ITK_USE_SYSTEM_GDCM:BOOL=ON
""")
if sys.platform != "win32" :
    command = ["gcc", "-dumpversion"]
    process = subprocess.Popen(command, stdout=subprocess.PIPE)
    stdout, _ = process.communicate()
    elements = stdout.split(".")
    major, minor = int(elements[0]), int(elements[1])
    if major > 4 or (major == 4 and minor > 4) :
        cmake_cache.write("CMAKE_C_COMPILER:FILEPATH=gcc-4.4\n")
        cmake_cache.write("CMAKE_CXX_COMPILER:FILEPATH=g++-4.4\n")

    cmake_cache.write("CMAKE_INSTALL_PREFIX:STRING=/usr\n")
else :
    cmake_cache.write("CMAKE_INSTALL_PREFIX:STRING=C:/Program Files/ITK\n")
cmake_cache.close()

# Run cmake
subprocess.call(
    ["cmake", os.path.join(os.getcwd(), "InsightToolkit-3.14.0")], 
    cwd=build_directory)

# Compile
if sys.platform == "linux2" :
    subprocess.call(["make"]+sys.argv[1:], cwd=build_directory)
elif sys.platform == "win32" :
    subprocess.call(
        ["devenv", "/build", "RelWithDebInfo", 
         "/project", "ALL_BUILD", "ITK.sln"],
        cwd=build_directory)
else : 
    raise Exception("Unknown platform {0}".format(sys.platform))

# Install to local dir
if sys.platform == "win32" :
    local_install_dir = os.path.join(build_directory, "local_install")
else :
    local_install_dir = os.path.join(build_directory, "local_install", "usr")
os.makedirs(local_install_dir)
subprocess.call(
    ["cmake", "-D", "CMAKE_INSTALL_PREFIX={0}".format(local_install_dir),
     "-D", "BUILD_TYPE=RelWithDebInfo", "-P", "cmake_install.cmake"], 
    cwd=build_directory)

# Package
if sys.platform == "win32" :
    data = open(os.path.join(os.path.dirname(__file__), "itk.nsi.in")).read()
    data = data.format(source_dir=local_install_dir)
    nsi_filename = os.path.join(build_directory, "itk.nsi")
    open(nsi_filename, "w").write(data)
    subprocess.call(["makensis", "/V2", nsi_filename], cwd=build_directory)
elif sys.platform == "linux2" :
    if "MAINTAINER" not in os.environ :
        raise Exception("MAINTAINER environment variable must be defined "
                        "(e.g. 'John Doe <john@example.com>'")
    
    maintainer = os.environ["MAINTAINER"]
    architecture = subprocess.Popen(
        ["dpkg", "--print-architecture"], 
        stdout=subprocess.PIPE).communicate()[0].strip()
    
    debian_dir = os.path.join(build_directory, "local_install", "DEBIAN")
    
    # Create binary package control file
    os.makedirs(debian_dir)
    data = open(os.path.join(os.path.dirname(__file__), "DEBIAN", "control.in")).read()
    data = data.format(architecture=architecture, maintainer=maintainer)
    open(os.path.join(debian_dir, "control"), "w").write(data)
    
    # Add library dir to ld
    os.makedirs(os.path.join(build_directory, "local_install", "etc", 
                             "ld.so.conf.d"))
    shutil.copy(os.path.join(os.path.dirname(__file__), "DEBIAN", "itk.conf"),
                os.path.join(build_directory, "local_install", "etc", "ld.so.conf.d"))
    for file in ["postinst", "prerm"] :
        shutil.copy(os.path.join(os.path.dirname(__file__), "DEBIAN", file),
                    os.path.join(debian_dir, file))
        os.chmod(os.path.join(debian_dir, file), 0755)
    
    vendor = subprocess.Popen(
        ["lsb_release", "-is"], stdout=subprocess.PIPE).communicate()[0].strip()
    release = subprocess.Popen(
        ["lsb_release", "-rs"], stdout=subprocess.PIPE).communicate()[0].strip()
    package_name = "insighttoolkit-3.14.0_{0}_{1}_{2}.deb".format(vendor, release, architecture)
    
    subprocess.call(
        ["fakeroot", "dpkg-deb", "-b", "local_install", package_name],
        cwd=build_directory)
else :
    raise Exception("Unknown platform {0}".format(sys.platform))
