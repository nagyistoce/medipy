#!/usr/bin/env python

import os
import shutil
import subprocess
import sys
import tarfile
import urllib
    
# Get the sources, unzip them, delete the archive
filename,_ = urllib.urlretrieve("http://wrapitk.googlecode.com/files/wrapitk-0.3.0.tar.bz2")
archive = tarfile.open(filename)
archive.extractall()
archive.close()
os.remove(filename)

# Prepare the build directory
build_directory = os.path.join(os.getcwd(), "wrapitk-0.3.0-build")
os.makedirs(build_directory)

# Prepare the CMakeCache.txt
cmake_cache = open(os.path.join(build_directory, "CMakeCache.txt"), "w")
cmake_cache.write("""BUILD_TESTING:BOOL=OFF

CMAKE_BUILD_TYPE:STRING=RELEASE
CMAKE_INSTALL_PREFIX:PATH=/usr

INSTALL_WRAP_ITK_COMPATIBILITY:BOOL=ON

WRAP_ITK_DIMS:STRING=2;3
WRAP_ITK_DOC:BOOL=ON
WRAP_ITK_GCCXML:BOOL=ON
WRAP_ITK_PYTHON:BOOL=ON
WRAP_ITK_SWIGINTERFACE:BOOL=ON

WRAP_complex_double:BOOL=OFF
WRAP_complex_float:BOOL=ON

WRAP_covariant_vector_double:BOOL=OFF
WRAP_covariant_vector_float:BOOL=ON

WRAP_double:BOOL=OFF
WRAP_float:BOOL=ON

WRAP_rgb_unsigned_char:BOOL=OFF
WRAP_rgb_unsigned_short:BOOL=ON
WRAP_rgba_unsigned_char:BOOL=OFF
WRAP_rgba_unsigned_short:BOOL=ON

WRAP_signed_char:BOOL=ON
WRAP_signed_long:BOOL=ON
WRAP_signed_short:BOOL=ON
WRAP_unsigned_char:BOOL=ON
WRAP_unsigned_long:BOOL=ON
WRAP_unsigned_short:BOOL=ON

WRAP_vector_double:BOOL=OFF
WRAP_vector_float:BOOL=ON
""")
if sys.platform != "win32" :
    cmake_cache.write("CMAKE_INSTALL_PREFIX:STRING=/usr\n")
    # TODO bigobj
else :
    cmake_cache.write("CMAKE_INSTALL_PREFIX:STRING=C:/Program Files/ITK\n")
cmake_cache.close()

# Run cmake
subprocess.call(
    ["cmake", os.path.join(os.getcwd(), "wrapitk-0.3.0")], 
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
