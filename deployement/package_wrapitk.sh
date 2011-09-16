#!/bin/sh

##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

MAINTAINER=${MAINTAINER?"environment variable must be defined (e.g. 'John Doe <john@example.com>'"}
VENDOR=`lsb_release -is`
RELEASE=`lsb_release -rs`
ARCHITECTURE=`dpkg --print-architecture`

PACKAGE_NAME="wrapitk-0.3.0_${VENDOR}_${RELEASE}_${ARCHITECTURE}.deb"

if test "x$JOBS" = "x"
then
    MAKE_OPTIONS=""
else
    MAKE_OPTIONS="-j $JOBS"
fi

# Get the sources, unzip them, delete the archive
wget http://wrapitk.googlecode.com/files/wrapitk-0.3.0.tar.bz2
tar jxf wrapitk-0.3.0.tar.bz2
rm wrapitk-0.3.0.tar.bz2

# Prepare the build directory
mkdir wrapitk-0.3.0-build
cd wrapitk-0.3.0-build

# Prepare the CMakeCache.txt
cat <<EOF > CMakeCache.txt
BUILD_TESTING:BOOL=OFF

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
EOF

# Run cmake
cmake ../wrapitk-0.3.0

# Compile
make $MAKE_OPTIONS

# Install to local dir
mkdir -p local_install/usr
sed -e "s/^FILE/#FILE/" -i install_wrapitk_compatibility.cmake
cmake -D CMAKE_INSTALL_PREFIX=local_install/usr -P cmake_install.cmake
sed -e "s/^#FILE/FILE/" -i install_wrapitk_compatibility.cmake

# Install additional files
mkdir -p local_install/usr/lib/python2.6/dist-packages
mkdir -p local_install/usr/share/cmake-2.6/Modules
cp Languages/Python/InstallOnly/WrapITK.pth local_install/usr/lib/python2.6/dist-packages
cp InstallOnly/FindWrapITK.cmake local_install/usr/share/cmake-2.6/Modules

# Create binary package control file
mkdir -p local_install/DEBIAN
cat <<EOF > local_install/DEBIAN/control
Package: wrapitk
Version: 0.3.0
Section: libs
Priority: extra
Architecture: $ARCHITECTURE
Maintainer: $MAINTAINER
Description: Automated dynamic language binding for Insight Toolkit (ITK)
 WrapITK is an effort to automate the language binding process of one of
 the largest highly template-oriented c++ libraries, the Insight Toolkit
 image processing library.
 .
 Currently Python, Java and Tcl language bindings are implemented, but
 only Python is fully supported. For ITK .NET languages wrapper you may
 refer to ManagedITK.
EOF

# Build the package
fakeroot dpkg-deb -b local_install ./$PACKAGE_NAME

cd ..
