#!/bin/sh

MAINTAINER=${MAINTAINER?"environment variable must be defined (e.g. 'John Doe <john@example.com>'"}
VENDOR=`lsb_release -is`
RELEASE=`lsb_release -rs`
ARCHITECTURE=`dpkg --print-architecture`

PACKAGE_NAME="itk-3.14.0_${VENDOR}_${RELEASE}_${ARCHITECTURE}.deb"

if test "x$JOBS" = "x"
then
    MAKE_OPTIONS=""
else
    MAKE_OPTIONS="-j $JOBS"
fi

# Get the sources, unzip them, delete the archive
wget http://voxel.dl.sourceforge.net/sourceforge/itk/InsightToolkit-3.14.0.tar.gz
tar zxf InsightToolkit-3.14.0.tar.gz
rm InsightToolkit-3.14.0.tar.gz

# Prepare the build directory
mkdir InsightToolkit-3.14.0-build
cd InsightToolkit-3.14.0-build

# Prepare the CMakeCache.txt
cat <<EOF > CMakeCache.txt
BUILD_TESTING:BOOL=OFF
CMAKE_BUILD_TYPE:STRING=Release
CMAKE_INSTALL_PREFIX:STRING=/usr
BUILD_SHARED_LIBS:BOOL=ON
BUILD_DOXYGEN:BOOL=ON
BUILD_DOCUMENTATION:BOOL=ON
BUILD_EXAMPLES:BOOL=OFF
ITK_USE_REVIEW:BOOL=ON
EOF

# Run cmake
cmake ../InsightToolkit-3.14.0

# Compile
make $MAKE_OPTIONS

# Install to local dir
mkdir -p local_install/usr
cmake -D CMAKE_INSTALL_PREFIX=local_install/usr -P cmake_install.cmake

# Create binary package control file
mkdir -p local_install/DEBIAN
cat <<EOF > local_install/DEBIAN/control
Package: InsightToolkit
Version: 3.14.0
Section: libs
Priority: extra
Architecture: $ARCHITECTURE
Maintainer: $MAINTAINER
Description: Image processing toolkit for registration and segmentation
 ITK is an open-source software toolkit for performing registration and 
 segmentation. Segmentation is the process of identifying and classifying 
 data found in a digitally sampled representation. Typically the sampled 
 representation is an image acquired from such medical instrumentation as 
 CT or MRI scanners. Registration is the task of aligning or developing 
 correspondences between data. For example, in the medical environment, a
 CT scan may be aligned with a MRI scan in order to combine the 
 information contained in both.
EOF

# Add library dir to ld
mkdir -p local_install/etc/ld.so.conf.d
echo /usr/lib/InsightToolkit > local_install/etc/ld.so.conf.d/itk.conf
cat <<EOF > local_install/DEBIAN/postinstall
#!/bin/sh
/sbin/ldconfig
EOF
chmod 755 local_install/DEBIAN/postinstall
cp -a local_install/DEBIAN/postinstall local_install/DEBIAN/prerm

# Build the package
fakeroot dpkg-deb -b local_install ./$PACKAGE_NAME

cd ..
