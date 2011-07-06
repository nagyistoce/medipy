Installation d'ITK et de WrapITK
================================

Paquets
-------

Des paquets pour `ITK <http://www.itk.org/>`_ sont disponibles pour Debian 5.0,
Ubuntu 8.10, Ubuntu 9.04 et Ubuntu 9.10. Ces paquets existent pour les
architectures ``i386`` et ``amd64``. Ils sont disponibles dans un dépôt non
officiel, qui peut être ajouté au système en ajoutant une ligne au fichier
``/etc/apt/sources.list`` :
  
* Debian 5.0 : ``deb http://download.opensuse.org/repositories/home:/tonddast:/itk/Debian_5.0/ ./``
* Ubuntu 8.10 : ``deb http://download.opensuse.org/repositories/home:/tonddast:/itk/xUbuntu_8.10/ ./``
* Ubuntu 9.04 : ``deb http://download.opensuse.org/repositories/home:/tonddast:/itk/xUbuntu_9.04/ ./``
* Ubuntu 9.10 : ``deb http://download.opensuse.org/repositories/home:/tonddast:/itk/xUbuntu_9.10/ ./``

ITK peut ensuite être installé par ``aptitude install libitk libitk-dev``.

Pour Ubuntu 10.04 LTS, un paquet est disponible dans ``/dskcommun/lamy/ITK``.

Des paquets existent également pour `WrapITK <http://code.google.com/p/wrapitk/>`_
pour Debian 5.0 (i386), Ubuntu 9.04 (i386 et amd64) et Ubuntu 9.10 (i386).
Ils sont disponibles sur notre `serveur FTP <ftp://alsace.u-strasbg.fr/pub/tmp/wrapitk/>`_
et dans ``/dskcommun/lamy/wrapitk``.

Construction d'un paquet
------------------------

Si aucun paquet n'existe pour votre distribution, il peut être intéressant d'en
fabriquer un et de l'héberger sur le `serveur FTP <ftp://alsace.u-strasbg.fr/pub/tmp/wrapitk/>`_.
`CMake <http://www.cmake.org/>`_ et `fakeroot <http://fakeroot.alioth.debian.org/>`_
sont nécessaires pour la construction du paquet.

ITK
^^^

* Génération du Makefile : ::

    mkdir InsightToolkit-3.14.0-build
    cd InsightToolkit-3.14.0-build
    ccmake ../InsightToolkit-3.14.0

* Options de cmake : 

  * BUILD_TESTING=OFF
  * CMAKE_BUILD_TYPE=Release
  * CMAKE_INSTALL_PREFIX=/usr
  * BUILD_SHARED_LIBS=ON
  * BUILD_DOXYGEN=ON
  * BUILD_DOCUMENTATION=ON
  * BUILD_EXAMPLES=OFF
  * ITK_USE_REVIEW=ON

* Compilation : la compilation parallèle (``make -j n``) fonctionne sans problème ::

    make -j 3

* Installation locale : ::

    mkdir -p local_install/usr
    cmake -D CMAKE_INSTALL_PREFIX=local_install/usr -P cmake_install.cmake

* Création du fichier de contrôle du paquet (adapter les champs ``Architecture``
  et ``Maintainer``) : ::

    mkdir -p local_install/DEBIAN
    cat <<EOF > local_install/DEBIAN/control
    Package: InsightToolkit
    Version: 3.14.0-lamy0
    Section: libs
    Priority: extra
    Architecture: i386
    Maintainer: Julien Lamy <lamy@unistra.fr>
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

* Intégration du répertoire ``/usr/local/InsightToolkit`` au chemin de recherche
  de ``ld`` ::

    mkdir -p local_install/etc/ld.so.conf.d
    echo /usr/lib/InsightToolkit > local_install/etc/ld.so.conf.d/itk.conf
    cat <<EOF > local_install/DEBIAN/postinstall
    #!/bin/sh
    /sbin/ldconfig
    EOF
    chmod 755 local_install/DEBIAN/postinstall
    cp -a local_install/DEBIAN/postinstall local_install/DEBIAN/prerm

* Création du paquet : ::

    fakeroot dpkg-deb -b local_install ./itk-3.14.0-lamy0.deb


WrapITK
^^^^^^^

* Récupération des sources : ::

    cd ~/src
    wget http://wrapitk.googlecode.com/files/wrapitk-0.3.0.tar.bz2
    tar jxf wrapitk-0.3.0.tar.bz2

* Génération du Makefile : ::

    mkdir wrapitk-0.3.0-build
    cd wrapitk-0.3.0-build
    ccmake ../wrapitk-0.3.0

* Options de cmake : 

  * BUILD_TESTING=OFF
  * CMAKE_BUILD_TYPE=Release
  * CMAKE_INSTALL_PREFIX=/usr
  * WRAP_ITK_DOC=ON
  * WRAP_ITK_PYTHON=ON
  * INSTALL_WRAP_ITK_COMPATIBILITY=ON
  * compiler pour tous les types entiers
  
* Compilation : la compilation parallèle (``make -j n``) fonctionne sans 
  problème ::

    make -j 3

* Installation locale : ::

    mkdir -p local_install/usr
    sed -e "s/^FILE/#FILE/" -i install_wrapitk_compatibility.cmake
    cmake -D CMAKE_INSTALL_PREFIX=local_install/usr -P cmake_install.cmake
    sed -e "s/^#FILE/FILE/" -i install_wrapitk_compatibility.cmake

* Copie des fichiers « InstallOnly » ::

    mkdir -p local_install/usr/lib/python2.6/dist-packages 
    mkdir -p local_install/usr/share/cmake-2.6/Modules
    cp Languages/Python/InstallOnly/WrapITK.pth local_install/usr/lib/python2.6/dist-packages
    cp InstallOnly/FindWrapITK.cmake local_install/usr/share/cmake-2.6/Modules

* Création du fichier de contrôle du paquet (adapter les champs ``Architecture``
  et ``Maintainer``) ::

    mkdir -p local_install/DEBIAN
    cat <<EOF > local_install/DEBIAN/control
    Package: wrapitk
    Version: 0.3.0
    Section: libs
    Priority: extra
    Architecture: i386
    Maintainer: Julien Lamy <lamy@unistra.fr>
    Description: Automated dynamic language binding for Insight Toolkit (ITK)
     WrapITK is an effort to automate the language binding process of one of 
     the largest highly template-oriented c++ libraries, the Insight Toolkit 
     image processing library.
     .
     Currently Python, Java and Tcl language bindings are implemented, but 
     only Python is fully supported. For ITK .NET languages wrapper you may 
     refer to ManagedITK.
    EOF

* Création du paquet : ::

    fakeroot dpkg-deb -b local_install ./wrapitk-0.3.0.deb

