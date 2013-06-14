Paquets pour Debian
===================

Debian 7.0 (Wheezy)
-------------------

Attention, la version de Swig fournie avec Debian 7.0 ne permet plus, à 
l'heure actuelle, de compiler les wrappers ITK. Il faut donc installer la 
version de Debian 6.0, et marquer le paquet comme étant à conserver. Les
commandes suivantes réalisent cette procédure.

**Outils** ::

    echo "deb http://ftp.u-strasbg.fr/debian/ squeeze main" >> /etc/apt/sources.list
    apt-get update
    apt-get install mercurial python-setuptools swig/squeeze ipython python-tk \
        cmake doxygen gccxml cableswig g++-4.4 cython
    echo "swig hold" | dpkg --set-selections

**Bilbiothèques** ::

    sudo apt-get install python-vtk libvtk5-dev python-wxgtk2.8 python-scipy \
        python-nifti python-docutils python-wxtools uuid-dev dcmtk \
        libdcmtk2-dev libwrap0-dev libgdcm2-dev

ITK et WrapITK sont aussi nécessaires : voir :doc:`la page les concernant <itk>`.

Debian 6.0 (Squeeze)
--------------------

**Outils** ::

    sudo aptitude install mercurial python-setuptools swig ipython cmake \
        doxygen gccxml cableswig make cython

**Bilbiothèques** ::

    sudo aptitude install python-vtk libvtk5-dev python-wxgtk2.8 python-scipy \
        python-nifti python-docutils python-wxtools uuid-dev dcmtk \
        libdcmtk1-dev libwrap0-dev libgdcm2-dev


ITK et WrapITK sont aussi nécessaires : voir :doc:`la page les concernant <itk>`.
