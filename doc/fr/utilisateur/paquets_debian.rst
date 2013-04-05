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
    apt-get install scons mercurial python-setuptools swig/squeeze ipython python-tk \
        cmake doxygen gccxml cableswig g++-4.4
    echo "swig hold" | dpkg --set-selections

**Bilbiothèques** ::

    sudo apt-get install python-vtk libvtk5-dev python-wxgtk2.8 python-scipy \
        libboost-python-dev python-nifti python-docutils python-wxtools

ITK et WrapITK sont aussi nécessaires : voir :doc:`la page les concernant <itk>`.

Debian 6.0 (Squeeze)
--------------------

**Outils** ::

    sudo aptitude install scons mercurial python-setuptools swig ipython cmake doxygen gccxml cableswig make

**Bilbiothèques** ::

    sudo aptitude install python-vtk libvtk5-dev python-wxgtk2.8 python-scipy \
        libboost-python-dev python-nifti python-docutils python-wxtools


ITK et WrapITK sont aussi nécessaires : voir :doc:`la page les concernant <itk>`.

Debian 5.0 (Lenny)
------------------

L'installation de la distribution a été réalisée avec les options par défaut
pour tasksel (*Desktop environment* et *Base system*).

**Outils** ::

    sudo aptitude install build-essential scons subversion python-setuptools swig ipython 

**Bilbiothèques** ::

    sudo aptitude install python-vtk libvtk5-dev python-wxgtk2.8 python-scipy \
        libboost-python-dev python-nifti python-docutils</code>

ITK et WrapITK sont aussi nécessaires : voir :doc:`la page les concernant <itk>`.
