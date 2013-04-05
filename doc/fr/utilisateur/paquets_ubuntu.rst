Paquets pour Ubuntu
===================

Ubuntu 12.04 LTS (Precise Pangolin)
-----------------------------------

Attention, la version de Swig fournie avec Ubuntu 12.04 ne permet plus, à 
l'heure actuelle, de compiler les wrappers ITK. Il faut donc installer la 
version d'Ubuntu 10.04, et marquer le paquet comme étant à conserver. Les
commandes suivantes réalisent cette procédure.

**Outils** ::

    echo "deb http://fr.archive.ubuntu.com/ubuntu lucid main restricted universe multiverse" >> /etc/apt/sources.list
    apt-get update
    apt-get install scons mercurial python-setuptools swig/lucid ipython python-tk \
        cmake doxygen gccxml cableswig g++-4.4
    echo "swig hold" | dpkg --set-selections

**Bilbiothèques** ::

    sudo apt-get install python-vtk libvtk5-dev python-wxgtk2.8 python-scipy \
        libboost-python-dev python-nifti python-docutils python-wxtools

ITK et WrapITK sont aussi nécessaires : voir :doc:`la page les concernant <itk>`.

Ubuntu 10.04 LTS (Lucid Lynx), 9.10 (Karmic Koala) et 9.04 (Jaunty Jackalope)
-----------------------------------------------------------------------------

**Outils** ::

    sudo aptitude install scons mercurial python-setuptools swig ipython cmake doxygen gccxml cableswig

**Bilbiothèques** ::

    sudo aptitude install python-vtk libvtk5-dev python-wxgtk2.8 python-scipy \
        libboost-python-dev python-nifti python-docutils python-wxtools


ITK et WrapITK sont aussi nécessaires : voir :doc:`la page les concernant <itk>`.

Ubuntu 8.04 LTS (Hardy Heron)
-----------------------------

Paquets de la distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    sudo aptitude install python-setuptools python-docutils python-numpy \
        python-scipy python-wxgtk2.8 python-vtk subversion \
        libstdc++6-dev g++ python-dev swig libboost-python-dev \
        libvtk5-dev zlib1g-dev


SCons
^^^^^

::

    cd ~/src
    wget http://prdownloads.sourceforge.net/scons/scons-1.2.0.tar.gz
    tar zxf scons-1.2.0.tar.gz
    rm scons-1.2.0.tar.gz
    cd scons-1.2.0
    python setup.py build
    python setup.py install --prefix=$HOME/local


PyNifti
^^^^^^^

::

    cd ~/src
    wget http://downloads.sourceforge.net/niftilib/pynifti_0.20090303.1.tar.gz
    tar zxf pynifti_0.20090303.1.tar.gz
    rm pynifti_0.20090303.1.tar.gz
    cd pynifti-0.20090303.1
    make
    python setup.py install --prefix=$HOME/local

ITK et WrapITK sont aussi nécessaires : voir :doc:`la page les concernant <itk>`.
