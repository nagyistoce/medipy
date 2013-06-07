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
    apt-get install mercurial python-setuptools swig/lucid ipython python-tk \
        cmake doxygen gccxml cableswig g++-4.4
    echo "swig hold" | dpkg --set-selections

**Bilbiothèques** ::

    sudo apt-get install python-vtk libvtk5-dev python-wxgtk2.8 python-scipy \
        python-nifti python-docutils python-wxtools uuid-dev libdcmtk2-dev \
        libwrap0-dev

ITK et WrapITK sont aussi nécessaires : voir :doc:`la page les concernant <itk>`.

Ubuntu 10.04 LTS (Lucid Lynx)
-----------------------------

**Outils** ::

    sudo aptitude install mercurial python-setuptools swig ipython cmake doxygen gccxml cableswig

**Bilbiothèques** ::

    sudo aptitude install python-vtk libvtk5-dev python-wxgtk2.8 python-scipy \
        python-nifti python-docutils python-wxtools uuid-dev libxml2-dev \
        libdcmtk1-dev libwrap0-dev


ITK et WrapITK sont aussi nécessaires : voir :doc:`la page les concernant <itk>`.


