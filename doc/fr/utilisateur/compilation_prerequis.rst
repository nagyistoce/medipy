.. MediPy - Copyright (C) Universite de Strasbourg, 2011
   Distributed under the terms of the CeCILL-B license, as published by
   the CEA-CNRS-INRIA. Refer to the LICENSE file or to
   http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
   for details.

Compilation des pré-requis
==========================

Compilation des pré-requis
--------------------------

Les instructions d'ITK et de WrapITK sont décrites sur une :doc:`page propre <itk>`.

Python 2.5
^^^^^^^^^^

Le code de MediPy utilise Python 2.5, la version 2.6 est également acceptée.

* Télécharger `les sources de Python <http://www.python.org/download/>`_
* Les décompresser et les configurer : ``./configure --prefix=$HOME/local``
* Compiler : ``make``
* Installer : ``make install``
* Vérifier que l'appel à ``python -V`` donne bien une version 2.5

setuptools
^^^^^^^^^^

* Télécharger `le script d'installation <http://peak.telecommunity.com/dist/ez_setup.py>`_
* L'exécuter avec l'option ``--prefix=$HOME/local``
* Vérifier que l'appel à ``easy_install --help`` vous donne bien un message d'aide


SCons
^^^^^

* Récupérer les `sources <http://www.scons.org/download.php>`_
* Les décompresser et lancer ``python setup.py install --prefix=$HOME/local``

Docutils
^^^^^^^^

Installation : ``easy_install --prefix=$HOME/local docutils``

Numpy
^^^^^

Installation : ``easy_install --prefix=$HOME/local numpy``

Scipy
^^^^^

Installation : ``easy_install --prefix=$HOME/local scipy``

wxPython
^^^^^^^^

* récupérer les `sources <http://wxpython.org/download.php#sources>`_, les
  décompresser et se placer dans le répertoire nouvellement créé. Puis : ::

    mkdir bld
    cd bld
    ../configure --prefix=$HOME/local \
                 --with-gtk \
                 --with-gnomeprint \
                 --with-opengl \
                 --enable-optimise \
                 --enable-geometry \
                 --enable-graphics_ctx \
                 --enable-sound \
                 --with-sdl \
                 --enable-mediactrl \
                 --enable-display \
                 --enable-unicode 
    make
    make -C contrib/src/gizmos
    make -C contrib/src/stc
    make install
    make -C contrib/src/gizmos install
    make -C contrib/src/stc install
    cd ../wxPython
    python setup.py build_ext --inplace 
    python setup.py install --prefix=$HOME/local

L'installation peut être testée grâce au répertoire ``demo`` : ::

    cd demo
    python demo.py

VTK
^^^

* Récupérer les `sources de VTK <http://www.vtk.org/get-software.php#latest>`_
* Les décompresser (par exemple dans ``$HOME/src/VTK``)
* Dans un *autre répertoire* (par exemple ``$HOME/src/VTK-build``), lancer
  ``ccmake $HOME/src/VTK``
* Modifiez les options de configuration suivantes (certaines options sont dans
  les options avancées, accessibles par ``t`` : 

    * ``BUILD_DOCUMENTATION``, ``BUILD_EXAMPLES`` et ``BUILD_TESTING`` peuvent
      être à ``OFF`` pour accélérer la compilation.
    * ``CMAKE_BUILD_TYPE`` vaudra soit ``Release`` soit ``RelWithDebInfo`` selon
      que vous souhaitiez la présence des informations de débuggage dans les
      bibliothèques générées.
    * ``CMAKE_INSTALL_PREFIX`` vaudra ``$HOME/local``
    * ``PYTHON_INCLUDE_PATH`` et ``PYTHON_LIBRARY`` pointeront sur votre
      installation de Python 2.5.
    * ``VTK_WRAP_PYTHON`` devra valoir ``ON``
    * ``VTK_USE_PARALLEL`` devra valoir ``ON``

* Une fois le ``Makefile`` généré, compilez (``make``) et installez (``make install``) VTK.


PyNifti
^^^^^^^

Récupérer `les sources <http://sourceforge.net/project/showfiles.php?group_id=126549>`_
et les décompresser, puis : ::

    make
    python setup.py install --prefix=$HOME/local


