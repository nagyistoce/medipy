Installation
============

Depuis les sources
------------------

Les outils suivants sont nécessaires pour installer MediPy depuis les sources :

* `Mercurial <http://mercurial.selenic.com/>`_
* `CMake <www.cmake.org/>`_
* `SWIG <http://www.swig.org/>`_
* `GCC-XML <http://www.gccxml.org/HTML/Index.html>`_
* `CableSwig <http://www.itk.org/ITK/resources/CableSwig.html>`_

Les bibliothèques suivantes sont nécessaires à la compilation et à l'exécution :

* `wxPython <http://www.wxpython.org/>`_
* `NumPy <http://numpy.scipy.org/>`_ et `SciPy <http://www.scipy.org/>`_
* `PyNIfTI <http://niftilib.sourceforge.net/pynifti/>`_
* :doc:`ITK et WrapITK <itk>` (seule la version 3.14 d'ITK a été testée)
* `VTK  <http://www.vtk.org/>`_ (seules les versions 5.0, 5.2 et 5.8 de VTK ont
  été testées)

Des listes de paquets sont disponibles pour certaines distributions de Linux :

* :doc:`Paquets pour Ubuntu <paquets_ubuntu>`
* :doc:`Paquets pour Debian <paquets_debian>`

Préparation de l'environnement
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sous Linux, certaines variables d'environnement doivent être mises à jour 
(``~/.bashrc`` ou ``~/.profile``) : ::

    export MEDIPY_HOME=$HOME/src/medipy
    export PYTHONPATH=$MEDIPY_HOME/lib${PYTHONPATH:+:$PYTHONPATH}
    export MEDIPY_PLUGINS_PATH=$MEDIPY_HOME/plugins/

.. important::
    
    Démarrez un nouveau terminal pour que ces modifications soient prises en compte.

Récupération des sources
^^^^^^^^^^^^^^^^^^^^^^^^

L'ensemble des sources de MediPy est stocké sur un `serveur 
Mercurial <http://fr.wikipedia.org/wiki/Mercurial>`_, à l'adresse
https://medipy.googlecode.com/hg/. Les sources peuvent donc être récupérées par : ::

    hg clone https://medipy.googlecode.com/hg/ $MEDIPY_HOME

Compilation de MediPy
^^^^^^^^^^^^^^^^^^^^^

`CMake <www.cmake.org/>`_ ne compile pas directement le code, mais permet de
générer un ``Makefile`` (ou équivalent, selon la plate-forme). Pour générer le
``Makefile``, on lancera (à la racine de MediPy) : ::

    ccmake .

L'interface permet alors de choisir des version spécifiques des outils et 
bibliothèques nécessaires à la compilation de MediPy, mais aussi de choisir les
options de compilation : l'option ``CMAKE_BUILD_TYPE`` permet de choisir un
type de compilation (``Debug`` ou ``Release`` pour les valeurs les plus 
courantes).

Cette étape de génération du ``Makefile`` peut également être réalisée sans 
passer par l'interface. Pour créer un ``Makefile`` en version optimisée en 
utilisant les valeurs de base des options, il suffit de lancer (à la racine de 
MediPy) : ::

    cmake -D CMAKE_BUILD_TYPE=Release .

La compilation peut ensuite être réalisée en lançant ``make`` si on se trouve
à la racine de MediPy, ou ``make -C $MEDIPY_HOME`` dans le cas contraire. La 
compilation en parallèle (option ``-j N`` de ``make`` où ``N`` est le nombre de
cœurs de la machine) fonctionne correctement ; il faut cependant noter que les
éventuels messages d'erreur seront moins lisibles.

Le nettoyage des sources est effectué par ``make clean``. Sous Windows, il est
nécessaire d'utiliser le script ``win32-build.bat``, qui assigne correctement 
les variables d'environnement.

MediPy peut ensuite être lancé par : ::

    $MEDIPY_HOME/apps/medipy/medipy

Options en ligne de commande
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``-d <LEVEL>``, ``--debug=LEVEL``
    Active les informations de débuggage. ``LEVEL`` peut valoir (du moins précis
    au plus précis) ``critical``, ``error``, ``warning`` ou ``debug``

``-f <FILE>``, ``--file=FILE``
    Exécute le contenu du fichier donné en paramètre après le lancement de 
    l'application

``-m <FILE>``, ``--menu-file=FILE``
    Utilise le fichier de menu donné en paramètre. Si cette option est absente,
    les fichiers ``api.py`` seront utilisés pour générer le menu.

``-r`` 
    Redirige les messages d'erreur vers la console. Par défault, ils seront
    affichés dans une fenêtre wx
