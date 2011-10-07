Installation
============

Depuis les sources
------------------

Les outils suivants sont nécessaires pour installer MediPy depuis les sources :

* `Mercurial <http://mercurial.selenic.com/>`_
* `SCons <http://www.scons.org/>`_
* `SWIG <http://www.swig.org/>`_
* `GCC-XML <http://www.gccxml.org/HTML/Index.html>`_
* `CableSwig <http://www.itk.org/ITK/resources/CableSwig.html>`_

Les bibliothèques suivantes sont nécessaires à la compilation et à l'exécution :

* `wxPython <http://www.wxpython.org/>`_
* `NumPy <http://numpy.scipy.org/>`_ et `SciPy <http://www.scipy.org/>`_
* `Boost.Python <http://www.boost.org/doc/libs/1_42_0/libs/python/doc/index.html>`_
* `PyNIfTI <http://niftilib.sourceforge.net/pynifti/>`_
* :doc:`ITK et WrapITK <itk>` (seule la version 3.14 d'ITK a été testée)
* `VTK  <http://www.vtk.org/>`_ (seules les versions 5.0 et 5.2 de VTK ont été testées)

Ces outils et bibliothèques sont tous compilables depuis les sources, mais il
est préférable d'utiliser les paquets d'une distribution Linux ou les binaires
Windows que de :doc:`tout recompiler <compilation_prerequis>`. Dans le cas des
paquets Python, l'utilisation de 
`setuptools <http://pypi.python.org/pypi/setuptools>`_ permet également d'éviter
la recompilation pour des paquets dont les binaires ne seraient pas disponibles
facilement.

Des listes de paquets sont disponibles pour certaines distributions de Linux :

* :doc:`Paquets pour Ubuntu <paquets_ubuntu>`
* :doc:`Paquets pour Debian <paquets_debian>`

Préparation de l'environnement
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sous Linux il est conseillé de créer une hiérarchie locale dans votre répertoire
personnel. Les exemples suivants supposent l'existence de la hiérarchie 
``$HOME/local``, contenant les répertoires ``bin``, ``lib`` et ``include`` : ::

    mkdir -p $HOME/local/{bin,include,lib/python`python -c "import sys; print sys.version[:3]"`/site-packages}


Certaines variables d'environnement doivent ensuite être mises à jour 
(``~/.bashrc`` ou ``~/.profile``) : ::

    export PATH=$HOME/local/bin:$PATH
    export LD_LIBRARY_PATH=$HOME/local/lib:$HOME/src/medipy/build/lib:$LD_LIBRARY_PATH
    export PYTHONPATH=$HOME/local/lib/python`python -c "import sys; print sys.version[:3]"`/site-packages/:$PYTHONPATH
    export MEDIPY_PLUGINS_PATH=$HOME/src/medipy/build/medipy/plugins/

Démarrez un nouveau terminal pour que ces modifications soient prises en compte.

Récupération des sources
^^^^^^^^^^^^^^^^^^^^^^^^

L'ensemble des sources de MediPy est stocké sur un serveur 
`Mercurial <http://fr.wikipedia.org/wiki/Mercurial>`_, à l'adresse
https://medipy.googlecode.com/hg/. Les sources peuvent donc être récupérées par : ::

    hg clone https://medipy.googlecode.com/hg/ $HOME/src/medipy

Compilation de MediPy
^^^^^^^^^^^^^^^^^^^^^

Le script de compilation accepte les options suivantes (la liste complète est
visible par un appel à ``scons --help``) : 

* ``optimized`` : active les options de compilation pour la version optimisée
  ou pour la version debug
* ``verbose`` : affiche plus de messages lors de la compilation

Ainsi, pour compiler MediPy en version optimisée, on lancera (à la racine de
MediPy) : ::

    scons optimized=1

ou, de façon équivalente : ::

    scons optimized=True

Le nettoyage des sources est effectué par ``scons -c`` (équivalent à 
``make clean``). Sous Windows, il est nécessaire d'utiliser le 
``win32-build.bat``, qui assigne correctement les variables d'environnement.

Sous Linux, afin d'intégrer MediPy aux paquets Python, il faut créer un lien
symbolique depuis le répertoire contenant les fichiers compilés
(``build/medipy``) vers la hiérarchie locale : ::

    ln -s $HOME/src/medipy/build/medipy/lib $HOME/local/lib/python`python -c "import sys; print sys.version[:3]"`/site-packages/medipy

MediPy peut ensuite être lancé par : ::

    $HOME/src/medipy/apps/medipy/medipy

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
