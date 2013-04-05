Installation d'ITK et de WrapITK
================================

Paquets
-------

Nous mettons à disposition des paquets pour ITK 3.14 et WrapITK/Python 0.3.0 
pour Ubuntu (9.04, 10.04 et 12.04, 32 bits et 64 bits) et Debian (6.0 et 7.0, 
64 bits). Ces paquets sont disponibles sur le 
`site Google Code de MediPy <http://code.google.com/p/medipy/downloads/list>`_.

Construction d'un paquet
------------------------

Si aucun paquet n'existe pour votre distribution, nous avons publié deux scripts
permettant de créer des paquets ``.deb`` pour ITK et WrapITK.
`CMake <http://www.cmake.org/>`_ et 
`fakeroot <http://fakeroot.alioth.debian.org/>`_ sont nécessaires pour la 
construction des paquets.

Pour ces deux scripts, la variable d'environnement ``MAINTAINER`` doit être 
définie ; cette variable permet de renseigner le champ correspondant dans le 
paquet. Elle doit être au format ``"Prénom Nom <adresse électronique>"``. La 
variable ``JOBS`` peut également être définie, elle permet de contrôler le 
nombre de tâches de compilation lancées simultannément (option ``-j`` de la 
commande ``make``).

Pour compiler et créer les paquets dans le répertoire courant en utilisant trois
tâches de compilation simultannées, on aura donc : ::

    export MAINTAINER="John Doe <john@example.com>"
    export JOBS=3
    $HOME/src/medipy/deployement/itk/package.py -j $JOBS 
    $HOME/src/medipy/deployement/package_wrapitk.sh
