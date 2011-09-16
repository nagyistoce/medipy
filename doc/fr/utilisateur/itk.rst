Installation d'ITK et de WrapITK
================================

Paquets
-------

Nous mettons à disposition des paquets pour ITK 3.14 et WrapITK/Python 0.3.0 pour Ubuntu (versions 9.04 et 10.04), à la fois en 32 bits et en 64 bits.
Ces paquets sont disponibles sur le `site Google Code de MediPy <http://code.google.com/p/medipy/downloads/list>`_.

Construction d'un paquet
------------------------

Si aucun paquet n'existe pour votre distribution, nous avons publié deux scripts permettant de créer des paquets ``.deb`` pour ITK et WrapITK.
`CMake <http://www.cmake.org/>`_ et `fakeroot <http://fakeroot.alioth.debian.org/>`_ sont nécessaires pour la construction du paquet.

Pour ces deux scripts, la variable d'environnement ``MAINTAINER`` doit être définie ; cette variable permet de renseigner le champ correspondant dans le paquet.
Elle doit être au format ``"Prénom Nom <adresse électronique>"``.
La variable ``JOBS`` peut également être définie, elle permet de contrôler le nombre de tâches de compilation lancées simultannément (option ``-j`` de la commande ``make``).

Pour compiler et créer les paquets dans le répertoire courant en utilisant trois tâches de compilation simultannées, on aura donc :
::

    MAINTAINER="John Doe <john@example.com>"
    JOBS=3
    $HOME/src/medipy/deployement/package_itk.sh
    $HOME/src/medipy/deployement/package_wrapitk.sh
