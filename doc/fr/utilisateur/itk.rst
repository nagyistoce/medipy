Installation d'ITK
==================

Paquets
-------

Nous mettons à disposition des paquets de ITK 4.6.1 pour Ubuntu (12.04 et 
14.04, 32 bits et 64 bits) et Debian (7.0, 32 bits et 64 bits). Ces paquets sont 
disponibles sur le dépôt de paquets de MediPy, à ajouter par les commandes
suivantes : ::

    wget -O - https://medipy.u-strasbg.fr/packages/medipy.gpg.key | sudo apt-key add -
    add-apt-repository http://medipy.u-strasbg.fr/packages
    apt-get update

ITK est alors installable comme tout autre paquet : ::

    apt-get install insighttoolkit
