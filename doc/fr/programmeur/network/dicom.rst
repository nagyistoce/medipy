Fonctionnalités réseau de DICOM
===============================

Afin d'exploiter les fonctionnalités réseau de DICOM, MediPy s'appuie sur DCMTK_.
Des paquets sont disponibles dans les distributions Linux supportées par MediPy,
et des binaires Windows sont disponibles sur le site du projet.

Informations de connexion
-------------------------

La classe ``medipy.network.dicom.scu.Connection`` contient les différentes 
informations pour se connecter à un nœud DICOM sur le réseau. Il s'agit pour la
couche TCP/IP de l'adresse de la machine et du port, et, pour la couche DICOM,
des noms des deux entités applicatives (Application Entity title) concernées
par la connexion.

    * :class:`~medipy.network.dicom.Connection`
    * :class:`~medipy.network.dicom.SSHTunnelConnection`

SCU : clients des services DICOM
--------------------------------

Le package ``medipy.network.dicom.scu`` contient, outre la classe ``Connection``
des classes encapsulant les différents SCU (Service Class User, équivalents au
terme « client » dans le sens réseau du terme). Ces classes sont toutes des 
objets-fonction. 

    * :class:`~medipy.network.dicom.scu.Echo`
    * :class:`~medipy.network.dicom.scu.Find`
    * :class:`~medipy.network.dicom.scu.Store`


Requêtes de haut niveau
-----------------------

Les fonctions contenues dans le module ``medipy.network.dicom.query`` permettent
de faire des requêtes de façon plus simple qu'en passant par le SCU Find. 
Cependant, ces fonctions sont potentiellement plus lentes.

    * :func:`medipy.network.dicom.query.relational`

Accès web aux objets DICOM (WADO)
---------------------------------

Hors des services permettant de copier des données (Get et Move), DICOM propose
aussi un accès aux données par HTTP : WADO pour *Web Access to DICOM Objects*
(PS 3.18-2011).

    * :func:`medipy.network.dicom.wado.get`

.. _DCMTK: http://dicom.offis.de/dcmtk.php.en