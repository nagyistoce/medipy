Menu
====

Deux méthodes de génération du menu sont possibles, l'une automatique, l'autre 
manuelle :

* la méthode automatique permet de générer un menu complet de l'ensemble des 
  fonctions rendues publiques, mais n'est pas flexible (hiérarchie et ordre de 
  présentation fixes, peu de possibilité de renommage d'une fonction)
* la méthode manuelle se base sur un fichier texte contenant une description du
  menu ; cette méthode est plus flexible, mais nécessite plus d'intervention de
  la part de l'utilisateur

Génération automatique
----------------------

La méthode automatique de génération du menu se base sur des fichiers nommés 
``api.py``. Ces fichiers, présents dans chaque paquet (i.e. à chaque niveau de
l'arborescence), fonctionnent à la manière des fichiers ``__init__.py`` en
important certains noms des différents modules du paquet. Chacun de ces noms
donnera un item du menu. Soit la hiérarchie de fichiers suivante : ::

    components
    |- arithmetic
    |  |- __init__.py
    |  |- api.py
    |  \- basic.py (contient les fonctions add et subtract)
    |
    \- segmentation
       |- ...

Le fichier ``api.py`` contient l'ensemble des noms qui apparaîtront dans le 
menu, sous l'entrée ``arithmetic``. Pour y ajouter les fonctions ``add`` et
``subtract``, ce fichier contiendra donc : ::

    from basic import add, subtract

La forme ``from ... import ... as ...`` permet de changer le nom que la fonction
aura dans le menu. L'exemple suivant fera apparaître *Addition* et 
*Soustraction* dans le menu, ces fonctions étant respectivement reliées à 
``arithmetic.add`` et ``arithmetic.subtract``. ::

    from basic import add as Addition
    from basic import subtract as Soustraction

La génération automatique est gérée par la fonction
`medipy.gui.menu_builder.from_api.build_menu`

Génération manuelle
-------------------

La génération manuelle passe par un fichier texte décrivant le menu. Chaque 
ligne du fichier décrit une entrée du menu par la syntaxe suivante : ::

    [indentation]<nom_python>[nom_à_afficher]

L'indentation permet de définir la hiérarchie des entrées, elle doit être 
cohérente d'une entrée à l'autre. Il est donc recommandé d'utiliser la même
convention que pour les codes sources Python, soit quatre espaces par niveau
d'indentation.

Le nom Python donne un nom de module ou de fonction défini soit de façon absolue
(si ce nom commence par ``medipy``), soit de façon relative à la 
hiérarchie des noms de entrées précédentes. Dans l'exemple suivant, l'entrée 
``bet`` correspond à la fonction ``medipy.segmentation.brain.bet``, 
et l'entrée ``medipy.other_module.skull_segmentation`` donne 
explicitement le module de la fonction ``skull_segmentation`` : ::
    
    segmentation
        brain
            bet
        skull
            medipy.other_module.skull_segmentation

Si le nom à afficher est vide, les caractères ``_`` sont remplacés par des 
espaces dans le composant le plus à droite du nom Python.

L'utilisation d'un nom de fonction ou de module inexistant permet une grande 
flexibilité : le menu suivant donne l'image ci-dessous. ::

    segmentation
        brain
            bet Brain Extraction Tool
    none --- (work in progress) ---
        belghith
        grigis
        irmd_work

.. image:: menu.png

Structure de données
--------------------

Chaque item de l'arbre est représenté par un triplet (nom, fonction, enfants) :

* nom : nom à afficher
* fonction : fonction à exécuter ou None ; la fonction doit contenir dans sa
  docstring une :doc:`description de l'interface graphique <function_gui>` 
* enfants : liste de triplets
