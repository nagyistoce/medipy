Description de l'interface graphique d'une fonction
===================================================

Les fonctions présentes dans le menu de MediPy vont avoir une interface 
graphique associée. L'ensemble des boîtes de dialogue de Medimax 3 (créees par 
les macros GET_INT, GET_FLOAT, GET_QCM, ...) sont ici remplacées par une unique
interface graphique.

Cette interface est décrite dans la docstring des fonctions python, selon la 
syntaxe précisée ci-dessous. Pour une fonction de seuillage prenant comme 
paramètres une image d'entrée, les bornes du seuillage et une image de sortie, 
on obtiendra l'interface suivante :

.. image:: threshold2_gui.png


Syntaxe
-------

La syntaxe à utiliser dans la docstring de la fonction est celle utilisée 
notamment par `Sphinx <http://sphinx.pocoo.org/>`_ : il s'agit du format RST, 
simple et relativement proche d'une syntaxe de Wiki. Les éléments de l'interface
graphique sont recensés dans un `champ
<http://docutils.sourceforge.net/docs/ref/rst/restructuredtext.html#field-lists>`_,
nommé ``gui``. Chaque élément est une `définition
<http://docutils.sourceforge.net/docs/ref/rst/restructuredtext.html#definition-lists>`_,
avec :

* terme : nom du paramètre correspondant de la fonction
* premier classificateur : type du contrôle relié à ce paramètre
* second classificateur (optionnel) : paramètre d'initialisation du contrôle
* premier paragraphe de la définition : label associé
* second paragraphe de la définition (optionnel) : tooltip associé

D'autres éléments peuvent bien sur être présents dans la docstring, seul le 
champ ``gui`` est pris en compte par le constructeur d'interface.

**Attention** : la définition ne peut pas être vide selon la syntaxe de RST. 
Avec une définition vide, l'élément ne sera pas affiché.

Exemples
^^^^^^^^

Des exemples valant mieux qu'un long discours, voici l'exemple d'une fonction de
seuillage : ::

    def threshold(input=None, min=None, max=None, output=None) :
        """
        Threshold the input image and store the result in the output image.
        
        :gui: 
            input : Image
                Input image
            
            min : Float
                Minimum value
                
                Minimum value in threshold (included)
                
            max : Float
                Maximum value
                
                Maximum value in threshold (included)
                
            output : Image
                Output image
                
                Output image. If input and output are the same, in-place processing
                is used
        """
        pass

Le champ ``input`` crée un contrôle de type ``Image``, qui va permettre de
choisir entre plusieurs images grâce à des radio-buttons. Les champs ``min`` et 
``max`` ont une syntaxe similaire, mais ajoutent un tooltip. L'interface générée
est la suivante : 

.. image:: threshold2_gui.png

La fonction suivante montre l'utilisation du paramètre d'initialisation du
contrôle : ::

    def binary_erosion(input, seType, seSize, output) :
        """
        Erode the input image and store the result in output.

        :gui:
            input : Image
                Input

            seType : Enum : ("box", "ball")
                Shape

                Type of structuring element

            seSize : Int : 1
                Size

                Must be greater than 1

            output : Image : output=True
                Output
        """
        if input.shape != output.shape :
           output.data = numpy.ndarray(shape=input.shape, dtype=input.dtype)
        se = structuring_element(seType, seSize, len(input.shape))
        scipy.ndimage.binary_erosion(input, se, output=output.data)

Le paramètre d'initialisation de ``seSize`` est optionnel : un contrôle de type
``int`` peut être initialisé sans paramètre, mais on souhaite un élément
structurant de taille au moins 1. Le paramètre d'initialisation de seType est en
revanche obligatoire : on ne peut pas initialiser une liste de choix sans lui
donner la liste des valeurs qu'elle peut prendre. L'interface générée par cette
docstring est la suivante : 

.. image:: binary_erosion.png

On peut également spécifier un intervalle pour les valeurs prises pour un
contrôle de type ``Int`` grâce au paramètre ``range``. A la création de 
l'interface graphique, si le paramètre ``range`` est renseigné, on ajoute alors
une barre de défilement au contrôle. En règle générale, on peut également
choisir de stocker le résultat du traitement dans une nouvelle image en cochant
le bouton ``new`` du paramètre ``Output`` de l'interface graphique. On peut
également spécifier l'intervalle en fonction d'un des autres éléments de
l'interface graphique comme dans l'exemple ci dessous : ::

    def something(input, value, output) :
        """ Do something
            :gui:
                input : Image
                    Input
                value : Int : range = (${input}.data.min(), ${input}.data.max())
                    Value
                output : Image : output = True
                    Output
        """
    
        pass

.. image:: something.png
