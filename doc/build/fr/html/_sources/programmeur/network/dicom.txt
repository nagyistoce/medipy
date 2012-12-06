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

Exemple : ::
    
    connection = medipy.network.dicom.scu.Connection(
        "pacs.example.com", 11112, 
        "MY_MACHINE", "REMOTE_PACS")

SCU : clients des services DICOM
--------------------------------

Le package ``medipy.network.dicom.scu`` contient, outre la classe ``Connection``
des classes encapsulant les différents SCU (Service Class User, équivalents au
terme « client » dans le sens réseau du terme). Ces classes sont toutes des 
objets-fonction. 

Echo
^^^^

Le SCU le plus simple est Echo, équivalent du niveau DICOM d'un ping. Aucun
paramètre n'est à fournir, puisque les AE concernés sont déjà spécifiés dans la
connexion. Ce SCU n'a pas de valeur de retour, mais lance une exception si la
connexion n'a pas pu être établie. ::

    echo = medipy.network.dicom.scu.Echo(connection)
    try :
        echo()
    except medipy.base.Exception, e:
        print "Service injoignable : {0}".format(e)
    else :
        print "OK"

Find
^^^^

Le SCU Find permet de faire des requêtes simples sur un nœud DICOM. Ce SCU a
trois paramètres : l'objet de plus haut niveau (``patient`` ou ``study``), le
niveau auquel se fait la requête (``patient``, ``study``, ``series`` ou ``image``)
et les critères de recherche, contenus dans un DataSet. Les attributs présents
dans ce DataSet peuvent contenir les valeurs suivantes :

* une valeur simple (pas de joker) : ce critère sélectionne les DataSet ayant
  un attribut avec *exactement* cette valeur. Pour les chaînes de caractères,
  la recherche dépend de la casse, à l'exception du type PN.
* une liste d'UID : il s'agit d'une chaîne de caractères contenant les 
  différents UID, séparés par un backslash (``\``). Ce critère sélectionne les
  DataSet ayant un attribut dont la valeur est *contenue* dans la liste.
* une valeur vide (``None`` ou ``""``). Ce critère sélectionne *tous* les 
  DataSet, il est utilisé afin de spécifier un attribut qu'on souhaite récupérer.
* une chaîne de caractères avec des jokers (``*`` et ``?``). Ces deux caractères
  gardent leur signification habituelle (``*`` pour signifier « n'importe quelle
  chaîne » et ``?`` pour signifier « n'importe quel caractère »).
* une séquence : la séquence devra contenir un seul élément, et les critères
  seront appliqué à cet élément selon les règles ci-dessus.

L'exemple d'utilisation suivant retourne un DataSet pour chaque patient dont le
nom commence par L. Chaque DataSet résultat contiendra l'attribut Patient Id. ::

    query = medipy.io.dicom.DataSet()
    query.patients_name = "L*"
    query.patient_id = None
    find = medipy.network.dicom.scu.Find(connection, "patient", "patient", query)
    results = find()   


Ce type de requête n'est pas récursif : on ne peut pas faire de requète sur 
l'attribut Study Description au niveau ``patient``. De plus, les attributs-clés 
des niveaux supérieurs doivent être renseignés :

* Patient ID (0010,0020) pour le niveau ``patient``
* Study Instance UID (0020,000D) pour le niveau ``study``
* Study Instance UID (0020,000E) pour le niveau ``series``
* SOP Instance UID (0008,0018) pour le niveau ``image``
    

Requêtes de haut niveau
-----------------------

Les fonctions contenues dans le module ``medipy.network.dicom.query`` permettent
de faire des requêtes de façon plus simple qu'en passant par le SCU Find. 
Cependant, ces fonctions sont potentiellement plus lentes.

Requête relationnelle
^^^^^^^^^^^^^^^^^^^^^

La fonction ``medipy.network.dicom.query.relational`` permet de faire une 
requête relationnelle, selon la définition de DICOM (PS 3.4-2011, C.4.1.2.2.1 et
C.4.1.3.2.2). Cette fonction permet donc de faire une recherche à un niveau de
départ, et d'étendre la recherche aux niveaux inférieurs jusqu'à que tous les
critères soient remplis. Si des critères ne sont pas remplis au niveau le plus
bas, une liste vide est retournée.

Dans l'exemple suivant, les résultats de la requête auront un attribut 
Query/Retrieve Level (0008,0052) valant "STUDY", puisqu'un des attributs de la
requête porte sur l'entité Study. Les attributs-clés du niveau STUDY et des
niveaux supérieurs (Study Instance UID (0020,000D) et Patient ID (0010,0020))
seront présents dans les résultats. ::

    query = medipy.io.dicom.DataSet()
    query.study_description = "PROTOCOLES^MON_PROTOCOLE"
    
    datasets = medipy.network.dicom.query.relational(connection, "patient", "patient", query) 

Accès web aux objets DICOM (WADO)
---------------------------------

Hors des services permettant de copier des données (Get et Move), DICOM propose
aussi un accès aux données par HTTP : WADO pour *Web Access to DICOM Objects*
(PS 3.18-2011).
Ce service permet de récupérer des données, au format DICOM, par une requète
HTTP de type GET qui spécifie la clé de l'objet à récupérer. Cette clé est formée
par les mêmes attributs pour ceux utilisés par le SCU Find. Le module 
``medipy.network.dicom.wado`` encapsule ce service : ::

    wado_url = "http://pacs.example.com/wado"
    
    query = medipy.io.dicom.DataSet()
    query.patient_id = "12345"
    query.study_instance_uid = "..."
    query.series_instance_uid = "..."
    query.sop_instance_uid = "..."
    
    # Stockage du résultat dans un fichier
    medipy.network.dicom.wado.get(wado_url, query, "result.dcm")
    
    # Retour du DataSet résultat
    dataset = medipy.network.dicom.wado.get(wado_url, query)

.. _DCMTK: http://dicom.offis.de/dcmtk.php.en