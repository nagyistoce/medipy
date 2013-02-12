:numbered:

Utilisation de Torque
=====================

`Torque <http://http://www.adaptivecomputing.com/products/open-source/torque/>`_
est un gestionnaire de ressources permettant de distribuer des calculs sur 
plusieurs nœuds : au lieu de choisir un serveur de calcul peu chargé, 
l'utilisateur soumet une tâche à un nœud maître. Ce nœud maître maintient une
liste de nœuds de calcul, et répartit les tâches en fonction de l'utilisation de
ces nœuds de calcul.

Le plugin :mod:`medipy.torque` encapsule l'opération de soumission des tâches et
facilite certains aspects de cette soumission. Les différentes fonctions de ce
plugin supposent que la tâche à effectuer est un fichier exécutable.

Tâche unique, sans argument
---------------------------

::

    import medipy.io
    image = medipy.io.load("/some/where/image.nii.gz")
    medipy.io.save(image, "/some/where/other_image.nii.gz")

Tâche unique, avec arguments
----------------------------

Tâches multiples
----------------

Il est bien entendu possible de soumettre plusieurs tâches, chacune avec ses 
arguments propres : ::

    values = range(10)
    for value in values :
        medipy.torque.submit(job_script, {"value":value}, working_directory)

Cette méthode n'est cependant pas optimale, car les différentes tâches ne sont
pas reliées entre elles. Torque permet de soumettre un tableau de tâches
(`job array <http://www.clusterresources.com/torquedocs21/2.1jobsubmission.shtml#jobarrays>`_),
qui donne l'avantage de pouvoir référencer un ensemble de tâches. Cette 
fonctionalité est également disponible avec la fonction ``medipy.torque.submit_array`` : ::

    values = range(10)
    medipy.torque.submit_array(job_script, [{"value" : values} for value in values], working_directory)