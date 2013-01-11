Encapsulation d'un document au format DICOM
===========================================

Le standard DICOM permet de stocker tout type de document dans un Data Set [#]_.
Le module ``medipy.io.dicom.encapsulated_document`` permet d'encapsuler un 
fichier dans un Data Set, et d'exposer un fichier précédemment encapsulé. Les
fonctions contenues dans ce module ne génèrent pas un Data Set directement 
envoyable sur un nœud DICOM, seuls les modules DICOM concernant l'objet
Encapsulated Document sont présents. L'exemple suivant permet d'obtenir un
Data Set « correct » pour le stockage. 

::

    dataset = medipy.io.dicom.encapsulated_document.encapsulate(filename)
    
    # Use Raw Data Storage as a SOP class since Encapsulated Document SOP classes
    # are limited to PDF and HL7 documents
    dataset.header.media_storage_sop_class_uid = medipy.io.dicom.UI(
        medipy.io.dicom.dictionary.uid_name_dictionary["raw_data_storage"])
    dataset.sop_class_uid = medipy.io.dicom.UI(
        medipy.io.dicom.dictionary.uid_name_dictionary["raw_data_storage"])
    
    # New file, hence new SOP Instance UID
    dataset.sop_instance_uid = medipy.io.dicom.generate_uid()
    
    # Most DICOM nodes will complain without Patient information
    dataset.patient_id = medipy.io.dicom.LO(medipy.io.dicom.generate_uid(False))
    dataset.patients_name = medipy.io.dicom.PN(dataset.patient_id.value)
    
    # Generate a new Study Instance UID
    dataset.study_instance_uid = medipy.io.dicom.generate_uid()
    
    connection = medipy.network.dicom.Connection(host, port, aet, aec)
    store = medipy.network.dicom.scu.Store(connection, dataset)
    store()

.. [#] PS 3.3-2011, A.1.2.16, p. 121 et PS 3.3-2011, A.45, p. 278 