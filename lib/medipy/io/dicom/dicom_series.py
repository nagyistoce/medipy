##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import csv
import logging
import medipy.base

class DicomSeries(object) :
    """ A dictionary of DICOM series, mapping a Series Instance UID to the
        Series Description, and optionally to a custom name given by the user.
    """
    
    def __init__(self, root=None, series=None) :
        self.root = root
        self.series = []
        for serie in series :
            if len(serie) == 3 :
                self.series.append(serie)
            elif len(serie) == 2 :
                self.series.append(list(serie)+[""])
            elif len(serie) == 1 :
                self.series.append(list(serie)+["", ""])
    
    def has_uid(self, uid) :
        """ Test whether the Series Instance UID is in the dictionary.
        """
        
        matches = [x for x in self.series if x[0] == uid]
        return (len(matches) != 0)

    def custom_name_from_uid(self, uid) :
        """ Return the custom name associated with the Series Instance UID
        """
        
        if not self.has_uid(uid) :
            raise medipy.base.Exception("No such UID {0}".format(repr(uid)))
        matches = [x for x in self.series if x[0] == uid]
        if len(matches) > 1 :
            logging.warning("UID is not unique : {0}".format(uid))
        return matches[0][2]

    def url_from_uid(self, uid) :
        """ Return the URL associated with the Series Instance UID
        """
        
        if not self.has_uid(uid) :
            raise medipy.base.Exception("No such UID {0}".format(repr(uid)))
        matches = [x for x in self.series if x[0] == uid]
        if len(matches) > 1 :
            logging.warning("UID is not unique : {0}".format(uid))
        return "{0}#series_instance_uid={1}".format(self.root, matches[0][0])

    def has_description(self, description) :
        """ Test whether the Series Description is in the dictionary.
        """
        
        matches = [x for x in self.series if x[1] == description]
        return (len(matches) != 0)

    def url_from_description(self, description) :
        """ Return the URL associated with the Series Description.
        """
        
        if not self.has_description(description) :
            raise medipy.base.Exception("No such description {0}".format(repr(description)))
        matches = [x for x in self.series if x[1] == description]
        if len(matches) > 1 :
            logging.warning("Description is not unique : {0}".format(description))
        return "{0}#series_instance_uid={1}".format(self.root, matches[0][0])

    def has_custom_name(self, custom_name) :
        """ Test whether the custom name is in the dictionary.
        """
        
        matches = [x for x in self.series if x[2] == custom_name]
        return (len(matches) != 0)

    def url_from_custom_name(self, custom_name) :
        """ Return the URL associated with the custom name.
        """
        
        if not self.has_custom_name(custom_name) :
            raise medipy.base.Exception("No such custom name {0}".format(repr(custom_name)))
        matches = [x for x in self.series if x[2] == custom_name]
        if len(matches) > 1 :
            logging.warning("Custom name is not unique : {0}".format(custom_name))
        return "{0}#series_instance_uid={1}".format(self.root, matches[0][0])

    def write(self, filename):
        """ Write the dictionary to a file.
        
            The dictionary is saved in CSV, using the following format :
            
            * The first line contains only the root URL
            * Each successive line contain the UID, description and custom
              name of each series.
        """
        
        fd = open(filename, "w")
        
        fd.write("{0}\n".format(self.root))
        # TODO : delimiter & quoting
        writer = csv.writer(fd)
        for uid, description, custom_name in self.series :
            writer.writerow(
                [uid, description.encode("ascii", "replace"), custom_name])
        fd.close()

    @staticmethod
    def read(filename) :
        """ Read a dictionary from a file.
            
            See :func:`write` for the file format.
        """
        
        fd = open(filename)
        root = fd.readline().strip()
        reader = csv.reader(fd)
        series = [x for x in reader]
        fd.close()

        return DicomSeries(root, series)
