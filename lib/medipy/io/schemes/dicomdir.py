##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

""" Load images from a DICOMDIR.

    The path is interpreted as a local filesystem path to a DICOMDIR.
    
    The fragment is of form ``[ tag "=" value { "&" tag "=" value } ]``, where ``tag``
    is a DICOM tag in one of the following forms :
      
      * Numerical (accepted forms : ``"(0020,000e)"``, ``"(0x0020,0x000e)"``, 
        ``"0020000e"``, ``"0x0020000e"``)
      * Named, in which case it must be a key of
        :attr:`medipy.io.dicom.dictionary.name_dictionary`
    
    A record in the Directory Record Sequence is considered a match for the 
    fragment if it (or one of its ancestors) contains all the tags in the 
    filters, and matches all the corresponding values. 
    
    See :mod:`medipy.io.schemes.dicom` for fragment examples.
"""

import re
import urlparse

import medipy.base
import medipy.io.dicom

def load_serie(path, fragment=None) :
    """ Load a serie of images
    """

    datasets = _get_matching_datasets(path, fragment)
    datasets = medipy.io.dicom.load_dicomdir_records(datasets)

    if not datasets :
        return None
    else :
        image_datasets = medipy.io.dicom.split.images(datasets)
        normalized_datasets = medipy.io.dicom.normalize.normalize(image_datasets)
        stacks = medipy.io.dicom.split.stacks(normalized_datasets)
        images = [medipy.io.dicom.image(stack) for stack in stacks]
        
        # Make sure the images are in their acquisition order 
        images.sort(key = lambda x:x.metadata.get("acquisition_time", ""))

        return images

def load(path, fragment=None) :
    """ Load an image.
    """
    
    datasets = _get_matching_datasets(path, fragment)
    
    if not datasets :
        return None
    else :
        datasets = medipy.io.dicom.load_dicomdir_records(datasets)
        return medipy.io.dicom.image(datasets)

def number_of_images(path, fragment=None) :
    """ Return the number of series in given DICOMDIR.
    """
    
    datasets = _get_matching_datasets(path, fragment)
    return len(medipy.io.dicom.series(datasets))

def _get_filters(fragment) :
    """ Return a list of filters from the URL fragment.
    """
    
    filters = []
    for tag, value in urlparse.parse_qsl(fragment) :
        matches = [re.match(r"^\((?:0x)?([\da-fA-F]{4}),(?:0x)?([\da-fA-F]{4})\)$", tag),
                   re.match(r"^(?:0x)?([\da-fA-F]+)$", tag)]
        if matches[0] :
            tag = medipy.io.dicom.Tag(int(matches[0].group(1), 16),
                                      int(matches[0].group(2), 16))
        elif matches[1] :
            tag = medipy.io.dicom.Tag(int(matches[1].group(1), 16))
        else :
            try :
                tag = medipy.io.dicom.Tag(medipy.io.dicom.dictionary.name_dictionary[tag])
            except KeyError :
                raise medipy.base.Exception("No such DICOM tag : \"{0}\"".format(tag))
        filters.append((tag, value))
    
    return filters

def _get_matching_datasets(path, fragment) :
    """ Return the datasets matching the filters defined in fragment
    """
    
    filters = _get_filters(fragment)
    
    dicomdir = medipy.io.dicom.read(path)
    datasets = []
    
    for record in dicomdir.directory_record_sequence.value :
        match = True
        for tag, value in filters :
            if tag not in record or record[tag].value != value :
                match = False
                break
        
        if match :
            queue = [record]
            while queue :
                dataset = queue.pop(0)
                datasets.append(dataset)
                for child in dataset.children :
                    queue.append(child)
    
    return datasets
