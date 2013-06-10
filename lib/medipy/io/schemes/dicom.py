##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

""" Load images from a local filesystem path.

    The path is interpreted as a local filesystem path to either a DICOM file
    or to a directory. If the path is a directory, DICOM files will be searched
    in that path, using recursion unless specified otherwise.
    
    The fragment is of form ``[ tag "=" value { "&" tag "=" value } ]``, where ``tag``
    is a DICOM tag in one of the following forms :
      
      * Numerical (accepted forms : ``"(0020,000e)"``, ``"(0x0020,0x000e)"``, 
        ``"0020000e"``, ``"0x0020000e"``)
      * Named, in which case it must be a key of
        :attr:`medipy.io.dicom.dictionary.name_dictionary`
    
    The tag may also have the following special values :
      
      * ``recursive`` : the value must be ``True`` (default) or ``False``. If ``True``, all 
        files below the given path will be used. If ``False``, only files directly
        below the given path will be used.
    
    A dataset is considered a match for the fragment if it contains all the 
    tags in the filters, and matches all the corresponding values.
    
    The following examples are correct fragments, matching datasets where the
    Patient's Name is equal to ``"Doe^John"`` and the Patient ID is equal to
    ``"1234"``. Since ``recursive`` is not specified for the first three, it is 
    assumed to be ``True``. ::
    
        numerical_fragment = "(0010,0010)=Doe^John&(0010,0020)=1234"
        named_fragment = "patients_name=Doe^John&patient_id=1234"
        mixed_fragment = "00100010=Doe^John&patient_id=1234"
        non_recursive_fragment = "recursive=False&patients_name=Doe^John&patient_id=1234"
"""

import os
import re
import urlparse
import numpy as np

import medipy.base
import medipy.io.dicom

import dicomdir

def load_serie(path, fragment=None) :
    """ Load a serie of images
    """

    datasets = _get_matching_datasets(path, fragment)

    if not datasets :
        return None
    else :
        datasets = medipy.io.dicom.split.images(datasets)
        datasets = medipy.io.dicom.normalize.normalize(datasets)
        stacks = medipy.io.dicom.split.stacks(datasets)
        limages = [medipy.io.dicom.image(stack) for stack in stacks]
        
        # Make sure the images are in their acquisition order 
        limages.sort(key = lambda x:x.metadata.get("acquisition_time", ""))

        return limages

def load(path, fragment=None) :
    """ Load an image.
    """
    
    datasets = _get_matching_datasets(path, fragment)
    
    if not datasets :
        return None
    else :
        return medipy.io.dicom.image(datasets)

def number_of_images(path, fragment=None) :
    """ Return the number of series in given directory.
    """
    
    datasets = _get_matching_datasets(path, fragment)
    return len(medipy.io.dicom.series(datasets))

def _get_filters(fragment) :
    """ Return a list of filters from the URL fragment.
    """
    
    filters = []
    recursive = True
    
    for tag, value in urlparse.parse_qsl(fragment) :
        matches = [re.match(r"^\((?:0x)?([\da-fA-F]{4}),(?:0x)?([\da-fA-F]{4})\)$", tag),
                   re.match(r"^(?:0x)?([\da-fA-F]+)$", tag)]
        if matches[0] :
            tag = medipy.io.dicom.Tag(int(matches[0].group(1), 16),
                                      int(matches[0].group(2), 16))
        elif matches[1] :
            tag = medipy.io.dicom.Tag(int(matches[1].group(1), 16))
        elif tag == "recursive" :
            try :
                recursive = bool(value)
            except ValueError :
                raise medipy.base.Exception("Invalid recursive value : \"{0}\"".format(value))
        else :
            try :
                tag = medipy.io.dicom.Tag(medipy.io.dicom.dictionary.name_dictionary[tag])
            except KeyError :
                raise medipy.base.Exception("No such DICOM tag : \"{0}\"".format(tag))
        filters.append((tag, value))
    
    return filters, recursive

def _get_matching_datasets(path, fragment) :
    """ Return the datasets matching the filters defined in fragment
    """
    
    filters, recursive = _get_filters(fragment)
    
    datasets = []
    
    filenames = []
    if os.path.isdir(path) :
        if recursive :
            for dirpath, dirnames, local_filenames in os.walk(path) :
                filenames.extend([os.path.join(dirpath, x) for x in local_filenames])
        else :
            filenames = [os.path.join(path, x) for x in os.listdir(path)]
    elif os.path.isfile(path) :
        filenames = [path]
    else :
        raise medipy.base.Exception("Cannot find any dataset in \"{0}\"".format(path))
    
    for filename in filenames :
        if not medipy.io.dicom.can_read(str(filename)) :
            continue
        dataset = medipy.io.dicom.read(str(filename))
        
        match = True
        for tag, value in filters :
            if tag not in dataset or dataset[tag].value != value :
                match = False
                break
        
        if match :
            datasets.append(dataset)
    
    return datasets
